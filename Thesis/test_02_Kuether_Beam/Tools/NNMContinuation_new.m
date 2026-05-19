classdef NNMContinuation_new < handle
    properties
        ReducedAssembly
        spring
        P
        SystemNLvib

        X_out
        energies
        frequencies
        rms_displacements 
        modal_amplitudes 
        inorm
    end

    methods
        function obj = NNMContinuation_new(reduced_assembly, spring_obj, P_matrix)
            obj.ReducedAssembly = reduced_assembly;
            obj.spring = spring_obj;
            obj.P = P_matrix;

            % use wrapper for nonlinear force
            obj.ReducedAssembly.DATA.fnl_CUSTOM = @(q) obj.nl_force_wrapper(q);

            n_dofs = size(P_matrix, 2);
            F_ext_null = zeros(n_dofs, 1);
            obj.SystemNLvib = FE_system(obj.ReducedAssembly, F_ext_null, 'custom');
        end

        function [Kt, f_nl] = nl_force_wrapper(obj, q)
            u_full = obj.P * q;
            u_nl = u_full(obj.spring.dof_idx);

            f_full = zeros(size(u_full, 1), 1);
            K_full = sparse(size(u_full, 1), size(u_full, 1));

            % --- REGULARIZED CONTACT (replaces hard if/else) ---
            % The smoothing width eps_reg is ~1% of the gap. This avoids
            % the discontinuous Jacobian that destabilizes Newton-Raphson
            % near the tongue bifurcation points.
            g = u_nl - obj.spring.a;
            eps_reg = obj.spring.a * 0.01;
            if eps_reg == 0, eps_reg = 1e-8; end

            % Smooth ramp: f = k/2 * (g + sqrt(g^2 + eps^2))
            % Recovers the exact piecewise-linear law for |g| >> eps_reg
            sqrt_term = sqrt(g^2 + eps_reg^2);
            f_val = (obj.spring.k / 2) * (g + sqrt_term);
            k_val = (obj.spring.k / 2) * (1 + g / sqrt_term);

            f_full(obj.spring.dof_idx) = f_val;
            K_full(obj.spring.dof_idx, obj.spring.dof_idx) = k_val;

            f_nl = obj.P' * f_full;
            Kt   = obj.P' * K_full * obj.P;
        end


        function solve(obj, mode_idx, log10a_start, log10a_end, solver_options)
            if nargin < 5
                solver_options = struct();
            end

            % --- Default solver options ---
            if ~isfield(solver_options, 'Ntd'),           solver_options.Ntd = 1000; end
            if ~isfield(solver_options, 'ds'),            solver_options.ds = 1e-3;  end
            if ~isfield(solver_options, 'stepmax'),       solver_options.stepmax = 30000; end
            if ~isfield(solver_options, 'dynamicDscale'), solver_options.dynamicDscale = 1; end

            % FIX 1: dsmax is now capped at 5*ds, not 100*ds.
            % A large dsmax is the primary cause of the solver jumping across
            % tongue bifurcations and then getting stuck.
            if ~isfield(solver_options, 'dsmax')
                solver_options.dsmax = solver_options.ds * 5;
            end
            if ~isfield(solver_options, 'dsmin')
                solver_options.dsmin = solver_options.ds / 1e5;
            end

            Ntd = solver_options.Ntd;
            ds  = solver_options.ds;

            M_r = obj.ReducedAssembly.DATA.M;
            K_r = obj.ReducedAssembly.DATA.K;
            [Phi_red, Om2_red] = eig(K_r, M_r);
            [om_lin_red, sort_idx] = sort(sqrt(diag(Om2_red)));
            Phi_red = Phi_red(:, sort_idx);

            om_start = om_lin_red(mode_idx);
            phi_start = Phi_red(:, mode_idx);

            % --- Setting correct mode orientation (+a) ---
            phi_phys = obj.P * phi_start;
            if phi_phys(obj.spring.dof_idx) < 0
                phi_start = -phi_start;
            end

            % Normalization
            [~, inorm] = max(abs(phi_start));
            phi_start = phi_start / phi_start(inorm);

            m = size(M_r, 1);
            nnorm = setdiff(1:2*m, [inorm, inorm+m]);
            ys0 = [phi_start; zeros(m, 1)];
            a0 = ys0(inorm);
            x0 = [ys0(nnorm)/a0; om_start; 0];

            % --- Shooting parameters ---
            Np   = 1;
            qscl = obj.spring.a;
            fscl = mean(diag(K_r)) * qscl;

            % FIX 2: Physics-informed Dscale.
            % Instead of all-ones, we scale the state coordinates by the
            % expected modal amplitude ratio (M-orthonormal eigenvectors),
            % and the frequency coordinate by om_start.
            % This keeps the arclength metric isotropic near tongue points,
            % where subsidiary modes grow large relative to the primary one.
            state_scale = zeros(length(x0) - 2, 1);
            % Displacement part of nnorm entries (indices 1..m in ys, minus inorm)
            % All start at amplitude ~1 (normalized), so use qscl as baseline
            state_scale(:) = qscl;
            dscale = [state_scale; om_start; qscl; 1e0];

            Sopt = struct(...
                'Dscale',          dscale, ...
                'dynamicDscale',   solver_options.dynamicDscale, ...
                'stepmax',         solver_options.stepmax, ...
                'dsmin',           solver_options.dsmin, ...
                'dsmax',           solver_options.dsmax);

            fprintf('Starting Computation NNM %d (Shooting, Ntd = %d, ds = %g, dsmax = %g)...\n', ...
                mode_idx, Ntd, ds, solver_options.dsmax);

            obj.X_out = solve_and_continue(x0, ...
                @(X) shooting_residual(X, obj.SystemNLvib, Ntd, Np, 'NMA', qscl, fscl, inorm), ...
                log10a_start, log10a_end, ds, Sopt);
            obj.inorm = inorm;
            obj.compute_energies(inorm);
        end

        function compute_energies(obj, inorm)
            m = size(obj.P, 2);
            om_sh = obj.X_out(end-2, :);
            log10a_sh = obj.X_out(end, :);
            a_sh = 10.^log10a_sh;

            nnorm = setdiff(1:2*m, [inorm, inorm+m]);

            obj.frequencies = om_sh / (2*pi);
            obj.energies = zeros(1, size(obj.X_out, 2));
            obj.modal_amplitudes = zeros(m, size(obj.X_out, 2));
            obj.rms_displacements = zeros(1, size(obj.X_out, 2));

            M_r = obj.ReducedAssembly.DATA.M;
            K_r = obj.ReducedAssembly.DATA.K;

            dir_mult = sign(obj.P(obj.spring.dof_idx, 1));
            if dir_mult == 0, dir_mult = 1; end

            for i = 1:size(obj.X_out, 2)
                ai = a_sh(i);
                ys_nnorm_i = obj.X_out(1:end-3, i) * ai;

                ys = zeros(2*m, 1);
                ys(inorm) = ai;
                ys(nnorm) = ys_nnorm_i;

                qi = ys(1:m);
                vi = ys(m+1:end);

                obj.modal_amplitudes(:, i) = abs(qi);
                x_full = obj.P * qi;
                obj.rms_displacements(i) = sqrt(mean(x_full.^2));

                E_lin = 0.5 * vi' * M_r * vi + 0.5 * qi' * K_r * qi;

                u_nl = x_full(obj.spring.dof_idx) * dir_mult;
                E_nl = 0;
                if u_nl > obj.spring.a
                    E_nl = 0.5 * obj.spring.k * (u_nl - obj.spring.a)^2;
                end

                obj.energies(i) = E_lin + E_nl;
            end
        end

        function plot_backbone(obj, style_str, display_name)
            if nargin < 2, style_str = 'b-'; end
            if nargin < 3, display_name = 'NNM Shooting'; end

            semilogx(obj.energies, obj.frequencies, style_str, 'LineWidth', 2, 'DisplayName', display_name);
            grid on;
            xlabel('Energy [J]');
            ylabel('Frequency [Hz]');
        end

        function plot_frequency_vs_amplitude(obj, style_str, display_name)
            if nargin < 2, style_str = 'r-'; end
            if nargin < 3, display_name = 'Freq vs Amp'; end
            
            log10a_sh = obj.X_out(end, :);
            a_sh = 10.^log10a_sh;
            
            semilogx(a_sh, obj.frequencies, style_str, 'LineWidth', 2, 'DisplayName', display_name);
            grid on;
            xlabel('Continuation Parameter (Modal Amplitude) [m]');
            ylabel('Frequency [Hz]');
            title('Frequency vs Modal Amplitude');
        end

        function plot_time_history_and_shape(obj, target_energy, x_coords)
            [~, idx] = min(abs(obj.energies - target_energy));
            actual_energy = obj.energies(idx);
            om_sh = obj.frequencies(idx) * 2 * pi;
            
            fprintf('NNM at E = %.2e J (Target: %.2e J)\n', actual_energy, target_energy);
            
            m = size(obj.P, 2);
            log10a_sh = obj.X_out(end, idx);
            ai = 10^log10a_sh;
            ys_nnorm_i = obj.X_out(1:end-3, idx) * ai;
            
            nnorm = setdiff(1:2*m, [obj.inorm, obj.inorm+m]);
            ys = zeros(2*m, 1);
            ys(obj.inorm) = ai;
            ys(nnorm) = ys_nnorm_i;
            
            q0 = ys(1:m);
            v0 = ys(m+1:end);
            
            T = 2 * pi / om_sh;
            t_span = linspace(0, T, 500);
            
            M_r = obj.ReducedAssembly.DATA.M;
            K_r = obj.ReducedAssembly.DATA.K;
            
            odefun = @(t, Y) obj.rom_ode_system(t, Y, M_r, K_r);
            options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
            [t_out, Y_out] = ode45(odefun, t_span, [q0; v0], options);
            
            q_history = Y_out(:, 1:m)';
            u_full_history = obj.P * q_history;
            u_mid_history = u_full_history(obj.spring.dof_idx, :);
            
            u_full_0 = u_full_history(:, 1);
            uy_indices = 2:3:size(u_full_0, 1);
            uy_0 = u_full_0(uy_indices);
            
            figure('Name', sprintf('NNM Dynamics (E = %.2e J)', actual_energy), 'Color', 'w', 'Position', [100, 100, 1200, 500]);
            sgtitle(sprintf('NNM Response at E = %.2e J', actual_energy), 'FontSize', 14, 'FontWeight', 'bold');
            
            subplot(1, 2, 1);
            plot(t_out, u_mid_history * 1000, 'b-', 'LineWidth', 2);
            hold on;
            yline(obj.spring.a * 1000, 'r--', 'Gap Clearance', 'LabelHorizontalAlignment', 'left', 'LineWidth', 1.5);
            grid on;
            xlabel('Time [s]'); ylabel('Midpoint Displacement [mm]');
            title(sprintf('Time History (T = %.4f s)', T));
            legend('Midpoint Trajectory', 'Contact Limit', 'Location', 'best');
            
            subplot(1, 2, 2);
            plot(x_coords, uy_0 * 1000, '-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
            grid on;
            xlabel('Beam X [m]'); ylabel('Deflection [mm]');
            title('Deformed Shape at t = 0');
            hold off;
        end
        
        function dYdt = rom_ode_system(obj, ~, Y, M_r, K_r)
            m = size(obj.P, 2);
            q = Y(1:m);
            v = Y(m+1:end);
            
            u_full = obj.P * q;
            u_nl = u_full(obj.spring.dof_idx);
            
            % Use same regularized contact as nl_force_wrapper for consistency
            g = u_nl - obj.spring.a;
            eps_reg = obj.spring.a * 0.01;
            if eps_reg == 0, eps_reg = 1e-8; end
            f_val = (obj.spring.k / 2) * (g + sqrt(g^2 + eps_reg^2));
            
            f_nl_full = zeros(size(u_full, 1), 1);
            f_nl_full(obj.spring.dof_idx) = f_val;
            f_nl_red = obj.P' * f_nl_full;
            
            qdd = M_r \ (-K_r * q - f_nl_red);
            dYdt = [v; qdd];
        end
    end
end
