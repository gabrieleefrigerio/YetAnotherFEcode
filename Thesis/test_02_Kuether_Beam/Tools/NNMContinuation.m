classdef NNMContinuation < handle
% NNMContinuation  Computes Nonlinear Normal Modes via shooting + continuation.
%
% Replicates the methodology of Kuether, Brake & Allen (2014):
%   "Evaluating Convergence of Reduced Order Models Using Nonlinear Normal Modes"
%
% Two continuation strategies are available:
%   solve()       - uses log10(amplitude) as continuation parameter.
%                   Robust for monotone backbones (hardening systems).
%                   Recommended first choice.
%   solve_omega() - uses omega as continuation parameter with arc-length.
%                   Can follow folds and internal resonances tongues.
%                   Use if solve() gets stuck at a fold.
%
% Usage:
%   nnm = NNMContinuation(reduced_assembly, spring, P);
%   nnm.solve(1, -5, -1);          % mode 1, log10a from -5 to -1
%   nnm.plot_backbone('b-', 'ROM [1 MC]');
%
% Output properties after solve():
%   frequencies        [1 x Npts]  frequency in Hz at each solution point
%   energies           [1 x Npts]  total mechanical energy [J]
%   rms_displacements  [1 x Npts]  spatial RMS of physical displacement [m]
%   modal_amplitudes   [m x Npts]  absolute modal coordinate amplitudes
%   X_out              raw output matrix from solve_and_continue

    properties
        ReducedAssembly     % YAFEC ReducedAssembly object
        spring              % LocalNonlinearity object with fields:
                            %   .k         spring stiffness [N/m]
                            %   .a         clearance [m]
                            %   .dof_idx   physical DOF index of contact
        P                   % reduction basis matrix [N_full x m]
        SystemNLvib         % FE_system object used by NLvib shooting

        % --- Results ---
        X_out               % raw continuation output [n_vars x Npts]
        frequencies         % [1 x Npts]  Hz
        energies            % [1 x Npts]  J
        rms_displacements   % [1 x Npts]  m  (spatial RMS over all DOFs)
        modal_amplitudes    % [m x Npts]  absolute modal amplitudes
    end

    % =====================================================================
    methods
    % =====================================================================

        function obj = NNMContinuation(reduced_assembly, spring_obj, P_matrix)
        % Constructor.
        %   reduced_assembly : YAFEC ReducedAssembly (must have .DATA.M and .DATA.K)
        %   spring_obj       : struct/object with .k, .a, .dof_idx
        %   P_matrix         : reduction basis [N_full x m]

            obj.ReducedAssembly = reduced_assembly;
            obj.spring          = spring_obj;
            obj.P               = P_matrix;

            % Inject the smoothed unilateral spring as a custom nonlinear element
            obj.ReducedAssembly.DATA.fnl_CUSTOM = @(q) obj.nl_force_wrapper(q);

            n_dofs         = size(P_matrix, 2);
            F_ext_null     = zeros(n_dofs, 1);
            obj.SystemNLvib = FE_system(obj.ReducedAssembly, F_ext_null, 'custom');
        end

        % -----------------------------------------------------------------
        function [Kt, f_nl] = nl_force_wrapper(obj, q)
        % Smoothed unilateral spring force and tangent stiffness.
        %
        % Uses a hyperbolic regularisation to smooth the non-smooth contact:
        %   f = (k/2) * (g + sqrt(g^2 + eps^2))
        % where g = u - a  is the penetration.
        % The tangent stiffness is the exact analytical derivative.
        % eps_smooth controls the width of the transition zone.

            u_full = obj.P * q;
            u_nl   = u_full(obj.spring.dof_idx);

            f_full = zeros(size(u_full, 1), 1);
            K_full = sparse(size(u_full, 1), size(u_full, 1));

            g          = u_nl - obj.spring.a;
            eps_smooth = obj.spring.a * 0.10;
            if eps_smooth == 0, eps_smooth = 1e-6; end

            f_val = (obj.spring.k / 2) * (g + sqrt(g^2 + eps_smooth^2));
            k_val = (obj.spring.k / 2) * (1 + g / sqrt(g^2 + eps_smooth^2));

            f_full(obj.spring.dof_idx)                         = f_val;
            K_full(obj.spring.dof_idx, obj.spring.dof_idx)     = k_val;

            f_nl = obj.P' * f_full;
            Kt   = obj.P' * K_full * obj.P;
        end

        % =================================================================
        %  STRATEGY 1: log10(amplitude) as continuation parameter
        %  Robust for hardening backbones without folds.
        % =================================================================
        function solve(obj, mode_idx, log10a_start, log10a_end, varargin)
        % SOLVE  Compute NNM backbone using log10(amplitude) as parameter.
        %
        %   solve(mode_idx, log10a_start, log10a_end)
        %   solve(mode_idx, log10a_start, log10a_end, 'Ntd', 500, 'ds', 1e-5)
        %
        %   mode_idx      : which linear mode to continue from (1 = first)
        %   log10a_start  : starting log10(amplitude), e.g. -5
        %   log10a_end    : ending   log10(amplitude), e.g. -1
        %
        %   Optional name-value pairs:
        %     'Ntd'   number of time steps per period (default 500)
        %     'ds'    initial continuation step in log10a (default 1e-4)

            % --- Parse optional arguments ---
            p = inputParser;
            addParameter(p, 'Ntd', 500,   @isnumeric);
            addParameter(p, 'ds',  1e-4,  @isnumeric);
            parse(p, varargin{:});
            Ntd = p.Results.Ntd;
            ds  = p.Results.ds;

            % --- Linear eigenvalue problem on ROM ---
            [inorm, nnorm, phi_start, om_lin] = ...
                obj.get_linear_mode(mode_idx);
            m = size(obj.P, 2);

            % --- Initial state vector (NLvib NMA format) ---
            % X_nlvib = [ys_nnorm/a (2m-1);  om;  D;  log10a]
            % lambda  = log10a  (last element)
            a0        = 10^log10a_start;
            ys0       = [phi_start; zeros(m, 1)];
            ys_nnorm0 = ys0(nnorm) / a0;    % normalised (divided by a)

            % x0 = free variables (everything except lambda = log10a)
            x0    = [ys_nnorm0; om_lin; 0];   % [ys_nnorm/a; om; D=0]
            lam_s = log10a_start;
            lam_e = log10a_end;

            % --- Scaling ---
            qscl = obj.spring.a;
            fscl = mean(diag(obj.ReducedAssembly.DATA.K)) * qscl;

            % Dscale: scaled variables should all be O(1)
            % [ys_nnorm/a ~ 1;  om ~ om_lin;  D ~ 0 (use small scale);  log10a ~ 1]
            dscale = [ones(length(nnorm), 1); ...
                      om_lin; ...
                      1e-2; ...
                      1.0];

            Sopt = struct(...
                'Dscale',          dscale, ...
                'dynamicDscale',   1, ...
                'stepmax',         1000, ...
                'dsmin',           ds / 1e5, ...
                'dsmax',           ds * 100);

            fprintf('--- NNM Shooting: mode %d, log10a in [%.1f, %.1f] ---\n', ...
                mode_idx, log10a_start, log10a_end);

            obj.X_out = solve_and_continue(x0, ...
                @(X) shooting_residual(X, obj.SystemNLvib, Ntd, 1, ...
                                       'NMA', qscl, fscl, inorm), ...
                lam_s, lam_e, ds, Sopt);

            fprintf('Done: %d solution points.\n', size(obj.X_out, 2));
            obj.compute_energies_loga(inorm, nnorm);
        end

        % =================================================================
        %  STRATEGY 2: omega as continuation parameter (arc-length)
        %  Can follow folds and internal resonance tongues.
        % =================================================================
        function solve_omega(obj, mode_idx, om_start_factor, om_end_factor, varargin)
        % SOLVE_OMEGA  Compute NNM using omega as continuation parameter.
        %
        %   solve_omega(mode_idx, om_start_factor, om_end_factor)
        %   solve_omega(1, 1.0, 1.5)   -> from linear freq to 1.5x linear freq
        %
        %   Optional name-value pairs:
        %     'Ntd'         time steps per period (default 500)
        %     'ds'          initial arc-length step, fraction of om0 (default 0.005)
        %     'log10a_max'  stop when log10(a) exceeds this (default log10(10*a_gap))
        %     'log10a_min'  stop when log10(a) falls below this (default -8)

            % --- Parse optional arguments ---
            p = inputParser;
            addParameter(p, 'Ntd',       500,                      @isnumeric);
            addParameter(p, 'ds',        0.005,                    @isnumeric);
            addParameter(p, 'log10a_max', log10(obj.spring.a * 10), @isnumeric);
            addParameter(p, 'log10a_min', -8,                      @isnumeric);
            parse(p, varargin{:});
            Ntd        = p.Results.Ntd;
            ds_frac    = p.Results.ds;
            log10a_max = p.Results.log10a_max;
            log10a_min = p.Results.log10a_min;

            % --- Linear eigenvalue problem on ROM ---
            [inorm, nnorm, phi_start, om_lin] = ...
                obj.get_linear_mode(mode_idx);
            m = size(obj.P, 2);

            % --- Initial amplitude: well below clearance (linear regime) ---
            a0        = obj.spring.a * 1e-3;
            log10a0   = log10(a0);
            om0       = om_lin * om_start_factor;
            D0        = 0;

            ys0       = [phi_start; zeros(m, 1)];
            ys_nnorm0 = ys0(nnorm) / a0;

            % X_ext layout:   [ys_nnorm/a (2m-1);  D;  log10a;  om]
            %                  ^---free vars-----^              ^lam^
            % X_nlvib layout: [ys_nnorm/a (2m-1);  om;  D;  log10a]
            x0    = [ys_nnorm0; D0; log10a0];
            lam_s = om0;
            lam_e = om_lin * om_end_factor;

            % --- Scaling (fixed, no dynamic rescale to avoid loop) ---
            qscl = obj.spring.a;
            fscl = mean(diag(obj.ReducedAssembly.DATA.K)) * qscl;

            dscale = [ones(length(nnorm), 1); ...  % ys_nnorm/a  ~ 1
                      1e-2; ...                     % D           ~ 0
                      1.0; ...                      % log10a      ~ O(1)
                      om0];                         % omega (lam) ~ om0

            ds = ds_frac * om0;   % arc-length step in physical omega units

            % --- Termination criteria ---
            % Stop if amplitude goes too high or too low
            % In X_out, index of log10a = end-1, index of om = end
            term_high = @(X) X(end-1) >  log10a_max;
            term_low  = @(X) X(end-1) <  log10a_min;

            Sopt = struct(...
                'Dscale', dscale, ...
                'dynamicDscale', 1, ...
                'stepmax', 50000, ...               
                'dsmin', ds/10000, ...              
                'dsmax', 50*ds, ...                 
                'parametrization', 'pseudo_arc_length', ... 
                'pseudoxi', 0.5, ...   
                'predictor', 'secant', ... % <-- AGGIUNGI IL PREDITTORE SECANTE
                'noconv_stepcor', 'red', ...
                'reversaltolerance', 10.0);

            fprintf('--- NNM Shooting (omega param): mode %d ---\n', mode_idx);
            fprintf('    om_lin = %.2f rad/s = %.3f Hz\n', om_lin, om_lin/(2*pi));
            fprintf('    Range: %.2f -> %.2f rad/s\n', lam_s, lam_e);
            fprintf('    Stop: log10a in [%.1f, %.1f]\n', log10a_min, log10a_max);

            residuo = @(X_ext) obj.residual_om_as_lam( ...
                X_ext, Ntd, inorm, nnorm, qscl, fscl);

            obj.X_out = solve_and_continue(x0, residuo, lam_s, lam_e, ds, Sopt);

            fprintf('Done: %d solution points.\n', size(obj.X_out, 2));
            obj.compute_energies_omega(inorm, nnorm);
        end

        % -----------------------------------------------------------------
        function [R, dRdX_ext] = residual_om_as_lam(obj, X_ext, Ntd, ...
                                                      inorm, nnorm, qscl, fscl)
        % Wrapper that permutes X_ext -> X_nlvib before calling NLvib.
        %
        % X_ext   = [ys_nnorm/a (2m-1);  D;  log10a;  om]   <- om is lambda
        % X_nlvib = [ys_nnorm/a (2m-1);  om;  D;  log10a]   <- NLvib format
        %
        % Jacobian columns are permuted accordingly on the way back.

            n      = length(X_ext);
            D_val  = X_ext(n-2);
            log10a = X_ext(n-1);
            om     = X_ext(n);       % lambda

            % Build X in NLvib's expected format
            X_nlvib = [X_ext(1:n-3); om; D_val; log10a];

            [R, dRdX_nlvib] = shooting_residual(X_nlvib, obj.SystemNLvib, ...
                Ntd, 1, 'NMA', qscl, fscl, inorm);

            % Permute Jacobian columns: nlvib order -> ext order
            % nlvib: [1:n-3 | n-2=om | n-1=D | n=log10a]
            % ext:   [1:n-3 | n-2=D  | n-1=log10a | n=om]
            dRdX_ext = [dRdX_nlvib(:, 1:n-3), ...   % ys_nnorm/a
                        dRdX_nlvib(:, n-1), ...       % D
                        dRdX_nlvib(:, n), ...          % log10a
                        dRdX_nlvib(:, n-2)];           % om (was col n-2)
        end

        % =================================================================
        %  POST-PROCESSING: compute physical quantities from X_out
        % =================================================================
        function compute_energies_loga(obj, inorm, nnorm)
        % Post-process results from solve() where X_out format is:
        %   X_out = [ys_nnorm/a (2m-1);  om;  D;  log10a]
        %   col indices:               end-2  end-1  end

            m         = size(obj.P, 2);
            om_sh     = obj.X_out(end-2, :);
            log10a_sh = obj.X_out(end,   :);
            a_sh      = 10.^log10a_sh;

            obj.frequencies       = om_sh / (2*pi);
            obj.energies          = zeros(1, size(obj.X_out, 2));
            obj.modal_amplitudes  = zeros(m, size(obj.X_out, 2));
            obj.rms_displacements = zeros(1, size(obj.X_out, 2));

            obj.fill_energy_arrays(inorm, nnorm, a_sh, 1:end-3);
        end

        function compute_energies_omega(obj, inorm, nnorm)
        % Post-process results from solve_omega() where X_out format is:
        %   X_out = [ys_nnorm/a (2m-1);  D;  log10a;  om]
        %   col indices:               end-2  end-1  end

            m         = size(obj.P, 2);
            om_sh     = obj.X_out(end,   :);
            log10a_sh = obj.X_out(end-1, :);
            a_sh      = 10.^log10a_sh;

            obj.frequencies       = om_sh / (2*pi);
            obj.energies          = zeros(1, size(obj.X_out, 2));
            obj.modal_amplitudes  = zeros(m, size(obj.X_out, 2));
            obj.rms_displacements = zeros(1, size(obj.X_out, 2));

            obj.fill_energy_arrays(inorm, nnorm, a_sh, 1:end-3);
        end

        function fill_energy_arrays(obj, inorm, nnorm, a_sh, ys_rows)
        % Common loop for both post-processing paths.
        %   inorm    : scalar index of normalised DOF
        %   nnorm    : indices of free DOFs in ys vector
        %   a_sh     : [1 x Npts] amplitude values
        %   ys_rows  : row range of ys_nnorm/a in X_out (e.g. 1:end-3)

            M_r = obj.ReducedAssembly.DATA.M;
            K_r = obj.ReducedAssembly.DATA.K;
            m   = size(obj.P, 2);

            dir_mult = sign(obj.P(obj.spring.dof_idx, 1));
            if dir_mult == 0, dir_mult = 1; end

            for i = 1:size(obj.X_out, 2)
                ai = a_sh(i);

                % Reconstruct full state vector ys = [q; qdot]
                ys        = zeros(2*m, 1);
                ys(inorm) = ai;
                ys(nnorm) = obj.X_out(ys_rows, i) * ai;

                qi = ys(1:m);
                vi = ys(m+1:end);

                % Modal amplitudes (absolute value)
                obj.modal_amplitudes(:, i) = abs(qi);

                % Physical displacement vector
                x_full = obj.P * qi;

                % Spatial RMS over all physical DOFs
                obj.rms_displacements(i) = sqrt(mean(x_full.^2));

                % Total mechanical energy = kinetic + linear potential + NL spring
                E_kin  = 0.5 * vi' * M_r * vi;
                E_pot  = 0.5 * qi' * K_r * qi;
                u_nl   = x_full(obj.spring.dof_idx) * dir_mult;
                E_nl   = 0;
                if u_nl > obj.spring.a
                    E_nl = 0.5 * obj.spring.k * (u_nl - obj.spring.a)^2;
                end
                obj.energies(i) = E_kin + E_pot + E_nl;
            end
        end

        % =================================================================
        %  PLOTTING
        % =================================================================
        function plot_backbone(obj, style_str, display_name)
        % Plot the NNM backbone on a frequency-energy plot.
        % Call after solve() or solve_omega().
        %
        %   plot_backbone()
        %   plot_backbone('r--', 'ROM [1 3 5 MC]')
        %
        % Creates a semilogx plot (energy on log axis) matching Fig. 28.3
        % of Kuether et al. (2014).

            if nargin < 2 || isempty(style_str),   style_str    = 'b-';            end
            if nargin < 3 || isempty(display_name), display_name = 'NNM Shooting'; end

            semilogx(obj.energies, obj.frequencies, style_str, ...
                'LineWidth', 2, 'DisplayName', display_name);
            grid on;
            xlabel('Energy [J]');
            ylabel('Frequency [Hz]');
        end

        function plot_rms_backbone(obj, style_str, display_name)
        % Plot backbone with RMS displacement on x-axis instead of energy.
        % Useful for comparing with force-response curves.

            if nargin < 2 || isempty(style_str),   style_str    = 'b-';            end
            if nargin < 3 || isempty(display_name), display_name = 'NNM Shooting'; end

            semilogx(obj.rms_displacements, obj.frequencies, style_str, ...
                'LineWidth', 2, 'DisplayName', display_name);
            grid on;
            xlabel('RMS Displacement [m]');
            ylabel('Frequency [Hz]');
        end

        function plot_debug(obj)
        % Diagnostic plot to inspect continuation progress.
        % Shows frequency and log10(amplitude) vs. continuation step.
        % Useful to detect loops, reversals, or stagnation.

            figure('Name', 'NNM Continuation Debug', 'Color', 'w');

            subplot(2, 1, 1);
            plot(obj.frequencies, 'b-', 'LineWidth', 1.5);
            ylabel('Frequency [Hz]');
            xlabel('Continuation step');
            title('Frequency along continuation path');
            grid on;

            subplot(2, 1, 2);
            % Detect which format we have from size of X_out
            % Try to extract log10a heuristically
            if size(obj.X_out, 1) >= 2
                % Try both formats and pick the one with reasonable values
                loga_candidate1 = obj.X_out(end,   :);   % loga format
                loga_candidate2 = obj.X_out(end-1, :);   % omega format
                if abs(mean(loga_candidate1)) < abs(mean(loga_candidate2))
                    loga = loga_candidate1;
                else
                    loga = loga_candidate2;
                end
                plot(loga, 'r-', 'LineWidth', 1.5);
                hold on;
                yline(log10(obj.spring.a), 'k--', 'clearance', ...
                    'LabelHorizontalAlignment', 'left');
                ylabel('log_{10}(amplitude)');
                xlabel('Continuation step');
                title('Amplitude along continuation path');
                grid on;
            end
        end

    end % methods

    % =====================================================================
    methods (Access = private)
    % =====================================================================

        function [inorm, nnorm, phi_start, om_lin] = get_linear_mode(obj, mode_idx)
        % Solve the ROM eigenvalue problem and return mode data.
        %   inorm     : index of DOF used for normalisation (max abs component)
        %   nnorm     : all other DOF indices (2m-1 elements)
        %   phi_start : normalised mode shape (phi(inorm) = 1), pointing
        %               toward the contact spring (+a direction)
        %   om_lin    : natural frequency in rad/s

            M_r = obj.ReducedAssembly.DATA.M;
            K_r = obj.ReducedAssembly.DATA.K;
            m   = size(M_r, 1);

            [Phi_red, Om2_red] = eig(K_r, M_r);
            [om_vals, sort_idx] = sort(sqrt(abs(diag(Om2_red))));
            Phi_red = Phi_red(:, sort_idx);

            om_lin    = om_vals(mode_idx);
            phi_start = Phi_red(:, mode_idx);

            % Orient mode toward the contact spring (positive displacement)
            phi_phys = obj.P * phi_start;
            if phi_phys(obj.spring.dof_idx) < 0
                phi_start = -phi_start;
            end

            % Normalise so that the largest component equals 1
            [~, inorm] = max(abs(phi_start));
            phi_start  = phi_start / phi_start(inorm);

            nnorm = setdiff(1:2*m, [inorm, inorm+m]);
        end

    end % private methods

end % classdef