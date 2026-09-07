classdef NNMHarmonicBalance < handle
    % NNMHARMONICBALANCE Calculation of Nonlinear Normal Modes using the 
    % Harmonic Balance (HB) method and Numerical Continuation.
    
    properties
        ReducedAssembly % YAFEC ReducedAssembly object
        spring          % LocalNonlinearity object (the contact spring)
        P               % ROM transformation matrix (Phi)
        SystemNLvib     % NLvib FE_system object
        
        % Outputs
        X_out           % Raw state vector from the continuation solver
        energies        % Total energy array [J]
        frequencies     % Modal frequencies array [Hz]
        segments        % Struct array filled by solve_with_tongues: one entry
                        % per continuation run, with fields X, energies,
                        % frequencies and i_exit. Points 1:i_exit lie on the
                        % backbone, the rest (if any) on the internal resonance
                        % tongue the branch turned into.
    end
    
    methods
        function obj = NNMHarmonicBalance(reduced_assembly, spring_obj, P_matrix)
            obj.ReducedAssembly = reduced_assembly;
            obj.spring = spring_obj;
            obj.P = P_matrix;
            
            m = size(P_matrix, 2);
            
            % --- THE WORKAROUND ---
            % We create a "dummy" custom function to satisfy the FE_system 
            % initialization, which strictly requires a 'fnl_CUSTOM' field
            % when using the 'custom' flag.
            obj.ReducedAssembly.DATA.fnl_CUSTOM = @(q) obj.dummy_force(q);
            
            % Now the FE_system constructor will safely compile without crashing
            obj.SystemNLvib = FE_system(obj.ReducedAssembly, zeros(m, 1), 'custom');
            
            % --- NLVIB NATIVE SPRING INJECTION ---
            % Check the arbitrary direction of the MATLAB eigenvector
            dir_mult = sign(obj.P(obj.spring.dof_idx, 1));
            if dir_mult == 0, dir_mult = 1; end
            
            % Compute the projection column vector for the specific DOF
            w = obj.P(obj.spring.dof_idx, :)' * dir_mult;

            k_nl   = obj.spring.k;
            gap_nl = obj.spring.a;
            if m == 1
                % HB_residual ignores force_direction when the system has a
                % single DOF: it hardcodes w = 1 (see the size(Q,1)==H+1 test
                % in NLvib). Fold the projection into an equivalent stiffness
                % and gap, which give exactly the same force:
                %   w*k*(w*q - gap) = (k*w^2)*(q - gap/w)
                k_nl   = obj.spring.k * w^2;
                gap_nl = obj.spring.a / w;
                w      = 1;
            end

            % Define the NLvib built-in unilateral spring structure
            nl_elem = struct('type', 'unilateralspring', ...
                             'force_direction', w, ...
                             'stiffness', k_nl, ...
                             'gap', gap_nl, ...
                             'islocal', true, ...        % Required by NLvib AFT
                             'ishysteretic', false);     % Required by NLvib AFT
                         
            % Overwrite the dummy custom function with the native analytical spring
            obj.SystemNLvib.nonlinear_elements = {nl_elem};
        end
        
        function solve(obj, mode_idx, log10a_start, log10a_end, H, N)
            % SOLVE Runs the Harmonic Balance continuation
            % H = Number of Harmonics
            % N = Number of time samples for Alternating Frequency/Time (AFT)
            if nargin < 5, H = 5; end
            if nargin < 6, N = 128; end

            [x0, inorm, om_start, fscl] = obj.initial_guess(mode_idx, H);

            % Harmonic Balance is numerically smoother, allowing for a generous step size
            ds =  0.0004;

            % Dscale helps the fsolve algorithm to balance variable magnitudes
            dscale = ones(size(x0,1)+1, 1);
            dscale(end-2) = om_start;
            Sopt = struct('Dscale', dscale, 'dynamicDscale', 1, 'dsmin', ds/100000, 'dsmax', 10*ds, 'stepmax', 10000, 'reversaltolerance', 10);

            fprintf('Starting HB (H=%d, N=%d) for NNM %d...\n', H, N, mode_idx);

            % Launch the continuation solver using the HB_residual function
            obj.X_out = solve_and_continue(x0, ...
                @(X) HB_residual(X, obj.SystemNLvib, H, N, 'NMA', inorm, fscl), ...
                log10a_start, log10a_end, ds, Sopt);

            fprintf('HB Calculation completed successfully!\n');

            % Post-process the results to compute frequencies and energies
            obj.compute_energies(H);
        end

        function solve_with_tongues(obj, mode_idx, log10a_start, log10a_end, H_bb, H_hi, N, opts)
            % SOLVE_WITH_TONGUES Backbone plus internal resonance tongues.
            %
            % An internal resonance of order n with the j-th mode of the ROM
            % shows up where the fundamental frequency of the branch crosses
            % f_j/n, and the solution needs the n-th harmonic to represent it.
            % That gives a clean separation:
            %
            %   phase 1 - backbone. Run the continuation with H_bb harmonics,
            %       chosen below the lowest resonance order the branch can meet
            %       (H_bb = 5 here, since f_2/f_1 is about 9). No tongue can be
            %       represented, so the branch has nothing to turn into and the
            %       backbone comes out in a single sweep.
            %
            %   phase 2 - tongues. For every crossing f_j/n, restart a
            %       continuation with H_hi harmonics slightly below the crossing
            %       amplitude, using the backbone solution padded with zeros as
            %       initial guess. The branch now follows the backbone up to the
            %       crossing and turns into the tongue, which is what we want to
            %       trace. It is not followed back to the backbone: along a
            %       tongue the primary modal coordinate vanishes while the
            %       resonating one takes over, and the amplitude normalisation
            %       of NLvib degenerates with it.
            %
            % opts fields (all optional):
            %   ds_bb, ds_t   arclength step, backbone and tongues (4e-4, 1e-4)
            %   delta         how far below the crossing a tongue run starts,
            %                 in log10 of the amplitude                 (0.05)
            %   stepmax_t     maximum steps of a tongue run             (1500)
            %   df_tol        frequency drop marking the exit [Hz]       (1.0)

            if nargin < 5 || isempty(H_bb), H_bb = 5;    end
            if nargin < 6 || isempty(H_hi), H_hi = 15;   end
            if nargin < 7 || isempty(N),    N = 2048;    end
            if nargin < 8, opts = struct(); end
            if ~isfield(opts, 'ds_bb'),     opts.ds_bb = 4e-4;   end
            if ~isfield(opts, 'ds_t'),      opts.ds_t = 1e-4;    end
            if ~isfield(opts, 'delta'),     opts.delta = 0.05;   end
            if ~isfield(opts, 'stepmax_t'), opts.stepmax_t = 1500; end
            if ~isfield(opts, 'df_tol'),    opts.df_tol = 1.0;   end

            m = size(obj.P, 2);
            obj.segments = struct('X', {}, 'energies', {}, 'frequencies', {}, ...
                                  'i_exit', {}, 'H', {}, 'kind', {});

            % ---------- phase 1: backbone ----------
            [x0, inorm, om_start, fscl] = obj.initial_guess(mode_idx, H_bb);
            dscale = ones(size(x0, 1) + 1, 1);
            dscale(end-2) = om_start;
            Sopt = struct('Dscale', dscale, 'dynamicDscale', 1, ...
                          'dsmin', opts.ds_bb/1e5, 'dsmax', 10*opts.ds_bb, ...
                          'stepmax', 10000, 'reversaltolerance', 0.05);

            Xbb = solve_and_continue(x0, ...
                @(X) HB_residual(X, obj.SystemNLvib, H_bb, N, 'NMA', inorm, fscl), ...
                log10a_start, log10a_end, opts.ds_bb, Sopt);
            [fbb, Ebb] = obj.fep_from_X(Xbb, H_bb);

            obj.segments(end+1) = struct('X', Xbb, 'energies', Ebb, ...
                'frequencies', fbb, 'i_exit', size(Xbb, 2), 'H', H_bb, ...
                'kind', 'backbone');
            fprintf('  backbone (H = %d): %d points, %.2f -> %.2f Hz\n', ...
                    H_bb, size(Xbb, 2), fbb(1), fbb(end));

            % ---------- phase 2: tongues ----------
            f_rom = sort(sqrt(eig(full(obj.ReducedAssembly.DATA.K), ...
                                  full(obj.ReducedAssembly.DATA.M)))) / (2*pi);
            lam_bb = Xbb(end, :);

            for j = 2:numel(f_rom)
                for n = (H_bb+1):H_hi
                    f_t = f_rom(j) / n;
                    if f_t <= fbb(1) || f_t >= max(fbb), continue; end

                    % Amplitude at which the backbone crosses the tongue
                    icross = find(fbb >= f_t, 1, 'first');
                    lam_t0 = lam_bb(icross) - opts.delta;
                    if lam_t0 <= log10a_start, continue; end

                    % Backbone solution just below the crossing, padded with
                    % zeros on the harmonics the coarse run did not carry
                    [~, istart] = min(abs(lam_bb - lam_t0));
                    Psi_hi = zeros((2*H_hi+1)*m, 1);
                    Psi_hi(1:(2*H_bb+1)*m) = Xbb(1:end-3, istart);
                    x_t = [Psi_hi; Xbb(end-2, istart); Xbb(end-1, istart)];

                    dscale_t = ones(size(x_t, 1) + 1, 1);
                    dscale_t(end-2) = om_start;
                    Sopt_t = struct('Dscale', dscale_t, 'dynamicDscale', 1, ...
                                    'dsmin', opts.ds_t/1e5, 'dsmax', 10*opts.ds_t, ...
                                    'stepmax', opts.stepmax_t, ...
                                    'reversaltolerance', 0.05);
                    try
                        Xt = solve_and_continue(x_t, ...
                            @(X) HB_residual(X, obj.SystemNLvib, H_hi, N, 'NMA', inorm, fscl), ...
                            lam_bb(istart), log10a_end, opts.ds_t, Sopt_t);
                    catch
                        fprintf('  tongue mode %d, 1:%-2d (%.2f Hz): no convergence\n', ...
                                j, n, f_t);
                        continue;
                    end

                    [ft, Et] = obj.fep_from_X(Xt, H_hi);
                    iex = obj.last_backbone_point(Xt, ft, opts.df_tol);
                    obj.segments(end+1) = struct('X', Xt, 'energies', Et, ...
                        'frequencies', ft, 'i_exit', iex, 'H', H_hi, ...
                        'kind', sprintf('mode %d, 1:%d', j, n));
                    fprintf('  tongue mode %d, 1:%-2d (%6.2f Hz): %4d points, left the backbone at %.2f Hz, ends at %.2f Hz\n', ...
                            j, n, f_t, size(Xt, 2), ft(iex), ft(end));
                end
            end

            % Default outputs describe the backbone
            obj.X_out       = Xbb;
            obj.frequencies = fbb;
            obj.energies    = Ebb;
        end

        function compute_energies(obj, H)
            % COMPUTE_ENERGIES Calculates the total mechanical energy for each step
            [obj.frequencies, obj.energies] = obj.fep_from_X(obj.X_out, H);
        end

        function [x0, inorm, om_start, fscl] = initial_guess(obj, mode_idx, H)
            % INITIAL_GUESS Linear mode of the ROM written as an HB state vector

            M_r = obj.ReducedAssembly.DATA.M;
            K_r = obj.ReducedAssembly.DATA.K;
            [Phi_red, Om2_red] = eig(K_r, M_r);
            [om_lin_red, sort_idx] = sort(sqrt(diag(Om2_red)));
            phi_start = Phi_red(:, sort_idx(mode_idx));

            % Ensure the initial eigenvector points *towards* the spring
            u_full = obj.P * phi_start;
            if u_full(obj.spring.dof_idx) < 0
                phi_start = -phi_start;
            end

            % Normalize the eigenvector
            [~, inorm] = max(abs(phi_start));
            phi_start = phi_start / phi_start(inorm);

            m = size(M_r, 1);

            % --- SETUP HB STATE VECTOR ---
            % In NLvib, Psi = [Q0; Q1_c; Q1_s; Q2_c; Q2_s; ...]
            % We place the linear mode shape entirely in the first harmonic (cosine)
            Psi = zeros((2*H+1)*m, 1);
            Psi(m + (1:m)) = phi_start;

            om_start = om_lin_red(mode_idx);

            % Initial guess vector x0 = [Psi; Om; modal_damping]
            x0 = [Psi; om_start; 0];

            % Use the clearance gap as the scaling baseline for displacements
            qscl = obj.spring.a;
            fscl = mean(diag(K_r)) * qscl;
        end

        function igood = last_backbone_point(obj, X, freq, df_tol)         %#ok<INUSL>
            % LAST_BACKBONE_POINT Index of the last point still on the backbone.
            %
            % Along the backbone of this system both the amplitude and the
            % frequency grow monotonically. A tongue always shows up as a
            % turning point of the amplitude parameter, and usually also as a
            % drop in frequency; small wobbles of either quantity are tolerated
            % so that a bumpy but still rising branch is not cut in two.

            n = size(X, 2);
            lam = X(end, :);
            lam_tol = 0.002;         % 0.2 % of a decade of amplitude
            igood = n;
            fmax = freq(1);
            lmax = lam(1);
            for i = 2:n
                if lam(i) < lmax - lam_tol || freq(i) < fmax - df_tol
                    igood = i - 1;
                    return;
                end
                fmax = max(fmax, freq(i));
                lmax = max(lmax, lam(i));
            end
        end

        function [freq, E] = fep_from_X(obj, X, H)
            % FEP_FROM_X Frequency [Hz] and total mechanical energy [J] of every
            % periodic solution stored in the continuation output X.

            m = size(obj.P, 2);
            Psi_HB = X(1:end-3, :);
            om_HB = X(end-2, :);
            a_HB = 10.^X(end, :);

            freq = om_HB / (2*pi);
            E = zeros(1, size(X, 2));

            M_r = obj.ReducedAssembly.DATA.M;
            K_r = obj.ReducedAssembly.DATA.K;

            dir_mult = sign(obj.P(obj.spring.dof_idx, 1));
            if dir_mult == 0, dir_mult = 1; end
            w = obj.P(obj.spring.dof_idx, :)' * dir_mult;

            for i = 1:size(X, 2)
                % Reconstruct the scaled harmonics matrix [m x (2H+1)]
                Qi = reshape(Psi_HB(:, i) * a_HB(i), m, 2*H+1);
                
                % Evaluate displacement and velocity at t=0.
                % x(t)   = Q0 + sum_k [ Qc_k*cos(k*Om*t) + Qs_k*sin(k*Om*t) ]
                % xdot(t)= sum_k k*Om*[ -Qc_k*sin(k*Om*t) + Qs_k*cos(k*Om*t) ]
                % so at t = 0 every sine harmonic must be weighted by its own
                % order k: omitting it underestimates the kinetic energy, which
                % matters a lot here because the contact generates strong
                % higher harmonics.
                q0 = Qi(:, 1);
                u0 = zeros(m, 1);
                for k = 1:H
                    q0 = q0 + Qi(:, 2*k);
                    u0 = u0 + k * om_HB(i) * Qi(:, 2*k+1);
                end
                
                % Linear kinetic and strain energy
                E_lin = 0.5 * u0' * M_r * u0 + 0.5 * q0' * K_r * q0;
                
                % Nonlinear strain energy (contact spring)
                u_nl = w' * q0;
                E_nl = 0;
                if u_nl > obj.spring.a
                    E_nl = 0.5 * obj.spring.k * (u_nl - obj.spring.a)^2;
                end
                
                E(i) = E_lin + E_nl;
            end
        end
        
        function plot_backbone(obj, style_str, display_name)
            % PLOT_BACKBONE Plots the Frequency-Energy Curve
            if nargin < 2, style_str = 'b-'; end
            if nargin < 3, display_name = 'NNM HB'; end
            
            semilogx(obj.energies, obj.frequencies, style_str, 'LineWidth', 2, 'DisplayName', display_name);
            grid on;
            xlabel('Energy [J]');
            ylabel('Frequency [Hz]');
        end
        
        function [Kt, f_nl] = dummy_force(~, q)
            % DUMMY_FORCE A placeholder function that returns zero.
            % It is only used to prevent the FE_system constructor from throwing 
            % a "missing field fnl_CUSTOM" error during initialization.
            m = length(q);
            Kt = sparse(m, m);
            f_nl = zeros(m, 1);
        end
    end
end