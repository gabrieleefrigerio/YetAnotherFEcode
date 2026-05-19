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
            
            % Define the NLvib built-in unilateral spring structure
            nl_elem = struct('type', 'unilateralspring', ...
                             'force_direction', w, ...
                             'stiffness', obj.spring.k, ...
                             'gap', obj.spring.a, ...
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
            
            % Extract linear matrices to compute the initial guess
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
        
        function compute_energies(obj, H)
            % COMPUTE_ENERGIES Calculates the total mechanical energy for each step
            m = size(obj.P, 2);
            Psi_HB = obj.X_out(1:end-3, :);
            om_HB = obj.X_out(end-2, :);
            log10a_HB = obj.X_out(end, :);
            a_HB = 10.^log10a_HB;
            
            obj.frequencies = om_HB / (2*pi);
            obj.energies = zeros(1, size(obj.X_out, 2));
            
            M_r = obj.ReducedAssembly.DATA.M;
            K_r = obj.ReducedAssembly.DATA.K;
            
            dir_mult = sign(obj.P(obj.spring.dof_idx, 1));
            if dir_mult == 0, dir_mult = 1; end
            w = obj.P(obj.spring.dof_idx, :)' * dir_mult;
            
            for i = 1:size(obj.X_out, 2)
                % Reconstruct the scaled harmonics matrix [m x (2H+1)]
                Qi = reshape(Psi_HB(:, i) * a_HB(i), m, 2*H+1);
                
                % Evaluate displacement and velocity at t=0 (assuming peak oscillation)
                q0 = Qi(:, 1) + sum(Qi(:, 2:2:end), 2); 
                u0 = sum(Qi(:, 3:2:end), 2) .* om_HB(i); 
                
                % Linear kinetic and strain energy
                E_lin = 0.5 * u0' * M_r * u0 + 0.5 * q0' * K_r * q0;
                
                % Nonlinear strain energy (contact spring)
                u_nl = w' * q0;
                E_nl = 0;
                if u_nl > obj.spring.a
                    E_nl = 0.5 * obj.spring.k * (u_nl - obj.spring.a)^2;
                end
                
                obj.energies(i) = E_lin + E_nl;
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