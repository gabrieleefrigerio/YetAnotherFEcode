classdef MasslessCBBuilder < handle
    properties
        Assembly, P, numModes, nl_dof, M_r, K_r, alpha
    end
    
    methods
        function obj = MasslessCBBuilder(yafec_assembly, num_fixed_modes, nonlinear_dof)
            obj.Assembly = yafec_assembly;
            obj.numModes = num_fixed_modes;
            obj.nl_dof = nonlinear_dof;
        end
        
        function build(obj)
            M_full = obj.Assembly.mass_matrix();
            K_full = obj.Assembly.stiffness_matrix();
            
            Mc = obj.Assembly.constrain_matrix(M_full);
            Kc = obj.Assembly.constrain_matrix(K_full);
            n_dofs_c = size(Kc, 1);
            
            unit_full = zeros(size(K_full, 1), 1);
            unit_full(obj.nl_dof) = 1;
            unit_c = obj.Assembly.constrain_vector(unit_full);
            nl_dof_c = find(unit_c);
            inner_idx_c = setdiff(1:n_dofs_c, nl_dof_c);
            
            K_ii = Kc(inner_idx_c, inner_idx_c);
            K_ib = Kc(inner_idx_c, nl_dof_c);
            M_ii = Mc(inner_idx_c, inner_idx_c);
            M_ib = Mc(inner_idx_c, nl_dof_c);
            
            % 1. Fixed-interface modes (Standard CB)
            m = obj.numModes;
            [Phi_i, D] = eigs(K_ii, M_ii, m, 'smallestabs');
            [~, sort_idx] = sort(diag(D));
            Phi_i = Phi_i(:, sort_idx);
            
            for i = 1:m
                Phi_i(:,i) = Phi_i(:,i) / sqrt(Phi_i(:,i)' * M_ii * Phi_i(:,i));
            end
            
            % 2. Static modes (Standard CB)
            Psi_c = - (K_ii \ full(K_ib));
            
            % 3. Calculation of alpha to nullify M_bi (Inertial Decoupling)
            % alpha = Phi^T * (M_ib + M_ii * Psi) 
            alpha = Phi_i' * (full(M_ib) + M_ii * Psi_c);
            
            % ... [rest of the existing code] ...
            
            % 4. New Projection Matrix (Component Modes)
            P_alpha = zeros(n_dofs_c, 1 + m);
            P_alpha(nl_dof_c, 1) = 1;
            P_alpha(inner_idx_c, 1) = Psi_c - Phi_i * alpha; 
            P_alpha(inner_idx_c, 2:end) = Phi_i;
            
            obj.P = obj.Assembly.unconstrain_vector(P_alpha);
            
            % 5. Final projection of reduced matrices
            K_complete = P_alpha' * Kc * P_alpha;
            M_complete = P_alpha' * Mc * P_alpha;
            
            % Assignment and elimination of boundary inertia (M_bb = 0)
            obj.K_r = K_complete;
            obj.M_r = M_complete;
            obj.M_r(1,1) = 0; % Replacement of the boundary mass block
            % Save alpha in the object
            obj.alpha = alpha;
            
        end
        
        function display_frequencies(obj)
            omega2 = diag(obj.K_r(2:end, 2:end));
            omega = sqrt(omega2);
            f_hz = omega / (2 * pi);
            
            fprintf('\n=============================================================\n');
            fprintf('   ROM NATURAL FREQUENCIES (Fixed-Interface Modes)   \n');
            fprintf('=============================================================\n');
            for idx = 1:length(f_hz)
                fprintf('   Mode %2d   |   Angular Frequency: %10.3f rad/s   |   Frequency: %10.3f Hz\n', ...
                    idx, omega(idx), f_hz(idx));
            end
            fprintf('=============================================================\n\n');
        end
    end
end