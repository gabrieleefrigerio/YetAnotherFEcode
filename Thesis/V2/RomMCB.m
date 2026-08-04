classdef RomMCB < handle
    properties
        Structure, P, Pc, numModes, contactDofs, n_bnd, M_r, K_r, alpha
    end
    
    methods
        function obj = RomMCB(dummy_struct, num_fixed_modes, contact_dofs_constrained)
            obj.Structure = dummy_struct;
            obj.numModes = num_fixed_modes;
            % FORZATURA A COLONNA
            obj.contactDofs = contact_dofs_constrained(:);
            obj.n_bnd = length(obj.contactDofs);
        end
        
        function build(obj)
            Mc = obj.Structure.AssemblyObj.constrain_matrix(obj.Structure.M);
            Kc = obj.Structure.AssemblyObj.constrain_matrix(obj.Structure.K);
            n_dofs_c = size(Kc, 1);
            
            fprintf('\n--- Building Massless CB ROM Base ---\n');
            
            % Partizionamento che gestisce n-nodi aggregati in automatico
            nl_dof_c = obj.contactDofs;
            inner_idx_c = setdiff(1:n_dofs_c, nl_dof_c);
            
            K_ii = Kc(inner_idx_c, inner_idx_c);
            K_ib = Kc(inner_idx_c, nl_dof_c);
            M_ii = Mc(inner_idx_c, inner_idx_c);
            M_ib = Mc(inner_idx_c, nl_dof_c);
            
            m = obj.numModes;
            [Phi_i, D] = eigs(K_ii, M_ii, m, 'smallestabs');
            [~, sort_idx] = sort(diag(D));
            Phi_i = Phi_i(:, sort_idx);
            
            for i = 1:m
                Phi_i(:,i) = Phi_i(:,i) / sqrt(Phi_i(:,i)' * M_ii * Phi_i(:,i));
            end
            
            Psi_c = - (K_ii \ full(K_ib));
            alpha = Phi_i' * (full(M_ib) + M_ii * Psi_c);
            
            P_alpha = zeros(n_dofs_c, obj.n_bnd + m);
            P_alpha(nl_dof_c, 1:obj.n_bnd) = eye(obj.n_bnd);
            P_alpha(inner_idx_c, 1:obj.n_bnd) = Psi_c - Phi_i * alpha; 
            P_alpha(inner_idx_c, obj.n_bnd+1:end) = Phi_i;
            
            obj.Pc = P_alpha;
            obj.P = obj.Structure.AssemblyObj.unconstrain_vector(P_alpha);
            
            obj.K_r = P_alpha' * Kc * P_alpha;
            M_complete = P_alpha' * Mc * P_alpha;
            
            obj.M_r = M_complete;
            obj.M_r(1:obj.n_bnd, 1:obj.n_bnd) = 0; 
            obj.alpha = alpha;
            fprintf('Base costruita: %d Physical DOFs + %d Modal DOFs (Massless Boundary)\n', obj.n_bnd, obj.numModes);
        end
        
        function [Mr, Kr, Cr] = get_reduced_matrices(obj)
            Mr = obj.M_r; Kr = obj.K_r;
            Cc = obj.Structure.AssemblyObj.constrain_matrix(obj.Structure.C);
            if isscalar(Cc) && Cc == 0, Cr = sparse(size(Mr,1), size(Mr,2)); else, Cr = obj.Pc' * Cc * obj.Pc; end
            Mr = (Mr + Mr') / 2; Kr = (Kr + Kr') / 2; Cr = (Cr + Cr') / 2;
        end
    end
end