classdef RomRubin < handle
    properties
        Structure, P, Pc, numModes, contactDofs, n_bnd
    end
    
    methods
        function obj = RomRubin(dummy_struct, num_linear_modes, contact_dofs_constrained)
            obj.Structure = dummy_struct;
            obj.numModes = num_linear_modes;
            % FORZATURA A COLONNA
            obj.contactDofs = contact_dofs_constrained(:);
            obj.n_bnd = length(obj.contactDofs);
        end
        
        function build(obj)
            Mc = obj.Structure.AssemblyObj.constrain_matrix(obj.Structure.M);
            Kc = obj.Structure.AssemblyObj.constrain_matrix(obj.Structure.K);
            n_dofs_c = size(Kc, 1);
            
            fprintf('\n--- Building Classic Rubin ROM Base (CMS) ---\n');
            
            [Phi_k, D] = eigs(Kc, Mc, obj.numModes, 'smallestabs');
            [~, sort_idx] = sort(diag(D));
            Phi_k = Phi_k(:, sort_idx);
            D = D(sort_idx, sort_idx);
            omega2 = diag(D);
            
            for i = 1:size(Phi_k,2)
                Phi_k(:,i) = Phi_k(:,i) / sqrt(Phi_k(:,i)' * Mc * Phi_k(:,i));
            end
            
            % Ottimizzazione: Assegnazione sparsa diretta per tutti i nodi aggregati
            F_int = sparse(obj.contactDofs, 1:obj.n_bnd, 1, n_dofs_c, obj.n_bnd);
            
            Psi_a = Kc \ F_int;
            modal_participation = Phi_k' * F_int; 
            Lambda_inv = diag(1 ./ omega2);
            Psi_c = Phi_k * (Lambda_inv * modal_participation);
            Psi_res = Psi_a - Psi_c;
            
            T1 = [Psi_res, Phi_k];
            T1_b = T1(obj.contactDofs, :); 
            
            T1_bb = T1_b(:, 1:obj.n_bnd);
            T1_bi = T1_b(:, obj.n_bnd+1:end);
            
            inv_T1_bb = inv(T1_bb); 
            
            T2_top_left = inv_T1_bb;
            T2_top_right = -inv_T1_bb * T1_bi;
            T2_bot_left = zeros(obj.numModes, obj.n_bnd);
            T2_bot_right = eye(obj.numModes);
            
            T2 = [T2_top_left, T2_top_right; 
                  T2_bot_left, T2_bot_right];
              
            obj.Pc = T1 * T2;
            obj.P = obj.Structure.AssemblyObj.unconstrain_vector(obj.Pc);
            fprintf('--- Rubin ROM Build Complete ---\n');
        end
        
        function [Mr, Kr, Cr] = get_reduced_matrices(obj)
            Mc = obj.Structure.AssemblyObj.constrain_matrix(obj.Structure.M);
            Kc = obj.Structure.AssemblyObj.constrain_matrix(obj.Structure.K);
            Cc = obj.Structure.AssemblyObj.constrain_matrix(obj.Structure.C);
            Mr = obj.Pc' * Mc * obj.Pc; Kr = obj.Pc' * Kc * obj.Pc;
            if isscalar(Cc) && Cc == 0, Cr = sparse(size(Mr,1), size(Mr,2)); else, Cr = obj.Pc' * Cc * obj.Pc; end
            Mr = (Mr + Mr') / 2; Kr = (Kr + Kr') / 2; Cr = (Cr + Cr') / 2;
        end
    end
end