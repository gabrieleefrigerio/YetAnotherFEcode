classdef RomMN < handle
    properties
        Structure, P, Pc, numModes, contactDofs, n_bnd, M_r, K_r
    end
    
    methods
        function obj = RomMN(dummy_struct, num_linear_modes, contact_dofs_constrained)
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
            
            % --- Difesa numerica: piccole asimmetrie introdotte dall'eliminazione
            % dei vincoli possono far restituire a eigs autovalori/autovettori
            % con parte immaginaria spuria, che si propaga corrompendo l'intero
            % ROM. Si forza la simmetria prima di procedere.
            Kc = (Kc + Kc') / 2;
            Mc = (Mc + Mc') / 2;
            
            fprintf('\n--- Building MacNeal ROM Base (Free-Interface CMS) ---\n');
            
            nl_dof = obj.contactDofs;
            inner_idx = setdiff(1:n_dofs_c, nl_dof);
            
            m = obj.numModes;
            [Phi, D] = eigs(Kc, Mc, m, 'smallestabs');
            
            % --- Cast difensivo a reale: per matrici simmetriche/SPD gli
            % autovalori/autovettori devono essere reali; eventuali parti
            % immaginarie sono solo rumore numerico.
            if ~isreal(Phi) || ~isreal(D)
                max_imag = max(max(abs(imag(Phi(:)))), max(abs(imag(diag(D)))));
                if max_imag > 1e-8
                    warning('RomMN:ComplexEigs', ...
                        ['eigs ha restituito una parte immaginaria non trascurabile ' ...
                         '(max %.3e). Verificare la simmetria di M, K.'], max_imag);
                end
                Phi = real(Phi);
                D = real(D);
            end
            
            [~, sort_idx] = sort(diag(D));
            Phi = Phi(:, sort_idx);
            omega2 = diag(D(sort_idx, sort_idx));
            

            
            for i = 1:m
                Phi(:,i) = Phi(:,i) / sqrt(Phi(:,i)' * Mc * Phi(:,i));
            end
            
            Phi_b = Phi(nl_dof, :);
            Phi_i = Phi(inner_idx, :);
            
            F_int = sparse(nl_dof, 1:obj.n_bnd, 1, n_dofs_c, obj.n_bnd);
            Flex_cols = Kc \ F_int; 
            
            F_bb = Flex_cols(nl_dof, :);
            F_ib = Flex_cols(inner_idx, :);
            
            inv_omega2 = diag(1 ./ omega2);
            F_bb_res = F_bb - Phi_b * inv_omega2 * Phi_b';
            F_ib_res = F_ib - Phi_i * inv_omega2 * Phi_b';
            
            % --- Controllo condizionamento: F_bb_res quasi singolare indica
            % che i modi trattenuti non catturano abbastanza la flessibilita'
            % statica al boundary (numModes troppo basso rispetto a n_bnd).
            if rcond(F_bb_res) < 1e-12
                warning('RomMN:IllConditionedFbb', ...
                    ['F_bb_res mal condizionata (rcond = %.2e). Aumentare numModes ' ...
                     'o verificare la scelta dei contactDofs.'], rcond(F_bb_res));
            end
            
            T_ib = F_ib_res / F_bb_res;
            
            Pc_matrix = zeros(n_dofs_c, obj.n_bnd + m);
            Pc_matrix(nl_dof, 1:obj.n_bnd) = eye(obj.n_bnd);
            Pc_matrix(inner_idx, 1:obj.n_bnd) = T_ib;
            Pc_matrix(nl_dof, obj.n_bnd+1:end) = zeros(obj.n_bnd, m);
            Pc_matrix(inner_idx, obj.n_bnd+1:end) = Phi_i - T_ib * Phi_b;
            
            obj.Pc = Pc_matrix;
            obj.P = obj.Structure.AssemblyObj.unconstrain_vector(Pc_matrix);
            
            obj.K_r = Pc_matrix' * Kc * Pc_matrix;
            obj.M_r = zeros(obj.n_bnd + m);
            obj.M_r(obj.n_bnd+1:end, obj.n_bnd+1:end) = eye(m);
            
            fprintf('Base costruita: %d Physical DOFs + %d Modal DOFs (MacNeal Massless)\n', obj.n_bnd, obj.numModes);
        end
        
        function [Mr, Kr, Cr] = get_reduced_matrices(obj)
            Mr = obj.M_r; Kr = obj.K_r;
            Cc = obj.Structure.AssemblyObj.constrain_matrix(obj.Structure.C);
            if isscalar(Cc) && Cc == 0, Cr = sparse(size(Mr,1), size(Mr,2)); else, Cr = obj.Pc' * Cc * obj.Pc; end
            Mr = (Mr + Mr') / 2; Kr = (Kr + Kr') / 2; Cr = (Cr + Cr') / 2;
        end
    end
end