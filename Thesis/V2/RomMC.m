classdef RomMC < handle
    % ROMMC Builds the reduction basis (ROM) using Milman-Chu vectors 
    % Multi-vector approach: one enrichment vector per contact DOF.
    
    properties
        DummyStruct 
        P           
        numModes    
        contactDofs 
        spring_k    
        include_MC  
    end
    
    methods
        function obj = RomMC(dummy_struct, num_linear_modes, contact_dofs_constrained, penalty_k, include_MC)
            obj.DummyStruct = dummy_struct;
            obj.numModes = num_linear_modes;
            % FORZATURA A VETTORE COLONNA PER MULTI-INTERFACCIA
            obj.contactDofs = contact_dofs_constrained(:); 
            
            if nargin < 4 || isempty(penalty_k)
                obj.spring_k = 0;
            else
                obj.spring_k = penalty_k;
            end
            
            if nargin < 5
                obj.include_MC = true;
            else
                obj.include_MC = include_MC;
            end
        end
        
        function build(obj)
            Mc = obj.DummyStruct.AssemblyObj.constrain_matrix(obj.DummyStruct.M);
            Kc = obj.DummyStruct.AssemblyObj.constrain_matrix(obj.DummyStruct.K);
            
            n_dofs_c = size(Kc, 1);
            n_bnd = length(obj.contactDofs);
            
            fprintf('\n--- Building Multi-Vector Milman-Chu ROM Base ---\n');
            fprintf('Extracting %d linear modes...\n', obj.numModes);
            
            [Phi_c, D] = eigs(Kc, Mc, obj.numModes, 'smallestabs');
            [~, sort_idx] = sort(diag(D));
            Phi_c = Phi_c(:, sort_idx);
            
            for i = 1:size(Phi_c,2)
                Phi_c(:,i) = Phi_c(:,i) / sqrt(Phi_c(:,i)' * Mc * Phi_c(:,i));
            end
            
            if obj.include_MC
                fprintf('Calculating %d Milman-Chu static enrichment vectors for %d interfaces...\n', n_bnd, n_bnd);
                
                % Forza unitaria su ogni nodo dell'interfaccia aggregata
                F_int_c = sparse(obj.contactDofs, 1:n_bnd, 1, n_dofs_c, n_bnd);
                
                % Rigidezza applicata a tutti i nodi di contatto contemporaneamente
                K_spring_c = sparse(obj.contactDofs, obj.contactDofs, obj.spring_k, n_dofs_c, n_dofs_c);
                Kc_closed = Kc + K_spring_c;
                
                Psi_c = Kc_closed \ F_int_c;
                V_initial_c = [Phi_c, Psi_c];
                fprintf('Base constructed: %d Linear Modes + %d MC Vectors\n', obj.numModes, n_bnd);
            else
                V_initial_c = Phi_c;
                fprintf('Base constructed: %d Linear Modes (MC Skipped)\n', obj.numModes);
            end
            
            fprintf('Applying Mass-Orthogonal Gram-Schmidt...\n');
            P_c = obj.gram_schmidt(V_initial_c, Mc);
            
            obj.P = obj.DummyStruct.AssemblyObj.unconstrain_vector(P_c);
            fprintf('--- ROM Build Complete ---\n');
        end
        
        function P_ortho = gram_schmidt(obj, V, M)
            n_vecs = size(V, 2);
            P_ortho = zeros(size(V));
            for i = 1:n_vecs
                v_curr = V(:,i);
                for j = 1:i-1
                    proj_coeff = (P_ortho(:,j)' * M * v_curr) / (P_ortho(:,j)' * M * P_ortho(:,j));
                    v_curr = v_curr - proj_coeff * P_ortho(:,j);
                end
                norm_v = sqrt(v_curr' * M * v_curr);
                if norm_v > 1e-20
                    P_ortho(:,i) = v_curr / norm_v;
                else
                    P_ortho(:,i) = 0; 
                end
            end
            zero_cols = all(P_ortho == 0, 1);
            if any(zero_cols)
                fprintf('Warning: Removed %d linearly dependent MC modes during orthogonalization.\n', sum(zero_cols));
                P_ortho(:, zero_cols) = [];
            end
        end
        
        function [Mr, Kr, Cr] = get_reduced_matrices(obj)
            Mr = obj.P' * obj.DummyStruct.M * obj.P;
            Kr = obj.P' * obj.DummyStruct.K * obj.P;
            Cr = obj.P' * obj.DummyStruct.C * obj.P;
            Mr = (Mr + Mr') / 2;
            Kr = (Kr + Kr') / 2;
            Cr = (Cr + Cr') / 2;
        end
    end
end