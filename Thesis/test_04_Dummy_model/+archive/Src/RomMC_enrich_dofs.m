classdef RomMC_enrich_dofs < handle
    % ROMMC Builds the reduction basis (ROM) using Milman-Chu vectors 
    % Multi-vector approach with SELECTIVE ENRICHMENT DOFs.
    
    properties
        DummyStruct % Reference to the DummyStructureOG object
        P           % Transformation matrix P = [Phi, Psi] expanded to full space
        numModes    % Number of linear modes to include
        contactDofs % Constrained indices of the DOFs where the penalty acts
        spring_k    % Penalty stiffness of the contact wall [N/m]
        include_MC  % Boolean flag to include MC vector in P
        enrichDofs  % SOTTOINSIEME di DoF da perturbare per la base statica
    end
    
    methods
        function obj = RomMC_enrich_dofs(dummy_struct, num_linear_modes, contact_dofs_constrained, penalty_k, include_MC, enrich_dofs)
            % Constructor
            obj.DummyStruct = dummy_struct;
            obj.numModes = num_linear_modes;
            obj.contactDofs = contact_dofs_constrained;
            
            if nargin < 4 || isempty(penalty_k)
                obj.spring_k = 0;
            else
                obj.spring_k = penalty_k;
            end
            
            if nargin < 5 || isempty(include_MC)
                obj.include_MC = true;
            else
                obj.include_MC = include_MC;
            end
            
            % GESTIONE SOTTOCAMPIONAMENTO NODI
            if nargin < 6 || isempty(enrich_dofs)
                obj.enrichDofs = contact_dofs_constrained; % Default: tutti i nodi
            else
                obj.enrichDofs = enrich_dofs;
            end
        end
        
        function build(obj)
            % 1. Retrieve CONSTRAINED matrices
            Mc = obj.DummyStruct.AssemblyObj.constrain_matrix(obj.DummyStruct.M);
            Kc = obj.DummyStruct.AssemblyObj.constrain_matrix(obj.DummyStruct.K);
            
            n_dofs_c = size(Kc, 1);
            n_bnd = length(obj.contactDofs);
            n_enrich = length(obj.enrichDofs);
            
            fprintf('\n--- Building Selective Milman-Chu ROM Base ---\n');
            fprintf('Extracting %d linear modes...\n', obj.numModes);
            
            % 2. Calculation of CONSTRAINED linear modes (Phi_c)
            [Phi_c, D] = eigs(Kc, Mc, obj.numModes, 'smallestabs');
            [~, sort_idx] = sort(diag(D));
            Phi_c = Phi_c(:, sort_idx);
            
            for i = 1:size(Phi_c,2)
                Phi_c(:,i) = Phi_c(:,i) / sqrt(Phi_c(:,i)' * Mc * Phi_c(:,i));
            end
            
            if obj.include_MC
                fprintf('Calculating %d Milman-Chu static enrichment vectors (Out of %d Interface DoFs)...\n', n_enrich, n_bnd);
                
                % 3. Create Interface Force Matrix SOLO per i nodi scelti
                F_int_c = sparse(obj.enrichDofs, 1:n_enrich, 1, n_dofs_c, n_enrich);
                
                % 4. Create a constrained stiffness matrix for the closed gap
                % LA PENALITA' AGISCE SEMPRE SU TUTTI I NODI DI CONTATTO
                K_spring_c = sparse(obj.contactDofs, obj.contactDofs, obj.spring_k, n_dofs_c, n_dofs_c);
                Kc_closed = Kc + K_spring_c;
                
                % 5. Calculate the static response (Psi_c)
                Psi_c = Kc_closed \ F_int_c;
                
                V_initial_c = [Phi_c, Psi_c];
                fprintf('Base constructed: %d Linear Modes + %d MC Vectors\n', obj.numModes, n_enrich);
            else
                V_initial_c = Phi_c;
                fprintf('Base constructed: %d Linear Modes (MC Skipped)\n', obj.numModes);
            end
            
            % 6. Orthonormalization
            fprintf('Applying Mass-Orthogonal Gram-Schmidt...\n');
            P_c = obj.gram_schmidt(V_initial_c, Mc);
            
            % 7. EXPANSION TO THE FULL SPACE
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
                fprintf('Warning: Removed %d linearly dependent MC modes.\n', sum(zero_cols));
                P_ortho(:, zero_cols) = [];
            end
        end
        
        function display_rom_frequencies(obj)
            % Omissis (Stessa di prima)
        end
        
        function [Mr, Kr, Cr] = get_reduced_matrices(obj)
            if isempty(obj.P)
                error('The projection matrix P is empty. Call build() first.');
            end
            Mr = obj.P' * obj.DummyStruct.M * obj.P;
            Kr = obj.P' * obj.DummyStruct.K * obj.P;
            Cr = obj.P' * obj.DummyStruct.C * obj.P;
            Mr = (Mr + Mr') / 2; Kr = (Kr + Kr') / 2; Cr = (Cr + Cr') / 2;
        end
    end
end