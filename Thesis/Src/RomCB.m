classdef RomCB < handle
    % ROMCB Classic Craig-Bampton reduction basis.
    %
    % The interface degrees of freedom are kept as physical coordinates.
    % Unlike the massless Craig-Bampton (see RomMCB), the interface mass is
    % retained: there is no inertial decoupling through the alpha transform,
    % and M_r(bnd,bnd) is left untouched.
    %
    % This is the Hurty/Craig-Bampton model the interface reduction literature
    % is built on, so it is the natural target of interface_reduction(): its
    % M_r(bnd,bnd) is SPD, which the massless variants' is not.

    properties
        Structure, P, Pc, numModes, contactDofs, n_bnd, M_r, K_r
    end

    methods
        function obj = RomCB(dummy_struct, num_fixed_modes, contact_dofs_constrained)
            obj.Structure = dummy_struct;
            obj.numModes = num_fixed_modes;
            obj.contactDofs = contact_dofs_constrained(:);
            obj.n_bnd = length(contact_dofs_constrained);
        end

        function build(obj)
            Mc = obj.Structure.AssemblyObj.constrain_matrix(obj.Structure.M);
            Kc = obj.Structure.AssemblyObj.constrain_matrix(obj.Structure.K);
            n_dofs_c = size(Kc, 1);

            fprintf('\n--- Building Classic Craig-Bampton (CB) ROM Base ---\n');

            % Partitioning
            nl_dof_c = obj.contactDofs;
            inner_idx_c = setdiff(1:n_dofs_c, nl_dof_c);

            K_ii = Kc(inner_idx_c, inner_idx_c);
            K_ib = Kc(inner_idx_c, nl_dof_c);
            M_ii = Mc(inner_idx_c, inner_idx_c);
            % M_ib and M_bb are not needed explicitly to build the classic CB basis

            % 1. Fixed-interface modes
            m = obj.numModes;
            [Phi_i, D] = eigs(K_ii, M_ii, m, 'smallestabs');
            [~, sort_idx] = sort(diag(D));
            Phi_i = Phi_i(:, sort_idx);

            % Mass normalization
            for i = 1:m
                Phi_i(:,i) = Phi_i(:,i) / sqrt(Phi_i(:,i)' * M_ii * Phi_i(:,i));
            end

            % 2. Static constraint modes
            Psi_c = - (K_ii \ full(K_ib));

            % 3. Projection matrix (standard CB).
            % No alpha term: the basis is simply [I, 0; Psi_c, Phi_i].
            P_cb = zeros(n_dofs_c, obj.n_bnd + m);
            P_cb(nl_dof_c, 1:obj.n_bnd) = eye(obj.n_bnd);
            P_cb(inner_idx_c, 1:obj.n_bnd) = Psi_c;
            P_cb(inner_idx_c, obj.n_bnd+1:end) = Phi_i;

            obj.Pc = P_cb;
            obj.P = obj.Structure.AssemblyObj.unconstrain_vector(P_cb);

            % 4. Reduced matrices.
            % Note that M_r(1:n_bnd, 1:n_bnd) is deliberately NOT zeroed:
            % that is what separates this class from RomMCB.
            obj.K_r = P_cb' * Kc * P_cb;
            obj.M_r = P_cb' * Mc * P_cb;

            fprintf('Basis built: %d physical DOFs + %d modal DOFs (classic CB, mass retained)\n', ...
                obj.n_bnd, obj.numModes);
        end

        function [Mr, Kr, Cr] = get_reduced_matrices(obj)
            if isempty(obj.Pc)
                error('RomCB:NotBuilt', 'Projection matrix is empty. Call build() first.');
            end
            Mr = obj.M_r;
            Kr = obj.K_r;

            Cc = obj.Structure.AssemblyObj.constrain_matrix(obj.Structure.C);
            if isscalar(Cc) && Cc == 0
                Cr = sparse(size(Mr,1), size(Mr,2));
            else
                Cr = obj.Pc' * Cc * obj.Pc;
            end

            % Defensive symmetrization
            Mr = (Mr + Mr') / 2;
            Kr = (Kr + Kr') / 2;
            Cr = (Cr + Cr') / 2;
        end
    end
end
