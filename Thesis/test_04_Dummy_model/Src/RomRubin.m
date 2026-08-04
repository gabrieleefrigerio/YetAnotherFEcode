classdef RomRubin < handle
    % ROMRUBIN Classic Rubin free-interface component mode synthesis.
    %
    % The basis collects the residual attachment modes and the retained free
    % interface normal modes, then transforms them so that the first n_bnd
    % coordinates are the physical interface displacements. Unlike MacNeal,
    % the reduced matrices are consistent (Galerkin): Kr = Pc'*Kc*Pc.

    properties
        Structure, P, Pc, numModes, contactDofs, n_bnd
        scaleD                % Per-column scaling factors
        applyScaling = true   % Flag, kept for A/B comparison
    end

    methods
        function obj = RomRubin(dummy_struct, num_linear_modes, contact_dofs_constrained, apply_scaling)
            obj.Structure   = dummy_struct;
            obj.numModes    = num_linear_modes;
            obj.contactDofs = contact_dofs_constrained(:);
            obj.n_bnd       = length(obj.contactDofs);
            if nargin >= 4 && ~isempty(apply_scaling)
                obj.applyScaling = logical(apply_scaling);
            end
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

            % Direct sparse assignment for all aggregated interface DOFs
            F_int = sparse(obj.contactDofs, 1:obj.n_bnd, 1, n_dofs_c, obj.n_bnd);

            % Residual attachment modes: full flexibility minus the part
            % already spanned by the retained modes
            Psi_a = Kc \ F_int;
            modal_participation = Phi_k' * F_int;
            Lambda_inv = diag(1 ./ omega2);
            Psi_c = Phi_k * (Lambda_inv * modal_participation);
            Psi_res = Psi_a - Psi_c;

            T1 = [Psi_res, Phi_k];
            T1_b = T1(obj.contactDofs, :);

            T1_bb = T1_b(:, 1:obj.n_bnd);
            T1_bi = T1_b(:, obj.n_bnd+1:end);

            % Change of coordinates that makes the first n_bnd generalized
            % coordinates coincide with the physical interface displacements
            inv_T1_bb = inv(T1_bb);

            T2_top_left = inv_T1_bb;
            T2_top_right = -inv_T1_bb * T1_bi;
            T2_bot_left = zeros(obj.numModes, obj.n_bnd);
            T2_bot_right = eye(obj.numModes);

            T2 = [T2_top_left, T2_top_right;
                  T2_bot_left, T2_bot_right];

            obj.Pc = T1 * T2;
            if obj.applyScaling
                % Diagonal scaling in the mass norm. It preserves the span
                % (invertible change of basis) and makes interface
                % coordinates [m] and modal amplitudes [kg^-1/2 m]
                % commensurable. The mass metric is the right one because the
                % ode15s iteration matrix, W = M/(h*gamma) - J, is dominated
                % by the mass term as h -> 0.
                d = sqrt(full(sum(obj.Pc .* (Mc * obj.Pc), 1)));   % 1 x r
                obj.scaleD = d(:);
                obj.Pc = obj.Pc ./ d;
            else
                obj.scaleD = ones(size(obj.Pc,2), 1);
            end

            Mr_dbg = full(obj.Pc' * Mc * obj.Pc);
            fprintf('  spread d = %.2e | cond(Mr) = %.3e | cond(Kr) = %.3e\n', ...
                max(obj.scaleD)/min(obj.scaleD), cond(Mr_dbg), ...
                cond(full(obj.Pc' * Kc * obj.Pc)));

            obj.P = obj.Structure.AssemblyObj.unconstrain_vector(obj.Pc);
            fprintf('--- Rubin ROM Build Complete ---\n');
        end

        function [Mr, Kr, Cr] = get_reduced_matrices(obj)
            Mc = obj.Structure.AssemblyObj.constrain_matrix(obj.Structure.M);
            Kc = obj.Structure.AssemblyObj.constrain_matrix(obj.Structure.K);
            Cc = obj.Structure.AssemblyObj.constrain_matrix(obj.Structure.C);
            Mr = obj.Pc' * Mc * obj.Pc;
            Kr = obj.Pc' * Kc * obj.Pc;
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

        function [gap_run, k_run] = contact_params(obj, gap_phys, k_phys)
            % CONTACT_PARAMS Gap and penalty stiffness in reduced coordinates.
            %
            % The scaling above maps physical interface displacements to
            % generalized ones as q_phys_i = y_i / d_i, hence
            %   gap'_i = d_i * gap_i ,   k'_i = k / d_i^2
            % With applyScaling = false all d_i are 1 and the values are
            % returned unchanged.
            d = obj.scaleD(1:obj.n_bnd);
            gap_run = d .* gap_phys(:);
            k_run   = k_phys ./ d.^2;
        end
    end
end
