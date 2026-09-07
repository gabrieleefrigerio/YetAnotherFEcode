classdef RomMN < handle
    % ROMMN MacNeal free-interface CMS with a massless boundary.
    %
    % Reference: Monjaraz Tec et al., Eq. (5.3)-(5.5).
    %
    % IMPORTANT: MacNeal is an INCONSISTENT (non-Galerkin) method. K_r is NOT
    % Pc'*Kc*Pc, which would be Rubin. It has to be assembled explicitly from
    % the residual flexibility.

    properties
        Structure
        P, Pc
        numModes
        contactDofs
        n_bnd
        M_r, K_r, C_r
        omega2          % Eigenvalues of the retained free-interface modes
        Phi_b           % Modes restricted to the boundary DOFs
        Fbb_res         % Residual flexibility at the boundary
        zeta            % Modal damping ratios actually used
    end

    methods
        function obj = RomMN(dummy_struct, num_linear_modes, contact_dofs_constrained)
            obj.Structure   = dummy_struct;
            obj.numModes    = num_linear_modes;
            obj.contactDofs = contact_dofs_constrained(:);
            obj.n_bnd       = length(obj.contactDofs);
        end

        function build(obj, zeta_modal)
            % BUILD Assemble the MacNeal reduced model.
            %   zeta_modal : scalar | vector [m x 1] | struct('alpha',a,'beta',b)
            %                Defaults to 0 if omitted.
            if nargin < 2, zeta_modal = 0; end

            Mc = obj.Structure.AssemblyObj.constrain_matrix(obj.Structure.M);
            Kc = obj.Structure.AssemblyObj.constrain_matrix(obj.Structure.K);
            n_dofs_c = size(Kc, 1);

            % Defensive symmetrization
            Kc = (Kc + Kc') / 2;
            Mc = (Mc + Mc') / 2;

            fprintf('\n--- Building MacNeal ROM (free-interface, massless boundary) ---\n');

            nl_dof    = obj.contactDofs;
            inner_idx = setdiff(1:n_dofs_c, nl_dof)';
            m         = obj.numModes;

            % ---------- 1. Free-interface normal modes ----------
            [Phi, D] = eigs(Kc, Mc, m, 'smallestabs');

            if ~isreal(Phi) || ~isreal(D)
                max_imag = max(max(abs(imag(Phi(:)))), max(abs(imag(diag(D)))));
                if max_imag > 1e-8
                    warning('RomMN:ComplexEigs', ...
                        'eigs returned an imaginary part of %.3e. Check the symmetry of M and K.', ...
                        max_imag);
                end
                Phi = real(Phi);  D = real(D);
            end

            [w2, sort_idx] = sort(diag(D));
            Phi = Phi(:, sort_idx);

            % ---------- 2. Rigid body mode guard ----------
            % MacNeal needs a non-singular K (it uses F = inv(K)), so a rigid
            % body motion is fatal. A rigid body mode is one in the NULL SPACE
            % OF K, and that is what gets tested here.
            %
            % The previous test compared the lowest retained eigenvalue with the
            % highest, flagging w2 < 1e-8*max(w2). That is a comparison between
            % two ELASTIC modes whenever there is no rigid motion at all, so it
            % turns into a false positive as soon as the basis spans more than
            % four decades of frequency - which a few hundred modes routinely
            % do. It rejected a perfectly legitimate 4982 Hz mode on the 2D
            % model for no reason other than the basis being wide.
            %
            % ||K*phi|| / (||K||*||phi||) is dimensionless and scale free, so it
            % behaves the same in the SI and the um/MPa/kg systems: a genuine
            % rigid mode sits at round-off, an elastic one does not.
            % Threshold placed between the two populations, not at the edge of
            % one of them: measured on the 2D model with 300 modes, the LOWEST
            % legitimate elastic mode sits at 1.3e-10 while a true null-space
            % vector sits at round-off, about 1e-16. 1e-13 is three orders above
            % the second and three below the first.
            Kn   = normest(Kc, 1e-3);
            rres = vecnorm(Kc*Phi) ./ (Kn * vecnorm(Phi) + realmin);
            rb_idx = find(rres < 1e-13);
            if ~isempty(rb_idx)
                error('RomMN:RigidBodyModes', ...
                    ['Found %d rigid body mode(s): K*phi is zero to round-off ' ...
                     '(relative residual %.3e) at f = %.4g Hz.\n' ...
                     'MacNeal cannot be applied when rigid body motions are ' ...
                     'present: constrain the model, or add artificial boundary ' ...
                     'stiffness (Sec. 5.1.1 of the paper).'], ...
                    numel(rb_idx), min(rres), sqrt(max(w2(rb_idx),0))/(2*pi));
            end
            obj.omega2 = w2;

            % ---------- 3. Mass normalization ----------
            for i = 1:m
                Phi(:,i) = Phi(:,i) / sqrt(Phi(:,i)' * Mc * Phi(:,i));
            end

            % Named Phi_bnd locally so it does not shadow the obj.Phi_b property.
            Phi_bnd = Phi(nl_dof, :);    % [n_bnd x m]
            Phi_i   = Phi(inner_idx, :); % [n_inn x m]
            obj.Phi_b = Phi_bnd;

            % ---------- 4. Attachment modes and flexibility ----------
            % Columns of F = inv(K) corresponding to the boundary DOFs.
            F_int     = sparse(nl_dof, 1:obj.n_bnd, 1, n_dofs_c, obj.n_bnd);
            Flex_cols = Kc \ full(F_int);

            F_bb = Flex_cols(nl_dof, :);      % [n_bnd x n_bnd]
            F_ib = Flex_cols(inner_idx, :);   % [n_inn x n_bnd]

            % RESIDUAL flexibility, Eq. (5.3): subtract the contribution of
            % the modes already retained in the basis.
            inv_w2   = diag(1 ./ w2);
            F_bb_res = F_bb - Phi_bnd * inv_w2 * Phi_bnd';
            F_ib_res = F_ib - Phi_i   * inv_w2 * Phi_bnd';

            F_bb_res = (F_bb_res + F_bb_res') / 2;   % must be symmetric
            obj.Fbb_res = F_bb_res;

            rc = rcond(F_bb_res);
            fprintf('  rcond(F_bb_res) = %.3e\n', rc);
            if rc < 1e-12
                warning('RomMN:IllConditionedFbb', ...
                    ['F_bb_res is ill-conditioned (rcond = %.2e). This usually means that ' ...
                     'numModes is too high for the residual flexibility left (the retained ' ...
                     'modes have already saturated F_bb), or that n_bnd is too large.'], rc);
            end

            % ---------- 5. Component mode matrix, Eq. (5.4) ----------
            T_ib = F_ib_res / F_bb_res;      % F'_ib * inv(F'_bb)

            Pc_matrix = zeros(n_dofs_c, obj.n_bnd + m);
            Pc_matrix(nl_dof,    1:obj.n_bnd)     = eye(obj.n_bnd);
            Pc_matrix(inner_idx, 1:obj.n_bnd)     = T_ib;
            Pc_matrix(nl_dof,    obj.n_bnd+1:end) = zeros(obj.n_bnd, m);
            Pc_matrix(inner_idx, obj.n_bnd+1:end) = Phi_i - T_ib * Phi_bnd;

            obj.Pc = Pc_matrix;
            obj.P  = obj.Structure.AssemblyObj.unconstrain_vector(Pc_matrix);

            % ---------- 6. Reduced matrices, Eq. (5.5) ----------
            % Assembled EXPLICITLY, not as Pc'*Kc*Pc, which would be Rubin.
            Fbb_inv = inv(F_bb_res);
            Fbb_inv = (Fbb_inv + Fbb_inv') / 2;

            K_bb_r = Fbb_inv;
            K_bi_r = -Fbb_inv * Phi_bnd;
            K_ii_r = diag(w2) + Phi_bnd' * Fbb_inv * Phi_bnd;

            obj.K_r = [ K_bb_r , K_bi_r ;
                        K_bi_r', K_ii_r ];
            obj.K_r = (obj.K_r + obj.K_r') / 2;

            % Mass: identity on the modal coordinates, ZERO at the boundary
            obj.M_r = zeros(obj.n_bnd + m);
            obj.M_r(obj.n_bnd+1:end, obj.n_bnd+1:end) = eye(m);

            % Sanity check: statically condensing the boundary out (open
            % contact) must reproduce the FOM frequencies.
            nb = obj.n_bnd;
            Kbb_t = obj.K_r(1:nb, 1:nb);
            Kbe_t = obj.K_r(1:nb, nb+1:end);
            Kee_t = obj.K_r(nb+1:end, nb+1:end);

            K_cond = Kee_t - Kbe_t' * (Kbb_t \ Kbe_t);
            K_cond = (K_cond + K_cond')/2;

            f_cond = sort(sqrt(abs(eig(K_cond))))/(2*pi);   % M = I on the modal block
            % The reference list is only as long as the modal analysis that was
            % actually run, which is not guaranteed to reach five modes.
            n_show = min([5, numel(f_cond), numel(obj.Structure.frequencies)]);
            if n_show > 0
                fprintf('\n  --- Condensed ROM (open contact, M = I) ---\n');
                for i = 1:n_show
                    fprintf('  f%d: ROM = %.4e | FOM = %.4e | err = %+.1f%%\n', ...
                        i, f_cond(i), obj.Structure.frequencies(i), ...
                        100*(f_cond(i)/obj.Structure.frequencies(i) - 1));
                end
            end

            % ---------- 7. Modal damping ----------
            w = sqrt(w2);

            if isstruct(zeta_modal)
                % Rayleigh: zeta_k = 0.5*(alpha/w_k + beta*w_k)
                zvec = 0.5 * (zeta_modal.alpha ./ w + zeta_modal.beta .* w);
            elseif isscalar(zeta_modal)
                zvec = repmat(zeta_modal, m, 1);
            else
                zvec = zeta_modal(:);
                assert(numel(zvec) == m, ...
                    'zeta_modal must be a scalar, an [m x 1] vector, or a Rayleigh struct.');
            end
            obj.zeta = zvec;

            obj.C_r = zeros(obj.n_bnd + m);
            obj.C_r(obj.n_bnd+1:end, obj.n_bnd+1:end) = diag(2 * zvec .* w);

            fprintf('  zeta: min = %.3e (mode %d), max = %.3e (mode %d)\n', ...
                min(zvec), find(zvec==min(zvec),1), max(zvec), find(zvec==max(zvec),1));
        end

        function [Mr, Kr, Cr] = get_reduced_matrices(obj)
            Mr = obj.M_r;  Kr = obj.K_r;  Cr = obj.C_r;
        end

        function check(obj)
            % CHECK Verify the properties the MacNeal ROM must satisfy.
            fprintf('\n--- RomMN check ---\n');
            nb = obj.n_bnd;

            % (a) Mass must be EXACTLY [0 0; 0 I]
            fprintf('  ||M_r(bnd,:)||       = %.3e   (expected 0)\n', ...
                norm(obj.M_r(1:nb,:), 'fro'));
            fprintf('  ||M_r(inn,inn) - I|| = %.3e   (expected 0)\n', ...
                norm(obj.M_r(nb+1:end, nb+1:end) - eye(obj.numModes), 'fro'));

            % (b) Reduced K_bb must equal inv(F'_bb), hence be SPD
            eK = eig(obj.K_r(1:nb, 1:nb));
            fprintf('  min eig K_r(bnd,bnd) = %.3e   (expected > 0)\n', min(eK));

            % (c) Global K_r must be positive definite
            eKg = eig(obj.K_r);
            fprintf('  min eig K_r global   = %.3e   (expected > 0)\n', min(eKg));

            % (d) Boundary-interior elastic coupling must NOT vanish, unlike
            %     standard CB. A zero here means something went wrong.
            fprintf('  ||K_r(bnd,inn)||     = %.3e   (expected > 0)\n', ...
                norm(obj.K_r(1:nb, nb+1:end), 'fro'));
        end
    end
end
