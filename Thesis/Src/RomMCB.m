classdef RomMCB < handle
    % ROMMCB Massless Craig-Bampton: fixed-interface CMS with a massless boundary.
    %
    % Reference: Monjaraz Tec et al., Sec. 5.1.2, Eq. (5.6)-(5.10).
    %
    % Unlike MacNeal, here K_r = R'*K*R is consistent (Galerkin). The only
    % inconsistency introduced is zeroing the M*_bb block, which is what makes
    % the boundary massless and lets the contact be solved as an LCP.

    properties
        Structure, P, Pc, numModes, contactDofs, n_bnd
        M_r, K_r, C_r
        alpha           % Inertial decoupling transform, Eq. (5.9)
        omega2          % Eigenvalues of the retained fixed-interface modes
        zeta            % Modal damping ratios actually used
    end

    methods
        function obj = RomMCB(dummy_struct, num_fixed_modes, contact_dofs_constrained)
            obj.Structure   = dummy_struct;
            obj.numModes    = num_fixed_modes;
            obj.contactDofs = contact_dofs_constrained(:);
            obj.n_bnd       = length(obj.contactDofs);
        end

        function build(obj, zeta_modal)
            % BUILD Assemble the massless Craig-Bampton reduced model.
            %   zeta_modal : scalar | vector [m x 1] | struct('alpha',a,'beta',b)
            if nargin < 2, zeta_modal = 0; end

            Mc = obj.Structure.AssemblyObj.constrain_matrix(obj.Structure.M);
            Kc = obj.Structure.AssemblyObj.constrain_matrix(obj.Structure.K);
            n_dofs_c = size(Kc, 1);

            Kc = (Kc + Kc')/2;
            Mc = (Mc + Mc')/2;

            fprintf('\n--- Building Massless Craig-Bampton ROM ---\n');

            nl_dof_c    = obj.contactDofs;
            inner_idx_c = setdiff(1:n_dofs_c, nl_dof_c)';

            K_ii = Kc(inner_idx_c, inner_idx_c);
            K_ib = Kc(inner_idx_c, nl_dof_c);
            M_ii = Mc(inner_idx_c, inner_idx_c);
            M_ib = Mc(inner_idx_c, nl_dof_c);

            m = obj.numModes;

            % --- Fixed-interface modes ---
            [Phi_i, D] = eigs(K_ii, M_ii, m, 'smallestabs');
            if ~isreal(Phi_i) || ~isreal(D)
                Phi_i = real(Phi_i);  D = real(D);
            end
            [w2, sort_idx] = sort(diag(D));
            Phi_i = Phi_i(:, sort_idx);
            obj.omega2 = w2;

            % Mass normalization
            for i = 1:m
                Phi_i(:,i) = Phi_i(:,i) / sqrt(Phi_i(:,i)' * M_ii * Phi_i(:,i));
            end

            % --- Constraint modes ---
            Psi_c = -(K_ii \ full(K_ib));

            % --- Alpha transform: decouples boundary and interior inertially, Eq. (5.9) ---
            % Named alpha_t locally so it does not shadow the obj.alpha property.
            alpha_t = Phi_i' * (full(M_ib) + M_ii * Psi_c);
            obj.alpha = alpha_t;

            % --- Basis R_alpha, Eq. (5.8) ---
            P_alpha = zeros(n_dofs_c, obj.n_bnd + m);
            P_alpha(nl_dof_c,    1:obj.n_bnd)     = eye(obj.n_bnd);
            P_alpha(inner_idx_c, 1:obj.n_bnd)     = Psi_c - Phi_i * alpha_t;
            P_alpha(inner_idx_c, obj.n_bnd+1:end) = Phi_i;

            obj.Pc = P_alpha;
            obj.P  = obj.Structure.AssemblyObj.unconstrain_vector(P_alpha);

            % --- Reduced matrices ---
            % Stiffness is Galerkin here, unlike MacNeal.
            obj.K_r = P_alpha' * Kc * P_alpha;
            obj.K_r = (obj.K_r + obj.K_r')/2;

            M_complete = P_alpha' * Mc * P_alpha;

            % Check that the alpha transform really decoupled the boundary
            coupling = norm(M_complete(1:obj.n_bnd, obj.n_bnd+1:end), 'fro');
            scale    = norm(M_complete, 'fro');
            fprintf('  M_bi coupling residual: %.2e (rel: %.2e)\n', coupling, coupling/scale);
            if coupling/scale > 1e-8
                warning('RomMCB:CouplingNotZero', ...
                    ['The alpha transform did not cancel M_bi (rel = %.2e). ' ...
                     'Check the mass normalization of the modes.'], coupling/scale);
            end

            % Mass is set EXPLICITLY to [0 0 ; 0 I], Eq. (5.10), second step
            obj.M_r = zeros(obj.n_bnd + m);
            obj.M_r(obj.n_bnd+1:end, obj.n_bnd+1:end) = ...
                M_complete(obj.n_bnd+1:end, obj.n_bnd+1:end);

            % --- Modal damping, interior coordinates only ---
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

            fprintf('  Basis: %d boundary DOFs + %d modal DOFs = %d total\n', ...
                obj.n_bnd, m, obj.n_bnd + m);
            fprintf('  zeta: min = %.3e | max = %.3e\n', min(zvec), max(zvec));
            fprintf('  Fixed-interface frequency range: %.3e - %.3e Hz\n', ...
                w(1)/(2*pi), w(end)/(2*pi));
        end

        function [Mr, Kr, Cr] = get_reduced_matrices(obj)
            Mr = obj.M_r;  Kr = obj.K_r;  Cr = obj.C_r;
        end

        function check(obj)
            % CHECK Verify the properties the massless CB model must satisfy.
            fprintf('\n--- RomMCB check ---\n');
            nb = obj.n_bnd;
            fprintf('  ||M_r(bnd,:)||       = %.3e   (expected 0)\n', ...
                norm(obj.M_r(1:nb,:), 'fro'));
            fprintf('  ||M_r(inn,inn) - I|| = %.3e   (expected ~0)\n', ...
                norm(obj.M_r(nb+1:end,nb+1:end) - eye(obj.numModes), 'fro'));
            fprintf('  min eig K_r(bnd,bnd) = %.3e   (expected > 0)\n', ...
                min(eig(obj.K_r(1:nb,1:nb))));
            fprintf('  min eig K_r global   = %.3e   (expected > 0)\n', min(eig(obj.K_r)));
            % In the massless CB the ELASTIC boundary-interior coupling is
            % non-zero, unlike standard CB where K_bi = 0.
            fprintf('  ||K_r(bnd,inn)||     = %.3e   (expected > 0)\n', ...
                norm(obj.K_r(1:nb, nb+1:end), 'fro'));
        end
    end
end
