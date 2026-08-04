classdef RomMCB < handle
    % Massless Craig-Bampton (fixed-interface CMS con boundary massless).
    % Riferimento: Monjaraz Tec et al., Sez. 5.1.2, Eq. (5.6)-(5.10).
    %
    % A differenza di MacNeal, qui K_r = R'KR e' consistente (Galerkin).
    % L'inconsistenza e' introdotta SOLO azzerando M*_bb.

    properties
        Structure, P, Pc, numModes, contactDofs, n_bnd
        M_r, K_r, C_r
        alpha, omega2, zeta
    end

    methods
        function obj = RomMCB(dummy_struct, num_fixed_modes, contact_dofs_constrained)
            obj.Structure   = dummy_struct;
            obj.numModes    = num_fixed_modes;
            obj.contactDofs = contact_dofs_constrained(:);
            obj.n_bnd       = length(obj.contactDofs);
        end

        function build(obj, zeta_modal)
            % zeta_modal : scalare | vettore [m x 1] | struct('alpha',a,'beta',b)
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

            % --- Modi fixed-interface ---
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

            % --- Trasformazione alpha: disaccoppia inerzialmente bordo/interno (Eq. 5.9) ---
            alpha = Phi_i' * (full(M_ib) + M_ii * Psi_c);
            obj.alpha = alpha;

            % --- Base R_alpha (Eq. 5.8) ---
            P_alpha = zeros(n_dofs_c, obj.n_bnd + m);
            P_alpha(nl_dof_c,    1:obj.n_bnd)     = eye(obj.n_bnd);
            P_alpha(inner_idx_c, 1:obj.n_bnd)     = Psi_c - Phi_i * alpha;
            P_alpha(inner_idx_c, obj.n_bnd+1:end) = Phi_i;

            obj.Pc = P_alpha;
            obj.P  = obj.Structure.AssemblyObj.unconstrain_vector(P_alpha);

            % --- Matrici ridotte ---
            % K_r e' Galerkin (a differenza di MacNeal)
            obj.K_r = P_alpha' * Kc * P_alpha;
            obj.K_r = (obj.K_r + obj.K_r')/2;

            M_complete = P_alpha' * Mc * P_alpha;

            % Verifica che alpha abbia disaccoppiato il bordo
            coupling = norm(M_complete(1:obj.n_bnd, obj.n_bnd+1:end), 'fro');
            scale    = norm(M_complete, 'fro');
            fprintf('  Residuo accoppiamento M_bi: %.2e (rel: %.2e)\n', coupling, coupling/scale);
            if coupling/scale > 1e-8
                warning('RomMCB:CouplingNotZero', ...
                    ['La trasformazione alpha non ha annullato M_bi (rel = %.2e). ' ...
                     'Verificare la mass-normalization dei modi.'], coupling/scale);
            end

            % Massa: ESPLICITAMENTE [0 0 ; 0 I]  (Eq. 5.10, secondo passo)
            obj.M_r = zeros(obj.n_bnd + m);
            obj.M_r(obj.n_bnd+1:end, obj.n_bnd+1:end) = ...
                M_complete(obj.n_bnd+1:end, obj.n_bnd+1:end);

            % --- Damping modale (solo coordinate interne) ---
            w = sqrt(w2);
            if isstruct(zeta_modal)
                zvec = 0.5 * (zeta_modal.alpha ./ w + zeta_modal.beta .* w);
            elseif isscalar(zeta_modal)
                zvec = repmat(zeta_modal, m, 1);
            else
                zvec = zeta_modal(:);
                assert(numel(zvec) == m, 'zeta_modal: scalare, [m x 1], o struct Rayleigh.');
            end
            obj.zeta = zvec;

            obj.C_r = zeros(obj.n_bnd + m);
            obj.C_r(obj.n_bnd+1:end, obj.n_bnd+1:end) = diag(2 * zvec .* w);

            fprintf('  Base: %d boundary DOFs + %d modal DOFs = %d totali\n', ...
                obj.n_bnd, m, obj.n_bnd + m);
            fprintf('  zeta: min = %.3e | max = %.3e\n', min(zvec), max(zvec));
            fprintf('  Range freq. fixed-interface: %.3e - %.3e Hz\n', ...
                w(1)/(2*pi), w(end)/(2*pi));
        end

        function [Mr, Kr, Cr] = get_reduced_matrices(obj)
            Mr = obj.M_r;  Kr = obj.K_r;  Cr = obj.C_r;
        end

        function check(obj)
            fprintf('\n--- RomMCB check ---\n');
            nb = obj.n_bnd;
            fprintf('  ||M_r(bnd,:)||       = %.3e   (atteso 0)\n', norm(obj.M_r(1:nb,:), 'fro'));
            fprintf('  ||M_r(inn,inn) - I|| = %.3e   (atteso ~0)\n', ...
                norm(obj.M_r(nb+1:end,nb+1:end) - eye(obj.numModes), 'fro'));
            fprintf('  min eig K_r(bnd,bnd) = %.3e   (atteso > 0)\n', min(eig(obj.K_r(1:nb,1:nb))));
            fprintf('  min eig K_r globale  = %.3e   (atteso > 0)\n', min(eig(obj.K_r)));
            % Nel CB massless l'accoppiamento ELASTICO bordo-interno e' NON nullo
            % (al contrario del CB standard, dove K_bi = 0).
            fprintf('  ||K_r(bnd,inn)||     = %.3e   (atteso > 0)\n', ...
                norm(obj.K_r(1:nb, nb+1:end), 'fro'));
        end
    end
end