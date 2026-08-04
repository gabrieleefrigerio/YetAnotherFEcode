classdef RomMN < handle
    % MACNEAL free-interface CMS con boundary massless.
    % Riferimento: Monjaraz Tec et al., Eq. (5.3)-(5.5).
    %
    % IMPORTANTE: MacNeal e' un metodo INCONSISTENTE (non-Galerkin).
    % K_r NON e' Pc'*Kc*Pc (quello sarebbe Rubin). Va costruita
    % esplicitamente dalla flessibilita' residua.

    properties
        Structure
        P, Pc
        numModes
        contactDofs
        n_bnd
        M_r, K_r, C_r
        omega2          % autovalori dei modi free-interface ritenuti
        Phi_b           % modi ristretti al boundary
        Fbb_res         % flessibilita' residua al boundary
        zeta            % damping ratio modale usato
    end

    methods
        function obj = RomMN(dummy_struct, num_linear_modes, contact_dofs_constrained)
            obj.Structure   = dummy_struct;
            obj.numModes    = num_linear_modes;
            obj.contactDofs = contact_dofs_constrained(:);
            obj.n_bnd       = length(obj.contactDofs);
        end

        function build(obj, zeta_modal)
            % zeta_modal : scalare (damping ratio uniforme) oppure vettore [m x 1].
            %              Se omesso -> 0.
            if nargin < 2, zeta_modal = 0; end

            Mc = obj.Structure.AssemblyObj.constrain_matrix(obj.Structure.M);
            Kc = obj.Structure.AssemblyObj.constrain_matrix(obj.Structure.K);
            n_dofs_c = size(Kc, 1);

            % Simmetrizzazione difensiva
            Kc = (Kc + Kc') / 2;
            Mc = (Mc + Mc') / 2;

            fprintf('\n--- Building MacNeal ROM (free-interface, massless boundary) ---\n');

            nl_dof    = obj.contactDofs;
            inner_idx = setdiff(1:n_dofs_c, nl_dof)';
            m         = obj.numModes;

            % ---------- 1. Modi normali free-interface ----------
            [Phi, D] = eigs(Kc, Mc, m, 'smallestabs');

            if ~isreal(Phi) || ~isreal(D)
                max_imag = max(max(abs(imag(Phi(:)))), max(abs(imag(diag(D)))));
                if max_imag > 1e-8
                    warning('RomMN:ComplexEigs', ...
                        'eigs ha restituito parte immaginaria %.3e. Verificare simmetria M,K.', max_imag);
                end
                Phi = real(Phi);  D = real(D);
            end

            [w2, sort_idx] = sort(diag(D));
            Phi = Phi(:, sort_idx);

            % ---------- 2. Guardia sui modi rigidi ----------
            % MacNeal richiede K non singolare (serve F = inv(K)) e omega != 0.
            w2_scale = max(abs(w2));
            rb_idx = find(w2 < 1e-8 * w2_scale);
            if ~isempty(rb_idx)
                error('RomMN:RigidBodyModes', ...
                    ['Rilevati %d modi a frequenza ~nulla (omega^2 = %.3e). ' ...
                     'MacNeal non e'' applicabile direttamente in presenza di moti rigidi. ' ...
                     'Rimuovere i moti rigidi o aggiungere rigidezza artificiale al boundary ' ...
                     '(vedi Sez. 5.1.1 del paper).'], numel(rb_idx), w2(1));
            end
            obj.omega2 = w2;

            % ---------- 3. Mass normalization ----------
            for i = 1:m
                Phi(:,i) = Phi(:,i) / sqrt(Phi(:,i)' * Mc * Phi(:,i));
            end

            Phi_b = Phi(nl_dof, :);      % [n_bnd x m]
            Phi_i = Phi(inner_idx, :);   % [n_inn x m]
            obj.Phi_b = Phi_b;

            % ---------- 4. Attachment modes / flessibilita' ----------
            % Colonne di F = inv(K) corrispondenti ai boundary DOF.
            F_int     = sparse(nl_dof, 1:obj.n_bnd, 1, n_dofs_c, obj.n_bnd);
            Flex_cols = Kc \ full(F_int);

            F_bb = Flex_cols(nl_dof, :);      % [n_bnd x n_bnd]
            F_ib = Flex_cols(inner_idx, :);   % [n_inn x n_bnd]

            % Flessibilita' RESIDUA (Eq. 5.3): sottrai il contributo dei modi ritenuti
            inv_w2   = diag(1 ./ w2);
            F_bb_res = F_bb - Phi_b * inv_w2 * Phi_b';
            F_ib_res = F_ib - Phi_i * inv_w2 * Phi_b';

            F_bb_res = (F_bb_res + F_bb_res') / 2;   % dev'essere simmetrica
            obj.Fbb_res = F_bb_res;

            rc = rcond(F_bb_res);
            fprintf('  rcond(F_bb_res) = %.3e\n', rc);
            if rc < 1e-12
                warning('RomMN:IllConditionedFbb', ...
                    ['F_bb_res mal condizionata (rcond = %.2e). Tipicamente significa che ' ...
                     'numModes e'' troppo alto rispetto alla flessibilita'' residua disponibile ' ...
                     '(i modi ritenuti hanno gia'' saturato F_bb), oppure n_bnd e'' troppo grande.'], rc);
            end

            % ---------- 5. Matrice dei component modes (Eq. 5.4) ----------
            T_ib = F_ib_res / F_bb_res;      % F'_ib * inv(F'_bb)

            Pc_matrix = zeros(n_dofs_c, obj.n_bnd + m);
            Pc_matrix(nl_dof,    1:obj.n_bnd)     = eye(obj.n_bnd);
            Pc_matrix(inner_idx, 1:obj.n_bnd)     = T_ib;
            Pc_matrix(nl_dof,    obj.n_bnd+1:end) = zeros(obj.n_bnd, m);
            Pc_matrix(inner_idx, obj.n_bnd+1:end) = Phi_i - T_ib * Phi_b;

            obj.Pc = Pc_matrix;
            obj.P  = obj.Structure.AssemblyObj.unconstrain_vector(Pc_matrix);

            % ---------- 6. Matrici ridotte (Eq. 5.5) ----------
            % NB: costruzione ESPLICITA, non Pc'*Kc*Pc (quello sarebbe Rubin).
            Fbb_inv = inv(F_bb_res);
            Fbb_inv = (Fbb_inv + Fbb_inv') / 2;

            K_bb_r = Fbb_inv;
            K_bi_r = -Fbb_inv * Phi_b;
            K_ii_r = diag(w2) + Phi_b' * Fbb_inv * Phi_b;

            obj.K_r = [ K_bb_r , K_bi_r ;
                        K_bi_r', K_ii_r ];
            obj.K_r = (obj.K_r + obj.K_r') / 2;

            % Massa: identita' sui modali, ZERO al boundary
            obj.M_r = zeros(obj.n_bnd + m);
            obj.M_r(obj.n_bnd+1:end, obj.n_bnd+1:end) = eye(m);
nb = obj.n_bnd;
Kbb_t = obj.K_r(1:nb, 1:nb);
Kbe_t = obj.K_r(1:nb, nb+1:end);
Kee_t = obj.K_r(nb+1:end, nb+1:end);

K_cond = Kee_t - Kbe_t' * (Kbb_t \ Kbe_t);
K_cond = (K_cond + K_cond')/2;

f_cond = sort(sqrt(abs(eig(K_cond))))/(2*pi);   % M = I sui modali
fprintf('\n  --- ROM condensato (contatto aperto, M=I) ---\n');
for i = 1:min(5, numel(f_cond))
    fprintf('  f%d: ROM = %.4e | FOM = %.4e | err = %+.1f%%\n', ...
        i, f_cond(i), obj.Structure.frequencies(i), ...
        100*(f_cond(i)/obj.Structure.frequencies(i) - 1));
end
            % ---------- 7. Damping modale ----------
            % zeta_modal puo' essere:
            %   - scalare            -> zeta uniforme
            %   - vettore [m x 1]    -> zeta per modo
            %   - struct('alpha',a,'beta',b) -> Rayleigh: zeta_k = 0.5*(a/w_k + b*w_k)
            w = sqrt(w2);

            if isstruct(zeta_modal)
                zvec = 0.5 * (zeta_modal.alpha ./ w + zeta_modal.beta .* w);
            elseif isscalar(zeta_modal)
                zvec = repmat(zeta_modal, m, 1);
            else
                zvec = zeta_modal(:);
                assert(numel(zvec) == m, 'zeta_modal deve essere scalare, vettore [m x 1], o struct Rayleigh.');
            end
            obj.zeta = zvec;

            obj.C_r = zeros(obj.n_bnd + m);
            obj.C_r(obj.n_bnd+1:end, obj.n_bnd+1:end) = diag(2 * zvec .* w);

            fprintf('  zeta: min = %.3e (modo %d), max = %.3e (modo %d)\n', ...
                min(zvec), find(zvec==min(zvec),1), max(zvec), find(zvec==max(zvec),1));
        end

        function [Mr, Kr, Cr] = get_reduced_matrices(obj)
            Mr = obj.M_r;  Kr = obj.K_r;  Cr = obj.C_r;
        end

        function check(obj)
            % Diagnostica: verifica le proprieta' attese del ROM MacNeal.
            fprintf('\n--- RomMN check ---\n');
            nb = obj.n_bnd;

            % (a) La massa deve essere ESATTAMENTE [0 0; 0 I]
            fprintf('  ||M_r(bnd,:)||       = %.3e   (atteso 0)\n', ...
                norm(obj.M_r(1:nb,:), 'fro'));
            fprintf('  ||M_r(inn,inn) - I|| = %.3e   (atteso 0)\n', ...
                norm(obj.M_r(nb+1:end, nb+1:end) - eye(obj.numModes), 'fro'));

            % (b) K_bb ridotta deve coincidere con inv(F'_bb): SPD
            eK = eig(obj.K_r(1:nb, 1:nb));
            fprintf('  min eig K_r(bnd,bnd) = %.3e   (atteso > 0)\n', min(eK));

            % (c) K_r globale deve essere definita positiva
            eKg = eig(obj.K_r);
            fprintf('  min eig K_r globale  = %.3e   (atteso > 0)\n', min(eKg));

            % (d) Accoppiamento elastico bordo-interno NON deve essere nullo
            %     (a differenza del CB standard). Se e' 0, c'e' un errore.
            fprintf('  ||K_r(bnd,inn)||     = %.3e   (atteso > 0)\n', ...
                norm(obj.K_r(1:nb, nb+1:end), 'fro'));
        end
    end
end