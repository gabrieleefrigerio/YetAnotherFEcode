classdef RomRubinTest < handle
    % ROMRUBINTEST  Rubin CMS ROM con scaling diagonale opzionale della base.
    %
    %   Lo scaling divide ogni colonna di Pc per la propria norma in massa.
    %   E' un cambio base diagonale invertibile: span, spettro del pencil
    %   ridotto (Kr,Mr) e quindi la soluzione fisica sono INVARIANTI.
    %   Cambia solo il condizionamento della parametrizzazione, che governa
    %   il costo online tramite la matrice di iterazione di ode15s
    %       W = M_state/(h*gamma_k) - J        con   kappa(W) -> kappa(Mr) per h->0.
    %
    %   Dopo lo scaling la relazione interfaccia <-> coordinate resta DIAGONALE:
    %       Pc(contactDofs, 1:n_bnd) = diag(1 ./ scaleD(1:n_bnd))
    %   quindi  q_fisico_i = y_i / scaleD(i).  L'interfaccia resta uno-a-uno
    %   (non diventa una proiezione densa come in Milman-Chu), ma NON e' piu'
    %   numericamente uguale allo spostamento fisico: con applyScaling = true
    %   il solver va lanciato sul ramo proiettato ('ModelType','MC').
    %
    %   Uso tipico:
    %       rom = RomRubinTest(DummyStruct, phi, contact_dofs);
    %       rom.build();                  % costruisce Pc_raw e applica il default
    %       rom.set_scaling(false);       % configurazione grezza
    %       rom.set_scaling(true);        % configurazione scalata (nessun rebuild)

    properties
        Structure
        P               % base nello spazio unconstrained
        Pc              % base nello spazio constrained (scalata o grezza)
        Pc_raw          % base constrained SEMPRE grezza (riferimento)
        numModes
        contactDofs
        n_bnd
        scaleD          % fattori di scala per colonna (1 se scaling disattivo)
        applyScaling = true
        diag_info = struct()
    end

    properties (Access = private)
        Mc_             % cache matrice di massa constrained
        Kc_             % cache matrice di rigidezza constrained
    end

    methods
        function obj = RomRubinTest(dummy_struct, num_linear_modes, ...
                                    contact_dofs_constrained, apply_scaling)
            obj.Structure   = dummy_struct;
            obj.numModes    = num_linear_modes;
            obj.contactDofs = contact_dofs_constrained(:);   % forzatura a colonna
            obj.n_bnd       = length(obj.contactDofs);

            if nargin >= 4 && ~isempty(apply_scaling)
                obj.applyScaling = logical(apply_scaling);
            end

            % DOF di contatto duplicati generano colonne identiche in F_int
            % e quindi una T1_bb singolare: meglio accorgersene subito.
            assert(numel(unique(obj.contactDofs)) == obj.n_bnd, ...
                   'RomRubinTest: contactDofs contiene DOF duplicati.');
        end

        % -----------------------------------------------------------------
        function build(obj)
            obj.Mc_ = obj.Structure.AssemblyObj.constrain_matrix(obj.Structure.M);
            obj.Kc_ = obj.Structure.AssemblyObj.constrain_matrix(obj.Structure.K);
            Mc = obj.Mc_;  Kc = obj.Kc_;
            n_dofs_c = size(Kc, 1);

            fprintf('\n--- Building Classic Rubin ROM Base (CMS) ---\n');

            % --- modi normali a interfaccia libera -------------------------
            [Phi_k, D] = eigs(Kc, Mc, obj.numModes, 'smallestabs');
            [~, sort_idx] = sort(diag(D));
            Phi_k = Phi_k(:, sort_idx);
            D     = D(sort_idx, sort_idx);
            omega2 = diag(D);

            for i = 1:size(Phi_k, 2)
                Phi_k(:, i) = Phi_k(:, i) / sqrt(Phi_k(:, i)' * Mc * Phi_k(:, i));
            end

            % --- flessibilita' residua ------------------------------------
            % La colonna j di F_int corrisponde a contactDofs(j): l'ordine
            % dell'interfaccia e' fissato qui e propagato a Psi_res, T1_bb,
            % alle coordinate ridotte e a scaleD.
            F_int = sparse(obj.contactDofs, 1:obj.n_bnd, 1, n_dofs_c, obj.n_bnd);

            Psi_a   = Kc \ F_int;
            Psi_c   = Phi_k * (diag(1 ./ omega2) * (Phi_k' * F_int));
            Psi_res = Psi_a - Psi_c;

            T1   = [Psi_res, Phi_k];
            T1_b = T1(obj.contactDofs, :);

            T1_bb = T1_b(:, 1:obj.n_bnd);            % flessibilita' residua d'interfaccia
            T1_bi = T1_b(:, obj.n_bnd+1:end);

            % T2 senza inversa esplicita (stesso risultato, meglio condizionato)
            W  = T1_bb \ [eye(obj.n_bnd), -T1_bi];   % n_bnd x (n_bnd + numModes)
            T2 = [W; zeros(obj.numModes, obj.n_bnd), eye(obj.numModes)];

            obj.Pc_raw = T1 * T2;

            % diagnostica offline dell'interfaccia (serve per interface reduction)
            Gres = full(T1_bb);  Gres = (Gres + Gres') / 2;
            s_res = sort(abs(eig(Gres)), 'descend');
            obj.diag_info.cond_T1bb   = cond(Gres);
            obj.diag_info.sigma_iface = s_res / s_res(1);

            fprintf('  r = %d  (%d interfaccia + %d modi) | cond(T1_bb) = %.3e\n', ...
                    size(obj.Pc_raw, 2), obj.n_bnd, obj.numModes, obj.diag_info.cond_T1bb);

            obj.set_scaling(obj.applyScaling);
            fprintf('--- Rubin ROM Build Complete ---\n');
        end

        % -----------------------------------------------------------------
        function set_scaling(obj, tf)
            % Attiva/disattiva lo scaling senza ricostruire la base.
            if nargin > 1, obj.applyScaling = logical(tf); end
            assert(~isempty(obj.Pc_raw), 'Chiamare build() prima di set_scaling().');

            Mc = obj.Mc_;  Kc = obj.Kc_;
            r  = size(obj.Pc_raw, 2);

            if obj.applyScaling
                % norma in massa di ogni colonna
                d = sqrt(full(sum(obj.Pc_raw .* (Mc * obj.Pc_raw), 1)));
                obj.scaleD = d(:);
                obj.Pc = obj.Pc_raw ./ d;
            else
                obj.scaleD = ones(r, 1);
                obj.Pc = obj.Pc_raw;
            end

            obj.P = obj.Structure.AssemblyObj.unconstrain_vector(obj.Pc);

            % --- diagnostica di condizionamento ---------------------------
            Mr_d = full(obj.Pc' * Mc * obj.Pc);  Mr_d = (Mr_d + Mr_d') / 2;
            Kr_d = full(obj.Pc' * Kc * obj.Pc);  Kr_d = (Kr_d + Kr_d') / 2;

            obj.diag_info.applyScaling = obj.applyScaling;
            obj.diag_info.r            = r;
            obj.diag_info.spread_d     = max(obj.scaleD) / min(obj.scaleD);
            obj.diag_info.cond_Mr      = cond(Mr_d);
            obj.diag_info.cond_Kr      = cond(Kr_d);

            % --- check di ordinamento dell'interfaccia --------------------
            % Per costruzione di T2 il blocco Pc(contactDofs, 1:n_bnd) deve
            % essere diagonale (= I senza scaling, = diag(1./d) con scaling).
            % Se l'ordine di contactDofs si e' desincronizzato in qualunque
            % punto della catena, qui esce una permutazione e l'assert scatta.
            Bchk = full(obj.Pc(obj.contactDofs, 1:obj.n_bnd));
            err_ord = norm(Bchk - diag(diag(Bchk)), 'fro') / norm(diag(Bchk));
            obj.diag_info.err_ordering = err_ord;

            fprintf('  [scaling %d] spread d = %.2e | cond(Mr) = %.3e | cond(Kr) = %.3e | ord %.1e\n', ...
                    obj.applyScaling, obj.diag_info.spread_d, ...
                    obj.diag_info.cond_Mr, obj.diag_info.cond_Kr, err_ord);

            assert(err_ord < 1e-8, ...
                'RomRubinTest: ordinamento interfaccia incoerente (err = %.2e).', err_ord);
        end

        % -----------------------------------------------------------------
        function [Mr, Kr, Cr] = get_reduced_matrices(obj)
            Mc = obj.Mc_;  Kc = obj.Kc_;
            Cc = obj.Structure.AssemblyObj.constrain_matrix(obj.Structure.C);

            Mr = obj.Pc' * Mc * obj.Pc;
            Kr = obj.Pc' * Kc * obj.Pc;
            if isscalar(Cc) && Cc == 0
                Cr = sparse(size(Mr, 1), size(Mr, 2));
            else
                Cr = obj.Pc' * Cc * obj.Pc;
            end

            Mr = (Mr + Mr') / 2;
            Kr = (Kr + Kr') / 2;
            Cr = (Cr + Cr') / 2;
        end

        % -----------------------------------------------------------------
        function [q0_r, qd0_r] = project_ic(obj, q0_c, qd0_c)
            % Proiezione M-ortogonale delle condizioni iniziali.
            % Valida per QUALUNQUE base, anche quando Mr ~= I.
            Mr   = obj.Pc' * obj.Mc_ * obj.Pc;  Mr = (Mr + Mr') / 2;
            rhs  = obj.Pc' * (obj.Mc_ * q0_c);
            q0_r = Mr \ rhs;
            if nargout > 1
                qd0_r = Mr \ (obj.Pc' * (obj.Mc_ * qd0_c));
            end
        end
    end
end
