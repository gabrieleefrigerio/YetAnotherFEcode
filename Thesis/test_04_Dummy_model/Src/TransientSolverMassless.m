classdef TransientSolverMassless < handle
    % Integratore semi-esplicito (leapfrog/Verlet) per modelli a boundary
    % massless con contatto unilaterale FRICTIONLESS set-valued.
    %
    % Riferimento: Monjaraz Tec et al., "A massless boundary component mode
    % synthesis method for elastodynamic contact problems", Comput. Struct. 260 (2022).
    % Cap. 4 (time stepping), Cap. 6 (algoritmo frictionless), App. C (aug. Lagrangian).
    %
    % Struttura del modello richiesta (MacNeal o massless CB):
    %   coordinate q = [q_b ; eta],  q_b = boundary (n_bnd),  eta = modali (m)
    %   M = [0 0 ; 0 I]        (massa NULLA al boundary)
    %   K = [K_bb K_be ; K_eb K_ee]
    %   C = [0 0 ; 0 D_ee]     (nessun damping al boundary)
    %
    % Equazioni risolte:
    %   K_bb q_b + K_be eta - W lambda = f_b(t)             (statica, al boundary)
    %   eta_dd + D_ee eta_d + K_ee eta + K_eb q_b = f_e(t)  (dinamica, interna)
    %   g = g0 + W' q_b ,   0 <= g  _|_  lambda >= 0
    %
    % Convenzione gap: per muro a destra, con q_b positivo verso il muro:
    %   W = -I ,  g0 = gap_wall   =>   g = gap_wall - q_b

    properties
        % --- matrici ridotte ---
        Kbb, Kbe, Kee, Dee
        n_bnd, n_mod

        % --- contatto ---
        W          % [n_bnd x n_c] matrice direzioni di contatto
        g0         % [n_c x 1] gap iniziale
        n_c

        % --- precomputazioni ---
        Kbb_fact   % fattorizzazione Cholesky di K_bb
        G_full     % W' * inv(K_bb) * W   [n_c x n_c]  (matrice di Delassus)
        Ainv       % inv( (1/dt)*I + 0.5*Dee )
        Bmat       % (1/dt)*I - 0.5*Dee

        % --- parametri augmented Lagrangian ---
        eps_AL              % passo di rilassamento
        max_iter_AL = 1000
        tol_AL      = 1e-8  % RELATIVA a ||c||_inf (residuo KKT normalizzato)

        % --- diagnostica ---
        stats
        warned_AL = false
    end

    methods
        function obj = TransientSolverMassless(Mr, Kr, Cr, n_bnd, W, g0)
            % Mr, Kr, Cr : matrici ridotte del ROM massless
            % n_bnd      : numero di DOF di boundary (primi n_bnd della base)
            % W          : [n_bnd x n_c] direzioni di contatto
            % g0         : [n_c x 1] gap iniziali

            n_tot = size(Kr, 1);
            obj.n_bnd = n_bnd;
            obj.n_mod = n_tot - n_bnd;

            ib = 1:n_bnd;
            ie = n_bnd+1 : n_tot;

            % ---------- validazione della struttura massless ----------
            if norm(Mr(ib, :), 'fro') > 1e-10 * (norm(Mr, 'fro') + eps)
                error('TSM:NotMassless', ...
                    ['La matrice di massa ha termini non nulli sulle righe di boundary ' ...
                     '(||M(bnd,:)|| = %.3e). Il modello NON e'' massless.'], ...
                     norm(Mr(ib,:), 'fro'));
            end

            Mee = full(Mr(ie, ie));
            devI = norm(Mee - eye(obj.n_mod), 'fro');
            if devI > 1e-8 * obj.n_mod
                warning('TSM:MassNotIdentity', ...
                    ['M(inn,inn) non e'' l''identita'' (dev = %.3e). Lo schema assume ' ...
                     'modi normalizzati in massa.'], devI);
            end

            if norm(Cr(ib, :), 'fro') > 1e-10 * (norm(Cr, 'fro') + eps)
                warning('TSM:BoundaryDamping', ...
                    ['La matrice di damping ha termini al boundary: verranno IGNORATI ' ...
                     '(la formulazione massless non li prevede).']);
            end

            % ---------- partizioni ----------
            obj.Kbb = full(Kr(ib, ib));
            obj.Kbe = full(Kr(ib, ie));
            obj.Kee = full(Kr(ie, ie));
            obj.Dee = full(Cr(ie, ie));

            obj.Kbb = (obj.Kbb + obj.Kbb') / 2;

            obj.W   = W;
            obj.g0  = g0(:);
            obj.n_c = size(W, 2);
            assert(size(W,1) == n_bnd,      'W deve avere n_bnd righe.');
            assert(numel(obj.g0) == obj.n_c,'g0 deve avere n_c elementi.');

            % ---------- precomputazioni ----------
            obj.Kbb_fact = decomposition(obj.Kbb, 'chol');

            Kbb_inv_W  = obj.Kbb_fact \ obj.W;        % inv(Kbb)*W
            obj.G_full = obj.W' * Kbb_inv_W;          % matrice di Delassus
            obj.G_full = (obj.G_full + obj.G_full') / 2;

            % ---------- eps_AL ----------
            % Jacobi proiettato converge per  0 < eps < 2/lambda_max(G).
            % Basarsi sullo SPETTRO di G (non sulla sola diagonale) e' robusto
            % anche quando K_r e' mal condizionata.
            if obj.n_c == 1
                lam_max = obj.G_full;
            else
                lam_max = eigs(obj.G_full, 1, 'largestabs', ...
                               'Tolerance', 1e-6, 'MaxIterations', 500);
            end
            obj.eps_AL = 1.0 / lam_max;

            fprintf('  [massless] n_bnd = %d | n_mod = %d | n_c = %d\n', ...
                obj.n_bnd, obj.n_mod, obj.n_c);
            fprintf('  [massless] lambda_max(G) = %.3e | eps_AL = %.3e\n', ...
                lam_max, obj.eps_AL);
        end

        % =================================================================
        function dt_crit = critical_timestep(obj)
            % Limite di stabilita' dello schema esplicito sulle coord. interne.
            % NOTA: il boundary NON contribuisce (risolto quasi-staticamente).
            %
            % Il condensamento statico del boundary modifica la rigidezza vista
            % dai modi:
            %   contatto APERTO : K_eff = Kee - Keb*inv(Kbb)*Kbe   (bordo libero)
            %   contatto CHIUSO : K_eff = Kee                      (bordo bloccato)
            % Si prende il caso piu' restrittivo.

            H = obj.Kbb_fact \ obj.Kbe;                 % inv(Kbb)*Kbe
            Keff_open = obj.Kee - obj.Kbe' * H;
            Keff_open = (Keff_open + Keff_open') / 2;
            Kee_sym   = (obj.Kee + obj.Kee') / 2;

            w2 = max([ max(eig(Keff_open)), max(eig(Kee_sym)) ]);
            dt_crit = 2 / sqrt(w2);
        end

        % =================================================================
        function [t, q, lambda_hist, info] = solve(obj, tmax, dt, q0_full, qd0_full, F_handle)
            % q0_full, qd0_full : CI sull'intero vettore ridotto [q_b ; eta].
            %                     La parte di boundary di qd0 e' ignorata.
            % F_handle          : @(t) -> forza ridotta [n_bnd+n_mod x 1]

            nb = obj.n_bnd;  nm = obj.n_mod;
            n_steps = round(tmax/dt);
            t = (0:n_steps) * dt;

            obj.warned_AL = false;

            % ---------- check stabilita' ----------
            dtc = obj.critical_timestep();
            fprintf('  [massless] dt = %.3e | dt_crit ~ %.3e | ratio = %.3f\n', ...
                dt, dtc, dt/dtc);
            if dt > dtc
                warning('TSM:Unstable', ...
                    ['dt = %.3e SUPERA il limite di stabilita'' stimato %.3e. ' ...
                     'Lo schema divergera''. Ridurre dt o numModes.'], dt, dtc);
            end

            % ---------- operatori dell'update esplicito ----------
            Im = eye(nm);
            A  = (1/dt)*Im + 0.5*obj.Dee;
            obj.Bmat = (1/dt)*Im - 0.5*obj.Dee;
            obj.Ainv = A \ Im;

            % ---------- storage ----------
            q           = zeros(nb + nm, n_steps+1);
            lambda_hist = zeros(obj.n_c, n_steps+1);
            n_active    = zeros(1, n_steps+1);
            iters_AL    = zeros(1, n_steps+1);
            kkt_rel     = zeros(1, n_steps+1);

            % ---------- inizializzazione ----------
            % Leapfrog: eta sulla griglia intera, eta_dot sulla semi-intera.
            % Approssimazione eta_dot^{1/2} = eta_dot(t0)  (Sez. 4.4 del paper).
            eta  = q0_full(nb+1:end);
            etad = qd0_full(nb+1:end);         % = eta_dot^{k-1/2}
            lam  = zeros(obj.n_c, 1);
% ---------- DEBUG: ampiezza della forzante proiettata ----------
            Fp = F_handle(5e-7);              % ~ meta' dello shock (t_shock = 1e-6)
            fprintf('  [debug] ||f_b|| = %.3e | ||f_e|| = %.3e\n', ...
                norm(Fp(1:nb)), norm(Fp(nb+1:end)));
            for k = 0:n_steps
                tk = k*dt;
                Fk = F_handle(tk);
                fb = Fk(1:nb);
                fe = Fk(nb+1:end);

                % ---------- 1. predizione del gap (lambda = 0) ----------
                rhs    = fb - obj.Kbe * eta;
                qb_pre = obj.Kbb_fact \ rhs;              % inv(Kbb)*(fb - Kbe*eta)
                g_pre  = obj.g0 + obj.W' * qb_pre;        % = c del paper

                Ia = find(g_pre <= 0);                    % set attivo
                n_active(k+1) = numel(Ia);

                % ---------- 2. soluzione del contatto ----------
                if isempty(Ia)
                    lam = zeros(obj.n_c, 1);
                    qb  = qb_pre;
                    iters_AL(k+1) = 0;
                    kkt_rel(k+1)  = 0;
                else
                    Ga = obj.G_full(Ia, Ia);
                    ca = g_pre(Ia);

                    lam_a = lam(Ia);                      % warm start

                    [lam_a, nit, res] = obj.solve_lcp(Ga, ca, lam_a);
                    iters_AL(k+1) = nit;
                    kkt_rel(k+1)  = res;

                    lam = zeros(obj.n_c, 1);
                    lam(Ia) = lam_a;

                    % q_b = inv(Kbb)*( fb - Kbe*eta + W*lambda )
                    qb = obj.Kbb_fact \ (rhs + obj.W * lam);
                end

                q(1:nb,     k+1) = qb;
                q(nb+1:end, k+1) = eta;
                lambda_hist(:,k+1) = lam;

                if ~all(isfinite(qb)) || ~all(isfinite(eta))
                    error('TSM:Diverged', ...
                        'Soluzione divergente al passo %d (t = %.3e s).', k, tk);
                end

                if k == n_steps, break; end

                % ---------- 3. update esplicito coordinate interne ----------
                % A*etad^{k+1/2} = fe - Kee*eta - Keb*qb + B*etad^{k-1/2}
                rhs_e = fe - obj.Kee * eta - obj.Kbe' * qb + obj.Bmat * etad;
                etad  = obj.Ainv * rhs_e;

                % ---------- 4. update posizione ----------
                eta = eta + etad * dt;
            end

            % ---------- diagnostica ----------
            act = n_active > 0;
            % ---------- DEBUG: ampiezze della risposta ----------
            fprintf('  [debug] max|q_b| = %.3e | max|eta| = %.3e\n', ...
                max(abs(q(1:nb,:)), [], 'all'), max(abs(q(nb+1:end,:)), [], 'all'));
            obj.stats = struct( ...
                'n_active', n_active, ...
                'iters_AL', iters_AL, ...
                'kkt_rel',  kkt_rel, ...
                'dt_crit',  dtc, ...
                'eps_AL',   obj.eps_AL);
            info = obj.stats;

            fprintf('  [massless] passi con contatto attivo: %d/%d (%.1f%%)\n', ...
                nnz(act), n_steps+1, 100*nnz(act)/(n_steps+1));
            if any(act)
                fprintf('  [massless] iter AL: media %.1f | max %d\n', ...
                    mean(iters_AL(act)), max(iters_AL));
                fprintf('  [massless] residuo KKT rel: max %.3e (tol %.1e)\n', ...
                    max(kkt_rel), obj.tol_AL);
                fprintf('  [massless] contatti attivi: max %d / %d\n', ...
                    max(n_active), obj.n_c);
            end
            fprintf('  [massless] penetrazione max: %.3e  (gap = %.3e)\n', ...
                obj.max_penetration(q), max(obj.g0));
        end

        % =================================================================
        function pen = max_penetration(obj, q)
            qb = q(1:obj.n_bnd, :);
            g  = obj.g0 + obj.W' * qb;
            pen = max(0, -min(g(:)));
        end
    end

    % =====================================================================
    methods (Access = private)
        function [lam, nit, res_rel] = solve_lcp(obj, G, c, lam0)
            % Risolve  0 <= (G*lam + c)  _|_  lam >= 0
            % via augmented Lagrangian + Jacobi proiettato (App. C):
            %   lam <- proj_{R+}( lam - eps_AL*(G*lam + c) )
            %
            % Arresto sul RESIDUO KKT NORMALIZZATO:
            %   r = || min(lam, G*lam + c) ||_inf / ||c||_inf
            % Un criterio sull'incremento di lam e' fragile: vicino alla
            % soluzione l'incremento e' dominato dal round-off e la soglia
            % assoluta non viene mai raggiunta.

            lam = max(lam0, 0);
            e   = obj.eps_AL;

            % scala di riferimento = grandezza del termine noto (gap predetto)
            scale = norm(c, inf);
            if scale < 1e-16
                scale = 1;      % floor: evita divisione per ~0
            end

            res_rel = Inf;

            for nit = 1:obj.max_iter_AL
                r   = G*lam + c;
                res_rel = norm(min(lam, r), inf) / scale;

                if res_rel <= obj.tol_AL
                    return;
                end

                lam = max(lam - e*r, 0);
            end

            % Non convergenza vera: warning UNA SOLA VOLTA per simulazione
            if ~obj.warned_AL
                warning('TSM:ALNoConv', ...
                    ['Aug. Lagrangian: %d iter senza convergenza. ' ...
                     'Residuo KKT relativo = %.3e (tol = %.1e). ' ...
                     'Warning emesso una sola volta per simulazione.'], ...
                    obj.max_iter_AL, res_rel, obj.tol_AL);
                obj.warned_AL = true;
            end
        end
    end
end