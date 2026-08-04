classdef TransientSolverMasslessCB < handle
    % TransientSolverMasslessCB
    %
    % Integratore semi-esplicito leapfrog (Verlet) per modelli ridotti a
    % contorno privo di massa (massless boundary), secondo lo schema del
    % Cap. 4 della tesi (Monjaraz Tec, 2025), specializzato per basi
    % ottenute con il metodo Massless Craig-Bampton (Sez. 5.1.2), ma
    % applicabile a qualunque ROM con la stessa struttura a blocchi:
    %
    %   stato = [ qb (n_bnd, coordinate fisiche di contatto) ;
    %             eta (m, coordinate modali) ]
    %
    %   K = [ K_bb  K_bi ;  K_ib  K_ii ]     (accoppiata, generalmente piena)
    %   M = [ 0     0    ;  0     M_ii ]     (massa nulla al contorno)
    %
    % Equazioni risolte (frictionless, Eq. 4.1-4.2, Fig. 4.2a):
    %
    %   K_bb*qb^k + K_bi*eta^k - lambda^k = f_b(t^k)      (statico, contorno)
    %   M_ii*u_i_dot + C_ii*u_i + K_ii*eta + K_ib*qb = f_i(t)   (dinamico, interno)
    %
    %   0 <= (gap - qb) _|_ lambda >= 0      (contatto unilaterale, Signorini)
    %
    % Il contatto e' risolto con un active-set ESATTO per l'LCP accoppiato
    % (necessario perche' K_bb non e' diagonale quando ci sono piu' nodi di
    % contatto elasticamente accoppiati - cfr. Sez. 4.3 della tesi, e
    % Appendice C per il caso generale).
    
    properties
        M, K, C
        MasslessTol  % tolleranza relativa per il controllo massless-boundary
    end
    
    methods
        function obj = TransientSolverMasslessCB(M, K, C, varargin)
            p = inputParser;
            addParameter(p, 'MasslessTol', 1e-8);
            parse(p, varargin{:});
            
            obj.K = K;
            if nargin < 3 || isempty(C)
                obj.C = zeros(size(K));
            else
                obj.C = C;
            end
            obj.MasslessTol = p.Results.MasslessTol;
            obj.M = M;
        end
        
        function [t_out, q_out, lambda_out] = solve(obj, tmax, dt, q0, qd0, F_handle, varargin)
            % 1. Parsing parametri
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p, 'ContactTargetDOF', 1);
            addParameter(p, 'ContactGap', 0);
            addParameter(p, 'LCPTol', 1e-10);
            addParameter(p, 'LCPMaxIter', 100);
            addParameter(p, 'RecordEvery', 1);
            parse(p, varargin{:});
            
            b   = p.Results.ContactTargetDOF;
            gap = p.Results.ContactGap;
            lcp_tol     = p.Results.LCPTol;
            lcp_maxiter = p.Results.LCPMaxIter;
            rec_every   = max(1, round(p.Results.RecordEvery));
            
            all_dofs = 1:size(obj.K, 1);
            i = setdiff(all_dofs, b);
            nb = length(b);
            
            if isscalar(gap)
                gap = gap * ones(nb, 1);
            else
                gap = gap(:);
                if length(gap) ~= nb
                    error('TransientSolverMasslessCB:GapSize', ...
                        'ContactGap deve essere scalare oppure avere la stessa lunghezza di ContactTargetDOF.');
                end
            end
            
            % 2. Partizione delle matrici
            K_bb = full(obj.K(b, b));
            K_bi = obj.K(b, i);
            K_ib = obj.K(i, b);
            K_ii = obj.K(i, i);
            
            M_bb = obj.M(b, b);
            M_bi = obj.M(b, i);
            M_ii = obj.M(i, i);
            C_ii = obj.C(i, i);
            
            % --- Controllo dell'ipotesi di contorno privo di massa ---
            % Lo schema semi-esplicito richiede M_bb = 0 e M_bi = 0
            % ESATTAMENTE (Eq. 4.1 non contiene alcun termine di massa al
            % contorno). Per una base Massless-CB questi blocchi sono nulli
            % solo ANALITICAMENTE (Eq. 5.9); numericamente restano residui.
            % Se sono trascurabili li azzeriamo esplicitamente (coerente con
            % "neglecting the inertia forces...", Sez. 5.1.2); se sono
            % significativi, la base non e' massless e l'integratore non e'
            % applicabile: si segnala con errore.
            scale_M = max(norm(full(M_ii), 'fro'), eps);
            res_bb = norm(full(M_bb), 'fro');
            res_bi = norm(full(M_bi), 'fro');
            rel_res = max(res_bb, res_bi) / scale_M;
            
            if rel_res > obj.MasslessTol
                if rel_res > 1e-2
                    error('TransientSolverMasslessCB:NotMassless', ...
                        ['La matrice di massa ridotta NON e'' a contorno privo di massa ' ...
                         '(||M_bb||=%.3e, ||M_bi||=%.3e, rispetto a ||M_ii||=%.3e, residuo relativo %.3e). ' ...
                         'Verificare la costruzione della base CMS (RomMCB.build): il blocco M_bi ' ...
                         'dovrebbe annullarsi analiticamente per costruzione (Eq. 5.9 della tesi).'], ...
                        res_bb, res_bi, scale_M, rel_res);
                else
                    warning('TransientSolverMasslessCB:SmallMassResidual', ...
                        ['Residuo di massa al contorno non nullo ma piccolo (residuo relativo %.3e). ' ...
                         'Azzerato esplicitamente prima dell''integrazione, come previsto dal metodo ' ...
                         '(Sez. 5.1.2 della tesi).'], rel_res);
                end
            end
            % Azzeramento esplicito (sicuro entro la tolleranza verificata sopra)
            M_bb = zeros(nb, nb); %#ok<NASGU> % non serve piu', il blocco non e' usato nell'algoritmo
            
            % 3. Preallocazione OUTPUT (sottocampionata)
            % --- Nota memoria: l'integrazione avviene sempre al passo dt fine
            % (necessario per stabilita'), ma si SALVA in memoria solo un
            % campione ogni RecordEvery passi. Con dt molto piccolo (es. 1e-9)
            % e tmax non altrettanto piccolo, salvare OGNI passo produce un
            % numero di colonne enorme, che esplode ulteriormente quando
            % espanso ai DOF fisici completi (Pc*q_rom). Impostare RecordEvery
            % (es. 100-1000) per tenere la memoria sotto controllo mantenendo
            % la risoluzione temporale interna fine.
            n_steps_full = floor(tmax/dt) + 1;
            n_out = floor((n_steps_full - 1) / rec_every) + 1;
            
            qb  = zeros(nb, n_out);
            eta = zeros(length(i), n_out);
            lambda_hist = zeros(nb, n_out);
            t_rec = zeros(n_out, 1);
            
            % 4. Condizioni iniziali
            qb_curr  = q0(b);
            eta_curr = q0(i);
            qb(:, 1)  = qb_curr;
            eta(:, 1) = eta_curr;
            t_rec(1) = 0;
            rec_idx = 1;
            
            F_ext_0 = F_handle(0);
            F_i_0   = F_ext_0(i);
            accel_0 = M_ii \ (F_i_0 - C_ii * qd0(i) - K_ib * qb_curr - K_ii * eta_curr);
            u_eta_half = qd0(i) + 0.5 * dt * accel_0;
            
            % --- Operatori impliciti per lo smorzamento (Fig. 4.2a) ---
            % u_i^(k+1/2) = A_impl^-1 [ f_i - K_ii q_i - K_ib q_b + A_expl*u_i^(k-1/2) ]
            A_impl = M_ii / dt + C_ii / 2;
            A_expl = M_ii / dt - C_ii / 2;
            [L_A, U_A, P_A] = lu(A_impl);
            solveInner = @(rhs) U_A \ (L_A \ (P_A * rhs));
            
            % Warm-start dell'active-set in base alla condizione iniziale
            active_set = qb_curr > gap;
            
            % 5. Ciclo di integrazione (avanza sempre a dt fine, salva ogni rec_every)
            lambda_curr = zeros(nb, 1);
            for k = 1:(n_steps_full - 1)
                t = (k - 1) * dt;
                F_ext = F_handle(t);
                F_b = F_ext(b);
                F_i = F_ext(i);
                
                % --- Step algebrico: contatto (LCP accoppiato, esatto) ---
                rhs_contact = F_b - K_bi * eta_curr;
                [qb_next, lambda_k, active_set] = obj.solveContactLCP( ...
                    K_bb, rhs_contact, gap, active_set, lcp_maxiter, lcp_tol);
                
                % --- Step dinamico: coordinate interne (implicito su smorzamento) ---
                rhs_i = F_i - K_ii * eta_curr - K_ib * qb_next + A_expl * u_eta_half;
                u_eta_half = solveInner(rhs_i);
                eta_next = eta_curr + dt * u_eta_half;
                
                qb_curr  = qb_next;
                eta_curr = eta_next;
                lambda_curr = lambda_k;
                
                if mod(k, rec_every) == 0
                    rec_idx = rec_idx + 1;
                    qb(:, rec_idx)  = qb_curr;
                    eta(:, rec_idx) = eta_curr;
                    lambda_hist(:, rec_idx) = lambda_curr;
                    t_rec(rec_idx) = k * dt;
                end
            end
            if rec_idx < n_out
                rec_idx = n_out;
                qb(:, rec_idx)  = qb_curr;
                eta(:, rec_idx) = eta_curr;
                lambda_hist(:, rec_idx) = lambda_curr;
                t_rec(rec_idx) = (n_steps_full - 1) * dt;
            end
            
            % 6. Ricomposizione output (solo campioni salvati, non ogni passo)
            t_out = t_rec;
            q_out = zeros(length(all_dofs), n_out);
            q_out(b, :) = qb;
            q_out(i, :) = eta;
            lambda_out = lambda_hist;
        end
    end
    
    methods (Access = private, Static)
        function [qb, lambda, active] = solveContactLCP(Kbb, rhs, gap, active_guess, maxiter, tol)
            % Risolve il problema (contatto unilaterale, vincolo SUPERIORE qb<=gap):
            %   Kbb*qb + lambda = rhs      (lambda = reazione del muro, >=0 se attivo)
            %   qb <= gap ,  lambda >= 0 ,  (gap-qb)'*lambda = 0
            %
            % NOTA SUL SEGNO: con vincolo superiore (il nodo si muove VERSO il
            % muro, che lo respinge), la reazione fisica lambda si SOTTRAE alla
            % forza applicata: Kbb*qb = rhs - lambda, cioe' lambda = rhs - Kbb*qb.
            % (Il segno opposto, lambda = Kbb*qb - rhs, corrisponderebbe invece
            % a un vincolo INFERIORE qb>=gap come nell'esempio canonico della
            % barra che rimbalza della tesi, Sez. 7.1 - qui la fisica e' speculare).
            %
            % Pivoting a UN indice per volta (il piu' violato) ad ogni
            % iterazione: necessario per la garanzia di terminazione finita
            % del metodo (Kbb e le sue sottomatrici principali sono definite
            % positive per costruzione CMS => P-matrix => l'active-set converge
            % in un numero finito di passi SOLO se si flippa un indice alla
            % volta; flippando piu' violatori insieme si puo' ciclare
            % indefinitamente tra due active-set incompatibili).
            n = length(rhs);
            active = active_guess(:);
            qb = zeros(n, 1);
            lambda = zeros(n, 1);
            
            for iter = 1:maxiter
                A = active;
                I = ~active;
                
                qb(A) = gap(A);
                if any(I)
                    if any(A)
                        rhsI = rhs(I) - Kbb(I, A) * gap(A);
                    else
                        rhsI = rhs(I);
                    end
                    qb(I) = Kbb(I, I) \ rhsI;
                end
                
                lambda = rhs - Kbb * qb;
                
                viol_open  = I & (qb > gap + tol);   % penetra: deve chiudersi
                viol_close = A & (lambda < -tol);     % reazione negativa (trazione): deve aprirsi
                
                if ~any(viol_open) && ~any(viol_close)
                    return;
                end
                
                % Flip di un solo indice: il piu' violato tra i due gruppi
                open_mag  = -inf(n, 1); open_mag(viol_open)   = qb(viol_open) - gap(viol_open);
                close_mag = -inf(n, 1); close_mag(viol_close) = -lambda(viol_close);
                [mo, io] = max(open_mag);
                [mc, ic] = max(close_mag);
                if mo >= mc
                    active(io) = true;
                else
                    active(ic) = false;
                end
            end
            
            warning('TransientSolverMasslessCB:LCPNotConverged', ...
                'Active-set del contatto non convergente in %d iterazioni.', maxiter);
        end
    end
end