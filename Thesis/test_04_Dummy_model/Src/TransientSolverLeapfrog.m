classdef TransientSolverLeapfrog < handle
    properties
        M, K, C
    end
    
    methods
        function obj = TransientSolverLeapfrog(M, K, C)
            obj.M = M;
            obj.K = K;
            if nargin < 3 || isempty(C)
                obj.C = zeros(size(K));
            else
                obj.C = C;
            end
        end
        
        function [t_out, q_out, lambda_out] = solve(obj, tmax, dt, q0, qd0, F_handle, varargin)
            % 1. Parsing dei parametri opzionali
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
            
            % 2. Setup degli indici e partizione matriciale
            all_dofs = 1:size(obj.K, 1);
            i = setdiff(all_dofs, b);
            nb = length(b);
            
            if isscalar(gap)
                gap = gap * ones(nb, 1);
            else
                gap = gap(:);
                if length(gap) ~= nb
                    error('TransientSolverLeapfrog:GapSize', ...
                        'ContactGap must be scalar or have the same length as ContactTargetDOF.');
                end
            end
            
            K_bb = full(obj.K(b, b));
            K_bi = obj.K(b, i);
            K_ib = obj.K(i, b);
            K_ii = obj.K(i, i);
            
            M_ii = obj.M(i, i);
            C_ii = obj.C(i, i);
            
            % 3. Vettori di tempo (integrazione interna) e preallocazione OUTPUT
            % --- Nota memoria: l'integrazione avviene sempre al passo dt fine
            % (necessario per stabilita'), ma si SALVA in memoria solo un
            % campione ogni RecordEvery passi. Con dt molto piccolo (es. 1e-9)
            % e tmax non altrettanto piccolo, salvare OGNI passo produce un
            % numero di colonne enorme (n_steps ~ tmax/dt), che poi esplode
            % ulteriormente quando espanso ai DOF fisici completi (Pc*q_rom).
            % Impostare RecordEvery (es. 100-1000) per tenere la memoria sotto
            % controllo mantenendo la risoluzione temporale interna fine.
            n_steps_full = floor(tmax/dt) + 1;
            n_out = floor((n_steps_full - 1) / rec_every) + 1;
            
            qb  = zeros(nb, n_out);
            eta = zeros(length(i), n_out);
            lambda_hist = zeros(nb, n_out);
            t_rec = zeros(n_out, 1);
            
            % 4. Condizioni Iniziali
            qb_curr  = q0(b);
            eta_curr0 = q0(i);
            qb(:, 1)  = qb_curr;
            eta(:, 1) = eta_curr0;
            t_rec(1) = 0;
            rec_idx = 1;
            
            % Inizializzazione Leapfrog: calcolo accelerazione a t=0
            F_ext_0 = F_handle(0);
            F_i_0   = F_ext_0(i);
            accel_0 = M_ii \ (F_i_0 - C_ii * qd0(i) - K_ib * qb(:, 1) - K_ii * eta(:, 1));
            u_eta_half = qd0(i) + 0.5 * dt * accel_0;
            
            % --- Precompute implicit damping operators (constant in time) ---
            % Fig. 4.2:  u_i^(k+1/2) = A_impl^-1 [ f_i - K_ii q_i - K_ib q_b + A_expl*u_i^(k-1/2) ]
            % with A_impl = M_ii/dt + D_ii/2,  A_expl = M_ii/dt - D_ii/2
            A_impl = M_ii / dt + C_ii / 2;
            A_expl = M_ii / dt - C_ii / 2;
            [L_A, U_A, P_A] = lu(A_impl);
            solveInner = @(rhs) U_A \ (L_A \ (P_A * rhs));
            
            % Warm-start active set based on initial condition feasibility
            active_set = qb_curr > gap;
            
            % 5. Ciclo di integrazione Leapfrog Semi-Implicito
            % --- L'integrazione avanza sempre al passo fine dt (stabilita');
            % il salvataggio in qb/eta/lambda_hist avviene solo ogni
            % rec_every passi (controllo memoria).
            lambda_curr = zeros(nb, 1);
            for k = 1:(n_steps_full - 1)
                t = (k - 1) * dt;
                F_ext = F_handle(t);
                F_b = F_ext(b);
                F_i = F_ext(i);
                
                % --- Step Algebrico (Contatto) ---
                % Kbb*qb + lambda = rhs_contact ,  0<=(gap-qb) _|_ lambda>=0
                rhs_contact = F_b - K_bi * eta_curr0;
                [qb_next, lambda_k, active_set] = obj.solveContactLCP( ...
                    K_bb, rhs_contact, gap, active_set, lcp_maxiter, lcp_tol);
                
                % --- Step Esplicito/Implicito (Dinamica Interna) ---
                rhs_i = F_i - K_ii * eta_curr0 - K_ib * qb_next + A_expl * u_eta_half;
                u_eta_half = solveInner(rhs_i);
                eta_next = eta_curr0 + dt * u_eta_half;
                
                qb_curr  = qb_next;
                eta_curr0 = eta_next;
                lambda_curr = lambda_k;
                
                % Salvataggio ogni rec_every passi (k qui e' il passo appena concluso,
                % cioe' lo stato al tempo t_{k+1} = k*dt)
                if mod(k, rec_every) == 0
                    rec_idx = rec_idx + 1;
                    qb(:, rec_idx)  = qb_curr;
                    eta(:, rec_idx) = eta_curr0;
                    lambda_hist(:, rec_idx) = lambda_curr;
                    t_rec(rec_idx) = k * dt;
                end
            end
            % Se l'ultimo passo non e' caduto su un multiplo di rec_every,
            % ci si assicura comunque di avere l'ultimo stato disponibile.
            if rec_idx < n_out
                rec_idx = n_out;
                qb(:, rec_idx)  = qb_curr;
                eta(:, rec_idx) = eta_curr0;
                lambda_hist(:, rec_idx) = lambda_curr;
                t_rec(rec_idx) = (n_steps_full - 1) * dt;
            end
            
            % 6. Ricomposizione dell'output (solo campioni salvati, non ogni passo)
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
            % muro, che lo respinge), la reazione fisica si SOTTRAE alla forza
            % applicata: Kbb*qb = rhs - lambda, cioe' lambda = rhs - Kbb*qb.
            %
            % Pivoting a UN indice per volta (il piu' violato) ad ogni
            % iterazione: necessario per la garanzia di terminazione finita
            % (Kbb e le sue sottomatrici principali sono definite positive per
            % costruzione CMS => P-matrix => convergenza garantita SOLO se si
            % flippa un indice alla volta; flippando piu' violatori insieme si
            % puo' ciclare indefinitamente tra due active-set incompatibili -
            % cfr. Appendice C della tesi, generalizzato al caso accoppiato).
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
                viol_close = A & (lambda < -tol);     % reazione negativa: deve aprirsi
                
                if ~any(viol_open) && ~any(viol_close)
                    return;
                end
                
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
            
            warning('TransientSolverLeapfrog:LCPNotConverged', ...
                ['Contact active-set non convergente in %d iterazioni. ' ...
                 'Risultato potrebbe violare le condizioni KKT.'], maxiter);
        end
    end
end