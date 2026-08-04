classdef MasslessVerletSolver < handle
    properties
        M, K, C
    end
    
    methods
        function obj = MasslessVerletSolver(M, K, C)
            obj.M = M;
            obj.K = K;
            obj.C = C;
        end
        
        function [t, q_hist] = solve(obj, tmax, dt, q0_reduced, qd0_reduced, F_reduced_handle, n_bnd, gap_wall)
            % Implementazione Verlet Leapfrog per ROM Massless (Craig-Bampton / MacNeal)
            % Include solutore LCP Active-Set con Single Pivoting e convenzione di segno esatta
            
            t = 0:dt:tmax;
            Nt = length(t);
            n_total = size(obj.K, 1);
            
            % Allocazione della cronologia degli spostamenti
            q_hist = zeros(n_total, Nt);
            q_hist(:, 1) = q0_reduced;
            
            % Indici per il partizionamento interno/bordo
            b_idx = 1:n_bnd;
            i_idx = (n_bnd + 1):n_total;
            
            % Estrazione sottomatrici 
            K_bb = obj.K(b_idx, b_idx);
            K_bi = obj.K(b_idx, i_idx);
            K_ib = obj.K(i_idx, b_idx);
            K_ii = obj.K(i_idx, i_idx);
            M_ii = obj.M(i_idx, i_idx);
            D_ii = obj.C(i_idx, i_idx); 
            
            % Inizializzazione delle velocità Leapfrog
            F_init = F_reduced_handle(0);
            f_i_0 = F_init(i_idx);
            q_b_0 = q0_reduced(b_idx);
            q_i_0 = q0_reduced(i_idx);
            
            ud_i_0 = M_ii \ (f_i_0 - K_ii*q_i_0 - K_ib*q_b_0 - D_ii*qd0_reduced(i_idx));
            u_i_half = qd0_reduced(i_idx) + (dt/2) * ud_i_0;
            
            % Pre-calcolo matrici dinamiche
            LHS_vel = (1/dt * M_ii + 0.5 * D_ii);
            RHS_vel_mat = (1/dt * M_ii - 0.5 * D_ii);
            K_bb_diag = diag(K_bb); % Utile per scalare le violazioni
            
            % Inizializzazione guess per l'active set (Cold start)
            is_active = false(n_bnd, 1);
            
            % --- LOOP TEMPORALE ---
            for j = 1:(Nt-1)
                t_j = t(j);
                q_i_j = q_hist(i_idx, j);
                F_j = F_reduced_handle(t_j);
                f_b_j = F_j(b_idx);
                f_i_j = F_j(i_idx);
                
                % --- STEP 1: SOLUTORE LCP (ACTIVE-SET SINGLE PIVOTING) ---
                % rhs = Forza esterna - Trascinamento elastico interno
                rhs = f_b_j - K_bi * q_i_j;
                
                q_b_real = zeros(n_bnd, 1);
                lambda = zeros(n_bnd, 1);
                
                % Garanzia di terminazione finita: numero di iterazioni max limitato
                max_iter = 5 * n_bnd; 
                for iter = 1:max_iter
                    % Imposizione vincoli set attuale
                    q_b_real(is_active) = gap_wall;
                    lambda(~is_active) = 0;
                    
                    % Risoluzione gradi di libertà inattivi (volo libero)
                    if any(~is_active)
                        q_b_real(~is_active) = K_bb(~is_active, ~is_active) \ ...
                            (rhs(~is_active) - K_bb(~is_active, is_active) * q_b_real(is_active));
                    end
                    
                    % Risoluzione moltiplicatori set attivo (forze di contatto R = rhs - Kbb*qb)
                    if any(is_active)
                        lambda(is_active) = rhs(is_active) - K_bb(is_active, :) * q_b_real;
                    end
                    
                    % Controllo Condizioni KKT
                    w = gap_wall - q_b_real; 
                    viol_pen = w < -1e-12 & ~is_active;    % Compenetrazione (il muro è stato superato)
                    viol_ten = lambda < -1e-12 & is_active; % Il muro sta tirando (fisicamente impossibile)
                    
                    if ~any(viol_pen) && ~any(viol_ten)
                        break; % Convergenza LCP raggiunta
                    end
                    
                    % SINGLE PIVOTING: Individua il singolo vincolo più violato
                    max_v = -1;
                    idx_to_flip = -1;
                    
                    for k = 1:n_bnd
                        if viol_pen(k)
                            % Scaliamo l'errore cinematico per la rigidezza locale (pseudo-forza)
                            val = -w(k) * K_bb_diag(k); 
                            if val > max_v
                                max_v = val; idx_to_flip = k;
                            end
                        elseif viol_ten(k)
                            % L'errore di forza non richiede scaling
                            val = -lambda(k);
                            if val > max_v
                                max_v = val; idx_to_flip = k;
                            end
                        end
                    end
                    
                    % Flip dello stato di attivazione per il singolo nodo peggiore
                    if idx_to_flip > 0
                        is_active(idx_to_flip) = ~is_active(idx_to_flip);
                    else
                        break; 
                    end
                end
                
                if iter == max_iter && max_iter > 0
                    warning('LCP non convergente a t = %g. Possibile instabilità numerica.', t_j);
                end
                
                % Salvataggio del bordo reale calcolato
                q_hist(b_idx, j) = q_b_real;
                
                % --- STEP 2: AGGIORNAMENTO VELOCITÀ INTERNE ELETTRO-DINAMICHE ---
                RHS_v = f_i_j - K_ii * q_i_j - K_ib * q_b_real + RHS_vel_mat * u_i_half;
                u_i_next_half = LHS_vel \ RHS_v;
                
                % --- STEP 3: AGGIORNAMENTO SPOSTAMENTI INTERNI ---
                q_hist(i_idx, j+1) = q_i_j + u_i_next_half * dt;
                
                u_i_half = u_i_next_half;
            end
            
            % Risoluzione approssimata all'ultimo istante per la coerenza dei grafici
            F_end = F_reduced_handle(t(end));
            rhs_end = F_end(b_idx) - K_bi * q_hist(i_idx, end);
            q_hist(b_idx, end) = K_bb \ rhs_end; 
        end
    end
end