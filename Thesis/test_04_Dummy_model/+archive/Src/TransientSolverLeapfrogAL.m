classdef TransientSolverLeapfrogAL < handle
    properties
        M, K, C
    end
    
    methods
        function obj = TransientSolverLeapfrogAL(M, K, C)
            obj.M = M;
            obj.K = K;
            if nargin < 3 || isempty(C)
                obj.C = zeros(size(K));
            else
                obj.C = C;
            end
        end
        
        function [t_out, q_out] = solve(obj, tmax, dt, q0, qd0, F_handle, varargin)
            % Parsing dei parametri
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p, 'ContactTargetDOF', 1);
            addParameter(p, 'ContactGap', 0);
            % Epsilon AL: Parametro di tuning cruciale del paper (es. 0.5 - 1.5)
            addParameter(p, 'EpsilonAL', 0.8); 
            addParameter(p, 'MaxIter', 150);
            addParameter(p, 'Tol', 1e-8);
            parse(p, varargin{:});
            
            b = p.Results.ContactTargetDOF;
            gap = p.Results.ContactGap;
            eps_AL = p.Results.EpsilonAL;
            max_iter = p.Results.MaxIter;
            tol = p.Results.Tol;
            
            % Setup indici partizione
            all_dofs = 1:size(obj.K, 1);
            i = setdiff(all_dofs, b);
            nb = length(b);
            
            K_bb = obj.K(b, b);
            K_bi = obj.K(b, i);
            K_ib = obj.K(i, b);
            K_ii = obj.K(i, i);
            M_ii = obj.M(i, i);
            C_ii = obj.C(i, i);
            
            % PRE-CALCOLO OFFLINE DELLA MATRICE DI DELASSUS STATICA (G)
            % Poiché K_bb è piccola (solo GdL di contatto X), l'inversione è istantanea
            G = inv(full(K_bb)); 
            
            time = 0:dt:tmax;
            n_steps = length(time);
            
            qb = zeros(nb, n_steps);
            eta = zeros(length(i), n_steps);
            lambda_hist = zeros(nb, n_steps);
            
            % Condizioni Iniziali
            qb(:, 1) = q0(b);
            eta(:, 1) = q0(i);
            
            % Inizializzazione Leapfrog (Accelerazione Modi Interni)
            F_ext_0 = F_handle(0);
            F_i_0 = F_ext_0(i);
            accel_0 = M_ii \ (F_i_0 - C_ii * qd0(i) - K_ib * qb(:, 1) - K_ii * eta(:, 1));
            u_eta_half = qd0(i) + 0.5 * dt * accel_0; 
            
            % Matrici Implicite per Smorzamento (Flowchart Fig 4.2 del paper)
            A_impl = M_ii / dt + C_ii / 2;
            A_expl = M_ii / dt - C_ii / 2;
            [L_A, U_A, P_A] = lu(A_impl);
            solveInner = @(rhs) U_A \ (L_A \ (P_A * rhs));

            for k = 1:(n_steps - 1)
                t = time(k);
                F_ext = F_handle(t);
                F_b = F_ext(b);
                
                eta_curr = eta(:, k);
                
                % 1. Calcolo del termine di Volo Libero
                rhs_b = F_b - K_bi * eta_curr;
                q_free = G * rhs_b; 
                
                % 2. AUGMENTED LAGRANGIAN (Jacobi Proiettato)
                lambda = zeros(nb, 1);
                if k > 1
                   lambda = lambda_hist(:, k-1); % Warm start per convergenza rapida
                end
                
                for iter = 1:max_iter
                    % Calcolo della compenetrazione corrente
                    % w >= 0 significa che non stiamo penetrando il muro
                    w = G * lambda - q_free + gap; 
                    
                    % Proiezione
                    lambda_new = max(0, lambda - eps_AL * w);
                    
                    % Controllo tolleranza
                    if norm(lambda_new - lambda) < tol
                        lambda = lambda_new;
                        break;
                    end
                    lambda = lambda_new;
                end
                
                % Se iter arriva a max_iter, significa che eps_AL va tarato
                if iter == max_iter && k == 2 
                    warning('LCP algebrico lento a convergere. Considera di tarare EpsilonAL.');
                end
                
                % 3. Calcolo coordinata di contatto effettiva
                qb_next = q_free - G * lambda;
                
                % 4. STEP SEMI-IMPLICITO (Dinamica Interna)
                F_i = F_ext(i);
                rhs_i = F_i - K_ii * eta_curr - K_ib * qb_next + A_expl * u_eta_half;
                u_eta_half = solveInner(rhs_i);
                eta_next = eta_curr + dt * u_eta_half;
                
                % Salvataggio
                qb(:, k+1) = qb_next;
                eta(:, k+1) = eta_next;
                lambda_hist(:, k) = lambda;
            end
            
            % Ricomposizione formattata [Step temporali x GdL]
            t_out = time(:);
            q_out = zeros(size(obj.K, 1), n_steps);
            q_out(b, :) = qb;
            q_out(i, :) = eta;
            q_out = q_out'; 
        end
    end
end