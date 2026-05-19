classdef TransientSimulatorFull < handle
    properties
        time_rom, qb_rom, lambda_rom
    end
    
    methods
        function obj = TransientSimulatorFull()
        end
        
        function solve_drop(obj, initial_height, gravity, tmax, dt, Mc, Kc, nl_dof_c)
            % IMPLEMENTAZIONE ESATTA: Discretization 7.2 (Doyen et al. 2011)
            % Metodo a Massa Modificata sul Modello FEM Completo (Semi-Esplicito)
            
            obj.time_rom = 0:dt:tmax;
            n_steps = length(obj.time_rom);
            n_dofs = size(Kc, 1);
            
            % Partizionamento Nodi Interni (*) e Nodo di Contatto (c)
            idx_c = nl_dof_c;
            idx_star = setdiff(1:n_dofs, idx_c);
            
            K_ss = Kc(idx_star, idx_star);
            K_sc = Kc(idx_star, idx_c);
            K_cs = Kc(idx_c, idx_star);
            K_cc = Kc(idx_c, idx_c);
            
            % La massa del nodo di contatto viene ignorata (Modified Mass)
            M_ss = Mc(idx_star, idx_star);
            
            % Vettore Forze (Gravità)
            g_vec = zeros(n_dofs, 1);
            g_vec(1:3:end) = gravity; % Applica g solo lungo l'asse longitudinale
            F_ext = Mc * g_vec;       % Calcolato con la massa originale
            
            F_s = F_ext(idx_star);
            F_c = 0; % Il nodo è massless, quindi il suo peso locale è 0
            
            % --- Inizializzazione Storico (Verlet) ---
            u_star = zeros(length(idx_star), n_steps);
            u_c    = zeros(1, n_steps);
            obj.lambda_rom = zeros(1, n_steps);
            
            u_star(:, 1) = initial_height;
            u_c(1)       = initial_height;
            
            % Calcolo accelerazione a t=0 (Caduta Libera)
            a_star_0 = M_ss \ (F_s - K_ss * u_star(:,1) - K_sc * u_c(1));
            
            % Step fittizio al tempo t = -dt
            u_star_prev = u_star(:,1) + 0.5 * dt^2 * a_star_0;
            
            % --- LOOP TEMPORALE (Eq. 7.16 - 7.20 del paper) ---
            for k = 1:(n_steps - 1)
                u_s_curr = u_star(:, k);
                u_c_curr = u_c(k);
                
                % 1. Step Esplicito per i nodi interni (Eq. 7.18)
                % Si usa la differenza centrale per l'accelerazione
                u_s_next = 2 * u_s_curr - u_star_prev + ...
                           dt^2 * (M_ss \ (F_s - K_ss * u_s_curr - K_sc * u_c_curr));
                
                % 2. Step Statico per il nodo di contatto (Eq. 7.19)
                % Il nodo massless si adatta istantaneamente al resto della struttura
                u_c_trial = (F_c - K_cs * u_s_next) / K_cc;
                
                % 3. Condizione LCP di Signorini (Eq. 7.20)
                if u_c_trial >= 0
                    % Volo: nessuna reazione
                    u_c_next = u_c_trial;
                    lambda_next = 0;
                else
                    % Impatto: penetrazione bloccata rigidamente a 0
                    u_c_next = 0;
                    lambda_next = K_cs * u_s_next - F_c;
                end
                
                % Salvataggio
                u_star(:, k+1) = u_s_next;
                u_c(k+1)       = u_c_next;
                obj.lambda_rom(k) = lambda_next;
                
                % Shift Memoria
                u_star_prev = u_s_curr;
            end
            
            obj.qb_rom = u_c;
            obj.lambda_rom(end) = 0;
        end
    end
end