classdef TransientSimulatorMoreau < handle
    properties
        M_r, K_r, P
        time_rom, qb_rom, lambda_rom, eta_rom
    end
    
    methods
        function obj = TransientSimulatorMoreau(rom_builder)
            obj.M_r = rom_builder.M_r;
            obj.K_r = rom_builder.K_r;
            obj.P   = rom_builder.P;
        end
        
        function solve_drop(obj, initial_height, gravity, tmax, dt, M_full)
            % IMPLEMENTAZIONE: Symmetric Moreau-like time integrator
            % Appendix B - Monjaraz Tec Thesis (Mass-carrying boundary)
            
            obj.time_rom = 0:dt:tmax;
            n_steps = length(obj.time_rom);
            n_dofs_rom = size(obj.M_r, 1);
            
            % Controllo di Sicurezza: Il modello DEVE avere massa al bordo
            if abs(obj.M_r(1,1)) < 1e-10
                error('Errore: M_r(1,1) è zero. Usa MacNeal o CB Standard, NON il Massless CB!');
            end
            
            obj.qb_rom = zeros(1, n_steps);
            obj.lambda_rom = zeros(1, n_steps);
            obj.eta_rom = zeros(n_dofs_rom - 1, n_steps);
            
            % 1. Setup Vettore Forze (Gravità)
            n_dofs_full = size(obj.P, 1);
            g_vec = zeros(n_dofs_full, 1);
            g_vec(1:3:end) = gravity;
            f_tk = obj.P' * (M_full * g_vec); % f(t^k)
            
            % Matrici fisse per l'algoritmo (Assumendo smorzamento D = 0)
            M_inv = inv(obj.M_r);
            
            % Eq. (B.4): G = W^T * M^-1 * W
            % Poiché W = [1; 0; ...; 0], G è semplicemente il primo elemento dell'inversa
            G = M_inv(1, 1);
            
            % --- STEP 1: Inizializzazione ---
            q_curr = zeros(n_dofs_rom, 1);
            q_curr(1) = initial_height;   % q^1
            u_prev = zeros(n_dofs_rom, 1); % u^{1/2}
            
            obj.qb_rom(1) = q_curr(1);
            
            % Coefficiente di restituzione normale (0 per Signorini/Doyen)
            e_rest = 0.0; 
            
            % --- STEP 4: Loop Temporale ---
            for k = 1:(n_steps - 1)
                
                % Integrazione delle forze nel passo dt (Impulsi)
                f_dt = dt * f_tk;
                Kq_dt = dt * (obj.K_r * q_curr);
                
                % --- STEP 2: Calcolo LCP (Eq. B.1, B.2, B.3) ---
                
                % Calcolo la velocità "libera" (senza contatto)
                % NOTA: Qui correggo il refuso di (B.5). M*u_prev ha il segno '+'
                u_free = M_inv * (f_dt - Kq_dt + obj.M_r * u_prev);
                
                % Eq. (B.5) Calcolo di 'c'
                c = e_rest * u_prev(1) + u_free(1);
                
                % Valutazione Set Attivo ("Contacts having non-positive normal gap")
                if q_curr(1) <= 0
                    % Contatto attivo: risolvo l'inclusione algebrica (Eq B.3)
                    DeltaP = max(0, -c / G);
                else
                    % Volo libero
                    DeltaP = 0;
                end
                
                % Aggiornamento Velocità u^{k+1/2} (Eq B.1)
                W_P = zeros(n_dofs_rom, 1);
                W_P(1) = DeltaP;
                
                u_next = u_free + M_inv * W_P;
                
                % --- STEP 3: Aggiornamento Posizioni ---
                % Eq: q^{k+1} = q^k + u^{k+1/2} * dt
                q_next = q_curr + dt * u_next;
                
                % Salvataggio dati per il plot
                obj.qb_rom(k+1) = q_next(1);
                obj.eta_rom(:, k+1) = q_next(2:end);
                obj.lambda_rom(k) = DeltaP; % Percussione al contatto
                
                % Avanzamento passo
                q_curr = q_next;
                u_prev = u_next;
            end
        end
    end
end