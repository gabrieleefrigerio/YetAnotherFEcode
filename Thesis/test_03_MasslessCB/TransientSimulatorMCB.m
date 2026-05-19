classdef TransientSimulatorMCB < handle
    properties
        M_r, K_r, P, alpha 
        K_bb, K_bi, K_ib, K_ii
        time_rom, qb_rom, lambda_rom, eta_rom
    end
    methods
        function obj = TransientSimulatorMCB(rom_builder)
            obj.M_r = rom_builder.M_r;
            obj.K_r = rom_builder.K_r;
            obj.P   = rom_builder.P;
            if isprop(rom_builder, 'alpha')
                obj.alpha = rom_builder.alpha;
            end
            obj.K_bb = obj.K_r(1, 1);
            obj.K_bi = obj.K_r(1, 2:end);
            obj.K_ib = obj.K_r(2:end, 1);
            obj.K_ii = obj.K_r(2:end, 2:end);
        end
        
        function solve_drop(obj, initial_height, gravity, tmax, dt, M_full)
            obj.time_rom = 0:dt:tmax;
            n_steps = length(obj.time_rom);
            num_modes = size(obj.K_ii, 1);
            
            obj.qb_rom = zeros(1, n_steps);
            obj.lambda_rom = zeros(1, n_steps);
            obj.eta_rom = zeros(num_modes, n_steps);
            
            % External forces projected via the new decoupled basis
            n_dofs_full = size(obj.P, 1);
            g_vec = zeros(n_dofs_full, 1);
            g_vec(1:3:end) = gravity;
            F_ext = obj.P' * (M_full * g_vec);
            F_b = F_ext(1);
            F_i = F_ext(2:end);
            
            % Initial Conditions
            qb_curr = initial_height;
            
            % Correction: Initialize eta_curr compensating for the base deformation
            if ~isempty(obj.alpha)
                eta_curr = obj.alpha * qb_curr;
            else
                eta_curr = zeros(num_modes, 1); % Return to default for other ROMs
            end
            
            obj.qb_rom(1) = qb_curr;
            obj.eta_rom(:, 1) = eta_curr;
            
            % Initial acceleration of flexible modes (M_ii = I)
            qdd_eta = F_i - obj.K_ib * qb_curr - obj.K_ii * eta_curr;
            u_eta_half = -0.5 * dt * qdd_eta; % Half step backward for Leapfrog
            
            for k = 1:(n_steps - 1)
                % 1. Algebraic calculation of the contact coordinate (Free-flight trial)
                qb_trial = (F_b - obj.K_bi * eta_curr) / obj.K_bb;
                
                % 2. Signorini condition
                if qb_trial >= 0
                    qb_curr = qb_trial;
                    lambda_curr = 0;
                else
                    qb_curr = 0;
                    lambda_curr = obj.K_bi * eta_curr - F_b;
                end
                
                % 3. Explicit dynamics of internal modes (pure ODE)
                qdd_eta = F_i - obj.K_ib * qb_curr - obj.K_ii * eta_curr;
                
                % 4. Verlet Leapfrog update
                u_eta_half = u_eta_half + dt * qdd_eta;
                eta_next = eta_curr + dt * u_eta_half;
                
                % Save data
                obj.qb_rom(k) = qb_curr;
                obj.lambda_rom(k) = lambda_curr;
                obj.qb_rom(k+1) = qb_curr; % Graphical approximation for the future step
                obj.eta_rom(:, k+1) = eta_next;
                
                eta_curr = eta_next;
            end
            obj.lambda_rom(end) = 0;
        end
    end
end