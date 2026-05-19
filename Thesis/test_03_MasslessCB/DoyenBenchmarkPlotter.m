classdef DoyenBenchmarkPlotter < handle
    % DOYENBENCHMARKPLOTTER Computes the exact analytical solution and compares it with Massless CB
    
    properties
        sim_cb          % Simulator with Massless Craig-Bampton ROM
        h0              % Initial height
        g0              % Gravity acceleration (magnitude)
        % Fixed physical parameters of the Doyen benchmark
        L = 10;
        E = 900;
        rho = 1;
        Area = 1;       
    end
    
    methods
        function obj = DoyenBenchmarkPlotter(sim_cb, h0, gravity)
            obj.sim_cb = sim_cb;
            obj.h0 = h0;
            obj.g0 = abs(gravity);
        end
        
        function plot_kinematics(obj)
            % Time vector
            time_ex = obj.sim_cb.time_rom; 
            u_ex = zeros(size(time_ex));
            
            % Derived parameters
            c0 = sqrt(obj.E / obj.rho);
            tau_w = obj.L / c0;           
            tau_f = sqrt(2 * obj.h0 / obj.g0);  
            vf = obj.g0 * tau_f;        
            T_cycle = 16 * tau_w;     
            
            % Exact Analytical Trajectory
            for i = 1:length(time_ex)
                t = time_ex(i);
                if t <= tau_f
                    u_ex(i) = obj.h0 - 0.5 * obj.g0 * t^2;
                else
                    t_local = mod(t - tau_f, T_cycle);
                    if t_local <= 2*tau_w
                        u_ex(i) = 0;
                    elseif t_local > 2*tau_w && t_local <= 8*tau_w
                        dt_flight = t_local - 2*tau_w;
                        P = obj.h0 - 0.5 * obj.g0 * (dt_flight - tau_f)^2; 
                        S2 = -(2 * obj.g0 * obj.L^2)/(3 * c0^2);
                        for n = 1:50 
                            lam_n = n * pi / obj.L;
                            bn = (4 * obj.g0) / (c0^2 * lam_n^2);
                            S2 = S2 + bn * cos(c0 * lam_n * dt_flight);
                        end
                        u_ex(i) = P + S2;
                    elseif t_local > 8*tau_w && t_local <= 10*tau_w
                        u_ex(i) = 0;
                    else
                        dt_flight = t_local - 10*tau_w;
                        u_ex(i) = vf * dt_flight - 0.5 * obj.g0 * dt_flight^2;
                    end
                end
            end
            
            % Normalization
            max_u = max(u_ex);
            u_ex_norm = u_ex / max_u;
            qb_cb_norm = obj.sim_cb.qb_rom / max_u;
            
            figure('Name', 'ROM vs Analytical (Kinematics)', 'Color', 'w', 'Position', [150 150 900 500]);
            plot(time_ex, u_ex_norm, 'k:', 'LineWidth', 2.5, 'DisplayName', 'Analytical Solution');
            hold on;
            plot(obj.sim_cb.time_rom, qb_cb_norm, 'r-', 'LineWidth', 1.0, 'DisplayName', 'Massless Craig-Bampton');
            
            yline(0, 'k-', 'HandleVisibility', 'off');
            xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 12);
            ylabel('Normalized displacement', 'Interpreter', 'latex', 'FontSize', 12);
            title('ROM vs Analytical (Kinematics)', 'Interpreter', 'latex', 'FontSize', 14);
            legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 11);
            xlim([0 max(time_ex)]); ylim([-0.05 1.1]); 
            grid on;
        end
        
        function plot_energies(obj)
            % --- ENERGY PARAMETERS ---
            M_tot = obj.L * obj.Area * obj.rho;
            E0 = M_tot * obj.g0 * obj.h0; 
            
            c0 = sqrt(obj.E / obj.rho);
            tau_w = obj.L / c0;           
            tau_f = sqrt(2 * obj.h0 / obj.g0); 
            T_cycle = 16 * tau_w;
            
            time_ex = obj.sim_cb.time_rom;
            periods = time_ex / T_cycle;
            
            % 1. EXACT ANALYTICAL ENERGY CALCULATION
            E_el_ex = zeros(size(time_ex));
            E_tot_ex = ones(size(time_ex)); % Normalized to 1
            
            for i = 1:length(time_ex)
                t = time_ex(i);
                if t <= tau_f
                    E_el_ex(i) = 0; % Initial undeformed flight
                else
                    t_local = mod(t - tau_f, T_cycle);
                    
                    if t_local <= 2*tau_w
                        % Contact 1: Shock wave (grows and decreases linearly)
                        E_res = (2 * obj.rho * obj.g0^2 * obj.L^3) / (3 * c0^2);
                        if t_local <= tau_w
                            E_el_ex(i) = (E0 / tau_w) * t_local; 
                        else
                            E_el_ex(i) = E0 - ((E0 - E_res) / tau_w) * (t_local - tau_w);
                        end
                        
                    elseif t_local > 2*tau_w && t_local <= 8*tau_w
                        % Flight 1: Analytical modal oscillation of elastic energy
                        dt_f = t_local - 2*tau_w;
                        E_el = 0;
                        for n = 1:50
                            lam_n = n * pi / obj.L;
                            % Exact integral of strain for the n-th mode
                            coeff = (4 * obj.rho * obj.L * obj.g0^2) / (c0^2 * lam_n^2);
                            E_el = E_el + coeff * cos(c0 * lam_n * dt_f)^2;
                        end
                        E_el_ex(i) = E_el;
                        
                    elseif t_local > 8*tau_w && t_local <= 10*tau_w
                        % Contact 2: Secondary bounce (absorbs E_res, goes to max, returns to 0)
                        dt_c = t_local - 8*tau_w;
                        E_res = (2 * obj.rho * obj.g0^2 * obj.L^3) / (3 * c0^2);
                        if dt_c <= tau_w
                            E_el_ex(i) = E_res + ((E0 - E_res) / tau_w) * dt_c;
                        else
                            E_el_ex(i) = E0 - (E0 / tau_w) * (dt_c - tau_w);
                        end
                    else
                        % Flight 2: Perfectly rigid and undeformed bounce
                        E_el_ex(i) = 0;
                    end
                end
            end
            
            E_el_norm_ex = E_el_ex / E0;
            E_rb_norm_ex = E_tot_ex - E_el_norm_ex;
            
            % 2. ENERGIES CALCULATION FOR MASSLESS CRAIG-BAMPTON
            q_cb = [obj.sim_cb.qb_rom; obj.sim_cb.eta_rom];
            E_el_cb = 0.5 * sum(q_cb .* (obj.sim_cb.K_r * q_cb), 1);
            
            E_el_norm_cb = E_el_cb / E0;
            E_tot_norm_cb = ones(size(E_el_norm_cb));
            E_rb_norm_cb = E_tot_norm_cb - E_el_norm_cb;
            
            % --- PLOT CREATION ---
            figure('Name', 'Energy Comparison', 'Color', 'w', 'Position', [100, 100, 900, 500]);
            
            c_el = '#d62728';   % Red
            c_rb = '#ff7f0e';   % Orange
            c_tot = '#98df8a';  % Light green
            
            % a) Analytical Reference
            subplot(2,1,1);
            hold on; box on; grid on;
            plot(periods, E_rb_norm_ex, 'Color', c_rb, 'LineWidth', 0.8, 'DisplayName', '$E_{\mathrm{rb}}$');
            plot(periods, E_el_norm_ex, 'Color', c_el, 'LineWidth', 0.8, 'DisplayName', '$E_{\mathrm{el}}$');
            plot(periods, E_tot_ex, 'Color', c_tot, 'LineWidth', 1.5, 'DisplayName', '$E_{\mathrm{tot}}$');
            title('a) analytical reference', 'Interpreter', 'latex', 'FontSize', 13);
            ylabel('normalized energy', 'Interpreter', 'latex', 'FontSize', 12);
            xlabel('number of periods', 'Interpreter', 'latex', 'FontSize', 12);
            xlim([0 max(periods)]); ylim([0 1.1]);
            set(gca, 'TickLabelInterpreter', 'latex');
            legend({'$E_{\mathrm{el}}$', '$E_{\mathrm{rb}}$', '$E_{\mathrm{tot}}$'}, ...
                   'Location', 'east', 'Interpreter', 'latex', 'FontSize', 11);
               
            % b) Massless Craig-Bampton
            subplot(2,1,2);
            hold on; box on; grid on;
            plot(periods, E_rb_norm_cb, 'Color', c_rb, 'LineWidth', 0.8);
            plot(periods, E_el_norm_cb, 'Color', c_el, 'LineWidth', 0.8);
            plot(periods, E_tot_norm_cb, 'Color', c_tot, 'LineWidth', 1.5);
            title('b) massless', 'Interpreter', 'latex', 'FontSize', 13);
            ylabel('normalized energy', 'Interpreter', 'latex', 'FontSize', 12);
            xlabel('number of periods', 'Interpreter', 'latex', 'FontSize', 12);
            xlim([0 max(periods)]); ylim([0 1.1]);
            set(gca, 'TickLabelInterpreter', 'latex');
        end
    end
end