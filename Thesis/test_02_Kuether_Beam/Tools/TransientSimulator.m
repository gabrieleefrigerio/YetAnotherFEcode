classdef TransientSimulator < handle
    properties
        Assembly        % Full YAFEC Assembly (per FOM)
        ReducedAssembly % Reduced YAFEC Assembly (per ROM)
        spring          % LocalNonlinearity object
        P               % Reduction basis matrix
        M_full          % Full mass matrix
        nl_dof          % Nonlinear DOF index (physical)
        
        % Risultati ROM
        time_rom
        disp_center_rom
        energy_rom
        max_disp_rom
        rms_disp_rom
        
        % Risultati FOM
        time_fom
        disp_center_fom
        energy_fom
        max_disp_fom
        rms_disp_fom

        global_rel_error
    end
    
    methods
        function obj = TransientSimulator(full_assembly, reduced_assembly, spring_obj, P_matrix, M_full, nl_dof)
            obj.Assembly = full_assembly;
            obj.ReducedAssembly = reduced_assembly;
            obj.spring = spring_obj;
            obj.P = P_matrix;
            obj.M_full = M_full;
            obj.nl_dof = nl_dof;
        end
        
        %% --- SOLVER ROM (REDUCED ORDER MODEL) ---
        function solve_rom(obj, A_mag, duration_pulse, tmax, h)
            % 1. Proiezione spaziale della forza sul ROM
            phi_1 = obj.P(:, 1);
            spatial_force_red = obj.P' * (obj.M_full * phi_1 * A_mag);
            F_ext_red = @(t) spatial_force_red * (sin(pi * t / duration_pulse) * (t <= duration_pulse));
            
            % 2. Condizioni iniziali ridotte
            m_red = size(obj.P, 2);
            q0 = zeros(m_red, 1);
            qd0 = zeros(m_red, 1);
            qdd0 = zeros(m_red, 1);
            
            % 3. Integratore temporale YAFEC
            TI_NL = ImplicitNewmark('timestep', h, 'alpha', 0.0);
            
            % 4. Definizione del residuo non lineare del ROM
            Residual_NL_red = @(q, qd, qdd, t) residual_rom_contact(q, qd, qdd, t, ...
                                               obj.ReducedAssembly, obj.spring, obj.P, F_ext_red);
            
            % 5. Integrazione
            fprintf('Starting ROM integration (A = %.2f N)...\n', A_mag);
            TI_NL.Integrate(q0, qd0, qdd0, tmax, Residual_NL_red);
            
            % 6. Post-processing ROM
            obj.time_rom = TI_NL.Solution.time;
            q_history = TI_NL.Solution.q;
            u_history_full = obj.P * q_history;
            obj.disp_center_rom = u_history_full(obj.nl_dof, :);
            
            % Calcolo energia finale ROM
            q_end = q_history(:, end);
            q_penultimate = q_history(:, end-1);
            dt = obj.time_rom(end) - obj.time_rom(end-1);
            qd_end = (q_end - q_penultimate) / dt;
            
            Mr = obj.ReducedAssembly.DATA.M;
            Kr = obj.ReducedAssembly.DATA.K;
            
            E_kin = 0.5 * qd_end' * Mr * qd_end;
            E_pot_beam = 0.5 * q_end' * Kr * q_end;
            E_pot_spring = 0;
            if obj.disp_center_rom(end) > obj.spring.a
                E_pot_spring = 0.5 * obj.spring.k * (obj.disp_center_rom(end) - obj.spring.a)^2;
            end
            obj.energy_rom = E_kin + E_pot_beam + E_pot_spring;
            
            % Metriche spostamento post-impulso
            post_pulse_idx = obj.time_rom > duration_pulse;
            disp_free = obj.disp_center_rom(post_pulse_idx);
            obj.max_disp_rom = max(abs(disp_free));
            obj.rms_disp_rom = sqrt(mean(disp_free.^2));
        end
        
        %% --- SOLVER FOM (FULL ORDER MODEL) ---
        function solve_fom(obj, A_mag, duration_pulse, tmax, h)
            % 1. Forza spaziale nel dominio completo (FOM)
            M_full_c = obj.Assembly.constrain_matrix(obj.M_full);
            K_full_c = obj.Assembly.constrain_matrix(obj.Assembly.stiffness_matrix());
            [Phi_fom_all, ~] = eigs(K_full_c, M_full_c, 1, 'smallestabs');
            Phi_fom_1 = obj.Assembly.unconstrain_vector(Phi_fom_all(:, 1));
            
            % Normalizzazione rispetto alla massa completa per coerenza
            Phi_fom_1 = Phi_fom_1 / sqrt(Phi_fom_1' * obj.M_full * Phi_fom_1);
            
            spatial_force_fom = obj.M_full * Phi_fom_1 * A_mag;
            F_ext_fom = @(t) spatial_force_fom * (sin(pi * t / duration_pulse) * (t <= duration_pulse));
            
            % 2. Condizioni iniziali nel dominio completo (vincolato)
            n_dofs_constrained = size(M_full_c, 1);
            u0 = zeros(n_dofs_constrained, 1);
            ud0 = zeros(n_dofs_constrained, 1);
            udd0 = zeros(n_dofs_constrained, 1);
            
            % 3. Integratore temporale YAFEC
            TI_NL = ImplicitNewmark('timestep', h, 'alpha', 0.0);
            
            % Calcoliamo M_c e K_c UNA SOLA VOLTA fuori dal ciclo!
            M_full_c = obj.Assembly.constrain_matrix(obj.M_full);
            K_full_c = obj.Assembly.constrain_matrix(obj.Assembly.stiffness_matrix());
            
            % 4. Definizione del residuo FOM lineare
            F_ext_fom_constrained = @(t) obj.Assembly.constrain_vector(F_ext_fom(t));
            
            Residual_NL_fom = @(u, ud, udd, t) residual_fom_contact(u, ud, udd, t, ...
                                               F_ext_fom_constrained, M_full_c, K_full_c, obj.Assembly, obj.spring);
            
            % 5. Integrazione FOM
            fprintf('Starting FOM integration (A = %.2f N)...\n', A_mag);
            TI_NL.Integrate(u0, ud0, udd0, tmax, Residual_NL_fom);
            
            % 6. Post-processing FOM
            obj.time_fom = TI_NL.Solution.time;
            u_constrained_history = TI_NL.Solution.q;
            u_full_history = obj.Assembly.unconstrain_vector(u_constrained_history);
            obj.disp_center_fom = u_full_history(obj.nl_dof, :);
            
            % Calcolo energia finale FOM
            u_end_c = u_constrained_history(:, end);
            u_penultimate_c = u_constrained_history(:, end-1);
            dt = obj.time_fom(end) - obj.time_fom(end-1);
            ud_end_c = (u_end_c - u_penultimate_c) / dt;
            
            E_kin = 0.5 * ud_end_c' * M_full_c * ud_end_c;
            E_pot_beam = 0.5 * u_end_c' * K_full_c * u_end_c;
            E_pot_spring = 0;
            if obj.disp_center_fom(end) > obj.spring.a
                E_pot_spring = 0.5 * obj.spring.k * (obj.disp_center_fom(end) - obj.spring.a)^2;
            end
            obj.energy_fom = E_kin + E_pot_beam + E_pot_spring;
            
            % Metriche spostamento post-impulso
            post_pulse_idx = obj.time_fom > duration_pulse;
            disp_free = obj.disp_center_fom(post_pulse_idx);
            obj.max_disp_fom = max(abs(disp_free));
            obj.rms_disp_fom = sqrt(mean(disp_free.^2));
        end
        
        
        %% --- UTILITY: FORZA NON LINEARE PER FOM ---
        function [Kt, f_nl] = f_nonlinear_fom(obj, u_constrained)
            u_full = obj.Assembly.unconstrain_vector(u_constrained);
            u_nl = u_full(obj.spring.dof_idx);
            
            f_full = zeros(size(u_full, 1), 1);
            K_full = sparse(size(u_full, 1), size(u_full, 1));
            
            g = u_nl - obj.spring.a;
            if g > 0
                f_val = obj.spring.k * g;
                k_val = obj.spring.k;
            else
                f_val = 0;
                k_val = 0;
            end
            
            f_full(obj.spring.dof_idx) = f_val;
            K_full(obj.spring.dof_idx, obj.spring.dof_idx) = k_val;
            
            f_nl = obj.Assembly.constrain_vector(f_full);
            Kt = obj.Assembly.constrain_matrix(K_full);
        end

        %% --- POST-PROCESSING E VISUALIZZAZIONE ---
        
        function compute_error(obj)
            % Calcola l'errore relativo globale (GRE L2 Norm) percentuale
            if isempty(obj.disp_center_fom) || isempty(obj.disp_center_rom)
                error('Devi prima eseguire solve_fom() e solve_rom() per calcolare l''errore.');
            end
            
            diff_disp = obj.disp_center_fom - obj.disp_center_rom;
            norm_fom_disp = norm(obj.disp_center_fom);
            
            % Protezione per divisione per zero
            if norm_fom_disp == 0
                obj.global_rel_error = 0;
            else
                obj.global_rel_error = (norm(diff_disp) / norm_fom_disp) * 100;
            end
        end
        
        function print_summary(obj)
            % Stampa la tabella comparativa sulla Command Window
            if isempty(obj.global_rel_error)
                obj.compute_error();
            end
            
            fprintf('\n=======================================================\n');
            fprintf('             TRANSIENT COMPARISON (FOM vs ROM)         \n');
            fprintf('=======================================================\n');
            fprintf('Metric                 FOM               ROM           \n');
            fprintf('-------------------------------------------------------\n');
            fprintf('Total Energy [J]       %.6e      %.6e\n', obj.energy_fom, obj.energy_rom);
            fprintf('Max Displacement [mm]  %.6f          %.6f\n', obj.max_disp_fom*1000, obj.max_disp_rom*1000);
            fprintf('RMS Displacement [mm]  %.6f          %.6f\n', obj.rms_disp_fom*1000, obj.rms_disp_rom*1000);
            fprintf('-------------------------------------------------------\n');
            fprintf('Global Relative Error (L2 Norm)        %.4f %%\n', obj.global_rel_error);
            fprintf('=======================================================\n');
        end
        
        function plot_comparison(obj, A_mag, tmax)
            % Genera il plot sovrapposto con il riquadro delle metriche
            if isempty(obj.global_rel_error)
                obj.compute_error();
            end
            
            figure('Name', 'Transient Response Comparison: FOM vs ROM', 'Color', 'w');
            hold on;
            grid on;
            
            % Risposta FOM e ROM
            plot(obj.time_fom, obj.disp_center_fom * 1000, 'b-', 'LineWidth', 2);
            plot(obj.time_rom, obj.disp_center_rom * 1000, 'r--', 'LineWidth', 2);
            
            % Linea del Gap (legge direttamente le proprietà interne dell'oggetto)
            yline(obj.spring.a * 1000, 'k-.', 'Gap Clearance (a)', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
            
            % Estetica
            title('Transient Response at Midpoint', 'FontSize', 12);
            xlabel('Time [s]');
            ylabel('Midpoint Displacement [mm]');
            legend({'FOM', 'ROM [1 MC]'}, 'Location', 'best');
            xlim([0, tmax]);
            
            % Box testuale
            info_text = {
                sprintf('\\bfGRE (L_2):\\rm %.4f%%', obj.global_rel_error), ...
                sprintf('\\bfEnergy:\\rm %.6e J', obj.energy_fom), ...
                sprintf('\\bfAmplitude A:\\rm %.2f', A_mag),
            };
            
            text(0.97, 0.05, info_text, ...
                'Units', 'normalized', ...
                'VerticalAlignment', 'bottom', ...
                'HorizontalAlignment', 'right', ...
                'BackgroundColor', [0.97 0.97 0.97], ... 
                'EdgeColor', [0.6 0.6 0.6], ...          
                'LineWidth', 1, ...
                'FontSize', 9.5, ...
                'Margin', 8);
            hold off;
        end

        function save_results(obj, A_mag)
            % SAVE_RESULTS Salva le time-histories in un file .mat leggero.
            % Il nome del file viene generato in base all'ampiezza dell'impulso.
            
            if isempty(obj.time_fom) || isempty(obj.time_rom)
                error('Nessun dato da salvare. Esegui prima le simulazioni solve_fom e solve_rom.');
            end
            
            % Generazione del nome file dinamico
            % Sostituiamo il punto decimale con un underscore per non corrompere l'estensione del file
            A_str = strrep(sprintf('%.2f', A_mag), '.', '_');
            tmax_val = obj.time_fom(end);
            tmax_str = strrep(sprintf('%.2f', tmax_val), '.', '_');
            
            filename = sprintf('TransientResults_[1_3_5_7_9_11_MC]_Amag%s_Tmax%s.mat', A_str, tmax_str);
            
            % Estrazione dei dati puri (evitiamo di salvare le enormi matrici del FOM)
            time_fom = obj.time_fom;
            disp_center_fom = obj.disp_center_fom;
            time_rom = obj.time_rom;
            disp_center_rom = obj.disp_center_rom;
            energy_fom = obj.energy_fom;
            energy_rom = obj.energy_rom;
            
            % Assicuriamoci che l'errore sia calcolato prima di salvare
            if isempty(obj.global_rel_error)
                obj.compute_error();
            end
            global_rel_error = obj.global_rel_error;
            
            fprintf('Salvataggio dati in corso...\n');
            save(filename, 'time_fom', 'disp_center_fom', 'time_rom', 'disp_center_rom', ...
                           'energy_fom', 'energy_rom', 'global_rel_error', 'A_mag');
            fprintf('Salvataggio completato: %s\n', filename);
        end
    end
end