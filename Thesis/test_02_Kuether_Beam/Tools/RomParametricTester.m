classdef RomParametricTester < handle
    properties
        % Modello Completo (FOM) e parametri base
        BeamAssembly
        MassFull
        NlDof
        Spring
        
        % Parametri di Simulazione
        PulseDuration = 0.001;
        Tmax = 0.2;
        TimeStep = 5e-6;
        
        % USIAMO UN CELL ARRAY {} INVECE DI UNO STRUCT ARRAY.
        % Questo impedisce a MATLAB di impacchettare gli oggetti per errore.
        RomConfigs = {}; 
    end
    
    methods
        % Costruttore: inizializza il tester con il modello completo (FOM)
        function obj = RomParametricTester(beam_asm, m_full, nl_dof, spring)
            obj.BeamAssembly = beam_asm;
            obj.MassFull = m_full;
            obj.NlDof = nl_dof;
            obj.Spring = spring;
        end
        
        % Metodo per aggiungere progressivamente i modelli ridotti
        function add_rom(obj, rom_asm, rom_p, label)
            % Creiamo una singola struct e la ficchiamo dritta nella cella successiva
            new_rom = struct('Assembly', rom_asm, 'P', rom_p, 'Label', label);
            obj.RomConfigs{end+1} = new_rom;
        end
        
        % Metodo principale per eseguire lo sweep sulle ampiezze
        function run_sweep(obj, amplitudes)
            for i = 1:length(amplitudes)
                A = amplitudes(i);
                obj.print_header(A);
                obj.setup_figure(A);
                
                % Inizializziamo flag per il FOM
                fom_solved = false;
                fom_disp = [];
                
                for j = 1:length(obj.RomConfigs)
                    % Estraiamo la struct corrente dal cell array usando le GRAFFE
                    cur_rom = obj.RomConfigs{j};
                    
                    % Inizializza il simulatore per la coppia FOM/ROM corrente
                    sim = TransientSimulator(obj.BeamAssembly, cur_rom.Assembly, ...
                                             obj.Spring, cur_rom.P, obj.MassFull, obj.NlDof);
                    
                    % Risolve il ROM
                    sim.solve_rom(A, obj.PulseDuration, obj.Tmax, obj.TimeStep);
                    
                    % Risolve e plotta il FOM solo una volta per ampiezza
                    if ~fom_solved
                        sim.solve_fom(A, obj.PulseDuration, obj.Tmax, obj.TimeStep);
                        fom_disp = sim.disp_center_fom;
                        
                        plot(sim.time_fom, fom_disp * 1000, 'k-', 'LineWidth', 2, 'DisplayName', 'FOM (Reference)');
                        obj.print_row('FOM', sim.energy_fom, sim.max_disp_fom, sim.rms_disp_fom, NaN);
                        
                        fom_solved = true;
                    end
                    
                    % Calcola Errore L2
                    err_perc = norm(fom_disp - sim.disp_center_rom) / norm(fom_disp) * 100;
                    
                    % Plot ROM
                    plot(sim.time_rom, sim.disp_center_rom * 1000, '--', 'LineWidth', 1.5, 'DisplayName', cur_rom.Label);
                    
                    % Stampa riga ROM
                    obj.print_row(cur_rom.Label, sim.energy_rom, sim.max_disp_rom, sim.rms_disp_rom, err_perc);
                end
                
                obj.finalize_figure(A);
                fprintf(repmat('=', 1, 75)); fprintf('\n');
            end
        end
        
    end
    
    % Metodi privati per l'estetica
    methods (Access = private)
        function print_header(~, A)
            fprintf('\n=== IMPATTO: %.2f N ===\n', A);
            fprintf('%-14s | %-12s | %-13s | %-13s | %-12s\n', 'Modello', 'Energy [J]', 'Max Disp [mm]', 'RMS Disp [mm]', 'L2 Error [%]');
            fprintf(repmat('-', 1, 75)); fprintf('\n');
        end
        
        function print_row(~, label, energy, max_d, rms_d, err)
            if isnan(err)
                err_str = '-';
            else
                err_str = sprintf('%.4f %%', err);
            end
            fprintf('%-14s | %.4e  | %-13.4f | %-13.4f | %-12s\n', ...
                label, energy, max_d * 1000, rms_d * 1000, err_str);
        end
        
        function setup_figure(obj, A)
            figure('Name', sprintf('Transient A=%.2fN', A), 'Color', 'w');
            hold on; grid on;
            yline(obj.Spring.a * 1000, 'k-.', 'Gap Clearance', 'LineWidth', 1.5, 'HandleVisibility', 'off');
        end
        
        function finalize_figure(~, A)
            title(sprintf('Midpoint Displacement (Amplitude = %.2f N)', A), 'FontSize', 12);
            xlabel('Time [s]'); ylabel('Displacement [mm]');
            legend('Location', 'best');
            xlim('tight');
            hold off;
        end
    end
end