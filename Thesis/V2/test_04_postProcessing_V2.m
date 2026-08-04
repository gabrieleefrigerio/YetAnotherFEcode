%% =====================================================================
%  POST-PROCESSING: FOM vs ROM (MC, Rubin, MCB & MN) for Shock V2
% =====================================================================
clear; close all; clc;

% --- Mappatura facce, direzioni e gap ---
faces = {'T', 'B', 'L', 'R'};
dirs  = {'Y', 'Y', 'X', 'X'};
gaps  = [5e-6, -1.5e-6, -5e-6, 1.5e-6];

methods_id = {'MC_', 'Rubin', 'MCB', 'MN_'};
method_names = {'Milman-Chu', 'Rubin', 'Massless CB', 'MacNeal'};

% 1. Folder Selection
results_dir = uigetdir(pwd, 'Select the results folder (e.g., Shock_...)');
if results_dir == 0
    error('No folder selected. Exiting.');
end
fprintf('Selected directory: %s\n', results_dir);

% 2. Find all FOM files
fom_files = dir(fullfile(results_dir, 'FOM_*.mat'));
if isempty(fom_files)
    error('No FOM files found in the selected directory.');
end

% Initialize log file for GRE results
log_file = fopen(fullfile(results_dir, 'GRE_Results_V2.txt'), 'w');
fprintf(log_file, '======================================================\n');
fprintf(log_file, ' GLOBAL RELATIVE ERROR (GRE) & TIME REPORT\n');
fprintf(log_file, '======================================================\n\n');

% 3. Loop over FOM files
for i_fom = 1:length(fom_files)
    fom_name = fom_files(i_fom).name;
    
    tokens = regexp(fom_name, 'FOM_Q(\d+)_K([\d\.]+)\.mat', 'tokens');
    if isempty(tokens), continue; end
    
    Q_val = str2double(tokens{1}{1});
    K_val = str2double(tokens{1}{2});
    
    fprintf('\n======================================================\n');
    fprintf('Analyzing case: Q = %d, K_mult = %g\n', Q_val, K_val);
    fprintf('======================================================\n');
    fprintf(log_file, '>>> CASE: Q = %d | K_mult = %g <<<\n', Q_val, K_val);
    
    % Load FOM data
    fom_data = load(fullfile(results_dir, fom_name));
    t_fom = fom_data.t_fom(:); % Forza colonna per interpolazioni sicure
    
    if isfield(fom_data, 'cpu_time')
        fom_cpu = fom_data.cpu_time;
        fom_legend = sprintf('FOM (On: %.2fs)', fom_cpu);
    else
        fom_cpu = NaN;
        fom_legend = 'FOM';
    end
    
    % Cerca tutti i file ROM per questa coppia Q, K
    rom_files = dir(fullfile(results_dir, sprintf('ROM_*_Q%04d_K%g.mat', Q_val, K_val)));
    
    % Loop sui metodi per creare una figura per metodo
    for m = 1:4
        current_method_id = methods_id{m};
        current_method_name = method_names{m};
        
        % Trova i file ROM che appartengono a questo metodo
        method_rom_files = {};
        for r = 1:length(rom_files)
            if contains(rom_files(r).name, ['ROM_' current_method_id])
                method_rom_files{end+1} = rom_files(r).name; %#ok<SAGROW>
            end
        end
        
        % Se non ci sono ROM per questo metodo, salta
        if isempty(method_rom_files)
            continue;
        end
        
        % Inizializza la figura per il metodo corrente
        % Altezza incrementata (1400) per gestire comodamente 5 subplot
        fig = figure('Name', sprintf('%s - Q%d - K%g', current_method_name, Q_val, K_val), ...
                     'NumberTitle', 'off', 'Position', [100, 50, 1000, 1400], 'Color', 'w');
        
        axs = gobjects(5,1);
        
        % Setup dei primi 4 subplot in colonna (Facce)
        for f = 1:4
            axs(f) = subplot(5, 1, f);
            hold(axs(f), 'on'); grid(axs(f), 'on');
            
            face = faces{f};
            dir_plot = dirs{f};
            gap_plot = gaps(f);
            
            if isfield(fom_data.y_contact, face)
                y_fom = fom_data.y_contact.(face).(dir_plot)(1, :);
                
                if f == 1
                    plot(axs(f), t_fom, y_fom(:), 'k-', 'LineWidth', 2, 'DisplayName', fom_legend);
                    yline(axs(f), gap_plot, 'r-.', 'LineWidth', 1.5, 'DisplayName', 'Wall Gap');
                else
                    plot(axs(f), t_fom, y_fom(:), 'k-', 'LineWidth', 2, 'HandleVisibility', 'off');
                    yline(axs(f), gap_plot, 'r-.', 'LineWidth', 1.5, 'HandleVisibility', 'off');
                end
                
                title(axs(f), sprintf('Face %s (Dir: %s)', face, dir_plot));
            else
                title(axs(f), sprintf('Face %s non presente nei dati', face));
            end
            
            xlabel(axs(f), 'Time [s]'); 
            ylabel(axs(f), 'Displacement [m]');
        end
        
        % Setup del 5° subplot (GRE nel tempo)
        axs(5) = subplot(5, 1, 5);
        hold(axs(5), 'on'); grid(axs(5), 'on');
        title(axs(5), 'Global Relative Error (GRE) Over Time');
        xlabel(axs(5), 'Time [s]');
        ylabel(axs(5), 'GRE [%]');
        
        % Plot delle curve ROM per il metodo corrente
        colors = lines(length(method_rom_files));
        
        for i_rom = 1:length(method_rom_files)
            rom_name = method_rom_files{i_rom};
            rom_data = load(fullfile(results_dir, rom_name));
            t_rom = rom_data.t_rom(:); % Forza colonna
            
            if isfield(rom_data, 'cpu_time'), rom_cpu = rom_data.cpu_time; else, rom_cpu = NaN; end
            if isfield(rom_data, 'offline_time'), rom_offline = rom_data.offline_time; else, rom_offline = NaN; end
            
            phi_tokens = regexp(rom_name, 'Phi(\d+)', 'tokens');
            phi_val = str2double(phi_tokens{1}{1});
            
            % Estrazione matrici complete di tutti i nodi di contatto per calcolare il GRE
            y_fom_cat = []; y_rom_cat = [];
            for f_idx = 1:4
                fc = faces{f_idx};
                dc = dirs{f_idx};
                if isfield(fom_data.y_contact, fc) && isfield(rom_data.y_contact, fc)
                    y_f = fom_data.y_contact.(fc).(dc);
                    y_r = rom_data.y_contact.(fc).(dc);
                    
                    y_fom_cat = [y_fom_cat; y_f]; %#ok<AGROW>
                    y_rom_interp_all = interp1(t_rom, y_r', t_fom, 'linear', 'extrap')';
                    y_rom_cat = [y_rom_cat; y_rom_interp_all]; %#ok<AGROW>
                end
            end
            % --- CALCOLO GRE GLOBALE SCALARE (su tutta la simulazione) ---
            gre_all = (norm(y_fom_cat - y_rom_cat, 'fro') / norm(y_fom_cat, 'fro')) * 100;
            

            % --- CALCOLO GRE NEL TEMPO (Istantaneo Normalizzato) ---
            norm_diff_t = sqrt(sum((y_fom_cat - y_rom_cat).^2, 1));
            max_fom_global  = max(sqrt(sum(y_fom_cat.^2, 1))); % Massimo assoluto
            gre_t = (norm_diff_t ./ (max_fom_global + eps)) * 100;
            
            % Costruzione Legenda ROM
            time_info = '';
            if ~isnan(rom_offline), time_info = sprintf('Off: %.2fs', rom_offline); end
            if ~isnan(rom_cpu)
                if isempty(time_info), time_info = sprintf('On: %.2fs', rom_cpu);
                else, time_info = sprintf('%s, On: %.2fs', time_info, rom_cpu); end
            end
            
            if isempty(time_info)
                legend_str = sprintf('ROM \\phi=%d (Avg GRE: %.2f%%)', phi_val, gre_all);
            else
                legend_str = sprintf('ROM \\phi=%d (Avg GRE: %.2f%%, %s)', phi_val, gre_all, time_info);
            end
            
            % Plot spostamenti nei 4 subplot superiori
            for f = 1:4
                face = faces{f};
                dir_plot = dirs{f};
                
                if isfield(rom_data.y_contact, face)
                    y_rom_face = rom_data.y_contact.(face).(dir_plot)(1, :);
                    y_rom_node1_interp = interp1(t_rom, y_rom_face(:), t_fom, 'linear', 'extrap');
                    
                    if f == 1
                        plot(axs(f), t_fom, y_rom_node1_interp, '--', 'Color', colors(i_rom, :), ...
                             'LineWidth', 1.5, 'DisplayName', legend_str);
                    else
                        plot(axs(f), t_fom, y_rom_node1_interp, '--', 'Color', colors(i_rom, :), ...
                             'LineWidth', 1.5, 'HandleVisibility', 'off');
                    end
                end
            end
            
            % Plot dell'andamento dell'errore (GRE nel tempo) nel 5° subplot
            plot(axs(5), t_fom, gre_t(:), '-', 'Color', colors(i_rom, :), ...
                 'LineWidth', 1.5, 'HandleVisibility', 'off');
            
            % Logging su file e console
            if isnan(rom_offline), off_str = 'N/A'; else, off_str = sprintf('%5.2fs', rom_offline); end
            if isnan(rom_cpu), on_str = 'N/A'; else, on_str = sprintf('%5.2fs', rom_cpu); end
            
            fprintf('  [%-5s] Phi: %3d | Avg GRE: %6.3f%% | Off: %s | On: %s\n', ...
                    current_method_name, phi_val, gre_all, off_str, on_str);
            fprintf(log_file, '  %-5s Phi: %03d | Avg GRE: %6.3f%% | Off: %s | On: %s\n', ...
                    current_method_name, phi_val, gre_all, off_str, on_str);
        end
        
        % Aggiungi la legenda solo al primo subplot
        legend(axs(1), 'Location', 'best');
        sgtitle(fig, sprintf('%s Method (Q=%d, K_{mult}=%g)', current_method_name, Q_val, K_val), 'FontSize', 16, 'FontWeight', 'bold');
        
        % Salvataggio figure
        clean_method_name = strrep(current_method_name, ' ', '');
        fig_filename = fullfile(results_dir, sprintf('Compare_%s_Q%d_K%g.png', clean_method_name, Q_val, K_val));
        exportgraphics(fig, fig_filename, 'Resolution', 300);
        savefig(fig, fullfile(results_dir, sprintf('Compare_%s_Q%d_K%g.fig', clean_method_name, Q_val, K_val)));
        
        fprintf('  -> Plot saved as %s\n', fig_filename);
    end
    fprintf(log_file, '\n');
end

fclose(log_file);
fprintf('\nPost-processing complete! Results are saved in %s\n', results_dir);