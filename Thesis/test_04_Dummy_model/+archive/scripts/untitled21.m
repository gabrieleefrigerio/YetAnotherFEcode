%% =====================================================================
%  POST-PROCESSING: FOM vs ALL ROMs (Only 200 Modes) for Shock V2
% =====================================================================
clear; close all; clc;

% --- Configurazione Analisi ---
target_phi = 200; % Numero di modi fisso da confrontare

% --- Mappatura facce, direzioni e gap ---
faces = {'T', 'B', 'L', 'R'};
dirs  = {'Y', 'Y', 'X', 'X'};
gaps  = [5e-6, -1.5e-6, -5e-6, 1.5e-6];

% --- Aggiunto il nuovo metodo MT_ ---
methods_id = {'MC_', 'Rubin', 'MCB', 'MN_', 'MT_'};
method_names = {'Milman-Chu', 'Rubin', 'Massless CB', 'MacNeal', 'Modal Truncation'};

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
log_file = fopen(fullfile(results_dir, sprintf('GRE_Results_Comparison_Phi%d.txt', target_phi)), 'w');
fprintf(log_file, '======================================================\n');
fprintf(log_file, ' GLOBAL RELATIVE ERROR (GRE) & TIME REPORT - PHI = %d\n', target_phi);
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
    t_fom = fom_data.t_fom(:); % Forza colonna
    
    if isfield(fom_data, 'cpu_time')
        fom_cpu = fom_data.cpu_time;
        fom_legend = sprintf('FOM (On: %.2fs)', fom_cpu);
    else
        fom_cpu = NaN;
        fom_legend = 'FOM';
    end
    
    % --- Inizializza la singola figura per questo FOM ---
    fig = figure('Name', sprintf('All Methods - Q%d - K%g - Phi%d', Q_val, K_val, target_phi), ...
                 'NumberTitle', 'off', 'Position', [100, 50, 1000, 1400], 'Color', 'w');
    axs = gobjects(5,1);
    
    % Setup dei primi 4 subplot in colonna (Facce) con i dati FOM
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
            title(axs(f), sprintf('Face %s non presente nei dati FOM', face));
        end
        
        xlabel(axs(f), 'Time [s]'); 
        ylabel(axs(f), 'Displacement [m]');
    end
    
    % Setup del 5° subplot (GRE nel tempo)
    axs(5) = subplot(5, 1, 5);
    hold(axs(5), 'on'); grid(axs(5), 'on');
    title(axs(5), sprintf('Global Relative Error (GRE) Over Time - \\phi = %d', target_phi));
    xlabel(axs(5), 'Time [s]');
    ylabel(axs(5), 'GRE [%]');
    
    % Cerca tutti i file ROM per questa coppia Q, K
    rom_files = dir(fullfile(results_dir, sprintf('ROM_*_Q%04d_K%g.mat', Q_val, K_val)));
    
    % Colori per i diversi metodi
    colors = lines(length(methods_id));
    
    % --- Loop sui metodi per sovrapporre i risultati a 200 modi ---
    for m = 1:length(methods_id)
        current_method_id = methods_id{m};
        current_method_name = method_names{m};
        
        % Trova il file ROM specifico per questo metodo e per target_phi
        target_rom_file = '';
        for r = 1:length(rom_files)
            rom_name_tmp = rom_files(r).name;
            if contains(rom_name_tmp, ['ROM_' current_method_id])
                phi_tokens = regexp(rom_name_tmp, 'Phi(\d+)', 'tokens');
                if ~isempty(phi_tokens) && str2double(phi_tokens{1}{1}) == target_phi
                    target_rom_file = rom_name_tmp;
                    break; % Trovato il file corretto per questo metodo
                end
            end
        end
        
        % Se non c'è un ROM a 200 modi per questo metodo, salta
        if isempty(target_rom_file)
            fprintf('  [%-5s] Nessun file trovato per Phi = %d\n', current_method_name, target_phi);
            continue;
        end
        
        % Load ROM data
        rom_data = load(fullfile(results_dir, target_rom_file));
        t_rom = rom_data.t_rom(:); % Forza colonna
        
        if isfield(rom_data, 'cpu_time'), rom_cpu = rom_data.cpu_time; else, rom_cpu = NaN; end
        if isfield(rom_data, 'offline_time'), rom_offline = rom_data.offline_time; else, rom_offline = NaN; end
        
        % Estrazione matrici complete per calcolare il GRE
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
        
        % --- CALCOLO GRE GLOBALE SCALARE ---
        gre_all = (norm(y_fom_cat - y_rom_cat, 'fro') / norm(y_fom_cat, 'fro')) * 100;
        
        % --- CALCOLO GRE NEL TEMPO ---
        norm_diff_t = sqrt(sum((y_fom_cat - y_rom_cat).^2, 1));
        max_fom_global  = max(sqrt(sum(y_fom_cat.^2, 1))); 
        gre_t = (norm_diff_t ./ (max_fom_global + eps)) * 100;
        
        % Costruzione Legenda ROM
        time_info = '';
        if ~isnan(rom_offline), time_info = sprintf('Off: %.2fs', rom_offline); end
        if ~isnan(rom_cpu)
            if isempty(time_info), time_info = sprintf('On: %.2fs', rom_cpu);
            else, time_info = sprintf('%s, On: %.2fs', time_info, rom_cpu); end
        end
        
        if isempty(time_info)
            legend_str = sprintf('%s (GRE: %.2f%%)', current_method_name, gre_all);
        else
            legend_str = sprintf('%s (GRE: %.2f%%, %s)', current_method_name, gre_all, time_info);
        end
        
        % Plot spostamenti nei 4 subplot superiori
        for f = 1:4
            face = faces{f};
            dir_plot = dirs{f};
            
            if isfield(rom_data.y_contact, face)
                y_rom_face = rom_data.y_contact.(face).(dir_plot)(1, :);
                y_rom_node1_interp = interp1(t_rom, y_rom_face(:), t_fom, 'linear', 'extrap');
                
                if f == 1
                    plot(axs(f), t_fom, y_rom_node1_interp, '--', 'Color', colors(m, :), ...
                         'LineWidth', 1.5, 'DisplayName', legend_str);
                else
                    plot(axs(f), t_fom, y_rom_node1_interp, '--', 'Color', colors(m, :), ...
                         'LineWidth', 1.5, 'HandleVisibility', 'off');
                end
            end
        end
        
        % Plot dell'andamento dell'errore (GRE nel tempo) nel 5° subplot
        plot(axs(5), t_fom, gre_t(:), '-', 'Color', colors(m, :), ...
             'LineWidth', 1.5, 'DisplayName', current_method_name);
        
        % Logging su file e console
        if isnan(rom_offline), off_str = 'N/A'; else, off_str = sprintf('%5.2fs', rom_offline); end
        if isnan(rom_cpu), on_str = 'N/A'; else, on_str = sprintf('%5.2fs', rom_cpu); end
        
        fprintf('  [%-16s] Avg GRE: %6.3f%% | Off: %s | On: %s\n', ...
                current_method_name, gre_all, off_str, on_str);
        fprintf(log_file, '  %-16s | Avg GRE: %6.3f%% | Off: %s | On: %s\n', ...
                current_method_name, gre_all, off_str, on_str);
    end
    
    % Aggiungi la legenda al primo subplot (Displacements) e al quinto (GRE)
    legend(axs(1), 'Location', 'best');
    legend(axs(5), 'Location', 'best');
    sgtitle(fig, sprintf('Methods Comparison (\\phi = %d) - Q=%d, K_{mult}=%g', target_phi, Q_val, K_val), 'FontSize', 16, 'FontWeight', 'bold');
    
    % Salvataggio figure
    fig_filename = fullfile(results_dir, sprintf('Compare_ALL_Phi%d_Q%d_K%g.png', target_phi, Q_val, K_val));
    exportgraphics(fig, fig_filename, 'Resolution', 300);
    savefig(fig, fullfile(results_dir, sprintf('Compare_ALL_Phi%d_Q%d_K%g.fig', target_phi, Q_val, K_val)));
    
    fprintf('  -> Plot saved as %s\n', fig_filename);
    fprintf(log_file, '\n');
end

fclose(log_file);
fprintf('\nPost-processing complete! Results are saved in %s\n', results_dir);