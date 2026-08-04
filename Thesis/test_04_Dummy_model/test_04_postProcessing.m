%% =====================================================================
%  POST-PROCESSING: FOM vs ROM (MC, Rubin, MCB & MN) for Shock Simulations
% =====================================================================
clear; close all; clc;

% 1. Folder Selection
results_dir = uigetdir(pwd, 'Select the results folder (e.g., Shock_...)');
if results_dir == 0
    error('No folder selected. Exiting.');
end
fprintf('Selected directory: %s\n', results_dir);

% 2. Find all FOM files to extract Q and K combinations
fom_files = dir(fullfile(results_dir, 'FOM_*.mat'));
if isempty(fom_files)
    error('No FOM files found in the selected directory.');
end

% Initialize log file for GRE results
log_file = fopen(fullfile(results_dir, 'GRE_Results.txt'), 'w');
fprintf(log_file, '======================================================\n');
fprintf(log_file, ' GLOBAL RELATIVE ERROR (GRE) & TIME REPORT\n');
fprintf(log_file, '======================================================\n\n');

% 3. Loop over FOM files (which define the base Q and K combinations)
for i_fom = 1:length(fom_files)
    fom_name = fom_files(i_fom).name;
    
    % Extract Q and K from the filename
    tokens = regexp(fom_name, 'FOM_Q(\d+)_K([\d\.]+)\.mat', 'tokens');
    if isempty(tokens), continue; end
    
    Q_val = str2double(tokens{1}{1});
    K_val = str2double(tokens{1}{2});
    
    fprintf('\nAnalyzing case: Q = %d, K_mult = %g\n', Q_val, K_val);
    fprintf(log_file, '>>> CASE: Q = %d | K_mult = %g <<<\n', Q_val, K_val);
    
    % Load FOM data
    fom_data = load(fullfile(results_dir, fom_name));
    t_fom = fom_data.t_fom;
    y_X_fom = fom_data.y_contact_nodes_X_fom; % Dimensions: [n_nodes, n_time]
    y_X_fom_node1 = y_X_fom(1, :);
    
    % Estrai la CPU time del FOM (con check di sicurezza)
    if isfield(fom_data, 'cpu_time')
        fom_cpu = fom_data.cpu_time;
        fom_legend = sprintf('FOM (On: %.2fs)', fom_cpu);
    else
        fom_cpu = NaN;
        fom_legend = 'FOM';
    end
    
    % Setup Figure (Height incrementata per accomodare 4 subplot)
    fig = figure('Name', sprintf('Contact Node 1 - Q%d - K%g', Q_val, K_val), ...
                 'NumberTitle', 'off', 'Position', [100, 50, 1200, 1200], 'Color', 'w');
             
    % --- Subplot 1: Milman-Chu ---
    ax1 = subplot(4,1,1);
    plot(t_fom, y_X_fom_node1, 'k-', 'LineWidth', 2, 'DisplayName', fom_legend);
    hold on; grid on;
    yline(1.5e-6, 'r-.', 'LineWidth', 1.5, 'DisplayName', 'Wall (1.5 \mum)');
    title(sprintf('Milman-Chu Method (Q=%d, K_{mult}=%g)', Q_val, K_val));
    xlabel('Time [s]'); ylabel('X Displacement [m]');
    
    % --- Subplot 2: Rubin ---
    ax2 = subplot(4,1,2);
    plot(t_fom, y_X_fom_node1, 'k-', 'LineWidth', 2, 'DisplayName', fom_legend);
    hold on; grid on;
    yline(1.5e-6, 'r-.', 'LineWidth', 1.5, 'DisplayName', 'Wall (1.5 \mum)');
    title(sprintf('Rubin Method (Q=%d, K_{mult}=%g)', Q_val, K_val));
    xlabel('Time [s]'); ylabel('X Displacement [m]');
    
    % --- Subplot 3: Massless CB (MCB) ---
    ax3 = subplot(4,1,3);
    plot(t_fom, y_X_fom_node1, 'k-', 'LineWidth', 2, 'DisplayName', fom_legend);
    hold on; grid on;
    yline(1.5e-6, 'r-.', 'LineWidth', 1.5, 'DisplayName', 'Wall (1.5 \mum)');
    title(sprintf('Massless CB Method (Q=%d, K_{mult}=%g)', Q_val, K_val));
    xlabel('Time [s]'); ylabel('X Displacement [m]');
    
    % --- Subplot 4: MacNeal (MN) ---
    ax4 = subplot(4,1,4);
    plot(t_fom, y_X_fom_node1, 'k-', 'LineWidth', 2, 'DisplayName', fom_legend);
    hold on; grid on;
    yline(1.5e-6, 'r-.', 'LineWidth', 1.5, 'DisplayName', 'Wall (1.5 \mum)');
    title(sprintf('MacNeal Method (Q=%d, K_{mult}=%g)', Q_val, K_val));
    xlabel('Time [s]'); ylabel('X Displacement [m]');
    
    % Find and Plot corresponding ROM files for this (Q,K) combination
    rom_files = dir(fullfile(results_dir, sprintf('ROM_*_Q%04d_K%g.mat', Q_val, K_val)));
    
    % Dynamic color array based on the number of files found (approximate)
    colors = lines(15); 
    color_idx_mc = 1; color_idx_rubin = 1; color_idx_mcb = 1; color_idx_mn = 1;
    
    for i_rom = 1:length(rom_files)
        rom_name = rom_files(i_rom).name;
        rom_data = load(fullfile(results_dir, rom_name));
        t_rom = rom_data.t_rom;
        
        % Estrai la CPU time (Online) del ROM
        if isfield(rom_data, 'cpu_time')
            rom_cpu = rom_data.cpu_time;
        else
            rom_cpu = NaN;
        end
        
        % Estrai l'Offline time del ROM
        if isfield(rom_data, 'offline_time')
            rom_offline = rom_data.offline_time;
        else
            rom_offline = NaN;
        end
        
        % Identify Method and extract Phi modes
        phi_tokens = regexp(rom_name, 'Phi(\d+)', 'tokens');
        phi_val = str2double(phi_tokens{1}{1});
        
        % LOGICA DI RICONOSCIMENTO
        is_mcb = contains(rom_name, 'ROM_MCB');
        is_mc = contains(rom_name, 'ROM_MC_'); 
        is_rubin = contains(rom_name, 'ROM_Rubin');
        is_mn = contains(rom_name, 'ROM_MN_');
        
        if is_mcb
            y_X_rom = rom_data.y_contact_nodes_X_romMCB;
            current_ax = ax3;
            col = colors(color_idx_mcb, :);
            color_idx_mcb = color_idx_mcb + 1;
            method_str = 'MCB';
        elseif is_mc
            y_X_rom = rom_data.y_contact_nodes_X_romMC;
            current_ax = ax1;
            col = colors(color_idx_mc, :);
            color_idx_mc = color_idx_mc + 1;
            method_str = 'MC';
        elseif is_rubin
            y_X_rom = rom_data.y_contact_nodes_X_romRubin;
            current_ax = ax2;
            col = colors(color_idx_rubin, :);
            color_idx_rubin = color_idx_rubin + 1;
            method_str = 'Rubin';
        elseif is_mn
            y_X_rom = rom_data.y_contact_nodes_X_romMN;
            current_ax = ax4;
            col = colors(color_idx_mn, :);
            color_idx_mn = color_idx_mn + 1;
            method_str = 'MN';
        else
            continue; % File non riconosciuto
        end
        
        % INTERPOLATION: ode15s time steps differ. Interpolate ROM onto t_fom
        y_X_rom_interp = interp1(t_rom, y_X_rom', t_fom, 'linear', 'extrap')';
        y_X_rom_node1_interp = y_X_rom_interp(1, :);
        
        % GRE CALCULATION (Global Relative Error)
        % GRE Node 1
        num_node1 = norm(y_X_fom_node1 - y_X_rom_node1_interp);
        den_node1 = norm(y_X_fom_node1);
        gre_node1 = (num_node1 / den_node1) * 100; % Percentage
        
        % GRE All contact nodes
        num_all = norm(y_X_fom(:) - y_X_rom_interp(:));
        den_all = norm(y_X_fom(:));
        gre_all = (num_all / den_all) * 100; % Percentage
        
        % Costruzione dinamica delle informazioni sui tempi per la legenda
        time_info = '';
        if ~isnan(rom_offline)
            time_info = sprintf('Off: %.2fs', rom_offline);
        end
        if ~isnan(rom_cpu)
            if isempty(time_info)
                time_info = sprintf('On: %.2fs', rom_cpu);
            else
                time_info = sprintf('%s, On: %.2fs', time_info, rom_cpu);
            end
        end
        
        % PLOT with GRE and Time in Legend
        if isempty(time_info)
            legend_str = sprintf('ROM \\phi=%d (GRE: %.2f%%)', phi_val, gre_all);
        else
            legend_str = sprintf('ROM \\phi=%d (GRE: %.2f%%, %s)', phi_val, gre_all, time_info);
        end
        
        plot(current_ax, t_fom, y_X_rom_node1_interp, '--', 'Color', col, ...
             'LineWidth', 1.5, 'DisplayName', legend_str);
         
        % Preparazione stringhe per i log
        if isnan(rom_offline), off_str = 'N/A'; else, off_str = sprintf('%5.2fs', rom_offline); end
        if isnan(rom_cpu), on_str = 'N/A'; else, on_str = sprintf('%5.2fs', rom_cpu); end
        
        % Console and File Logging
        fprintf('  [%-5s] Phi: %3d | GRE Node 1: %6.3f%% | Global GRE: %6.3f%% | Off: %s | On: %s\n', ...
                method_str, phi_val, gre_node1, gre_all, off_str, on_str);
        fprintf(log_file, '  %-5s Phi: %03d | GRE Node 1: %6.3f%% | GRE All Nodes: %6.3f%% | Off: %s | On: %s\n', ...
                method_str, phi_val, gre_node1, gre_all, off_str, on_str);
    end
    
    % Finalize Figures
    legend(ax1, 'Location', 'best');
    legend(ax2, 'Location', 'best');
    legend(ax3, 'Location', 'best');
    legend(ax4, 'Location', 'best');
    
    % Save Figures
    fig_filename = fullfile(results_dir, sprintf('Compare_Q%d_K%g.png', Q_val, K_val));
    exportgraphics(fig, fig_filename, 'Resolution', 300);
    savefig(fig, fullfile(results_dir, sprintf('Compare_Q%d_K%g.fig', Q_val, K_val)));
    
    fprintf('  -> Plot saved as %s\n', fig_filename);
    fprintf(log_file, '\n');
end

fclose(log_file);
fprintf('\nPost-processing complete! Results are saved in %s\n', results_dir);