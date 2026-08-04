%% =====================================================================
%  POST-PROCESSING: SELECTIVE MILMAN-CHU (Time-History Analysis)
% =====================================================================
clear; close all; clc;

% 1. Folder Selection
results_dir = uigetdir(pwd, 'Select the Selective MC Grid Search results folder');
if results_dir == 0
    error('No folder selected. Exiting.');
end
fprintf('Selected directory: %s\n', results_dir);

% 2. Find all FOM files to extract K combinations
fom_files = dir(fullfile(results_dir, 'FOM_Ref_K*.mat'));
if isempty(fom_files)
    error('No FOM files found in the selected directory.');
end

% Initialize log file for GRE results
log_file = fopen(fullfile(results_dir, 'GRE_Results.txt'), 'w');
fprintf(log_file, '======================================================\n');
fprintf(log_file, ' GLOBAL RELATIVE ERROR (GRE) & TIME REPORT\n');
fprintf(log_file, '======================================================\n\n');

% 3. Loop over FOM files (which define the base K combination)
for i_fom = 1:length(fom_files)
    fom_name = fom_files(i_fom).name;
    
    % Extract K from the filename
    tokens = regexp(fom_name, 'FOM_Ref_K([\d\.]+)\.mat', 'tokens');
    if isempty(tokens), continue; end
    
    k_mult = str2double(tokens{1}{1});
    
    fprintf('\n======================================================\n');
    fprintf('Analyzing case: K_mult = %g\n', k_mult);
    fprintf('======================================================\n');
    fprintf(log_file, '>>> CASE: K_mult = %g <<<\n\n', k_mult);
    
    % Matrix to collect data for THIS K_mult summary table
    % Columns: [Phi, Skip, N_Vecs, GRE_All, Off_Time, On_Time]
    % (We use Phi = -1 for FOM so it sorts to the top)
    summary_data_k = [];
    
    % Load FOM data
    fom_data = load(fullfile(results_dir, fom_name));
    t_fom = fom_data.t_fom;
    y_X_fom = fom_data.y_contact_nodes_X_fom; 
    y_X_fom_node1 = y_X_fom(1, :);
    
    % Extract FOM CPU time
    if isfield(fom_data, 'cpu_time')
        fom_cpu = fom_data.cpu_time;
        fom_legend = sprintf('FOM (On: %.2fs)', fom_cpu);
    else
        fom_cpu = NaN;
        fom_legend = 'FOM';
    end
    
    % Add FOM row to summary (Phi = -1)
    summary_data_k = [summary_data_k; -1, 0, NaN, 0.0, NaN, fom_cpu];
    
    % 4. Find all UNIQUE Phi values for this K_mult
    pattern_all = sprintf('ROM_MC_K%g_Phi*_Skip*_NVec*.mat', k_mult);
    rom_files_all = dir(fullfile(results_dir, pattern_all));
    unique_phis = [];
    for r = 1:length(rom_files_all)
        phi_tokens = regexp(rom_files_all(r).name, 'Phi(\d+)', 'tokens');
        if ~isempty(phi_tokens)
            unique_phis(end+1) = str2double(phi_tokens{1}{1});
        end
    end
    unique_phis = unique(unique_phis);
    
    if isempty(unique_phis)
        fprintf('No ROM files found for K_mult = %g\n', k_mult);
        continue;
    end
    
    % 5. Loop over each unique Phi to create a SEPARATE figure and collect data
    for p_idx = 1:length(unique_phis)
        curr_phi = unique_phis(p_idx);
        fprintf('\n  --- Generating plot for Phi = %d ---\n', curr_phi);
        
        % Setup Figure (Single axis per Phi)
        fig_name = sprintf('Contact Node 1 - K%g - Phi%d', k_mult, curr_phi);
        fig = figure('Name', fig_name, 'NumberTitle', 'off', 'Position', [100, 50, 1200, 800]);
                 
        ax1 = axes(fig);
        plot(ax1, t_fom, y_X_fom_node1, 'k-', 'LineWidth', 2.5, 'DisplayName', fom_legend);
        hold on; grid on;
        yline(ax1, 1.5e-6, 'r-.', 'LineWidth', 2, 'DisplayName', 'Wall (1.5 \mum)');
        title(ax1, sprintf('Selective Milman-Chu Method (K_{mult}=%g, \\Phi=%d)', k_mult, curr_phi), 'FontSize', 14);
        xlabel(ax1, 'Time [s]', 'FontWeight', 'bold'); 
        ylabel(ax1, 'X Displacement [m]', 'FontWeight', 'bold');
        
        % Find only ROMs with this K_mult and this Phi
        pattern_phi = sprintf('ROM_MC_K%g_Phi%03d_Skip*_NVec*.mat', k_mult, curr_phi);
        rom_files_phi = dir(fullfile(results_dir, pattern_phi));
        
        % Dynamic color array for Skip variants
        colors = lines(length(rom_files_phi)); 
        
        for i_rom = 1:length(rom_files_phi)
            rom_name = rom_files_phi(i_rom).name;
            rom_data = load(fullfile(results_dir, rom_name));
            t_rom = rom_data.t_rom;
            
            % Times
            rom_cpu = NaN;
            if isfield(rom_data, 'cpu_time'), rom_cpu = rom_data.cpu_time; end
            rom_offline = NaN;
            if isfield(rom_data, 'offline_time'), rom_offline = rom_data.offline_time; end
            
            % Parameters
            skip_val = rom_data.skip;
            nvec_val = rom_data.num_enrich;
            y_X_rom = rom_data.y_contact_nodes_X_romMC;
            
            % INTERPOLATION
            y_X_rom_interp = interp1(t_rom, y_X_rom', t_fom, 'linear', 'extrap')';
            y_X_rom_node1_interp = y_X_rom_interp(1, :);
            
            % CALCULATION: Only Global Relative Error
            num_all = norm(y_X_fom(:) - y_X_rom_interp(:));
            den_all = norm(y_X_fom(:));
            gre_all = (num_all / den_all) * 100; % All contact nodes
            
            % Save data for THIS K_mult summary table
            summary_data_k = [summary_data_k; curr_phi, skip_val, nvec_val, gre_all, rom_offline, rom_cpu];
            
            % Dynamic time info for legend
            time_info = '';
            if ~isnan(rom_offline), time_info = sprintf('Off: %.2fs', rom_offline); end
            if ~isnan(rom_cpu)
                if isempty(time_info), time_info = sprintf('On: %.2fs', rom_cpu);
                else, time_info = sprintf('%s, On: %.2fs', time_info, rom_cpu); end
            end
            
            % Plotting
            legend_str = sprintf('Skip=%d (%dV) | GRE: %.2f%% | %s', skip_val, nvec_val, gre_all, time_info);
            plot(ax1, t_fom, y_X_rom_node1_interp, '--', 'Color', colors(i_rom, :), 'LineWidth', 1.5, 'DisplayName', legend_str);
        end
        
        % Finalize Figure
        if length(rom_files_phi) > 6
            legend(ax1, 'Location', 'eastoutside');
        else
            legend(ax1, 'Location', 'best');
        end
        
        fig_filename = fullfile(results_dir, sprintf('TimeHistory_K%g_Phi%03d.png', k_mult, curr_phi));
        exportgraphics(fig, fig_filename, 'Resolution', 300);
        savefig(fig, fullfile(results_dir, sprintf('TimeHistory_K%g_Phi%03d.fig', k_mult, curr_phi)));
        fprintf('    -> Plot saved as %s\n', fig_filename);
    end
    
    %% --- 6. TABLE GENERATION FOR CURRENT K_MULT ---
    if ~isempty(summary_data_k)
        % Sort matrix by: Phi (Col 1), then Skip (Col 2)
        summary_data_k = sortrows(summary_data_k, [1, 2]);
        
        % Pre-allocate cell array for the GUI Figure Table
        table_cell = {};
        current_phi_tracker = -999;
        
        fprintf(log_file, '---------------------------------------------------------------------------------------\n');
        fprintf(log_file, '|  Method            | Phi Modes | Skip | MC Vecs | Global GRE [%%] | Off Time | On Time |\n');
        fprintf(log_file, '---------------------------------------------------------------------------------------\n');
        
        for i = 1:size(summary_data_k, 1)
            phi_val   = summary_data_k(i, 1);
            skip_val  = summary_data_k(i, 2);
            nvec_val  = summary_data_k(i, 3);
            gre_all   = summary_data_k(i, 4);
            off_time  = summary_data_k(i, 5);
            on_time   = summary_data_k(i, 6);
            
            % INSERT SECTION HEADERS WHEN PHI CHANGES
            if phi_val ~= current_phi_tracker
                if phi_val == -1
                    sec_title = '=== FOM REFERENCE ===';
                else
                    sec_title = sprintf('=== ROMs: PHI = %d ===', phi_val);
                end
                
                % Add separator to GUI cell array
                table_cell(end+1,:) = {sec_title, '---', '---', '---', '---', '---', '---'};
                % Add separator to Text Log
                fprintf(log_file, '| %-83s |\n', sec_title);
                
                current_phi_tracker = phi_val;
            end
            
            % Formats
            if phi_val == -1
                method_str = 'FOM';
                phi_str = 'N/A';
                skip_str = 'N/A';
                nvec_str = 'N/A';
                gre_str = '0.000';
            else
                method_str = 'ROM MC';
                phi_str = num2str(phi_val);
                skip_str = num2str(skip_val);
                nvec_str = num2str(nvec_val);
                gre_str = sprintf('%.3f', gre_all);
            end
            
            if isnan(off_time), off_str = 'N/A'; else, off_str = sprintf('%.2f', off_time); end
            if isnan(on_time), on_str = 'N/A'; else, on_str = sprintf('%.2f', on_time); end
            
            % Append to GUI cell array
            table_cell(end+1,:) = {method_str, phi_str, skip_str, nvec_str, gre_str, off_str, on_str};
            
            % Append to Text Log
            row_str = sprintf('| %-18s | %9s | %4s | %7s | %14s | %8s | %7s |\n', ...
                method_str, phi_str, skip_str, nvec_str, gre_str, off_str, on_str);
            fprintf(log_file, '%s', row_str);
        end
        fprintf(log_file, '---------------------------------------------------------------------------------------\n\n');
        
        % --- CREATE THE UI TABLE FIGURE ---
        fprintf('\nGenerating visual summary table for K_mult = %g...\n', k_mult);
        fig_table = figure('Name', sprintf('Summary Results Table - K%g', k_mult), ...
                           'Position', [150, 150, 900, 500], 'MenuBar', 'none', ...
                           'ToolBar', 'none', 'NumberTitle', 'off', 'Color', 'w');
                       
        col_names = {'Method', 'Phi Modes', 'Skip', 'MC Vecs', 'Global GRE [%]', 'Offline Time [s]', 'Online Time [s]'};
        
        t = uitable('Parent', fig_table, 'Data', table_cell, 'ColumnName', col_names, ...
                    'Units', 'normalized', 'Position', [0.02 0.02 0.96 0.96], ...
                    'RowStriping', 'on', 'FontSize', 11);
                
        % Adjust column widths (Method column is wider to fit section headers)
        t.ColumnWidth = {180, 80, 50, 70, 120, 120, 120};
        
        % Save the table as an image
        table_filename = fullfile(results_dir, sprintf('Summary_Table_K%g.png', k_mult));
        exportgraphics(fig_table, table_filename, 'Resolution', 300);
        fprintf('  -> Visual Table saved as %s\n', table_filename);
    end
end

fclose(log_file);
fprintf('\nPost-processing complete! All plots and tables are saved in %s\n', results_dir);