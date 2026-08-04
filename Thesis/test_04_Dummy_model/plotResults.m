%% =====================================================================
%  POST-PROCESSING: Compare ALL Methods (FOM vs MT, MC, Rubin, MCB, MN)
%  Figura unica: 4 facce + 1 subplot andamento GRE temporale + singola legenda
% =====================================================================
clear; close all; clc;

% --- 1. Parametri e Selettori ---
target_phi = 200; % <--- Specifica qui il numero di modi da confrontare

% Definizione delle facce, direzioni e gap associati
faces = {'T', 'B', 'L', 'R'};
dirs  = {'Y', 'Y', 'X', 'X'};
gaps  = [5e-6, -1.5e-6, -5e-6, 1.5e-6];

% Identificativi dei metodi nei nomi dei file e nomi puliti per la legenda
methods_id = {'MT_', 'MC_', 'Rubin', 'MCB', 'MN_'};
method_names = {'Modal Trunc. (MT)', 'Milman-Chu (MC)', 'Rubin', 'Massless CB (MCB)', 'MacNeal (MN)'};

% Colori per i 5 metodi
line_colors = lines(length(methods_id));

% --- 2. Selezione Cartella ---
results_dir = uigetdir(pwd, 'Select the results folder (V2)');
if results_dir == 0, error('Nessuna cartella selezionata.'); end
fprintf('Cartella selezionata: %s\n', results_dir);

% --- 3. Trova i file FOM per estrarre Q e K ---
fom_files = dir(fullfile(results_dir, 'FOM_*.mat'));
if isempty(fom_files), error('Nessun file FOM trovato nella cartella.'); end

for i_fom = 1:length(fom_files)
    fom_name = fom_files(i_fom).name;
    tokens = regexp(fom_name, 'FOM_Q(\d+)_K([\d\.]+)\.mat', 'tokens');
    if isempty(tokens), continue; end
    
    Q_val = str2double(tokens{1}{1});
    K_val = str2double(tokens{1}{2});
    fprintf('\n=================================================\n');
    fprintf('Analizzando caso: Q = %d, K_mult = %g (Phi = %d)\n', Q_val, K_val, target_phi);
    
    % Carica Dati FOM
    fom_data = load(fullfile(results_dir, fom_name));
    t_fom = fom_data.t_fom;
    
    % Costruzione della matrice FOM Globale (tutti i nodi di contatto)
    Y_fom_all = [];
    for f = 1:4
        loc = faces{f}; dir_str = dirs{f};
        if isfield(fom_data.y_contact, loc) && ~isempty(fom_data.y_contact.(loc).(dir_str))
            Y_fom_all = [Y_fom_all; fom_data.y_contact.(loc).(dir_str)];
        end
    end
    % Calcola la norma massima raggiunta dal FOM per normalizzare il GRE nel tempo
    fom_norms = vecnorm(Y_fom_all, 2, 1);
    max_fom_norm = max(fom_norms);
    
    % --- 4. Setup della Figura (Griglia 3x2) ---
    fig_name = sprintf('Compare All - Phi %d - Q%d K%g', target_phi, Q_val, K_val);
    fig = figure('Name', fig_name, 'Position', [50, 50, 1400, 950], 'Color', 'w');
    axs = zeros(1, 5);
    
    % Pre-traccia il FOM sulle prime 4 facce
    for f = 1:4
        axs(f) = subplot(3, 2, f);
        hold(axs(f), 'on'); grid(axs(f), 'on');
        
        loc = faces{f}; dir_str = dirs{f};
        
        if isfield(fom_data.y_contact, loc) && ~isempty(fom_data.y_contact.(loc).(dir_str))
            y_fom_node1 = fom_data.y_contact.(loc).(dir_str)(1, :);
            % Spegniamo l'HandleVisibility in modo che il FOM non sporchi la legenda finale unica
            plot(axs(f), t_fom, y_fom_node1, 'k-', 'LineWidth', 2.5, 'HandleVisibility', 'off');
yline(gaps(f), 'r-.', 'LineWidth', 1.5, 'Parent', axs(f), 'HandleVisibility', 'off');
        end
        
        title(axs(f), sprintf('Interface %s (%s)', loc, dir_str), 'FontSize', 11, 'FontWeight', 'bold');
        xlabel(axs(f), 'Time [s]'); ylabel(axs(f), sprintf('Disp. %s [m]', dir_str));
    end
    
    % Setup 5° Subplot: Andamento GRE Globale (Occupa l'intera terza riga)
    axs(5) = subplot(3, 2, [5, 6]);
    hold(axs(5), 'on'); grid(axs(5), 'on');
    title(axs(5), 'Time Evolution of Global Relative Error (GRE)', 'FontSize', 12, 'FontWeight', 'bold');
    xlabel(axs(5), 'Time [s]', 'FontWeight', 'bold'); 
    ylabel(axs(5), 'Instantaneous GRE [%]', 'FontWeight', 'bold');
    
    % --- 5. Estrazione e Plot dei ROM ---
    for m = 1:length(methods_id)
        pattern = sprintf('ROM_%s*Phi%03d_Q%04d_K%g*.mat', methods_id{m}, target_phi, Q_val, K_val);
        rom_files = dir(fullfile(results_dir, pattern));
        
        if isempty(rom_files)
            fprintf('  [SKIP] Nessun file %s trovato\n', method_names{m});
            continue;
        end
        
        rom_name = rom_files(1).name;
        rom_data = load(fullfile(results_dir, rom_name));
        t_rom = rom_data.t_rom;
        
        % Matrice ROM globale per questo metodo
        Y_rom_all = [];
        for f = 1:4
            loc = faces{f}; dir_str = dirs{f};
            if isfield(rom_data.y_contact, loc) && ~isempty(rom_data.y_contact.(loc).(dir_str))
                Y_rom_all = [Y_rom_all; rom_data.y_contact.(loc).(dir_str)];
            end
        end
        
        % Interpolazione sui tempi del FOM
        Y_rom_all_interp = interp1(t_rom, Y_rom_all', t_fom, 'linear', 'extrap')';
        
        % Calcolo GRE TOTALE (Scalare)
        gre_tot = norm(Y_fom_all(:) - Y_rom_all_interp(:)) / norm(Y_fom_all(:)) * 100;
        
        % Calcolo GRE NEL TEMPO (Vettore)
        diff_matrix = Y_fom_all - Y_rom_all_interp;
        diff_norms = vecnorm(diff_matrix, 2, 1);
        gre_t = (diff_norms / max_fom_norm) * 100;
        
        % Stringa dei tempi
        if isfield(rom_data, 'cpu_time'), rom_cpu = rom_data.cpu_time; else, rom_cpu = NaN; end
        if isfield(rom_data, 'offline_time'), rom_off = rom_data.offline_time; else, rom_off = NaN; end
        
        time_str = '';
        if ~isnan(rom_off), time_str = sprintf('Off: %.1fs', rom_off); end
        if ~isnan(rom_cpu)
            if isempty(time_str), time_str = sprintf('On: %.1fs', rom_cpu);
            else, time_str = sprintf('%s, On: %.1fs', time_str, rom_cpu); end
        end
        
        % Stringa COMPLETA per l'unica legenda
        legend_str = sprintf('%s | Total GRE: %6.2f%% | %s', method_names{m}, gre_tot, time_str);
        
        % 1. Plot sulle singole facce (Invisibile alla legenda)
        for f = 1:4
            loc = faces{f}; dir_str = dirs{f};
            if isfield(rom_data.y_contact, loc) && ~isempty(rom_data.y_contact.(loc).(dir_str))
                % Prende il primo nodo interpolato per tracciarlo
                y_rom_interp = interp1(t_rom, rom_data.y_contact.(loc).(dir_str)(1,:)', t_fom, 'linear', 'extrap')';
                plot(axs(f), t_fom, y_rom_interp, '--', 'LineWidth', 1.5, ...
                     'Color', line_colors(m, :), 'HandleVisibility', 'off');
            end
        end
        
        % 2. Plot del GRE nel quinto subplot (VISIBILE ALLA LEGENDA)
        plot(axs(5), t_fom, gre_t, '-', 'LineWidth', 2, ...
             'Color', line_colors(m, :), 'DisplayName', legend_str);
         
        fprintf('  [OK] %-18s | GRE Tot: %6.2f%%\n', method_names{m}, gre_tot);
    end
    
    % Aggiunta della linea nera per il FOM nella legenda (manualmente)
    if isfield(fom_data, 'cpu_time')
        fom_leg_str = sprintf('FOM Reference (On: %.1fs)', fom_data.cpu_time);
    else
        fom_leg_str = 'FOM Reference';
    end
    plot(axs(5), NaN, NaN, 'k-', 'LineWidth', 2.5, 'DisplayName', fom_leg_str);
    
    % --- 6. Finalizzazione Figura ---
    % Mette la legenda sull'ultimo subplot posizionandola in modo intelligente
    legend(axs(5), 'Location', 'eastoutside', 'FontSize', 10, 'Interpreter', 'none');
    
    sgtitle(fig, sprintf('Global Simulation Benchmarks (\\Phi = %d) | Q=%d, K_{mult}=%g', target_phi, Q_val, K_val), ...
            'FontSize', 16, 'FontWeight', 'bold');
    
    % Salvataggio
    fig_filename = fullfile(results_dir, sprintf('Compare_AllMethods_Phi%03d_Q%d_K%g.png', target_phi, Q_val, K_val));
    exportgraphics(fig, fig_filename, 'Resolution', 300);
    savefig(fig, fullfile(results_dir, sprintf('Compare_AllMethods_Phi%03d_Q%d_K%g.fig', target_phi, Q_val, K_val)));
    
    fprintf('  -> Grafico salvato come: %s\n', fig_filename);
end
fprintf('\nPost-processing completato con successo!\n');