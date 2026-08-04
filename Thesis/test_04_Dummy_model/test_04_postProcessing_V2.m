%% =====================================================================
%  POST-PROCESSING V2: FOM vs ROM (MT, MC, Rubin, MCB & MN)
%
%  Adattato agli export del main V3. Modifiche marcate [V3]:
%    - interp1 RIMOSSA: t_rom e t_fom sono ora la stessa griglia (t_common)
%      -> sostituita da un assert di coerenza (fallisce subito se un run
%         e' stato lanciato senza 'OutputTimes')
%    - aggiunto Modal Truncation ai metodi
%    - GRE calcolato su DUE finestre: intera + primi impatti
%      (oltre i primi urti il vibro-impatto diverge in modo caotico,
%       quindi il confronto puntuale perde significato)
%    - asse GRE in scala logaritmica + linea del floor di integrazione
%    - lettura di Eref / RelTol per la riproducibilita'
%    - due figure di sintesi finali: GRE vs phi e Pareto accuratezza-costo
%    - tabella riassuntiva esportata in CSV
%
%  [FLAG] gre_mode: seleziona la metrica GRE di riferimento
%    'window' -> comportamento V3 (primi impatti)
%    'full'   -> comportamento dello script vecchio (intera simulazione)
%    In entrambi i casi NESSUNA interpolazione: griglie gia' allineate.
% =====================================================================
clear; close all; clc;

%% --- Mappatura facce, direzioni e gap ---
faces = {'T', 'B', 'L', 'R'};
dirs  = {'Y', 'Y', 'X', 'X'};
gaps  = [5e-6, -1.5e-6, -5e-6, 1.5e-6];

% [V3] Modal Truncation aggiunto. MCB e MN restano in lista: se in futuro
% arrivano dal solutore dedicato vengono raccolti automaticamente.
methods_id   = {'MT_', 'MC_', 'Rubin', 'MCB', 'MN_'};
method_names = {'Modal Truncation', 'Milman-Chu', 'Rubin', 'Massless CB', 'MacNeal'};
n_methods    = numel(methods_id);

%% --- [V3] Parametri di validita' del confronto ---
% Frazione della simulazione su cui il confronto puntuale e' significativo.
% Oltre i primi impatti l'esponente di Lyapunov positivo fa divergere le
% traiettorie: l'errore satura per ragioni fisiche, non numeriche.
win_frac = 0.25;

% Floor di integrazione misurato con lo studio di convergenza in tolleranza
% (MC 3.0e-7, Rubin 2.82e-7 -> indipendente dalla base, origine nel contatto).
% Sotto questa soglia il GRE non e' piu' errore di riduzione.
gre_floor_pct = 3.0e-7 * 100;

%% --- [FLAG] Metrica GRE di riferimento -------------------------------
%  'window' -> [V3]  GRE valutato sui primi win_frac della simulazione
%  'full'   -> [old] GRE valutato sull'intera time history
%  Entrambi i valori vengono comunque sempre calcolati e loggati; il flag
%  decide quale propaga in legenda, figure di sintesi e fit di convergenza.
gre_mode = 'full';

% [FLAG] scala dell'asse Y del subplot GRE(t):
%   false -> lineare, come nello script vecchio (default)
%   true  -> logaritmica + linea del floor di integrazione
% In scala lineare il floor (%.2e %%) collassa sullo zero e la sua linea
% non e' leggibile, quindi viene disegnata solo in modalita' log.
gre_plot_log = false;

switch lower(gre_mode)
    case 'window'
        gre_tag     = 'GRE_win';        % piatto: console, log, CSV
        gre_tag_tex = 'GRE_{win}';      % TeX: legende e assi
        gre_desc    = sprintf('GRE on first %.0f%% of the simulation [%%]', 100*win_frac);
    case 'full'
        gre_tag     = 'GRE_full';
        gre_tag_tex = 'GRE_{full}';
        gre_desc    = 'GRE on the whole time history [%]';
    otherwise
        error('gre_mode non valido: usare ''window'' oppure ''full''.');
end
fprintf('Metrica GRE di riferimento: %s (modalita'' ''%s'')\n', gre_tag, gre_mode);

%% 1. Folder Selection
results_dir = uigetdir(pwd, 'Select the results folder (e.g., Shock_...)');
if results_dir == 0
    error('No folder selected. Exiting.');
end
fprintf('Selected directory: %s\n', results_dir);

%% 2. Find all FOM files
fom_files = dir(fullfile(results_dir, 'FOM_*.mat'));
if isempty(fom_files)
    error('No FOM files found in the selected directory.');
end

log_file = fopen(fullfile(results_dir, 'GRE_Results_V3.txt'), 'w');
fprintf(log_file, '======================================================\n');
fprintf(log_file, ' GLOBAL RELATIVE ERROR (GRE) & TIME REPORT - V3\n');
fprintf(log_file, '======================================================\n');
fprintf(log_file, ' Finestra di validita'': primo %.0f%% della simulazione\n', 100*win_frac);
fprintf(log_file, ' Floor di integrazione: %.2e %%\n', gre_floor_pct);
fprintf(log_file, ' Metrica di riferimento: %s (modalita'' ''%s'')\n\n', gre_tag, gre_mode);

% [V3] accumulatore per le figure di sintesi
summary = struct('method', {}, 'phi', {}, 'Q', {}, 'K', {}, ...
                 'gre_full', {}, 'gre_win', {}, 'gre_ref', {}, ...
                 'cpu', {}, 'offline', {});

%% 3. Loop over FOM files
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

    fom_data = load(fullfile(results_dir, fom_name));
    t_fom = fom_data.t_fom(:);

    % [V3] indici della finestra di validita'
    nT    = numel(t_fom);
    i_win = 1 : max(2, round(win_frac * nT));

    % [FLAG] indici effettivamente usati per la metrica di riferimento
    if strcmpi(gre_mode, 'window'), i_gre = i_win; else, i_gre = 1:nT; end

    % [V3] parametri di riproducibilita'
    if isfield(fom_data, 'Eref')
        fprintf('  Eref = %.4e J\n', fom_data.Eref);
        fprintf(log_file, '  Eref = %.4e J\n', fom_data.Eref);
    end
    if isfield(fom_data, 'fom_RelTol')
        fprintf(log_file, '  RelTol FOM = %.1e\n', fom_data.fom_RelTol);
    elseif isfield(fom_data, 'bench_RelTol')
        fprintf(log_file, '  RelTol FOM = %.1e\n', fom_data.bench_RelTol);
    end

    if isfield(fom_data, 'cpu_time')
        fom_cpu    = fom_data.cpu_time;
        fom_legend = sprintf('FOM (On: %.2fs)', fom_cpu);
    else
        fom_cpu    = NaN;
        fom_legend = 'FOM';
    end

    rom_files = dir(fullfile(results_dir, sprintf('ROM_*_Q%04d_K%g.mat', Q_val, K_val)));

    %% Loop sui metodi
    for m = 1:n_methods
        current_method_id   = methods_id{m};
        current_method_name = method_names{m};

        method_rom_files = {};
        for r = 1:length(rom_files)
            if contains(rom_files(r).name, ['ROM_' current_method_id])
                method_rom_files{end+1} = rom_files(r).name; %#ok<SAGROW>
            end
        end
        if isempty(method_rom_files), continue; end

        fig = figure('Name', sprintf('%s - Q%d - K%g', current_method_name, Q_val, K_val), ...
                     'NumberTitle', 'off', 'Position', [100, 50, 1000, 1400], 'Color', 'w');

        axs = gobjects(5,1);

        for f = 1:4
            axs(f) = subplot(5, 1, f);
            hold(axs(f), 'on'); grid(axs(f), 'on');

            face     = faces{f};
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

        axs(5) = subplot(5, 1, 5);
        hold(axs(5), 'on'); grid(axs(5), 'on');
        title(axs(5), 'Global Relative Error (GRE) Over Time');
        xlabel(axs(5), 'Time [s]');
        ylabel(axs(5), 'GRE [%]');

        % [FLAG] scala e floor: in lineare il floor e' indistinguibile da 0
        if gre_plot_log
            set(axs(5), 'YScale', 'log');
            yline(axs(5), gre_floor_pct, 'r--', 'LineWidth', 1.2, ...
                  'HandleVisibility', 'off');
        end

        % [V3] confine della finestra di validita': ha senso solo se la
        % metrica di riferimento e' quella windowed
        if strcmpi(gre_mode, 'window')
            xline(axs(5), t_fom(i_win(end)), 'k:', 'LineWidth', 1.2, ...
                  'HandleVisibility', 'off');
        end

        colors = lines(length(method_rom_files));

        for i_rom = 1:length(method_rom_files)
            rom_name = method_rom_files{i_rom};
            rom_data = load(fullfile(results_dir, rom_name));
            t_rom    = rom_data.t_rom(:);

            % ============ [V3] COERENZA DELLA GRIGLIA ============
            % Con 'OutputTimes' i due vettori sono identici: l'interpolazione
            % non serve piu'. L'assert e' la rete di sicurezza: se un run e'
            % stato lanciato senza 'OutputTimes', fallisce qui invece di
            % interpolare silenziosamente dati non confrontabili.
            if numel(t_rom) ~= numel(t_fom) || ...
               max(abs(t_rom - t_fom)) > 1e-9 * (t_fom(end) - t_fom(1))
                error(['Griglia temporale non coerente in %s.\n' ...
                       'FOM: %d punti, ROM: %d punti.\n' ...
                       'Rilanciare il run passando ''OutputTimes'', t_common.'], ...
                       rom_name, numel(t_fom), numel(t_rom));
            end
            % =====================================================

            if isfield(rom_data, 'cpu_time'),     rom_cpu = rom_data.cpu_time;         else, rom_cpu = NaN; end
            if isfield(rom_data, 'offline_time'), rom_offline = rom_data.offline_time; else, rom_offline = NaN; end

            phi_tokens = regexp(rom_name, 'Phi(\d+)', 'tokens');
            phi_val    = str2double(phi_tokens{1}{1});

            % --- concatenazione di tutti i nodi di contatto ---
            y_fom_cat = []; y_rom_cat = [];
            for f_idx = 1:4
                fc = faces{f_idx};
                dc = dirs{f_idx};
                if isfield(fom_data.y_contact, fc) && isfield(rom_data.y_contact, fc)
                    y_fom_cat = [y_fom_cat; fom_data.y_contact.(fc).(dc)]; %#ok<AGROW>
                    y_rom_cat = [y_rom_cat; rom_data.y_contact.(fc).(dc)]; %#ok<AGROW>
                    % [V3] nessuna interpolazione: stessa griglia
                end
            end

            % --- GRE globale, intera simulazione ---
            gre_all = (norm(y_fom_cat - y_rom_cat, 'fro') / ...
                       norm(y_fom_cat, 'fro')) * 100;

            % --- [V3] GRE sulla finestra di validita' (primi impatti) ---
            gre_win = (norm(y_fom_cat(:,i_win) - y_rom_cat(:,i_win), 'fro') / ...
                       norm(y_fom_cat(:,i_win), 'fro')) * 100;

            % --- [FLAG] metrica di riferimento selezionata ---
            gre_ref = (norm(y_fom_cat(:,i_gre) - y_rom_cat(:,i_gre), 'fro') / ...
                       norm(y_fom_cat(:,i_gre), 'fro')) * 100;

            % --- GRE nel tempo (istantaneo normalizzato) ---
            % NOTA: la normalizzazione resta sul massimo dell'INTERA storia in
            % entrambe le modalita', cosi' la curva plottata e' invariante
            % rispetto al flag e resta identica a quella dello script vecchio.
            norm_diff_t    = sqrt(sum((y_fom_cat - y_rom_cat).^2, 1));
            max_fom_global = max(sqrt(sum(y_fom_cat.^2, 1)));
            gre_t          = (norm_diff_t ./ (max_fom_global + eps)) * 100;

            % --- legenda ---
            time_info = '';
            if ~isnan(rom_offline), time_info = sprintf('Off: %.2fs', rom_offline); end
            if ~isnan(rom_cpu)
                if isempty(time_info), time_info = sprintf('On: %.2fs', rom_cpu);
                else,                  time_info = sprintf('%s, On: %.2fs', time_info, rom_cpu); end
            end

            % [FLAG] la legenda riporta la metrica selezionata; l'altra
            % finisce comunque nel log e nel CSV.
            if isempty(time_info)
                legend_str = sprintf('ROM \\phi=%d (%s: %.3f%%)', ...
                                     phi_val, gre_tag_tex, gre_ref);
            else
                legend_str = sprintf('ROM \\phi=%d (%s: %.3f%%, %s)', ...
                                     phi_val, gre_tag_tex, gre_ref, time_info);
            end

            for f = 1:4
                face     = faces{f};
                dir_plot = dirs{f};
                if isfield(rom_data.y_contact, face)
                    y_rom_face = rom_data.y_contact.(face).(dir_plot)(1, :);
                    if f == 1
                        plot(axs(f), t_fom, y_rom_face(:), '--', 'Color', colors(i_rom, :), ...
                             'LineWidth', 1.5, 'DisplayName', legend_str);
                    else
                        plot(axs(f), t_fom, y_rom_face(:), '--', 'Color', colors(i_rom, :), ...
                             'LineWidth', 1.5, 'HandleVisibility', 'off');
                    end
                end
            end

            plot(axs(5), t_fom, gre_t(:), '-', 'Color', colors(i_rom, :), ...
                 'LineWidth', 1.5, 'HandleVisibility', 'off');

            % --- [V3] accumulo per le figure di sintesi ---
            summary(end+1) = struct('method', current_method_name, ...
                                    'phi', phi_val, 'Q', Q_val, 'K', K_val, ...
                                    'gre_full', gre_all, 'gre_win', gre_win, ...
                                    'gre_ref', gre_ref, ...
                                    'cpu', rom_cpu, 'offline', rom_offline); %#ok<SAGROW>

            % --- logging ---
            if isnan(rom_offline), off_str = 'N/A'; else, off_str = sprintf('%6.2fs', rom_offline); end
            if isnan(rom_cpu),      on_str = 'N/A'; else,  on_str = sprintf('%6.2fs', rom_cpu);     end

            % [V3] segnala se il GRE si avvicina al floor di integrazione
            if gre_ref < 100 * gre_floor_pct
                flag = '  [!] vicino al floor di integrazione';
            else
                flag = '';
            end

            fprintf('  [%-16s] Phi: %3d | %s: %8.4f%% | GRE_full: %8.4f%% | Off: %s | On: %s%s\n', ...
                    current_method_name, phi_val, gre_tag, gre_ref, gre_all, off_str, on_str, flag);
            fprintf(log_file, '  %-16s Phi: %03d | %s: %8.4f%% | GRE_full: %8.4f%% | Off: %s | On: %s%s\n', ...
                    current_method_name, phi_val, gre_tag, gre_ref, gre_all, off_str, on_str, flag);
        end

        legend(axs(1), 'Location', 'best');
        sgtitle(fig, sprintf('%s Method (Q=%d, K_{mult}=%g)', ...
                current_method_name, Q_val, K_val), 'FontSize', 16, 'FontWeight', 'bold');

        clean_method_name = strrep(current_method_name, ' ', '');
        fig_filename = fullfile(results_dir, sprintf('Compare_%s_Q%d_K%g.png', clean_method_name, Q_val, K_val));
        exportgraphics(fig, fig_filename, 'Resolution', 300);
        savefig(fig, fullfile(results_dir, sprintf('Compare_%s_Q%d_K%g.fig', clean_method_name, Q_val, K_val)));

        fprintf('  -> Plot saved as %s\n', fig_filename);
    end
    fprintf(log_file, '\n');
end

fclose(log_file);

%% =====================================================================
%  [V3] FIGURE DI SINTESI - le due che servono in tesi
% =====================================================================
if isempty(summary)
    fprintf('\nNessun ROM trovato: figure di sintesi saltate.\n');
    return
end

T_summary = struct2table(summary);
writetable(T_summary, fullfile(results_dir, ...
           sprintf('summary_V3_%s.csv', lower(gre_mode))));

uniq_methods = unique(T_summary.method, 'stable');
mk = {'o-','s-','^-','d-','v-'};

% --- Figura A: convergenza GRE vs phi ---------------------------------
figA = figure('Name','Convergenza in phi','Color','w','Position',[100 100 800 600]);
hold on; grid on;
for m = 1:numel(uniq_methods)
    sel = strcmp(T_summary.method, uniq_methods{m});
    [phis, iord] = sort(T_summary.phi(sel));
    g = T_summary.gre_ref(sel);
    loglog(phis, g(iord), mk{min(m,numel(mk))}, 'LineWidth', 1.8, ...
           'MarkerSize', 7, 'DisplayName', uniq_methods{m});
end
yline(gre_floor_pct, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Integration floor');
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Number of retained modes \phi');
ylabel(gre_desc);
title('ROM convergence: reduction error vs basis size');
legend('Location','southwest'); box on;
exportgraphics(figA, fullfile(results_dir, ...
    sprintf('Summary_GRE_vs_phi_%s.png', lower(gre_mode))), 'Resolution', 300);
savefig(figA, fullfile(results_dir, ...
    sprintf('Summary_GRE_vs_phi_%s.fig', lower(gre_mode))));

% --- Figura B: Pareto accuratezza-costo -------------------------------
figB = figure('Name','Pareto accuratezza-costo','Color','w','Position',[100 100 800 600]);
hold on; grid on;
for m = 1:numel(uniq_methods)
    sel = strcmp(T_summary.method, uniq_methods{m});
    [cpus, iord] = sort(T_summary.cpu(sel));
    g = T_summary.gre_ref(sel);
    loglog(cpus, g(iord), mk{min(m,numel(mk))}, 'LineWidth', 1.8, ...
           'MarkerSize', 7, 'DisplayName', uniq_methods{m});
end
yline(gre_floor_pct, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Integration floor');
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Online CPU time [s]');
ylabel(gre_desc);
title('Accuracy-cost trade-off (online)');
legend('Location','southwest'); box on;
exportgraphics(figB, fullfile(results_dir, ...
    sprintf('Summary_Pareto_%s.png', lower(gre_mode))), 'Resolution', 300);
savefig(figB, fullfile(results_dir, ...
    sprintf('Summary_Pareto_%s.fig', lower(gre_mode))));

%% --- [V3] Ordine di convergenza osservato ----------------------------
fprintf('\n--- Ordine di convergenza osservato (%s ~ phi^-p) ---\n', gre_tag);
for m = 1:numel(uniq_methods)
    sel  = strcmp(T_summary.method, uniq_methods{m});
    phis = T_summary.phi(sel);
    g    = T_summary.gre_ref(sel);
    ok   = phis > 0 & g > 0;
    if nnz(ok) >= 2
        pfit = polyfit(log(phis(ok)), log(g(ok)), 1);
        fprintf('  %-16s : p = %.2f   (%d punti)\n', uniq_methods{m}, -pfit(1), nnz(ok));
    end
end

fprintf('\nPost-processing complete! Results are saved in %s\n', results_dir);