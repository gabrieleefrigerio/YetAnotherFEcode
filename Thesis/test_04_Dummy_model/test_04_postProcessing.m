%% =====================================================================
%  POST-PROCESSING - FOM vs ROM, confronto e convergenza
%
%  Main unico, adattivo nel numero di interfacce di contatto: facce,
%  direzioni e gap NON sono scritti qui, ma letti da run_config.mat nella
%  cartella dei risultati. Il post-processing si adatta quindi da solo sia
%  al modello a quattro interfacce sia a quello a una sola, e non puo'
%  andare fuori sincrono con i gap effettivamente usati nella simulazione.
%
%  Prodotti:
%    - una figura per metodo, con un subplot per interfaccia + GRE(t)
%    - GRE_Results.txt      log testuale
%    - summary_<mode>.csv   tabella riassuntiva
%    - Summary_GRE_vs_phi   convergenza dell'errore al crescere dei modi
%    - Summary_Pareto       compromesso accuratezza / costo online
%
%  Nessuna interpolazione: tutti i modelli sono salvati sulla stessa griglia
%  t_common. Se un file non la rispetta lo script si ferma, invece di
%  confrontare silenziosamente dati non confrontabili.
% =====================================================================
clear; close all; clc;

%% --- Opzioni ----------------------------------------------------------

% Metrica GRE di riferimento (quella che finisce in legenda e nelle sintesi):
%   'full'   -> intera time history
%   'window' -> primi win_frac della simulazione
% Oltre i primi impatti l'esponente di Lyapunov positivo fa divergere le
% traiettorie: l'errore satura per ragioni fisiche, non numeriche. Entrambe
% le metriche vengono comunque sempre calcolate e loggate.
gre_mode = 'full';
win_frac = 0.25;

% Floor di integrazione misurato con lo studio di convergenza in tolleranza
% (MC 3.0e-7, Rubin 2.82e-7 -> indipendente dalla base, origine nel contatto).
% Sotto questa soglia il GRE non e' piu' errore di riduzione.
gre_floor_pct = 3.0e-7 * 100;

% Scala dell'asse Y del subplot GRE(t). In lineare il floor collassa sullo
% zero e la sua linea non e' leggibile, quindi viene disegnata solo in log.
gre_plot_log = false;

% Nodo di ciascuna interfaccia da disegnare nelle time history.
node_idx = 1;

switch lower(gre_mode)
    case 'window'
        gre_tag = 'GRE_win';  gre_tag_tex = 'GRE_{win}';
        gre_desc = sprintf('GRE on first %.0f%% of the simulation [%%]', 100*win_frac);
    case 'full'
        gre_tag = 'GRE_full'; gre_tag_tex = 'GRE_{full}';
        gre_desc = 'GRE on the whole time history [%]';
    otherwise
        error('gre_mode non valido: usare ''full'' oppure ''window''.');
end
fprintf('Metrica di riferimento: %s\n', gre_tag);

%% --- 1. Cartella dei risultati ----------------------------------------
results_dir = uigetdir(pwd, 'Selezionare la cartella dei risultati');
if results_dir == 0
    error('Nessuna cartella selezionata.');
end
fprintf('Cartella: %s\n', results_dir);

cfg_file = fullfile(results_dir, 'run_config.mat');
if ~exist(cfg_file, 'file')
    error('PP:NoConfig', ...
        ['run_config.mat non trovato in\n  %s\n' ...
         'Questa cartella e'' stata prodotta da una versione precedente del main, ' ...
         'che non salvava la configurazione. Rilanciare la simulazione con test_04_main.m.'], ...
        results_dir);
end
R = load(cfg_file);

%% --- 2. Tabella delle interfacce (letta dal file, non riscritta a mano) ---
faces    = R.active_labels;
n_faces  = numel(faces);
face_dir = cell(1, n_faces);
face_gap = zeros(1, n_faces);
for i = 1:n_faces
    if R.Interfaces.(faces{i}).dir == 1
        face_dir{i} = 'X';
    else
        face_dir{i} = 'Y';
    end
    face_gap(i) = R.Interfaces.(faces{i}).gap;
end

fprintf('Interfacce: %d\n', n_faces);
for i = 1:n_faces
    fprintf('  %-4s dir %s  gap %+9.3e m  (%d nodi)\n', ...
        faces{i}, face_dir{i}, face_gap(i), numel(R.Interfaces.(faces{i}).nodes));
end
fprintf('Eref = %.4e J\n', R.Eref);

%% --- 3. Log e accumulatore --------------------------------------------
fom_files = dir(fullfile(results_dir, 'FOM_*.mat'));
if isempty(fom_files)
    error('PP:NoFOM', 'Nessun file FOM_*.mat: serve un riferimento per il confronto.');
end

% NB: chiusura esplicita in fondo allo script. onCleanup qui non servirebbe:
% in uno script le variabili restano nel workspace base a fine esecuzione,
% quindi l'oggetto non viene distrutto e il file resterebbe aperto.
log_file = fopen(fullfile(results_dir, 'GRE_Results.txt'), 'w');
fprintf(log_file, '======================================================\n');
fprintf(log_file, ' GLOBAL RELATIVE ERROR (GRE) & TIME REPORT\n');
fprintf(log_file, '======================================================\n');
fprintf(log_file, ' Interfacce: %s\n', strjoin(faces, ', '));
fprintf(log_file, ' Metrica di riferimento: %s (finestra %.0f%%)\n', gre_tag, 100*win_frac);
fprintf(log_file, ' Floor di integrazione: %.2e %%\n', gre_floor_pct);
fprintf(log_file, ' Eref = %.4e J\n\n', R.Eref);

summary = struct('method', {}, 'phi', {}, 'Q', {}, 'K', {}, ...
                 'gre_full', {}, 'gre_win', {}, 'gre_ref', {}, ...
                 'cpu', {}, 'offline', {});

%% --- 4. Ciclo sui casi (Q, k_mult) ------------------------------------
for i_fom = 1:numel(fom_files)
    tok = regexp(fom_files(i_fom).name, 'FOM_Q(\d+)_K([\d\.]+)\.mat', 'tokens');
    if isempty(tok), continue; end
    Q_val = str2double(tok{1}{1});
    K_val = str2double(tok{1}{2});

    fprintf('\n======================================================\n');
    fprintf('Caso: Q = %d | k_mult = %g\n', Q_val, K_val);
    fprintf('======================================================\n');
    fprintf(log_file, '>>> CASO: Q = %d | k_mult = %g <<<\n', Q_val, K_val);

    fom   = load(fullfile(results_dir, fom_files(i_fom).name));
    t_ref = fom.t(:);
    nT    = numel(t_ref);
    i_win = 1 : max(2, round(win_frac * nT));
    if strcmpi(gre_mode, 'window'), i_gre = i_win; else, i_gre = 1:nT; end

    if isfield(fom, 'cpu_time')
        fom_legend = sprintf('FOM (On: %.2fs)', fom.cpu_time);
    else
        fom_legend = 'FOM';
    end

    % --- ROM disponibili per questo caso, raggruppati per metodo ---
    rom_files = dir(fullfile(results_dir, sprintf('ROM_*_Q%04d_K%g.mat', Q_val, K_val)));
    if isempty(rom_files)
        fprintf('  Nessun ROM per questo caso.\n');
        continue;
    end
    rom_models = cell(1, numel(rom_files));
    for r = 1:numel(rom_files)
        mt = regexp(rom_files(r).name, '^ROM_([A-Za-z]+)_Phi', 'tokens', 'once');
        rom_models{r} = mt{1};
    end
    methods_here = unique(rom_models, 'stable');

    %% --- Ciclo sui metodi ---
    for im = 1:numel(methods_here)
        method    = methods_here{im};
        sel_files = rom_files(strcmp(rom_models, method));

        fig = figure('Name', sprintf('%s - Q%d - K%g', method, Q_val, K_val), ...
                     'NumberTitle', 'off', 'Color', 'w', ...
                     'Position', [100, 50, 1000, 250*(n_faces+1)]);
        n_sub = n_faces + 1;
        axs   = gobjects(n_sub, 1);

        % --- un subplot per interfaccia: risposta del FOM e posizione del muro ---
        for f = 1:n_faces
            axs(f) = subplot(n_sub, 1, f);
            hold(axs(f), 'on'); grid(axs(f), 'on');

            y_fom_face = fom.y_contact.(faces{f}).(face_dir{f})(node_idx, :);
            if f == 1
                plot(axs(f), t_ref, y_fom_face(:), 'k-', 'LineWidth', 2, ...
                    'DisplayName', fom_legend);
                yline(axs(f), face_gap(f), 'r-.', 'LineWidth', 1.5, ...
                    'DisplayName', 'Wall gap');
            else
                plot(axs(f), t_ref, y_fom_face(:), 'k-', 'LineWidth', 2, ...
                    'HandleVisibility', 'off');
                yline(axs(f), face_gap(f), 'r-.', 'LineWidth', 1.5, ...
                    'HandleVisibility', 'off');
            end
            title(axs(f), sprintf('Interfaccia %s (dir %s, gap %+.2e m) - nodo %d', ...
                faces{f}, face_dir{f}, face_gap(f), node_idx));
            xlabel(axs(f), 'Tempo [s]');
            ylabel(axs(f), 'Spostamento [m]');
        end

        % --- subplot GRE(t) ---
        axs(n_sub) = subplot(n_sub, 1, n_sub);
        hold(axs(n_sub), 'on'); grid(axs(n_sub), 'on');
        title(axs(n_sub), 'Global Relative Error nel tempo');
        xlabel(axs(n_sub), 'Tempo [s]');
        ylabel(axs(n_sub), 'GRE [%]');
        if gre_plot_log
            set(axs(n_sub), 'YScale', 'log');
            yline(axs(n_sub), gre_floor_pct, 'r--', 'LineWidth', 1.2, ...
                'HandleVisibility', 'off');
        end
        if strcmpi(gre_mode, 'window')
            xline(axs(n_sub), t_ref(i_win(end)), 'k:', 'LineWidth', 1.2, ...
                'HandleVisibility', 'off');
        end

        colors = lines(numel(sel_files));

        for i_rom = 1:numel(sel_files)
            rom   = load(fullfile(results_dir, sel_files(i_rom).name));
            t_rom = rom.t(:);

            % Griglia comune: nessuna interpolazione. Se un run non la
            % rispetta ci si ferma qui invece di confrontare dati diversi.
            if numel(t_rom) ~= nT || max(abs(t_rom - t_ref)) > 1e-9*(t_ref(end)-t_ref(1))
                error('PP:GridMismatch', ...
                    ['Griglia temporale incoerente in %s.\n' ...
                     'FOM: %d punti, ROM: %d punti.\n' ...
                     'Rilanciare quel run sulla griglia comune t_common.'], ...
                    sel_files(i_rom).name, nT, numel(t_rom));
            end

            if isfield(rom, 'n_modes')
                phi_val = rom.n_modes;
            else
                pt = regexp(sel_files(i_rom).name, 'Phi(\d+)', 'tokens', 'once');
                phi_val = str2double(pt{1});
            end
            if isfield(rom, 'cpu_time'),     rom_cpu = rom.cpu_time;     else, rom_cpu = NaN; end
            if isfield(rom, 'offline_time'), rom_off = rom.offline_time; else, rom_off = NaN; end

            % --- tutti i GdL di contatto concatenati, su tutte le interfacce ---
            y_fom_cat = [];
            y_rom_cat = [];
            for f = 1:n_faces
                y_fom_cat = [y_fom_cat; fom.y_contact.(faces{f}).(face_dir{f})]; %#ok<AGROW>
                y_rom_cat = [y_rom_cat; rom.y_contact.(faces{f}).(face_dir{f})]; %#ok<AGROW>
            end

            gre_full = norm(y_fom_cat - y_rom_cat, 'fro') / norm(y_fom_cat, 'fro') * 100;
            gre_win  = norm(y_fom_cat(:,i_win) - y_rom_cat(:,i_win), 'fro') / ...
                       norm(y_fom_cat(:,i_win), 'fro') * 100;
            gre_ref  = norm(y_fom_cat(:,i_gre) - y_rom_cat(:,i_gre), 'fro') / ...
                       norm(y_fom_cat(:,i_gre), 'fro') * 100;

            % GRE istantaneo, normalizzato sul massimo dell'intera storia,
            % cosi' la curva disegnata non dipende dal flag gre_mode.
            norm_diff_t = sqrt(sum((y_fom_cat - y_rom_cat).^2, 1));
            gre_t = norm_diff_t ./ (max(sqrt(sum(y_fom_cat.^2, 1))) + eps) * 100;

            % --- legenda ---
            time_info = '';
            if ~isnan(rom_off), time_info = sprintf('Off: %.2fs', rom_off); end
            if ~isnan(rom_cpu)
                if isempty(time_info)
                    time_info = sprintf('On: %.2fs', rom_cpu);
                else
                    time_info = sprintf('%s, On: %.2fs', time_info, rom_cpu);
                end
            end
            if isempty(time_info)
                legend_str = sprintf('ROM \\phi=%d (%s: %.3f%%)', phi_val, gre_tag_tex, gre_ref);
            else
                legend_str = sprintf('ROM \\phi=%d (%s: %.3f%%, %s)', ...
                    phi_val, gre_tag_tex, gre_ref, time_info);
            end

            for f = 1:n_faces
                y_rom_face = rom.y_contact.(faces{f}).(face_dir{f})(node_idx, :);
                if f == 1
                    plot(axs(f), t_ref, y_rom_face(:), '--', 'Color', colors(i_rom,:), ...
                        'LineWidth', 1.5, 'DisplayName', legend_str);
                else
                    plot(axs(f), t_ref, y_rom_face(:), '--', 'Color', colors(i_rom,:), ...
                        'LineWidth', 1.5, 'HandleVisibility', 'off');
                end
            end
            plot(axs(n_sub), t_ref, gre_t(:), '-', 'Color', colors(i_rom,:), ...
                'LineWidth', 1.5, 'HandleVisibility', 'off');

            summary(end+1) = struct('method', method, 'phi', phi_val, ...
                'Q', Q_val, 'K', K_val, 'gre_full', gre_full, 'gre_win', gre_win, ...
                'gre_ref', gre_ref, 'cpu', rom_cpu, 'offline', rom_off); %#ok<SAGROW>

            % --- log ---
            if isnan(rom_off), off_str = 'N/A'; else, off_str = sprintf('%6.2fs', rom_off); end
            if isnan(rom_cpu), on_str  = 'N/A'; else, on_str  = sprintf('%6.2fs', rom_cpu); end
            if gre_ref < 100 * gre_floor_pct
                flag = '  [!] vicino al floor di integrazione';
            else
                flag = '';
            end
            % Entrambe le metriche sempre in chiaro; quale sia il riferimento
            % e' scritto nell'intestazione ed e' quella usata nelle sintesi.
            fprintf('  [%-6s] Phi %3d | GRE_full %9.4f%% | GRE_win %9.4f%% | Off %s | On %s%s\n', ...
                method, phi_val, gre_full, gre_win, off_str, on_str, flag);
            fprintf(log_file, '  %-6s Phi %03d | GRE_full %9.4f%% | GRE_win %9.4f%% | Off %s | On %s%s\n', ...
                method, phi_val, gre_full, gre_win, off_str, on_str, flag);
        end

        legend(axs(1), 'Location', 'best');
        sgtitle(fig, sprintf('Metodo %s (Q = %d, k_{mult} = %g)', method, Q_val, K_val), ...
            'FontSize', 16, 'FontWeight', 'bold');

        base_name = fullfile(results_dir, sprintf('Compare_%s_Q%d_K%g', method, Q_val, K_val));
        exportgraphics(fig, [base_name '.png'], 'Resolution', 300);
        savefig(fig, [base_name '.fig']);
        fprintf('  -> figura salvata: %s.png\n', base_name);
    end
    fprintf(log_file, '\n');
end

%% --- 5. Figure di sintesi ---------------------------------------------
if isempty(summary)
    fprintf('\nNessun ROM confrontato: figure di sintesi saltate.\n');
    fclose(log_file);
    return
end

T_summary = struct2table(summary);
writetable(T_summary, fullfile(results_dir, sprintf('summary_%s.csv', lower(gre_mode))));

uniq_methods = unique(T_summary.method, 'stable');
mk = {'o-','s-','^-','d-','v-','>-'};

% --- Figura A: convergenza del GRE al crescere dei modi ---
figA = figure('Name','Convergenza in phi','Color','w','Position',[100 100 800 600]);
hold on; grid on;
for m = 1:numel(uniq_methods)
    sel = strcmp(T_summary.method, uniq_methods{m});
    [phis, iord] = sort(T_summary.phi(sel));
    g = T_summary.gre_ref(sel);
    plot(phis, g(iord), mk{min(m,numel(mk))}, 'LineWidth', 1.8, ...
        'MarkerSize', 7, 'DisplayName', uniq_methods{m});
end
yline(gre_floor_pct, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Integration floor');
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Numero di modi ritenuti \phi');
ylabel(gre_desc);
title('Convergenza del ROM: errore di riduzione vs dimensione della base');
legend('Location','southwest'); box on;
exportgraphics(figA, fullfile(results_dir, ...
    sprintf('Summary_GRE_vs_phi_%s.png', lower(gre_mode))), 'Resolution', 300);
savefig(figA, fullfile(results_dir, sprintf('Summary_GRE_vs_phi_%s.fig', lower(gre_mode))));

% --- Figura B: compromesso accuratezza / costo online ---
figB = figure('Name','Pareto accuratezza-costo','Color','w','Position',[100 100 800 600]);
hold on; grid on;
for m = 1:numel(uniq_methods)
    sel = strcmp(T_summary.method, uniq_methods{m});
    [cpus, iord] = sort(T_summary.cpu(sel));
    g = T_summary.gre_ref(sel);
    plot(cpus, g(iord), mk{min(m,numel(mk))}, 'LineWidth', 1.8, ...
        'MarkerSize', 7, 'DisplayName', uniq_methods{m});
end
yline(gre_floor_pct, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Integration floor');
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Tempo CPU online [s]');
ylabel(gre_desc);
title('Compromesso accuratezza-costo (online)');
legend('Location','southwest'); box on;
exportgraphics(figB, fullfile(results_dir, ...
    sprintf('Summary_Pareto_%s.png', lower(gre_mode))), 'Resolution', 300);
savefig(figB, fullfile(results_dir, sprintf('Summary_Pareto_%s.fig', lower(gre_mode))));

%% --- 6. Ordine di convergenza osservato -------------------------------
fprintf('\n--- Ordine di convergenza osservato (%s ~ phi^-p) ---\n', gre_tag);
fprintf(log_file, '\n--- Ordine di convergenza osservato (%s ~ phi^-p) ---\n', gre_tag);
for m = 1:numel(uniq_methods)
    sel  = strcmp(T_summary.method, uniq_methods{m});
    phis = T_summary.phi(sel);
    g    = T_summary.gre_ref(sel);
    ok   = phis > 0 & g > 0;
    if nnz(ok) >= 2
        pfit = polyfit(log(phis(ok)), log(g(ok)), 1);
        fprintf('  %-8s : p = %.2f   (%d punti)\n', uniq_methods{m}, -pfit(1), nnz(ok));
        fprintf(log_file, '  %-8s : p = %.2f   (%d punti)\n', uniq_methods{m}, -pfit(1), nnz(ok));
    end
end

fclose(log_file);
fprintf('\nPost-processing completato. Risultati in %s\n', results_dir);
