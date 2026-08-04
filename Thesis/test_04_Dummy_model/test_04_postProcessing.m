%% =====================================================================
%  POST-PROCESSING - FOM vs ROM comparison and convergence
%
%  Single main, adaptive in the number of contact interfaces: faces,
%  directions and gaps are NOT written here but read from run_config.mat in
%  the results folder. The post-processing therefore adapts by itself to the
%  four-interface model as well as to a single-interface one, and cannot
%  drift out of sync with the gaps actually used in the simulation.
%
%  Outputs:
%    - one figure per method, with a subplot per interface plus GRE(t)
%    - GRE_Results.txt      text log
%    - summary_<mode>.csv   summary table
%    - Summary_GRE_vs_phi   error convergence against the number of modes
%    - Summary_Pareto       accuracy / online cost trade-off
%
%  No interpolation: every model is stored on the same t_common grid. If a
%  file does not respect it the script stops, rather than silently comparing
%  data that is not comparable.
% =====================================================================
clear; close all; clc;

%% --- Options ----------------------------------------------------------

% Reference GRE metric, the one shown in the legends and the summaries:
%   'full'   -> whole time history
%   'window' -> first win_frac of the simulation
% Past the first impacts the positive Lyapunov exponent makes the
% trajectories diverge: the error saturates for physical reasons, not
% numerical ones. Both metrics are always computed and logged anyway.
gre_mode = 'full';
win_frac = 0.25;

% Integration floor measured with the tolerance convergence study
% (MC 3.0e-7, Rubin 2.82e-7 -> independent of the basis, originating in the
% contact). Below this threshold the GRE is no longer a reduction error.
gre_floor_pct = 3.0e-7 * 100;

% Y-axis scale of the GRE(t) subplot. On a linear scale the floor collapses
% onto zero and its line is unreadable, so it is drawn only in log scale.
gre_plot_log = false;

% Node of each interface to draw in the time histories.
node_idx = 1;

switch lower(gre_mode)
    case 'window'
        gre_tag = 'GRE_win';  gre_tag_tex = 'GRE_{win}';
        gre_desc = sprintf('GRE on first %.0f%% of the simulation [%%]', 100*win_frac);
    case 'full'
        gre_tag = 'GRE_full'; gre_tag_tex = 'GRE_{full}';
        gre_desc = 'GRE on the whole time history [%]';
    otherwise
        error('PP:BadMode', 'Invalid gre_mode: use ''full'' or ''window''.');
end
fprintf('Reference metric: %s\n', gre_tag);

%% --- 1. Results folder ------------------------------------------------
results_dir = uigetdir(pwd, 'Select the results folder');
if results_dir == 0
    error('PP:NoFolder', 'No folder selected.');
end
fprintf('Folder: %s\n', results_dir);

cfg_file = fullfile(results_dir, 'run_config.mat');
if ~exist(cfg_file, 'file')
    error('PP:NoConfig', ...
        ['run_config.mat not found in\n  %s\n' ...
         'This folder was produced by an earlier version of the main, which did ' ...
         'not save the configuration. Re-run the simulation with test_04_main.m.'], ...
        results_dir);
end
R = load(cfg_file);

%% --- 2. Interface table (read from the file, not rewritten by hand) ----
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

fprintf('Interfaces: %d\n', n_faces);
for i = 1:n_faces
    fprintf('  %-4s dir %s  gap %+9.3e m  (%d nodes)\n', ...
        faces{i}, face_dir{i}, face_gap(i), numel(R.Interfaces.(faces{i}).nodes));
end
fprintf('Eref = %.4e J\n', R.Eref);

%% --- 3. Log and accumulator -------------------------------------------
fom_files = dir(fullfile(results_dir, 'FOM_*.mat'));
if isempty(fom_files)
    error('PP:NoFOM', 'No FOM_*.mat file: a reference is needed for the comparison.');
end

% Note: closed explicitly at the end of the script. onCleanup would not help
% here, because in a script the variables stay in the base workspace after
% execution, so the object is never destroyed and the file would stay open.
log_file = fopen(fullfile(results_dir, 'GRE_Results.txt'), 'w');
fprintf(log_file, '======================================================\n');
fprintf(log_file, ' GLOBAL RELATIVE ERROR (GRE) & TIME REPORT\n');
fprintf(log_file, '======================================================\n');
fprintf(log_file, ' Interfaces: %s\n', strjoin(faces, ', '));
fprintf(log_file, ' Reference metric: %s (window %.0f%%)\n', gre_tag, 100*win_frac);
fprintf(log_file, ' Integration floor: %.2e %%\n', gre_floor_pct);
fprintf(log_file, ' Eref = %.4e J\n\n', R.Eref);

summary = struct('method', {}, 'phi', {}, 'Q', {}, 'K', {}, ...
                 'gre_full', {}, 'gre_win', {}, 'gre_ref', {}, ...
                 'cpu', {}, 'offline', {});

%% --- 4. Loop over the cases (Q, k_mult) -------------------------------
for i_fom = 1:numel(fom_files)
    tok = regexp(fom_files(i_fom).name, 'FOM_Q(\d+)_K([\d\.]+)\.mat', 'tokens');
    if isempty(tok), continue; end
    Q_val = str2double(tok{1}{1});
    K_val = str2double(tok{1}{2});

    fprintf('\n======================================================\n');
    fprintf('Case: Q = %d | k_mult = %g\n', Q_val, K_val);
    fprintf('======================================================\n');
    fprintf(log_file, '>>> CASE: Q = %d | k_mult = %g <<<\n', Q_val, K_val);

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

    % --- ROMs available for this case, grouped by method ---
    rom_files = dir(fullfile(results_dir, sprintf('ROM_*_Q%04d_K%g.mat', Q_val, K_val)));
    if isempty(rom_files)
        fprintf('  No ROM for this case.\n');
        continue;
    end
    rom_models = cell(1, numel(rom_files));
    for r = 1:numel(rom_files)
        mt = regexp(rom_files(r).name, '^ROM_([A-Za-z]+)_Phi', 'tokens', 'once');
        rom_models{r} = mt{1};
    end
    methods_here = unique(rom_models, 'stable');

    %% --- Loop over the methods ---
    for im = 1:numel(methods_here)
        method    = methods_here{im};
        sel_files = rom_files(strcmp(rom_models, method));

        fig = figure('Name', sprintf('%s - Q%d - K%g', method, Q_val, K_val), ...
                     'NumberTitle', 'off', 'Color', 'w', ...
                     'Position', [100, 50, 1000, 250*(n_faces+1)]);
        n_sub = n_faces + 1;
        axs   = gobjects(n_sub, 1);

        % --- one subplot per interface: FOM response and wall position ---
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
            title(axs(f), sprintf('Interface %s (dir %s, gap %+.2e m) - node %d', ...
                faces{f}, face_dir{f}, face_gap(f), node_idx));
            xlabel(axs(f), 'Time [s]');
            ylabel(axs(f), 'Displacement [m]');
        end

        % --- GRE(t) subplot ---
        axs(n_sub) = subplot(n_sub, 1, n_sub);
        hold(axs(n_sub), 'on'); grid(axs(n_sub), 'on');
        title(axs(n_sub), 'Global Relative Error over time');
        xlabel(axs(n_sub), 'Time [s]');
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

            % Common grid, no interpolation. If a run does not respect it we
            % stop here rather than comparing different things.
            if numel(t_rom) ~= nT || max(abs(t_rom - t_ref)) > 1e-9*(t_ref(end)-t_ref(1))
                error('PP:GridMismatch', ...
                    ['Inconsistent time grid in %s.\n' ...
                     'FOM: %d points, ROM: %d points.\n' ...
                     'Re-run that case on the common grid t_common.'], ...
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

            % --- all contact DOFs concatenated, across every interface ---
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

            % Instantaneous GRE, normalized on the maximum of the whole
            % history, so that the plotted curve does not depend on gre_mode.
            norm_diff_t = sqrt(sum((y_fom_cat - y_rom_cat).^2, 1));
            gre_t = norm_diff_t ./ (max(sqrt(sum(y_fom_cat.^2, 1))) + eps) * 100;

            % --- legend ---
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
                flag = '  [!] close to the integration floor';
            else
                flag = '';
            end
            % Both metrics are always shown; which one is the reference is
            % stated in the header and is the one used in the summaries.
            fprintf('  [%-6s] Phi %3d | GRE_full %9.4f%% | GRE_win %9.4f%% | Off %s | On %s%s\n', ...
                method, phi_val, gre_full, gre_win, off_str, on_str, flag);
            fprintf(log_file, '  %-6s Phi %03d | GRE_full %9.4f%% | GRE_win %9.4f%% | Off %s | On %s%s\n', ...
                method, phi_val, gre_full, gre_win, off_str, on_str, flag);
        end

        legend(axs(1), 'Location', 'best');
        sgtitle(fig, sprintf('%s method (Q = %d, k_{mult} = %g)', method, Q_val, K_val), ...
            'FontSize', 16, 'FontWeight', 'bold');

        base_name = fullfile(results_dir, sprintf('Compare_%s_Q%d_K%g', method, Q_val, K_val));
        exportgraphics(fig, [base_name '.png'], 'Resolution', 300);
        savefig(fig, [base_name '.fig']);
        fprintf('  -> figure saved: %s.png\n', base_name);
    end
    fprintf(log_file, '\n');
end

%% --- 5. Summary figures -----------------------------------------------
if isempty(summary)
    fprintf('\nNo ROM compared: summary figures skipped.\n');
    fclose(log_file);
    return
end

T_summary = struct2table(summary);
writetable(T_summary, fullfile(results_dir, sprintf('summary_%s.csv', lower(gre_mode))));

uniq_methods = unique(T_summary.method, 'stable');
mk = {'o-','s-','^-','d-','v-','>-'};

% --- Figure A: GRE convergence against the number of modes ---
figA = figure('Name','Convergence in phi','Color','w','Position',[100 100 800 600]);
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
xlabel('Number of retained modes \phi');
ylabel(gre_desc);
title('ROM convergence: reduction error vs basis size');
legend('Location','southwest'); box on;
exportgraphics(figA, fullfile(results_dir, ...
    sprintf('Summary_GRE_vs_phi_%s.png', lower(gre_mode))), 'Resolution', 300);
savefig(figA, fullfile(results_dir, sprintf('Summary_GRE_vs_phi_%s.fig', lower(gre_mode))));

% --- Figure B: accuracy / online cost trade-off ---
figB = figure('Name','Accuracy-cost Pareto','Color','w','Position',[100 100 800 600]);
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
xlabel('Online CPU time [s]');
ylabel(gre_desc);
title('Accuracy-cost trade-off (online)');
legend('Location','southwest'); box on;
exportgraphics(figB, fullfile(results_dir, ...
    sprintf('Summary_Pareto_%s.png', lower(gre_mode))), 'Resolution', 300);
savefig(figB, fullfile(results_dir, sprintf('Summary_Pareto_%s.fig', lower(gre_mode))));

%% --- 6. Observed convergence order ------------------------------------
fprintf('\n--- Observed convergence order (%s ~ phi^-p) ---\n', gre_tag);
fprintf(log_file, '\n--- Observed convergence order (%s ~ phi^-p) ---\n', gre_tag);
for m = 1:numel(uniq_methods)
    sel  = strcmp(T_summary.method, uniq_methods{m});
    phis = T_summary.phi(sel);
    g    = T_summary.gre_ref(sel);
    ok   = phis > 0 & g > 0;
    if nnz(ok) >= 2
        pfit = polyfit(log(phis(ok)), log(g(ok)), 1);
        fprintf('  %-8s : p = %.2f   (%d points)\n', uniq_methods{m}, -pfit(1), nnz(ok));
        fprintf(log_file, '  %-8s : p = %.2f   (%d points)\n', uniq_methods{m}, -pfit(1), nnz(ok));
    end
end

fclose(log_file);
fprintf('\nPost-processing complete. Results in %s\n', results_dir);
