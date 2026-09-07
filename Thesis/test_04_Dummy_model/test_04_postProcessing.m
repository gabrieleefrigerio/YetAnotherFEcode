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
%    - Summary_GRE_vs_ncc   error vs interface basis size, when interface
%                           reduction was used (Kuether et al. 2017, Table 2)
%    - Summary_Pareto       accuracy / online cost trade-off
%
%  No interpolation: every model is stored on the same t_common grid. If a
%  file does not respect it the script stops, rather than silently comparing
%  data that is not comparable.
% =====================================================================
% Set pp_results_dir before calling to drive this from another script without
% the folder dialog. It is consumed here and cleared with everything else, so
% running the script by hand ALWAYS asks: an override that survived would
% silently reuse the previous folder, which is worse than no override at all.
if exist('pp_results_dir', 'var') && ~isempty(pp_results_dir)
    pp_dir = pp_results_dir;
else
    pp_dir = '';
end
clearvars -except pp_dir; close all; clc;

addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'Src'));

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

% Tracking time: for how long does the ROM stay within a given error of the
% FOM. On a chaotic problem this discriminates better than a windowed GRE,
% because it does not saturate and needs no window chosen by hand. In the
% exponential regime one decade of threshold buys ln(10)/lambda of extra
% horizon, so the ranking is insensitive to the choice but the numbers are
% not: keep a set, not a single value.
%   headline = index of the threshold carried into the summary table.
track_thresholds = [0.1 1 10];    % [%]
track_headline   = 2;             % -> 1 %
track_hold_frac  = 0.01;          % a crossing counts after 1% of the run above

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
if isempty(pp_dir)
    results_dir = uigetdir(pwd, 'Select the results folder');
    if results_dir == 0
        error('PP:NoFolder', 'No folder selected.');
    end
else
    results_dir = pp_dir;
    if ~isfolder(results_dir)
        error('PP:NoFolder', 'pp_results_dir does not exist: %s', results_dir);
    end
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

% n_cc = 0 marks a run without interface reduction; ir_mode says which CC
% variant produced the reduced ones ('global' or 'per_interface').
summary = struct('method', {}, 'phi', {}, 'n_cc', {}, 'ir_mode', {}, ...
                 'Q', {}, 'K', {}, ...
                 'gre_full', {}, 'gre_win', {}, 'gre_ref', {}, ...
                 't_track_us', {}, 't_track_cens', {}, 'frac_track', {}, ...
                 'cpu', {}, 'offline', {});

% Tracking time at EVERY threshold, one row per ROM, kept apart from the
% summary table because it is a matrix: the table carries only the headline.
TRK.T    = [];      % [n_rom x k] horizon [s]
TRK.cens = [];      % [n_rom x k] logical, threshold never crossed
TRK.frac = [];      % [n_rom x k] fraction of the run spent inside the band
TRK.name = {};      % [n_rom x 1] label

%% --- 4. Loop over the cases (Q, k_mult) -------------------------------
for i_fom = 1:numel(fom_files)
    tok = regexp(fom_files(i_fom).name, 'FOM_Q([\d\.eE+-]+)_K([\d\.eE+-]+)\.mat', 'tokens');
    if isempty(tok), continue; end
    Q_val = str2double(tok{1}{1});
    K_val = str2double(tok{1}{2});

    fprintf('\n======================================================\n');
    fprintf('Case: Q = %g | k_mult = %g\n', Q_val, K_val);
    fprintf('======================================================\n');
    fprintf(log_file, '>>> CASE: Q = %g | k_mult = %g <<<\n', Q_val, K_val);

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
    rom_files = dir(fullfile(results_dir, sprintf('ROM_*_Q%g_K%g.mat', Q_val, K_val)));
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

        fig = figure('Name', sprintf('%s - Q%g - K%g', method, Q_val, K_val), ...
                     'NumberTitle', 'off', 'Color', 'w', ...
                     'Position', [100, 50, 1000, 250*(n_faces+1)]);
        n_sub = n_faces + 1;
        axs   = gobjects(n_sub, 1);

        % tiledlayout instead of subplot: it is what lets the legend get a
        % tile of its OWN ('north', below) rather than sitting inside an
        % axes and covering data. TileSpacing stays 'compact' as subplot
        % effectively was; the legend's row is added on top of that.
        tl = tiledlayout(fig, n_sub, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

        % --- one subplot per interface: FOM response and wall position ---
        for f = 1:n_faces
            axs(f) = nexttile(tl);
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
        axs(n_sub) = nexttile(tl);
        hold(axs(n_sub), 'on'); grid(axs(n_sub), 'on');
        title(axs(n_sub), 'Global Relative Error over time');
        % Tracking thresholds: the horizon in the summary is the abscissa at
        % which a curve leaves its band, so drawing the bands makes that
        % number readable off the figure instead of taken on trust.
        for t_thr = track_thresholds
            yline(axs(n_sub), t_thr, ':', sprintf('%g%%', t_thr), ...
                'Color', [0.45 0.45 0.45], 'LineWidth', 1.0, ...
                'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
        end
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

            % Interface reduction: absent from the file or from the name means
            % this run had none.
            if isfield(rom, 'n_cc')
                n_cc_val = rom.n_cc;
            else
                ct = regexp(sel_files(i_rom).name, 'CC(\d+)', 'tokens', 'once');
                if isempty(ct), n_cc_val = 0; else, n_cc_val = str2double(ct{1}); end
            end
            if isfield(rom, 'ir_mode'), ir_mode_val = rom.ir_mode; else, ir_mode_val = 'none'; end

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

            % Tracking time, fed the SAME curve that gets plotted so that the
            % number in the table and the crossing visible in the figure can
            % never disagree.
            rom_label = sprintf('%s Phi%d', method, phi_val);
            if n_cc_val > 0
                rom_label = sprintf('%s CC%d', rom_label, n_cc_val);
            end
            trk = rom_tracking_time(t_ref, gre_t, track_thresholds, ...
                'HoldFraction', track_hold_frac, 'Label', rom_label);
            TRK.T(end+1, :)    = trk.T;
            TRK.cens(end+1, :) = trk.censored;
            TRK.frac(end+1, :) = trk.frac;
            TRK.name{end+1, 1} = rom_label;

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
            if n_cc_val > 0
                cc_tag = sprintf(', n_{cc}=%d', n_cc_val);
            else
                cc_tag = '';
            end
            if isempty(time_info)
                legend_str = sprintf('ROM \\phi=%d%s (%s: %.3f%%)', ...
                    phi_val, cc_tag, gre_tag_tex, gre_ref);
            else
                legend_str = sprintf('ROM \\phi=%d%s (%s: %.3f%%, %s)', ...
                    phi_val, cc_tag, gre_tag_tex, gre_ref, time_info);
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
                'n_cc', n_cc_val, 'ir_mode', ir_mode_val, ...
                'Q', Q_val, 'K', K_val, 'gre_full', gre_full, 'gre_win', gre_win, ...
                'gre_ref', gre_ref, ...
                't_track_us', 1e6*trk.T(track_headline), ...
                't_track_cens', trk.censored(track_headline), ...
                'frac_track', trk.frac(track_headline), ...
                'cpu', rom_cpu, 'offline', rom_off); %#ok<SAGROW>

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
            if n_cc_val > 0, cc_str = sprintf('%3d', n_cc_val); else, cc_str = '  -'; end
            fprintf('  [%-9s] Phi %3d | CC %s | GRE_full %9.4f%% | GRE_win %9.4f%% | Off %s | On %s%s\n', ...
                method, phi_val, cc_str, gre_full, gre_win, off_str, on_str, flag);
            fprintf(log_file, '  %-9s Phi %03d | CC %s | GRE_full %9.4f%% | GRE_win %9.4f%% | Off %s | On %s%s\n', ...
                method, phi_val, cc_str, gre_full, gre_win, off_str, on_str, flag);
            if trk.censored(track_headline), pre = '>'; else, pre = ' '; end
            fprintf('                 within %g%%: %s%.3f us  (inside the band %.0f%% of the run)\n', ...
                track_thresholds(track_headline), pre, 1e6*trk.T(track_headline), ...
                100*trk.frac(track_headline));
            fprintf(log_file, '                 within %g%%: %s%.3f us  (inside %.0f%%)\n', ...
                track_thresholds(track_headline), pre, 1e6*trk.T(track_headline), ...
                100*trk.frac(track_headline));
        end

        % A tile of its own ('north') rather than 'Location','best' inside
        % axs(1): with many ROM curves 'best' has nowhere left to hide the
        % legend without covering a trace. The north tile reserves exactly
        % the height the legend needs and pushes every plot down by that
        % much, so nothing is ever covered.
        lg = legend(axs(1));
        lg.Layout.Tile = 'north';
        lg.Orientation = 'horizontal';
        lg.NumColumns  = min(numel(sel_files) + 2, 4);

        sgtitle(fig, sprintf('%s method (Q = %g, k_{mult} = %g)', method, Q_val, K_val), ...
            'FontSize', 16, 'FontWeight', 'bold');

        base_name = fullfile(results_dir, sprintf('Compare_%s_Q%g_K%g', method, Q_val, K_val));
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

% --- Tracking time at every threshold ----------------------------------
% Its own file because it is a matrix, not a column of the summary. A ">" in
% the console means the threshold was never crossed: that row is CENSORED and
% its number is a lower bound, so it must not be averaged in with the others
% nor read as a horizon.
if ~isempty(TRK.name)
    T_track = table(TRK.name, 'VariableNames', {'rom'});
    for i_thr = 1:numel(track_thresholds)
        tag = strrep(strrep(sprintf('%g', track_thresholds(i_thr)), '.', 'p'), '-', 'm');
        T_track.(sprintf('T_%s_us',   tag)) = 1e6 * TRK.T(:, i_thr);
        T_track.(sprintf('cens_%s',   tag)) = TRK.cens(:, i_thr);
        T_track.(sprintf('inside_%s', tag)) = TRK.frac(:, i_thr);
    end
    writetable(T_track, fullfile(results_dir, 'tracking_time.csv'));

    fprintf('\n--- Tracking time: how long the ROM stays within a threshold ---\n');
    fprintf(log_file, '\n--- Tracking time ---\n');
    hdr = sprintf('%-24s', 'ROM');
    for thr = track_thresholds
        hdr = [hdr sprintf('%12s', sprintf('< %g%%', thr))]; %#ok<AGROW>
    end
    fprintf('%s\n', hdr); fprintf(log_file, '%s\n', hdr);
    for r = 1:numel(TRK.name)
        row = sprintf('%-24s', TRK.name{r});
        for i_thr = 1:numel(track_thresholds)
            if TRK.cens(r, i_thr), pre = '>'; else, pre = ' '; end
            row = [row sprintf('%11s', sprintf('%s%.3fus', pre, 1e6*TRK.T(r, i_thr)))]; %#ok<AGROW>
        end
        fprintf('%s\n', row); fprintf(log_file, '%s\n', row);
    end

    figT = figure('Name','Tracking time','Color','w','Position',[100 100 1000 560]);
    hold on; grid on; box on
    bh = bar(1e6 * TRK.T, 'grouped');
    for i_thr = 1:numel(track_thresholds)
        cs = logical(TRK.cens(:, i_thr));   % assignment into [] makes it double
        if any(cs)
            plot(bh(i_thr).XEndPoints(cs), 1e6*TRK.T(cs, i_thr), '^k', ...
                'MarkerFaceColor', 'w', 'MarkerSize', 7, 'HandleVisibility', 'off');
        end
        bh(i_thr).DisplayName = sprintf('within %g%%', track_thresholds(i_thr));
    end
    set(gca, 'XTick', 1:numel(TRK.name), 'XTickLabel', TRK.name, ...
        'TickLabelInterpreter', 'none', 'XTickLabelRotation', 45);
    ylabel('time inside the error band [\mus]');
    title(['How long each ROM tracks the FOM    ' ...
        '(open marker = never left the band, so the bar is a lower bound)']);
    legend('Location', 'northeast');
    exportgraphics(figT, fullfile(results_dir, 'Summary_TrackingTime.png'), 'Resolution', 300);
    savefig(figT, fullfile(results_dir, 'Summary_TrackingTime.fig'));
end

uniq_methods = unique(T_summary.method, 'stable');
mk = {'o-','s-','^-','d-','v-','>-'};

% --- Figure A: GRE convergence against the number of modes ---
% Restricted to runs without interface reduction, so that phi is the only
% parameter varying along each curve. The reduced ones get their own figure.
base_rows = T_summary.n_cc == 0;
figA = figure('Name','Convergence in phi','Color','w','Position',[100 100 800 600]);
hold on; grid on;
for m = 1:numel(uniq_methods)
    sel = strcmp(T_summary.method, uniq_methods{m}) & base_rows;
    if ~any(sel), continue; end
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
figB = plot_pareto(T_summary, gre_desc, gre_floor_pct);
exportgraphics(figB, fullfile(results_dir, ...
    sprintf('Summary_Pareto_%s.png', lower(gre_mode))), 'Resolution', 300);
savefig(figB, fullfile(results_dir, sprintf('Summary_Pareto_%s.fig', lower(gre_mode))));

% --- Figure C: GRE against the number of interface (CC) modes ---
% The equivalent of Table 2 in Kuether et al. 2017: how far the interface can
% be truncated before accuracy degrades. One curve per (method, phi) pair, so
% that the two CC variants (global vs per_interface) can be read side by side.
ir_rows = T_summary.n_cc > 0;
if any(ir_rows)
    figC = figure('Name','Convergence in n_cc','Color','w','Position',[100 100 800 600]);
    hold on; grid on;

    ir_methods = unique(T_summary.method(ir_rows), 'stable');
    ir_phis    = unique(T_summary.phi(ir_rows));
    ic = 0;
    for m = 1:numel(ir_methods)
        for p = 1:numel(ir_phis)
            sel = ir_rows & strcmp(T_summary.method, ir_methods{m}) & ...
                  T_summary.phi == ir_phis(p);
            if nnz(sel) < 2, continue; end
            ic = ic + 1;
            [ccs, iord] = sort(T_summary.n_cc(sel));
            g = T_summary.gre_ref(sel);
            plot(ccs, g(iord), mk{mod(ic-1,numel(mk))+1}, 'LineWidth', 1.8, ...
                'MarkerSize', 7, 'DisplayName', ...
                sprintf('%s, \\phi=%d', ir_methods{m}, ir_phis(p)));
        end
    end

    % Reference: the same methods without interface reduction. Any CC curve
    % that reaches this level has lost nothing to the interface truncation.
    for m = 1:numel(ir_methods)
        base_name = regexprep(ir_methods{m}, 'IR[GP]$', '');
        sel = base_rows & strcmp(T_summary.method, base_name);
        if any(sel)
            yline(min(T_summary.gre_ref(sel)), 'k--', 'LineWidth', 1.2, ...
                'DisplayName', sprintf('%s, no IR (best)', base_name));
        end
    end

    yline(gre_floor_pct, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Integration floor');
    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlabel('Number of retained interface (CC) modes n_{cc}');
    ylabel(gre_desc);
    title('Interface reduction: error vs interface basis size');
    legend('Location','best'); box on;
    exportgraphics(figC, fullfile(results_dir, ...
        sprintf('Summary_GRE_vs_ncc_%s.png', lower(gre_mode))), 'Resolution', 300);
    savefig(figC, fullfile(results_dir, sprintf('Summary_GRE_vs_ncc_%s.fig', lower(gre_mode))));
end

%% --- 6. Observed convergence order ------------------------------------
fprintf('\n--- Observed convergence order (%s ~ phi^-p) ---\n', gre_tag);
fprintf(log_file, '\n--- Observed convergence order (%s ~ phi^-p) ---\n', gre_tag);
for m = 1:numel(uniq_methods)
    % Only the runs without interface reduction: with IR, phi is fixed along a
    % curve and n_cc is what varies, so a fit in phi would be meaningless.
    sel  = strcmp(T_summary.method, uniq_methods{m}) & base_rows;
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
