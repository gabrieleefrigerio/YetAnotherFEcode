%% =====================================================================
%  POST-PROCESSING 3D - FOM vs ROM comparison and convergence
%
%  Same contract as the 2D one: faces, gaps and normals are NOT written here
%  but read from run_config.mat in the results folder, so this script cannot
%  drift out of sync with the model that was actually run.
%
%  What differs from the 2D version, and why:
%
%  * The 2D script picked the component to plot from Interfaces.<face>.dir
%    (1 = X, 2 = Y). In 3D there is no such thing: a wall has a NORMAL, and
%    extract_contact_response already stores the projection along it in
%    y_contact.<face>.normal. That is the component the gap is measured on,
%    whatever the orientation of the wall, so everything here works on it.
%
%  * One subplot per bumpstop, but each shows the whole face rather than one
%    arbitrary node:
%       - a shaded band spanning min..max over the 29 nodes of the FOM. Its
%         WIDTH is the flat/corner verdict at a glance: a face landing
%         parallel collapses to a thin ribbon, one touching on an edge fans
%         out.
%       - a solid line for one reference node, chosen by physics and not by
%         index: the deepest-penetrating node of the FOM, i.e. the one that
%         governs the contact force. If the stop never engages, the one with
%         the largest normal excursion.
%       - the ROMs, dashed, on that SAME node index, so the comparison is
%         like for like.
%    One line per ROM keeps the legend readable, which 29 lines would not.
%
%  * Layout: the eight stoppers in a 4x2 grid, one tab per row with its top
%    and bottom face side by side, and GRE(t) spanning the width underneath.
%    Pairing top and bottom is the comparison that matters: did the mass go
%    up, down, or rattle between the two.
%
%  Outputs: one figure per method, GRE_Results.txt, summary_<mode>.csv, and
%  the three summary figures (convergence in phi, Pareto, convergence in n_cc).
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
% Same single-use contract for the animation: a batch that only wants the
% numbers should not spend minutes rendering a GIF it will never open.
if exist('pp_make_animation', 'var') && ~isempty(pp_make_animation)
    pp_anim = logical(pp_make_animation);
else
    pp_anim = [];
end
clearvars -except pp_dir pp_anim; close all; clc;

addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'Src'));

%% --- Options ----------------------------------------------------------

% Reference GRE metric shown in the legends and the summaries:
%   'full'   -> whole time history
%   'window' -> first win_frac of the simulation
% Past the first impacts the positive Lyapunov exponent makes trajectories
% diverge and the error saturates for physical reasons. Both are always
% computed and logged; see rom_divergence for where the boundary sits.
gre_mode = 'full';
win_frac = 0.25;

% Integration floor, MEASURED on this model by comparing the two 100 us FOM
% runs at RelTolFOM 1e-7 and 1e-8: they differ by 4.9e-7 % at 1 us, 6.1e-6 %
% at 10 us and 2.5e-2 % at 100 us. Below that curve a GRE measures the
% integrator, not the reduction. The single number kept here is the value at
% 10 us: conservative on short windows, optimistic on long ones, which is one
% more reason to read it together with the tracking time below.
gre_floor_pct = 6.1e-6;

gre_plot_log = false;      % log scale on the GRE(t) subplot

% Zoom on one tile of the transient figures drives every other tile along the
% time axis. Set this to true to share the vertical scale between the face
% tiles as well, which makes the bumpstops directly comparable but hides the
% detail of the quieter ones.
link_y_faces = false;

% Tracking time: for how long does the ROM stay within a given error of the
% FOM. On a chaotic problem this discriminates better than a windowed GRE,
% because it does not saturate and needs no window chosen by hand. In the
% exponential regime one decade of threshold buys ln(10)/lambda of extra
% horizon, measured at 24 us on this model, so the ranking is insensitive to
% the choice but the numbers are not: keep a set, not a single value.
%   headline = index of the threshold carried into the summary table.
track_thresholds = [0.1 1 10];    % [%]
track_headline   = 2;             % -> 1 %
track_hold_frac  = 0.01;          % a crossing counts after 1% of the run above

% What the per-stopper subplots show:
%   'physical'  -> displacement along a single reference direction, the same
%                  for every panel (the shock direction), with each wall drawn
%                  at its own SIGNED position. A proof mass driven along +Z
%                  rises in every panel; the upper walls sit at +gap and get
%                  hit, the lower ones at -gap and are left behind.
%   'clearance' -> distance still to travel before touching, gap - u_n
%   'normal'    -> the raw normal component u_n, with the wall at +gap
%
% 'physical' is the default because the other two both fight intuition, in
% opposite directions. With 'normal', a stopper underneath has normal
% [0 0 -1], so a mass moving UP plots as going DOWN there. With 'clearance',
% every curve falls towards the wall, so the upper stoppers plot as going
% down while the mass rises. Only the common reference direction lets the
% eight panels be read as one picture of the same motion.
plot_quantity = 'physical';

% Spatial view of the contact nodes on the mesh (needs the mesh to be
% rebuilt, a few seconds). Off by default for a quick pass.
show_spatial_view = false;

% Animation of the stoppers moving against their walls, written as a GIF next
% to the other figures. Costs a minute or two of rendering, so set it to false
% for a quick pass over the numbers.
%
% The whole deformed mesh is animated only if the FOM was run with
% cfg.save_full_field, which makes run_fom store the full displacement field.
% Without it the animation still works, but the mesh stays static and only the
% contact nodes move: the result file holds the response at those nodes alone.
make_animation = false;
if ~isempty(pp_anim), make_animation = pp_anim; end
anim_amplify   = 8;      % displacements are microns on a structure 400 um wide
anim_stride    = 1;      % one frame every N output samples

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
    results_dir = uigetdir(fullfile(fileparts(mfilename('fullpath')), 'results'), ...
        'Select the results folder');
else
    results_dir = pp_dir;
end
if isequal(results_dir, 0)
    error('PP:NoFolder', 'No folder selected.');
end
fprintf('Folder: %s\n', results_dir);

cfg_file = fullfile(results_dir, 'run_config.mat');
if ~exist(cfg_file, 'file')
    error('PP:NoConfig', 'run_config.mat not found in\n  %s', results_dir);
end
R = load(cfg_file);

%% --- 2. Interface table (read from the file, not rewritten by hand) ----
[faces, face_gap, face_nrm, d_ref] = interface_table(R);
n_faces = numel(faces);
fprintf('Eref = %.4e\n', R.Eref);

%% --- 3. Log and accumulator -------------------------------------------
fom_files = dir(fullfile(results_dir, 'FOM_*.mat'));
if isempty(fom_files)
    error('PP:NoFOM', 'No FOM_*.mat file: a reference is needed for the comparison.');
end

% Closed explicitly at the end: in a script the variables stay in the base
% workspace, so an onCleanup object would never be destroyed.
log_file = fopen(fullfile(results_dir, 'GRE_Results.txt'), 'w');
fprintf(log_file, '======================================================\n');
fprintf(log_file, ' GLOBAL RELATIVE ERROR (GRE) & TIME REPORT - 3D\n');
fprintf(log_file, '======================================================\n');
fprintf(log_file, ' Interfaces: %s\n', strjoin(faces, ', '));
fprintf(log_file, ' Reference metric: %s (window %.0f%%)\n', gre_tag, 100*win_frac);
fprintf(log_file, ' Integration floor: %.2e %%\n', gre_floor_pct);
fprintf(log_file, ' Eref = %.4e\n\n', R.Eref);

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

    % Named variables only. With cfg.save_full_field the FOM file also holds
    % the full displacement field, which is hundreds of MB and is needed by
    % the animation alone: loading it here would cost that on every pass.
    fom_file = fullfile(results_dir, fom_files(i_fom).name);
    fom   = load(fom_file, 't', 'y_contact', 'cpu_time');
    t_ref = fom.t(:);
    t_us  = 1e6 * t_ref;
    nT    = numel(t_ref);
    i_win = 1 : max(2, round(win_frac * nT));
    if strcmpi(gre_mode, 'window'), i_gre = i_win; else, i_gre = 1:nT; end

    if isfield(fom, 'cpu_time')
        fom_legend = sprintf('FOM (On: %.0fs)', fom.cpu_time);
    else
        fom_legend = 'FOM';
    end

    % --- reference node per face, chosen from the FOM ------------------
    % The deepest-penetrating node governs the contact force, so it is the
    % one worth following. With no contact at all, the largest excursion.
    ref_node = ones(1, n_faces);
    for f = 1:n_faces
        un  = fom.y_contact.(faces{f}).normal;
        pen = max(un - face_gap(f), [], 2);
        if any(pen > 0)
            [~, ref_node(f)] = max(pen);
        else
            [~, ref_node(f)] = max(max(abs(un), [], 2));
        end
    end

    % --- ROMs available for this case, grouped by method ---------------
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

        n_rows = ceil(n_faces/2) + 1;          % faces in a 2-column grid, GRE below
        fig = figure('Name', sprintf('%s - Q%g - K%g', method, Q_val, K_val), ...
                     'NumberTitle', 'off', 'Color', 'w', ...
                     'Position', [60, 40, 1200, 210*n_rows]);
        tl  = tiledlayout(fig, n_rows, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
        axs = gobjects(n_faces + 1, 1);

        % Handles kept so that every curve of one model, across ALL the tiles,
        % can be switched as a single object. Without this the plot browser
        % only ever hides the copy that carries the legend entry, which is the
        % one in the first tile, and the other faces keep drawing it.
        h_fom  = gobjects(n_faces, 1);
        h_wall = gobjects(n_faces, 1);
        vis_links = {};        % linkprop handles: they must stay alive, so
                               % they end up in the figure appdata below

        % --- one tile per bumpstop: FOM band, reference node, wall ------
        for f = 1:n_faces
            axs(f) = nexttile(tl);
            hold(axs(f), 'on'); grid(axs(f), 'on'); box(axs(f), 'on');

            [un, wall_level, ylab] = contact_plot_quantity(fom.y_contact.(faces{f}), ...
                face_nrm{f}, face_gap(f), plot_quantity, d_ref);
            u_lo  = min(un, [], 1);
            u_hi  = max(un, [], 1);
            u_ref = un(ref_node(f), :);

            % Spread over the whole face. A narrow ribbon means the face
            % lands parallel; a wide one means it touches on an edge.
            fill(axs(f), [t_us; flipud(t_us)], [u_hi(:); flipud(u_lo(:))], ...
                [0.25 0.25 0.25], 'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
                'HandleVisibility', 'off');

            h_fom(f)  = plot(axs(f), t_us, u_ref(:), 'k-', 'LineWidth', 2);
            h_wall(f) = yline(axs(f), wall_level, 'r-.', 'LineWidth', 1.4);
            if f == 1
                % Only the first tile feeds the legend, otherwise every entry
                % would appear once per face.
                h_fom(f).DisplayName  = fom_legend;
                h_wall(f).DisplayName = 'Wall';
            else
                h_fom(f).HandleVisibility  = 'off';
                h_wall(f).HandleVisibility = 'off';
            end

            spread = max(u_hi - u_lo);   % um: the mesh is in the um/MPa/kg system
            % Interpreter stays 'none': the face labels (e.g. stop_A_top) have
            % underscores that TeX would read as subscripts and mangle.
            title(axs(f), sprintf('%s  -  node %d  -  face spread %.3g um', ...
                faces{f}, ref_node(f), spread), 'Interpreter', 'none');
            if f > n_faces - 2, xlabel(axs(f), 't [\mus]'); end
            ylabel(axs(f), ylab);
        end

        % --- GRE(t), spanning the width -------------------------------
        axs(end) = nexttile(tl, [1 2]);
        hold(axs(end), 'on'); grid(axs(end), 'on'); box(axs(end), 'on');
        title(axs(end), 'Global Relative Error over time');
        % Tracking thresholds: the horizon in the summary is the abscissa at
        % which a curve leaves its band, so drawing the bands makes that
        % number readable off the figure instead of taken on trust.
        for t_thr = track_thresholds
            yline(axs(end), t_thr, ':', sprintf('%g%%', t_thr), ...
                'Color', [0.45 0.45 0.45], 'LineWidth', 1.0, ...
                'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
        end
        xlabel(axs(end), 't [\mus]'); ylabel(axs(end), 'GRE [%]');
        if gre_plot_log
            set(axs(end), 'YScale', 'log');
            yline(axs(end), gre_floor_pct, 'r--', 'LineWidth', 1.2, ...
                'HandleVisibility', 'off');
        end
        if strcmpi(gre_mode, 'window')
            xline(axs(end), t_us(i_win(end)), 'k:', 'LineWidth', 1.2, ...
                'HandleVisibility', 'off');
        end

        colors = lines(numel(sel_files));

        for i_rom = 1:numel(sel_files)
            rom   = load(fullfile(results_dir, sel_files(i_rom).name));
            t_rom = rom.t(:);

            % Common grid, no interpolation.
            if numel(t_rom) ~= nT || max(abs(t_rom - t_ref)) > 1e-9*(t_ref(end)-t_ref(1))
                error('PP:GridMismatch', ...
                    ['Inconsistent time grid in %s.\n' ...
                     'FOM: %d points, ROM: %d points.'], ...
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

            if isfield(rom, 'n_cc')
                n_cc_val = rom.n_cc;
            else
                ct = regexp(sel_files(i_rom).name, 'CC(\d+)', 'tokens', 'once');
                if isempty(ct), n_cc_val = 0; else, n_cc_val = str2double(ct{1}); end
            end
            if isfield(rom, 'ir_mode'), ir_mode_val = rom.ir_mode; else, ir_mode_val = 'none'; end

            % --- all contact DOFs concatenated, normal component -------
            y_fom_cat = [];
            y_rom_cat = [];
            for f = 1:n_faces
                y_fom_cat = [y_fom_cat; fom.y_contact.(faces{f}).normal]; %#ok<AGROW>
                y_rom_cat = [y_rom_cat; rom.y_contact.(faces{f}).normal]; %#ok<AGROW>
            end

            gre_full = norm(y_fom_cat - y_rom_cat, 'fro') / norm(y_fom_cat, 'fro') * 100;
            gre_win  = norm(y_fom_cat(:,i_win) - y_rom_cat(:,i_win), 'fro') / ...
                       norm(y_fom_cat(:,i_win), 'fro') * 100;
            gre_ref  = norm(y_fom_cat(:,i_gre) - y_rom_cat(:,i_gre), 'fro') / ...
                       norm(y_fom_cat(:,i_gre), 'fro') * 100;

            % Instantaneous GRE, normalized on the maximum of the whole
            % history so the curve does not depend on gre_mode.
            norm_diff_t = sqrt(sum((y_fom_cat - y_rom_cat).^2, 1));
            gre_t = norm_diff_t ./ (max(sqrt(sum(y_fom_cat.^2, 1))) + eps) * 100;

            % Tracking time. It is fed the SAME curve that gets plotted, so
            % the number in the table and the crossing visible in the figure
            % can never disagree.
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
            if ~isnan(rom_off), time_info = sprintf('Off: %.1fs', rom_off); end
            if ~isnan(rom_cpu)
                if isempty(time_info)
                    time_info = sprintf('On: %.1fs', rom_cpu);
                else
                    time_info = sprintf('%s, On: %.1fs', time_info, rom_cpu);
                end
            end
            if n_cc_val > 0, cc_tag = sprintf(', n_{cc}=%d', n_cc_val); else, cc_tag = ''; end
            if isempty(time_info)
                legend_str = sprintf('ROM \\phi=%d%s (%s: %.3f%%)', ...
                    phi_val, cc_tag, gre_tag_tex, gre_ref);
            else
                legend_str = sprintf('ROM \\phi=%d%s (%s: %.3f%%, %s)', ...
                    phi_val, cc_tag, gre_tag_tex, gre_ref, time_info);
            end

            % One handle per tile plus the GRE curve, linked below so that
            % hiding this model hides its error curve too.
            h_this = gobjects(n_faces + 1, 1);

            % Same reference node as the FOM, so the curves are comparable.
            for f = 1:n_faces
                y_rom_all  = contact_plot_quantity(rom.y_contact.(faces{f}), face_nrm{f}, ...
                                     face_gap(f), plot_quantity, d_ref);
                y_rom_face = y_rom_all(ref_node(f), :);
                h_this(f) = plot(axs(f), t_us, y_rom_face(:), '--', ...
                    'Color', colors(i_rom,:), 'LineWidth', 1.4);
                if f == 1
                    h_this(f).DisplayName = legend_str;
                else
                    h_this(f).HandleVisibility = 'off';
                end
            end
            h_this(end) = plot(axs(end), t_us, gre_t(:), '-', ...
                'Color', colors(i_rom,:), 'LineWidth', 1.4, ...
                'HandleVisibility', 'off');
            vis_links{end+1} = linkprop(h_this, 'Visible'); %#ok<SAGROW>

            summary(end+1) = struct('method', method, 'phi', phi_val, ...
                'n_cc', n_cc_val, 'ir_mode', ir_mode_val, ...
                'Q', Q_val, 'K', K_val, 'gre_full', gre_full, 'gre_win', gre_win, ...
                'gre_ref', gre_ref, ...
                't_track_us', 1e6*trk.T(track_headline), ...
                't_track_cens', trk.censored(track_headline), ...
                'frac_track', trk.frac(track_headline), ...
                'cpu', rom_cpu, 'offline', rom_off); %#ok<SAGROW>

            % --- log ---
            if isnan(rom_off), off_str = 'N/A'; else, off_str = sprintf('%7.1fs', rom_off); end
            if isnan(rom_cpu), on_str  = 'N/A'; else, on_str  = sprintf('%7.1fs', rom_cpu); end
            if gre_ref < 100 * gre_floor_pct
                flag = '  [!] close to the integration floor';
            else
                flag = '';
            end
            if n_cc_val > 0, cc_str = sprintf('%3d', n_cc_val); else, cc_str = '  -'; end
            fprintf('  [%-9s] Phi %3d | CC %s | GRE_full %9.4f%% | GRE_win %9.4f%% | Off %s | On %s%s\n', ...
                method, phi_val, cc_str, gre_full, gre_win, off_str, on_str, flag);
            fprintf(log_file, '  %-9s Phi %03d | CC %s | GRE_full %9.4f%% | GRE_win %9.4f%% | Off %s | On %s%s\n', ...
                method, phi_val, cc_str, gre_full, gre_win, off_str, on_str, flag);
            if trk.censored(track_headline), pre = '>'; else, pre = ' '; end
            fprintf('               within %g%%: %s%.2f us  (inside the band %.0f%% of the run)\n', ...
                track_thresholds(track_headline), pre, 1e6*trk.T(track_headline), ...
                100*trk.frac(track_headline));
            fprintf(log_file, '               within %g%%: %s%.2f us  (inside %.0f%%)\n', ...
                track_thresholds(track_headline), pre, 1e6*trk.T(track_headline), ...
                100*trk.frac(track_headline));
        end

        % The FOM and the walls are single objects across the tiles too.
        vis_links{end+1} = linkprop(h_fom,  'Visible'); %#ok<SAGROW>
        vis_links{end+1} = linkprop(h_wall, 'Visible'); %#ok<SAGROW>
        % linkprop dies with its variable, so the links have to be parked on
        % the figure or the browser stops propagating as soon as this scope ends.
        setappdata(fig, 'VisibilityLinks', vis_links);

        % Zooming or panning one tile moves them all along the time axis. Y is
        % deliberately NOT linked: each bumpstop keeps its own scale, and the
        % GRE tile is a per cent, not a displacement. Set link_y_faces to true
        % to put the face tiles on a common vertical scale as well.
        linkaxes(axs, 'x');
        if link_y_faces
            setappdata(fig, 'YLink', linkprop(axs(1:n_faces), 'YLim'));
        end

        legend(axs(1), 'Location', 'best');
        switch lower(plot_quantity)
            case 'physical'
                qdesc = sprintf('displacement along [%g %g %g], walls at their signed gap', d_ref);
            case 'clearance'
                qdesc = 'clearance to the wall (0 = touching)';
            otherwise
                qdesc = 'displacement along each wall normal';
        end
        title(tl, sprintf(['%s method (Q = %g, k_{mult} = %g)   -   %s' ...
            '   -   shaded band = spread over the face'], ...
            method, Q_val, K_val, qdesc), 'FontSize', 13, 'FontWeight', 'bold');

        base_name = fullfile(results_dir, sprintf('Compare_%s_Q%g_K%g', method, Q_val, K_val));
        exportgraphics(fig, [base_name '.png'], 'Resolution', 200);
        savefig(fig, [base_name '.fig']);
        fprintf('  -> figure saved: %s.png\n', base_name);
    end

    % --- spatial view and animation of the contact, once per case ------
    % Both draw on the mesh, so the structure is rebuilt once for the two.
    if show_spatial_view || make_animation
        try
            Struct = rebuild_structure(R.cfg, results_dir);

            if show_spatial_view
                plot_contact_nodes(Struct, R.Interfaces, faces, fom.y_contact);
                title(sprintf('FOM contact nodes - Q %g, k_{mult} %g', Q_val, K_val));
            end

            if make_animation
                % The full field is optional: whos tells whether the run saved
                % it, so a run without it animates the contact nodes on a
                % static mesh instead of failing.
                u_full = [];
                if ismember('u_full', {whos('-file', fom_file).name})
                    fprintf('  Loading the full displacement field...\n');
                    S = load(fom_file, 'u_full');
                    u_full = S.u_full;
                    clear S;
                else
                    fprintf(['  [note] The run has no full field ' ...
                        '(cfg.save_full_field was off): the mesh stays static.\n']);
                end

                gif_name = fullfile(results_dir, ...
                    sprintf('Animation_FOM_Q%g_K%g.gif', Q_val, K_val));
                fprintf('  Rendering the animation (stride %d)...\n', anim_stride);
                animate_contact_3d(Struct, R.Interfaces, faces, fom.y_contact, ...
                    t_ref, 'Amplify', anim_amplify, 'Stride', anim_stride, ...
                    'Displacement', u_full, 'GifFile', gif_name);
                clear u_full;
                fprintf('  -> animation saved: %s\n', gif_name);
            end
        catch ME
            warning('PP:SpatialView', 'Spatial view / animation skipped: %s', ME.message);
        end
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
% Kept in its own file because it is a matrix, not a column of the summary.
% A ">" in the console means the threshold was never crossed: that row is
% CENSORED and its number is a lower bound, so it must not be averaged in
% with the others or read as a horizon.
if ~isempty(TRK.name)
    T_track = table(TRK.name, 'VariableNames', {'rom'});
    for i_thr = 1:numel(track_thresholds)
        tag = sprintf('%g', track_thresholds(i_thr));
        tag = strrep(strrep(tag, '.', 'p'), '-', 'm');
        T_track.(sprintf('T_%s_us',   tag)) = 1e6 * TRK.T(:, i_thr);
        T_track.(sprintf('cens_%s',   tag)) = TRK.cens(:, i_thr);
        T_track.(sprintf('inside_%s', tag)) = TRK.frac(:, i_thr);
    end
    writetable(T_track, fullfile(results_dir, 'tracking_time.csv'));

    fprintf('\n--- Tracking time: how long the ROM stays within a threshold ---\n');
    fprintf(log_file, '\n--- Tracking time ---\n');
    hdr = sprintf('%-22s', 'ROM');
    for thr = track_thresholds, hdr = [hdr sprintf('%12s', sprintf('< %g%%', thr))]; end %#ok<AGROW>
    fprintf('%s\n', hdr); fprintf(log_file, '%s\n', hdr);
    for r = 1:numel(TRK.name)
        line = sprintf('%-22s', TRK.name{r});
        for i_thr = 1:numel(track_thresholds)
            if TRK.cens(r, i_thr), pre = '>'; else, pre = ' '; end
            line = [line sprintf('%11s', sprintf('%s%.2fus', pre, 1e6*TRK.T(r, i_thr)))]; %#ok<AGROW>
        end
        fprintf('%s\n', line); fprintf(log_file, '%s\n', line);
    end

    figT = figure('Name','Tracking time','Color','w','Position',[100 100 900 520]);
    hold on; grid on; box on
    bh = bar(1e6 * TRK.T, 'grouped');
    % Censored bars are hatched by an open marker on top: the bar is a lower
    % bound, and without the mark a saturated run reads as the best result.
    for i_thr = 1:numel(track_thresholds)
        cs = logical(TRK.cens(:, i_thr));   % assignment into [] makes it double
        if any(cs)
            plot(bh(i_thr).XEndPoints(cs), 1e6*TRK.T(cs, i_thr), '^k', ...
                'MarkerFaceColor', 'w', 'MarkerSize', 7, 'HandleVisibility', 'off');
        end
        bh(i_thr).DisplayName = sprintf('within %g%%', track_thresholds(i_thr));
    end
    set(gca, 'XTick', 1:numel(TRK.name), 'XTickLabel', TRK.name, ...
        'TickLabelInterpreter', 'none', 'XTickLabelRotation', 30);
    ylabel('time inside the error band [\mus]');
    title(['How long each ROM tracks the FOM    ' ...
        '(open marker = never left the band, so the bar is a lower bound)']);
    legend('Location', 'northwest');
    exportgraphics(figT, fullfile(results_dir, 'Summary_TrackingTime.png'), 'Resolution', 300);
    savefig(figT, fullfile(results_dir, 'Summary_TrackingTime.fig'));
end

uniq_methods = unique(T_summary.method, 'stable');
mk = {'o-','s-','^-','d-','v-','>-'};

% --- Figure A: GRE convergence against the number of modes ---
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
xlabel('Number of retained modes \phi'); ylabel(gre_desc);
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
    for m = 1:numel(ir_methods)
        bn = regexprep(ir_methods{m}, 'IR[GP]$', '');
        sel = base_rows & strcmp(T_summary.method, bn);
        if any(sel)
            yline(min(T_summary.gre_ref(sel)), 'k--', 'LineWidth', 1.2, ...
                'DisplayName', sprintf('%s, no IR (best)', bn));
        end
    end
    yline(gre_floor_pct, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Integration floor');
    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlabel('Number of retained interface (CC) modes n_{cc}'); ylabel(gre_desc);
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

%% =====================================================================
% =====================================================================
function Struct = rebuild_structure(cfg, results_dir) %#ok<INUSD>
% REBUILD_STRUCTURE The model behind a results folder, for the spatial view.
% Only needed by plot_contact_nodes, which draws on the mesh.
    here = fileparts(mfilename('fullpath'));
    Struct = FeStructure();
    Struct.mesh_file    = fullfile(here, 'mesh', cfg.mesh_file);
    Struct.element_type = cfg.element_type;
    Struct.E = cfg.E; Struct.nu = cfg.nu; Struct.rho = cfg.rho;
    Struct.set_specs = cfg.node_sets;
    Struct.bc_sets   = cfg.bc_sets;
    Struct.build();
end
