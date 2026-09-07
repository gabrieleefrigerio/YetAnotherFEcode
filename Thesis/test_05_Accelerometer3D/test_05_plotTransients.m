%% =====================================================================
%  TRANSIENTS ONLY - no reference model needed
%
%  Draws the contact response of whatever result files a folder holds, FOM or
%  ROM or both, one tile per bumpstop. It asks no question about accuracy:
%  there is no GRE, no tracking time and no summary figure, because all three
%  need a reference. Use test_05_postProcessing for that.
%
%  The point of a separate script is that a run does not need a FOM to be
%  worth looking at. A ROM sweep alone, or a single exploratory run, can be
%  inspected here while the eight-hour reference is still going.
%
%  Everything about the model is read from run_config.mat, so this needs no
%  prior knowledge of the mesh, the stoppers or the shock direction.
% =====================================================================
% Set pp_results_dir before calling to drive this from another script without
% the folder dialog. It is consumed here and cleared with everything else, so
% running the script by hand ALWAYS asks.
if exist('pp_results_dir', 'var') && ~isempty(pp_results_dir)
    pp_dir = pp_results_dir;
else
    pp_dir = '';
end
if exist('pp_select_files', 'var') && ~isempty(pp_select_files)
    pp_sel = pp_select_files;
else
    pp_sel = [];
end
clearvars -except pp_dir pp_sel; close all; clc;

addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'Src'));

%% --- Options ----------------------------------------------------------

% What the vertical axis shows. 'physical' projects on the shock direction and
% puts each wall at its signed position, which is the only one of the three
% that does not send half the panels the wrong way. See contact_plot_quantity.
plot_quantity = 'physical';    % 'physical' | 'clearance' | 'normal'

% Spread of the whole face, as a shaded band behind the curves. First selected
% file only: one band per model would hide the curves.
show_band = true;

% Empty -> a dialog lists the files. Otherwise a cell of names or a pattern,
% e.g. {'ROM_MT_Phi500_Q1000_K0.003125.mat'} or 'ROM_MT_*'.
select_files = '';
if ~isempty(pp_sel), select_files = pp_sel; end

% Time window [us]. Empty -> the whole run.
t_window = [];

%% --- 1. Results folder ------------------------------------------------
if isempty(pp_dir)
    results_dir = uigetdir(fullfile(fileparts(mfilename('fullpath')), 'results'), ...
        'Select the results folder');
    if results_dir == 0, error('PT:NoFolder', 'No folder selected.'); end
else
    results_dir = pp_dir;
    if ~isfolder(results_dir)
        error('PT:NoFolder', 'pp_results_dir does not exist: %s', results_dir);
    end
end
fprintf('Folder: %s\n', results_dir);

cfg_file = fullfile(results_dir, 'run_config.mat');
if ~isfile(cfg_file)
    error('PT:NoConfig', 'run_config.mat is missing from %s.', results_dir);
end
R = load(cfg_file);

[faces, face_gap, face_nrm, d_ref] = interface_table(R);
n_faces = numel(faces);

%% --- 2. Which files to draw -------------------------------------------
listing = [dir(fullfile(results_dir, 'FOM_*.mat')); ...
           dir(fullfile(results_dir, 'ROM_*.mat'))];
if isempty(listing)
    error('PT:NoResults', 'No FOM_*.mat or ROM_*.mat in %s.', results_dir);
end
names = {listing.name};

if isempty(select_files)
    [pick, ok] = listdlg('ListString', names, 'SelectionMode', 'multiple', ...
        'Name', 'Transients', 'PromptString', 'Files to draw:', ...
        'ListSize', [420 320], 'InitialValue', 1:min(4, numel(names)));
    if ~ok, fprintf('Nothing selected.\n'); return; end
    sel = names(pick);
elseif ischar(select_files) || isstring(select_files)
    hit = dir(fullfile(results_dir, char(select_files)));
    sel = {hit.name};
    if isempty(sel), error('PT:NoMatch', 'Nothing matches %s.', char(select_files)); end
else
    sel = select_files;
end
fprintf('Drawing %d file(s).\n', numel(sel));

%% --- 3. Load ----------------------------------------------------------
D = struct('name', {}, 'label', {}, 't', {}, 'y', {});
for i = 1:numel(sel)
    S = load(fullfile(results_dir, sel{i}), 't', 'y_contact', 'model_tag', ...
        'n_modes', 'n_cc', 'model');
    if ~isfield(S, 'y_contact')
        warning('PT:NoResponse', '%s holds no y_contact, skipped.', sel{i});
        continue
    end
    % A label built from the file content rather than from its name, so a
    % renamed file still says what it is.
    if isfield(S, 'model_tag'),  lbl = S.model_tag;
    elseif isfield(S, 'model'),  lbl = S.model;
    else,                        lbl = strrep(sel{i}, '.mat', '');
    end
    if isfield(S, 'n_modes') && ~isempty(S.n_modes) && ~strcmpi(lbl, 'FOM')
        lbl = sprintf('%s phi=%d', lbl, S.n_modes);
    end
    if isfield(S, 'n_cc') && ~isempty(S.n_cc) && S.n_cc > 0
        lbl = sprintf('%s, n_cc=%d', lbl, S.n_cc);
    end
    D(end+1) = struct('name', sel{i}, 'label', lbl, ...
        't', S.t(:)', 'y', S.y_contact); %#ok<SAGROW>
end
if isempty(D), error('PT:Empty', 'None of the selected files holds a response.'); end

%% --- 4. Reference node per face ---------------------------------------
% The deepest-penetrating node governs the contact force, so it is the one
% worth following. Taken from the FIRST selected file: with no reference model
% there is nothing more authoritative, and a different node per model would
% make the curves incomparable.
ref_node = ones(1, n_faces);
for f = 1:n_faces
    un  = D(1).y.(faces{f}).normal;
    pen = max(un - face_gap(f), [], 2);
    if any(pen > 0), [~, ref_node(f)] = max(pen);
    else,            [~, ref_node(f)] = max(max(abs(un), [], 2));
    end
end

%% --- 5. Figure --------------------------------------------------------
n_rows = ceil(n_faces/2);
fig = figure('Name', 'Transients', 'NumberTitle', 'off', 'Color', 'w', ...
             'Position', [60, 40, 1200, 230*n_rows]);
tl  = tiledlayout(fig, n_rows, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
axs = gobjects(n_faces, 1);
colors = lines(numel(D));
vis_links = {};
h_wall = gobjects(n_faces, 1);

for f = 1:n_faces
    axs(f) = nexttile(tl);
    hold(axs(f), 'on'); grid(axs(f), 'on'); box(axs(f), 'on');

    [v1, wall_level, ylab] = contact_plot_quantity(D(1).y.(faces{f}), ...
        face_nrm{f}, face_gap(f), plot_quantity, d_ref);

    if show_band
        t_us = 1e6 * D(1).t(:);
        lo = min(v1, [], 1);  hi = max(v1, [], 1);
        fill(axs(f), [t_us; flipud(t_us)], [hi(:); flipud(lo(:))], ...
            [0.25 0.25 0.25], 'FaceAlpha', 0.12, 'EdgeColor', 'none', ...
            'HandleVisibility', 'off');
    end

    h_wall(f) = yline(axs(f), wall_level, 'r-.', 'LineWidth', 1.4);
    if f == 1, h_wall(f).DisplayName = 'Wall';
    else,      h_wall(f).HandleVisibility = 'off';
    end

    title(axs(f), sprintf('%s  -  node %d', faces{f}, ref_node(f)), ...
        'Interpreter', 'none');
    if f > n_faces - 2, xlabel(axs(f), 't [\mus]'); end
    ylabel(axs(f), ylab);
end

for i = 1:numel(D)
    h_this = gobjects(n_faces, 1);
    for f = 1:n_faces
        v = contact_plot_quantity(D(i).y.(faces{f}), face_nrm{f}, ...
            face_gap(f), plot_quantity, d_ref);
        if i == 1, sty = '-'; lw = 2; else, sty = '--'; lw = 1.4; end
        h_this(f) = plot(axs(f), 1e6*D(i).t, v(ref_node(f), :), sty, ...
            'Color', colors(i,:), 'LineWidth', lw);
        if f == 1
            h_this(f).DisplayName = D(i).label;
        else
            h_this(f).HandleVisibility = 'off';
        end
    end
    % Every copy of one model across the tiles switches as a single object,
    % otherwise the plot browser only hides the one carrying the legend entry.
    vis_links{end+1} = linkprop(h_this, 'Visible'); %#ok<SAGROW>
end
vis_links{end+1} = linkprop(h_wall, 'Visible');
% linkprop dies with its variable, so the links have to live on the figure.
setappdata(fig, 'VisibilityLinks', vis_links);

linkaxes(axs, 'x');            % zoom on one tile moves them all in time
if ~isempty(t_window), xlim(axs(1), t_window); end

legend(axs(1), 'Location', 'best', 'Interpreter', 'none');
if show_band
    band_note = sprintf(', band = face spread of %s', D(1).label);
else
    band_note = '';
end
sgtitle(fig, sprintf('Contact transients (%s%s)', plot_quantity, band_note), ...
    'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'none');

fprintf('Done: %d tiles, %d curves each.\n', n_faces, numel(D));
