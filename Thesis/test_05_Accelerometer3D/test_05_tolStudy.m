%% =====================================================================
%  TOLERANCE STUDY - is RelTol 1e-7 as good as 1e-8?
%
%  Runs ONE reduced model twice, changing nothing but the integration
%  tolerance, and measures how far the two transients drift apart. The point
%  is to decide whether the FOM reference can be run at 1e-7 instead of 1e-8,
%  which on this model is the difference between a day and several.
%
%  WHY A ROM AND NOT THE FOM. The question is about the FOM, but asking it
%  directly costs the very days the study is meant to save. A ROM integrated
%  by the same solver, on the same contact, with the same damping, is a proxy:
%  it shares the stiff contact events that set the step size, and differs only
%  in how much high-frequency content it carries. That makes it a
%  NON-CONSERVATIVE proxy - the FOM keeps more high-frequency content than any
%  ROM, so it can be MORE tolerance-sensitive, not less. Read a pass here as
%  "probably safe", never as proof, and keep a margin.
%
%  WHY IT MATTERS NOW. The earlier 1e-7 vs 1e-8 comparison was made with the
%  old damping, which put zeta = 1 already at 30 MHz and so annihilated the
%  high-frequency content the contact generates. The current damping moves
%  that threshold to about 930 MHz, so that content now survives and has to be
%  integrated. The old equivalence does not automatically carry over, which is
%  exactly what this script re-checks.
%
%  The configuration is INHERITED from a run_config.mat rather than retyped,
%  so the study cannot silently drift from the run it is meant to inform.
%
%  See also RUN_ROM_SWEEP, TRANSIENTSOLVERODE, ROM_TRACKING_TIME.
% =====================================================================
clear; close all; clc;

%% --- 0. PATHS ---------------------------------------------------------
test_root = fileparts(mfilename('fullpath'));
addpath(fullfile(test_root, '..', 'Src'));
addpath(fullfile(test_root, 'mesh'));
if isempty(which('KirchoffMaterial'))
    run(fullfile(test_root, '..', '..', 'startup.m'));
end
cd(test_root);

%% --- 1. STUDY CONFIGURATION -------------------------------------------

% Where the model configuration comes from. '' picks the NEWEST
% results/*/run_config.mat, which is the run currently being questioned.
% Set it by hand to study a specific past run instead.
study.config_from = '';

% The proxy model. MT is the cheapest basis that still sees the contact, so
% it is the natural probe; raise phi if you want the proxy to carry more
% high-frequency content and therefore be a harsher test.
study.model = 'MT';
study.phi   = 500;

% Tolerances to compare. The TIGHTEST is used as the reference.
study.reltols = [1e-8, 1e-7];

% [] inherits cfg.tmax. Shorten it for a first look: the drift between two
% tolerances grows with time, so a short window is optimistic by construction
% and only ever proves a problem, never the absence of one.
study.tmax = [];

% Thresholds for the tracking time, in per cent. Much tighter than the ones
% used for ROM error, because here we are measuring the integrator against
% itself and expect the difference to be small.
study.thresholds = [1e-3, 1e-2, 1e-1, 1];
study.hold_frac  = 0.01;

% The smallest ROM error the thesis actually cares about resolving. The
% verdict below asks whether the tolerance noise sits at least a decade under
% it; if it does not, the looser tolerance is measuring the integrator.
study.gre_of_interest = 0.1;   % [%]

% THE WINDOW THE VERDICT IS JUDGED ON, and the single most important knob here.
% A vibro-impact response is chaotic: two runs of the SAME model at different
% tolerances diverge exponentially and the error SATURATES around 100%. Judging
% a tolerance on the full-window GRE therefore measures the decorrelation, not
% the integrator, and condemns every tolerance including the good ones.
%
% Measured on this model with the current damping: the error grows by one
% decade every ~4 us (lambda = 6e5 1/s), so 1e-7 against 1e-8 sits at 4e-3% at
% 10 us, 0.9% at 20 us and saturates past 30 us. Under the old damping the same
% rate was 9.6e4 1/s, one decade per 24 us: the high-frequency content the new
% damping no longer annihilates is what accelerated it.
%
% The default matches the 25% window the post-processing already uses as its
% reference metric. Shorten it to be stricter about what the reference is
% trusted for; lengthening it past the decorrelation horizon is meaningless.
study.window = 25e-6;          % [s]

%% --- 2. INHERIT THE MODEL CONFIGURATION -------------------------------
if isempty(study.config_from)
    d = dir(fullfile('results', '*', 'run_config.mat'));
    if isempty(d)
        error('TOL:NoConfig', ...
            ['No results/*/run_config.mat found. Run the main at least once, ' ...
             'or set study.config_from by hand.']);
    end
    [~, newest] = max([d.datenum]);
    study.config_from = fullfile(d(newest).folder, d(newest).name);
end
fprintf('Inheriting the configuration from:\n  %s\n\n', study.config_from);

loaded = load(study.config_from, 'cfg');
cfg    = loaded.cfg;

if ~isempty(study.tmax), cfg.tmax = study.tmax; end

fprintf('Damping: Q_freq = %s Hz | Q = %s\n', ...
    mat2str(cfg.Q_freq), mat2str(cfg.array_QFactor));
fprintf('tmax = %.3e s | dt = %.3e s | k_mult = %g\n\n', ...
    cfg.tmax, cfg.dt, cfg.array_k_mult(1));

%% --- 3. MODEL ---------------------------------------------------------
Struct = FeStructure();
Struct.mesh_file    = cfg.mesh_file;
Struct.element_type = cfg.element_type;
Struct.E            = cfg.E;
Struct.nu           = cfg.nu;
Struct.rho          = cfg.rho;
Struct.set_specs    = cfg.node_sets;
Struct.bc_sets      = cfg.bc_sets;
Struct.build();
Struct.compute_eigenmodes(study.phi);

contact = build_interfaces(Struct, cfg);
shock   = shock_forcing(Struct, cfg);
labels  = fieldnames(contact.Interfaces)';

%% --- 4. ONE BASIS, BUILT ONCE -----------------------------------------
% Built once and reused, so the ONLY difference between the two runs is the
% tolerance. Rebuilding per run would be deterministic anyway, but this makes
% that independent of any future change inside build_rom.
[Q_pairs, f_anchor] = damping_spec(cfg);
if size(Q_pairs, 2) > 1
    warning('TOL:QSweep', ...
        'The inherited cfg defines %d damping cases; the study uses the first.', ...
        size(Q_pairs, 2));
end
[~, alpha_ray, beta_ray] = Struct.compute_rayleigh_damping( ...
    Q_pairs(1,1), Q_pairs(2,1), f_anchor);
rayleigh = struct('alpha', alpha_ray, 'beta', beta_ray);

k_contact = contact.k_base * cfg.array_k_mult(1);

fprintf('\nBuilding the %s basis (phi = %d)...\n', study.model, study.phi);
[~, Pc, Mr, Kr, Cr] = build_rom(Struct, study.model, study.phi, ...
    contact, k_contact, rayleigh);

N_run    = contact.N * Pc;
F_r      = Pc' * shock.F_spatial;
F_handle = @(tt) F_r * shock.profile(tt);
z0       = zeros(size(Kr, 1), 1);

%% --- 5. ONE INTEGRATION PER TOLERANCE ---------------------------------
n_tol = numel(study.reltols);
res(n_tol) = struct('reltol', [], 't', [], 'y', [], 'cpu', []);

for i = 1:n_tol
    fprintf('\n========== RelTol %.0e ==========\n', study.reltols(i));
    solver = TransientSolverOde(Mr, Kr, Cr);
    tic;
    [t_i, q_i] = solver.solve(cfg.tmax, cfg.dt, z0, z0, F_handle, ...
        'ContactOperator', N_run, ...
        'ContactGap',      contact.gaps, ...
        'ContactPenalty',  k_contact, ...
        'Label',           sprintf('%s%d @ %.0e', study.model, study.phi, study.reltols(i)), ...
        'Eref',            shock.Eref, ...
        'RelTol',          study.reltols(i), ...
        'OutputTimes',     shock.t_out);
    res(i).cpu    = toc;
    res(i).reltol = study.reltols(i);
    res(i).t      = t_i;
    res(i).y      = extract_contact_response(Struct, contact.Interfaces, ...
                                             labels, Pc * q_i);
end

%% --- 6. COMPARISON ----------------------------------------------------
% All contact DOFs concatenated, normal component: the same quantity the
% post-processing compares ROMs on, so the numbers here are on the same scale
% as the GRE in the summary tables.
cat_normal = @(y) cell2mat(cellfun(@(f) y.(f).normal, labels(:), ...
                                   'UniformOutput', false));

[~, i_ref] = min([res.reltol]);
Y    = arrayfun(@(r) cat_normal(r.y), res, 'UniformOutput', false);
Yref = Y{i_ref};
tref = res(i_ref).t;
ref_scale = max(sqrt(sum(Yref.^2, 1)));

fprintf('\n=========================================\n');
fprintf('  Tolerance study: %s phi %d, tmax %.3e s\n', ...
    study.model, study.phi, cfg.tmax);
fprintf('  Reference: RelTol %.0e\n', res(i_ref).reltol);
fprintf('=========================================\n\n');

% Both are reported, because they answer different questions. gre_full is the
% honest description of the whole run and SATURATES on the chaotic
% decorrelation; gre_win is the one the verdict uses, restricted to the window
% where a pointwise comparison still carries information.
i_win = tref <= study.window;

fprintf('%-10s %9s %9s %14s %14s\n', ...
    'RelTol', 'cpu [s]', 'speedup', 'GRE full', sprintf('GRE %.0f us', 1e6*study.window));
for i = 1:n_tol
    D = Y{i} - Yref;
    res(i).gre      = norm(D, 'fro') / norm(Yref, 'fro') * 100;
    res(i).gre_win  = norm(D(:,i_win), 'fro') / norm(Yref(:,i_win), 'fro') * 100;
    res(i).err_t    = sqrt(sum(D.^2, 1)) ./ (ref_scale + eps) * 100;
    fprintf('%-10.0e %9.1f %8.2fx %12.3e %% %12.3e %%\n', res(i).reltol, res(i).cpu, ...
        res(i_ref).cpu / res(i).cpu, res(i).gre, res(i).gre_win);
end

fprintf('\n--- Tracking time of each looser tolerance against the reference ---\n');
for i = 1:n_tol
    if i == i_ref, continue; end
    % Verbose, or the horizons are computed and silently kept in the struct -
    % which is exactly what happened on the first run of this study, leaving
    % the only informative number out of the log.
    res(i).trk = rom_tracking_time(tref, res(i).err_t, study.thresholds, ...
        'HoldFraction', study.hold_frac, 'Verbose', true, ...
        'Label', sprintf('RelTol %.0e', res(i).reltol));
end

%% --- 7. VERDICT -------------------------------------------------------
fprintf('\n=========================================\n');
fprintf('  VERDICT\n');
fprintf('=========================================\n');
% The verdict is judged on gre_win, NOT on the full-window GRE. The latter
% saturates on the chaotic decorrelation of two runs of the same model, so it
% would condemn every tolerance and tell us nothing about the integrator.
for i = 1:n_tol
    if i == i_ref, continue; end
    margin = study.gre_of_interest / res(i).gre_win;
    fprintf('\nRelTol %.0e vs %.0e, judged over the first %.1f us:\n', ...
        res(i).reltol, res(i_ref).reltol, 1e6*study.window);
    fprintf('  tolerance noise in window  %.3e %%\n', res(i).gre_win);
    fprintf('  ROM error of interest      %.3e %%\n', study.gre_of_interest);
    fprintf('  margin                     %.1fx\n', margin);
    fprintf('  (full-window GRE           %.3e %% - saturated, not a verdict)\n', res(i).gre);
    fprintf('  cpu saved                  %.1f%% (%.0f s of %.0f s)\n', ...
        100*(1 - res(i).cpu/res(i_ref).cpu), ...
        res(i_ref).cpu - res(i).cpu, res(i_ref).cpu);
    if margin >= 10
        fprintf(['  => PROBABLY SAFE over this window. The noise sits more than a\n' ...
                 '     decade under the smallest error of interest. The FOM carries\n' ...
                 '     more high-frequency content than this proxy and so diverges\n' ...
                 '     FASTER: keep the margin rather than spending it.\n']);
    elseif margin >= 3
        fprintf(['  => MARGINAL. The noise is within a decade of the errors being\n' ...
                 '     measured. Usable for ranking methods, NOT for quoting an\n' ...
                 '     absolute accuracy near %.2g%%. Shorten study.window or keep %.0e.\n'], ...
                 study.gre_of_interest, res(i_ref).reltol);
    else
        fprintf(['  => NOT SAFE over this window. Either keep %.0e, or shorten\n' ...
                 '     study.window to the horizon reported above where the error\n' ...
                 '     is still %.2g%%.\n'], res(i_ref).reltol, study.gre_of_interest);
    end
end

% The decorrelation rate itself, which is a property of the MODEL and not of
% the tolerances: it bounds how long ANY pointwise FOM-to-ROM comparison can
% carry information, for every method in the sweep.
for i = 1:n_tol
    if i == i_ref, continue; end
    e   = res(i).err_t;
    sel = e > 0 & tref > 0 & e < 10;      % exponential stretch, before saturation
    if nnz(sel) > 10
        p      = polyfit(tref(sel), log(e(sel)), 1);
        lambda = p(1);
        % The saturation instant is MEASURED, not derived from lambda: the time
        % to decorrelate depends on where the error starts, not on the rate
        % alone, so any closed-form "horizon" built from lambda by itself is
        % wrong. Here it would have read 7 us against a measured 30.
        j_sat = find(e >= 50, 1);
        if isempty(j_sat)
            sat = sprintf('not reached within %.0f us', 1e6*tref(end));
        else
            sat = sprintf('%.0f us', 1e6*tref(j_sat));
        end
        fprintf(['\nDecorrelation of two runs of the SAME model:\n' ...
                 '  lambda = %.2e 1/s, i.e. one decade of error every %.1f us\n' ...
                 '  50%% error reached at %s.\n' ...
                 '  Past that instant no pointwise FOM-to-ROM comparison carries\n' ...
                 '  information, for ANY method: the window, not the tolerance,\n' ...
                 '  is what bounds the whole study.\n'], ...
            lambda, 1e6*log(10)/lambda, sat);
    end
end

%% --- 8. FIGURES -------------------------------------------------------
% The busiest contact DOF: where a tolerance difference would show first.
[~, i_node] = max(max(abs(Yref), [], 2));

figure('Name', 'Tolerance study', 'Position', [100 100 900 700]);

subplot(2,1,1); hold on; grid on;
for i = 1:n_tol
    plot(1e6*res(i).t, Y{i}(i_node,:), 'LineWidth', 1.2, ...
        'DisplayName', sprintf('RelTol %.0e', res(i).reltol));
end
xlabel('t [\mus]'); ylabel('normal displacement [\mum]');
title(sprintf('%s \\phi%d - busiest contact DOF (#%d)', ...
    study.model, study.phi, i_node));
legend('Location', 'best');

subplot(2,1,2); hold on; grid on;
for i = 1:n_tol
    if i == i_ref, continue; end
    plot(1e6*tref, max(res(i).err_t, eps), 'LineWidth', 1.2, ...
        'DisplayName', sprintf('RelTol %.0e', res(i).reltol));
end
for th = study.thresholds
    yline(th, '--', sprintf('%g%%', th), 'Color', [0.6 0.6 0.6], ...
        'HandleVisibility', 'off');
end
set(gca, 'YScale', 'log');
xlabel('t [\mus]'); ylabel('difference vs reference [%]');
title('Tolerance-induced difference (the noise floor of the reference)');
legend('Location', 'southeast');

%% --- 9. SAVE ----------------------------------------------------------
out_dir = fileparts(study.config_from);
stem    = fullfile(out_dir, sprintf('tolStudy_%s%d', study.model, study.phi));

% Exported explicitly: this study is meant to be launched unattended, and when
% MATLAB runs headless (matlab -batch) an unsaved figure is simply lost.
savefig(gcf, [stem '.fig']);
exportgraphics(gcf, [stem '.png'], 'Resolution', 150);

save([stem '.mat'], 'study', 'res', 'cfg', '-v7.3');
fprintf('\nSaved %s .mat / .png / .fig\n', stem);
