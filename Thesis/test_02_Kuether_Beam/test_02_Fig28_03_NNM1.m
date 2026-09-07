%% Reproduction of Fig. 28.3 - Kuether, Brake & Allen, IMAC 2014
%  "Evaluating Convergence of Reduced Order Models Using Nonlinear Normal Modes"
%
%  First NNM of the simply supported beam with a unilateral contact spring at
%  midspan, for the four reduction bases of the paper:
%      ROM [1], ROM [1 MC], ROM [1 3 5 MC], ROM [1 3 5 7 9 11 MC]
%
%  Method: harmonic balance (alternating frequency/time) plus pseudo-arclength
%  continuation. The paper uses shooting instead, but a fixed step Newmark
%  shooting residual is attracted to spurious grazing orbits as soon as the
%  spring engages (midspan amplitude equal to the gap, contact lasting less
%  than 1 % of the period): those orbits satisfy the discrete residual while
%  their true periodicity error is of order 100 %. Harmonic balance does not
%  suffer from this, and reproduces the frequencies quoted in the paper to
%  better than 1 %.
%
%  Number of harmonics: the internal resonances of these ROMs are 1:n with a
%  higher ROM mode (339 Hz for [1 MC], 314 Hz and 871 Hz for the larger bases),
%  and the solution can only represent one once H reaches n. That drives the two
%  phases of NNMHarmonicBalance.solve_with_tongues: the backbone is computed
%  with few harmonics, so that no tongue exists to be caught, and each tongue is
%  then traced by a dedicated run with many harmonics started just below its
%  crossing (see that method for details).
%
%  Each basis is saved to Results/Fig28_03_<tag>.mat as soon as it is done, so
%  the figure can be redrawn from disk without recomputing.

clear; close all; clc;

here = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(here, 'Tools')));
resdir = fullfile(here, 'Results');

IN_LBF = 0.1129848;        % 1 in-lbf = 0.1129848 J

% --- Paper data (Table 28.1) converted to SI ---
k_spring  = 35025;         % N/m, 200 lbf/in
clearance = 0.394e-3;      % m,   0.0155 in

% --- Reduction bases of Fig. 28.3 ---
% nModes is the number of eigenvectors extracted; RomBuilder keeps the odd
% ones only (the even modes have a node at the impact location).
% H_bb is the number of harmonics used for the backbone: it must stay below the
% lowest resonance order the branch can meet, otherwise the continuation turns
% into a tongue (the second ROM mode sits at about 9 times the first, so H = 5
% is safe). H_hi is used for the tongue runs; the highest order it can resolve
% is H_hi itself. ROM [1] has a single DOF, hence no internal resonance at all,
% and can use more harmonics for the backbone.
cfg(1) = struct('tag','1',         'label','ROM [1]',               'nModes', 1, 'useMC', false, 'color',[0.8 0 0.8],  'ls',':',  'H_bb', 9, 'H_hi', 9,  'N', 2048);
cfg(2) = struct('tag','1MC',       'label','ROM [1 MC]',            'nModes', 1, 'useMC', true,  'color',[0 0.6 0.2],  'ls','--', 'H_bb', 5, 'H_hi', 16, 'N', 2048);
cfg(3) = struct('tag','135MC',     'label','ROM [1 3 5 MC]',        'nModes', 5, 'useMC', true,  'color',[0.85 0.2 0], 'ls','-.', 'H_bb', 5, 'H_hi', 16, 'N', 2048);
cfg(4) = struct('tag','1357911MC', 'label','ROM [1 3 5 7 9 11 MC]', 'nModes',11, 'useMC', true,  'color',[0 0.3 0.9],  'ls','-',  'H_bb', 5, 'H_hi', 16, 'N', 2048);

% Continuation range, as log10 of the modal amplitude used by NLvib
log10a_start = -6;
log10a_end   =  0;

% Set to a subset (e.g. 1:2) to recompute only part of the bases
run_list = 1:numel(cfg);

%% Common model
beam   = BeamModel();
spring = LocalNonlinearity(k_spring, clearance, beam.nl_dof);

%% Continuation, one basis at a time
for ii = run_list
    c = cfg(ii);
    fprintf('\n================ %s (H = %d / %d) ================\n', ...
            c.label, c.H_bb, c.H_hi);

    rom = RomBuilder(beam.Assembly, c.nModes, beam.nl_dof, k_spring, c.useMC);
    rom.build();
    m = size(rom.P, 2);

    RA = ReducedAssembly(beam.Mesh, rom.P);
    RA.DATA.M = RA.mass_matrix();
    RA.DATA.K = RA.stiffness_matrix();
    RA.DATA.C = zeros(m);
    RA.DATA.D = zeros(m);

    % Frequencies of the ROM: the tongues sit at f_j/n
    f_rom = sort(sqrt(eig(full(RA.DATA.K), full(RA.DATA.M)))) / (2*pi);

    hb = NNMHarmonicBalance(RA, spring, rom.P);

    t0 = tic;
    % solve_and_continue prints one line per step: swallow it, keep ours
    txt = evalc(['hb.solve_with_tongues(1, log10a_start, log10a_end, ' ...
                 'c.H_bb, c.H_hi, c.N)']);
    disp(strjoin(regexp(txt, '^\s+(backbone|tongue).*$', 'match', 'lineanchors'), newline));
    fprintf('%d branches in %.1f min\n', numel(hb.segments), toc(t0)/60);

    seg = hb.segments;
    for is = 1:numel(seg)
        seg(is).energy_in_lbf = seg(is).energies / IN_LBF;
    end
    S = struct('label', c.label, 'tag', c.tag, 'color', c.color, 'ls', c.ls, ...
               'H_bb', c.H_bb, 'H_hi', c.H_hi, 'segments', seg, 'f_rom', f_rom); %#ok<NASGU>
    save(fullfile(resdir, ['Fig28_03_' c.tag '.mat']), '-struct', 'S');

    % Backbone frequencies at the energy levels quoted in the paper
    Ebb = hb.energies / IN_LBF;
    for Etar = [3e-4 2e-3 5e-2 0.76 5]
        [d, idx] = min(abs(log10(Ebb) - log10(Etar)));
        if d < 0.1
            fprintf('   E = %8.1e in-lbf -> f = %6.2f Hz\n', Ebb(idx), hb.frequencies(idx));
        end
    end
end

%% Frequency-energy plot (same axes as Fig. 28.3)
figure('Name', 'Fig. 28.3 - First NNM', 'Color', 'w');
hold on; grid on; box on;

% Richest basis first, so that the coarser ones stay visible on top of it: the
% 4 and 7 mode backbones lie on each other, exactly the point of the figure.
hleg = gobjects(1, numel(cfg));
for ii = numel(cfg):-1:1
    f = fullfile(resdir, ['Fig28_03_' cfg(ii).tag '.mat']);
    if ~isfile(f), continue; end
    S = load(f);
    % Backbone: thick line, the only one in the legend. Tongues: thin lines of
    % the same colour, drawn from a few points before they leave the backbone so
    % that they read as spikes attached to it, as in the paper.
    hleg(ii) = semilogx(S.segments(1).energy_in_lbf, S.segments(1).frequencies, ...
                        S.ls, 'Color', S.color, 'LineWidth', 1.8, ...
                        'DisplayName', S.label);
    for is = 2:numel(S.segments)
        s = S.segments(is);
        if ~tongue_is_useful(s, S.segments(1).frequencies(1)), continue; end
        k = max(1, s.i_exit - 20);
        semilogx(s.energy_in_lbf(k:end), s.frequencies(k:end), '-', ...
                 'Color', S.color, 'LineWidth', 0.7, 'HandleVisibility', 'off');
    end
end

set(gca, 'XScale', 'log');
xlim([1e-6, 1e2]);
ylim([34, 65]);
xlabel('Energy [in-lbf]');
ylabel('Frequency [Hz]');
title('First NNM of the beam with contact - Kuether et al., Fig. 28.3');
legend(hleg(isgraphics(hleg)), 'Location', 'northwest');

% Onset of contact quoted in the paper
xline(3.0e-4, 'k:', 'linear \rightarrow contact (paper: 3.0e-4)', ...
      'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');

exportgraphics(gcf, fullfile(resdir, 'Fig28_03_reproduction.png'), 'Resolution', 200);

%% Identification of the tongues of the highest fidelity ROM
%  Every tongue sits at f_j/n, where f_j is a ROM natural frequency and n the
%  order of the internal resonance. Drawn on top of the branch, these lines
%  say which mode each tongue belongs to.
ftag = fullfile(resdir, 'Fig28_03_1357911MC.mat');
if isfile(ftag)
    S = load(ftag);
    figure('Name', 'Internal resonances', 'Color', 'w');
    hold on; grid on; box on;
    semilogx(S.segments(1).energy_in_lbf, S.segments(1).frequencies, '-', ...
             'Color', S.color, 'LineWidth', 1.6, 'HandleVisibility', 'off');
    for is = 2:numel(S.segments)
        s = S.segments(is);
        if ~tongue_is_useful(s, S.segments(1).frequencies(1)), continue; end
        k = max(1, s.i_exit - 20);
        semilogx(s.energy_in_lbf(k:end), s.frequencies(k:end), '-', ...
                 'Color', S.color, 'LineWidth', 0.8, 'HandleVisibility', 'off');
    end
    for j = 2:numel(S.f_rom)
        for n = 2:S.H_hi
            fr = S.f_rom(j) / n;
            if fr > 35 && fr < 65
                yline(fr, 'k:', sprintf('mode %d, 1:%d', j, n), ...
                      'FontSize', 7, 'HandleVisibility', 'off');
            end
        end
    end
    set(gca, 'XScale', 'log');
    xlim([1e-4, 1e2]); ylim([34, 65]);
    xlabel('Energy [in-lbf]'); ylabel('Frequency [Hz]');
    title([S.label ': tongues at f_j/n']);
    exportgraphics(gcf, fullfile(resdir, 'Fig28_03_tongues.png'), 'Resolution', 200);
end

%% Helpers
function ok = tongue_is_useful(s, f_linear)
% A tongue run is worth drawing only if it actually left the backbone and
% travelled along the resonance. Runs that turn immediately at the very onset
% of contact, where the crossing sits exactly on the linear frequency, carry no
% visible branch.
n = numel(s.frequencies);
ok = (n - s.i_exit) > 50 && s.frequencies(s.i_exit) > f_linear + 0.5;
end
