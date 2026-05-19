%% --- Linear Beam Model  ---
clear; close all; clc;

beam = BeamModel();

% Spring Setup
k_spring = 35025;       % N/m
clearance = 0.394e-3;   % m
spring = LocalNonlinearity(k_spring, clearance, beam.nl_dof);

% ROM Construction
num_linear_modes = 11;
use_MC = 1;
rom = RomBuilder(beam.Assembly, num_linear_modes, beam.nl_dof, k_spring, use_MC);
rom.build();

m = size(rom.P, 2);
BeamReducedAssembly = ReducedAssembly(beam.Mesh, rom.P);
BeamReducedAssembly.DATA.M = BeamReducedAssembly.mass_matrix();
BeamReducedAssembly.DATA.K = BeamReducedAssembly.stiffness_matrix();
BeamReducedAssembly.DATA.C = zeros(m, m);
BeamReducedAssembly.DATA.D = zeros(m, m);

%% --- NONLINEAR NORMAL MODES (Shooting) ---

log10a_start = -4.5;
log10a_end   = -1;

figure('Name', 'NNM with Shooting', 'Color', 'w');
hold on; grid on;

nnm_solver = NNMContinuation_new(BeamReducedAssembly, spring, rom.P);

opts = struct();

% --- KEY PARAMETERS ---
% ds controls the base arclength step.
% A value of 1e-3 is 10x smaller than before (was 1e-2), which gives the
% corrector more room to converge near tongue bifurcation points.
opts.ds      = 1e-3;

% dsmax is now 5*ds instead of 100*ds (the old implicit default).
% This prevents the predictor from jumping across tongue loops.
opts.dsmax   = 5e-3;

% dsmin: allow very small steps so the solver squeezes through
% the tight turns at the tip of each tongue instead of giving up.
opts.dsmin   = 1e-9;

% Ntd = 1000: double the integration points per period.
% Impact events have steep velocity gradients; 500 pts is often
% insufficient to accurately evaluate the monodromy matrix,
% which corrupts the Newton-Raphson Jacobian near resonances.
opts.Ntd     = 1000;

opts.stepmax = 40000;

nnm_solver.solve(1, log10a_start, log10a_end, opts);

nnm_solver.plot_backbone('b-', 'ROM Shooting (Mode 1)');

set(gca, 'XScale', 'log');
title('NNM Backbone — Shooting Method');
xlabel('Energy [J]');
ylabel('Frequency [Hz]');
legend('Location', 'northwest');
hold off;

%% --- Frequency vs Amplitude diagnostic plot ---
% Useful to see WHERE the solver went: if it followed a tongue,
% the frequency will oscillate up/down as a function of the
% continuation parameter (log10a). This is physically expected.
figure('Name', 'Freq vs Continuation Amplitude', 'Color', 'w');
hold on; grid on;
nnm_solver.plot_frequency_vs_amplitude('b-', 'Mode 1');
hold off;

%% --- Time history at a chosen energy level ---
energy_level = 1e-4; % [J] — adjust to a point on the backbone
x_coords = beam.Assembly.Mesh.nodes(:, 1);
nnm_solver.plot_time_history_and_shape(energy_level, x_coords);

%% --- Transient verification ---
A_mag          = 165.94;   % [N]
duration_pulse = 0.001;    % [s]
tmax           = 1;        % [s]
h              = 5e-6;     % [s]

transient_sim = TransientSimulator(beam.Assembly, BeamReducedAssembly, ...
                    spring, rom.P, beam.M_full, beam.nl_dof);

transient_sim.solve_rom(A_mag, duration_pulse, tmax, h);
transient_sim.solve_fom(A_mag, duration_pulse, tmax, h);
transient_sim.print_summary();
transient_sim.plot_comparison(A_mag, tmax);
transient_sim.save_results(A_mag);

%% --- Frequency Analysis ---
rom.display_rom_frequencies();
freq_analyzer = FrequencyAnalyzer(transient_sim, A_mag, duration_pulse);
freq_analyzer.plot_results(1450);
