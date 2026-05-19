%% --- REPLICATE MONJARAZ-TEC PAPER RESULTS ---
clear; close all; clc;

%% 1. Physical Model Setup
beam = FallingBarModel();
beam.display_linear_modes();

%% 2. Craig-Bampton ROM Setup (Massless, fixed-interface base)
num_linear_modes = 20; 
rom_cb = MasslessCBBuilder(beam.Assembly, num_linear_modes, beam.nl_dof);
rom_cb.build();
rom_cb.display_frequencies();

%% 3. Integrators (LCP)
sim_cb  = TransientSimulatorMCB(rom_cb);
% --- BENCHMARK PARAMETERS (Doyen et al.) ---
initial_height = 5.0;   % [m] 
gravity = -10.0;        % [m/s^2] 
tmax = 500;            % [s] 
dt = 1e-4;             % [s]

% Integration
fprintf('Running Massless CB simulation...\n');
sim_cb.solve_drop(initial_height, gravity, tmax, dt, beam.M_full);


%% 4. PLOT RESULTS
plotter = DoyenBenchmarkPlotter(sim_cb, initial_height, gravity);
plotter.plot_kinematics();
% plotter.plot_energies();