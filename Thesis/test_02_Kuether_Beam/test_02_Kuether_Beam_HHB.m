
%% --- Linear Beam Model  ---
clear; close all; clc;
% 1. Create the Linear Beam Model
beam = BeamModel();
% beam.plot_mode(5);
% beam.display_linear_modes();

% 2. Spring Setup
k_spring = 35025;       % N/m
clearance = 0.394e-3;   % m
spring = LocalNonlinearity(k_spring, clearance, beam.nl_dof);

% 3. Reduced Order Model (ROM) Construction
num_linear_modes =1;
use_MC = 1;
rom = RomBuilder(beam.Assembly, num_linear_modes, beam.nl_dof, k_spring, use_MC);
rom.build();
% rom.plot_milman_chu_comparison();

m = size(rom.P, 2); % Number of degrees of freedom in the ROM
BeamReducedAssembly = ReducedAssembly(beam.Mesh, rom.P); % YaFEc class
BeamReducedAssembly.DATA.M = BeamReducedAssembly.mass_matrix();
BeamReducedAssembly.DATA.K = BeamReducedAssembly.stiffness_matrix();
BeamReducedAssembly.DATA.C = zeros(m, m);
BeamReducedAssembly.DATA.D = zeros(m, m);
%% --- NONLINEAR NORMAL MODES  ---

hb_solver = NNMHarmonicBalance(BeamReducedAssembly, spring, rom.P);
hb_solver.solve(1, -5, 10, 100, 1024);

figure('Name', 'NNM con Harmonic Balance', 'Color', 'w');
hb_solver.plot_backbone('b-', 'ROM (Harmonic Balance)');

ylim([30, 65]); 
title('Convergenza del Nonlinear Normal Mode (HB)');
legend('Location', 'northwest');

%% --- HARMONIC BALANCE CONVERGENCE (NUMBER OF HARMONICS) ---
% Define the harmonic values to test
H_values = [3, 7, 15, 31];
colors = lines(length(H_values));

% Conversion factor from in-lbf to Joules
in_lbf_to_joules = 0.1129848;

% Create the main figure
main_fig_hb = figure('Name', 'Harmonic Balance Convergence', 'Color', 'w');
hold on; grid on; box on;

% Loop over the number of harmonics
for i = 1:length(H_values)
    current_H = H_values(i);
    fprintf('\n--- RUN HB %d/%d (Harmonics = %d) ---\n', i, length(H_values), current_H);
    
    % Initialize the solver for this iteration
    hb_solver = NNMHarmonicBalance(BeamReducedAssembly, spring, rom.P);
    
    % solve(mode_idx, log10a_start, log10a_end, H, Ntd)
    % We use 1024 time samples (Ntd) for a highly accurate Alternating Frequency/Time (AFT) scheme
    hb_solver.solve(1, -5, 2, current_H, 1024);
    
    % Convert the extracted energies to Joules
    energies_J = hb_solver.energies * in_lbf_to_joules;
    
    % Plot the data using a logarithmic scale for the X-axis
    display_name = sprintf('HB (H = %d)', current_H);
    semilogx(energies_J, hb_solver.frequencies, ...
        '-', 'LineWidth', 2, 'Color', colors(i,:), 'DisplayName', display_name);
    
    drawnow; % Update the plot in real-time
end

% Final plot formatting
set(gca, 'XScale', 'log'); % Explicitly enforce logarithmic scale on the X-axis
ylim([30, 65]); 
title('NNM Convergence by Varying Harmonics (HB)', 'FontSize', 12);
xlabel('Energy [J]', 'FontSize', 11);
ylabel('Frequency [Hz]', 'FontSize', 11);
legend('Location', 'northwest', 'FontSize', 11);

fprintf('\nHarmonic Balance convergence analysis successfully completed!\n');