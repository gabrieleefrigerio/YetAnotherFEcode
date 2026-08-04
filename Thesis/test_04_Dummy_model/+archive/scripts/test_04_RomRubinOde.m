%% =====================================================================
%  DUMMY MODEL for ROM BENCHMARKING
%  =====================================================================
clear; close all; clc;

%% =====================================================================
%  SETUP DIRECTORY
%  =====================================================================
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
save_dir = fullfile('results', timestamp);
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end
fprintf('Directory di salvataggio creata: %s\n', save_dir);

%% --- 1. FULL ORDER MODEL (FOM) ---
DummyStruct = DummyStructure('elementType', 'QUAD4');
DummyStruct.build();
% DummyStruct.plot_undeformed();

num_modi_plot = 5;
DummyStruct.compute_eigenmodes(num_modi_plot);
% for i = 1:num_modi_plot
%     DummyStruct.plot_mode(i);
% end

Mc = DummyStruct.AssemblyObj.constrain_matrix(DummyStruct.M);
Kc = DummyStruct.AssemblyObj.constrain_matrix(DummyStruct.K);
Cc = DummyStruct.AssemblyObj.constrain_matrix(DummyStruct.C);
n_dofs_fom = size(Mc, 1);

% Contact Parameters
gap_wall = 1.5e-6; 
k_contact = max(diag(Kc)) * 10; 
contact_dofs = DummyStruct.get_right_wall_dofs(); 

% --- Identificazione Nodi e GdL per Estrazione Dati ---
nodo_spia = length(DummyStruct.nodes);
nodi_parete = find(abs(DummyStruct.nodes(:,1) - max(DummyStruct.nodes(:,1))) < 1e-9);

dof_spia_global = (nodo_spia - 1) * DummyStruct.MeshObj.nDOFPerNode + 1;
dofs_parete_global = (nodi_parete - 1) * DummyStruct.MeshObj.nDOFPerNode + 1;

%% --- 2. TRANSIENT SETUP ---
% Simularion parameters
dt = 0.5e-8; 
tmax = (1e-3)/2; 

% Initial Conditions
q0 = zeros(n_dofs_fom, 1);
qd0 = zeros(n_dofs_fom, 1);

% Shock Parameters
shock_g = 1e5;                 % Acceleration [g]
shock_amp = shock_g * 9.81;    % Acceleration [m/s^2]
t_shock = 10e-7;               % Impulse duratrion [s]

% Building Shock
n_dof_per_node = DummyStruct.MeshObj.nDOFPerNode;
dir_vector = zeros(n_dofs_fom, 1);
dir_vector(1:n_dof_per_node:n_dofs_fom) = 1; 

F_spatial_fom = Mc * dir_vector;
F_fom_handle = @(t) F_spatial_fom * shock_amp * sin(pi * t / t_shock) * (t <= t_shock);

%% --- 3. FOM TRANSIENT ---
solverFOM = TransientSolver(Mc, Kc, Cc);

% 1. Newmark (Default)
[t_fom_newmark, q_fom_newmark] = solverFOM.solve(tmax, dt, q0, qd0, F_fom_handle, ...
    'ContactTargetDOF', contact_dofs, ...
    'ContactGap', gap_wall, ...
    'ContactPenalty', k_contact, ...
    'ModelType', 'FOM');
y_fom_newmark = DummyStruct.AssemblyObj.unconstrain_vector(q_fom_newmark);

% --- Data Extraction and Saving ---
y_spia_fom_newmark = y_fom_newmark(dof_spia_global, :);
y_parete_fom_newmark = y_fom_newmark(dofs_parete_global, :);

save(fullfile(save_dir, 'FOM_Newmark.mat'), 't_fom_newmark', 'y_spia_fom_newmark', 'y_parete_fom_newmark');
fprintf('Salvato: FOM_Newmark.mat (Dimensioni ridotte)\n');
clear y_fom_newmark q_fom_newmark; 

%% --- 4. RUBIN REDUCED ORDER MODEL (ROM) ---
num_linear_modes = 20; % Number of linear modes

% ROM Building
rom_builder = RomRubin(DummyStruct, num_linear_modes, contact_dofs);
rom_builder.build();
rom_builder.display_rom_frequencies();

% Reduced Matrices
[Mr_rubin, Kr_rubin, Cr_rubin] = rom_builder.get_reduced_matrices();

Pc_rubin = rom_builder.Pc; 
n_rom_rubin = size(Pc_rubin, 2);

% Contact DoF Mapping
n_dofs_contatto = length(contact_dofs);
contact_dofs_rubin = 1:n_dofs_contatto; 

% Display Dimensions
fprintf('\n--- ROM Dimensions ---\n');
fprintf('GdL FOM: %d\n', size(Kc, 1));
fprintf('GdL ROM (Rubin): %d\n', n_rom_rubin);
fprintf('---------------------------\n');
%% --- 5. ROM TRANSIENT (ODE15s) ---
fprintf('\n=== Risoluzione ROM (Milman-Chu) ===\n');

% Initial Conditions
q0_rubin = zeros(n_rom_rubin, 1);
qd0_rubin = zeros(n_rom_rubin, 1);

% Forcing
F_rubin_handle = @(t) Pc_rubin' * F_fom_handle(t);

% 1. Istanziamo la classe passando le matrici RIDOTTE (Mr, Kr, Cr)
solverROM = TransientSolverOde(Mr_rubin, Kr_rubin, Cr_rubin);

% 2. Lanciamo il solve. Notare ModelType = 'MC' e l'aggiunta di 'ProjectionMatrix'
[t_rubin, q_rom_rubin] = solverROM.solve(tmax, dt, q0_rubin, qd0_rubin, F_rubin_handle, ...
    'ContactTargetDOF', contact_dofs_rubin, ...
    'ContactGap', gap_wall, ...
    'ContactPenalty', k_contact, ...
    'ModelType', 'Rubin', ...
    'ProjectionMatrix', Pc_rubin);

% ROM Solution Expansion
q_rubin_c = Pc_rubin * q_rom_rubin; 
y_rom_rubin = DummyStruct.AssemblyObj.unconstrain_vector(q_rubin_c);

% --- Data Extraction and Saving ---
y_spia_rom_rubin = y_rom_rubin(dof_spia_global, :);
y_parete_rom_rubin = y_rom_rubin(dofs_parete_global, :);

save(fullfile(save_dir, 'ROM_Rubin.mat'), 't_rubin', 'y_spia_rom_rubin', 'y_parete_rom_rubin');
fprintf('Salvato: ROM_Rubin_Light.mat (Dimensioni ridotte)\n');
clear y_rom_rubin q_rom_rubin q_rubin_c;
%% --- 6. FOM vs ROM ---
fprintf('\n=== Post-Processing and Benchmarking ===\n');

nodo_spia = length(DummyStruct.nodes);
nodi_parete = find(abs(DummyStruct.nodes(:,1) - max(DummyStruct.nodes(:,1))) < 1e-9);

% GdL globali per l'estrazione
dof_spia_global = (nodo_spia - 1) * DummyStruct.MeshObj.nDOFPerNode + 1;
dofs_parete_global = (nodi_parete - 1) * DummyStruct.MeshObj.nDOFPerNode + 1;

% =====================================================================
% CALCOLO DEL GLOBAL RELATIVE ERROR (GRE) SU NODO SPIA
% =====================================================================
% Poiché t_fom_newmark e t_rubin sono diversi a causa dell'ATS,
% interpoliamo il ROM sulla griglia temporale del FOM
y_spia_rom_interp = interp1(t_rubin, y_spia_rom_rubin, t_fom_ode', 'linear', 'extrap');

% Calcolo dell'errore globale in norma L2 relativa
errore_assoluto = y_spia_fom_ode - y_spia_rom_interp;
gre_spostamento = norm(errore_assoluto) / norm(y_spia_fom_ode);
gre_percentage = gre_spostamento * 100;

fprintf('Global Relative Error (Spostamento Nodo %d): %.4f%%\n', nodo_spia, gre_percentage);

% =====================================================================
% --- A. Confronto Cinematica (Nodo Spia) con GRE ---
% =====================================================================
figure('Name', 'Impact Transient FOM vs ROM', 'Color', 'w', 'Units', 'normalized', 'Position', [0.1 0.4 0.4 0.4]);
plot(t_fom_ode', y_spia_fom_ode * 1e6, 'b', 'LineWidth', 2, 'DisplayName', 'FOM (Newmark)');
hold on;

% Inseriamo il valore del GRE direttamente nella legenda del ROM
rom_legend_str = sprintf('ROM (Rubin Method) - GRE: %.3f%%', gre_percentage);
plot(t_rubin, y_spia_rom_rubin * 1e6, 'r--', 'LineWidth', 1.5, 'DisplayName', rom_legend_str);

yline(gap_wall * 1e6, 'k:', 'LineWidth', 2, 'DisplayName', 'Wall');
title(sprintf('Impact Transient - Node #%d', nodo_spia));
xlabel('Time [s]'); ylabel('Displacement X [\mum]');
legend('Location', 'best'); grid on;

% =====================================================================
% --- B. Confronto Forza di Contatto Totale ---
% =====================================================================
% pen_fom = max(0, y_parete_fom_newmark - gap_wall);
% forza_fom = sum(k_contact * pen_fom, 1);
% 
% pen_rom = max(0, y_parete_rom_rubin - gap_wall);
% forza_rom = sum(k_contact * pen_rom, 1);
% 
% figure('Name', 'Total Contact Force FOM vs ROM', 'Color', 'w', 'Units', 'normalized', 'Position', [0.5 0.4 0.4 0.4]);
% plot(t_fom_newmark, forza_fom * 1e6, 'b', 'LineWidth', 2, 'DisplayName', 'FOM (Newmark)');
% hold on;
% plot(t_rubin, forza_rom * 1e6, 'r--', 'LineWidth', 1.5, 'DisplayName', 'ROM (Rubin)');
% title('Total Contact Force');
% xlabel('Time [s]'); ylabel('Force [\muN]');
% legend('Location', 'best'); grid on;