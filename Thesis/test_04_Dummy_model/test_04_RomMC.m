%% =====================================================================
%  DUMMY MODEL for ROM BENCHMARKING (MILMAN-CHU VERSION)
% =====================================================================
clear; close all; clc;

%% =====================================================================
%  SETUP DIRECTORY DI SALVATAGGIO
% =====================================================================
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

% --- Identificazione Nodi e GdL Globali per Estrazione Dati ---
target_node = length(DummyStruct.nodes);
nodo_spia = target_node; 
nodi_parete = find(abs(DummyStruct.nodes(:,1) - max(DummyStruct.nodes(:,1))) < 1e-9);

dof_spia_global = (nodo_spia - 1) * DummyStruct.MeshObj.nDOFPerNode + 1;
dofs_parete_global = (nodi_parete - 1) * DummyStruct.MeshObj.nDOFPerNode + 1;

%% --- 2. TRANSIENT SETUP ---
% Simulation parameters
dt = 1e-8; 
tmax = 1e-3; 

% Initial Conditions (Struttura pre-deflessa sul 1° modo)
q0 = 2e-6 * DummyStruct.AssemblyObj.constrain_vector(DummyStruct.mode_shapes(:,1));
qd0 = zeros(size(Kc, 1), 1);

% Forcing
F_spatial_c = DummyStruct.create_constrained_force_vector(target_node, 1);
F_fom_handle = @(t) F_spatial_c * (0 * sin(2 * pi * 2000 * t));

%% --- 3. FOM TRANSIENT ---
solverFOM = YafecTransientSolver(Mc, Kc, Cc);
[t_nl, q_nl] = solverFOM.solve(tmax, dt, q0, qd0, F_fom_handle, ...
    'ContactTargetDOF', contact_dofs, ...
    'ContactGap', gap_wall, ...
    'ContactPenalty', k_contact);

% Espansione nello spazio fisico completo
y_nl_full = DummyStruct.AssemblyObj.unconstrain_vector(q_nl);

% --- Estrazione Intelligente FOM ---
y_spia_fom = y_nl_full(dof_spia_global, :);
y_parete_fom = y_nl_full(dofs_parete_global, :);

% Salvataggio leggero su disco
save(fullfile(save_dir, 'FOM_MilmanChu_Light.mat'), 't_nl', 'y_spia_fom', 'y_parete_fom');
fprintf('Salvato: FOM_MilmanChu_Light.mat\n');

% Pulizia immediata della RAM
clear y_nl_full q_nl;

%% --- 3. REDUCED ORDER MODEL (ROM) ---
num_linear_modes = 20; % Numero di modi da includere
include_MC = 1;        % Arricchimento statico

% Creazione Base ROM
rom_builder = RomMC(DummyStruct, num_linear_modes, contact_dofs, k_contact, include_MC);
rom_builder.build();
rom_builder.display_rom_frequencies();

% Estrazione Matrici Ridotte
[Mr, Kr, Cr] = rom_builder.get_reduced_matrices();

% Estrazione Matrice di Proiezione Vincolata (Pc) per le condizioni iniziali
n_rom = size(rom_builder.P, 2);
Pc = zeros(size(Kc, 1), n_rom);
for i = 1:n_rom
    Pc(:, i) = DummyStruct.AssemblyObj.constrain_vector(rom_builder.P(:, i));
end

% Verifica Dimensionale
fprintf('--- Verifica Dimensionale ---\n');
fprintf('GdL FOM: %d\n', size(Kc, 1));
fprintf('GdL ROM: %d\n', size(Kr, 1));
fprintf('---------------------------\n');

%% --- 4. ROM TRANSIENT ---
fprintf('\n=== SEZIONE 4: Simulazione ROM (Benchmarking) ===\n');
fprintf('\n--- Avvio Metodo Milman-Chu ---\n');

% Proiezione Condizioni Iniziali (Corretto: usiamo q0_fom)
q0_mc = Pc \ q0; 
qd0_mc = zeros(size(Pc, 2), 1);

% Proiezione Forzante nel sottospazio MC
F_mc_handle = @(t) Pc' * F_fom_handle(t);

% Lancio Solutore ROM (Milman-Chu)
solverMC = TransientSolver(Mr, Kr, Cr);
[t_rom_mc, q_rom_mc] = solverMC.solve(tmax, dt, q0_mc, qd0_mc, F_mc_handle, ...
    'ContactTargetDOF', contact_dofs, ... 
    'ContactGap', gap_wall, ...
    'ContactPenalty', k_contact, ...
    'ProjectionMatrix', Pc, ...
    'ModelType', 'MC');                   

% Espansione nello spazio fisico completo
q_mc_c = Pc * q_rom_mc; 
y_rom_mc = DummyStruct.AssemblyObj.unconstrain_vector(q_mc_c);

% --- Estrazione Intelligente ROM ---
y_spia_rom_mc = y_rom_mc(dof_spia_global, :);
y_parete_rom_mc = y_rom_mc(dofs_parete_global, :);

% Salvataggio leggero su disco
save(fullfile(save_dir, 'ROM_MilmanChu_Light.mat'), 't_rom_mc', 'y_spia_rom_mc', 'y_parete_rom_mc');
fprintf('Salvato: ROM_MilmanChu_Light.mat\n');

% Pulizia della RAM
clear y_rom_mc q_rom_mc q_mc_c;
%% --- 5. FOM vs ROM ---
fprintf('\n=== SEZIONE 6: Post-Processing e Confronto ===\n');

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
y_spia_rom_interp = interp1(t_rubin, y_spia_rom_mc, t_fom_newmark, 'linear', 'extrap');

% Calcolo dell'errore globale in norma L2 relativa
errore_assoluto = y_spia_fom_newmark - y_spia_rom_interp;
gre_spostamento = norm(errore_assoluto) / norm(y_spia_fom_newmark);
gre_percentage = gre_spostamento * 100;

fprintf('Global Relative Error (Spostamento Nodo %d): %.4f%%\n', nodo_spia, gre_percentage);

% =====================================================================
% --- A. Confronto Cinematica (Nodo Spia) con GRE ---
% =====================================================================
figure('Name', 'Confronto Cinematica FOM vs ROM', 'Color', 'w', 'Units', 'normalized', 'Position', [0.1 0.4 0.4 0.4]);
plot(t_nl, y_spia_fom * 1e6, 'b', 'LineWidth', 2, 'DisplayName', 'FOM (Completo)');
hold on;

% Inseriamo il valore del GRE direttamente nella legenda del ROM Milman-Chu
rom_legend_str = sprintf('ROM (Milman-Chu) - GRE: %.3f%%', gre_percentage);
plot(t_rom_mc, y_spia_rom_mc * 1e6, 'r--', 'LineWidth', 1.5, 'DisplayName', rom_legend_str);

yline(gap_wall * 1e6, 'k:', 'LineWidth', 2, 'DisplayName', 'Parete');
title(sprintf('Cinematica di Impatto - Nodo %d', nodo_spia));
xlabel('Tempo [s]'); ylabel('Spostamento X [\mum]');
legend('Location', 'best'); grid on;
hold off;

% =====================================================================
% --- B. Confronto Forza di Contatto Totale ---
% =====================================================================
pen_fom = max(0, y_parete_fom - gap_wall);
forza_fom = sum(k_contact * pen_fom, 1);

pen_rom = max(0, y_parete_rom_mc - gap_wall);
forza_rom = sum(k_contact * pen_rom, 1);

figure('Name', 'Confronto Forza FOM vs ROM', 'Color', 'w', 'Units', 'normalized', 'Position', [0.5 0.4 0.4 0.4]);
plot(t_nl, forza_fom * 1e6, 'b', 'LineWidth', 2, 'DisplayName', 'FOM (Completo)');
hold on;
plot(t_rom_mc, forza_rom * 1e6, 'r--', 'LineWidth', 1.5, 'DisplayName', 'ROM (Milman-Chu)');

title('Forza Impulsiva Totale');
xlabel('Tempo [s]'); ylabel('Forza di Reazione [\muN]');
legend('Location', 'best'); grid on;
hold off;