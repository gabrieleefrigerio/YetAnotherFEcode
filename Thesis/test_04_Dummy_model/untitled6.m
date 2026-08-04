%% =====================================================================
%  SENSITIVITY ANALYSIS: ode15s Tolerances (RelTol & AbsTol)
% =====================================================================
clear; close all; clc;

%% --- 1. SETUP PARAMETRI FISSI ---
fixed_phi = 50;
fixed_Q = 1000;
fixed_k_mult = 10;

shock_g = 1e5;
shock_amp = shock_g * 9.81;    
t_shock = 10e-7;               
dt = 0.5e-8;
tmax = (1e-3)/2;
gap_wall = 1.5e-6;

% Arrays per lo Sweep delle Tolleranze
array_RelTol = [1e-3, 1e-4, 1e-5, 1e-6];

array_AbsTol = [1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11];

%% --- 2. INIZIALIZZAZIONE STRUTTURA ---
fprintf('Setting up the model for Sensitivity Analysis...\n');
DummyStruct = AbaqusStructure();
DummyStruct.filename = 'DummyStructureAbaqus.inp'; 
DummyStruct.elementType = 'TRI3';           
DummyStruct.build();
DummyStruct.compute_eigenmodes(fixed_phi);
DummyStruct.compute_rayleigh_damping(fixed_Q, fixed_Q);

Mc = DummyStruct.AssemblyObj.constrain_matrix(DummyStruct.M);
Kc = DummyStruct.AssemblyObj.constrain_matrix(DummyStruct.K);
Cc = DummyStruct.AssemblyObj.constrain_matrix(DummyStruct.C);

n_dofs_fom = size(Mc, 1);
k_base = max(diag(Kc));
k_contact = k_base * fixed_k_mult;

contact_dofs = DummyStruct.get_contact_dofs(1); 
contact_nodes_global_X = (DummyStruct.contact_nodes - 1) * DummyStruct.MeshObj.nDOFPerNode + 1;

n_dof_per_node = DummyStruct.MeshObj.nDOFPerNode;
dir_vector = zeros(n_dofs_fom, 1);
dir_vector(1:n_dof_per_node:n_dofs_fom) = 1; 
F_spatial_fom = Mc * dir_vector;
F_fom_handle = @(t) F_spatial_fom * shock_amp * sin(pi * t / t_shock) * (t <= t_shock);

q0 = zeros(n_dofs_fom, 1); qd0 = zeros(n_dofs_fom, 1);

%% --- 3. GROUND TRUTH (FOM SUPER PRECISO) ---
fprintf('\n--- Computing Ground Truth (FOM con RelTol=1e-5, AbsTol=1e-10) ---\n');
solverTruth = TransientSolverOde_NEW(Mc, Kc, Cc);
[t_ref, q_truth] = solverTruth.solve(tmax, dt, q0, qd0, F_fom_handle, ...
    'ContactTargetDOF', contact_dofs, 'ContactGap', gap_wall, 'ContactPenalty', k_contact, ...
    'ModelType', 'FOM', 'RelTol', 1e-5, 'AbsTol', 1e-10);

y_truth_full = DummyStruct.AssemblyObj.unconstrain_vector(q_truth);
y_X_truth = y_truth_full(contact_nodes_global_X, :);

%% --- 4. COSTRUZIONE BASI ROM (OFFLINE) ---
fprintf('\n--- Pre-computing ROM Bases (Offline) ---\n');
% MC
rom_mc = RomMC(DummyStruct, fixed_phi, contact_dofs, k_contact, 1); rom_mc.build();
[Mr_mc, Kr_mc, Cr_mc] = rom_mc.get_reduced_matrices();
Pc_mc = sparse(size(Kc, 1), size(rom_mc.P, 2));
for i = 1:size(Pc_mc, 2), Pc_mc(:, i) = DummyStruct.AssemblyObj.constrain_vector(rom_mc.P(:, i)); end

% Rubin
rom_rubin = RomRubin(DummyStruct, fixed_phi, contact_dofs); rom_rubin.build();
[Mr_rubin, Kr_rubin, Cr_rubin] = rom_rubin.get_reduced_matrices(); Pc_rubin = rom_rubin.Pc;

% MCB
rom_mcb = RomMCB(DummyStruct, fixed_phi, contact_dofs); rom_mcb.build();
[Mr_mcb, Kr_mcb, Cr_mcb] = rom_mcb.get_reduced_matrices(); Pc_mcb = rom_mcb.Pc;

% MN
rom_mn = RomMN(DummyStruct, fixed_phi, contact_dofs); rom_mn.build();
[Mr_mn, Kr_mn, Cr_mn] = rom_mn.get_reduced_matrices(); Pc_mn = rom_mn.Pc;

rom_contact_dofs = 1:length(contact_dofs);

%% --- 5. SWEEP DELLE TOLLERANZE ---
% Inizializza array di salvataggio
results = struct();
models = {'FOM', 'MC', 'Rubin', 'MCB', 'MN'};
for m = 1:length(models)
    results.(models{m}).cpu = zeros(length(array_RelTol), length(array_AbsTol));
    results.(models{m}).gre = zeros(length(array_RelTol), length(array_AbsTol));
end

fprintf('\n=== INIZIO SWEEP TOLLERANZE ===\n');
for i = 1:length(array_RelTol)
    RT = array_RelTol(i);
    for j = 1:length(array_AbsTol)
        AT = array_AbsTol(j);
        fprintf('\n--- Test: RelTol = %g | AbsTol = %g ---\n', RT, AT);
        
        % 1. FOM Sweep
        solverFOM = TransientSolverOde_NEW(Mc, Kc, Cc);
        tic;
        [t_out, q_fom] = solverFOM.solve(tmax, dt, q0, qd0, F_fom_handle, ...
            'ContactTargetDOF', contact_dofs, 'ContactGap', gap_wall, 'ContactPenalty', k_contact, ...
            'ModelType', 'FOM', 'RelTol', RT, 'AbsTol', AT);
        results.FOM.cpu(i,j) = toc;
        y_fom_full = DummyStruct.AssemblyObj.unconstrain_vector(q_fom);
        y_X_fom = y_fom_full(contact_nodes_global_X, :);
        % Interpolazione necessaria se ode15s cambia step
        y_X_fom_int = interp1(t_out, y_X_fom', t_ref, 'linear', 'extrap')';
        results.FOM.gre(i,j) = (norm(y_X_truth(:) - y_X_fom_int(:)) / norm(y_X_truth(:))) * 100;
        
        % 2. MC
        q0_mc = Pc_mc \ q0; F_mc_handle = @(t) Pc_mc' * F_fom_handle(t);
        solverMC = TransientSolverOde_NEW(Mr_mc, Kr_mc, Cr_mc); tic;
        [t_out, q_rom] = solverMC.solve(tmax, dt, q0_mc, zeros(size(q0_mc)), F_mc_handle, ...
            'ContactTargetDOF', contact_dofs, 'ContactGap', gap_wall, 'ContactPenalty', k_contact, ...
            'ModelType', 'MC', 'ProjectionMatrix', Pc_mc, 'RelTol', RT, 'AbsTol', AT);
        results.MC.cpu(i,j) = toc;
        y_X_rom = DummyStruct.AssemblyObj.unconstrain_vector(Pc_mc * q_rom);
        y_X_rom = y_X_rom(contact_nodes_global_X, :);
        y_X_rom_int = interp1(t_out, y_X_rom', t_ref, 'linear', 'extrap')';
        results.MC.gre(i,j) = (norm(y_X_truth(:) - y_X_rom_int(:)) / norm(y_X_truth(:))) * 100;
        
        % 3. Rubin
        q0_r = Pc_rubin \ q0; F_r_handle = @(t) Pc_rubin' * F_fom_handle(t);
        solverRubin = TransientSolverOde_NEW(Mr_rubin, Kr_rubin, Cr_rubin); tic;
        [t_out, q_rom] = solverRubin.solve(tmax, dt, q0_r, zeros(size(q0_r)), F_r_handle, ...
            'ContactTargetDOF', rom_contact_dofs, 'ContactGap', gap_wall, 'ContactPenalty', k_contact, ...
            'ModelType', 'Rubin', 'RelTol', RT, 'AbsTol', AT);
        results.Rubin.cpu(i,j) = toc;
        y_X_rom = DummyStruct.AssemblyObj.unconstrain_vector(Pc_rubin * q_rom);
        y_X_rom = y_X_rom(contact_nodes_global_X, :);
        y_X_rom_int = interp1(t_out, y_X_rom', t_ref, 'linear', 'extrap')';
        results.Rubin.gre(i,j) = (norm(y_X_truth(:) - y_X_rom_int(:)) / norm(y_X_truth(:))) * 100;
        
        % 4. MCB
        q0_mcb = Pc_mcb \ q0; F_mcb_handle = @(t) Pc_mcb' * F_fom_handle(t);
        solverMCB = TransientSolverOde_NEW(Mr_mcb, Kr_mcb, Cr_mcb); tic;
        [t_out, q_rom] = solverMCB.solve(tmax, dt, q0_mcb, zeros(size(q0_mcb)), F_mcb_handle, ...
            'ContactTargetDOF', rom_contact_dofs, 'ContactGap', gap_wall, 'ContactPenalty', k_contact, ...
            'ModelType', 'MCB', 'RelTol', RT, 'AbsTol', AT);
        results.MCB.cpu(i,j) = toc;
        y_X_rom = DummyStruct.AssemblyObj.unconstrain_vector(Pc_mcb * q_rom);
        y_X_rom = y_X_rom(contact_nodes_global_X, :);
        y_X_rom_int = interp1(t_out, y_X_rom', t_ref, 'linear', 'extrap')';
        results.MCB.gre(i,j) = (norm(y_X_truth(:) - y_X_rom_int(:)) / norm(y_X_truth(:))) * 100;
        
        % 5. MN
        q0_mn = Pc_mn \ q0; F_mn_handle = @(t) Pc_mn' * F_fom_handle(t);
        solverMN = TransientSolverOde_NEW(Mr_mn, Kr_mn, Cr_mn); tic;
        [t_out, q_rom] = solverMN.solve(tmax, dt, q0_mn, zeros(size(q0_mn)), F_mn_handle, ...
            'ContactTargetDOF', rom_contact_dofs, 'ContactGap', gap_wall, 'ContactPenalty', k_contact, ...
            'ModelType', 'MN', 'RelTol', RT, 'AbsTol', AT);
        results.MN.cpu(i,j) = toc;
        y_X_rom = DummyStruct.AssemblyObj.unconstrain_vector(Pc_mn * q_rom);
        y_X_rom = y_X_rom(contact_nodes_global_X, :);
        y_X_rom_int = interp1(t_out, y_X_rom', t_ref, 'linear', 'extrap')';
        results.MN.gre(i,j) = (norm(y_X_truth(:) - y_X_rom_int(:)) / norm(y_X_truth(:))) * 100;
        
    end
end

%% --- 6. PLOT FRONTIERA DI PARETO ---
figure('Name', 'Pareto Front: Accuracy vs CPU Time', 'Position', [100 100 1200 800]);
colors = lines(5);
markers = {'o', 's', 'd', '^', 'v'};

hold on; grid on;
for m = 1:length(models)
    cpu_vec = results.(models{m}).cpu(:);
    gre_vec = results.(models{m}).gre(:);
    
    % Plot punti base
    scatter(cpu_vec, gre_vec, 150, markers{m}, 'MarkerFaceColor', colors(m,:), ...
        'MarkerEdgeColor', 'k', 'DisplayName', models{m});
    
    % Aggiungi etichette ai punti (RT e AT)
    count = 1;
    for i = 1:length(array_RelTol)
        for j = 1:length(array_AbsTol)
            label = sprintf(' R%g, A%g', log10(array_RelTol(i)), log10(array_AbsTol(j)));
            text(results.(models{m}).cpu(i,j), results.(models{m}).gre(i,j), label, ...
                'FontSize', 8, 'Color', [0.4 0.4 0.4]);
        end
    end
end

xlabel('Online CPU Time [s]', 'FontWeight', 'bold');
ylabel('Global Relative Error (GRE) [%]', 'FontWeight', 'bold');
title(sprintf('Sensitivity Analysis (Phi = %d, Q = %d, K_{mult} = %g)', fixed_phi, fixed_Q, fixed_k_mult));
legend('Location', 'northeast');
set(gca, 'XScale', 'log'); % Scala logaritmica aiuta a visualizzare bene i tempi