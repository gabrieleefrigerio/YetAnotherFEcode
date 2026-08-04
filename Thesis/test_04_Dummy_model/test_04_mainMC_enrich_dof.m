%% =====================================================================
%  MAIN BENCHMARKING: SELECTIVE MILMAN-CHU (Grid Search: K, Phi, Skip)
% =====================================================================
clear; close all; clc;

%% --- 2. INPUTS & SELECTORS ---
run_FOM    = 1; 
run_ROM_MC = 1;

% --- PARAMETRI DEL GRID SEARCH ---
array_linModes = [10, 20, 50, 70, 100]; % Number of linear modes for ROMs
array_k_mult    = [0.1, 1, 10];            % Moltiplicatori della rigidezza di contatto
array_node_skip = [1, 3, 5, 10];      % 1 = tutti i nodi, 3 = 1 ogni 3, ecc.

fixed_Q = 1000;   

shock_g = 1e5; shock_amp = shock_g * 9.81;    
t_shock = 10e-7; dt = 0.5e-8; tmax = (1e-3)/2; gap_wall = 1.5e-6;

%% --- 1. SETUP DIRECTORY ---
test_name = sprintf('SelectiveMC_GridSearch_%dg', shock_g);
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM');
save_dir = fullfile('results', sprintf('%s_%s', test_name, timestamp));
if ~exist(save_dir, 'dir'), mkdir(save_dir); end

%% --- 3. INIZIALIZZAZIONE STRUTTURA BASE ---
fprintf('Setting up the model...\n');
DummyStruct = AbaqusStructure();
DummyStruct.filename = 'DummyStructureAbaqus.inp'; DummyStruct.elementType = 'TRI3'; DummyStruct.build();

max_phi = max(array_linModes);
fprintf('Estrazione modale massima (%d modi)...\n', max_phi);
DummyStruct.compute_eigenmodes(max_phi);
DummyStruct.compute_rayleigh_damping(fixed_Q, fixed_Q);

Mc = DummyStruct.AssemblyObj.constrain_matrix(DummyStruct.M);
Kc = DummyStruct.AssemblyObj.constrain_matrix(DummyStruct.K);
Cc = DummyStruct.AssemblyObj.constrain_matrix(DummyStruct.C);
n_dofs_fom = size(Mc, 1);
k_base = max(diag(Kc));

contact_dofs = DummyStruct.get_contact_dofs(1); 
contact_nodes_global_X = (DummyStruct.contact_nodes - 1) * DummyStruct.MeshObj.nDOFPerNode + 1;

n_dof_per_node = DummyStruct.MeshObj.nDOFPerNode;
dir_vector = zeros(n_dofs_fom, 1); dir_vector(1:n_dof_per_node:n_dofs_fom) = 1; 
F_spatial_fom = Mc * dir_vector;
F_fom_handle = @(t) F_spatial_fom * shock_amp * sin(pi * t / t_shock) * (t <= t_shock);

q0 = zeros(n_dofs_fom, 1); qd0 = zeros(n_dofs_fom, 1);

%% --- 4. CICLO SUI PARAMETRI ---
for k_mult = array_k_mult
    k_contact = k_base * k_mult;
    fprintf('\n=======================================================\n');
    fprintf('   INIZIO CICLO CON K_MULT = %g\n', k_mult);
    fprintf('=======================================================\n');
    
    % ---------------------------------------------------------
    % ESECUZIONE FOM (Ground Truth per questo K)
    % ---------------------------------------------------------
    if run_FOM
        fprintf('\n[FOM] Esecuzione modello completo per K_mult = %g...\n', k_mult);
        tic; 
        solverFOM = TransientSolverOde(Mc, Kc, Cc);
        [t_fom, q_fom] = solverFOM.solve(tmax, dt, q0, qd0, F_fom_handle, ...
            'ContactTargetDOF', contact_dofs, 'ContactGap', gap_wall, ...
            'ContactPenalty', k_contact, 'ModelType', 'FOM');
        cpu_time = toc; 
        
        y_fom_full = DummyStruct.AssemblyObj.unconstrain_vector(q_fom);
        y_contact_nodes_X_fom = y_fom_full(contact_nodes_global_X, :);
        
        file_name_fom = sprintf('FOM_Ref_K%g.mat', k_mult);
        save(fullfile(save_dir, file_name_fom), 't_fom', 'y_contact_nodes_X_fom', 'cpu_time', 'k_mult');
    end

    % ---------------------------------------------------------
    % ESECUZIONE ROM SELECTIVE MILMAN-CHU
    % ---------------------------------------------------------
    if run_ROM_MC
        for phi = array_linModes
            for skip = array_node_skip
                fprintf('\n--- ROM MC | K: %g | Phi: %d | Skip: %d ---\n', k_mult, phi, skip);
                
                % SELEZIONE SOTTOINSIEME DI NODI DI INTERFACCIA
                idx_to_keep = 1:skip:length(contact_dofs);
                if idx_to_keep(end) ~= length(contact_dofs)
                    idx_to_keep = [idx_to_keep, length(contact_dofs)]; % Forza l'inclusione dell'ultimo nodo
                end
                enrich_dofs_selected = contact_dofs(idx_to_keep);
                num_enrich = length(enrich_dofs_selected);
                
                % 1. Costruzione Base Offline
                tic_offline = tic; 
                rom_mc = RomMC_enrich_dofs(DummyStruct, phi, contact_dofs, k_contact, true, enrich_dofs_selected);
                rom_mc.build();
                [Mr_mc, Kr_mc, Cr_mc] = rom_mc.get_reduced_matrices();
                offline_time = toc(tic_offline);
                
                Pc_mc = zeros(size(Kc, 1), size(rom_mc.P, 2));
                for i = 1:size(Pc_mc, 2), Pc_mc(:, i) = DummyStruct.AssemblyObj.constrain_vector(rom_mc.P(:, i)); end
                
                q0_mc = Pc_mc \ q0; qd0_mc = zeros(size(Pc_mc, 2), 1);
                F_mc_handle = @(t) Pc_mc' * F_fom_handle(t);
                
                % 2. Integrazione Temporale Online
                tic; 
                solverMC = TransientSolverOde(Mr_mc, Kr_mc, Cr_mc);
                [t_rom, q_rom_mc] = solverMC.solve(tmax, dt, q0_mc, qd0_mc, F_mc_handle, ...
                    'ContactTargetDOF', contact_dofs, 'ContactGap', gap_wall, ...
                    'ContactPenalty', k_contact, 'ModelType', 'MC', 'ProjectionMatrix', Pc_mc);
                cpu_time = toc; 
                
                % 3. Espansione
                q_mc_c = Pc_mc * q_rom_mc; 
                y_rom_full = DummyStruct.AssemblyObj.unconstrain_vector(q_mc_c);
                y_contact_nodes_X_romMC = y_rom_full(contact_nodes_global_X, :);
                
                file_name_rom = sprintf('ROM_MC_K%g_Phi%03d_Skip%02d_NVec%02d.mat', k_mult, phi, skip, num_enrich);
                save(fullfile(save_dir, file_name_rom), 't_rom', 'y_contact_nodes_X_romMC', 'cpu_time', 'offline_time', 'num_enrich', 'k_mult', 'phi', 'skip');
            end
        end
    end
end
fprintf('\nGrid Search Completato!\n');