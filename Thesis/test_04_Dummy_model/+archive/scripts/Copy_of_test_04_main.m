%% =====================================================================
%  MAIN BENCHMARKING
% =====================================================================
clear; close all; clc;
%% --- 2. INPUTS & SELECTORS ---
% Select Models
run_FOM       = 0; 
run_ROM_MT    = 0; 
run_ROM_MC    = 0;
run_ROM_Rubin = 0;
run_ROM_MCB   = 1; 
run_ROM_MN    = 1;
% Vettori dei Parametri per il Loop
array_linModes = [10, 20, 50, 70, 100]; % Number of linear modes for ROMs
array_QFactor   = [1000]; % Q factors
array_k_mult    = [10]; % Moltiplicatori per la rigidezza di contatto
% Transient Parameters
shock_g = 1e5;
shock_amp = shock_g * 9.81;    
t_shock = 10e-7;               
dt = 0.5e-8;
tmax = (1e-3)/2;
gap_wall = 1.5e-6;
%% --- 1. SETUP DIRECTORY ---
% Set up custom name for Results Folder
test_name = sprintf('Shock_%dg_%fs_%dxspring_Ode15', shock_g, tmax);
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM');
dir_name = sprintf('%s_%s', test_name, timestamp);
save_dir = fullfile('results', dir_name);
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end
fprintf('Directory creata: %s\n\n', save_dir);
%% --- 3. INIZIALIZZAZIONE STRUTTURA BASE ---
fprintf('Setting up the model...\n');
DummyStruct = AbaqusStructure();
DummyStruct.filename = 'DummyStructureAbaqus.inp'; 
DummyStruct.elementType = 'TRI3';           
DummyStruct.build();
max_phi = max(array_linModes);
fprintf('Estrazione modale di %d modi...\n', max_phi);
DummyStruct.compute_eigenmodes(max_phi);
% Extracting Structure's data
Mc = DummyStruct.AssemblyObj.constrain_matrix(DummyStruct.M);
Kc = DummyStruct.AssemblyObj.constrain_matrix(DummyStruct.K);
n_dofs_fom = size(Mc, 1);
k_base = max(diag(Kc)); % Rigidezza base da moltiplicare nel ciclo
% Extracting contact DoFs
contact_dofs = DummyStruct.get_contact_dofs(1); 
contact_nodes = DummyStruct.contact_nodes;
coord_Y_contact_nodes = DummyStruct.nodes(contact_nodes, 2); 
% GdL Global for X and Y
contact_nodes_global_X = (contact_nodes - 1) * DummyStruct.MeshObj.nDOFPerNode + 1;
contact_nodes_global_Y = (contact_nodes - 1) * DummyStruct.MeshObj.nDOFPerNode + 2;
% Forcing
n_dof_per_node = DummyStruct.MeshObj.nDOFPerNode;
dir_vector = zeros(n_dofs_fom, 1);
dir_vector(1:n_dof_per_node:n_dofs_fom) = 1; 
F_spatial_fom = Mc * dir_vector;
F_fom_handle = @(t) F_spatial_fom * shock_amp * sin(pi * t / t_shock) * (t <= t_shock);
% Initial Conditions
q0 = zeros(n_dofs_fom, 1);
qd0 = zeros(n_dofs_fom, 1);
%% --- 4. Run FOM ---
if run_FOM
    fprintf('\n=========================================\n');
    fprintf('           Run FOM                \n');
    fprintf('=========================================\n');
    for Q = array_QFactor
        fprintf('\n[FOM] FOM with Q = %d...\n', Q);
        
        % 1. Update C matrix based on current Q factor
        DummyStruct.compute_rayleigh_damping(Q, Q);
        Cc = DummyStruct.AssemblyObj.constrain_matrix(DummyStruct.C);
        
        for k_mult = array_k_mult
            k_contact = k_base * k_mult;
            fprintf('      -> Testing Contact Penalty Mult = %g\n', k_mult);
            
            % 2. Run FOM solver
            tic; 
            solverFOM = TransientSolverOde(Mc, Kc, Cc);
            [t_fom, q_fom] = solverFOM.solve(tmax, dt, q0, qd0, F_fom_handle, ...
                'ContactTargetDOF', contact_dofs, ...
                'ContactGap', gap_wall, ...
                'ContactPenalty', k_contact, ...
                'ModelType', 'FOM');
            cpu_time = toc; 
                
            % 3. Extract FOM results
            y_fom_full = DummyStruct.AssemblyObj.unconstrain_vector(q_fom);
            y_contact_nodes_X_fom = y_fom_full(contact_nodes_global_X, :);
            y_contact_nodes_Y_fom = y_fom_full(contact_nodes_global_Y, :);
            
            % Save
            file_name = sprintf('FOM_Q%04d_K%g.mat', Q, k_mult);
            save(fullfile(save_dir, file_name), 't_fom', 'y_contact_nodes_X_fom', 'y_contact_nodes_Y_fom', 'coord_Y_contact_nodes', 'cpu_time');
        end
    end
end
%% --- 5. Run ROM ---
fprintf('\n=========================================\n');
fprintf('           Run ROM                \n');
fprintf('=========================================\n');
for Q = array_QFactor
    % 1. Update C matrix based on current Q factor
    DummyStruct.compute_rayleigh_damping(Q, Q);
    
    for k_mult = array_k_mult
        k_contact = k_base * k_mult;
        
        for phi = array_linModes
            fprintf('\n--- ROM Sweep | Phi: %d | Q: %d | K_mult: %g ---\n', phi, Q, k_mult);
            
            % ==============================================
            % ROM 0: MODAL TRUNCATION (MT)
            % ==============================================
            if run_ROM_MT
                tic_offline = tic; % <--- INIZIO CRONOMETRO OFFLINE
                fprintf(' Building ROM Modal Truncation (MT) Basis...\n');
                include_MC = 0; % <--- FLAG SPENTO
                rom_mt = RomMC(DummyStruct, phi, contact_dofs, k_contact, include_MC);
                rom_mt.build();
                [Mr_mt, Kr_mt, Cr_mt] = rom_mt.get_reduced_matrices();
                offline_time = toc(tic_offline); % <--- FINE CRONOMETRO OFFLINE
                
                % IC and Forcing Projection
                Pc_mt = zeros(size(Kc, 1), size(rom_mt.P, 2));
                for i = 1:size(Pc_mt, 2)
                    Pc_mt(:, i) = DummyStruct.AssemblyObj.constrain_vector(rom_mt.P(:, i));
                end
                q0_mt = Pc_mt \ q0; 
                qd0_mt = zeros(size(Pc_mt, 2), 1);
                F_mt_handle = @(t) Pc_mt' * F_fom_handle(t);
                
                % Solve MT
                tic; 
                solverMT = TransientSolverOde(Mr_mt, Kr_mt, Cr_mt);
                [t_rom, q_rom_mt] = solverMT.solve(tmax, dt, q0_mt, qd0_mt, F_mt_handle, ...
                    'ContactTargetDOF', contact_dofs, ...
                    'ContactGap', gap_wall, ...
                    'ContactPenalty', k_contact, ...
                    'ModelType', 'MC', ...  % <--- Come per MC, serve la ProjectionMatrix
                    'ProjectionMatrix', Pc_mt);
                cpu_time = toc; 
                    
                % Expansion and Extraction
                q_mt_c = Pc_mt * q_rom_mt; 
                y_rom_full = DummyStruct.AssemblyObj.unconstrain_vector(q_mt_c);
                y_contact_nodes_X_romMT = y_rom_full(contact_nodes_global_X, :);
                y_contact_nodes_Y_romMT = y_rom_full(contact_nodes_global_Y, :);
                
                file_name = sprintf('ROM_MT_Phi%03d_Q%04d_K%g.mat', phi, Q, k_mult);
                save(fullfile(save_dir, file_name), 't_rom', 'y_contact_nodes_X_romMT', 'y_contact_nodes_Y_romMT', 'coord_Y_contact_nodes', 'cpu_time', 'offline_time');
            end

            % ==============================================
            % ROM 1: MILMAN-CHU (MC)
            % ==============================================
            if run_ROM_MC
                tic_offline = tic; % <--- INIZIO CRONOMETRO OFFLINE
                fprintf(' Building ROM MC Basis...\n');
                include_MC = 1;
                rom_mc = RomMC(DummyStruct, phi, contact_dofs, k_contact, include_MC);
                rom_mc.build();
                [Mr_mc, Kr_mc, Cr_mc] = rom_mc.get_reduced_matrices();
                offline_time = toc(tic_offline); % <--- FINE CRONOMETRO OFFLINE
                
                % IC and Forcing Projection
                Pc_mc = zeros(size(Kc, 1), size(rom_mc.P, 2));
                for i = 1:size(Pc_mc, 2)
                    Pc_mc(:, i) = DummyStruct.AssemblyObj.constrain_vector(rom_mc.P(:, i));
                end
                q0_mc = Pc_mc \ q0; 
                qd0_mc = zeros(size(Pc_mc, 2), 1);
                F_mc_handle = @(t) Pc_mc' * F_fom_handle(t);
                
                % Solve MC
                tic; 
                solverMC = TransientSolverOde(Mr_mc, Kr_mc, Cr_mc);
                [t_rom, q_rom_mc] = solverMC.solve(tmax, dt, q0_mc, qd0_mc, F_mc_handle, ...
                    'ContactTargetDOF', contact_dofs, ...
                    'ContactGap', gap_wall, ...
                    'ContactPenalty', k_contact, ...
                    'ModelType', 'MC', ...
                    'ProjectionMatrix', Pc_mc);
                cpu_time = toc; 
                    
                % Expansion and Extraction
                q_mc_c = Pc_mc * q_rom_mc; 
                y_rom_full = DummyStruct.AssemblyObj.unconstrain_vector(q_mc_c);
                y_contact_nodes_X_romMC = y_rom_full(contact_nodes_global_X, :);
                y_contact_nodes_Y_romMC = y_rom_full(contact_nodes_global_Y, :);
                
                file_name = sprintf('ROM_MC_Phi%03d_Q%04d_K%g.mat', phi, Q, k_mult);
                save(fullfile(save_dir, file_name), 't_rom', 'y_contact_nodes_X_romMC', 'y_contact_nodes_Y_romMC', 'coord_Y_contact_nodes', 'cpu_time', 'offline_time');
            end
            
            % ==============================================
            % ROM 2: RUBIN (CMS)
            % ==============================================
            if run_ROM_Rubin
                tic_offline = tic; % <--- INIZIO CRONOMETRO OFFLINE
                fprintf(' Building ROM Rubin Basis...\n');
                rom_rubin = RomRubin(DummyStruct, phi, contact_dofs);
                rom_rubin.build();
                [Mr_rubin, Kr_rubin, Cr_rubin] = rom_rubin.get_reduced_matrices();
                Pc_rubin = rom_rubin.Pc;
                offline_time = toc(tic_offline); % <--- FINE CRONOMETRO OFFLINE
               
                % Interface DoFs are on top (1:n_bnd)
                rom_contact_dofs = 1:length(contact_dofs);
                
                % IC and Forcing Projection
                q0_rubin = Pc_rubin \ q0; 
                qd0_rubin = zeros(size(Pc_rubin, 2), 1);
                F_rubin_handle = @(t) Pc_rubin' * F_fom_handle(t);
                
                % Solve Rubin
                tic; 
                solverRubin = TransientSolverOde(Mr_rubin, Kr_rubin, Cr_rubin);
                [t_rom, q_rom_rubin] = solverRubin.solve(tmax, dt, q0_rubin, qd0_rubin, F_rubin_handle, ...
                    'ContactTargetDOF', rom_contact_dofs, ...
                    'ContactGap', gap_wall, ...
                    'ContactPenalty', k_contact, ...
                    'ModelType', 'Rubin'); 
                cpu_time = toc; 
                    
                % Expansion and Extraction
                q_rubin_c = Pc_rubin * q_rom_rubin; 
                y_rom_full = DummyStruct.AssemblyObj.unconstrain_vector(q_rubin_c);
                y_contact_nodes_X_romRubin = y_rom_full(contact_nodes_global_X, :);
                y_contact_nodes_Y_romRubin = y_rom_full(contact_nodes_global_Y, :);
                
                file_name = sprintf('ROM_Rubin_Phi%03d_Q%04d_K%g.mat', phi, Q, k_mult);
                save(fullfile(save_dir, file_name), 't_rom', 'y_contact_nodes_X_romRubin', 'y_contact_nodes_Y_romRubin', 'coord_Y_contact_nodes', 'cpu_time', 'offline_time'); 
            end
            
            % ==============================================
            % ROM 3: MASSLESS CRAIG-BAMPTON (MCB)
            % ==============================================
            if run_ROM_MCB
tic_offline = tic; % <--- INIZIO CRONOMETRO OFFLINE
                fprintf(' Building ROM Massless Craig-Bampton (MCB) Basis...\n');
                
                % Utilizziamo la classe che hai scritto tu
                rom_mcb = RomMCB(DummyStruct, phi, contact_dofs);
                rom_mcb.build();
                [Mr_mcb, Kr_mcb, Cr_mcb] = rom_mcb.get_reduced_matrices();
                Pc_mcb = rom_mcb.Pc;
                offline_time = toc(tic_offline); % <--- FINE CRONOMETRO OFFLINE
               
                % I GDL fisici di contatto (il bordo) sono ordinati in testa (1:n_bnd)
                n_bnd_mcb = length(contact_dofs);
                
                % Proiezione delle Condizioni Iniziali e delle Forze
                q0_mcb = Pc_mcb \ q0; 
                qd0_mcb = zeros(size(Pc_mcb, 2), 1);
                F_mcb_handle = @(t) Pc_mcb' * F_fom_handle(t);
                
                % Risoluzione con l'integratore dedicato Massless Verlet
                tic; 
                solverMCB = MasslessVerletSolver(Mr_mcb, Kr_mcb, Cr_mcb);
                [t_rom, q_rom_mcb] = solverMCB.solve(tmax, dt, q0_mcb, qd0_mcb, F_mcb_handle, ...
                    n_bnd_mcb, gap_wall); 
                cpu_time = toc; 
                    
                % Espansione allo spazio vincolato completo e recupero spostamenti fisici
                q_mcb_c = Pc_mcb * q_rom_mcb; 
                y_rom_full = DummyStruct.AssemblyObj.unconstrain_vector(q_mcb_c);
                y_contact_nodes_X_romMCB = y_rom_full(contact_nodes_global_X, :);
                y_contact_nodes_Y_romMCB = y_rom_full(contact_nodes_global_Y, :);
                
                % Salvataggio dei Risultati
                file_name = sprintf('ROM_MCB_Phi%03d_Q%04d_K%g.mat', phi, Q, k_mult);
                save(fullfile(save_dir, file_name), 't_rom', 'y_contact_nodes_X_romMCB', 'y_contact_nodes_Y_romMCB', 'coord_Y_contact_nodes', 'cpu_time', 'offline_time'); 
                fprintf('Done MCB solver sweep.\n');
            end
% ==============================================
            % ROM 4: MACNEAL (MN) - Free-Interface CMS
            % ==============================================
            if run_ROM_MN
                tic_offline = tic; 
                fprintf(' Building ROM MacNeal (MN) Basis...\n');
                
                rom_mn = RomMN(DummyStruct, phi, contact_dofs);
                rom_mn.build();
                [Mr_mn, Kr_mn, Cr_mn] = rom_mn.get_reduced_matrices();
                Pc_mn = rom_mn.Pc;
                offline_time = toc(tic_offline); 
               
                % I GDL fisici di contatto sono in testa (1:n_bnd)
                n_bnd_mn = length(contact_dofs);
                
                % Proiezione Condizioni Iniziali e Forze
                q0_mn = Pc_mn \ q0; 
                qd0_mn = zeros(size(Pc_mn, 2), 1);
                F_mn_handle = @(t) Pc_mn' * F_fom_handle(t);
                
                % Risoluzione con lo STESSO integratore Massless Verlet
                tic; 
                solverMN = MasslessVerletSolver(Mr_mn, Kr_mn, Cr_mn);
                [t_rom, q_rom_mn] = solverMN.solve(tmax, dt, q0_mn, qd0_mn, F_mn_handle, ...
                    n_bnd_mn, gap_wall); 
                cpu_time = toc; 
                    
                % Espansione
                q_mn_c = Pc_mn * q_rom_mn; 
                y_rom_full = DummyStruct.AssemblyObj.unconstrain_vector(q_mn_c);
                y_contact_nodes_X_romMN = y_rom_full(contact_nodes_global_X, :);
                y_contact_nodes_Y_romMN = y_rom_full(contact_nodes_global_Y, :);
                
                file_name = sprintf('ROM_MN_Phi%03d_Q%04d_K%g.mat', phi, Q, k_mult);
                save(fullfile(save_dir, file_name), 't_rom', 'y_contact_nodes_X_romMN', 'y_contact_nodes_Y_romMN', 'coord_Y_contact_nodes', 'cpu_time', 'offline_time'); 
                fprintf('Done MacNeal solver sweep.\n');
            end
        end
    end
end