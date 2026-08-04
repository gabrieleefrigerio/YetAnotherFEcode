%% =====================================================================
%  MAIN BENCHMARKING (Updated for V2 - 4 Contact Interfaces + MT)
% =====================================================================
clear; close all; clc;

%% --- 1. INPUTS & SELECTORS ---
% Select Models
run_FOM       = 0; 
run_ROM_MT    = 0; 
run_ROM_MC    = 1;
run_ROM_Rubin = 0;
run_ROM_MCB   = 0; 
run_ROM_MN    = 0;

% Vettori dei Parametri per il Loop (dal file originale)
array_linModes = [10, 50, 100, 150, 200]; 
array_QFactor   = [1000];               
array_k_mult    = [10];              

% Transient Parameters
shock_g = 1e5;
shock_amp = shock_g * 9.81;    
t_shock = 10e-7;               
dt = 0.5e-8;
tmax = (1e-3)/2;

% --- Gaps FIRMATI specifici per ciascuna interfaccia ---
% Il segno indica in quale direzione si trova il muro rispetto all'origine
gap_T =  5e-6;   % Muro in Y positivo
gap_B = -1.5e-6; % Muro in Y negativo
gap_L = -5e-6;   % Muro in X negativo
gap_R =  1.5e-6; % Muro in X positivo

%% --- 2. SETUP DIRECTORY ---
test_name = sprintf('Shock_%dg_%fs_%dxspring_Ode15_V2', shock_g, tmax);
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM');
dir_name = sprintf('%s_%s', test_name, timestamp);
save_dir = fullfile('results_V2', dir_name);
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end
fprintf('Directory creata: %s\n\n', save_dir);

%% --- 3. INIZIALIZZAZIONE STRUTTURA BASE ---
fprintf('Setting up the model...\n');
DummyStruct = AbaqusStructure_V2();                  % <--- Usiamo la classe V2
DummyStruct.filename = 'DummyStructureAbaqus_V4.inp'; % <--- Mesh V2
DummyStruct.elementType = 'TRI3';           
DummyStruct.build();

max_phi = max(array_linModes);
fprintf('Estrazione modale di %d modi...\n', max_phi);
DummyStruct.compute_eigenmodes(max_phi);

% Extracting Structure's data
Mc = DummyStruct.AssemblyObj.constrain_matrix(DummyStruct.M);
Kc = DummyStruct.AssemblyObj.constrain_matrix(DummyStruct.K);
n_dofs_fom = size(Mc, 1);
k_base = max(diag(Kc)); 

%% --- 3.1 SETUP INTERFACCE DI CONTATTO (T, B, L, R) ---
% Estrazione GdL per le 4 facce (2=Y, 1=X)
dofs_T = DummyStruct.get_specific_contact_dofs('T', 2);
dofs_B = DummyStruct.get_specific_contact_dofs('B', 2);
dofs_L = DummyStruct.get_specific_contact_dofs('L', 1);
dofs_R = DummyStruct.get_specific_contact_dofs('R', 1);

% Vettori concatenati per il solutore ODE
contact_dofs = [dofs_T; dofs_B; dofs_L; dofs_R];
gaps_array   = [repmat(gap_T, length(dofs_T), 1); ...
                repmat(gap_B, length(dofs_B), 1); ...
                repmat(gap_L, length(dofs_L), 1); ...
                repmat(gap_R, length(dofs_R), 1)];

% Struct organizzata per il Post-Processing
Interfaces = struct();
labels = {'T', 'B', 'L', 'R'};
nodes_list = {DummyStruct.contact_nodes_T, DummyStruct.contact_nodes_B, ...
              DummyStruct.contact_nodes_L, DummyStruct.contact_nodes_R};

for i = 1:length(labels)
    loc = labels{i};
    n = nodes_list{i};
    Interfaces.(loc).nodes = n;
    if ~isempty(n)
        Interfaces.(loc).global_X = (n - 1) * DummyStruct.MeshObj.nDOFPerNode + 1;
        Interfaces.(loc).global_Y = (n - 1) * DummyStruct.MeshObj.nDOFPerNode + 2;
        Interfaces.(loc).coord_X  = DummyStruct.nodes(n, 1);
        Interfaces.(loc).coord_Y  = DummyStruct.nodes(n, 2);
    end
end

%% --- 3.2 Forcing & Initial Conditions ---
n_dof_per_node = DummyStruct.MeshObj.nDOFPerNode;
dir_vector = zeros(n_dofs_fom, 1);
dir_vector(1:n_dof_per_node:n_dofs_fom) = 1; 
F_spatial_fom = Mc * dir_vector;
F_fom_handle = @(t) F_spatial_fom * shock_amp * sin(pi * t / t_shock) * (t <= t_shock);

q0 = zeros(n_dofs_fom, 1);
qd0 = zeros(n_dofs_fom, 1);

%% --- 4. Run FOM ---
if run_FOM
    fprintf('\n=========================================\n');
    fprintf('           Run FOM                \n');
    fprintf('=========================================\n');
    for Q = array_QFactor
        fprintf('\n[FOM] FOM with Q = %d...\n', Q);
        DummyStruct.compute_rayleigh_damping(Q, Q);
        Cc = DummyStruct.AssemblyObj.constrain_matrix(DummyStruct.C);
        
        for k_mult = array_k_mult
            k_contact = k_base * k_mult;
            fprintf('      -> Testing Contact Penalty Mult = %g\n', k_mult);
            
            tic; 
            solverFOM = TransientSolverOde_V2(Mc, Kc, Cc);
            [t_fom, q_fom] = solverFOM.solve(tmax, dt, q0, qd0, F_fom_handle, ...
                'ContactTargetDOF', contact_dofs, ...
                'ContactGap', gaps_array, ...
                'ContactPenalty', k_contact, ...
                'ModelType', 'FOM');
            cpu_time = toc; 
                
            % Estrazione Risultati Dinamica (suddivisi per interfaccia)
            y_fom_full = DummyStruct.AssemblyObj.unconstrain_vector(q_fom);
            y_contact = struct();
            for loc = labels
                l = loc{1};
                if ~isempty(Interfaces.(l).nodes)
                    y_contact.(l).X = y_fom_full(Interfaces.(l).global_X, :);
                    y_contact.(l).Y = y_fom_full(Interfaces.(l).global_Y, :);
                end
            end
            
            file_name = sprintf('FOM_Q%04d_K%g.mat', Q, k_mult);
            save(fullfile(save_dir, file_name), 't_fom', 'y_contact', 'Interfaces', 'cpu_time');
        end
    end
end

%% --- 5. Run ROM ---
fprintf('\n=========================================\n');
fprintf('           Run ROM                \n');
fprintf('=========================================\n');
for Q = array_QFactor
    DummyStruct.compute_rayleigh_damping(Q, Q);
    
    for k_mult = array_k_mult
        k_contact = k_base * k_mult;
        
        for phi = array_linModes
            fprintf('\n--- ROM Sweep | Phi: %d | Q: %d | K_mult: %g ---\n', phi, Q, k_mult);
            
            % ==============================================
            % ROM 0: MODAL TRUNCATION (MT) 
            % (RomMC with include_MC = 0)
            % ==============================================
            if run_ROM_MT
                tic_offline = tic; 
                fprintf(' Building ROM Modal Truncation (MT) Basis...\n');
                include_MC = 0; % FLAG SPENTO
                rom_mt = RomMC(DummyStruct, phi, contact_dofs, k_contact, include_MC);
                rom_mt.build();
                [Mr_mt, Kr_mt, Cr_mt] = rom_mt.get_reduced_matrices();
                offline_time = toc(tic_offline); 
                
                Pc_mt = zeros(size(Kc, 1), size(rom_mt.P, 2));
                for i = 1:size(Pc_mt, 2), Pc_mt(:, i) = DummyStruct.AssemblyObj.constrain_vector(rom_mt.P(:, i)); end
                q0_mt = Pc_mt \ q0; qd0_mt = zeros(size(Pc_mt, 2), 1);
                F_mt_handle = @(t) Pc_mt' * F_fom_handle(t);
                
                % Utilizziamo 'ModelType', 'MC' perché il solutore deve comunque usare 
                % la ProjectionMatrix per recuperare i GdL fisici per il contatto.
                tic; solverMT = TransientSolverOde_V2(Mr_mt, Kr_mt, Cr_mt);
                [t_rom, q_rom_mt] = solverMT.solve(tmax, dt, q0_mt, qd0_mt, F_mt_handle, ...
                    'ContactTargetDOF', contact_dofs, 'ContactGap', gaps_array, ...
                    'ContactPenalty', k_contact, 'ModelType', 'MC', 'ProjectionMatrix', Pc_mt);
                cpu_time = toc; 
                    
                q_mt_c = Pc_mt * q_rom_mt; 
                y_rom_full = DummyStruct.AssemblyObj.unconstrain_vector(q_mt_c);
                
                y_contact = struct();
                for loc = labels, l=loc{1}; if ~isempty(Interfaces.(l).nodes)
                    y_contact.(l).X = y_rom_full(Interfaces.(l).global_X, :);
                    y_contact.(l).Y = y_rom_full(Interfaces.(l).global_Y, :);
                end; end
                
                file_name = sprintf('ROM_MT_Phi%03d_Q%04d_K%g.mat', phi, Q, k_mult);
                save(fullfile(save_dir, file_name), 't_rom', 'y_contact', 'Interfaces', 'cpu_time', 'offline_time');
            end

            % ==============================================
            % ROM 1: MILMAN-CHU (MC)
            % ==============================================
            if run_ROM_MC
                tic_offline = tic; 
                fprintf(' Building ROM MC Basis...\n');
                include_MC = 1; % FLAG ACCESO
                rom_mc = RomMC(DummyStruct, phi, contact_dofs, k_contact, include_MC);
                rom_mc.build();
                [Mr_mc, Kr_mc, Cr_mc] = rom_mc.get_reduced_matrices();
                offline_time = toc(tic_offline); 
                
                Pc_mc = zeros(size(Kc, 1), size(rom_mc.P, 2));
                for i = 1:size(Pc_mc, 2), Pc_mc(:, i) = DummyStruct.AssemblyObj.constrain_vector(rom_mc.P(:, i)); end
                q0_mc = Pc_mc \ q0; qd0_mc = zeros(size(Pc_mc, 2), 1);
                F_mc_handle = @(t) Pc_mc' * F_fom_handle(t);
                
                tic; solverMC = TransientSolverOde_V2(Mr_mc, Kr_mc, Cr_mc);
                [t_rom, q_rom_mc] = solverMC.solve(tmax, dt, q0_mc, qd0_mc, F_mc_handle, ...
                    'ContactTargetDOF', contact_dofs, 'ContactGap', gaps_array, ...
                    'ContactPenalty', k_contact, 'ModelType', 'MC', 'ProjectionMatrix', Pc_mc);
                cpu_time = toc; 
                    
                q_mc_c = Pc_mc * q_rom_mc; 
                y_rom_full = DummyStruct.AssemblyObj.unconstrain_vector(q_mc_c);
                
                y_contact = struct();
                for loc = labels, l=loc{1}; if ~isempty(Interfaces.(l).nodes)
                    y_contact.(l).X = y_rom_full(Interfaces.(l).global_X, :);
                    y_contact.(l).Y = y_rom_full(Interfaces.(l).global_Y, :);
                end; end
                
                file_name = sprintf('ROM_MC_Phi%03d_Q%04d_K%g.mat', phi, Q, k_mult);
                save(fullfile(save_dir, file_name), 't_rom', 'y_contact', 'Interfaces', 'cpu_time', 'offline_time');
            end
            
            % ==============================================
            % ROM 2: RUBIN / ROM 3: MCB / ROM 4: MN
            % ==============================================
            rom_types = {'Rubin', 'MCB', 'MN'};
            run_flags = [run_ROM_Rubin, run_ROM_MCB, run_ROM_MN];
            
            for r_idx = 1:3
                if run_flags(r_idx)
                    tic_offline = tic; 
                    fprintf(' Building ROM %s Basis...\n', rom_types{r_idx});
                    
                    switch rom_types{r_idx}
                        case 'Rubin', rom_obj = RomRubin(DummyStruct, phi, contact_dofs);
                        case 'MCB',   rom_obj = RomMCB(DummyStruct, phi, contact_dofs);
                        case 'MN',    rom_obj = RomMN(DummyStruct, phi, contact_dofs);
                    end
                    rom_obj.build();
                    [Mr_r, Kr_r, Cr_r] = rom_obj.get_reduced_matrices();
                    Pc_r = rom_obj.Pc;
                    offline_time = toc(tic_offline); 
                   
                    rom_contact_dofs = 1:length(contact_dofs);
                    
                    q0_r = Pc_r \ q0; qd0_r = zeros(size(Pc_r, 2), 1);
                    F_r_handle = @(t) Pc_r' * F_fom_handle(t);
                    
                    tic; solverROM = TransientSolverOde_V2(Mr_r, Kr_r, Cr_r);
                    [t_rom, q_rom] = solverROM.solve(tmax, dt, q0_r, qd0_r, F_r_handle, ...
                        'ContactTargetDOF', rom_contact_dofs, 'ContactGap', gaps_array, ...
                        'ContactPenalty', k_contact, 'ModelType', rom_types{r_idx}); 
                    cpu_time = toc; 
                        
                    q_c = Pc_r * q_rom; 
                    y_rom_full = DummyStruct.AssemblyObj.unconstrain_vector(q_c);
                    
                    y_contact = struct();
                    for loc = labels, l=loc{1}; if ~isempty(Interfaces.(l).nodes)
                        y_contact.(l).X = y_rom_full(Interfaces.(l).global_X, :);
                        y_contact.(l).Y = y_rom_full(Interfaces.(l).global_Y, :);
                    end; end
                    
                    file_name = sprintf('ROM_%s_Phi%03d_Q%04d_K%g.mat', rom_types{r_idx}, phi, Q, k_mult);
                    save(fullfile(save_dir, file_name), 't_rom', 'y_contact', 'Interfaces', 'cpu_time', 'offline_time'); 
                end
            end
        end
    end
end