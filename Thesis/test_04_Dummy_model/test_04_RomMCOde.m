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

%% --- 1. FULL ORDER MODEL (FOM) SETUP ---
DummyStruct = AbaqusStructure();
DummyStruct.filename = 'DummyStructureAbaqus.inp'; % Sostituisci col nome esatto del tuo .inp
DummyStruct.elementType = 'TRI3';             % Allinea al tipo di mesh (es. TRI3, QUAD8)
DummyStruct.thickness = 10e-6;

DummyStruct.build();

num_modi_plot = 5;
DummyStruct.compute_eigenmodes(num_modi_plot);

Mc = DummyStruct.AssemblyObj.constrain_matrix(DummyStruct.M);
Kc = DummyStruct.AssemblyObj.constrain_matrix(DummyStruct.K);
Cc = DummyStruct.AssemblyObj.constrain_matrix(DummyStruct.C);
n_dofs_fom = size(Mc, 1);

% Contact Parameters
gap_wall = 1.5e-6; 
k_contact = max(diag(Kc)) * 10; 

% --- Estrazione GdL e Nodi tramite Node Set (NUOVA STRUTTURA) ---
contact_dofs = DummyStruct.get_contact_dofs(1); 

nodi_parete = DummyStruct.contact_nodes;
coordinate_Y_parete = DummyStruct.nodes(nodi_parete, 2); 
[~, idx_top] = max(coordinate_Y_parete);
nodo_spia = nodi_parete(idx_top); 

dof_spia_global = (nodo_spia - 1) * DummyStruct.MeshObj.nDOFPerNode + 1;
dofs_parete_global = (nodi_parete - 1) * DummyStruct.MeshObj.nDOFPerNode + 1;

%% --- 3. REDUCED ORDER MODEL (ROM) BUILD ---
num_linear_modes = 20; % Numero di modi da includere
include_MC = 1;         % Arricchimento statico

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

%% --- 2. TRANSIENT SETUP ---
% Simulation parameters
dt = 0.5e-8; 
tmax = (1e-3)/2; 

% Initial Conditions
q0 = zeros(n_dofs_fom, 1);
qd0 = zeros(n_dofs_fom, 1);

% Shock Parameters
shock_g = 1e5;                 % Acceleration [g]
shock_amp = shock_g * 9.81;    % Acceleration [m/s^2]
t_shock = 10e-7;               % Impulse duration [s]

% Building Shock
n_dof_per_node = DummyStruct.MeshObj.nDOFPerNode;
dir_vector = zeros(n_dofs_fom, 1);
dir_vector(1:n_dof_per_node:n_dofs_fom) = 1; 
F_spatial_fom = Mc * dir_vector;
F_fom_handle = @(t) F_spatial_fom * shock_amp * sin(pi * t / t_shock) * (t <= t_shock);

%% --- 4. ROM TRANSIENT (ODE15s) ---
fprintf('\n=== Risoluzione ROM (Milman-Chu) ===\n');

% Proiezione Condizioni Iniziali 
q0_mc = Pc \ q0; 
qd0_mc = zeros(size(Pc, 2), 1);

% Proiezione Forzante nel sottospazio MC
F_mc_handle = @(t) Pc' * F_fom_handle(t);

% 1. Istanziamo la classe passando le matrici RIDOTTE (Mr, Kr, Cr)
solverROM = TransientSolverOde(Mr, Kr, Cr);

% 2. Lanciamo il solve. Notare ModelType = 'MC' e l'aggiunta di 'ProjectionMatrix'
[t_rom_mc, q_rom_mc] = solverROM.solve(tmax, dt, q0_mc, qd0_mc, F_mc_handle, ...
    'ContactTargetDOF', contact_dofs, ...
    'ContactGap', gap_wall, ...
    'ContactPenalty', k_contact, ...
    'ModelType', 'MC', ...
    'ProjectionMatrix', Pc);

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
fprintf('\n=== SEZIONE 5: Post-Processing e Confronto ===\n');

% =====================================================================
% CARICAMENTO DATI FOM (NECESSARIO PER IL CONFRONTO)
% Poiché la memoria è stata pulita all'inizio, recuperiamo i risultati FOM
% =====================================================================
try
    % CHIEDI ALL'UTENTE IL FILE FOM SE NON E' NELLA DIRECTORY CORRENTE
    [fom_file, fom_path] = uigetfile('*.mat', 'Seleziona il file FOM_ode15s.mat generato in precedenza');
    if isequal(fom_file,0)
        disp('Caricamento FOM annullato. Il plot di confronto non includerà i dati FOM.');
        t_fom_ode = []; y_spia_fom_ode = [];
    else
        load(fullfile(fom_path, fom_file), 't_fom_ode', 'y_spia_fom_ode');
        fprintf('Dati FOM caricati con successo da: %s\n', fom_file);
    end
catch
    warning('Errore nel caricamento del file FOM.');
end

% =====================================================================
% CALCOLO DEL GLOBAL RELATIVE ERROR (GRE) SU NODO SPIA
% =====================================================================
if ~isempty(t_fom_ode) && ~isempty(y_spia_fom_ode)
    % Interpoliamo il ROM sulla griglia temporale del FOM
    y_spia_rom_interp = interp1(t_rom_mc, y_spia_rom_mc, t_fom_ode', 'linear', 'extrap');
    
    % Calcolo dell'errore globale in norma L2 relativa
    errore_assoluto = y_spia_fom_ode - y_spia_rom_interp;
    gre_spostamento = norm(errore_assoluto) / norm(y_spia_fom_ode);
    gre_percentage = gre_spostamento * 100;
    fprintf('Global Relative Error (Spostamento Nodo %d): %.4f%%\n', nodo_spia, gre_percentage);
else
    gre_percentage = NaN; % Se il FOM non è caricato, saltiamo il GRE
end

% =====================================================================
% --- A. Confronto Cinematica (Nodo Spia) con GRE ---
% =====================================================================
figure('Name', 'Confronto Cinematica FOM vs ROM', 'Color', 'w', 'Units', 'normalized', 'Position', [0.1 0.4 0.4 0.4]);

if ~isempty(t_fom_ode)
    plot(t_fom_ode, y_spia_fom_ode * 1e6, 'b-', 'LineWidth', 2, 'DisplayName', 'FOM (Completo)');
    hold on;
end

% Inseriamo il valore del GRE direttamente nella legenda del ROM Milman-Chu
if ~isnan(gre_percentage)
    rom_legend_str = sprintf('ROM (Milman-Chu) - GRE: %.3f%%', gre_percentage);
else
    rom_legend_str = 'ROM (Milman-Chu)';
end

plot(t_rom_mc, y_spia_rom_mc * 1e6, 'r--', 'LineWidth', 1.5, 'DisplayName', rom_legend_str);
yline(gap_wall * 1e6, 'k:', 'LineWidth', 2, 'DisplayName', 'Parete');

title(sprintf('Cinematica di Impatto - Nodo %d', nodo_spia));
xlabel('Tempo [s]'); ylabel('Spostamento X [\mum]');
legend('Location', 'best'); grid on;