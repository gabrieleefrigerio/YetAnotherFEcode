%% =====================================================================
%  SCRIPT FOR SOLVER COMPARISON
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
% Simulation parameters
dt = 0.5e-8; 
tmax = (1e-3)/2; 

% Initial Conditions (Spostamenti e Velocità nulle a t=0)
q0 = zeros(n_dofs_fom, 1);
qd0 = zeros(n_dofs_fom, 1);

% Shock Parameters
shock_g = 1e5;                 % Accelerazione di picco [g]
shock_amp = shock_g * 9.81;    % Accelerazione [m/s^2]
t_shock = 10e-7;               % Durata dell'impulso [s]

% Building Shock / Direzione di accelerazione
n_dof_per_node = DummyStruct.MeshObj.nDOFPerNode;
dir_vector = zeros(n_dofs_fom, 1);
dir_vector(1:n_dof_per_node:n_dofs_fom) = 1; % Forziamo 1 solo sulle componenti X

% --- Definizione del Vettore di Accelerazione Iniziale qdd0 ---
qdd0 = zeros(n_dofs_fom, 1); 

% Forzante 
F_spatial_fom = Mc * dir_vector;
F_fom_handle = @(t) F_spatial_fom * shock_amp * sin(pi * t / t_shock) * (t <= t_shock);

%% --- 3. FOM TRANSIENT (GENERALIZED-ALPHA / NEWMARK) ---
fprintf('\n=== Risoluzione FOM con Newmark/Gen-Alpha ===\n');
solverFOM = TransientSolver(Mc, Kc, Cc);

% Lancio con Generalized-alpha
tic;
[t_fom_newmark, q_fom_newmark] = solverFOM.solve(tmax, dt, q0, qd0, F_fom_handle, ...
    'ContactTargetDOF', contact_dofs, ...
    'ContactGap', gap_wall, ...
    'ContactPenalty', k_contact, ...
    'ModelType', 'FOM', ...
    'qdd0', qdd0); 
fprintf('Tempo Newmark: %.2f s\n', toc);

y_fom_newmark = DummyStruct.AssemblyObj.unconstrain_vector(q_fom_newmark);

% --- Estrazione intelligente e Salvataggio (Zero sprechi) ---
y_spia_fom_newmark = y_fom_newmark(dof_spia_global, :);
y_parete_fom_newmark = y_fom_newmark(dofs_parete_global, :);

save(fullfile(save_dir, 'FOM_Newmark.mat'), 't_fom_newmark', 'y_spia_fom_newmark', 'y_parete_fom_newmark');
fprintf('Salvato: FOM_Newmark.mat\n');

clear y_fom_newmark q_fom_newmark; 

%% --- 4. FOM TRANSIENT (ODE15S OTTIMIZZATO) ---
fprintf('\n=== Risoluzione FOM con ode15s (Versione Ottimizzata) ===\n');

% 1. Vettore di stato iniziale
y0_fom = [q0; qd0];

% 2. MODIF. 1: State-Space Mass Matrix
% Evitiamo di fare M \ (...) dentro la funzione. 
% Diciamo a MATLAB che l'equazione è M_state * y' = f(y)
M_state = blkdiag(speye(n_dofs_fom), Mc);

% 3. MODIF. 2: Jacobian Topology
I_sparse = speye(n_dofs_fom);
Zero_sparse = sparse(n_dofs_fom, n_dofs_fom);
J_pattern = [Zero_sparse, I_sparse; Kc ~= 0, Cc ~= 0];

% Inseriamo le nuove opzioni (Tol un po' più alte per velocizzare ulteriormente)
options = odeset('RelTol', 1e-4, 'AbsTol', 1e-5, 'MaxStep', t_shock/20, ...
                 'Mass', M_state, 'JPattern', J_pattern);

tic;
% NOTA: Passiamo K e C alla funzione, NON più M!
[t_fom_ode, y_fom_ode_raw] = ode15s(@(t, y) fom_state_space_fast(t, y, Kc, Cc, contact_dofs, gap_wall, k_contact, F_fom_handle), [0 tmax], y0_fom, options);
fprintf('Tempo ode15s (Ottimizzato): %.2f s\n', toc);

% Post-Processing ODE
q_fom_ode = y_fom_ode_raw(:, 1:n_dofs_fom)'; 
y_fom_ode_full = DummyStruct.AssemblyObj.unconstrain_vector(q_fom_ode);
y_spia_fom_ode = y_fom_ode_full(dof_spia_global, :);

save(fullfile(save_dir, 'FOM_ode15s.mat'), 't_fom_ode', 'y_spia_fom_ode');
fprintf('Salvato: FOM_ode15s.mat\n');

clear y_fom_ode_raw y_fom_ode_full q_fom_ode M_state J_pattern;

%% --- 6. FOM PLOT (CONFRONTO SOLUTORI) ---
fprintf('\n=== SEZIONE 6: Post-Processing e Confronto ===\n');

figure('Name', 'Integrator Comparison: Newmark vs ODE15s', 'Color', 'w', 'Units', 'normalized', 'Position', [0.1 0.3 0.5 0.5]);
plot(t_fom_newmark, y_spia_fom_newmark * 1e6, 'b-', 'LineWidth', 2, 'DisplayName', 'FOM (Newmark)');
hold on;
plot(t_fom_ode, y_spia_fom_ode * 1e6, 'r--', 'LineWidth', 2, 'DisplayName', 'FOM (ode15s)');

yline(gap_wall * 1e6, 'k:', 'LineWidth', 2, 'DisplayName', 'Wall');
title(sprintf('Impact Transient Comparison - Node #%d', nodo_spia));
xlabel('Time [s]'); 
ylabel('Displacement X [\mum]');
legend('Location', 'best'); 
grid on;




%% =====================================================================
%  FUNZIONI LOCALI (Da inserire in fondo al file)
%  ====================================================================

function f = fom_state_space_fast(t, y, K, C, contact_dofs, gap, k_pen, F_ext_handle)
    % Estrazione coordinate
    n = size(K, 1);
    q = y(1:n);
    qd = y(n+1:end);

    % Controllo Compenetrazione
    penetration = q(contact_dofs) - gap;
    is_pen = penetration > 0;

    % Vettore Forze di Penalità
    F_pen = zeros(n, 1);
    if any(is_pen)
        F_pen(contact_dofs(is_pen)) = k_pen * penetration(is_pen);
    end

    % Richiamo forzante al tempo t
    F_ext = F_ext_handle(t);

    % Costruzione del lato destro dell'equazione (Senza invertire M!)
    % Dato che MATLAB ora usa M_state * y' = f
    f = [qd; F_ext - K*q - C*qd - F_pen];
end