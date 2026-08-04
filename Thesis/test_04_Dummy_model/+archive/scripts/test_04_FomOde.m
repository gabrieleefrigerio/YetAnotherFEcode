%% =====================================================================
%  SCRIPT FOR SOLVER COMPARISON (ABSTOL EFFECT)
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
DummyStruct = AbaqusStructure();
% Assegna il file (assicurati che sia nella stessa cartella o metti il path completo)
DummyStruct.filename = 'DummyStructureAbaqus.inp'; 
% ASSICURATI che il tipo di elemento coincida con quello esportato.
DummyStruct.elementType = 'TRI3'; 
% Seleziona lo spessore fuori piano desiderato
DummyStruct.thickness = 10e-6; 

% Costruisci la geometria, le matrici K, M, C e applica i vincoli
DummyStruct.build();
% Plotta la struttura indeformata per un sanity check visivo
DummyStruct.plot_undeformed();
% Calcola e plotta i primi 5 modi propri
DummyStruct.compute_eigenmodes(5);

Mc = DummyStruct.AssemblyObj.constrain_matrix(DummyStruct.M);
Kc = DummyStruct.AssemblyObj.constrain_matrix(DummyStruct.K);
DummyStruct.compute_rayleigh_damping(1000, 1000);
Cc = DummyStruct.AssemblyObj.constrain_matrix(DummyStruct.C);

n_dofs_fom = size(Mc, 1);
n_of_per_node = DummyStruct.MeshObj.nDOFPerNode;

% Contact Parameters
gap_wall = 1.5e-6; 
k_contact = max(diag(Kc)) * 10; 
% Estrazione sicura dei GdL per il contatto tramite Node Set
contact_dofs = DummyStruct.get_contact_dofs(1); 

%% --- Identificazione Nodi e GdL per Estrazione Dati ---
nodi_parete = DummyStruct.contact_nodes;
coordinate_Y_parete = DummyStruct.nodes(nodi_parete, 2); 
[~, idx_top] = max(coordinate_Y_parete);
nodo_spia = nodi_parete(idx_top); 

dof_spia_global = (nodo_spia - 1) * n_of_per_node + 1;
dofs_parete_global = (nodi_parete - 1) * n_of_per_node + 1;

%% --- 2. TRANSIENT SETUP ---
tmax = (1e-3)/2; 
% Initial Conditions
q0 = zeros(n_dofs_fom, 1);
qd0 = zeros(n_dofs_fom, 1);

% Shock Parameters
shock_g = 1e5;                 
shock_amp = shock_g * 9.81;    
t_shock = 10e-7;               

% Building Shock
dir_vector = zeros(n_dofs_fom, 1);
dir_vector(1:n_of_per_node:n_dofs_fom) = 1; 

F_spatial_fom = Mc * dir_vector;
F_fom_handle = @(t) F_spatial_fom * shock_amp * sin(pi * t / t_shock) * (t <= t_shock);

%% --- 4. FOM TRANSIENT (CONFRONTO ABSTOL) ---
fprintf('\n=== Risoluzione FOM con ode15s (Confronto AbsTol) ===\n');
y0_fom = [q0; qd0];

% State-Space Mass Matrix e Jacobian Topology
M_state = blkdiag(speye(n_dofs_fom), Mc);
I_sparse = speye(n_dofs_fom);
Zero_sparse = sparse(n_dofs_fom, n_dofs_fom);
J_pattern = [Zero_sparse, I_sparse; Kc ~= 0, Cc ~= 0];

% Array delle AbsTol da testare (mantenendo RelTol a 1e-4 e MaxStep fisso)
array_abs_tol = [1e-5, 1e-9];
nomi_legenda = {'FOM: AbsTol = 1e-5 (Larga)', 'FOM: AbsTol = 1e-9 (Stretta)'};
colori = {'r--', 'b-'};

% Celle per salvare i risultati da plottare
t_results = cell(length(array_abs_tol), 1);
y_results = cell(length(array_abs_tol), 1);
cpu_times = zeros(length(array_abs_tol), 1);

for i = 1:length(array_abs_tol)
    current_abs_tol = array_abs_tol(i);
    fprintf('\n-> Esecuzione Run %d: AbsTol = %g...\n', i, current_abs_tol);
    
    options = odeset('RelTol', 1e-4, 'AbsTol', current_abs_tol, 'MaxStep', 0.5e-8, ...
                     'Mass', M_state, 'JPattern', J_pattern);
    tic;
    
    % Integrazione
    [t_fom_ode, y_fom_ode_raw] = ode15s(@(t, y) fom_state_space_fast(t, y, Kc, Cc, contact_dofs, gap_wall, k_contact, F_fom_handle), [0 tmax], y0_fom, options);
    cpu_times(i) = toc;
    fprintf('   Tempo ode15s: %.2f s\n', cpu_times(i));
    
    % Post-Processing ODE
    q_fom_ode = y_fom_ode_raw(:, 1:n_dofs_fom)'; 
    y_fom_ode_full = DummyStruct.AssemblyObj.unconstrain_vector(q_fom_ode);
    
    % Estrazione Dati
    y_parete_fom_ode = y_fom_ode_full(dofs_parete_global, :);
    y_spia_fom_ode = y_fom_ode_full(dof_spia_global, :);
    
    % Salvataggio in memoria per il plot
    t_results{i} = t_fom_ode;
    y_results{i} = y_spia_fom_ode;
    
    % Salvataggio file
    file_name = sprintf('FOM_ode15s_AbsTol_%e.mat', current_abs_tol);
    save(fullfile(save_dir, file_name), 't_fom_ode', 'y_parete_fom_ode', 'y_spia_fom_ode', 'current_abs_tol');
end
clear y_fom_ode_raw y_fom_ode_full q_fom_ode M_state J_pattern;

%% --- 6. FOM PLOT (RISULTATI TRANSITORIO) ---
fprintf('\n=== SEZIONE 6: Post-Processing Spostamenti ===\n');
figure('Name', 'Confronto Effetto AbsTol', 'Color', 'w', 'Units', 'normalized', 'Position', [0.1 0.2 0.6 0.6]);
hold on; grid on;

% Plotta il limite fisico della parete
yline(gap_wall * 1e6, 'k:', 'LineWidth', 2, 'DisplayName', 'Wall Gap (1.5 \mum)');

% Plotta i risultati ciclando sulle AbsTol effettive
for i = 1:length(array_abs_tol)
    disp_name = sprintf('%s - CPU: %.1fs', nomi_legenda{i}, cpu_times(i));
    plot(t_results{i}, y_results{i} * 1e6, colori{i}, 'LineWidth', 1.5, 'DisplayName', disp_name);
end

% Formattazione del grafico
title(sprintf('Impact Transient Response - Nodo #%d', nodo_spia), 'FontSize', 14);
xlabel('Time [s]', 'FontSize', 12); 
ylabel('Displacement X [\mum]', 'FontSize', 12);
legend('Location', 'best', 'FontSize', 10); 

%% --- 7. PLOT DELLA FORZA DI CONTATTO ---
fprintf('\n=== SEZIONE 7: Calcolo Forza di Contatto ===\n');
figure('Name', 'Confronto Forza di Contatto', 'Color', 'w', 'Units', 'normalized', 'Position', [0.6 0.2 0.4 0.6]);
hold on; grid on;

% Plotta la forza calcolata per i due run
for i = 1:length(array_abs_tol)
    % 1. Calcolo la compenetrazione del nodo spia (Spostamento - Gap)
    compenetrazione = y_results{i} - gap_wall;
    
    % 2. La molla non "tira" mai, quindi azzero tutte le compenetrazioni negative
    compenetrazione(compenetrazione < 0) = 0;
    
    % 3. Calcolo la Forza (k * delta_x)
    forza_contatto = k_contact .* compenetrazione;
    
    % 4. Plot
    disp_name = sprintf('%s - CPU: %.1fs', nomi_legenda{i}, cpu_times(i));
    plot(t_results{i}, forza_contatto, colori{i}, 'LineWidth', 1.5, 'DisplayName', disp_name);
end

% Formattazione del grafico
title(sprintf('Contact Force Transient - Nodo #%d', nodo_spia), 'FontSize', 14);
xlabel('Time [s]', 'FontSize', 12); 
ylabel('Contact Force [N]', 'FontSize', 12);
legend('Location', 'best', 'FontSize', 10);

%% --- 8. KINEMATICS ANALYSIS (VELOCITY & ACCELERATION) ---
fprintf('\n=== SEZIONE 8: Analisi Cinematica al Primo Impatto ===\n');
% Prendiamo i risultati del Run più accurato (il secondo, AbsTol = 1e-9)
t_plot = t_results{2};
q_spia = y_results{2}; 

% 1. Calcolo numerico Velocità e Accelerazione tramite gradiente centrale
v_spia = gradient(q_spia, t_plot);
a_spia = gradient(v_spia, t_plot);

% 2. Trova l'indice del primo istante in cui si tocca la parete
idx_impatto = find(q_spia >= gap_wall, 1);
if ~isempty(idx_impatto)
    t_impatto = t_plot(idx_impatto);
    
    % Per sicurezza, prendiamo i valori un istante prima della compenetrazione
    idx_pre = idx_impatto - 1; 
    if idx_pre < 1; idx_pre = 1; end
    
    v_impatto = v_spia(idx_pre);
    a_impatto_g = a_spia(idx_pre) / 9.81;
    
    % 3. Calcolo dell'accelerazione massima subita in volo
    [max_acc_volo, idx_max_a] = max(a_spia(1:idx_pre));
    max_acc_volo_g = max_acc_volo / 9.81;
    
    fprintf('--- DIAGNOSTICA PRIMO URTO (Dati da AbsTol = 1e-9) ---\n');
    fprintf('Tempo del primo impatto       : %.3e s\n', t_impatto);
    fprintf('Velocita'' di schianto         : %.3f m/s\n', v_impatto);
    fprintf('Accelerazione istante urto    : %.1f g\n', a_impatto_g);
    fprintf('Picco di Accelerazione subito : %.1f g (al tempo %.3e s)\n', max_acc_volo_g, t_plot(idx_max_a));
    fprintf('------------------------------------------------------\n');
else
    fprintf('ATTENZIONE: La massa non ha mai toccato la parete nel tempo simulato!\n');
end

% Plot per controllo visivo (riferito ad AbsTol = 1e-9)
figure('Name', 'Velocità e Accelerazione (Fini)', 'Color', 'w', 'Units', 'normalized', 'Position', [0.3 0.2 0.4 0.6]);
subplot(2,1,1);
plot(t_plot, v_spia, 'b', 'LineWidth', 1.5);
if ~isempty(idx_impatto), xline(t_impatto, 'r--', 'Impatto'); end
title('Velocità Nodo Spia (AbsTol 1e-9)'); ylabel('Velocità [m/s]'); grid on;

subplot(2,1,2);
plot(t_plot, a_spia / 9.81, 'k', 'LineWidth', 1.5);
if ~isempty(idx_impatto), xline(t_impatto, 'r--', 'Impatto'); end
title('Accelerazione Nodo Spia (AbsTol 1e-9)'); ylabel('Accelerazione [g]'); xlabel('Tempo [s]'); grid on;

%% =====================================================================
%  FUNZIONI LOCALI
%  ====================================================================
function f = fom_state_space_fast(t, y, K, C, contact_dofs, gap, k_pen, F_ext_handle)
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
    
    % Costruzione della derivata di stato
    f = [qd; F_ext - K*q - C*qd - F_pen];
end