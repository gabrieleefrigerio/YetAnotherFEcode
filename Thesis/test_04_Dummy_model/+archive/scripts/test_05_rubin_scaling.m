%% =====================================================================
%  TEST SCALING BASE DI RUBIN  (A/B/C)
% ---------------------------------------------------------------------
%  Domanda: il divario di CPU Rubin vs Milman-Chu e' dovuto alla
%  parametrizzazione (lunghezze delle colonne incommensurabili) o al
%  metodo?
%
%  check_ortho ha gia' misurato, sulla base di Rubin:
%      spread lunghezze = 1.57e+06     -> problema di SCALA
%      cond(G) direzioni = 2.84e+01    -> direzioni SANE
%  quindi kappa(Mr) ~ (1.57e6)^2 * 28.4 ~ 7e13 e' interamente artefatto
%  delle unita' di misura (interfaccia in [m], modi in [kg^-1/2 m]).
%
%  Tre configurazioni, per separare l'effetto scaling dall'effetto
%  "cambio ramo del solutore":
%     A  grezza  + ramo 'Rubin' (contatto su coordinate ridotte)  -> baseline
%     B  grezza  + ramo 'MC'    (contatto proiettato)             -> costo ramo
%     C  scalata + ramo 'MC'                                      -> effetto scaling
%  L'effetto scaling e' C vs B. Il costo del ramo e' B vs A.
%
%  ATTESO su C:  cond(Mr) -> ~28 | nfailed% 4.6 -> ~1 | LU/Jac 9.7 -> ~40
%  INVARIANTE:   GRE identico in A, B, C (cambio base invertibile).
%                Se il GRE cambia, e' un BUG (indicizzazione o CI).
% =====================================================================
clear; close all; clc;

%% --- 1. INPUTS ---
array_linModes = 200;
array_QFactor  = 1000;
array_k_mult   = 10;

% Se true, tutte le configurazioni girano con AbsTol scalare identica.
% Serve perche' la AbsTol energetica usa diag(Mr) e diag(Kr), che NON sono
% invarianti per cambio base: A/B e C riceverebbero criteri d'errore diversi.
% Richiede il parametro 'AbsTolOverride' nel solutore (patch di 4 righe).
use_fixed_abstol = false;

impulse_g         = 1e5;
impulse_amp       = impulse_g * 9.81;
impulse_angle_deg = 0;
impulse_sign      = 1;
t_shock = 10e-7;
dt      = 0.5e-8;
tmax    = (1e-3)/2;

gap_T =  5e-6;    gap_B = -1.5e-6;
gap_L = -5e-6;    gap_R =  1.5e-6;

bench_RelTol = 1e-9;

%% --- 2. SETUP DIRECTORY ---
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM');
save_dir  = fullfile('results_V2', sprintf('RubinScaling_%s', timestamp));
if ~exist(save_dir, 'dir'), mkdir(save_dir); end
fprintf('Directory creata: %s\n\n', save_dir);

%% --- 3. STRUTTURA BASE ---
fprintf('Setting up the model...\n');
DummyStruct = AbaqusStructure_V2();
DummyStruct.filename    = 'DummyStructureAbaqus_V4.inp';
DummyStruct.elementType = 'TRI3';
DummyStruct.build();

max_phi = max(array_linModes);
fprintf('Estrazione modale di %d modi...\n', max_phi);
DummyStruct.compute_eigenmodes(max_phi);

Mc = DummyStruct.AssemblyObj.constrain_matrix(DummyStruct.M);
Kc = DummyStruct.AssemblyObj.constrain_matrix(DummyStruct.K);
n_dofs_fom = size(Mc, 1);
k_base = max(diag(Kc));

%% --- 3.1 INTERFACCE DI CONTATTO ---
dofs_T = DummyStruct.get_specific_contact_dofs('T', 2);
dofs_B = DummyStruct.get_specific_contact_dofs('B', 2);
dofs_L = DummyStruct.get_specific_contact_dofs('L', 1);
dofs_R = DummyStruct.get_specific_contact_dofs('R', 1);

% L'ordine fissato qui e' l'UNICA verita': RomRubinTest lo eredita in
% contactDofs, quindi nelle colonne di Psi_res, nelle coordinate ridotte
% e in scaleD. gaps_array e' costruito nello stesso ordine per costruzione.
contact_dofs = [dofs_T; dofs_B; dofs_L; dofs_R];
gaps_array   = [repmat(gap_T, length(dofs_T), 1); ...
                repmat(gap_B, length(dofs_B), 1); ...
                repmat(gap_L, length(dofs_L), 1); ...
                repmat(gap_R, length(dofs_R), 1)];
n_bnd = length(contact_dofs);

assert(numel(unique(contact_dofs)) == n_bnd, 'contact_dofs contiene duplicati.');
assert(numel(gaps_array) == n_bnd, 'gaps_array disallineato da contact_dofs.');

Interfaces = struct();
labels     = {'T', 'B', 'L', 'R'};
nodes_list = {DummyStruct.contact_nodes_T, DummyStruct.contact_nodes_B, ...
              DummyStruct.contact_nodes_L, DummyStruct.contact_nodes_R};
for i = 1:length(labels)
    l = labels{i};  n = nodes_list{i};
    Interfaces.(l).nodes = n;
    if ~isempty(n)
        Interfaces.(l).global_X = (n - 1) * DummyStruct.MeshObj.nDOFPerNode + 1;
        Interfaces.(l).global_Y = (n - 1) * DummyStruct.MeshObj.nDOFPerNode + 2;
        Interfaces.(l).coord_X  = DummyStruct.nodes(n, 1);
        Interfaces.(l).coord_Y  = DummyStruct.nodes(n, 2);
    end
end

%% --- 3.2 FORZANTE E CONDIZIONI INIZIALI ---
n_dof_per_node = DummyStruct.MeshObj.nDOFPerNode;
impulse_dir    = impulse_sign * [cosd(impulse_angle_deg); sind(impulse_angle_deg)];

dir_vector = zeros(n_dofs_fom, 1);
dir_vector(1:n_dof_per_node:n_dofs_fom) = impulse_dir(1);
dir_vector(2:n_dof_per_node:n_dofs_fom) = impulse_dir(2);

F_spatial_fom = Mc * dir_vector;
F_fom_handle  = @(t) F_spatial_fom * impulse_amp * sin(pi * t / t_shock) * (t <= t_shock);

q0  = zeros(n_dofs_fom, 1);
qd0 = zeros(n_dofs_fom, 1);

v_max = impulse_amp * 2 * t_shock / pi;
m_eff = dir_vector' * Mc * dir_vector;
Eref  = 0.5 * m_eff * v_max^2;          % scalare fisico: NON dipende dalla base

t_common = 0 : 10*dt : tmax;

fprintf('Impulso: %.1e g @ %.1f deg (sign %+d)\n', impulse_g, impulse_angle_deg, impulse_sign);
fprintf('Eref = %.4e J   (v_max = %.4f m/s)\n', Eref, v_max);
fprintf('Interfaccia: %d DOF di contatto\n\n', n_bnd);

%% --- 4. CONFIGURAZIONI DEL TEST ---
cfgs = struct( ...
    'name',    {'A_raw_std',  'B_raw_proj', 'C_scaled_proj'}, ...
    'scaling', {false,        false,        true          }, ...
    'branch',  {'std',        'proj',       'proj'        }, ...
    'enable',  {true,         true,         true          });

%% --- 5. RUN ---
for Q = array_QFactor
    DummyStruct.compute_rayleigh_damping(Q, Q);

    for k_mult = array_k_mult
        k_contact = k_base * k_mult;

        % AbsTol scalare comune (opzionale): calibrata sull'energia di
        % riferimento e sulla rigidezza di contatto, uguale per tutte le cfg.
        epsE = 1e-2 * bench_RelTol;
        abstol_fixed = epsE * min(sqrt(2*Eref / k_contact), sqrt(2*Eref / m_eff));

        for phi = array_linModes
            fprintf('\n========= Rubin scaling test | Phi %d | Q %d | Kmult %g =========\n', ...
                    phi, Q, k_mult);

            % ---- build UNA sola volta (eigs domina il costo offline) ----
            tic_offline = tic;
            rom = RomRubinTest(DummyStruct, phi, contact_dofs, false);
            rom.build();
            offline_time = toc(tic_offline);
            fprintf(' Offline build: %.2f s\n', offline_time);

            % spettro della flessibilita' residua d'interfaccia:
            % serve gia' ora come input per l'interface reduction
            s_if = rom.diag_info.sigma_iface;
            fprintf(' sigma interfaccia (norm.): ');
            fprintf('%.2e ', s_if(1:min(end,10)));  fprintf('...\n');
            fprintf(' rango effettivo interfaccia @1e-6: %d su %d\n', ...
                    sum(s_if > 1e-6), numel(s_if));

            for c = 1:numel(cfgs)
                cfg = cfgs(c);
                if ~cfg.enable, continue; end

                fprintf('\n>>> CONFIG %s  (scaling %d, ramo %s)\n', ...
                        cfg.name, cfg.scaling, cfg.branch);

                % lo scaling rende le coordinate d'interfaccia pari a
                % d_i * q_fisico: il ramo standard le leggerebbe come
                % spostamenti fisici, in silenzio e con risultati sbagliati.
                assert(~(cfg.scaling && strcmpi(cfg.branch, 'std')), ...
                    'Config %s: scaling attivo richiede il ramo proiettato.', cfg.name);

                rom.set_scaling(cfg.scaling);
                [Mr, Kr, Cr] = rom.get_reduced_matrices();
                Pc_r = rom.Pc;

                % --- CI proiettate: forma valida anche quando Mr ~= I ---
                [q0_r, qd0_r] = rom.project_ic(q0, qd0);
                F_r_handle = @(t) Pc_r' * F_fom_handle(t);

                % --- selezione ramo e semantica di ContactTargetDOF -----
                %   'std'  -> indici di COORDINATE RIDOTTE (1:n_bnd)
                %   'proj' -> indici di DOF FISICI constrained (contact_dofs)
                % In entrambi i casi gaps_array e' nell'ordine di contact_dofs,
                % che per costruzione di T2 coincide con l'ordine delle
                % prime n_bnd coordinate ridotte.
                switch lower(cfg.branch)
                    case 'std'
                        model_type  = 'Rubin';
                        target_dofs = (1:n_bnd)';
                        extra_args  = {};
                    case 'proj'
                        model_type  = 'MC';
                        target_dofs = contact_dofs;
                        extra_args  = {'ProjectionMatrix', Pc_r};
                end

                if use_fixed_abstol
                    extra_args = [extra_args, {'AbsTolOverride', abstol_fixed}]; %#ok<AGROW>
                end

                tic;
                solverROM = TransientSolverOde_V2(Mr, Kr, Cr);
                [t_rom, q_rom] = solverROM.solve(tmax, dt, q0_r, qd0_r, F_r_handle, ...
                    'ContactTargetDOF', target_dofs, ...
                    'ContactGap',       gaps_array, ...
                    'ContactPenalty',   k_contact, ...
                    'ModelType',        model_type, ...
                    'Eref',             Eref, ...
                    'RelTol',           bench_RelTol, ...
                    'OutputTimes',      t_common, ...
                    extra_args{:});
                cpu_time = toc;

                % --- ricostruzione fisica (valida per qualunque scaling) ---
                q_c        = Pc_r * q_rom;
                y_rom_full = DummyStruct.AssemblyObj.unconstrain_vector(q_c);

                y_contact = struct();
                for loc = labels
                    l = loc{1};
                    if ~isempty(Interfaces.(l).nodes)
                        y_contact.(l).X = y_rom_full(Interfaces.(l).global_X, :);
                        y_contact.(l).Y = y_rom_full(Interfaces.(l).global_Y, :);
                    end
                end

                % --- amplificazione delle coordinate ridotte -------------
                % ||y|| / ||Pc*y||_M : misura la cancellazione. Se e' grande,
                % il controllo d'errore su y non protegge la ricostruzione q.
                nrm_y  = max(vecnorm(q_rom));
                nrm_Py = max(sqrt(abs(sum(q_c .* (Mc * q_c), 1))));
                amp    = nrm_y / max(nrm_Py, realmin);

                diag_info = rom.diag_info;
                diag_info.cfg_name   = cfg.name;
                diag_info.branch     = cfg.branch;
                diag_info.amp_coords = amp;
                diag_info.cpu_time   = cpu_time;

                fprintf('  CPU %.2f s | cond(Mr) %.3e | ||y||/||Py||_M %.3e\n', ...
                        cpu_time, diag_info.cond_Mr, amp);

                file_name = sprintf('ROM_Rubin_%s_Phi%03d_Q%04d_K%g.mat', ...
                                    cfg.name, phi, Q, k_mult);
                save(fullfile(save_dir, file_name), 't_rom', 'y_contact', ...
                     'Interfaces', 'cpu_time', 'offline_time', 'diag_info');
            end
        end
    end
end

fprintf('\n=== Test completato. Risultati in %s ===\n', save_dir);
fprintf('Da confrontare:  C vs B = effetto scaling | B vs A = costo del ramo\n');
fprintf('Il GRE deve essere IDENTICO nelle tre configurazioni.\n');
