%% =====================================================================
%  MAIN SIMULAZIONI - benchmark FOM vs ROM con contatto unilaterale
%
%  Main unico, adattivo nel numero di interfacce di contatto. Le interfacce
%  vengono lette dal file .inp (node set 'ContactInterface[_<label>]') e
%  dichiarate qui sotto in cfg.interfaces con direzione e gap firmato.
%  Aggiungere o togliere una riga a quella tabella e' l'unica modifica
%  necessaria per passare da un modello all'altro.
%
%  Metodi disponibili:
%    FOM     modello completo, contatto a penalita' (ode15s)
%    MT      troncamento modale, penalita' proiettata (ode15s)
%    MC      Milman-Chu, penalita' proiettata (ode15s)
%    Rubin   CMS free-interface, penalita' (ode15s)
%    MCB     Craig-Bampton massless, contatto set-valued esatto (LCP+leapfrog)
%    MN      MacNeal massless,       contatto set-valued esatto (LCP+leapfrog)
%
%  I risultati vanno in results/<nome_test>_<timestamp>/, insieme a
%  run_config.mat che contiene la configurazione completa: il post-processing
%  la rilegge da li' e non ha bisogno di sapere nulla a priori sul modello.
% =====================================================================
clear; close all; clc;

%% --- 1. CONFIGURAZIONE ------------------------------------------------

% --- Modello ---
cfg.mesh_file    = 'DummyStructureAbaqus_V4.inp';
cfg.element_type = 'TRI3';

% --- Interfacce di contatto ---
%   etichetta | direzione (1 = X, 2 = Y) | gap firmato [m]
% Il SEGNO del gap indica da che parte sta il muro rispetto all'origine:
%   gap > 0  muro nella direzione positiva del GdL  (penetra se q >  gap)
%   gap < 0  muro nella direzione negativa del GdL  (penetra se q <  gap)
% Le etichette devono esistere nel .inp come *Nset, nset=ContactInterface_<label>.
% Per un modello a interfaccia unica (*Nset, nset=ContactInterface) l'etichetta
% e' 'C' e la tabella si riduce a una sola riga.
cfg.interfaces = { ...
    'T', 2,  5.0e-6 ; ...   % muro in Y positivo
    'B', 2, -1.5e-6 ; ...   % muro in Y negativo
    'L', 1, -5.0e-6 ; ...   % muro in X negativo
    'R', 1,  1.5e-6 };      % muro in X positivo

% --- Metodi da eseguire ---
cfg.run.FOM   = 1;
cfg.run.MT    = 1;
cfg.run.MC    = 1;
cfg.run.Rubin = 1;
cfg.run.MCB   = 0;
cfg.run.MN    = 0;

% --- Sweep dei parametri ---
cfg.array_linModes = [10, 50, 100, 150, 200];
cfg.array_QFactor  = [1000];
cfg.array_k_mult   = [10];        % moltiplicatore rigidezza di contatto
                                  % (ignorato dai modelli massless MCB/MN)

% --- Forzante impulsiva ---
cfg.impulse_g         = 1e5;      % ampiezza [g]
cfg.impulse_angle_deg = 0;        % direzione nel piano XY [deg]
cfg.impulse_sign      = 1;        % verso (+1 / -1)
cfg.t_shock           = 10e-7;    % durata del semiseno [s]

% --- Integrazione ---
cfg.dt      = 0.5e-8;
cfg.tmax    = (1e-3)/2;
cfg.RelTol  = 1e-9;               % ode15s, ROM
cfg.RelTolFOM = 1e-10;            % ode15s, FOM (riferimento piu' stretto)
cfg.output_stride = 10;           % output ogni N passi dt (griglia comune)

% --- Nome del test ---
cfg.test_name = sprintf('Shock_%dg_%.1es', cfg.impulse_g, cfg.tmax);

%% --- 2. DIRECTORY DEI RISULTATI ---------------------------------------
timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd_HH-mm'));
save_dir  = fullfile('results', sprintf('%s_%s', cfg.test_name, timestamp));
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end
fprintf('Directory dei risultati: %s\n\n', save_dir);

%% --- 3. MODELLO --------------------------------------------------------
fprintf('Costruzione del modello...\n');
Struct = AbaqusStructure();
Struct.filename    = cfg.mesh_file;
Struct.elementType = cfg.element_type;
Struct.build();
Struct.describe_interfaces();

max_phi = max(cfg.array_linModes);
fprintf('Estrazione di %d modi...\n', max_phi);
Struct.compute_eigenmodes(max_phi);

Mc = Struct.AssemblyObj.constrain_matrix(Struct.M);
Kc = Struct.AssemblyObj.constrain_matrix(Struct.K);
n_dofs_fom = size(Mc, 1);
k_base     = max(diag(Kc));

%% --- 4. INTERFACCE DI CONTATTO ----------------------------------------
% Da cfg.interfaces si ricavano, in un colpo solo:
%   contact_dofs  GdL vincolati di contatto, concatenati
%   gaps_array    gap firmato per ciascun GdL
%   Interfaces    metadati per il post-processing (nodi, GdL globali, coord.)
labels     = cfg.interfaces(:, 1)';
dirs       = cell2mat(cfg.interfaces(:, 2))';
gaps_iface = cell2mat(cfg.interfaces(:, 3))';

% Controllo: ogni etichetta dichiarata deve esistere nel file .inp
missing = setdiff(labels, Struct.contact_labels);
if ~isempty(missing)
    error('MAIN:MissingInterface', ...
        ['Le interfacce {%s} non esistono in %s.\n' ...
         'Interfacce disponibili nel file: {%s}'], ...
        strjoin(missing, ', '), cfg.mesh_file, strjoin(Struct.contact_labels, ', '));
end
unused = setdiff(Struct.contact_labels, labels);
if ~isempty(unused)
    fprintf('[nota] Interfacce presenti nel .inp ma non usate: %s\n', strjoin(unused, ', '));
end

contact_dofs = [];
gaps_array   = [];
Interfaces   = struct();
nDOFPerNode  = Struct.MeshObj.nDOFPerNode;

for i = 1:numel(labels)
    lbl = labels{i};
    d   = Struct.get_contact_dofs(lbl, dirs(i));
    if isempty(d)
        warning('MAIN:EmptyInterface', ...
            'Interfaccia %s: nessun GdL libero (tutti i nodi sono vincolati). Saltata.', lbl);
        continue;
    end

    contact_dofs = [contact_dofs; d];                       %#ok<AGROW>
    gaps_array   = [gaps_array;   gaps_iface(i)*ones(numel(d), 1)]; %#ok<AGROW>

    n = Struct.get_contact_nodes(lbl);
    Interfaces.(lbl).nodes    = n;
    Interfaces.(lbl).dir      = dirs(i);
    Interfaces.(lbl).gap      = gaps_iface(i);
    Interfaces.(lbl).global_X = (n - 1) * nDOFPerNode + 1;
    Interfaces.(lbl).global_Y = (n - 1) * nDOFPerNode + 2;
    Interfaces.(lbl).coord_X  = Struct.nodes(n, 1);
    Interfaces.(lbl).coord_Y  = Struct.nodes(n, 2);
end
active_labels = fieldnames(Interfaces)';

if isempty(contact_dofs)
    error('MAIN:NoContact', 'Nessun GdL di contatto attivo: controllare cfg.interfaces.');
end
fprintf('\nContatto: %d interfacce, %d GdL totali\n', numel(active_labels), numel(contact_dofs));

%% --- 5. FORZANTE E CONDIZIONI INIZIALI --------------------------------
if nDOFPerNode < 2
    error('MAIN:Not2D', 'Il modello non ha abbastanza GdL per un impulso nel piano.');
end

impulse_amp = cfg.impulse_g * 9.81;
impulse_dir = cfg.impulse_sign * [cosd(cfg.impulse_angle_deg); sind(cfg.impulse_angle_deg)];

dir_vector = zeros(n_dofs_fom, 1);
dir_vector(1:nDOFPerNode:n_dofs_fom) = impulse_dir(1);   % GdL X
dir_vector(2:nDOFPerNode:n_dofs_fom) = impulse_dir(2);   % GdL Y

F_spatial_fom = Mc * dir_vector;
F_fom_handle  = @(t) F_spatial_fom * impulse_amp * sin(pi*t/cfg.t_shock) * (t <= cfg.t_shock);

q0  = zeros(n_dofs_fom, 1);
qd0 = zeros(n_dofs_fom, 1);

% Energia di riferimento per la AbsTol pesata in energia
v_max = impulse_amp * 2 * cfg.t_shock / pi;
m_eff = dir_vector' * Mc * dir_vector;
Eref  = 0.5 * m_eff * v_max^2;

% Griglia di output comune a tutti i modelli: il post-processing puo' cosi'
% confrontare le time history senza interpolare.
t_common = 0 : cfg.output_stride*cfg.dt : cfg.tmax;

fprintf('Impulso: %.1e g @ %.1f deg (verso %+d)\n', ...
    cfg.impulse_g, cfg.impulse_angle_deg, cfg.impulse_sign);
fprintf('Eref = %.4e J   (v_max = %.4f m/s)\n', Eref, v_max);

% Configurazione salvata una sola volta: e' il contratto col post-processing
save(fullfile(save_dir, 'run_config.mat'), ...
    'cfg', 'Interfaces', 'active_labels', 'contact_dofs', 'gaps_array', ...
    'Eref', 't_common', 'n_dofs_fom', 'k_base');

%% --- 6. FOM ------------------------------------------------------------
if cfg.run.FOM
    fprintf('\n=========================================\n');
    fprintf('                 FOM\n');
    fprintf('=========================================\n');
    for Q = cfg.array_QFactor
        Struct.compute_rayleigh_damping(Q, Q);
        Cc = Struct.AssemblyObj.constrain_matrix(Struct.C);

        for k_mult = cfg.array_k_mult
            k_contact = k_base * k_mult;
            fprintf('\n[FOM] Q = %d | k_mult = %g\n', Q, k_mult);

            tic;
            solver = TransientSolverOde(Mc, Kc, Cc);
            [t, q] = solver.solve(cfg.tmax, cfg.dt, q0, qd0, F_fom_handle, ...
                'ContactTargetDOF', contact_dofs, ...
                'ContactGap',       gaps_array, ...
                'ContactPenalty',   k_contact, ...
                'ModelType',        'FOM', ...
                'Eref',             Eref, ...
                'RelTol',           cfg.RelTolFOM, ...
                'OutputTimes',      t_common);
            cpu_time = toc;

            y_contact = extract_contact_response(Struct, Interfaces, active_labels, q);

            model = 'FOM'; n_modes = n_dofs_fom; offline_time = 0;
            save(fullfile(save_dir, sprintf('FOM_Q%04d_K%g.mat', Q, k_mult)), ...
                't', 'y_contact', 'Interfaces', 'cpu_time', 'offline_time', ...
                'model', 'n_modes', 'Q', 'k_mult');
        end
    end
end

%% --- 7. ROM ------------------------------------------------------------
rom_list = {'MT', 'MC', 'Rubin', 'MCB', 'MN'};
rom_list = rom_list(cellfun(@(m) cfg.run.(m), rom_list) == 1);

if ~isempty(rom_list)
    fprintf('\n=========================================\n');
    fprintf('                 ROM\n');
    fprintf('=========================================\n');
end

for Q = cfg.array_QFactor

    % compute_rayleigh_damping fa due cose:
    %   (a) aggiorna Struct.C          -> serve ai ROM a penalita' (MT/MC/Rubin)
    %   (b) restituisce alpha e beta   -> servono ai ROM massless, che ne
    %       ricavano uno smorzamento modale diagonale equivalente
    % Usando gli stessi alpha/beta tutti i ROM hanno lo stesso smorzamento.
    [~, alpha_ray, beta_ray] = Struct.compute_rayleigh_damping(Q, Q);
    rayleigh = struct('alpha', alpha_ray, 'beta', beta_ray);

    % Diagnostica: il Rayleigh sovrasmorza i modi alti.
    w_hi = 2*pi * Struct.frequencies(max_phi);
    z_hi = 0.5*(alpha_ray/w_hi + beta_ray*w_hi);
    fprintf('  zeta(modo %d) = %.4f | zeta target (modi 1-2) = %.4f\n', ...
        max_phi, z_hi, 1/(2*Q));
    if z_hi > 1
        warning('MAIN:Overdamped', ...
            ['Rayleigh rende SOVRACRITICI i modi alti (zeta_%d = %.2f). ' ...
             'Questo affligge TUTTI i ROM, non solo i massless.'], max_phi, z_hi);
    end

    for k_mult = cfg.array_k_mult
        k_contact = k_base * k_mult;

        for phi = cfg.array_linModes
            for im = 1:numel(rom_list)
                model = rom_list{im};

                % I modelli massless usano contatto set-valued esatto: la
                % rigidezza di penalita' non li riguarda, quindi il ciclo su
                % k_mult sarebbe degenere. Girano solo al primo valore.
                is_massless = any(strcmp(model, {'MCB', 'MN'}));
                if is_massless && k_mult ~= cfg.array_k_mult(1)
                    continue;
                end

                fprintf('\n--- %s | Phi %d | Q %d | k_mult %g ---\n', model, phi, Q, k_mult);

                % ---------- costruzione della base ----------
                tic_offline = tic;
                switch model
                    case 'MT'
                        rom = RomMC(Struct, phi, contact_dofs, k_contact, 0);  % include_MC = 0
                        rom.build();
                    case 'MC'
                        rom = RomMC(Struct, phi, contact_dofs, k_contact, 1);
                        rom.build();
                    case 'Rubin'
                        rom = RomRubin(Struct, phi, contact_dofs);
                        rom.build();
                    case 'MCB'
                        rom = RomMCB(Struct, phi, contact_dofs);
                        rom.build(rayleigh);
                        rom.check();
                    case 'MN'
                        rom = RomMN(Struct, phi, contact_dofs);
                        rom.build(rayleigh);
                        rom.check();
                end
                [Mr, Kr, Cr] = rom.get_reduced_matrices();

                % Matrice di proiezione sui GdL vincolati.
                % MT e MC espongono P sui GdL globali, gli altri Pc gia' vincolata.
                if any(strcmp(model, {'MT', 'MC'}))
                    Pc = zeros(n_dofs_fom, size(rom.P, 2));
                    for ic = 1:size(Pc, 2)
                        Pc(:, ic) = Struct.AssemblyObj.constrain_vector(rom.P(:, ic));
                    end
                else
                    Pc = rom.Pc;
                end
                offline_time = toc(tic_offline);

                % ---------- condizioni iniziali e forzante ridotte ----------
                if any(q0)
                    q0_r = Pc \ q0;
                else
                    q0_r = zeros(size(Pc, 2), 1);   % q0 nullo: nessuna proiezione
                end
                qd0_r    = zeros(size(Pc, 2), 1);
                F_handle = @(t) Pc' * F_fom_handle(t);

                % ---------- integrazione ----------
                lambda = []; info = [];
                tic;
                if is_massless
                    % Contatto set-valued esatto (Monjaraz-Tec et al. 2022).
                    % Convenzione del solutore: g = g0 + W'*q_b, contatto se g <= 0.
                    % Con gap firmato s:  W = -diag(sign(s)),  g0 = |s|.
                    n_bnd = rom.n_bnd;
                    if n_bnd ~= numel(gaps_array)
                        error('MAIN:BndMismatch', ...
                            '%s: n_bnd = %d ma i GdL di contatto sono %d.', ...
                            model, n_bnd, numel(gaps_array));
                    end
                    W  = -diag(sign(gaps_array));
                    g0 = abs(gaps_array);

                    solver = TransientSolverMassless(Mr, Kr, Cr, n_bnd, W, g0);
                    [t, q_rom, lambda, info] = solver.solve( ...
                        cfg.tmax, cfg.dt, q0_r, qd0_r, F_handle);

                    % Il solutore massless integra a passo fisso dt e non
                    % conosce OutputTimes: riportiamo l'uscita sulla griglia
                    % comune. Essendo t_common a passo output_stride*dt, i suoi
                    % istanti sono un sottoinsieme ESATTO della griglia del
                    % solutore: il campionamento non introduce interpolazione.
                    idx   = 1 : cfg.output_stride : numel(t);
                    t     = t(idx);
                    q_rom = q_rom(:, idx);
                    lambda = lambda(:, idx);
                    if numel(t) ~= numel(t_common) || max(abs(t - t_common)) > 1e-12*cfg.tmax
                        warning('MAIN:GridMismatch', ...
                            ['%s: la griglia campionata (%d punti) non coincide con ' ...
                             't_common (%d punti). Il post-processing dovra'' interpolare.'], ...
                            model, numel(t), numel(t_common));
                    end
                else
                    % Contatto a penalita' con ode15s.
                    if strcmp(model, 'Rubin')
                        % Rubin scala la base: gap e penalita' vanno riportati
                        % nelle coordinate scalate.
                        [gaps_run, k_run] = rom.contact_params(gaps_array, k_contact);
                    else
                        gaps_run = gaps_array;
                        k_run    = k_contact;
                    end

                    if any(strcmp(model, {'MT', 'MC'}))
                        % Penalita' proiettata: i GdL di contatto restano fisici
                        % e il solutore li raggiunge tramite Pc.
                        solver_args = {'ContactTargetDOF', contact_dofs, ...
                                       'ModelType', 'MC', 'ProjectionMatrix', Pc};
                    else
                        % ROM CMS: l'interfaccia e' in testa al vettore ridotto.
                        solver_args = {'ContactTargetDOF', 1:numel(contact_dofs), ...
                                       'ModelType', model};
                    end

                    solver = TransientSolverOde(Mr, Kr, Cr);
                    [t, q_rom] = solver.solve(cfg.tmax, cfg.dt, q0_r, qd0_r, F_handle, ...
                        solver_args{:}, ...
                        'ContactGap',     gaps_run, ...
                        'ContactPenalty', k_run, ...
                        'Eref',           Eref, ...
                        'RelTol',         cfg.RelTol, ...
                        'OutputTimes',    t_common);
                end
                cpu_time = toc;

                % ---------- ricostruzione e salvataggio ----------
                y_contact = extract_contact_response(Struct, Interfaces, active_labels, Pc * q_rom);

                n_modes = phi;
                save(fullfile(save_dir, ...
                        sprintf('ROM_%s_Phi%03d_Q%04d_K%g.mat', model, phi, Q, k_mult)), ...
                    't', 'y_contact', 'Interfaces', 'cpu_time', 'offline_time', ...
                    'model', 'n_modes', 'Q', 'k_mult', 'lambda', 'info');
            end
        end
    end
end

fprintf('\n=========================================\n');
fprintf('  Benchmark completato.\n  Risultati in: %s\n', save_dir);
fprintf('=========================================\n');
