%% =====================================================================
%  SIMULATION MAIN - FOM vs ROM benchmark with unilateral contact
%
%  Single main, adaptive in the number of contact interfaces. The interfaces
%  are read from the .inp file (node sets 'ContactInterface[_<label>]') and
%  declared below in cfg.interfaces with their direction and signed gap.
%  Adding or removing a row of that table is the only change needed to move
%  from one model to another.
%
%  Available methods:
%    FOM     full model, penalty contact (ode15s)
%    MT      modal truncation, projected penalty (ode15s)
%    MC      Milman-Chu, projected penalty (ode15s)
%    Rubin   free-interface CMS, penalty (ode15s)
%    MCB     massless Craig-Bampton, exact set-valued contact (LCP + leapfrog)
%    MN      massless MacNeal,       exact set-valued contact (LCP + leapfrog)
%
%  Results go to results/<test_name>_<timestamp>/, together with
%  run_config.mat holding the full configuration: the post-processing reads
%  it from there and needs no prior knowledge of the model.
% =====================================================================
clear; close all; clc;

%% --- 1. CONFIGURATION -------------------------------------------------

% --- Model ---
cfg.mesh_file    = 'DummyStructureAbaqus_V4.inp';
cfg.element_type = 'TRI3';

% --- Contact interfaces ---
%   label | direction (1 = X, 2 = Y) | signed gap [m]
% The SIGN of the gap tells which side the wall is on:
%   gap > 0  wall along the positive direction of the DOF (penetrates if q > gap)
%   gap < 0  wall along the negative direction of the DOF (penetrates if q < gap)
% The labels must exist in the .inp as *Nset, nset=ContactInterface_<label>.
% For a single-interface model (*Nset, nset=ContactInterface) the label is
% 'C' and the table collapses to a single row.
cfg.interfaces = { ...
    'T', 2,  5.0e-6 ; ...   % wall on the positive Y side
    'B', 2, -1.5e-6 ; ...   % wall on the negative Y side
    'L', 1, -5.0e-6 ; ...   % wall on the negative X side
    'R', 1,  1.5e-6 };      % wall on the positive X side

% --- Methods to run ---
cfg.run.FOM   = 1;
cfg.run.MT    = 1;
cfg.run.MC    = 1;
cfg.run.Rubin = 1;
cfg.run.MCB   = 0;
cfg.run.MN    = 0;

% --- Parameter sweeps ---
cfg.array_linModes = [10, 50, 100, 150, 200];
cfg.array_QFactor  = [1000];
cfg.array_k_mult   = [10];        % contact stiffness multiplier
                                  % (ignored by the massless models MCB/MN)

% --- Impulsive forcing ---
cfg.impulse_g         = 1e5;      % amplitude [g]
cfg.impulse_angle_deg = 0;        % direction in the XY plane [deg]
cfg.impulse_sign      = 1;        % orientation (+1 / -1)
cfg.t_shock           = 10e-7;    % half-sine duration [s]

% --- Integration ---
cfg.dt      = 0.5e-8;
cfg.tmax    = (1e-3)/2;
cfg.RelTol  = 1e-9;               % ode15s, ROM
cfg.RelTolFOM = 1e-10;            % ode15s, FOM (tighter reference)
cfg.output_stride = 10;           % output every N steps of dt (common grid)

% --- Test name ---
cfg.test_name = sprintf('Shock_%dg_%.1es', cfg.impulse_g, cfg.tmax);

%% --- 2. RESULTS DIRECTORY ---------------------------------------------
timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd_HH-mm'));
save_dir  = fullfile('results', sprintf('%s_%s', cfg.test_name, timestamp));
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end
fprintf('Results directory: %s\n\n', save_dir);

%% --- 3. MODEL ---------------------------------------------------------
fprintf('Building the model...\n');
Struct = AbaqusStructure();
Struct.filename    = cfg.mesh_file;
Struct.elementType = cfg.element_type;
Struct.build();
Struct.describe_interfaces();

max_phi = max(cfg.array_linModes);
fprintf('Extracting %d modes...\n', max_phi);
Struct.compute_eigenmodes(max_phi);

Mc = Struct.AssemblyObj.constrain_matrix(Struct.M);
Kc = Struct.AssemblyObj.constrain_matrix(Struct.K);
n_dofs_fom = size(Mc, 1);
k_base     = max(diag(Kc));

%% --- 4. CONTACT INTERFACES --------------------------------------------
% From cfg.interfaces we derive, in one pass:
%   contact_dofs  constrained contact DOFs, concatenated
%   gaps_array    signed gap for each of those DOFs
%   Interfaces    metadata for the post-processing (nodes, global DOFs, coords)
labels     = cfg.interfaces(:, 1)';
dirs       = cell2mat(cfg.interfaces(:, 2))';
gaps_iface = cell2mat(cfg.interfaces(:, 3))';

% Every declared label must exist in the .inp file
missing = setdiff(labels, Struct.contact_labels);
if ~isempty(missing)
    error('MAIN:MissingInterface', ...
        ['Interfaces {%s} do not exist in %s.\n' ...
         'Interfaces available in the file: {%s}'], ...
        strjoin(missing, ', '), cfg.mesh_file, strjoin(Struct.contact_labels, ', '));
end
unused = setdiff(Struct.contact_labels, labels);
if ~isempty(unused)
    fprintf('[note] Interfaces present in the .inp but not used: %s\n', strjoin(unused, ', '));
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
            'Interface %s: no free DOF (all its nodes are constrained). Skipped.', lbl);
        continue;
    end

    contact_dofs = [contact_dofs; d];                               %#ok<AGROW>
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
    error('MAIN:NoContact', 'No active contact DOF: check cfg.interfaces.');
end
fprintf('\nContact: %d interfaces, %d DOFs in total\n', numel(active_labels), numel(contact_dofs));

%% --- 5. FORCING AND INITIAL CONDITIONS --------------------------------
if nDOFPerNode < 2
    error('MAIN:Not2D', 'The model does not have enough DOFs for an in-plane impulse.');
end

impulse_amp = cfg.impulse_g * 9.81;
impulse_dir = cfg.impulse_sign * [cosd(cfg.impulse_angle_deg); sind(cfg.impulse_angle_deg)];

dir_vector = zeros(n_dofs_fom, 1);
dir_vector(1:nDOFPerNode:n_dofs_fom) = impulse_dir(1);   % X DOFs
dir_vector(2:nDOFPerNode:n_dofs_fom) = impulse_dir(2);   % Y DOFs

F_spatial_fom = Mc * dir_vector;
F_fom_handle  = @(t) F_spatial_fom * impulse_amp * sin(pi*t/cfg.t_shock) * (t <= cfg.t_shock);

q0  = zeros(n_dofs_fom, 1);
qd0 = zeros(n_dofs_fom, 1);

% Reference energy for the energy-weighted AbsTol
v_max = impulse_amp * 2 * cfg.t_shock / pi;
m_eff = dir_vector' * Mc * dir_vector;
Eref  = 0.5 * m_eff * v_max^2;

% Output grid shared by every model, so that the post-processing can compare
% the time histories without interpolating.
t_common = 0 : cfg.output_stride*cfg.dt : cfg.tmax;

fprintf('Impulse: %.1e g @ %.1f deg (orientation %+d)\n', ...
    cfg.impulse_g, cfg.impulse_angle_deg, cfg.impulse_sign);
fprintf('Eref = %.4e J   (v_max = %.4f m/s)\n', Eref, v_max);

% Configuration saved once: this is the contract with the post-processing
save(fullfile(save_dir, 'run_config.mat'), ...
    'cfg', 'Interfaces', 'active_labels', 'contact_dofs', 'gaps_array', ...
    'Eref', 't_common', 'n_dofs_fom', 'k_base');

%% --- 6. FOM -----------------------------------------------------------
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

%% --- 7. ROM -----------------------------------------------------------
rom_list = {'MT', 'MC', 'Rubin', 'MCB', 'MN'};
rom_list = rom_list(cellfun(@(m) cfg.run.(m), rom_list) == 1);

if ~isempty(rom_list)
    fprintf('\n=========================================\n');
    fprintf('                 ROM\n');
    fprintf('=========================================\n');
end

for Q = cfg.array_QFactor

    % compute_rayleigh_damping does two things:
    %   (a) updates Struct.C            -> used by the penalty ROMs (MT/MC/Rubin)
    %   (b) returns alpha and beta      -> used by the massless ROMs, which
    %       build an equivalent diagonal modal damping from them
    % Using the same alpha/beta keeps the damping identical across all ROMs.
    [~, alpha_ray, beta_ray] = Struct.compute_rayleigh_damping(Q, Q);
    rayleigh = struct('alpha', alpha_ray, 'beta', beta_ray);

    % Diagnostic: Rayleigh damping overdamps the high modes.
    w_hi = 2*pi * Struct.frequencies(max_phi);
    z_hi = 0.5*(alpha_ray/w_hi + beta_ray*w_hi);
    fprintf('  zeta(mode %d) = %.4f | target zeta (modes 1-2) = %.4f\n', ...
        max_phi, z_hi, 1/(2*Q));
    if z_hi > 1
        warning('MAIN:Overdamped', ...
            ['Rayleigh makes the high modes OVERDAMPED (zeta_%d = %.2f). ' ...
             'This affects ALL ROMs, not just the massless ones.'], max_phi, z_hi);
    end

    for k_mult = cfg.array_k_mult
        k_contact = k_base * k_mult;

        for phi = cfg.array_linModes
            for im = 1:numel(rom_list)
                model = rom_list{im};

                % The massless models use exact set-valued contact: the
                % penalty stiffness does not apply to them, so the loop over
                % k_mult would be degenerate. They run at the first value only.
                is_massless = any(strcmp(model, {'MCB', 'MN'}));
                if is_massless && k_mult ~= cfg.array_k_mult(1)
                    continue;
                end

                fprintf('\n--- %s | Phi %d | Q %d | k_mult %g ---\n', model, phi, Q, k_mult);

                % ---------- build the basis ----------
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

                % Projection matrix on the constrained DOFs.
                % MT and MC expose P on the global DOFs, the others expose Pc
                % already constrained.
                if any(strcmp(model, {'MT', 'MC'}))
                    Pc = zeros(n_dofs_fom, size(rom.P, 2));
                    for ic = 1:size(Pc, 2)
                        Pc(:, ic) = Struct.AssemblyObj.constrain_vector(rom.P(:, ic));
                    end
                else
                    Pc = rom.Pc;
                end
                offline_time = toc(tic_offline);

                % ---------- reduced initial conditions and forcing ----------
                if any(q0)
                    q0_r = Pc \ q0;
                else
                    q0_r = zeros(size(Pc, 2), 1);   % q0 is zero: no projection needed
                end
                qd0_r    = zeros(size(Pc, 2), 1);
                F_handle = @(t) Pc' * F_fom_handle(t);

                % ---------- integration ----------
                lambda = []; info = [];
                tic;
                if is_massless
                    % Exact set-valued contact (Monjaraz-Tec et al. 2022).
                    % Solver convention: g = g0 + W'*q_b, contact when g <= 0.
                    % With a signed gap s:  W = -diag(sign(s)),  g0 = |s|.
                    n_bnd = rom.n_bnd;
                    if n_bnd ~= numel(gaps_array)
                        error('MAIN:BndMismatch', ...
                            '%s: n_bnd = %d but there are %d contact DOFs.', ...
                            model, n_bnd, numel(gaps_array));
                    end
                    W  = -diag(sign(gaps_array));
                    g0 = abs(gaps_array);

                    solver = TransientSolverMassless(Mr, Kr, Cr, n_bnd, W, g0);
                    [t, q_rom, lambda, info] = solver.solve( ...
                        cfg.tmax, cfg.dt, q0_r, qd0_r, F_handle);

                    % The massless solver integrates at fixed dt and knows
                    % nothing about OutputTimes, so its output is brought back
                    % onto the common grid. Since t_common has a step of
                    % output_stride*dt, its instants are an EXACT subset of the
                    % solver grid: sampling introduces no interpolation.
                    idx    = 1 : cfg.output_stride : numel(t);
                    t      = t(idx);
                    q_rom  = q_rom(:, idx);
                    lambda = lambda(:, idx);
                    if numel(t) ~= numel(t_common) || max(abs(t - t_common)) > 1e-12*cfg.tmax
                        warning('MAIN:GridMismatch', ...
                            ['%s: the sampled grid (%d points) does not match t_common ' ...
                             '(%d points). The post-processing will have to interpolate.'], ...
                            model, numel(t), numel(t_common));
                    end
                else
                    % Penalty contact with ode15s.
                    if strcmp(model, 'Rubin')
                        % Rubin scales the basis: gap and penalty must be
                        % expressed in the scaled coordinates.
                        [gaps_run, k_run] = rom.contact_params(gaps_array, k_contact);
                    else
                        gaps_run = gaps_array;
                        k_run    = k_contact;
                    end

                    if any(strcmp(model, {'MT', 'MC'}))
                        % Projected penalty: the contact DOFs stay physical and
                        % the solver reaches them through Pc.
                        solver_args = {'ContactTargetDOF', contact_dofs, ...
                                       'ModelType', 'MC', 'ProjectionMatrix', Pc};
                    else
                        % CMS ROM: the interface sits at the head of the
                        % reduced vector.
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

                % ---------- reconstruction and saving ----------
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
fprintf('  Benchmark complete.\n  Results in: %s\n', save_dir);
fprintf('=========================================\n');
