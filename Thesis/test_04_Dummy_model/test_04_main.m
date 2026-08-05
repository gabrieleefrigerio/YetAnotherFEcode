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
%    CB      Hurty/Craig-Bampton fixed-interface CMS, penalty (ode15s)
%    Rubin   free-interface CMS, penalty (ode15s)
%    MCB     massless Craig-Bampton, exact set-valued contact (LCP + leapfrog)
%    MN      massless MacNeal,       exact set-valued contact (LCP + leapfrog)
%
%  CB and Rubin keep the interface as physical coordinates at the head of the
%  reduced vector, so they can take a secondary interface reduction: see
%  cfg.interface_reduction below and Src/interface_reduction.m.
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
cfg.run.FOM   = 0;
cfg.run.MT    = 0;
cfg.run.MC    = 0;
cfg.run.CB    = 1;
cfg.run.Rubin = 1;
cfg.run.MCB   = 0;
cfg.run.MN    = 0;

% --- Interface reduction (CB and Rubin only) ---
% Secondary modal reduction of the interface partition: the n_bnd physical
% interface DOFs are replaced by a few characteristic constraint (CC) modes.
%   mode = 'global'         one eigenproblem on the coupled boundary partition
%                           (Kuether et al. 2017)
%   mode = 'per_interface'  one eigenproblem per contact face, modes pooled and
%                           sorted by frequency (Aoyama et al. / H-CC flavour).
%                           Localizes the basis per face, which matters when the
%                           contact itself is localized on one face at a time.
% When enabled, the methods eligible for it are run BOTH without reduction and
% at every value of array_ccModes, so the baseline is always available.
% Setting a value equal to n_bnd means no truncation: same model, but through
% the same code path as the reduced ones. Keep it in the sweep as the control
% point that separates the effect of the truncation from that of the change of
% contact evaluation path (direct -> projected).
% basis = where the CC modes are computed from:
%   'guyan' static condensation of the structure onto the contact DOFs. The
%           basis is then a property of the INTERFACE, not of the reduction
%           method. Identical to 'self' for CB, and the only usable choice for
%           Rubin (see below).
%   'self'  boundary partition of the ROM being reduced (Kuether et al. 2017
%           verbatim). For Rubin this is DEGENERATE: its interface block is
%           spanned by residual attachment modes, which carry no low-frequency
%           content, so the retained subspace is nearly orthogonal to the
%           interface motion the dynamics produces and the contact response
%           collapses to zero (measured: 94% of the FOM interface motion
%           unrepresentable at n_cc = 6, GRE = 100%). Kept because reproducing
%           that failure is itself a result worth reporting.
cfg.interface_reduction.enabled = 1;
cfg.interface_reduction.mode    = 'global';
cfg.interface_reduction.basis   = 'guyan';
cfg.array_ccModes               = [6, 14, 26];

% --- Parameter sweeps ---
cfg.array_linModes = [50, 150, 200];
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
%   iface_blocks  index ranges of each face within 1:n_bnd. Since contact_dofs
%                 is concatenated interface by interface these ranges are
%                 contiguous, and they are what interface_reduction() needs to
%                 build a per-face CC basis.
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
iface_blocks = {};
Interfaces   = struct();
nDOFPerNode  = Struct.MeshObj.nDOFPerNode;
offset       = 0;

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

    % Position of this face inside the concatenated interface vector
    block = offset + (1:numel(d));
    iface_blocks{end+1} = block;                                    %#ok<SAGROW>
    offset = offset + numel(d);

    n = Struct.get_contact_nodes(lbl);
    Interfaces.(lbl).rom_idx  = block;
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
n_bnd_total = numel(contact_dofs);
fprintf('\nContact: %d interfaces, %d DOFs in total\n', numel(active_labels), n_bnd_total);

if cfg.interface_reduction.enabled
    bad = cfg.array_ccModes(cfg.array_ccModes < 1 | cfg.array_ccModes > n_bnd_total);
    if ~isempty(bad)
        error('MAIN:BadCcModes', ...
            ['cfg.array_ccModes contains values out of range: %s. ' ...
             'They must lie between 1 and n_bnd = %d.'], mat2str(bad), n_bnd_total);
    end
    fprintf('Interface reduction: %s / %s basis | CC modes %s (n_bnd = %d)\n', ...
        cfg.interface_reduction.mode, cfg.interface_reduction.basis, ...
        mat2str(cfg.array_ccModes), n_bnd_total);
end

% Guyan condensation onto the contact DOFs. Depends only on the structure and
% on which DOFs are in contact, not on the ROM or on how many modes it keeps,
% so it is computed once and reused across the whole sweep.
Kbb_guyan = []; Mbb_guyan = [];
if cfg.interface_reduction.enabled && strcmpi(cfg.interface_reduction.basis, 'guyan')
    fprintf('Computing the Guyan interface pencil (%d constraint modes)...\n', n_bnd_total);
    tic_g = tic;
    [Kbb_guyan, Mbb_guyan] = guyan_interface_pencil(Struct, contact_dofs);
    fprintf('  done in %.1f s\n', toc(tic_g));
end

%% --- 5. FORCING AND INITIAL CONDITIONS --------------------------------
if nDOFPerNode < 2
    error('MAIN:Not2D', 'The model does not have enough DOFs for an in-plane impulse.');
end

impulse_amp = cfg.impulse_g * 9.81;
impulse_dir = cfg.impulse_sign * [cosd(cfg.impulse_angle_deg); sind(cfg.impulse_angle_deg)];

dir_vector = zeros(n_dofs_fom, 1);
dir_vector(1:nDOFPerNode:n_dofs_fom) = impulse_dir(1);   % X DOFs
dir_vector(2:nDOFPerNode:n_dofs_fom) = impulse_dir(2);   % Y DOFs

% The forcing is a fixed spatial vector times a scalar function of time. Keeping
% the two factors separate lets the ROMs project the spatial part ONCE, instead
% of redoing a dense (r x n_dofs) product at every ODE function evaluation.
F_spatial_fom = Mc * dir_vector;
shock_profile = @(t) impulse_amp * sin(pi*t/cfg.t_shock) * (t <= cfg.t_shock);
F_fom_handle  = @(t) F_spatial_fom * shock_profile(t);

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
    'iface_blocks', 'Eref', 't_common', 'n_dofs_fom', 'k_base');

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
rom_list = {'MT', 'MC', 'CB', 'Rubin', 'MCB', 'MN'};
rom_list = rom_list(cellfun(@(m) cfg.run.(m), rom_list) == 1);

% Methods whose basis keeps the interface as physical coordinates at the head
% of the reduced vector, and which can therefore take a secondary interface
% reduction. The massless ones are excluded because their M_bb is zero.
ir_eligible = {'CB', 'Rubin'};

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
                    case 'CB'
                        rom = RomCB(Struct, phi, contact_dofs);
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

                % Spatial part of the forcing projected once. The time
                % dependence is a scalar, so there is no need to redo this
                % product at every ODE function evaluation.
                F_r_spatial_full = Pc' * F_spatial_fom;

                % ---------- interface reduction variants ----------
                % n_cc = 0 is the baseline without reduction, always run so the
                % comparison term is present in every results folder.
                do_ir = cfg.interface_reduction.enabled && any(strcmp(model, ir_eligible));
                if do_ir
                    cc_list = [0, cfg.array_ccModes];
                else
                    cc_list = 0;
                end

                for n_cc = cc_list

                    if n_cc == 0
                        Mr_v = Mr;  Kr_v = Kr;  Cr_v = Cr;
                        Phi_CC = [];  ir_info = [];
                        ir_mode   = 'none';
                        model_tag = model;
                        F_r_spatial = F_r_spatial_full;
                    else
                        tic_ir = tic;

                        % Guyan basis: the pencil is known in PHYSICAL interface
                        % coordinates, while interface_reduction works in the
                        % ROM's own interface coordinates. With u_phys = A*u_rom
                        % (A = the contact rows of the boundary block of Pc, the
                        % identity for CB and the scaling for Rubin), the
                        % congruence A'*(.)*A moves the pencil across and leaves
                        % the eigenvalues unchanged, so Phi_CC comes back
                        % directly in ROM coordinates.
                        if strcmpi(cfg.interface_reduction.basis, 'guyan')
                            A = Pc(contact_dofs, 1:rom.n_bnd);
                            pencil = struct('K', A' * Kbb_guyan * A, ...
                                            'M', A' * Mbb_guyan * A);
                        else
                            pencil = [];
                        end

                        [Mr_v, Kr_v, Cr_v, Phi_CC, ir_info] = interface_reduction( ...
                            Mr, Kr, Cr, rom.n_bnd, n_cc, ...
                            cfg.interface_reduction.mode, iface_blocks, pencil);
                        ir_mode = ir_info.mode;
                        switch ir_mode
                            case 'global',        model_tag = [model 'IRG'];
                            case 'per_interface', model_tag = [model 'IRP'];
                        end

                        % T_CC' applied to the projected forcing, without forming
                        % T_CC: it is blkdiag(Phi_CC, I).
                        F_r_spatial = [Phi_CC' * F_r_spatial_full(1:rom.n_bnd); ...
                                       F_r_spatial_full(rom.n_bnd+1:end)];
                        offline_time = offline_time + toc(tic_ir);
                    end

                    fprintf('  [%s] reduced size %d\n', model_tag, size(Kr_v, 1));

                    % ---------- reduced initial conditions and forcing ----------
                    r_v = size(Kr_v, 1);
                    if any(q0)
                        % Only needed for a non-zero initial state; costs a full
                        % least-squares solve, hence the guard.
                        if n_cc == 0
                            q0_r = Pc \ q0;
                        else
                            q0_r = (Pc * blkdiag(Phi_CC, eye(r_v - n_cc))) \ q0;
                        end
                    else
                        q0_r = zeros(r_v, 1);
                    end
                    qd0_r    = zeros(r_v, 1);
                    F_handle = @(t) F_r_spatial * shock_profile(t);

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

                        solver = TransientSolverMassless(Mr_v, Kr_v, Cr_v, n_bnd, W, g0);
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
                        gaps_run = gaps_array;
                        k_run    = k_contact;

                        if n_cc > 0
                            % Interface reduced: the interface is now modal, so the
                            % contact must be evaluated in physical space through a
                            % projection. That projection is just Phi_CC: for CB and
                            % Rubin the contact rows of Pc are [D^-1, 0], so
                            %   Pc(contact,:) * T_CC = [D^-1 * Phi_CC, 0]
                            % and nothing of full size needs to be formed.
                            Pc_bnd = Pc(contact_dofs, 1:rom.n_bnd);
                            coupling = norm(Pc(contact_dofs, rom.n_bnd+1:end), 'fro');
                            % Rubin leaves round-off here (~1e-11 relative) from
                            % inverting T1_bb, so the threshold is loose enough to
                            % ignore that while still catching a genuinely coupled
                            % modal block, which would be O(1) relative.
                            if coupling > 1e-8 * max(norm(Pc_bnd, 'fro'), eps)
                                error('MAIN:InterfaceNotAtHead', ...
                                    ['%s: the modal block of Pc is not zero at the contact rows ' ...
                                     '(||.|| = %.3e). This ROM does not keep the interface as ' ...
                                     'physical coordinates at the head of the basis, so the ' ...
                                     'interface reduction cannot be applied in this form.'], ...
                                    model, coupling);
                            end
                            Pc_contact_ir = [Pc_bnd * Phi_CC, zeros(numel(contact_dofs), r_v - n_cc)];

                            solver_args = {'ContactTargetDOF', 1:numel(contact_dofs), ...
                                           'ModelType', 'MC', 'ProjectionMatrix', Pc_contact_ir};
                            % Gaps and penalty stay physical: the projected path
                            % works in physical coordinates, so Rubin's scaling must
                            % NOT be applied on top (it is already inside Pc).

                        elseif any(strcmp(model, {'MT', 'MC'}))
                            % Projected penalty: the contact DOFs stay physical and
                            % the solver reaches them through Pc.
                            solver_args = {'ContactTargetDOF', contact_dofs, ...
                                           'ModelType', 'MC', 'ProjectionMatrix', Pc};
                        else
                            % CMS ROM: the interface sits at the head of the
                            % reduced vector, so the contact DOFs are direct indices.
                            if strcmp(model, 'Rubin')
                                % Rubin scales the basis: gap and penalty must be
                                % expressed in the scaled coordinates.
                                [gaps_run, k_run] = rom.contact_params(gaps_array, k_contact);
                            end
                            solver_args = {'ContactTargetDOF', 1:numel(contact_dofs), ...
                                           'ModelType', model};
                        end

                        solver = TransientSolverOde(Mr_v, Kr_v, Cr_v);
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
                    % For the interface-reduced variants, go back through T_CC first
                    % and only then through Pc. Applying blkdiag(Phi_CC, I) to the
                    % time history is far cheaper than forming Pc*T_CC.
                    if n_cc == 0
                        q_rom_full = q_rom;
                    else
                        q_rom_full = [Phi_CC * q_rom(1:n_cc, :); q_rom(n_cc+1:end, :)];
                    end
                    y_contact = extract_contact_response(Struct, Interfaces, active_labels, ...
                                                         Pc * q_rom_full);

                    n_modes = phi;
                    if n_cc == 0
                        file_name = sprintf('ROM_%s_Phi%03d_Q%04d_K%g.mat', ...
                                            model_tag, phi, Q, k_mult);
                    else
                        file_name = sprintf('ROM_%s_Phi%03d_CC%03d_Q%04d_K%g.mat', ...
                                            model_tag, phi, n_cc, Q, k_mult);
                    end
                    save(fullfile(save_dir, file_name), ...
                        't', 'y_contact', 'Interfaces', 'cpu_time', 'offline_time', ...
                        'model', 'model_tag', 'n_modes', 'Q', 'k_mult', 'lambda', 'info', ...
                        'n_cc', 'ir_mode', 'ir_info');

                end   % n_cc
            end   % im (rom_list)
        end
    end
end

fprintf('\n=========================================\n');
fprintf('  Benchmark complete.\n  Results in: %s\n', save_dir);
fprintf('=========================================\n');
