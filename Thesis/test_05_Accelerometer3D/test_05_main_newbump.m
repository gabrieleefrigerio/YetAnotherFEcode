%% =====================================================================
%  SIMULATION MAIN - 3D accelerometer, FOM vs ROM with unilateral contact
%
%  Everything below the configuration is a call: the machinery lives in
%  Thesis/Src and is shared with the other test cases. To set up a run, edit
%  section 1.
%
%    FeStructure       mesh, material, node sets, boundary conditions, modes
%    build_interfaces  contact operator and interface bookkeeping
%    shock_forcing     separable impulsive forcing and reference energy
%    run_fom           reference full-order runs, with the contact diagnostic
%    run_rom_sweep     every reduction method over the requested sweep
%
%  Available methods:
%    FOM     full model, penalty contact (ode15s)
%    MT      modal truncation, projected penalty contact(ode15s)
%    MC      Milman-Chu, projected penalty contact(ode15s)
%    CB      Hurty/Craig-Bampton fixed-interface CMS, penalty contact(ode15s)
%    Rubin   free-interface CMS, penalty contact(ode15s)
%    MCB     massless Craig-Bampton, exact set-valued contact (LCP + leapfrog)
%    MN      massless MacNeal,       exact set-valued contact (LCP + leapfrog)
%
%  CB and Rubin keep the interface as physical coordinates at the head of the
%  reduced vector, so they can take a secondary interface reduction.
% =====================================================================
clear; close all; clc;

%% --- 0. PATHS ---------------------------------------------------------
test_root = fileparts(mfilename('fullpath'));
addpath(fullfile(test_root, '..', 'Src'));
addpath(fullfile(test_root, 'mesh'));

% YaFEc startup.
if isempty(which('KirchoffMaterial'))
    run(fullfile(test_root, '..', '..', 'startup.m'));
end

cd(test_root);

%% --- 1. CONFIGURATION -------------------------------------------------

% --- Mesh ---
cfg.mesh_file    = 'Accel3D_V1.mat';   % see mesh/README.md
cfg.element_type = 'WED15';

% --- Material ---
% UNITS: TDK mesh works in the consistent system um / MPa / kg.
cfg.E   = 168e3;      % [MPa]     (168 GPa)
cfg.nu  = 0.23;
cfg.rho = 2.33e-15;   % [kg/um^3] (2330 kg/m^3)

% --- Node sets, by geometric predicate ---
% Each entry is a box [n_dim x 2] of [min max], with -inf/inf on the free directions; 
% Check the Node sets with describe_node_sets() and plot_node_sets().
% The geometry, read off the mesh: the proof mass is a plate from z = 0 to
% z = 30.3. The upper stoppers act on the z = 30.3 face and the lower ones on z = 0.

z_top = 30.3;              % upper face of the proof mass
z_bot = 0.0;               % lower face
z_anc = [52.3, -2.18];     % the two anchoring levels
tol   = 1e-3;              % tolerance on the coordinates

tab = { 'A', [-202.875 -172.875], [-115.000 -105.000]
        'B', [ 172.875  202.875], [-115.000 -105.000]
        'C', [-151.350 -121.350], [ 139.825  149.825]
        'D', [ 162.875  192.875], [ 139.825  149.825] };

box = @(x, y, z) [x(1)-tol x(2)+tol; y(1)-tol y(2)+tol; z-tol z+tol];
for it = 1:size(tab, 1)
    cfg.node_sets.(['stop_' tab{it,1} '_top']) = ...
        struct('box', box(tab{it,2}, tab{it,3}, z_top));
    cfg.node_sets.(['stop_' tab{it,1} '_bot']) = ...
        struct('box', box(tab{it,2}, tab{it,3}, z_bot));
end

cfg.node_sets.anchor = [struct('box', [-inf inf; -inf inf; z_anc(1)-tol z_anc(1)+tol]), ...
                        struct('box', [-inf inf; -inf inf; z_anc(2)-tol z_anc(2)+tol])];
cfg.bc_sets = {'anchor'};

% --- Contact interfaces ---
% For a wall that is NOT parallel to the contacting surface, replace 'gap' with 
%     cfg.interfaces(1).plane = struct('point', [x0 y0 z0]);
% and the gap is computed node by node, varying along the surface exactly as
% the two planes diverge. An oblique normal also needs dofs = 'all', because
% the contact force then has components on all three DOFs of the node.

gap_top = 4.0;            % [um]
gap_bot = 2.0;            % [um]

k = 0;
for it = 1:size(tab, 1)
    for face = {'top', [0 0 1], gap_top; 'bot', [0 0 -1], gap_bot}'
        k = k + 1;
        cfg.interfaces(k).set    = ['stop_' tab{it,1} '_' face{1}];
        cfg.interfaces(k).normal = face{2};
        cfg.interfaces(k).gap    = face{3};
        cfg.interfaces(k).dofs   = 'normal';
    end
end

% --- Methods to run ---
cfg.run.FOM   = 1;
cfg.run.MT    = 0;
cfg.run.MC    = 0;
cfg.run.CB    = 0;
cfg.run.Rubin = 0;
cfg.run.MCB   = 0;  
cfg.run.MN    = 0;

% --- Interface reduction (CB and Rubin only) ---
cfg.interface_reduction.enabled = 1;
cfg.interface_reduction.mode    = 'per_interface';   % 'global' | 'per_interface'
% equal_per_face applies only to mode 'per_interface': when true, every contact
% face gets the SAME number of CC modes. When false, the modes are pooled
% across faces and chosen by frequency.
cfg.interface_reduction.equal_per_face = true;
cfg.interface_reduction.basis   = 'guyan';    % 'guyan'  | 'self'
% static_correction = put back what the truncation throws away. 
% It applies to the penalty methods with interface reduction (CB, Rubin) only. 
cfg.interface_reduction.static_correction = true;
cfg.array_ccModes  = [16, 32, 64, 128, 232]; % 232 = n_bnd, the control point

% --- Sweep ---
cfg.array_linModes = [200];   % 90 is what the previous thesis retained

% --- Damping ---
% Rayleigh, C = alpha*M + beta*K, so 1/Q(f) = alpha/(2*pi*f) + 2*pi*f*beta.
cfg.Q_freq         = [7e3, 2e6];    % anchor frequencies [Hz]
cfg.array_QFactor  = [5,  200];   % Q at 30 kHz, Q at 5 MHz
% cfg.array_QFactor  = 1000;
cfg.array_k_mult   = 0.003125;  % contact stiffness as a multiple of max(diag(K)).
                              % max(diag(K)) = 1.6008e9 N/m on this model, so
                              % 0.3125 reproduces the 5e8 N/m of the previous
                              % thesis. Note this is NOT the value the 2D
                              % benchmark used (10): the same multiplier here
                              % would give a contact 32 times stiffer.

% --- Impulsive forcing ---
cfg.impulse_g   = 1e6;        % amplitude [g]
cfg.impulse_dir = [0 0 -1];    % direction in space, normalized afterwards
cfg.g_value     = 9.81e6;     % gravity in um/s^2: the mesh is in um, not m
cfg.t_shock     = 10e-7;      % half-sine duration [s]

% --- Integration ---
cfg.dt        = 2e-9;
cfg.tmax      = 30e-6;
cfg.RelTol    = 1e-8;         % ROM
cfg.RelTolFOM = 1e-8;         % FOM
cfg.output_stride = 10;       % output every N steps of dt
cfg.integrator    = 'newmark';       % 'ode15s' | 'genalpha' | 'newmark'
cfg.h_int         = 1e-9;           % fixed step of the second-order schemes
cfg.rho_inf       = 0.7;            % genalpha only
cfg.newmark_alpha = 0;              % newmark only
% Store the full displacement field of the FOM.
cfg.save_full_field = true;

cfg.test_name = sprintf('Shock3D_%dg_%.1es', cfg.impulse_g, cfg.tmax);

%% --- 2. RESULTS DIRECTORY ---------------------------------------------
run_dir = fullfile('results', sprintf('%s_%s', cfg.test_name, ...
    char(datetime('now', 'Format', 'yyyy-MM-dd_HH-mm'))));
if ~exist(run_dir, 'dir'), mkdir(run_dir); end
fprintf('Results directory: %s\n\n', run_dir);

%% --- 3. MODEL ---------------------------------------------------------
Struct = FeStructure();
Struct.mesh_file    = cfg.mesh_file;
Struct.element_type = cfg.element_type;
Struct.E            = cfg.E;
Struct.nu           = cfg.nu;
Struct.rho          = cfg.rho;
Struct.set_specs    = cfg.node_sets;
Struct.bc_sets      = cfg.bc_sets;
Struct.build();
Struct.describe_node_sets();

if Struct.n_dim ~= 3
    error('MAIN:Not3D', ...
        'This main is for a 3D model, but the mesh has %d dimensions.', Struct.n_dim);
end

Struct.compute_eigenmodes(max(cfg.array_linModes));

% Reference check against the previous thesis.
fprintf('\nFirst frequencies [kHz]: ');
fprintf('%.2f ', Struct.frequencies(1:min(5, numel(Struct.frequencies)))/1e3);
fprintf('\n(previous thesis: 7.08 7.83 9.41 111.00 119.00)\n\n');

%% --- 4. CONTACT AND FORCING -------------------------------------------
contact = build_interfaces(Struct, cfg);
shock   = shock_forcing(Struct, cfg);

check_cc_modes(cfg, contact);

%% --- 5. RUN -----------------------------------------------------------
save_run_config(run_dir, cfg, Struct, contact, shock);

if cfg.run.FOM
    run_fom(Struct, contact, shock, cfg, run_dir);
end
run_rom_sweep(Struct, contact, shock, cfg, run_dir);

fprintf('\n=========================================\n');
fprintf('  Benchmark complete.\n  Results in: %s\n', run_dir);
fprintf('=========================================\n');

% -----------------------------------------------------------------------
% NOTE on cost. The Accel3D model has 22428 DOFs against the 6484 of the 2D one.
