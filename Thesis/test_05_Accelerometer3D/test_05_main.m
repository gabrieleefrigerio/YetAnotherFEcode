%% =====================================================================
%  SIMULATION MAIN - 3D accelerometer, FOM vs ROM with unilateral contact
%
%  Everything below the configuration is a call: the machinery lives in
%  Thesis/Src and is shared with the other test cases. To set up a run, edit
%  section 1 and nothing else.
%
%    FeStructure       mesh, material, node sets, boundary conditions, modes
%    build_interfaces  contact operator and interface bookkeeping
%    shock_forcing     separable impulsive forcing and reference energy
%    run_fom           reference full-order runs, with the contact diagnostic
%    run_rom_sweep     every reduction method over the requested sweep
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
%  reduced vector, so they can take a secondary interface reduction.
%
%  Reference values from the previous thesis, to be reproduced before trusting
%  anything else: f1..f5 = 7.08, 7.83, 9.41, 111, 119 kHz.
%
%  FIRST RUN THE FOM AND READ THE CONTACT ACTIVITY. If the contact is on a
%  corner rather than flat, an n_cc sweep only measures how fast the model
%  degrades. On the 2D benchmark the same mesh gave 81 active nodes out of 81
%  at one shock direction and 4 out of 81 at another, and the interface
%  reduction degraded 11800 times more in the second case.
% =====================================================================
clear; close all; clc;

%% --- 0. PATHS ---------------------------------------------------------
test_root = fileparts(mfilename('fullpath'));
addpath(fullfile(test_root, '..', 'Src'));
addpath(fullfile(test_root, 'mesh'));

% YaFEc itself, if the session has not been initialised yet. Without this the
% run only works after startup.m has been called by hand, which is fine in an
% interactive session and fails in a batch one (matlab -batch starts cold).
if isempty(which('KirchoffMaterial'))
    run(fullfile(test_root, '..', '..', 'startup.m'));
end

cd(test_root);

%% --- 1. CONFIGURATION -------------------------------------------------

% --- Mesh ---
cfg.mesh_file    = 'Accel3D_V1.mat';   % see mesh/README.md
cfg.element_type = 'WED15';

% --- Material ---
% UNITS: the TDK mesh is NOT in metres. Accelerometer_3D.m works in the
% consistent system um / MPa / kg, that is lengths in um, Young's modulus in
% MPa and density in kg/um^3. In that system the assembled K comes out in N/m
% and M in kg, so frequencies are in Hz and times in seconds, while
% DISPLACEMENTS are in um and forces in uN. Putting SI constants on a mesh
% expressed in um would give meaningless frequencies with nothing to flag it.
cfg.E   = 168e3;      % [MPa]     (168 GPa)
cfg.nu  = 0.23;
cfg.rho = 2.33e-15;   % [kg/um^3] (2330 kg/m^3)

% --- Node sets, by geometric predicate ---
% A .mat mesh carries no node sets, so they are declared here. Each entry is a
% box [n_dim x 2] of [min max], with -inf/inf on the free directions; a struct
% array of boxes is their union, which is what disjoint anchor pads need.
% The values come from Accelerometer_3D.m, so they are the ones that produced
% the previous thesis. Check them anyway with describe_node_sets() and
% plot_node_sets(): 'stop' must be a surface and not a thickness, and 'anchor'
% must catch the two levels and nothing in between. All coordinates in um.
% EIGHT stoppers: four tabs, each acting on both faces of the proof mass.
%
% The geometry, read off the mesh: the proof mass is a plate from z = 0 to
% z = 30.3, hanging between a lower post (z from -2.18 to -1.09) and an upper
% one (z from 35.8 to 52.3), both clamped. The upper stoppers act on the
% z = 30.3 face and the lower ones on z = 0.
%
% Accelerometer_3D.m declares one tab only, on the upper face, which is why it
% selects 29 nodes; the previous thesis reports 58, that is two of them. None
% of the four is modelled on the lower face, and without those the proof mass
% falls 139 um instead of the 2 um the device allows.
%
% They are kept as eight separate interfaces rather than one set of 232: with
% mode = 'global' the reduced model is identical either way, but eight blocks
% let contact_activity report each stopper on its own and make
% 'per_interface' available for comparison.
z_top = 30.3;              % upper face of the proof mass
z_bot = 0.0;               % lower face
z_anc = [52.3, -2.18];     % the two anchoring levels
tol   = 1e-3;              % tolerance on the coordinates

% The four tabs, read off the mesh. All of them have the same shape -- rows of
% 7, 4, 7, 4, 7 nodes over 30 um in x and 10 um in y, so 29 nodes each -- and
% each appears identically on both faces of the proof mass, the plan-view node
% pattern of z = 0 and z = 30.3 being the same.
%
% A tab is easy to crop by accident: leave out the last row in y and the count
% drops to 22, leave out a column in x and it drops further. If a set does not
% come out at 29, the box is too tight, not the mesh short of nodes.
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
% One wall above the proof mass, normal along +Z. For a wall that is NOT
% parallel to the contacting surface, replace 'gap' with a plane,
%     cfg.interfaces(1).plane = struct('point', [x0 y0 z0]);
% and the gap is computed node by node, varying along the surface exactly as
% the two planes diverge. An oblique normal also needs dofs = 'all', because
% the contact force then has components on all three DOFs of the node.
% The lower stoppers need no new code: the direction of a wall lives in the
% row of the contact operator, so they are the same declaration with the
% normal reversed and a positive gap.
gap_stop = 2.0;            % travel allowed either way [um]

% One interface per tab per face: eight in total, 232 contact DOFs. The upper
% face is stopped along +Z and the lower one along -Z; nothing else differs,
% and no code changes with the direction.
k = 0;
for it = 1:size(tab, 1)
    for face = {'top', [0 0 1]; 'bot', [0 0 -1]}'
        k = k + 1;
        cfg.interfaces(k).set    = ['stop_' tab{it,1} '_' face{1}];
        cfg.interfaces(k).normal = face{2};
        cfg.interfaces(k).gap    = gap_stop;
        cfg.interfaces(k).dofs   = 'normal';
    end
end

% --- Methods to run ---
cfg.run.FOM   = 1;
cfg.run.MT    = 0;
cfg.run.MC    = 0;
cfg.run.CB    = 0;
cfg.run.Rubin = 0;
cfg.run.MCB   = 0;    % fixed step: see the stability note at the bottom
cfg.run.MN    = 0;

% --- Interface reduction (CB and Rubin only) ---
% On the 2D model 'global' beat 'per_interface' by three to four orders of
% magnitude. With two stoppers the choice exists again, so it is worth one
% comparison run before settling.
cfg.interface_reduction.enabled = 1;
cfg.interface_reduction.mode    = 'per_interface';   % 'global' | 'per_interface'
% equal_per_face applies only to mode 'per_interface': when true, every contact
% face gets the SAME number of CC modes (array_ccModes distributed evenly, so
% k modes per interface means array_ccModes = k * number_of_faces). When false,
% the modes are pooled across faces and chosen by frequency, so a softer face
% takes more of them. The total, hence the reduced size and the file name, is
% array_ccModes either way, so the two are directly comparable in the plots.
cfg.interface_reduction.equal_per_face = true;
cfg.interface_reduction.basis   = 'guyan';    % 'guyan'  | 'self'
% static_correction = put back what the truncation throws away, instead of
% pretending it is not there. The CC modes above n_cc sit at 1e8-1e9 Hz against
% an excitation of about 1 MHz, so they carry no dynamics: they only deflect,
% quasi-statically, under the contact load. That deflection is the LOCAL
% COMPLIANCE of the interface, and a truncated model without it is artificially
% RIGID exactly where the contact is evaluated.
%
% Switching this on restores it in closed form: the contact spring is put in
% series with that compliance. It adds NO state to the ODE - the reduced size is
% unchanged - and only the contact force law changes. It costs 16-32% of CPU and
% is inactive by construction when n_cc = n_bnd, where there is nothing to
% restore.
%
% Measured on the 3D model, Rubin phi = 200 over 10 us (error against the FOM):
%   n_cc =  16    11.34% -> 6.39%
%   n_cc =  64     4.93% -> 1.03%
%   n_cc = 128     1.05% -> 0.204%
%
% Applies to the penalty methods with interface reduction (CB, Rubin). The
% massless models (MCB, MN) solve the contact by a different route and ignore it.
cfg.interface_reduction.static_correction = true;
cfg.array_ccModes  = [16, 32, 64, 128, 232]; % 232 = n_bnd, the control point
cfg.array_ccModes  = [64, 128]; % 232 = n_bnd, the control point
% --- Sweep ---
cfg.array_linModes = [200];   % 90 is what the previous thesis retained
% --- Damping ---
% Rayleigh, C = alpha*M + beta*K, so 1/Q(f) = alpha/(2*pi*f) + 2*pi*f*beta.
% The two anchors fix alpha and beta and EVERYTHING ELSE IS EXTRAPOLATION,
% which is why the anchor FREQUENCIES matter as much as the quality factors.
% They are given explicitly here rather than falling back on the first two
% natural frequencies, which on this model sit at 7.08 and 7.83 kHz, ten per
% cent apart: anchored there, asking for Q = 1000 quietly produces Q = 15 at
% 1 MHz and Q = 3 at 5 MHz, i.e. 17% of critical on exactly the frequencies
% the contact excites.
%
% This is the parametrisation used at TDK. One case only, no sweep: to compare
% two dampings, run this main once per case.
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
cfg.impulse_dir = [0 0 1];    % direction in space, normalized afterwards
cfg.g_value     = 9.81e6;     % gravity in um/s^2: the mesh is in um, not m
cfg.t_shock     = 10e-7;      % half-sine duration [s]

% --- Integration ---
cfg.dt        = 2e-9;
cfg.tmax      = 1e-4;
cfg.RelTol    = 1e-8;         % ROM
cfg.RelTolFOM = 1e-8;         % FOM
% Measured on this model, linear phase, 0.5 us window:
%
%   MaxStep = dt, RelTol 1e-8      236 steps,  36 LU,  132.5 s
%   MaxStep = dt, RelTol 1e-6      141 steps,  27 LU,   96.8 s
%   no MaxStep,   RelTol 1e-6      125 steps,  26 LU,   91.5 s
%   no MaxStep,   scalar AbsTol    502 steps,  70 LU,  287.5 s
%
% Three things follow. MaxStep is NOT what limits the step here: removing it
% saves 11 per cent, so it is kept as cheap insurance for resolving contact
% events, where it will matter. The energy-weighted AbsTol is worth a lot:
% replacing it with a scalar one triples the step count. And the cost is
% almost entirely the factorization, 26 LU at 3.05 s = 87 per cent of the
% total, because a 3D quadratic mesh gives 99.5 nonzeros per row of K and a
% 17.5x fill-in, against 12.0 and 3.7x for the 2D model. That is why the FOM
% costs 60 times more per factorization for only 3.4 times the DOFs.
cfg.output_stride = 10;       % output every N steps of dt

% Store the full displacement field of the FOM, not just the contact nodes.
% Needed by animate_contact_3d to deform the whole mesh; costs a few MB.
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
% NOTE on the massless models (MCB, MN). They integrate at a FIXED step with a
% leapfrog scheme, so they have a stability limit on the step that the penalty
% models do not: on the 2D model, at ten times the nominal step they diverged
% to 1e176 while every other method stayed correct. The 3D model has higher
% frequencies and therefore a tighter limit, so run a step convergence test
% before enabling them here.
%
% NOTE on cost. The TDK model has 22428 DOFs against the 6484 of the 2D one.
% Measure the cost of the FOM on a short window BEFORE setting up any sweep.
% The previous thesis sampled at 283 ns; in 2D the contact events lasted about
% 1 us, that is three or four steps, so that reference may be under-resolved:
% redo a time step convergence on the FOM before using their numbers.
