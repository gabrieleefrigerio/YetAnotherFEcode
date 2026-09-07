%% =====================================================================
%  SHOW SHAPES - look at the vectors any of the models is actually built on
%
%  Three families, all drawn by the same routine because they are all just
%  columns of a basis on the model's DOFs:
%
%    'eigen'       free vibration modes of the structure
%    'rom'         the basis of a ROM: constraint modes and fixed-interface
%                  modes for CB, attachment and normal modes for Rubin, modes
%                  plus static contact vectors for MC
%    'cc'          the interface (characteristic constraint) modes used by the
%                  interface reduction, extended to the whole body through the
%                  constraint modes - pick is the field they actually stand for
%
%  Set the block below and run. Nothing is written to disk.
% =====================================================================
clear; close all; clc;

%% --- CONFIGURATION ----------------------------------------------------
model  = '3D';            % '2D' | '3D'
family = 'cc';            % 'eigen' | 'rom' | 'cc'
pick   = 129:134;              % which shapes: indices into the chosen family
rom_kind = 'Rubin';          % for family 'rom': 'CB' | 'Rubin' | 'MC'
cc_mode  = 'global';          % for family 'cc': 'global' | 'per_interface'
n_modes  = 20;            % modes the ROM is built with (family 'rom' and 'cc')
scale    = [];            % [] = 20% of the model size

%% --- PATHS ------------------------------------------------------------
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'Src'));
run(fullfile(fileparts(here), 'startup.m'));

switch upper(model)
    case '2D'
        cfg_dir = fullfile(here, 'test_04_Dummy_model');
        cfg_run = 'results/Shock_100000g_5.0e-04s_2026-08-14_12-20/run_config.mat';
    case '3D'
        cfg_dir = fullfile(here, 'test_05_Accelerometer3D');
        cfg_run = 'results/Shock3D_1000000g_1.0e-04s_2026-09-01_12-22/run_config.mat';
    otherwise
        error('SHOW:BadModel', 'model must be ''2D'' or ''3D''.');
end
addpath(fullfile(cfg_dir, 'mesh'));

%% --- MODEL ------------------------------------------------------------
% The configuration is taken from a run that already exists, so the mesh, the
% material and above all the boundary conditions match the simulations. Typing
% them again here would be a second source of truth waiting to drift.
R   = load(fullfile(cfg_dir, cfg_run));
cfg = R.cfg;

Struct = FeStructure();
Struct.mesh_file    = cfg.mesh_file;
Struct.element_type = cfg.element_type;
if isfield(cfg, 'E'),  Struct.E   = cfg.E;   end
if isfield(cfg, 'nu'), Struct.nu  = cfg.nu;  end
if isfield(cfg, 'rho'), Struct.rho = cfg.rho; end
if isfield(cfg, 'node_sets'), Struct.set_specs = cfg.node_sets; end
if isfield(cfg, 'bc_sets'),   Struct.bc_sets   = cfg.bc_sets;   end
Struct.build();

%% --- SHAPES -----------------------------------------------------------
switch lower(family)

    case 'eigen'
        Struct.compute_eigenmodes(max(pick));
        V   = Struct.mode_shapes(:, pick);
        lab = arrayfun(@(i) sprintf('mode %d  -  %.4g Hz', i, Struct.frequencies(i)), ...
                       pick, 'Uni', 0);
        ttl = sprintf('%s: free vibration modes', model);

    case {'rom', 'cc'}
        Struct.compute_eigenmodes(n_modes);
        contact = build_interfaces(Struct, cfg);
        k_c = contact.k_base * cfg.array_k_mult(1);
        [~, al, be] = Struct.compute_rayleigh_damping(1000, 1000);
        rom = build_rom(Struct, rom_kind, n_modes, contact, k_c, ...
                        struct('alpha', al, 'beta', be));

        if strcmpi(family, 'rom')
            V   = rom.P(:, pick);
            nb  = numel(contact.dofs);
            lab = arrayfun(@(i) tag_column(i, nb, n_modes, rom_kind), pick, 'Uni', 0);
            ttl = sprintf('%s: %s basis (%d interface + %d modal)', ...
                          model, rom_kind, nb, n_modes);
        else
            % CC modes: the interface eigenproblem of the Guyan pencil, then
            % pushed out into the body by the PHYSICAL constraint modes Psi -
            % never by rom.P(:,1:n_bnd). That column block is what the ROM
            % happens to use as its own boundary parameterization, and for CB
            % it IS Psi, but for Rubin it is the residual attachment modes,
            % which carry a completely different, force-response scaling. Psi is
            % independent of the ROM, which is why guyan_interface_pencil
            % returns it separately (see its docstring).
            %
            % cc_modes is the SAME routine interface_reduction uses, so what is
            % drawn here is exactly the basis the reduction runs on, for either
            % variant:
            %   'global'         one eigenproblem on the whole interface
            %   'per_interface'  one per contact face, modes localized on a face
            [Kbb, Mbb, Psi] = guyan_interface_pencil(Struct, contact.dofs);
            n_show = max(pick);
            [Phi_CC, w2] = cc_modes(Kbb, Mbb, n_show, cc_mode, contact.blocks);
            V   = Struct.AssemblyObj.unconstrain_vector(Psi * Phi_CC(:, pick));
            lab = arrayfun(@(i) sprintf('CC %d  -  %.4g Hz', i, sqrt(max(w2(i),0))/(2*pi)), ...
                           pick, 'Uni', 0);
            ttl = sprintf('%s: CC modes (%s), %d contact DOFs', model, cc_mode, contact.n_bnd);
        end

    otherwise
        error('SHOW:BadFamily', 'family must be ''eigen'', ''rom'' or ''cc''.');
end

plot_shapes(Struct, V, 'Labels', lab, 'Scale', scale, 'Title', ttl);
fprintf('\nDrawn %d shapes: %s\n', numel(pick), ttl);

%% =====================================================================
function s = tag_column(i, n_bnd, n_modes, kind)
% Which half of the basis a column belongs to - the thing one actually wants
% to know when looking at a CMS basis. The two families order it differently:
% the CMS methods put the STATIC part first, Milman-Chu puts the MODES first,
% so the split point and its meaning both depend on the method.
switch upper(kind)
    case {'CB','MCB'}
        if i <= n_bnd, s = sprintf('constraint mode %d', i);
        else,          s = sprintf('fixed-interface mode %d', i - n_bnd); end
    case {'RUBIN','MN'}
        if i <= n_bnd, s = sprintf('attachment mode %d', i);
        else,          s = sprintf('normal mode %d', i - n_bnd); end
    case 'MC'
        if i <= n_modes, s = sprintf('normal mode %d', i);
        else,            s = sprintf('static contact vector %d', i - n_modes); end
    case 'MT'
        s = sprintf('normal mode %d', i);
    otherwise
        s = sprintf('column %d', i);
end
end
