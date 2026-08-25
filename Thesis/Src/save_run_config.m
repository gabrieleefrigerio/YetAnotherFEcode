function save_run_config(run_dir, cfg, Struct, contact, shock)
%SAVE_RUN_CONFIG Write run_config.mat, the contract with the post-processing.
%
%   SAVE_RUN_CONFIG(run_dir, cfg, Struct, contact, shock)
%
% Saved once at the start of a run, so that the post-processing needs no prior
% knowledge of the model: it reads the geometry of the interfaces, the contact
% operator, the output grid and the whole configuration from here.
%
% The variable names are flat rather than the structs used inside the run,
% because they are what the existing post-processing expects.

    Interfaces    = contact.Interfaces;   %#ok<NASGU>
    active_labels = contact.labels;       %#ok<NASGU>
    contact_dofs  = contact.dofs;         %#ok<NASGU>
    gaps_array    = contact.gaps;         %#ok<NASGU>
    iface_blocks  = contact.blocks;       %#ok<NASGU>
    N_fom         = contact.N;            %#ok<NASGU>
    cinfo         = contact.info;         %#ok<NASGU>
    k_base        = contact.k_base;       %#ok<NASGU>
    Eref          = shock.Eref;           %#ok<NASGU>
    t_common      = shock.t_out;          %#ok<NASGU>
    n_dofs_fom    = size(Struct.AssemblyObj.constrain_matrix(Struct.K), 1); %#ok<NASGU>
    frequencies   = Struct.frequencies;   %#ok<NASGU>

    save(fullfile(run_dir, 'run_config.mat'), ...
        'cfg', 'Interfaces', 'active_labels', 'contact_dofs', 'gaps_array', ...
        'iface_blocks', 'Eref', 't_common', 'n_dofs_fom', 'k_base', ...
        'N_fom', 'cinfo', 'frequencies');
end
