function [contact, Interfaces] = build_interfaces(Struct, cfg)
%BUILD_INTERFACES Contact operator and interface bookkeeping for one run.
%
%   [contact, Interfaces] = BUILD_INTERFACES(Struct, cfg)
%
% Turns cfg.interfaces into everything the rest of a run needs, so that a main
% does not have to assemble it by hand.
%
% OUTPUT
%   contact.N        [n_c x n_free] contact operator; penetration is N*u - gaps
%   contact.gaps     [n_c x 1] positive gaps
%   contact.dofs     interface partition, constrained numbering
%   contact.blocks   {1 x n_iface} index ranges of each face within 1:n_bnd,
%                    contiguous by construction, which is what
%                    interface_reduction needs for a per-face CC basis
%   contact.labels   {1 x n_iface} interface names
%   contact.n_bnd    total number of interface DOFs
%   contact.k_base   max(diag(Kc)), the reference for the contact penalty
%   contact.info     the raw output of contact_operator
%
%   Interfaces       per-interface metadata for the post-processing: node IDs,
%                    global DOF indices, nodal coordinates, wall normal and
%                    per-node gaps
%
% See also CONTACT_OPERATOR, INTERFACE_REDUCTION, EXTRACT_CONTACT_RESPONSE.

    [N, gaps, cinfo] = contact_operator(Struct, cfg.interfaces);

    Kc = Struct.AssemblyObj.constrain_matrix(Struct.K);
    nD = Struct.MeshObj.nDOFPerNode;
    axis_name = {'X', 'Y', 'Z'};

    contact = struct();
    contact.N      = N;
    contact.gaps   = gaps;
    contact.dofs   = cinfo.bnd_dofs;
    contact.labels = cinfo.labels;
    contact.n_bnd  = numel(cinfo.bnd_dofs);
    contact.k_base = full(max(diag(Kc)));
    contact.info   = cinfo;
    contact.blocks = cell(1, numel(cinfo.blocks));

    Interfaces = struct();
    offset = 0;
    for i = 1:numel(cinfo.labels)
        lbl   = cinfo.labels{i};
        block = offset + (1:cinfo.blocks(i));
        contact.blocks{i} = block;
        offset = offset + cinfo.blocks(i);

        ids = cinfo.nodes{i};
        Interfaces.(lbl).rom_idx   = block;
        Interfaces.(lbl).nodes     = ids;
        Interfaces.(lbl).normal    = cinfo.spec(i).normal;
        % Per-node gaps, which differ from one another as soon as the wall is
        % not parallel to the surface. The scalar is kept for the plots that
        % expect a single value.
        Interfaces.(lbl).gap_nodes = gaps(block);
        Interfaces.(lbl).gap       = gaps(block(1));
        for d = 1:nD
            Interfaces.(lbl).(['global_' axis_name{d}]) = (ids - 1)*nD + d;
            Interfaces.(lbl).(['coord_'  axis_name{d}]) = Struct.nodes(ids, d);
        end
    end

    contact.Interfaces = Interfaces;

    fprintf('Contact: %d interfaces, %d DOFs in total\n', ...
        numel(cinfo.labels), contact.n_bnd);
end
