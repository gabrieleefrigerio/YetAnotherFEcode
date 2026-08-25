function [N, g, info] = contact_operator(Struct, ifaces)
%CONTACT_OPERATOR Build the linear contact operator of a set of rigid walls.
%
%   [N, g, info] = CONTACT_OPERATOR(Struct, ifaces)
%
% Every contact node contributes one unilateral constraint, written as one row
% of N. With u the displacement vector in the constrained numbering,
%
%       p = N*u - g                                              (penetration)
%
% is positive where the node has crossed its wall, and the penalty force and
% its Jacobian follow from the same operator:
%
%       f = N(a,:)' * (k(a) .* p(a))        a = active rows, p(a) > 0
%       J = N(a,:)' * diag(k(a)) * N(a,:)
%
% Because the direction of the wall is carried by the row of N, the gap is a
% plain distance and is always positive. This replaces the earlier convention,
% where the gap was signed and its sign told which side the wall was on: that
% convention could only describe walls normal to a coordinate axis, and it
% needed the activation test sign(gap)*(u-gap) > 0 to undo the sign.
%
% Reduced models use the same operator. If u = Pc*q maps reduced coordinates to
% physical DOFs, then N*Pc is the operator in reduced coordinates and every
% formula above holds unchanged.
%
% INTERFACE SPECIFICATION
%   ifaces is a struct array, one entry per wall, with fields:
%
%     set      name of the node set of the contacting surface, as known to
%              Struct (from the mesh file or from a geometric selector).
%     normal   [1 x n_dim] direction along which penetration is measured,
%              pointing from the body towards the wall. Normalized internally.
%     gap      distance from the undeformed node to the wall, measured along
%              normal. Scalar (same for every node) or [n_nodes x 1].
%              Use this when the wall is parallel to the contacting surface.
%     plane    alternative to gap, for a wall that is NOT parallel to the
%              surface: struct('point', p0) with p0 a point of the wall plane.
%              The gap is then computed node by node,
%                   g_i = normal * (p0 - x_i)',
%              so it varies along the surface exactly as the two planes
%              diverge. This is the general case; a scalar gap is the
%              particular case of parallel planes.
%     dofs     'normal' (default) keeps one DOF per node, the one along the
%              normal, and requires the normal to be aligned with an axis.
%              'all' keeps every DOF of the node, which is what an oblique
%              normal needs, since the contact force then has components on
%              all of them. This choice only affects the interface partition
%              returned in info.bnd_dofs, not N itself.
%
%   A legacy cell array {label, direction, signed_gap; ...} is also accepted
%   and converted, so that configurations written for the signed-gap
%   convention keep working.
%
% OUTPUT
%   N      [n_c x n_free] sparse contact operator, one row per contact node
%   g      [n_c x 1] positive gaps
%   info   .nodes      {1 x n_iface} node IDs actually retained per interface
%          .bnd_dofs   interface partition, DOFs in the constrained numbering
%          .blocks     [1 x n_iface] size of each block of bnd_dofs
%          .rows       [1 x n_iface] rows of N belonging to each interface
%          .labels     {1 x n_iface} set names
%          .spec       the normalized interface specification
%
% See also FESTRUCTURE, TRANSIENTSOLVERODE, INTERFACE_REDUCTION.

    ifaces = normalize_interface_spec(ifaces);
    n_if   = numel(ifaces);
    nD     = Struct.MeshObj.nDOFPerNode;
    n_free = size(Struct.AssemblyObj.constrain_matrix(Struct.K), 1);

    rows = []; cols = []; vals = [];
    g = [];
    info = struct('nodes', {cell(1,n_if)}, 'bnd_dofs', [], ...
                  'blocks', zeros(1,n_if), 'rows', {cell(1,n_if)}, ...
                  'labels', {cell(1,n_if)}, 'spec', ifaces);
    bnd = [];
    r0 = 0;

    for i = 1:n_if
        s = ifaces(i);
        nrm = s.normal(:)' / norm(s.normal);
        ids = Struct.get_nodes(s.set);
        x   = Struct.nodes(ids, :);

        % --- gaps, uniform or from a wall plane --------------------------
        if ~isempty(s.plane)
            gi = (s.plane.point(:)' - x) * nrm(:);
        else
            gi = s.gap(:);
            if isscalar(gi), gi = gi * ones(numel(ids), 1); end
            if numel(gi) ~= numel(ids)
                error('contact_operator:GapSize', ...
                    ['Interface ''%s'': gap has %d entries but the set has %d ' ...
                     'nodes.'], s.set, numel(gi), numel(ids));
            end
        end
        if any(gi <= 0)
            error('contact_operator:NodeInsideWall', ...
                ['Interface ''%s'': %d nodes have a non-positive gap, i.e. they ' ...
                 'start on or beyond the wall. Check the sign of the normal: it ' ...
                 'must point from the body towards the wall.'], ...
                s.set, nnz(gi <= 0));
        end

        % --- which DOFs of a node take part -------------------------------
        [~, ax] = max(abs(nrm));
        aligned = abs(abs(nrm(ax)) - 1) < 1e-10;
        switch lower(s.dofs)
            case 'normal'
                if ~aligned
                    error('contact_operator:ObliqueNormal', ...
                        ['Interface ''%s'': dofs = ''normal'' keeps only the DOF ' ...
                         'along the normal, which requires the normal to be ' ...
                         'aligned with an axis. This one is [%s]. Use ' ...
                         'dofs = ''all''.'], s.set, num2str(nrm, ' %.4g'));
                end
                keep = ax;
            case 'all'
                keep = 1:nD;
            otherwise
                error('contact_operator:BadDofs', ...
                    'Interface ''%s'': dofs must be ''normal'' or ''all''.', s.set);
        end

        % --- assemble the rows, skipping clamped nodes --------------------
        kept = false(numel(ids), 1);
        n_blk0 = numel(bnd);
        kk = keep(:);
        nk = nrm(kk); nk = nk(:);
        for j = 1:numel(ids)
            cd = Struct.AssemblyObj.free2constrained_index((ids(j)-1)*nD + kk');
            cd = cd(:);

            % A row is created only for the DOFs that are free AND carry a
            % component of the normal. A node whose normal-bearing DOFs are
            % all clamped can never come into contact.
            use = cd > 0 & nk ~= 0;
            if ~any(use), continue; end
            kept(j) = true;
            r0 = r0 + 1;
            rows = [rows; repmat(r0, nnz(use), 1)];  %#ok<AGROW>
            cols = [cols; cd(use)];                  %#ok<AGROW>
            vals = [vals; nk(use)];                  %#ok<AGROW>

            % The interface partition keeps every free DOF selected by dofs,
            % including those with a zero normal component: the contact force
            % is transmitted to them through the stiffness of the substructure.
            bnd = [bnd; cd(cd > 0)];                 %#ok<AGROW>
        end

        n_skip = nnz(~kept);
        if n_skip > 0
            fprintf(['  Interface %-4s: %d of %d nodes are clamped and cannot ' ...
                'come into contact, dropped.\n'], s.set, n_skip, numel(ids));
        end
        g = [g; gi(kept)];                           %#ok<AGROW>

        info.nodes{i}  = ids(kept);
        info.labels{i} = s.set;
        info.blocks(i) = numel(bnd) - n_blk0;
        info.rows{i}   = (r0 - nnz(kept) + 1):r0;
    end

    N = sparse(rows, cols, vals, r0, n_free);

    [u_bnd, ia] = unique(bnd, 'stable');
    if numel(u_bnd) ~= numel(bnd)
        error('contact_operator:SharedDof', ...
            ['%d DOFs belong to more than one interface. Overlapping contact ' ...
             'sets would make the interface partition ambiguous; remove the ' ...
             'shared nodes from one of the sets.'], numel(bnd) - numel(u_bnd));
    end
    info.bnd_dofs = bnd(sort(ia));
end

% =====================================================================
function s = normalize_interface_spec(spec)
%NORMALIZE_INTERFACE_SPEC Bring both accepted formats to one struct array.
% The legacy format is the cell array {label, direction, signed_gap; ...},
% where the sign of the gap told which side of the DOF the wall was on. It
% maps to a normal along that axis with the sign of the gap, and a positive
% gap equal to its magnitude.

    if iscell(spec)
        n = size(spec, 1);
        s = repmat(struct('set', '', 'normal', [], 'gap', [], ...
                          'plane', [], 'dofs', 'normal'), 1, n);
        for i = 1:n
            d = spec{i,2}; gp = spec{i,3};
            if gp == 0
                error('contact_operator:ZeroGap', ...
                    'Interface ''%s'': a zero gap gives no side for the wall.', spec{i,1});
            end
            nrm = zeros(1, max(d, 2));
            nrm(d) = sign(gp);
            s(i).set = spec{i,1};
            s(i).normal = nrm;
            s(i).gap = abs(gp);
        end
        return
    end

    if ~isstruct(spec)
        error('contact_operator:BadSpec', ...
            'The interfaces must be a struct array or the legacy cell array.');
    end

    s = spec(:)';
    for i = 1:numel(s)
        if ~isfield(s, 'dofs') || isempty(s(i).dofs), s(i).dofs = 'normal'; end
        if ~isfield(s, 'plane'), s(i).plane = []; end
        if ~isfield(s, 'gap'),   s(i).gap = []; end
        if isempty(s(i).gap) && isempty(s(i).plane)
            error('contact_operator:NoGap', ...
                'Interface ''%s'': give either gap or plane.', s(i).set);
        end
        if ~isempty(s(i).gap) && ~isempty(s(i).plane)
            error('contact_operator:TwoGaps', ...
                'Interface ''%s'': gap and plane are alternatives, not both.', s(i).set);
        end
        if isempty(s(i).normal)
            error('contact_operator:NoNormal', ...
                'Interface ''%s'': the normal is missing.', s(i).set);
        end
    end
end
