function y_contact = extract_contact_response(Struct, Interfaces, labels, q_constrained)
%EXTRACT_CONTACT_RESPONSE Time histories of the contact nodes, per interface.
%
%   y_contact = EXTRACT_CONTACT_RESPONSE(Struct, Interfaces, labels, q)
%
%   Struct        FeStructure object, needed for unconstrain_vector
%   Interfaces    interface metadata struct produced by the main
%   labels        cell array of the interface labels to extract
%   q_constrained displacements on the free DOFs, [n_dofs_c x n_time]
%
%   Returns a struct with one field per interface:
%       .U        [n_nodes x n_dim x n_time]  every component, 2D or 3D
%       .normal   [n_nodes x n_time]          component along the wall normal
%       .X, .Y    [n_nodes x n_time]          kept for the 2D post-processing
%       .Z        only in 3D
%
%   The normal component is the one that defines the gap, and the only one
%   the contact condition sees: the penetration of a node is its normal
%   displacement minus the gap. The others describe the sliding along the
%   wall. Storing the normal explicitly is what lets the post-processing work
%   the same way for a wall aligned with an axis and for an oblique one, where
%   no single Cartesian component carries the whole story.

y_full    = Struct.AssemblyObj.unconstrain_vector(q_constrained);
n_dim     = Struct.MeshObj.nDOFPerNode;
n_time    = size(q_constrained, 2);
axis_name = {'X', 'Y', 'Z'};
y_contact = struct();

for i = 1:numel(labels)
    lbl = labels{i};
    if ~isfield(Interfaces, lbl) || isempty(Interfaces.(lbl).nodes)
        continue;
    end
    ids = Interfaces.(lbl).nodes(:);
    n_n = numel(ids);

    U = zeros(n_n, n_dim, n_time);
    for d = 1:n_dim
        Ud = y_full((ids-1)*n_dim + d, :);
        U(:, d, :) = reshape(Ud, n_n, 1, n_time);
        y_contact.(lbl).(axis_name{d}) = Ud;
    end
    y_contact.(lbl).U = U;

    % Projection on the wall normal. Older runs saved no normal, in which
    % case the contact direction alone identifies it.
    if isfield(Interfaces.(lbl), 'normal') && ~isempty(Interfaces.(lbl).normal)
        nrm = Interfaces.(lbl).normal(:);
        nrm = nrm / norm(nrm);
    else
        nrm = zeros(n_dim, 1);
        nrm(Interfaces.(lbl).dir) = sign(Interfaces.(lbl).gap);
    end
    if numel(nrm) < n_dim, nrm(n_dim) = 0; end

    Un = zeros(n_n, n_time);
    for d = 1:n_dim
        % reshape rather than squeeze: squeeze collapses the wrong dimension
        % when the interface has a single node or the history a single instant
        Un = Un + nrm(d) * reshape(U(:, d, :), n_n, n_time);
    end
    y_contact.(lbl).normal = Un;
end
end
