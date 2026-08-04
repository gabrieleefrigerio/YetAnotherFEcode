function y_contact = extract_contact_response(Struct, Interfaces, labels, q_constrained)
%EXTRACT_CONTACT_RESPONSE Time histories of the contact nodes, per interface.
%
%   y_contact = EXTRACT_CONTACT_RESPONSE(Struct, Interfaces, labels, q)
%
%   Struct        AbaqusStructure object, needed for unconstrain_vector
%   Interfaces    interface metadata struct produced by the main
%   labels        cell array of the interface labels to extract
%   q_constrained displacements on the free DOFs, [n_dofs_c x n_time]
%
%   Returns a struct with one field per interface:
%       y_contact.<label>.X   [n_nodes x n_time]  displacements along X
%       y_contact.<label>.Y   [n_nodes x n_time]  displacements along Y
%
%   Both components are stored regardless of the contact direction of the
%   interface: the normal one drives the gap, the tangential one shows the
%   sliding along the wall.

y_full    = Struct.AssemblyObj.unconstrain_vector(q_constrained);
y_contact = struct();

for i = 1:numel(labels)
    lbl = labels{i};
    if ~isfield(Interfaces, lbl) || isempty(Interfaces.(lbl).nodes)
        continue;
    end
    y_contact.(lbl).X = y_full(Interfaces.(lbl).global_X, :);
    y_contact.(lbl).Y = y_full(Interfaces.(lbl).global_Y, :);
end
end
