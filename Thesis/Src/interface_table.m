function [faces, face_gap, face_nrm, d_ref] = interface_table(R, verbose)
%INTERFACE_TABLE Per-interface plotting data, read from a run_config.
%
%   [faces, face_gap, face_nrm, d_ref] = INTERFACE_TABLE(R)
%   [...] = INTERFACE_TABLE(R, false)          % no printout
%
% R is the struct loaded from a run_config.mat. Everything is read from the
% file rather than restated by hand, so a post-processing script needs no
% prior knowledge of the model it is looking at.
%
% OUTPUT
%   faces     {1 x n} interface labels
%   face_gap  [1 x n] positive gap of each interface
%   face_nrm  {1 x n} wall normal of each interface, [] when unknown
%   d_ref     reference direction for the 'physical' view: the shock
%             direction, so the panels show the motion the loading drives
%
% See also CONTACT_PLOT_QUANTITY.

    if nargin < 2, verbose = true; end

    faces   = R.active_labels;
    n_faces = numel(faces);
    face_gap = zeros(1, n_faces);
    face_nrm = cell(1, n_faces);
    for i = 1:n_faces
        I = R.Interfaces.(faces{i});
        if isfield(I, 'gap_nodes') && ~isempty(I.gap_nodes)
            face_gap(i) = I.gap_nodes(1);
        else
            face_gap(i) = abs(I.gap);
        end
        if isfield(I, 'normal'), face_nrm{i} = I.normal(:)'; else, face_nrm{i} = []; end
    end

    if isfield(R, 'cfg') && isfield(R.cfg, 'impulse_dir') && ~isempty(R.cfg.impulse_dir)
        d_ref = R.cfg.impulse_dir(:)' / norm(R.cfg.impulse_dir);
    else
        d_ref = [0 0 1];
    end

    if verbose
        fprintf('Interfaces: %d\n', n_faces);
        for i = 1:n_faces
            if isempty(face_nrm{i}), nstr = 'n/a'; else, nstr = mat2str(face_nrm{i}); end
            fprintf('  %-12s n = %-12s gap %.3f  (%d nodes)\n', ...
                faces{i}, nstr, face_gap(i), numel(R.Interfaces.(faces{i}).nodes));
        end
    end
end
