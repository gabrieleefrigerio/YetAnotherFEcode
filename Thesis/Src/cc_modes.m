function [Phi_CC, w2, modes_per_face] = cc_modes(K_bb, M_bb, n_cc, mode, iface_blocks, dense_limit, alloc)
%CC_MODES Characteristic-constraint (interface) modes of a condensed pencil.
%
%   [Phi_CC, w2, modes_per_face] = CC_MODES(K_bb, M_bb, n_cc, mode, iface_blocks)
%   [...] = CC_MODES(..., dense_limit)
%   [...] = CC_MODES(..., dense_limit, alloc)          % per_interface only
%
% CC modes of the interface pencil (K_bb, M_bb), mass normalized, one per
% column. This is the single source of truth for the interface reduction
% basis: interface_reduction reduces the ROM with exactly these vectors, and
% the plotting/diagnostic scripts draw exactly these, so the two never drift.
%
%   mode = 'global'         one eigenproblem on the whole boundary partition
%                           (Kuether et al. 2017). The n_cc lowest-frequency
%                           modes are kept. modes_per_face returns [].
%   mode = 'per_interface'  one eigenproblem per contact face. HOW the modes
%                           are then selected depends on alloc:
%
%     alloc empty (default) POOL across faces, sort by frequency, keep the n_cc
%                           lowest. The allocation between faces is whatever the
%                           spectrum dictates, so a face with softer modes gets
%                           more of them. This is the Aoyama / H-CC flavour.
%
%     alloc scalar N        FIXED number per face: the N lowest-frequency modes
%                           of every face, total N*n_faces. n_cc is then only a
%                           sanity check (must equal N*n_faces). Use this to
%                           give every interface the same resolution regardless
%                           of its individual spectrum.
%
%     alloc vector [n1..nf] explicit per-face counts, total sum(alloc).
%
% Each per_interface mode is supported on a single face (Phi has block
% structure), so a deformation localized on one face lives in that face's
% subspace instead of being synthesized by cancellation between modes spread
% over all faces. modes_per_face returns how many kept modes came from each.
%
% iface_blocks is required for 'per_interface': a cell array with the index
% range of each contact face within 1:n_bnd, which must partition it exactly.
%
% See also INTERFACE_REDUCTION, GUYAN_INTERFACE_PENCIL.

    if nargin < 6 || isempty(dense_limit), dense_limit = 2000; end
    if nargin < 7, alloc = []; end
    K_bb = (K_bb + K_bb') / 2;
    M_bb = (M_bb + M_bb') / 2;
    n_bnd = size(K_bb, 1);

    switch lower(mode)
        case 'global'
            if ~isempty(alloc)
                error('CC_MODES:AllocGlobal', ...
                    'A per-face allocation only applies to mode ''per_interface''.');
            end
            [Phi_CC, w2]   = solve_cc(K_bb, M_bb, n_cc, dense_limit);
            modes_per_face = [];

        case 'per_interface'
            if nargin < 5 || isempty(iface_blocks)
                error('CC_MODES:NoBlocks', ...
                    'mode ''per_interface'' requires iface_blocks, one index vector per face.');
            end
            [Phi_CC, w2, modes_per_face] = ...
                solve_cc_per_face(K_bb, M_bb, n_cc, iface_blocks, n_bnd, dense_limit, alloc);

        otherwise
            error('CC_MODES:BadMode', ...
                'Unknown mode ''%s'': use ''global'' or ''per_interface''.', mode);
    end
end

% =====================================================================
function [Phi, w2] = solve_cc(K, M, n_keep, dense_limit)
% Lowest n_keep modes of (K - w^2 M) phi = 0, mass normalized.
n = size(K, 1);

if n <= dense_limit || n_keep >= n
    [V, D] = eig(K, M, 'chol');
    [w2_all, idx] = sort(real(diag(D)), 'ascend');
    V = real(V(:, idx));
    Phi = V(:, 1:n_keep);
    w2  = w2_all(1:n_keep);
else
    [V, D] = eigs(sparse(K), sparse(M), n_keep, 'smallestabs');
    [w2, idx] = sort(real(diag(D)), 'ascend');
    Phi = real(V(:, idx));
end

% Mass normalization with respect to the M actually used in the eigenproblem
for i = 1:size(Phi, 2)
    nrm = sqrt(Phi(:,i)' * M * Phi(:,i));
    if nrm > 0
        Phi(:,i) = Phi(:,i) / nrm;
    end
end
end

% =====================================================================
function [Phi, w2, modes_per_face] = solve_cc_per_face(K_bb, M_bb, n_cc, iface_blocks, n_bnd, dense_limit, alloc)
% One eigenproblem per contact face. Each resulting mode is supported on a
% single face, so Phi has block structure: a deformation localized on one face
% lives in that face's subspace instead of having to be synthesized by
% cancellation between modes spread over all the faces.
%
% The selection between faces follows alloc:
%   empty   pool every face's modes, sort by frequency, keep the n_cc lowest.
%   scalar  keep that many lowest-frequency modes from EACH face.
%   vector  keep alloc(f) lowest from face f.

n_faces = numel(iface_blocks);
if nargin < 7, alloc = []; end

% The blocks must partition 1:n_bnd exactly, otherwise the pooled basis would
% either miss interface DOFs or count some twice.
all_idx = sort([iface_blocks{:}]);
if ~isequal(all_idx(:)', 1:n_bnd)
    error('IR:BadBlocks', ...
        ['iface_blocks must be a partition of 1:%d (found %d indices, %d unique). ' ...
         'Check how the interface blocks were recorded in the main.'], ...
        n_bnd, numel(all_idx), numel(unique(all_idx)));
end

face_sizes = cellfun(@numel, iface_blocks);

% Resolve a per-face allocation from alloc, if given, and validate it.
if ~isempty(alloc)
    if isscalar(alloc), alloc = repmat(alloc, 1, n_faces); end
    alloc = alloc(:)';
    if numel(alloc) ~= n_faces
        error('CC_MODES:BadAlloc', ...
            'alloc has %d entries but there are %d contact faces.', numel(alloc), n_faces);
    end
    if any(alloc < 0) || any(alloc ~= round(alloc))
        error('CC_MODES:BadAlloc', 'alloc must hold non-negative integers.');
    end
    if any(alloc > face_sizes)
        bad = find(alloc > face_sizes, 1);
        error('CC_MODES:AllocTooLarge', ...
            ['Asked for %d modes on face %d but it has only %d DOFs. ' ...
             'A face cannot supply more CC modes than its own DOFs.'], ...
            alloc(bad), bad, face_sizes(bad));
    end
    if n_cc ~= sum(alloc)
        error('CC_MODES:AllocMismatch', ...
            ['n_cc = %d does not match the requested per-face allocation ' ...
             '(sum = %d). They must agree so the reduced size is unambiguous.'], ...
            n_cc, sum(alloc));
    end
end

% --- solve each face separately ---
Phi_pool  = zeros(n_bnd, n_bnd);   % at most n_bnd modes in total
w2_pool   = zeros(n_bnd, 1);
face_pool = zeros(n_bnd, 1);
filled    = 0;

for f = 1:n_faces
    idx = iface_blocks{f}(:)';
    nf  = numel(idx);
    if nf == 0, continue; end

    Kf = K_bb(idx, idx);  Kf = (Kf + Kf') / 2;
    Mf = M_bb(idx, idx);  Mf = (Mf + Mf') / 2;

    if isempty(alloc)
        % Keep every mode of the face; pooling below decides the split.
        keep = nf;
    else
        % Keep only the alloc(f) lowest of this face.
        keep = alloc(f);
    end
    if keep == 0, continue; end

    [Phi_f, w2_f] = solve_cc(Kf, Mf, keep, dense_limit);

    rows = filled + (1:keep);
    Phi_pool(idx, rows) = Phi_f;   % zero outside this face: block structure
    w2_pool(rows)       = w2_f;
    face_pool(rows)     = f;
    filled              = filled + keep;
end

% --- select ---
if isempty(alloc)
    % Pool across faces, sort by frequency, keep the n_cc lowest.
    [w2_sorted, ord] = sort(w2_pool(1:filled), 'ascend');
    ord = ord(1:n_cc);
else
    % Everything already collected is kept; still order by frequency so the
    % columns of Phi come out low-to-high like the global variant.
    [w2_sorted, ord] = sort(w2_pool(1:filled), 'ascend');
end

Phi = Phi_pool(:, ord);
w2  = w2_sorted;

modes_per_face = accumarray(face_pool(ord), 1, [n_faces, 1])';
end
