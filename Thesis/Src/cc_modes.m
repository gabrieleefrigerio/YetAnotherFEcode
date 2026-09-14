function [Phi_CC, w2, modes_per_face, cc_diag] = cc_modes(K_bb, M_bb, n_cc, mode, iface_blocks, dense_limit, alloc)
%CC_MODES Characteristic-constraint (interface) modes of a condensed pencil.
%
%   [Phi_CC, w2, modes_per_face] = CC_MODES(K_bb, M_bb, n_cc, mode, iface_blocks)
%   [...] = CC_MODES(..., dense_limit)
%   [...] = CC_MODES(..., dense_limit, alloc)          % per_interface only
%   [Phi_CC, w2, modes_per_face, cc_diag] = CC_MODES(...)
%
% CC modes of the interface pencil (K_bb, M_bb), M_bb-orthonormal over the
% WHOLE boundary partition, one per column. This is the single source of truth
% for the interface reduction basis: interface_reduction reduces the ROM with
% exactly these vectors, and the plotting/diagnostic scripts draw exactly
% these, so the two never drift.
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
% Each per_interface mode is SELECTED on a single face, so a deformation
% localized on one face lives in that face's subspace instead of being
% synthesized by cancellation between modes spread over all faces. That
% argument is about the SPAN, and the span is what the reduction sees.
%
% The selected per-face modes are however only orthonormal inside their own
% face: M_bb couples the faces, so Phi'*M_bb*Phi is nowhere near the identity
% (measured on the 3D accelerometer: defect of order 1, cond = 2e4, because
% the top and bottom faces of one tab sit across 30 um of proof mass and
% couple almost perfectly). They are therefore rotated onto the Ritz basis of
% the SELECTED subspace, which leaves the span - hence the reduced model -
% untouched while making Phi'*M_bb*Phi = I and Phi'*K_bb*Phi diagonal.
%
% Two consequences of that rotation, both intentional:
%   - the returned columns are no longer supported on a single face; the raw
%     per-face modes are handed back in cc_diag for plotting.
%   - w2 holds the RITZ values of the selected subspace, not the isolated
%     per-face eigenvalues. The Ritz values are the honest frequencies of the
%     subspace actually used; the per-face ones are not comparable across
%     faces. modes_per_face therefore describes the SELECTION, not the columns.
%
% cc_diag (optional 4th output) carries Phi_raw, w2_raw, the rotation C with
% Phi_CC = Phi_raw*C, and the measured orthogonality defect before and after.
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

    % Both variants are MEASURED, only per_interface is corrected. 'global'
    % already returns eigenvectors of this very pencil, so its defect is
    % roundoff and rotating it would only flip signs and stir degenerate
    % subspaces, changing saved bases and figures without changing physics.
    % Measuring it anyway costs one k-by-k product and turns a standing
    % assumption into a number - which matters here, where the pencil has been
    % measured at cond = 1e17 and an eigensolver can quietly lose orthogonality.
    [Phi_CC, w2, cc_diag] = orthonormalize_cc(Phi_CC, w2, K_bb, M_bb, ...
                                              strcmpi(mode, 'per_interface'));
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
[w2_sorted, ord] = sort(w2_pool(1:filled), 'ascend');
if isempty(alloc)
    % Pool across faces, sort by frequency, keep the n_cc lowest. BOTH the
    % index list and the frequencies have to be truncated: keeping the full
    % w2 while cutting ord left w2 longer than Phi has columns, so the caller
    % reported f_cc(end) - the top of the CC range, printed and saved into
    % info.f_cc - as the frequency of a mode that had been thrown away.
    ord       = ord(1:n_cc);
    w2_sorted = w2_sorted(1:n_cc);
end
% With a non-empty alloc everything collected is kept (filled = sum(alloc) =
% n_cc); the sort above only orders the columns low-to-high like the global
% variant does.

Phi = Phi_pool(:, ord);
w2  = w2_sorted;

modes_per_face = accumarray(face_pool(ord), 1, [n_faces, 1])';
end

% =====================================================================
function [Phi, w2, cc_diag] = orthonormalize_cc(Phi, w2, K_bb, M_bb, do_rotate)
% Measure the M_bb-orthonormality of the selected basis and, when asked,
% restore it by rotating onto the Ritz basis of the SAME subspace.
%
% Why a Ritz rotation and not a Cholesky/Gram-Schmidt factorisation of the
% Gram matrix, which would also orthonormalise while preserving the span:
%   - it is canonical. A triangular factor makes the result depend on the
%     order the modes happened to be stacked in, which is an implementation
%     detail, not a property of the model.
%   - it diagonalises Phi'*K_bb*Phi as well, so the frequencies handed back
%     are the Ritz values of the subspace instead of per-face eigenvalues
%     that are not comparable across faces.
%   - accuracy. One Cholesky-QR pass loses orthogonality like cond(G)^2*eps,
%     and cond(G) is 2e4 here, so a single pass would leave a defect around
%     1e-7 - better than the 1e0 we start from, but not good enough to call
%     the basis orthonormal.
    n_cc = size(Phi, 2);
    G0   = Phi' * M_bb * Phi;  G0 = (G0 + G0') / 2;
    e0   = norm(G0 - eye(n_cc), 'fro') / sqrt(n_cc);

    cc_diag = struct('Phi_raw', Phi, 'w2_raw', w2, 'C', eye(n_cc), ...
                     'defect_before', e0, 'defect_after', e0, ...
                     'cond_G', cond(G0), 'rotated', false);
    if ~do_rotate, return; end

    Ks = Phi' * K_bb * Phi;  Ks = (Ks + Ks') / 2;
    C  = ritz_rotation(Ks, G0);
    Phi = Phi * C;

    G1 = Phi' * M_bb * Phi;  G1 = (G1 + G1') / 2;
    e1 = norm(G1 - eye(n_cc), 'fro') / sqrt(n_cc);

    % One rotation leaves a defect of order cond(G)*eps, which at cond(G) = 2e4
    % lands around 1e-10 - already six orders better than the 1e0 we started
    % from, but not machine precision. A single Cholesky pass finishes the job,
    % and it is safe HERE precisely because the first rotation already brought
    % the Gram to cond ~ 1: the same pass applied to the raw basis would have
    % lost cond(G)^2*eps.
    if e1 > 1e-12
        [Rc, bad] = chol(G1);
        if bad == 0
            C   = C / Rc;
            Phi = Phi / Rc;
            G1  = Phi' * M_bb * Phi;  G1 = (G1 + G1') / 2;
            e1  = norm(G1 - eye(n_cc), 'fro') / sqrt(n_cc);
        end
    end

    w2 = real(diag(Phi' * K_bb * Phi));   % Ritz values of the selected subspace

    cc_diag.C = C;  cc_diag.defect_after = e1;  cc_diag.rotated = true;
    if e1 > 1e-10
        warning('CC_MODES:OrthoNotAchieved', ...
            ['The Ritz rotation left an M-orthonormality defect of %.2e ' ...
             '(was %.2e). The selected per-face modes are close to linearly ' ...
             'dependent in the M_bb metric, which means two faces couple so ' ...
             'strongly that their modes carry the same motion. Reduce n_cc, ' ...
             'or use mode = ''global''.'], e1, e0);
    end
end

% =====================================================================
function C = ritz_rotation(Ks, Ms)
% Rotation onto the Ritz basis of the subspace: C'*Ms*C = I, C'*Ks*C diagonal.
%
% Solved on the INVERTED, diagonally scaled pencil. The small pencil inherits
% the spread of the retained CC frequencies, and it is the LOWEST ones we care
% about - exactly the ones a direct eig resolves worst, since its relative
% accuracy there goes like eps*(lambda_max/lambda_min). Inverting turns them
% into the largest eigenvalues, which come back with full relative precision.
% Measured on the CB pencil of this project: noise 3.6e-06 direct, 4e-11
% inverted.
    d   = 1 ./ sqrt(abs(diag(Ks)));  d(~isfinite(d)) = 1;
    Ksd = (d .* Ks) .* d';  Ksd = (Ksd + Ksd') / 2;
    Msd = (d .* Ms) .* d';  Msd = (Msd + Msd') / 2;

    [Y, Mu] = eig(Msd, Ksd, 'chol');
    mu = real(diag(Mu));
    Y  = real(Y);

    % eig(...,'chol') normalises against Ksd; rescale to Y'*Msd*Y = I so the
    % rotation is M-orthonormal, which is the property the caller needs.
    nrm = sqrt(max(sum(Y .* (Msd * Y), 1), realmin));
    Y   = Y ./ nrm;

    lam = 1 ./ max(mu, realmin);
    [~, ord] = sort(lam, 'ascend');
    C = d .* Y(:, ord);
end
