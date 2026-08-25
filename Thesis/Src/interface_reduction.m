function [Mr2, Kr2, Cr2, Phi_CC, info] = interface_reduction(Mr, Kr, Cr, n_bnd, n_cc, mode, iface_blocks, pencil)
%INTERFACE_REDUCTION Secondary modal reduction of the interface partition of a CMS ROM.
%
%   [Mr2, Kr2, Cr2, Phi_CC, info] = INTERFACE_REDUCTION(Mr, Kr, Cr, n_bnd, n_cc, mode, iface_blocks, pencil)
%
%   Applies to ROMs whose reduced basis keeps the interface as physical
%   coordinates at the HEAD of the vector, i.e. q = [x_b ; q_modal] with
%   x_b of length n_bnd. That is the case for RomCB and RomRubin.
%
%   The characteristic constraint (CC) modes come from a secondary eigenvalue
%   analysis on the boundary partition of the already-reduced matrices
%   (Kuether et al., ASME IDETC 2017, Eq. 4):
%
%       ( K_bb - omega^2 M_bb ) phi_CC = 0
%
%   The n_cc lowest-frequency modes are collected in Phi_CC [n_bnd x n_cc],
%   giving the secondary transformation (Eq. 5)
%
%       T_CC = [ Phi_CC   0 ]      Mr2 = T_CC' * Mr * T_CC   (idem Kr, Cr)
%              [   0      I ]
%
%   Inputs
%     Mr, Kr, Cr   reduced matrices of the CMS ROM, size (n_bnd + m)
%     n_bnd        number of physical interface DOFs, at the head of the basis
%     n_cc         number of CC modes to retain (n_cc = n_bnd -> no truncation)
%     mode         'global'        one eigenproblem on the coupled K_bb, M_bb
%                  'per_interface' one eigenproblem per contact face; the modes
%                                  of all faces are pooled, sorted by frequency,
%                                  and the n_cc lowest are kept
%     iface_blocks cell array of index vectors into 1:n_bnd, one per contact
%                  face. Required by 'per_interface', ignored by 'global'.
%     pencil       optional struct('K',K_bb,'M',M_bb) giving the pencil the CC
%                  modes are computed from, EXPRESSED IN THE ROM's INTERFACE
%                  COORDINATES. Omit or leave empty to use the boundary
%                  partition of Mr/Kr themselves.
%
%                  The CC modes are only a Ritz basis for the interface
%                  displacement, so any full-rank choice is admissible and the
%                  choice affects accuracy alone. Supplying the Guyan pencil of
%                  the structure (see guyan_interface_pencil) instead of using
%                  the ROM's own boundary block reproduces Tran, Comput. Struct.
%                  79 (2001), Sec. 3.3: same Ritz subspace as his Eq. (19)-(20),
%                  up to a change of generalized coordinates. Without it, Rubin
%                  is degenerate - its interface block is spanned by residual
%                  attachment modes and carries no low-frequency content, so the
%                  truncated subspace is nearly orthogonal to the interface
%                  motion the dynamics produces and the response collapses.
%                  See guyan_interface_pencil for the derivation, the numbers
%                  and the validity condition.
%
%   Outputs
%     Mr2, Kr2, Cr2  interface-reduced matrices, size (n_cc + m)
%     Phi_CC         [n_bnd x n_cc] CC modes, mass normalized
%     info           struct with the CC frequencies in Hz, the mode and,
%                    for 'per_interface', how many modes each face contributed
%
%   On 'per_interface': solving one eigenproblem per face is mathematically
%   equivalent to zeroing the off-diagonal blocks of K_bb and M_bb and solving
%   a single eigenproblem, because the spectrum of a block-diagonal matrix is
%   the union of the spectra of its blocks. The per-face loop is used because
%   it scales: with a large interface, several small eigenproblems cost far
%   less than one big one. The equivalence is exercised by the verification
%   script rather than relied upon at runtime.
%
%   This reduction does NOT apply to the massless ROMs (RomMCB, RomMN), whose
%   M_bb is identically zero by construction; the eigenproblem would be
%   singular. Those methods use exact set-valued contact anyway.

% Dense eig is more robust and needs no convergence tuning, but is O(n^3) and
% dense in memory. Above this size switch to the sparse iterative solver.
DENSE_LIMIT = 500;

n_tot = size(Kr, 1);
m     = n_tot - n_bnd;

% ---------- input validation ----------
if n_bnd < 1 || n_bnd > n_tot
    error('IR:BadNbnd', 'n_bnd = %d is not compatible with a ROM of size %d.', n_bnd, n_tot);
end
if n_cc < 1 || n_cc > n_bnd
    error('IR:BadNcc', ...
        'n_cc = %d out of range: must be between 1 and n_bnd = %d.', n_cc, n_bnd);
end

ib = 1:n_bnd;
if nargin >= 8 && ~isempty(pencil)
    % Externally supplied pencil (e.g. the Guyan condensation of the structure),
    % already expressed in the ROM's interface coordinates by the caller.
    basis_src = 'guyan';
    K_bb = full(pencil.K);
    M_bb = full(pencil.M);
    if ~isequal(size(K_bb), [n_bnd n_bnd]) || ~isequal(size(M_bb), [n_bnd n_bnd])
        error('IR:BadPencil', ...
            'The supplied pencil must be %dx%d, got K %s and M %s.', ...
            n_bnd, n_bnd, mat2str(size(K_bb)), mat2str(size(M_bb)));
    end
else
    % Boundary partition of the ROM being reduced (Kuether et al. 2017).
    basis_src = 'self';
    K_bb = full(Kr(ib, ib));
    M_bb = full(Mr(ib, ib));
end
K_bb = (K_bb + K_bb') / 2;
M_bb = (M_bb + M_bb') / 2;

% The secondary eigenproblem needs M_bb positive definite. A Cholesky attempt
% is the direct test; comparing norms of M_bb against K_bb would be comparing
% masses with stiffnesses, which is dimensionally meaningless and gives false
% positives on models with small absolute masses (MEMS).
[~, not_spd] = chol(M_bb);
if not_spd > 0
    error('IR:BoundaryMassNotSPD', ...
        ['M_bb is not positive definite, so the secondary eigenproblem is singular. ' ...
         'This is the case for the massless-boundary ROMs (MCB / MacNeal), whose ' ...
         'M_bb is zero by construction: interface reduction in this form does not ' ...
         'apply to them.']);
end

% ---------- CC modes ----------
switch lower(mode)
    case 'global'
        [Phi_CC, w2] = solve_cc(K_bb, M_bb, n_cc, DENSE_LIMIT);
        modes_per_face = [];

    case 'per_interface'
        if nargin < 7 || isempty(iface_blocks)
            error('IR:NoBlocks', ...
                'mode ''per_interface'' requires iface_blocks, one index vector per contact face.');
        end
        [Phi_CC, w2, modes_per_face] = solve_cc_per_face(K_bb, M_bb, n_cc, iface_blocks, ...
                                                         n_bnd, DENSE_LIMIT);

    otherwise
        error('IR:BadMode', 'Unknown mode ''%s'': use ''global'' or ''per_interface''.', mode);
end

% ---------- secondary transformation ----------
T_CC = blkdiag(Phi_CC, eye(m));

Mr2 = T_CC' * Mr * T_CC;   Mr2 = (Mr2 + Mr2') / 2;
Kr2 = T_CC' * Kr * T_CC;   Kr2 = (Kr2 + Kr2') / 2;
Cr2 = T_CC' * Cr * T_CC;   Cr2 = (Cr2 + Cr2') / 2;

% ---------- report ----------
f_cc = sqrt(max(w2, 0)) / (2*pi);
info = struct('mode', lower(mode), 'basis', basis_src, 'n_cc', n_cc, 'n_bnd', n_bnd, ...
              'f_cc', f_cc, 'modes_per_face', modes_per_face);

fprintf('  [IR] %s / %s basis | %d/%d interface DOFs retained | CC freq %.3e - %.3e Hz\n', ...
    lower(mode), basis_src, n_cc, n_bnd, f_cc(1), f_cc(end));
if ~isempty(modes_per_face)
    fprintf('  [IR] modes per face: %s\n', mat2str(modes_per_face));
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
function [Phi, w2, modes_per_face] = solve_cc_per_face(K_bb, M_bb, n_cc, iface_blocks, n_bnd, dense_limit)
% One eigenproblem per contact face, then pool the modes across faces, sort by
% frequency and keep the n_cc lowest. Each resulting mode is supported on a
% single face, so Phi has block structure: a deformation localized on one face
% lives in that face's subspace instead of having to be synthesized by
% cancellation between modes spread over all the faces.

n_faces = numel(iface_blocks);

% The blocks must partition 1:n_bnd exactly, otherwise the pooled basis would
% either miss interface DOFs or count some twice.
all_idx = sort([iface_blocks{:}]);
if ~isequal(all_idx(:)', 1:n_bnd)
    error('IR:BadBlocks', ...
        ['iface_blocks must be a partition of 1:%d (found %d indices, %d unique). ' ...
         'Check how the interface blocks were recorded in the main.'], ...
        n_bnd, numel(all_idx), numel(unique(all_idx)));
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

    % Keep every mode of the face here; the truncation happens after pooling,
    % so that the frequency ordering decides the allocation between faces.
    [Phi_f, w2_f] = solve_cc(Kf, Mf, nf, dense_limit);

    rows = filled + (1:nf);
    Phi_pool(idx, rows) = Phi_f;   % zero outside this face: block structure
    w2_pool(rows)       = w2_f;
    face_pool(rows)     = f;
    filled              = filled + nf;
end

% --- pool, sort by frequency, truncate ---
[w2_sorted, ord] = sort(w2_pool(1:filled), 'ascend');
ord = ord(1:n_cc);

Phi = Phi_pool(:, ord);
w2  = w2_sorted(1:n_cc);

modes_per_face = accumarray(face_pool(ord), 1, [n_faces, 1])';
end
