function [Mr2, Kr2, Cr2, Phi_CC, info] = interface_reduction(Mr, Kr, Cr, n_bnd, n_cc, mode, iface_blocks, pencil, static_correction, cc_alloc, refine)
%INTERFACE_REDUCTION Secondary modal reduction of the interface partition of a CMS ROM.
%
%   [Mr2, Kr2, Cr2, Phi_CC, info] = INTERFACE_REDUCTION(Mr, Kr, Cr, n_bnd, n_cc, mode, iface_blocks, pencil, static_correction)
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
%     static_correction  logical (optional). Adds back the quasi-static
%                  flexibility of the truncated CC modes as info.R_res, a local
%                  contact compliance. Defined for the 'global' basis ONLY; for
%                  'per_interface' it is skipped and info.R_res stays empty.
%     cc_alloc     optional per-face mode allocation, passed to cc_modes;
%                  'per_interface' only. Empty pools the modes by frequency.
%     refine       logical (optional). Refines the CC modes with the residual
%                  modal flexibility, after Ahn et al., Mech. Syst. Signal
%                  Process. 178 (2022) 109265. Independent of, and composable
%                  with, static_correction: that one corrects the contact LAW,
%                  this one corrects the BASIS and so the reduced matrices
%                  themselves. The reduced size is IDENTICAL either way - the
%                  method buys accuracy without adding a single coordinate.
%                  Derived for the ROM's own boundary block ('self' pencil);
%                  with the Guyan pencil it warns and stays approximate.
%                  'global' basis only.
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
% Delegated to cc_modes so the reduction and the plotting scripts share one
% implementation and cannot drift apart. cc_alloc (optional) fixes how many CC
% modes each contact face contributes; empty leaves cc_modes to pool by
% frequency. n_cc is re-read from the result so everything downstream - the
% size of T_CC, the file name, the R_res guard - uses the count that was
% actually produced.
if nargin < 10, cc_alloc = []; end
[Phi_CC, w2, modes_per_face] = cc_modes(K_bb, M_bb, n_cc, mode, iface_blocks, DENSE_LIMIT, cc_alloc);
n_cc = size(Phi_CC, 2);
% ---------- residual flexibility of the truncated CC modes (GLOBAL only) ----------
% Truncation drops the CC modes above n_cc as if they did not exist. They do
% exist, and their static deflection is the ingredient BOTH corrections below
% are built from:
%
%       R_res = K_bb^-1 - Phi (Phi' K_bb Phi)^-1 Phi' = sum_{i>n_cc} phi_i phi_i'/w_i^2
%
% Computed in the deflated form K_bb\(I - M_bb Phi Phi'), which is the same
% matrix written so the cancellation happens inside a well-conditioned
% projector instead of between two large nearly-equal inverses.
%
% Defined for the 'global' basis only: there Phi is the M_bb-orthonormal set of
% eigenvectors of the whole boundary pencil. The 'per_interface' variant is
% orthonormal only within each face's block while M_bb couples the faces, so
% neither correction is defined for it and both are skipped - those runs use the
% plain penalty law with the plain basis, exactly like an uncorrected IR.
if nargin < 11 || isempty(refine), refine = false; end
want_static = nargin >= 9 && ~isempty(static_correction) && static_correction && n_cc < n_bnd;
want_refine = ~isempty(refine) && refine && n_cc < n_bnd;
is_global   = strcmpi(mode, 'global');

R_res = [];
if (want_static || want_refine) && is_global
    R_res = K_bb \ (eye(n_bnd) - M_bb * (Phi_CC * Phi_CC'));
    R_res = (R_res + R_res') / 2;
elseif want_static || want_refine
    fprintf('  [IR] corrections are only available for the global CC basis; skipped for ''%s''\n', ...
        lower(mode));
    want_static = false;  want_refine = false;
end

% ---------- secondary transformation ----------
T_CC = blkdiag(Phi_CC, eye(m));

% Refined CC modes, Ahn et al., MSSP 178 (2022) 109265.
%
% The truncated interface modes are not merely absent from the contact law -
% they also carry inertia loads that the substructure modes apply to the
% interface, and dropping them degrades the ROM's own dynamics. Projecting the
% interface equation onto the truncated modes and keeping the quasi-static term
% gives (their Eq. 20a)
%
%       du_b = w^2 R_res M_bi q_i .
%
% The interface's own inertia M_bb does NOT appear: it cancels through the
% M_bb-orthogonality between kept and truncated CC modes, which is why this is
% restricted to the global basis. O'Callahan (their Eq. 22) replaces w^2 times
% the state by (M\K) times the state, so the q_i rows of X = Mr2\Kr2 supply the
% missing w^2 q_i and the transformation becomes frequency independent.
%
% Two consequences, both visible in T_CC: the CC modes themselves change, and
% the zero block coupling interface to substructure modes becomes non-zero. The
% SIZE is untouched - this buys accuracy at n_cc + m coordinates, exactly as
% many as the plain reduction, which is the whole point of the method.
%
% Written with X = Mr2\Kr2 rather than the paper's Banachiewicz expansion (their
% Eq. 24). They are equivalent, but the expansion assumes the reduced M has unit
% diagonal blocks and the reduced K is block diagonal: true for CB, not for
% Rubin or MC, which this function also serves.
if want_refine
    if ~strcmpi(basis_src, 'self')
        warning('IR:RefineBasis', ...
            ['Refinement assumes the CC modes are M-orthogonal in the ROM''s own ' ...
             'boundary block (the ''self'' pencil, as derived by Ahn et al.). With ' ...
             'the ''%s'' pencil that orthogonality does not hold and the M_bb term ' ...
             'no longer cancels, so the correction is approximate.'], basis_src);
    end
    Mr2_0 = T_CC' * Mr * T_CC;   Mr2_0 = full((Mr2_0 + Mr2_0') / 2);
    Kr2_0 = T_CC' * Kr * T_CC;   Kr2_0 = full((Kr2_0 + Kr2_0') / 2);
    M_bi  = full(Mr(ib, n_bnd+1:end));
    X     = Mr2_0 \ Kr2_0;
    T_CC(ib, :) = T_CC(ib, :) + R_res * M_bi * X(n_cc+1:end, :);
    fprintf('  [IR] CC modes refined with the residual modal flexibility of the %d truncated modes\n', ...
        n_bnd - n_cc);
end

Mr2 = T_CC' * Mr * T_CC;   Mr2 = (Mr2 + Mr2') / 2;
Kr2 = T_CC' * Kr * T_CC;   Kr2 = (Kr2 + Kr2') / 2;
Cr2 = T_CC' * Cr * T_CC;   Cr2 = (Cr2 + Cr2') / 2;

% ---------- report ----------
f_cc = sqrt(max(w2, 0)) / (2*pi);
info = struct('mode', lower(mode), 'basis', basis_src, 'n_cc', n_cc, 'n_bnd', n_bnd, ...
              'f_cc', f_cc, 'modes_per_face', modes_per_face, 'refined', want_refine);

% The full transformation, so callers never have to assume it is block diagonal:
% once refined it is not. Reconstruction, the contact operator and the projected
% forcing all go through this one matrix.
info.T = T_CC;

% Handed to the solver only for the series-compliance correction, which puts the
% contact spring in series with the interface compliance. When only the
% refinement was asked for, the basis already carries the residual effect and the
% solver must use the plain penalty law.
info.R_res = [];
if want_static
    info.R_res = R_res;
    fprintf('  [IR] residual flexibility of the %d truncated modes retained (contact compliance)\n', ...
        n_bnd - n_cc);
end

fprintf('  [IR] %s / %s basis | %d/%d interface DOFs retained | CC freq %.3e - %.3e Hz\n', ...
    lower(mode), basis_src, n_cc, n_bnd, f_cc(1), f_cc(end));
if ~isempty(modes_per_face)
    fprintf('  [IR] modes per face: %s\n', mat2str(modes_per_face));
end

end
