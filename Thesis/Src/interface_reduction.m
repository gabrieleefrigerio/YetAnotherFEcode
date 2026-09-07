function [Mr2, Kr2, Cr2, Phi_CC, info] = interface_reduction(Mr, Kr, Cr, n_bnd, n_cc, mode, iface_blocks, pencil, static_correction, cc_alloc)
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
% ---------- secondary transformation ----------
T_CC = blkdiag(Phi_CC, eye(m));

Mr2 = T_CC' * Mr * T_CC;   Mr2 = (Mr2 + Mr2') / 2;
Kr2 = T_CC' * Kr * T_CC;   Kr2 = (Kr2 + Kr2') / 2;
Cr2 = T_CC' * Cr * T_CC;   Cr2 = (Cr2 + Cr2') / 2;

% ---------- report ----------
f_cc = sqrt(max(w2, 0)) / (2*pi);
info = struct('mode', lower(mode), 'basis', basis_src, 'n_cc', n_cc, 'n_bnd', n_bnd, ...
              'f_cc', f_cc, 'modes_per_face', modes_per_face);

% ---------- residual flexibility of the truncated CC modes (GLOBAL only) ----------
% Truncation drops the CC modes above n_cc as if they did not exist. They do
% exist, but at 1e8-1e9 Hz against an excitation of about 1 MHz, so they carry
% no dynamics: they deflect quasi-statically under the contact load. Adding
% their static contribution back as a local compliance keeps the contact from
% being artificially STIFF.
%
% This correction is defined for the 'global' basis only. There Phi is the
% M_bb-orthonormal set of eigenvectors of the whole boundary pencil, so
%
%       R_res = K_bb^-1 - Phi * (Phi' K_bb Phi)^-1 * Phi'
%              = sum_{i>n_cc} phi_i phi_i' / w_i^2
%
% is the residual MODAL flexibility of the modes above n_cc. The
% 'per_interface' variant is NOT M_bb-orthonormal (each face's modes are
% orthonormal only within their own block, while M_bb couples the faces), so
% this static correction is not defined for it and is skipped: those runs use
% the plain penalty law, exactly as an uncorrected interface reduction.
info.R_res = [];
want_static = nargin >= 9 && ~isempty(static_correction) && static_correction && n_cc < n_bnd;
if want_static && strcmpi(mode, 'global')
    Kinv = K_bb \ eye(n_bnd);
    info.R_res = Kinv - Phi_CC * ((Phi_CC' * K_bb * Phi_CC) \ Phi_CC');
    info.R_res = (info.R_res + info.R_res') / 2;
    fprintf('  [IR] residual flexibility of the %d truncated modes retained\n', ...
        n_bnd - n_cc);
elseif want_static
    fprintf('  [IR] static correction is only available for the global CC basis; skipped for ''%s''\n', ...
        lower(mode));
end

fprintf('  [IR] %s / %s basis | %d/%d interface DOFs retained | CC freq %.3e - %.3e Hz\n', ...
    lower(mode), basis_src, n_cc, n_bnd, f_cc(1), f_cc(end));
if ~isempty(modes_per_face)
    fprintf('  [IR] modes per face: %s\n', mat2str(modes_per_face));
end

end
