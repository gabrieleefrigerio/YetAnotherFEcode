function [K_bb, M_bb, Psi] = guyan_interface_pencil(Struct, contact_dofs)
%GUYAN_INTERFACE_PENCIL Static condensation of the structure onto the contact DOFs.
%
%   [K_bb, M_bb, Psi] = GUYAN_INTERFACE_PENCIL(Struct, contact_dofs)
%
%   Builds the constraint modes Psi, i.e. the static response of the whole
%   structure to a unit displacement of each contact DOF with all the others
%   held fixed, and returns the condensed pencil
%
%       K_bb = Psi' * Kc * Psi ,      M_bb = Psi' * Mc * Psi
%
%   in PHYSICAL interface coordinates.
%
%   The eigenvectors of this pencil are what Tran calls the INTERFACE MODES:
%   "the normal modes of the whole structure after performing a static
%   condensation to the interface". Reducing the interface by truncating them
%   was proposed for the fixed-interface method by Bourquin (doctoral thesis,
%   Paris VI, 1991; Math. Model. Numer. Anal. 26(3), 1992, where they appear as
%   the eigenmodes of the Poincare-Steklov operator) and extended to the free-
%   and hybrid-interface methods by Tran, Comput. Struct. 79 (2001) 209-222.
%
%   Why the basis must come from HERE and not from the ROM. The CC modes are
%   only a Ritz basis for the interface displacement, so any full-rank choice is
%   admissible and the choice affects accuracy alone. Taking them from the
%   boundary partition of the ROM being reduced ties them to the CMS
%   parameterization. That is harmless for Craig-Bampton but degenerate for
%   Rubin, whose interface block is spanned by RESIDUAL attachment modes and so
%   carries no low-frequency content by construction: its low CC modes come out
%   at 5.9e7 Hz against a response in the kHz range, nearly orthogonal to the
%   interface motion the dynamics produces. Measured on the dummy model, keeping
%   6 of 26 modes from the Rubin pencil leaves 94% of the FOM interface motion
%   unrepresentable and the contact response collapses to zero; from this
%   pencil, 0.01%.
%
%   Relation to Tran's procedure. Tran (Sec. 3.3, Eq. 19-20) reduces Rubin in
%   two steps: first replace the residual attachment block by the constraint
%   modes Psi, then substitute the truncated interface modes for Psi. Using this
%   pencil while leaving the ROM basis alone gives the SAME Ritz subspace, hence
%   the same reduced model up to a change of generalized coordinates. The reason
%   is Tran's Eq. (15): for a constrained substructure the attachment modes are
%   linear combinations of the constraint modes, so Psi already lies in Rubin's
%   space. Writing Psi = Pb*A1 + Pm*A2 and taking the interface rows gives
%   A1 = inv(Pb(bnd,:)), hence
%
%       Psi*Phi_CC = Pb*(Pb(bnd,:)\Phi_CC) + Pm*(A2*Phi_CC)
%
%   whose first term is the interface block obtained here and whose second term
%   already lies in the span of the modal block. Verified numerically on Rubin
%   at phi = 200: subspace angle between the two bases 6e-9 to 5e-8 over
%   n_cc = 6, 14, 26; the part of the difference not absorbed by the modal block
%   is 5e-10 to 3e-8; and Psi is reconstructed from Rubin's space with residual
%   1.1e-9, with the interface coefficient block equal to inv(Pb(bnd,:)) to
%   5.0e-10.
%
%   VALIDITY. That equivalence rests on Tran's Eq. (15), which holds for a
%   CONSTRAINED substructure. Tran states explicitly that it fails for a free
%   one, where the inertial forces balancing the applied unit loads act on the
%   whole body and not only at the interface; there his explicit construction
%   would be needed. The present model is constrained. A free substructure would
%   in any case make K_ii singular and this function fail.
%
%   For Craig-Bampton this pencil coincides with the boundary partition of the
%   ROM, because RomCB's first n_bnd basis columns ARE the constraint modes
%   (verified: 0.000e+00 relative difference on both K and M). Switching basis
%   therefore leaves CB untouched and only affects Rubin.
%
%   Psi does not depend on how many fixed-interface modes the ROM keeps, so this
%   can be computed once per model and reused across the whole sweep.

b = contact_dofs(:);

Mc = Struct.AssemblyObj.constrain_matrix(Struct.M);
Kc = Struct.AssemblyObj.constrain_matrix(Struct.K);

n_dofs_c = size(Kc, 1);
i_idx    = setdiff(1:n_dofs_c, b)';
n_b      = numel(b);

% Constraint modes: one linear solve per interface DOF
K_ii  = Kc(i_idx, i_idx);
K_ib  = Kc(i_idx, b);
Psi_i = -(K_ii \ full(K_ib));

Psi          = zeros(n_dofs_c, n_b);
Psi(b, :)    = eye(n_b);
Psi(i_idx,:) = Psi_i;

K_bb = Psi' * Kc * Psi;   K_bb = (K_bb + K_bb') / 2;
M_bb = Psi' * Mc * Psi;   M_bb = (M_bb + M_bb') / 2;
end
