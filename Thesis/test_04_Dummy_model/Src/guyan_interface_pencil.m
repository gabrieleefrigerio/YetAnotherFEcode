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
%   Why this exists. The characteristic constraint modes used by the interface
%   reduction are only a Ritz basis for the interface displacement, so any
%   full-rank choice is admissible and the choice only affects accuracy. Taking
%   them from the boundary partition of the ROM being reduced ties the basis to
%   the CMS parameterization, which is fine for Craig-Bampton but degenerate for
%   Rubin: there the interface block is spanned by RESIDUAL attachment modes,
%   which by construction carry no low-frequency content, so its low CC modes
%   are nearly orthogonal to the interface motion the dynamics actually
%   produces. Measured on the dummy model, keeping 6 of 26 modes from the Rubin
%   pencil leaves 94% of the FOM interface motion unrepresentable; from this
%   pencil, 0.01%.
%
%   Taking the basis from the Guyan condensation instead makes it a property of
%   the INTERFACE rather than of the reduction method, which is the direction
%   indicated by Tran (see Krattiger et al., MSSP 114 (2019), Sec. 1, on
%   extending S-CC to free- and hybrid-interface methods by first describing
%   the attachment modes in terms of constraint modes).
%
%   For Craig-Bampton this pencil coincides with the boundary partition of the
%   ROM, because RomCB's first n_bnd basis columns ARE the constraint modes.
%   Switching basis therefore leaves CB unchanged and only affects Rubin.
%
%   Note that Psi is independent of how many fixed-interface modes the ROM
%   keeps, so this can be computed once per model and reused across the sweep.

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
