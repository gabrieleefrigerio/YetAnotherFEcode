function [F_pen, K_eff, act, fa] = contact_tangent(q, K, Nop, Nop_t, gap, k_pen, Scomp)
%CONTACT_TANGENT Penalty contact force and the tangent consistent with it.
%
%   [F_pen, K_eff, act, fa] = CONTACT_TANGENT(q, K, Nop, Nop_t, gap, k_pen, Scomp)
%
% INPUT
%   q      [n x 1] displacement
%   K      [n x n] structural stiffness
%   Nop    [n_c x n] contact operator, penetration is Nop*q - gap
%   Nop_t  Nop' precomputed, since it is applied at every evaluation
%   gap    [n_c x 1] positive gaps
%   k_pen  [n_c x 1] penalty stiffness per constraint
%   Scomp  [n_c x n_c] contact compliance, or [] for the plain penalty law
%
% OUTPUT
%   F_pen  [n x 1] contact force on the structure
%   K_eff  [n x n] K plus the contact contribution, i.e. d(K*q + F_pen)/dq
%   act    [n_c x 1] logical active set
%   fa     force on the active constraints
%
% The force and the tangent come from ONE call to contact_solve, so they are
% always the exact derivative of one another. Deciding activity twice, once
% for each, is what gives an integrator an inconsistent pair and makes the
% Newton iteration stall on the step where a node enters the set.
%
% See also CONTACT_SOLVE, TRANSIENTSOLVERODE, CONTACT_RESIDUAL.

    n = size(K, 1);
    [fa, act] = contact_solve(Nop*q - gap, k_pen, Scomp);

    if any(act)
        F_pen = Nop_t(:, act) * fa;
        if isempty(Scomp)
            K_eff = K + Nop_t(:, act) * (k_pen(act) .* Nop(act, :));
        else
            Keff_a = (diag(1 ./ k_pen(act)) + Scomp(act, act)) \ Nop(act, :);
            K_eff  = K + Nop_t(:, act) * Keff_a;
        end
    else
        F_pen = zeros(n, 1);
        K_eff = K;
    end
end
