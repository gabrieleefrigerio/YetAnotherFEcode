function [fa, act] = contact_solve(pen, k_pen, Scomp)
%CONTACT_SOLVE Contact force and the active set consistent with it.
%
%   [fa, act] = CONTACT_SOLVE(pen, k_pen, Scomp)
%
% INPUT
%   pen    [n_c x 1] penetration N*q - gap, positive where a constraint is
%          violated
%   k_pen  [n_c x 1] penalty stiffness per constraint
%   Scomp  [n_c x n_c] contact compliance, or [] for the plain penalty law
%
% OUTPUT
%   fa     [nnz(act) x 1] force on the ACTIVE constraints only
%   act    [n_c x 1] logical active set, matching fa in length
%
% Shared by every integrator so that the force, the tangent and the residual
% can never disagree about which constraints are touching. It was a nested
% function of TransientSolverOde until the Newmark path needed the same
% answer; the body below is unchanged.
    % Without the compliance correction each constraint is
    % independent: it is active where the penetration is positive,
    % and the force is k*p. That is the plain penalty law.
    %
    % With the correction the constraints are COUPLED, because the
    % interface deforms under the contact load and that deformation
    % is felt by the neighbours. The penetration a spring actually
    % sees is then
    %
    %       p_true = p - S*f
    %
    % and activity has to be decided on p_true, not on p. Deciding
    % it on p and then solving the coupled system on that set is
    % what a first version of this did, and it is WRONG in a way
    % that is expensive rather than merely inaccurate: a node
    % entering the set changes the matrix being inverted, so every
    % activation makes the force on all the other nodes jump. On the
    % 3D model that was a 24% discontinuity at every event, and
    % ode15s restarted at each one - the run went from 10 seconds to
    % over twenty minutes without finishing.
    %
    % Iterating the set to consistency removes the discontinuity:
    % measured, the largest jump between adjacent samples falls from
    % 23.6% of the force scale to 0.34%, which is just the sampling
    % of a continuous curve. Since diag(1/k)+S is positive definite
    % the complementarity problem has a unique solution and the loop
    % converges in 1.23 iterations on average, 2 at worst, so the
    % correction costs barely more than the single solve it replaces.
    act = pen > 0;
    if isempty(Scomp)
        fa = k_pen(act) .* pen(act);
        return
    end

    for sweep = 1:40
        if any(act)
            fa = (diag(1 ./ k_pen(act)) + Scomp(act, act)) \ pen(act);
        else
            fa = zeros(0, 1);
        end

        f_full = zeros(numel(pen), 1);
        f_full(act) = fa;
        p_true = pen - Scomp * f_full;

        act_new = act;
        act_new(act  & f_full <= 0) = false;   % pulled, so not touching
        act_new(~act & p_true  > 0) = true;    % pushed in by a neighbour
        if isequal(act_new, act), return; end
        act = act_new;
    end

    % Did not settle in the sweep budget (a degenerate active set
    % can cycle). Return a force CONSISTENT with the final act -
    % solving on it and dropping any pulling node - so that fa and
    % act always match in length. Omitting this is what let a stale
    % fa reach the caller and crash the state-space product.
    if any(act)
        fa = (diag(1 ./ k_pen(act)) + Scomp(act, act)) \ pen(act);
        keep = fa > 0;
        if ~all(keep)
            idx = find(act);
            act(idx(~keep)) = false;
            fa = fa(keep);
        end
    else
        fa = zeros(0, 1);
    end
end
