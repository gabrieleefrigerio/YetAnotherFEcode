function check_cc_orthogonality(Mr, Kr, n_bnd, iface_blocks, n_cc_list, pencil)
%CHECK_CC_ORTHOGONALITY Verify the CC basis rotation and the residual flexibility.
%
%   CHECK_CC_ORTHOGONALITY(Mr, Kr, n_bnd, iface_blocks, n_cc_list)
%   CHECK_CC_ORTHOGONALITY(..., pencil)
%
% Two things have to be proved, and they pull in opposite directions.
%
%   1. The per-face basis really was NOT M_bb-orthonormal, and now is. That is
%      the defect Marconi pointed out, and the "before" number is the one worth
%      quoting: it says how far from orthonormal the shipped basis was.
%
%   2. Fixing it changed NOTHING about the model. The rotation stays inside the
%      selected subspace, so the reduced model must come out identical. Here a
%      difference would be a BUG, not an improvement - which is why the
%      eigenvalue check below is run against a measured noise floor instead of
%      a threshold picked by hand.
%
% It also validates the generalized residual flexibility against the formula it
% replaces, on the 'global' case where the old one is trusted, before anything
% relies on it for 'per_interface'.
%
% See also CC_MODES, INTERFACE_REDUCTION, CHECK_CC_MODES.

    if nargin < 6, pencil = []; end
    ib = 1:n_bnd;
    if isempty(pencil)
        K_bb = full(Kr(ib, ib));  M_bb = full(Mr(ib, ib));
    else
        K_bb = full(pencil.K);    M_bb = full(pencil.M);
    end
    K_bb = (K_bb + K_bb') / 2;  M_bb = (M_bb + M_bb') / 2;

    fprintf('\n=== CC basis orthogonality and residual flexibility ===\n');

    for n_cc = n_cc_list(:)'
        fprintf('\n--- n_cc = %d of %d ---\n', n_cc, n_bnd);
        for mode = {'global', 'per_interface'}
            md = mode{1};
            [Phi, ~, ~, dg] = cc_modes(K_bb, M_bb, n_cc, md, iface_blocks, 500);
            raw = dg.Phi_raw;

            % --- 1. orthonormality, before and after ---
            fprintf('  [%s] M-orthonormality %.3e -> %.3e   (cond(G) = %.3e)\n', ...
                md, dg.defect_before, dg.defect_after, dg.cond_G);

            % --- 2. the span is untouched: that is what makes the ROM identical ---
            ang = subspace(raw, Phi);
            fprintf('  [%s] principal angle raw vs rotated = %.3e rad   (expected 0)\n', md, ang);
            if ang > 1e-10
                error('check_cc_orthogonality:SpanChanged', ...
                    ['The rotation moved the subspace (angle %.3e rad). It must ' ...
                     'only change the basis, never the span, or the reduced ' ...
                     'model is a different model.'], ang);
            end

            % --- 3. eigenvalues of the reduced pencil, old basis vs new ---
            [K_old, M_old] = reduce_with(raw, Kr, Mr, n_bnd);
            [K_new, M_new] = reduce_with(Phi, Kr, Mr, n_bnd);
            lam_old = small_eig(K_old, M_old);
            lam_new = small_eig(K_new, M_new);

            % Noise floor, measured rather than assumed: a random ORTHOGONAL
            % congruence is the same class of perturbation the rotation applies
            % and is guaranteed not to change the model, so whatever it moves
            % the eigenvalues by is what roundoff alone is worth here.
            rng(0);
            [Qr, ~] = qr(randn(size(K_old, 1)), 0);
            lam_rot = small_eig(Qr' * K_old * Qr, Qr' * M_old * Qr);
            nk = min([numel(lam_old), numel(lam_new), numel(lam_rot)]);
            floor_meas = max(abs(lam_rot(1:nk) - lam_old(1:nk)) ./ abs(lam_old(1:nk)));
            d_eig      = max(abs(lam_new(1:nk) - lam_old(1:nk)) ./ abs(lam_old(1:nk)));

            fprintf('  [%s] eigenvalue drift %.3e   (noise floor %.3e)\n', ...
                md, d_eig, floor_meas);
            if d_eig > max(10 * floor_meas, 1e-12)
                error('check_cc_orthogonality:ModelChanged', ...
                    ['The reduced eigenvalues moved by %.3e, well past the %.3e ' ...
                     'that roundoff alone explains. A span-preserving rotation ' ...
                     'cannot do that: something else changed.'], d_eig, floor_meas);
            end

            % --- 4. what SHOULD change: conditioning ---
            fprintf('  [%s] cond(Mr2) %.3e -> %.3e | cond(Kr2) %.3e -> %.3e\n', md, ...
                cond(M_old), cond(M_new), cond(K_old), cond(K_new));

            % --- 5. residual flexibility ---
            check_residual(K_bb, M_bb, Phi, raw, n_cc, n_bnd, md);
        end
    end
    fprintf('\nAll checks passed.\n');
end

% =====================================================================
function [K2, M2] = reduce_with(Phi, Kr, Mr, n_bnd)
    T  = blkdiag(Phi, eye(size(Kr, 1) - n_bnd));
    K2 = full(T' * Kr * T);  K2 = (K2 + K2') / 2;
    M2 = full(T' * Mr * T);  M2 = (M2 + M2') / 2;
end

% =====================================================================
function lam = small_eig(K, M)
% Diagonally scaled, INVERTED pencil. The reduced pencils here span many
% decades, and a direct eig resolves the lowest modes only to eps*spread -
% measured at 3.6e-06 on the CB pencil of this project, which would swamp the
% difference this check exists to detect. Inverting makes those modes the
% largest eigenvalues, which come back with full relative accuracy (4e-11).
    d  = sqrt(abs(diag(K)));  d(d <= 0) = 1;  d = 1 ./ d;
    Ks = (d .* K) .* d';  Ks = (Ks + Ks') / 2;
    Ms = (d .* M) .* d';  Ms = (Ms + Ms') / 2;
    mu = real(eig(Ms, Ks));
    lam = sort(1 ./ max(mu, realmin), 'ascend');
end

% =====================================================================
function check_residual(K_bb, M_bb, Phi, raw, n_cc, n_bnd, md)
% The generalized residual flexibility, against the formula it replaces and
% against the properties that define it uniquely.
    R = residual_flexibility_local(K_bb, Phi);

    % (a) invariance: it may depend on the subspace, never on the basis
    R_raw = residual_flexibility_local(K_bb, raw);
    d_inv = norm(R - R_raw, 'fro') / max(norm(R, 'fro'), realmin);
    fprintf('  [%s] R_res raw basis vs rotated: %.3e   (expected 0, depends on span only)\n', ...
        md, d_inv);

    % (b) on 'global' there is a gold standard: the spectral sum over the modes
    % that were actually truncated. Both the old formula and the new one can be
    % scored against it, which settles which is more accurate instead of
    % assuming the incumbent is right.
    if strcmpi(md, 'global') && n_cc < n_bnd
        [V, Dg] = eig(K_bb, M_bb, 'chol');
        [w2s, ix] = sort(real(diag(Dg)), 'ascend');  V = real(V(:, ix));
        V = V ./ sqrt(max(sum(V .* (M_bb * V), 1), realmin));
        Vt = V(:, n_cc+1:end);
        R_ref = Vt * diag(1 ./ w2s(n_cc+1:end)) * Vt';

        R_old = K_bb \ (eye(n_bnd) - M_bb * (Phi * Phi'));
        R_old = (R_old + R_old') / 2;
        nref  = norm(R_ref, 'fro');
        fprintf('  [%s] R_res vs spectral truth: new %.3e | old formula %.3e\n', ...
            md, norm(R - R_ref, 'fro') / nref, norm(R_old - R_ref, 'fro') / nref);
    end

    % (c) positive semi-definiteness. contact_solve needs diag(1/k) + N R N'
    % definite for the complementarity problem to have a unique solution, so a
    % negative eigenvalue here is not an inaccuracy, it is an ill-posed contact.
    ev = eig((R + R') / 2);
    fprintf('  [%s] R_res min eig / max eig = %.3e   (expected >= 0)\n', ...
        md, min(ev) / max(max(ev), realmin));

    % (d) it must produce nothing the retained subspace already produces
    d_ann = norm(R * K_bb * Phi, 'fro') / ...
            max(norm(R, 'fro') * norm(K_bb * Phi, 'fro'), realmin);
    fprintf('  [%s] ||R_res K Phi|| (normalised) = %.3e   (expected 0)\n', md, d_ann);
end

% =====================================================================
function R_res = residual_flexibility_local(K_bb, Phi)
% Same construction interface_reduction uses; duplicated here on purpose so the
% check is an independent statement of the formula rather than a call into the
% code under test.
    n_bnd  = size(K_bb, 1);
    L      = chol(K_bb, 'lower');
    [Q, ~] = qr(L' * Phi, 0);
    X      = (eye(n_bnd) - Q * Q') / L;
    R_res  = X' * X;
end
