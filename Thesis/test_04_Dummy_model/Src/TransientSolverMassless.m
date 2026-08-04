classdef TransientSolverMassless < handle
    % TRANSIENTSOLVERMASSLESS Semi-explicit integrator for massless-boundary
    % models with frictionless set-valued unilateral contact.
    %
    % Reference: Monjaraz Tec et al., "A massless boundary component mode
    % synthesis method for elastodynamic contact problems", Comput. Struct. 260
    % (2022). Ch. 4 (time stepping), Ch. 6 (frictionless algorithm),
    % App. C (augmented Lagrangian).
    %
    % Required model structure (MacNeal or massless CB):
    %   coordinates q = [q_b ; eta],  q_b = boundary (n_bnd), eta = modal (m)
    %   M = [0 0 ; 0 I]        (ZERO mass at the boundary)
    %   K = [K_bb K_be ; K_eb K_ee]
    %   C = [0 0 ; 0 D_ee]     (no damping at the boundary)
    %
    % Equations solved:
    %   K_bb q_b + K_be eta - W lambda = f_b(t)             (static, boundary)
    %   eta_dd + D_ee eta_d + K_ee eta + K_eb q_b = f_e(t)  (dynamic, interior)
    %   g = g0 + W' q_b ,   0 <= g  _|_  lambda >= 0
    %
    % Gap convention: for a wall on the right, with q_b positive towards it,
    %   W = -I ,  g0 = gap_wall   =>   g = gap_wall - q_b
    % With signed gaps s the caller passes W = -diag(sign(s)) and g0 = |s|.
    %
    % The integrator advances on a FIXED step dt and returns its own uniform
    % time grid; it has no equivalent of the 'OutputTimes' option of ode15s.

    properties
        % --- Reduced matrices ---
        Kbb, Kbe, Kee, Dee
        n_bnd, n_mod

        % --- Contact ---
        W          % [n_bnd x n_c] contact direction matrix
        g0         % [n_c x 1] initial gaps
        n_c

        % --- Precomputed quantities ---
        Kbb_fact   % Cholesky decomposition of K_bb
        G_full     % W' * inv(K_bb) * W   [n_c x n_c]  (Delassus matrix)
        Ainv       % inv( (1/dt)*I + 0.5*Dee )
        Bmat       % (1/dt)*I - 0.5*Dee

        % --- Augmented Lagrangian parameters ---
        eps_AL              % Relaxation step
        max_iter_AL = 1000
        tol_AL      = 1e-8  % RELATIVE to ||c||_inf (normalized KKT residual)

        % --- Diagnostics ---
        stats
        warned_AL = false
    end

    methods
        function obj = TransientSolverMassless(Mr, Kr, Cr, n_bnd, W, g0)
            % Mr, Kr, Cr : reduced matrices of the massless ROM
            % n_bnd      : number of boundary DOFs (first n_bnd of the basis)
            % W          : [n_bnd x n_c] contact directions
            % g0         : [n_c x 1] initial gaps

            n_tot = size(Kr, 1);
            obj.n_bnd = n_bnd;
            obj.n_mod = n_tot - n_bnd;

            ib = 1:n_bnd;
            ie = n_bnd+1 : n_tot;

            % ---------- validate the massless structure ----------
            if norm(Mr(ib, :), 'fro') > 1e-10 * (norm(Mr, 'fro') + eps)
                error('TSM:NotMassless', ...
                    ['The mass matrix has non-zero terms on the boundary rows ' ...
                     '(||M(bnd,:)|| = %.3e). The model is NOT massless.'], ...
                     norm(Mr(ib,:), 'fro'));
            end

            Mee = full(Mr(ie, ie));
            devI = norm(Mee - eye(obj.n_mod), 'fro');
            if devI > 1e-8 * obj.n_mod
                warning('TSM:MassNotIdentity', ...
                    ['M(inn,inn) is not the identity (dev = %.3e). The scheme assumes ' ...
                     'mass-normalized modes.'], devI);
            end

            if norm(Cr(ib, :), 'fro') > 1e-10 * (norm(Cr, 'fro') + eps)
                warning('TSM:BoundaryDamping', ...
                    ['The damping matrix has boundary terms: they will be IGNORED ' ...
                     '(the massless formulation does not account for them).']);
            end

            % ---------- partitions ----------
            obj.Kbb = full(Kr(ib, ib));
            obj.Kbe = full(Kr(ib, ie));
            obj.Kee = full(Kr(ie, ie));
            obj.Dee = full(Cr(ie, ie));

            obj.Kbb = (obj.Kbb + obj.Kbb') / 2;

            obj.W   = W;
            obj.g0  = g0(:);
            obj.n_c = size(W, 2);
            assert(size(W,1) == n_bnd,       'W must have n_bnd rows.');
            assert(numel(obj.g0) == obj.n_c, 'g0 must have n_c elements.');

            % ---------- precomputations ----------
            obj.Kbb_fact = decomposition(obj.Kbb, 'chol');

            Kbb_inv_W  = obj.Kbb_fact \ obj.W;        % inv(Kbb)*W
            obj.G_full = obj.W' * Kbb_inv_W;          % Delassus matrix
            obj.G_full = (obj.G_full + obj.G_full') / 2;

            % ---------- eps_AL ----------
            % Projected Jacobi converges for 0 < eps < 2/lambda_max(G).
            % Basing it on the SPECTRUM of G, rather than on its diagonal
            % alone, stays robust even when K_r is ill-conditioned.
            if obj.n_c == 1
                lam_max = obj.G_full;
            else
                lam_max = eigs(obj.G_full, 1, 'largestabs', ...
                               'Tolerance', 1e-6, 'MaxIterations', 500);
            end
            obj.eps_AL = 1.0 / lam_max;

            fprintf('  [massless] n_bnd = %d | n_mod = %d | n_c = %d\n', ...
                obj.n_bnd, obj.n_mod, obj.n_c);
            fprintf('  [massless] lambda_max(G) = %.3e | eps_AL = %.3e\n', ...
                lam_max, obj.eps_AL);
        end

        % =================================================================
        function dt_crit = critical_timestep(obj)
            % CRITICAL_TIMESTEP Stability limit of the explicit scheme on the
            % interior coordinates.
            %
            % The boundary does NOT contribute: it is solved quasi-statically.
            % Static condensation of the boundary changes the stiffness seen
            % by the modes:
            %   contact OPEN   : K_eff = Kee - Keb*inv(Kbb)*Kbe   (free edge)
            %   contact CLOSED : K_eff = Kee                      (locked edge)
            % The more restrictive of the two is used.

            H = obj.Kbb_fact \ obj.Kbe;                 % inv(Kbb)*Kbe
            Keff_open = obj.Kee - obj.Kbe' * H;
            Keff_open = (Keff_open + Keff_open') / 2;
            Kee_sym   = (obj.Kee + obj.Kee') / 2;

            w2 = max([ max(eig(Keff_open)), max(eig(Kee_sym)) ]);
            dt_crit = 2 / sqrt(w2);
        end

        % =================================================================
        function [t, q, lambda_hist, info] = solve(obj, tmax, dt, q0_full, qd0_full, F_handle)
            % SOLVE Integrate the massless model with set-valued contact.
            %   q0_full, qd0_full : initial conditions on the whole reduced
            %                       vector [q_b ; eta]. The boundary part of
            %                       qd0 is ignored.
            %   F_handle          : @(t) -> reduced force [n_bnd+n_mod x 1]

            nb = obj.n_bnd;  nm = obj.n_mod;
            n_steps = round(tmax/dt);
            t = (0:n_steps) * dt;

            obj.warned_AL = false;

            % ---------- stability check ----------
            dtc = obj.critical_timestep();
            fprintf('  [massless] dt = %.3e | dt_crit ~ %.3e | ratio = %.3f\n', ...
                dt, dtc, dt/dtc);
            if dt > dtc
                warning('TSM:Unstable', ...
                    ['dt = %.3e EXCEEDS the estimated stability limit %.3e. ' ...
                     'The scheme will diverge. Reduce dt or numModes.'], dt, dtc);
            end

            % ---------- operators of the explicit update ----------
            Im = eye(nm);
            A  = (1/dt)*Im + 0.5*obj.Dee;
            obj.Bmat = (1/dt)*Im - 0.5*obj.Dee;
            obj.Ainv = A \ Im;

            % ---------- storage ----------
            q           = zeros(nb + nm, n_steps+1);
            lambda_hist = zeros(obj.n_c, n_steps+1);
            n_active    = zeros(1, n_steps+1);
            iters_AL    = zeros(1, n_steps+1);
            kkt_rel     = zeros(1, n_steps+1);

            % ---------- initialization ----------
            % Leapfrog: eta on the integer grid, eta_dot on the half grid.
            % Approximation eta_dot^{1/2} = eta_dot(t0), Sec. 4.4 of the paper.
            eta  = q0_full(nb+1:end);
            etad = qd0_full(nb+1:end);         % = eta_dot^{k-1/2}
            lam  = zeros(obj.n_c, 1);

            % Diagnostic: magnitude of the projected forcing at mid-shock
            Fp = F_handle(5e-7);
            fprintf('  [debug] ||f_b|| = %.3e | ||f_e|| = %.3e\n', ...
                norm(Fp(1:nb)), norm(Fp(nb+1:end)));

            for k = 0:n_steps
                tk = k*dt;
                Fk = F_handle(tk);
                fb = Fk(1:nb);
                fe = Fk(nb+1:end);

                % ---------- 1. gap prediction (lambda = 0) ----------
                rhs    = fb - obj.Kbe * eta;
                qb_pre = obj.Kbb_fact \ rhs;              % inv(Kbb)*(fb - Kbe*eta)
                g_pre  = obj.g0 + obj.W' * qb_pre;        % = c in the paper

                Ia = find(g_pre <= 0);                    % active set
                n_active(k+1) = numel(Ia);

                % ---------- 2. contact solution ----------
                if isempty(Ia)
                    lam = zeros(obj.n_c, 1);
                    qb  = qb_pre;
                    iters_AL(k+1) = 0;
                    kkt_rel(k+1)  = 0;
                else
                    Ga = obj.G_full(Ia, Ia);
                    ca = g_pre(Ia);

                    lam_a = lam(Ia);                      % warm start

                    [lam_a, nit, res] = obj.solve_lcp(Ga, ca, lam_a);
                    iters_AL(k+1) = nit;
                    kkt_rel(k+1)  = res;

                    lam = zeros(obj.n_c, 1);
                    lam(Ia) = lam_a;

                    % q_b = inv(Kbb)*( fb - Kbe*eta + W*lambda )
                    qb = obj.Kbb_fact \ (rhs + obj.W * lam);
                end

                q(1:nb,     k+1) = qb;
                q(nb+1:end, k+1) = eta;
                lambda_hist(:,k+1) = lam;

                if ~all(isfinite(qb)) || ~all(isfinite(eta))
                    error('TSM:Diverged', ...
                        'Solution diverged at step %d (t = %.3e s).', k, tk);
                end

                if k == n_steps, break; end

                % ---------- 3. explicit update of the interior coordinates ----------
                % A*etad^{k+1/2} = fe - Kee*eta - Keb*qb + B*etad^{k-1/2}
                rhs_e = fe - obj.Kee * eta - obj.Kbe' * qb + obj.Bmat * etad;
                etad  = obj.Ainv * rhs_e;

                % ---------- 4. position update ----------
                eta = eta + etad * dt;
            end

            % ---------- diagnostics ----------
            act = n_active > 0;
            fprintf('  [debug] max|q_b| = %.3e | max|eta| = %.3e\n', ...
                max(abs(q(1:nb,:)), [], 'all'), max(abs(q(nb+1:end,:)), [], 'all'));

            obj.stats = struct( ...
                'n_active', n_active, ...
                'iters_AL', iters_AL, ...
                'kkt_rel',  kkt_rel, ...
                'dt_crit',  dtc, ...
                'eps_AL',   obj.eps_AL);
            info = obj.stats;

            fprintf('  [massless] steps with active contact: %d/%d (%.1f%%)\n', ...
                nnz(act), n_steps+1, 100*nnz(act)/(n_steps+1));
            if any(act)
                fprintf('  [massless] AL iterations: mean %.1f | max %d\n', ...
                    mean(iters_AL(act)), max(iters_AL));
                fprintf('  [massless] relative KKT residual: max %.3e (tol %.1e)\n', ...
                    max(kkt_rel), obj.tol_AL);
                fprintf('  [massless] active contacts: max %d / %d\n', ...
                    max(n_active), obj.n_c);
            end
            fprintf('  [massless] max penetration: %.3e  (gap = %.3e)\n', ...
                obj.max_penetration(q), max(obj.g0));
        end

        % =================================================================
        function pen = max_penetration(obj, q)
            % MAX_PENETRATION Largest constraint violation over the history.
            qb = q(1:obj.n_bnd, :);
            g  = obj.g0 + obj.W' * qb;
            pen = max(0, -min(g(:)));
        end
    end

    % =====================================================================
    methods (Access = private)
        function [lam, nit, res_rel] = solve_lcp(obj, G, c, lam0)
            % SOLVE_LCP Solve  0 <= (G*lam + c)  _|_  lam >= 0
            % via augmented Lagrangian and projected Jacobi (App. C):
            %   lam <- proj_{R+}( lam - eps_AL*(G*lam + c) )
            %
            % Stopping on the NORMALIZED KKT RESIDUAL:
            %   r = || min(lam, G*lam + c) ||_inf / ||c||_inf
            % A criterion on the increment of lam would be fragile: near the
            % solution the increment is dominated by round-off and an absolute
            % threshold is never reached.

            lam = max(lam0, 0);
            e   = obj.eps_AL;

            % Reference scale = magnitude of the right-hand side (predicted gap)
            scale = norm(c, inf);
            if scale < 1e-16
                scale = 1;      % floor: avoid dividing by ~0
            end

            res_rel = Inf;

            for nit = 1:obj.max_iter_AL
                r   = G*lam + c;
                res_rel = norm(min(lam, r), inf) / scale;

                if res_rel <= obj.tol_AL
                    return;
                end

                lam = max(lam - e*r, 0);
            end

            % Genuine non-convergence: warn ONCE per simulation
            if ~obj.warned_AL
                warning('TSM:ALNoConv', ...
                    ['Augmented Lagrangian: %d iterations without convergence. ' ...
                     'Relative KKT residual = %.3e (tol = %.1e). ' ...
                     'This warning is issued only once per simulation.'], ...
                    obj.max_iter_AL, res_rel, obj.tol_AL);
                obj.warned_AL = true;
            end
        end
    end
end
