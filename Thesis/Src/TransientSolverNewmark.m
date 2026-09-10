classdef TransientSolverNewmark < handle
%TRANSIENTSOLVERNEWMARK Second-order time integration with unilateral contact.
%
%   solver = TransientSolverNewmark(M, K, C)
%   [t, q] = solver.solve(tmax, dt, q0, qd0, F_handle, ...)
%
% Same signature as TransientSolverOde and TransientSolverMassless, so a run
% can switch integrator without touching anything else.
%
% WHY NOT CALL THE YaFEc CLASSES DIRECTLY. src/TimeIntegration/ImplicitNewmark
% and GeneralizedAlpha implement exactly these schemes, and the coefficients
% and the Prediction/Correction formulas below are copied from them so the
% results are the same to machine precision. What they cannot do is run on a
% model this size: their solution grows by concatenation, q = [q q_new], at
% every step, which on 22302 DOFs and 10^4 steps is 1.8 GB reallocated
% quadratically. They also print one line per Newton iteration. This class
% keeps only the output grid and stays silent.
%
% WHY IT IS WORTH USING. ode15s integrates the first-order form and factorises
% a 2n x 2n state-space Jacobian; here the iteration matrix is n x n. Measured
% on the 3D accelerometer, 22302 DOFs, against the stored 24 h ode15s run at
% RelTol 1e-8, over a 3 us window with the same damping:
%
%   h = 1 ns   397 s/us   error 0.0052%   peak penetration 0.04832 vs 0.04831
%   h = 2 ns   184 s/us   error 0.0349%
%   h = 5 ns    88 s/us   error 0.4427%
%
% against 864 s/us for ode15s averaged over 100 us, so between 2.2x and 9.8x
% depending on the step. The usable step depends on the damping: with the
% stiffness-proportional term of the legacy anchors the response is smooth and
% 5 ns holds, while with the TDK anchors the high frequencies survive and the
% impact count already breaks at 5 ns. Check it, do not assume it.
%
% THE STEP IS FIXED. The adaptive stepping in the YaFEc classes is not usable
% and is not reproduced: it changes h inside the Newton loop rather than
% between steps, and hmax is declared but never applied, so h can grow without
% bound. Choose h from the contact duration: the events on this model last
% about 60 ns, and a step that resolves fewer than roughly twelve points
% across one starts merging impacts.
%
% See also TRANSIENTSOLVERODE, CONTACT_RESIDUAL, CONTACT_TANGENT.

    properties
        M, K, C
    end

    methods
        function obj = TransientSolverNewmark(M, K, C)
            obj.M = M;  obj.K = K;  obj.C = C;
        end

        function [t, q_history] = solve(obj, tmax, dt, q0, qd0, F_handle, varargin)
            p = inputParser;
            addParameter(p, 'ContactOperator', []);
            addParameter(p, 'ContactGap', []);
            addParameter(p, 'ContactPenalty', []);
            addParameter(p, 'ContactCompliance', []);
            addParameter(p, 'Scheme', 'genalpha');   % 'genalpha' | 'newmark'
            addParameter(p, 'RhoInf', 0.7);          % genalpha: 1 = no numerical damping
            addParameter(p, 'Alpha', 0);             % newmark: 0 = trapezoidal, conservative
            addParameter(p, 'TimeStep', []);         % h; defaults to dt
            addParameter(p, 'RelTol', 1e-6);         % on the Newton residual, NOT on the
                                                     % discretisation error: convergence
                                                     % here says the step is solved, not
                                                     % that h is small enough
            addParameter(p, 'MaxNRit', 25);
            addParameter(p, 'OutputTimes', []);
            addParameter(p, 'Label', 'model');
            addParameter(p, 'qdd0', []);
            parse(p, varargin{:});
            args = p.Results;

            n = size(obj.K, 1);
            N   = args.ContactOperator;
            gap = args.ContactGap(:);
            k   = args.ContactPenalty;
            if isempty(N) || isempty(gap) || isempty(k)
                error('TSN:NoContact', ...
                    'ContactOperator, ContactGap and ContactPenalty are required.');
            end
            if isscalar(k), k = k * ones(numel(gap), 1); end

            h = args.TimeStep;
            if isempty(h), h = dt; end

            % --- output grid ---------------------------------------------
            % The comparison contract of the post-processing is that every
            % model lands on the same instants. A step that does not divide
            % the output spacing would need interpolation, which would quietly
            % smooth exactly the impacts this is meant to resolve.
            if isempty(args.OutputTimes)
                t = 0 : h : tmax;
            else
                t = args.OutputTimes(:)';
                spacing = t(2) - t(1);
                ratio   = spacing / h;
                if abs(ratio - round(ratio)) > 1e-9
                    error('TSN:GridMismatch', ...
                        ['The time step %.4g s does not divide the output spacing ' ...
                         '%.4g s (ratio %.6f). Choose h so that it does: the ' ...
                         'alternative is interpolating the impacts away.'], ...
                        h, spacing, ratio);
                end
            end
            n_out = numel(t);
            q_history = zeros(n, n_out);

            % --- scheme coefficients -------------------------------------
            % Copied from src/TimeIntegration/GeneralizedAlpha.m (Arnold &
            % Bruels, optimal parameters of Chung & Hulbert) and from
            % ImplicitNewmark.m. Equivalence with both is checked in the
            % verification script rather than assumed.
            switch lower(args.Scheme)
                case 'genalpha'
                    ri = args.RhoInf;
                    if ~(ri > 0 && ri <= 1)
                        error('TSN:BadRhoInf', 'RhoInf must be in (0, 1], got %g.', ri);
                    end
                    alpha_m = (2*ri - 1) / (ri + 1);
                    alpha_f = ri / (ri + 1);
                    gam     = 0.5 + alpha_f - alpha_m;
                    bet     = 0.25 * (gam + 0.5)^2;
                    beta_p  = h^2 * bet * (1 - alpha_f) / (1 - alpha_m);
                    gamma_p = gam * h * (1 - alpha_f) / (1 - alpha_m);
                    use_aux = true;
                case 'newmark'
                    al      = args.Alpha;
                    bet     = (1 + al)^2 / 4;
                    gam     = 0.5 + al;
                    beta_p  = bet * h^2;
                    gamma_p = gam * h;
                    alpha_m = 0; alpha_f = 0;
                    use_aux = false;
                otherwise
                    error('TSN:BadScheme', ...
                        'Scheme must be ''genalpha'' or ''newmark'', got ''%s''.', ...
                        args.Scheme);
            end

            % --- initial acceleration ------------------------------------
            % M*qdd0 = F(0) - C*qd0 - K*q0 - F_contact(q0). At rest with a
            % forcing that starts from zero this is zero, but a run that starts
            % already loaded or already touching would need the real thing, so
            % it is solved rather than assumed.
            if isempty(args.qdd0)
                Nt = N';
                F_pen0 = contact_tangent(q0, obj.K, N, Nt, gap, k, args.ContactCompliance);
                rhs = F_handle(0) - obj.C*qd0 - obj.K*q0 - F_pen0;
                if norm(rhs) == 0
                    % Starting from rest under a forcing that begins at zero.
                    % Solving anyway warns about the conditioning of M, which
                    % here is scaling and not rank: in um/MPa/kg the nodal
                    % masses are around 1e-14, so rcond is tiny by construction.
                    qdd0 = zeros(n, 1);
                else
                    qdd0 = obj.M \ rhs;
                end
            else
                qdd0 = args.qdd0;
            end

            % --- time loop -----------------------------------------------
            Res = contact_residual(obj.M, obj.K, obj.C, N, gap, k, F_handle, ...
                                   args.ContactCompliance);
            q_old = q0;  qd_old = qd0;  qdd_old = qdd0;  a = qdd0;
            q_history(:, 1) = q0;
            i_out = 2;                       % next output instant to fill
            n_step = round(tmax / h);
            nr_tot = 0;  n_slow = 0;
            fprintf('Integrating %s with %s (h = %.4g s, %d steps)...\n', ...
                upper(args.Label), lower(args.Scheme), h, n_step);
            tic;

            for s = 1:n_step
                tc = s * h;

                % Prediction
                if use_aux
                    qd_new  = qd_old + h*(1 - gam)*a;
                    q_new   = q_old + h*qd_old + (0.5 - bet)*h^2*a;
                    a_new   = (alpha_f*qdd_old - alpha_m*a) / (1 - alpha_m);
                    q_new   = q_new  + h^2*bet*a_new;
                    qd_new  = qd_new + h*gam*a_new;
                else
                    qd_new  = qd_old + h*(1 - gam)*qdd_old;
                    q_new   = q_old + h*qd_old + (0.5 - bet)*h^2*qdd_old;
                end
                qdd_new = zeros(n, 1);

                % Newton-Raphson
                it = 0;
                while true
                    [r, Mt, Ct, Kt, c0] = Res(q_new, qd_new, qdd_new, tc);
                    if norm(r)/c0 < args.RelTol, break; end
                    it = it + 1;
                    if it > args.MaxNRit
                        n_slow = n_slow + 1;
                        break
                    end
                    S  = Mt + gamma_p*Ct + beta_p*Kt;
                    Da = -(S \ r);
                    q_new   = q_new   + beta_p*Da;
                    qd_new  = qd_new  + gamma_p*Da;
                    qdd_new = qdd_new + Da;
                end
                nr_tot = nr_tot + it;

                if use_aux
                    a = a_new + qdd_new*(1 - alpha_f)/(1 - alpha_m);
                end
                q_old = q_new;  qd_old = qd_new;  qdd_old = qdd_new;

                % Store only where the output grid asks for it.
                if i_out <= n_out && abs(tc - t(i_out)) < 0.5*h
                    q_history(:, i_out) = q_new;
                    i_out = i_out + 1;
                end
            end

            integration_time = toc;
            fprintf('%s integration time: %.2f s | %.2f Newton it/step\n', ...
                lower(args.Scheme), integration_time, nr_tot/max(n_step,1));
            if n_slow > 0
                % Not fatal, but it means the step was solved to a looser
                % residual than asked on those steps, so it is said out loud
                % rather than hidden in a return value nobody reads.
                warning('TSN:MaxNRit', ...
                    ['%d of %d steps hit the %d iteration cap without reaching ' ...
                     'RelTol %.1e. Reduce the time step.'], ...
                    n_slow, n_step, args.MaxNRit, args.RelTol);
            end
            if i_out <= n_out
                warning('TSN:ShortRun', ...
                    ['Only %d of %d output instants were filled: tmax is not an ' ...
                     'integer number of steps.'], i_out-1, n_out);
                t = t(1:i_out-1);
                q_history = q_history(:, 1:i_out-1);
            end
        end
    end
end
