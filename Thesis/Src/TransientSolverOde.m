classdef TransientSolverOde < handle
    % TRANSIENTSOLVERODE Implicit integrator (ode15s) with penalty contact.
    %
    % The contact enters through a single linear operator N, whatever the
    % model. With q the solved vector, the penetration, the penalty force and
    % its Jacobian are
    %
    %   p = N*q - g                 active where p > 0
    %   f = N(a,:)' * (k(a).*p(a))
    %   J = N(a,:)' * diag(k(a)) * N(a,:)
    %
    % and the gap g is a plain positive distance: the direction of each wall
    % is carried by the corresponding row of N. Build N with contact_operator.
    %
    % The same expressions serve the full model and every reduced one. If
    % q = q_rom and Pc maps the reduced coordinates to the physical DOFs, the
    % operator to pass is simply
    %
    %   N_rom = N * Pc
    %
    % which is why this solver no longer needs to know which model it is
    % integrating. Three consequences worth noting:
    %   - the CMS ROMs used to index their interface directly, relying on it
    %     sitting at the head of the reduced vector; N*Pc removes that
    %     assumption, so the ordering of the reduced coordinates is free;
    %   - a ROM whose interface has been reduced needs no special case, its
    %     operator being N*Pc*T_CC;
    %   - a scaled basis such as Rubin's carries its scaling inside Pc, so the
    %     gaps and the penalty stay in physical units and must NOT be
    %     rescaled on top.
    %
    % Options of solve():
    %   'ContactOperator'  [n_c x n_dofs] operator N, in the solved coordinates
    %   'ContactGap'       positive gaps, scalar or [n_c x 1]
    %   'ContactPenalty'   penalty stiffness, scalar or [n_c x 1]
    %   'MassSingular'     [] (default) detects a singular mass matrix from a
    %                      zero diagonal entry, which is what the massless
    %                      boundary of MCB and MacNeal produces, and switches
    %                      ode15s to DAE mode. Pass true or false to force it.
    %   'Label'            name of the model, used only in the printout
    %   'Eref'             reference energy -> energy-weighted AbsTol
    %   'AbsTol'           scalar AbsTol, used when 'Eref' is not supplied
    %   'RelTol'           relative tolerance (default 1e-8)
    %   'OutputTimes'      output time grid (default: internal steps)
    %
    % On AbsTol: a scalar AbsTol is not invariant with respect to the
    % reduction basis, so it applies different error criteria to different
    % ROMs. Passing 'Eref' weights the tolerance on the mechanical energy,
    % which is a physical scalar and therefore identical across bases:
    %   atol_q(i)  = epsE*sqrt(2*Eref/K_ii)
    %   atol_qd(i) = epsE*sqrt(2*Eref/M_ii)
    % Without 'Eref' the solver falls back to the scalar AbsTol.
    %
    % See also CONTACT_OPERATOR, TRANSIENTSOLVERMASSLESS.

    properties
        M % Mass matrix
        K % Stiffness matrix
        C % Damping matrix
    end

    methods
        function obj = TransientSolverOde(M, K, varargin)
            obj.M = M;
            obj.K = K;
            if nargin > 2 && ~isempty(varargin{1})
                obj.C = varargin{1};
            else
                obj.C = sparse(size(K,1), size(K,2));
            end
        end

        function [t, q_history] = solve(obj, tmax, dt, q0, qd0, F_handle, varargin)
            p = inputParser;
            addParameter(p, 'ContactOperator', []);
            addParameter(p, 'ContactGap', []);
            addParameter(p, 'ContactPenalty', []);
            addParameter(p, 'MassSingular', []);
            addParameter(p, 'Label', 'model');
            addParameter(p, 'qdd0', []);
            addParameter(p, 'Eref', []);
            addParameter(p, 'AbsTol', 1e-8);
            addParameter(p, 'RelTol', 1e-8);
            addParameter(p, 'OutputTimes', []);
            addParameter(p, 'Stats', 'on');
            parse(p, varargin{:});
            args = p.Results;

            N = args.ContactOperator;
            if isempty(N) || isempty(args.ContactGap) || isempty(args.ContactPenalty)
                error('TSO:LinearUnsupported', ...
                    ['Linear simulation is not supported: provide ' ...
                     'ContactOperator, ContactGap and ContactPenalty.']);
            end

            if isscalar(obj.C) && obj.C == 0
                obj.C = sparse(size(obj.K,1), size(obj.K,2));
            end

            n_dofs = size(obj.K, 1);
            n_c    = size(N, 1);
            y0     = [q0; qd0];

            if size(N, 2) ~= n_dofs
                error('TSO:OperatorSize', ...
                    ['ContactOperator has %d columns but the system has %d DOFs. ' ...
                     'For a reduced model the operator must be N*Pc, not N.'], ...
                    size(N, 2), n_dofs);
            end

            % --- expand gap and penalty to [n_c x 1] vectors ---
            gap = args.ContactGap(:);
            if isscalar(gap), gap = gap * ones(n_c, 1); end
            k_penalty = args.ContactPenalty(:);
            if isscalar(k_penalty), k_penalty = k_penalty * ones(n_c, 1); end

            if numel(gap) ~= n_c
                error('TSO:GapSize', 'ContactGap has %d elements, expected %d.', numel(gap), n_c);
            end
            if numel(k_penalty) ~= n_c
                error('TSO:PenaltySize', 'ContactPenalty has %d elements, expected %d.', numel(k_penalty), n_c);
            end
            if any(gap <= 0)
                error('TSO:NonPositiveGap', ...
                    ['%d gaps are not positive. The side of each wall is carried ' ...
                     'by the corresponding row of the contact operator, so the ' ...
                     'gap must be a plain distance.'], nnz(gap <= 0));
            end

            % --- mass matrix of the state-space system ---
            M_state = blkdiag(speye(n_dofs), obj.M);

            % A zero on the diagonal of M means the corresponding coordinate has
            % no inertia, as happens on the boundary of the massless ROMs, and
            % the second-order system degenerates into a DAE.
            if isempty(args.MassSingular)
                mass_singular = any(diag(obj.M) == 0);
            else
                mass_singular = logical(args.MassSingular);
            end
            if mass_singular
                sing_flag = 'yes';
                fprintf('  Singular mass matrix: ode15s runs in DAE mode\n');
            else
                sing_flag = 'no';
            end

            Nt = N';
            der_handle = @(tc, y) state_space(tc, y, obj.K, obj.C, N, Nt, ...
                gap, k_penalty, F_handle);
            jac_handle = @(tc, y) jacobian_contact(tc, y, obj.K, obj.C, N, Nt, ...
                gap, k_penalty);

            % --- tolerances ---
            reltol = args.RelTol;
            if isempty(args.Eref)
                abstol = args.AbsTol;
                fprintf('  Scalar AbsTol %.2e | RelTol %.1e\n', abstol, reltol);
            else
                epsE = 1e-2 * reltol;

                % Contact stiffness added to each coordinate, worst case with
                % every constraint active. For an operator that merely selects
                % a DOF this gives back the penalty itself.
                dK_c = full(sum((N.^2) .* k_penalty, 1)).';

                dM = abs(full(diag(obj.M)));
                dM = max(dM, 1e-6 * median(dM(dM > 0)));      % massless guard
                dK = abs(full(diag(obj.K))) + dK_c;
                dK = max(dK, (2*pi/tmax)^2 .* dM);            % floor on the slow modes

                abstol = epsE * [sqrt(2*args.Eref ./ dK); sqrt(2*args.Eref ./ dM)];
                fprintf('  Energy-weighted AbsTol [%.2e .. %.2e] | RelTol %.1e\n', ...
                    min(abstol), max(abstol), reltol);
            end

            options = odeset('RelTol', reltol, 'AbsTol', abstol, 'MaxStep', dt, ...
                'Mass', M_state, 'MassSingular', sing_flag, ...
                'Jacobian', jac_handle, 'Stats', args.Stats);

            % --- integration ---
            fprintf('Integrating %s model with analytical Jacobian (RelTol %g)...\n', ...
                upper(args.Label), reltol);
            tic;
            if isempty(args.OutputTimes)
                tspan_eval = [0 tmax];
            else
                tspan_eval = args.OutputTimes;
            end
            [t_out, y_out] = ode15s(der_handle, tspan_eval, y0, options);
            integration_time = toc;

            t = t_out';
            q_history = y_out(:, 1:n_dofs)';
            fprintf('ode15s integration time: %.2f s\n', integration_time);

            % ===============================================================
            % LOCAL FUNCTIONS: STATE SPACE AND ANALYTICAL JACOBIAN
            % ===============================================================
            function f = state_space(t_curr, y, K, C, Nop, Nop_t, g, k_pen, F_ext_handle)
                n  = size(K, 1);
                q  = y(1:n);
                qd = y(n+1:end);

                pen = Nop*q - g;
                act = pen > 0;

                if any(act)
                    F_pen = Nop_t(:, act) * (k_pen(act) .* pen(act));
                else
                    F_pen = zeros(n, 1);
                end

                f = [qd; F_ext_handle(t_curr) - K*q - C*qd - F_pen];
            end

            function J = jacobian_contact(~, y, K, C, Nop, Nop_t, g, k_pen)
                n = size(K, 1);

                pen = Nop*y(1:n) - g;
                act = pen > 0;

                if any(act)
                    K_eff = K + Nop_t(:, act) * (k_pen(act) .* Nop(act, :));
                else
                    K_eff = K;
                end

                Z = sparse(n, n);
                I = speye(n);
                J = [Z, I; -K_eff, -C];
            end
        end
    end
end
