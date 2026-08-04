classdef TransientSolverOde < handle
    % TRANSIENTSOLVERODE Implicit integrator (ode15s) with penalty contact.
    %
    % Handles any number of contact DOFs, each with its own SIGNED gap:
    %
    %   gap > 0  ->  wall along the positive direction of the DOF, penetrates if q > gap
    %   gap < 0  ->  wall along the negative direction of the DOF, penetrates if q < gap
    %
    % The activation test is sign(gap)*(q - gap) > 0, which reduces to q > gap
    % for a positive scalar gap. The single-interface behaviour is therefore a
    % particular case of this one, and no separate solver is needed.
    %
    % Options of solve():
    %   'ContactTargetDOF'  contact DOFs (indices in the solved vector)
    %   'ContactGap'        signed gap, scalar or vector [n_c x 1]
    %   'ContactPenalty'    penalty stiffness, scalar or vector
    %   'ModelType'         'FOM' | 'MC' | 'Rubin' | 'MCB' | 'MN'
    %   'ProjectionMatrix'  required by 'MC' only
    %   'Eref'              reference energy -> energy-weighted AbsTol
    %   'AbsTol'            scalar AbsTol, used when 'Eref' is not supplied
    %   'RelTol'            relative tolerance (default 1e-8)
    %   'OutputTimes'       output time grid (default: internal steps)
    %
    % On AbsTol: a scalar AbsTol is not invariant with respect to the
    % reduction basis, so it applies different error criteria to different
    % ROMs. Passing 'Eref' weights the tolerance on the mechanical energy,
    % which is a physical scalar and therefore identical across bases:
    %   atol_q(i)  = epsE*sqrt(2*Eref/K_ii)
    %   atol_qd(i) = epsE*sqrt(2*Eref/M_ii)
    % Without 'Eref' the solver falls back to the scalar AbsTol.

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
            addParameter(p, 'ContactTargetDOF', []);
            addParameter(p, 'ContactGap', []);
            addParameter(p, 'ContactPenalty', []);
            addParameter(p, 'ProjectionMatrix', []);
            addParameter(p, 'qdd0', []);
            addParameter(p, 'ModelType', 'FOM');
            addParameter(p, 'Eref', []);
            addParameter(p, 'AbsTol', 1e-8);
            addParameter(p, 'RelTol', 1e-8);
            addParameter(p, 'OutputTimes', []);
            addParameter(p, 'Stats', 'on');
            parse(p, varargin{:});
            args = p.Results;

            is_nonlinear = ~isempty(args.ContactTargetDOF) && ...
                           ~isempty(args.ContactGap) && ~isempty(args.ContactPenalty);
            if ~is_nonlinear
                error('TSO:LinearUnsupported', ...
                    'Linear simulation is not supported: provide ContactTargetDOF, ContactGap and ContactPenalty.');
            end

            if isscalar(obj.C) && obj.C == 0
                obj.C = sparse(size(obj.K,1), size(obj.K,2));
            end

            target_dofs = args.ContactTargetDOF(:);
            n_c         = numel(target_dofs);
            n_dofs      = size(obj.K, 1);
            y0          = [q0; qd0];

            % --- expand gap and penalty to [n_c x 1] vectors ---
            gap_wall = args.ContactGap(:);
            if isscalar(gap_wall), gap_wall = gap_wall * ones(n_c, 1); end
            k_penalty = args.ContactPenalty(:);
            if isscalar(k_penalty), k_penalty = k_penalty * ones(n_c, 1); end

            if numel(gap_wall) ~= n_c
                error('TSO:GapSize', 'ContactGap has %d elements, expected %d.', numel(gap_wall), n_c);
            end
            if numel(k_penalty) ~= n_c
                error('TSO:PenaltySize', 'ContactPenalty has %d elements, expected %d.', numel(k_penalty), n_c);
            end
            if any(gap_wall == 0)
                error('TSO:ZeroGap', ...
                    ['At least one gap is zero. The sign of the gap tells which side ' ...
                     'the wall is on, so a zero gap is ambiguous.']);
            end

            % --- mass matrix of the state-space system ---
            M_state = blkdiag(speye(n_dofs), obj.M);

            % --- model-specific setup ---
            mass_singular = 'no';
            Pc_contact = [];

            switch upper(args.ModelType)
                case {'FOM', 'RUBIN', 'MCB', 'MN'}
                    % The interface sits at the head of the reduced vector for
                    % the CMS ROMs, so the contact DOFs are direct indices.
                    if any(strcmpi(args.ModelType, {'MCB', 'MN'}))
                        % Zero mass at the boundary -> run ode15s in DAE mode.
                        mass_singular = 'yes';
                    end
                    der_handle = @(tc, y) state_space_standard(tc, y, obj.K, obj.C, ...
                        target_dofs, gap_wall, k_penalty, F_handle);
                    jac_handle = @(tc, y) jacobian_standard(tc, y, obj.K, obj.C, ...
                        target_dofs, gap_wall, k_penalty);

                case 'MC'
                    Pc = args.ProjectionMatrix;
                    if isempty(Pc)
                        error('TSO:NoProjection', ...
                            'The MC method requires ProjectionMatrix (Pc).');
                    end
                    Pc_contact = Pc(target_dofs, :);
                    der_handle = @(tc, y) state_space_projected(tc, y, obj.K, obj.C, ...
                        Pc_contact, gap_wall, k_penalty, F_handle);
                    jac_handle = @(tc, y) jacobian_projected(tc, y, obj.K, obj.C, ...
                        Pc_contact, gap_wall, k_penalty);

                otherwise
                    error('TSO:BadModelType', ...
                        'Unrecognized ModelType: use FOM, MC, Rubin, MCB or MN.');
            end

            % --- tolerances ---
            reltol = args.RelTol;
            if isempty(args.Eref)
                abstol = args.AbsTol;
                fprintf('  Scalar AbsTol %.2e | RelTol %.1e\n', abstol, reltol);
            else
                epsE = 1e-2 * reltol;

                % Contact stiffness, worst case with every DOF active
                dK_c = zeros(n_dofs, 1);
                if strcmpi(args.ModelType, 'MC')
                    dK_c = full(sum((Pc_contact.^2) .* k_penalty, 1)).';
                else
                    dK_c(target_dofs) = k_penalty;
                end

                dM = abs(full(diag(obj.M)));
                dM = max(dM, 1e-6 * median(dM(dM > 0)));      % massless guard
                dK = abs(full(diag(obj.K))) + dK_c;
                dK = max(dK, (2*pi/tmax)^2 .* dM);            % floor on the slow modes

                abstol = epsE * [sqrt(2*args.Eref ./ dK); sqrt(2*args.Eref ./ dM)];
                fprintf('  Energy-weighted AbsTol [%.2e .. %.2e] | RelTol %.1e\n', ...
                    min(abstol), max(abstol), reltol);
            end

            options = odeset('RelTol', reltol, 'AbsTol', abstol, 'MaxStep', dt, ...
                'Mass', M_state, 'MassSingular', mass_singular, ...
                'Jacobian', jac_handle, 'Stats', args.Stats);

            % --- integration ---
            fprintf('Integrating %s model with analytical Jacobian (RelTol %g)...\n', ...
                upper(args.ModelType), reltol);
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
            % LOCAL FUNCTIONS: STATE SPACE AND ANALYTICAL JACOBIANS
            % ===============================================================

            % --- 1. FOM, RUBIN, MCB, MN (direct physical DOFs) ---
            function f = state_space_standard(t_curr, y, K, C, contact_dofs, gap, k_pen, F_ext_handle)
                n  = size(K, 1);
                q  = y(1:n);
                qd = y(n+1:end);

                % Signed penetration
                penetration = q(contact_dofs) - gap;

                % Active only if the DOF crossed the wall in the right direction
                is_pen = (sign(gap) .* penetration) > 0;

                F_pen = zeros(n, 1);
                if any(is_pen)
                    F_pen(contact_dofs(is_pen)) = k_pen(is_pen) .* penetration(is_pen);
                end

                f = [qd; F_ext_handle(t_curr) - K*q - C*qd - F_pen];
            end

            function J = jacobian_standard(~, y, K, C, contact_dofs, gap, k_pen)
                n = size(K, 1);
                q = y(1:n);

                penetration = q(contact_dofs) - gap;
                is_pen = (sign(gap) .* penetration) > 0;

                if any(is_pen)
                    active_dofs = contact_dofs(is_pen);
                    K_pen = sparse(active_dofs, active_dofs, k_pen(is_pen), n, n);
                    K_eff = K + K_pen;
                else
                    K_eff = K;
                end

                Z = sparse(n, n);
                I = speye(n);
                J = [Z, I; -K_eff, -C];
            end

            % --- 2. MC (projected contact DOFs) ---
            function f = state_space_projected(t_curr, y, K, C, Pcc, gap, k_pen, F_ext_handle)
                n  = size(K, 1);
                q  = y(1:n);
                qd = y(n+1:end);

                penetration = Pcc * q - gap;
                is_pen = (sign(gap) .* penetration) > 0;

                F_pen_rom = zeros(n, 1);
                if any(is_pen)
                    F_pen_contact = k_pen(is_pen) .* penetration(is_pen);
                    Pc_active = Pcc(is_pen, :);
                    F_pen_rom = Pc_active' * F_pen_contact;
                end

                f = [qd; F_ext_handle(t_curr) - K*q - C*qd - F_pen_rom];
            end

            function J = jacobian_projected(~, y, K, C, Pcc, gap, k_pen)
                n = size(K, 1);
                q = y(1:n);

                penetration = Pcc * q - gap;
                is_pen = (sign(gap) .* penetration) > 0;

                if any(is_pen)
                    Pc_active = Pcc(is_pen, :);
                    K_pen_rom = Pc_active' * (k_pen(is_pen) .* Pc_active);
                    K_eff = K + K_pen_rom;
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
