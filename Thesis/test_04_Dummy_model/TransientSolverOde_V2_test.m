classdef TransientSolverOde_V2_test < handle
    % TRANSIENTSOLVERODE_V2_TEST
    %   Variante di test di TransientSolverOde_V2 con due sole aggiunte:
    %
    %   1) ContactPenalty puo' essere VETTORIALE (una rigidezza per DOF di
    %      contatto). Serve per la base di Rubin scalata sul ramo standard:
    %      con q_fisico_i = y_i / d_i la penetrazione e la forza generalizzata
    %      si riscalano come
    %           gap'_i = d_i * gap_i        k'_i = k_pen / d_i^2
    %      il che permette di conservare la risoluzione del contatto O(b)
    %      invece di passare al ramo proiettato, che costa O(b*r).
    %      Retro-compatibile: uno scalare viene espanso automaticamente.
    %
    %   2) 'AbsTolOverride' per imporre una AbsTol scalare identica a tutte
    %      le configurazioni. La AbsTol energetica usa diag(Mr) e diag(Kr),
    %      che NON sono invarianti per cambio base: senza override, basi
    %      diverse ricevono criteri d'errore diversi.
    %
    %   Il resto (jacobiano analitico, gestione DAE per MCB/MN, tolleranza
    %   energetica) e' invariato rispetto alla V2.

    properties
        M % Mass matrix
        K % Stiffness matrix
        C % Damping matrix
    end

    methods
        function obj = TransientSolverOde_V2_test(M, K, varargin)
            obj.M = M;
            obj.K = K;
            if nargin > 2 && ~isempty(varargin{1})
                obj.C = varargin{1};
            else
                obj.C = sparse(size(K,1), size(K,2));
            end
        end

        function [t, q_history, stats] = solve(obj, tmax, dt, q0, qd0, F_handle, varargin)
            p = inputParser;
            addParameter(p, 'ContactTargetDOF', []);
            addParameter(p, 'ContactGap', []);
            addParameter(p, 'ContactPenalty', []);
            addParameter(p, 'ProjectionMatrix', []);
            addParameter(p, 'qdd0', []);
            addParameter(p, 'ModelType', 'FOM');
            addParameter(p, 'Eref', []);
            addParameter(p, 'RelTol', 1e-8);
            addParameter(p, 'OutputTimes', []);
            addParameter(p, 'AbsTolOverride', []);      % <-- NUOVO
            parse(p, varargin{:});
            args = p.Results;

            is_nonlinear = ~isempty(args.ContactTargetDOF) && ...
                           ~isempty(args.ContactGap) && ~isempty(args.ContactPenalty);

            if isscalar(obj.C) && obj.C == 0
                obj.C = sparse(size(obj.K,1), size(obj.K,2));
            end

            if ~is_nonlinear
                error('Linear simulation currently unsupported in this version.');
            end

            target_dofs = args.ContactTargetDOF(:);
            gap_wall    = args.ContactGap(:);
            n_bnd_loc   = numel(target_dofs);
            n_dofs      = size(obj.K, 1);
            y0          = [q0; qd0];

            % --- penalita' vettoriale (retro-compatibile con scalare) -----
            k_penalty = args.ContactPenalty(:);
            if isscalar(k_penalty)
                k_penalty = k_penalty * ones(n_bnd_loc, 1);
            end
            assert(numel(k_penalty) == n_bnd_loc, ...
                'ContactPenalty deve essere scalare o lungo quanto ContactTargetDOF.');
            assert(numel(gap_wall) == n_bnd_loc, ...
                'ContactGap deve essere lungo quanto ContactTargetDOF.');

            M_state = blkdiag(speye(n_dofs), obj.M);
            mass_singular = 'no';

            % ==============================================================
            % SETUP SPECIFICO DEL MODELLO E JACOBIANO ANALITICO
            % ==============================================================
            if any(strcmpi(args.ModelType, {'FOM','Rubin','MCB','MN'}))
                if any(strcmpi(args.ModelType, {'MCB','MN'}))
                    mass_singular = 'yes';   % M_bb = 0 all'interfaccia
                end
                der_handle = @(t_curr, y) state_space_standard(t_curr, y, obj.K, obj.C, ...
                                target_dofs, gap_wall, k_penalty, F_handle);
                jac_handle = @(t_curr, y) jacobian_standard(t_curr, y, obj.K, obj.C, ...
                                target_dofs, gap_wall, k_penalty);
                Pc_contact = [];

            elseif strcmpi(args.ModelType, 'MC')
                Pc = args.ProjectionMatrix;
                if isempty(Pc)
                    error('For MC method, you must provide the ProjectionMatrix (Pc).');
                end
                Pc_contact = Pc(target_dofs, :);
                der_handle = @(t_curr, y) state_space_projected(t_curr, y, obj.K, obj.C, ...
                                Pc_contact, gap_wall, k_penalty, F_handle);
                jac_handle = @(t_curr, y) jacobian_projected(t_curr, y, obj.K, obj.C, ...
                                Pc_contact, gap_wall, k_penalty);
            else
                error('ModelType not recognized. Use ''FOM'', ''MC'', ''Rubin'', ''MCB'', or ''MN''.');
            end

            % ============ TOLLERANZE ============
            adaptive_reltol = args.RelTol;

            if ~isempty(args.AbsTolOverride)
                adaptive_abstol = args.AbsTolOverride;
                fprintf('  AbsTol OVERRIDE = %.2e | RelTol %.1e\n', ...
                        min(adaptive_abstol), adaptive_reltol);
            else
                if isempty(args.Eref)
                    error('Parametro ''Eref'' mancante: calcolarlo una volta dal FOM.');
                end
                epsE = 1e-2 * adaptive_reltol;

                % rigidezza di contatto, caso peggiore (tutti i DOF attivi)
                dK_c = zeros(n_dofs, 1);
                if strcmpi(args.ModelType, 'MC')
                    % somma pesata: diag(Pc_c' * diag(k) * Pc_c)
                    dK_c = full(sum((Pc_contact.^2) .* k_penalty, 1)).';
                else
                    dK_c(target_dofs) = k_penalty;
                end

                dM = abs(full(diag(obj.M)));
                dM = max(dM, 1e-6 * median(dM(dM > 0)));       % guardia massless
                dK = abs(full(diag(obj.K))) + dK_c;
                dK = max(dK, (2*pi/tmax)^2 .* dM);             % floor sui modi lenti

                adaptive_abstol = epsE * [sqrt(2*args.Eref ./ dK); ...
                                          sqrt(2*args.Eref ./ dM)];

                fprintf('  AbsTol energetica [%.2e .. %.2e] | RelTol %.1e\n', ...
                        min(adaptive_abstol), max(adaptive_abstol), adaptive_reltol);
            end

            % ==============================================================
            % OPZIONI E INTEGRAZIONE
            % ==============================================================
            options = odeset('RelTol', adaptive_reltol, 'AbsTol', adaptive_abstol, ...
                'MaxStep', dt, 'Mass', M_state, 'MassSingular', mass_singular, ...
                'Jacobian', jac_handle, 'Stats', 'on');

            fprintf('Integrating %s model with analytical Jacobian (RelTol: %g)...\n', ...
                    upper(args.ModelType), adaptive_reltol);
            tic;
            if isempty(args.OutputTimes)
                tspan_eval = [0 tmax];
            else
                tspan_eval = args.OutputTimes;
            end
            sol = ode15s(der_handle, tspan_eval, y0, options);
            integration_time = toc;

            if isempty(args.OutputTimes)
                t_out = sol.x(:);
                y_out = sol.y.';
            else
                t_out = args.OutputTimes(:);
                y_out = deval(sol, t_out).';
            end

            t = t_out';
            q_history = y_out(:, 1:n_dofs)';
            stats = sol.stats;
            fprintf('ode15s integration time: %.2f s\n', integration_time);

            % ===============================================================
            % LOCAL FUNCTIONS: STATE-SPACE & ANALYTICAL JACOBIANS
            % ===============================================================

            % --- 1. FOM, RUBIN, MCB, MN (coordinate d'interfaccia dirette) ---
            function f = state_space_standard(t_curr, y, K, C, contact_dofs, gap, k_pen, F_ext_handle)
                n  = size(K, 1);
                q  = y(1:n);
                qd = y(n+1:end);

                penetration = q(contact_dofs) - gap;
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

            % --- 2. MC e qualunque base con interfaccia proiettata ---------
            function f = state_space_projected(t_curr, y, K, C, Pc_c, gap, k_pen, F_ext_handle)
                n  = size(K, 1);
                q  = y(1:n);
                qd = y(n+1:end);

                penetration = Pc_c * q - gap;
                is_pen = (sign(gap) .* penetration) > 0;

                F_pen_rom = zeros(n, 1);
                if any(is_pen)
                    F_pen_contact = k_pen(is_pen) .* penetration(is_pen);
                    F_pen_rom = Pc_c(is_pen, :)' * F_pen_contact;
                end

                f = [qd; F_ext_handle(t_curr) - K*q - C*qd - F_pen_rom];
            end

            function J = jacobian_projected(~, y, K, C, Pc_c, gap, k_pen)
                n = size(K, 1);
                q = y(1:n);

                penetration = Pc_c * q - gap;
                is_pen = (sign(gap) .* penetration) > 0;

                if any(is_pen)
                    Pc_active = Pc_c(is_pen, :);
                    % Pc_a' * diag(k) * Pc_a
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
