classdef TransientSolverOde < handle
    % TRANSIENTSOLVERODE Integratore implicito (ode15s) con contatto a penalita'.
    %
    % Unifica le due versioni precedenti (TransientSolverOde e
    % TransientSolverOde_V2). Gestisce un numero qualsiasi di GdL di contatto,
    % ognuno con il proprio gap FIRMATO:
    %
    %   gap > 0  ->  muro nella direzione positiva del GdL, penetra se q > gap
    %   gap < 0  ->  muro nella direzione negativa del GdL, penetra se q < gap
    %
    % Il test e' sign(gap)*(q - gap) > 0, che con gap scalare positivo si
    % riduce a q > gap: il comportamento della versione a interfaccia singola
    % e' quindi contenuto in questo come caso particolare.
    %
    % Parametri di solve():
    %   'ContactTargetDOF'  GdL di contatto (indici nel vettore risolto)
    %   'ContactGap'        gap firmato, scalare o vettore [n_c x 1]
    %   'ContactPenalty'    rigidezza di penalita', scalare o vettore
    %   'ModelType'         'FOM' | 'MC' | 'Rubin' | 'MCB' | 'MN'
    %   'ProjectionMatrix'  richiesta solo da 'MC'
    %   'Eref'              energia di riferimento -> AbsTol pesata in energia
    %   'AbsTol'            AbsTol scalare, usata se 'Eref' non e' fornita
    %   'RelTol'            tolleranza relativa (default 1e-8)
    %   'OutputTimes'       griglia temporale di output (default: passi interni)
    %
    % AbsTol: una AbsTol scalare non e' invariante rispetto alla base di
    % riduzione, quindi applica criteri d'errore diversi a ROM diversi.
    % Passando 'Eref' la tolleranza viene pesata sull'energia meccanica, che e'
    % uno scalare fisico e quindi identico per tutte le basi:
    %   atol_q(i)  = epsE*sqrt(2*Eref/K_ii)
    %   atol_qd(i) = epsE*sqrt(2*Eref/M_ii)
    % Senza 'Eref' si ricade sulla AbsTol scalare (comportamento storico).

    properties
        M % Matrice di massa
        K % Matrice di rigidezza
        C % Matrice di smorzamento
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
                    'Simulazione lineare non supportata: fornire ContactTargetDOF, ContactGap e ContactPenalty.');
            end

            if isscalar(obj.C) && obj.C == 0
                obj.C = sparse(size(obj.K,1), size(obj.K,2));
            end

            target_dofs = args.ContactTargetDOF(:);
            n_c         = numel(target_dofs);
            n_dofs      = size(obj.K, 1);
            y0          = [q0; qd0];

            % --- gap e penalita' espansi a vettore [n_c x 1] ---
            gap_wall = args.ContactGap(:);
            if isscalar(gap_wall), gap_wall = gap_wall * ones(n_c, 1); end
            k_penalty = args.ContactPenalty(:);
            if isscalar(k_penalty), k_penalty = k_penalty * ones(n_c, 1); end

            if numel(gap_wall) ~= n_c
                error('TSO:GapSize', 'ContactGap ha %d elementi, attesi %d.', numel(gap_wall), n_c);
            end
            if numel(k_penalty) ~= n_c
                error('TSO:PenaltySize', 'ContactPenalty ha %d elementi, attesi %d.', numel(k_penalty), n_c);
            end
            if any(gap_wall == 0)
                error('TSO:ZeroGap', ...
                    ['Almeno un gap e'' nullo: il segno del gap definisce da che parte ' ...
                     'sta il muro, quindi un gap zero e'' ambiguo.']);
            end

            % --- matrice di massa del sistema in forma di stato ---
            M_state = blkdiag(speye(n_dofs), obj.M);

            % --- setup specifico del modello ---
            mass_singular = 'no';
            Pc_contact = [];

            switch upper(args.ModelType)
                case {'FOM', 'RUBIN', 'MCB', 'MN'}
                    % L'interfaccia e' in testa al vettore ridotto per i ROM CMS,
                    % quindi i GdL di contatto sono indici diretti.
                    if any(strcmpi(args.ModelType, {'MCB', 'MN'}))
                        % Massa nulla al boundary -> ode15s in modalita' DAE.
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
                            'Il metodo MC richiede ProjectionMatrix (Pc).');
                    end
                    Pc_contact = Pc(target_dofs, :);
                    der_handle = @(tc, y) state_space_projected(tc, y, obj.K, obj.C, ...
                        Pc_contact, gap_wall, k_penalty, F_handle);
                    jac_handle = @(tc, y) jacobian_projected(tc, y, obj.K, obj.C, ...
                        Pc_contact, gap_wall, k_penalty);

                otherwise
                    error('TSO:BadModelType', ...
                        'ModelType non riconosciuto: usare FOM, MC, Rubin, MCB o MN.');
            end

            % --- tolleranze ---
            reltol = args.RelTol;
            if isempty(args.Eref)
                abstol = args.AbsTol;
                fprintf('  AbsTol scalare %.2e | RelTol %.1e\n', abstol, reltol);
            else
                epsE = 1e-2 * reltol;

                % rigidezza di contatto, caso peggiore (tutti i GdL attivi)
                dK_c = zeros(n_dofs, 1);
                if strcmpi(args.ModelType, 'MC')
                    dK_c = full(sum((Pc_contact.^2) .* k_penalty, 1)).';
                else
                    dK_c(target_dofs) = k_penalty;
                end

                dM = abs(full(diag(obj.M)));
                dM = max(dM, 1e-6 * median(dM(dM > 0)));      % guardia massless
                dK = abs(full(diag(obj.K))) + dK_c;
                dK = max(dK, (2*pi/tmax)^2 .* dM);            % floor sui modi lenti

                abstol = epsE * [sqrt(2*args.Eref ./ dK); sqrt(2*args.Eref ./ dM)];
                fprintf('  AbsTol energetica [%.2e .. %.2e] | RelTol %.1e\n', ...
                    min(abstol), max(abstol), reltol);
            end

            options = odeset('RelTol', reltol, 'AbsTol', abstol, 'MaxStep', dt, ...
                'Mass', M_state, 'MassSingular', mass_singular, ...
                'Jacobian', jac_handle, 'Stats', args.Stats);

            % --- integrazione ---
            fprintf('Integrazione modello %s con Jacobiano analitico (RelTol %g)...\n', ...
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
            fprintf('Tempo di integrazione ode15s: %.2f s\n', integration_time);

            % ===============================================================
            % FUNZIONI LOCALI: SPAZIO DI STATO E JACOBIANI ANALITICI
            % ===============================================================

            % --- 1. FOM, RUBIN, MCB, MN (GdL fisici diretti) ---
            function f = state_space_standard(t_curr, y, K, C, contact_dofs, gap, k_pen, F_ext_handle)
                n  = size(K, 1);
                q  = y(1:n);
                qd = y(n+1:end);

                % Compenetrazione col segno
                penetration = q(contact_dofs) - gap;

                % Attivo solo se il GdL supera il muro nella direzione giusta
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

            % --- 2. MC (GdL di contatto proiettati) ---
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
