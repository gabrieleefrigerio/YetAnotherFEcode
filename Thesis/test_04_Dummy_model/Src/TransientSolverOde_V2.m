classdef TransientSolverOde_V2 < handle
    properties
        M % Mass matrix
        K % Stiffness matrix
        C % Damping matrix
    end
    
    methods
        function obj = TransientSolverOde_V2(M, K, varargin)
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
            addParameter(p, 'ModelType', 'FOM'); % Model selection: 'FOM', 'MC', 'Rubin', 'MCB', 'MN'
            addParameter(p, 'Eref', []);
            addParameter(p, 'RelTol', 1e-8);
            addParameter(p, 'OutputTimes', []);
            parse(p, varargin{:});
            args = p.Results;
            
            is_nonlinear = ~isempty(args.ContactTargetDOF) && ~isempty(args.ContactGap) && ~isempty(args.ContactPenalty);
            if isscalar(obj.C) && obj.C == 0
                obj.C = sparse(size(obj.K,1), size(obj.K,2));
            end
            
            if is_nonlinear
                target_dofs = args.ContactTargetDOF;
                gap_wall = args.ContactGap;
                k_penalty = args.ContactPenalty;
                n_dofs = size(obj.K, 1);
                y0 = [q0; qd0];
                if isscalar(k_penalty), k_penalty = k_penalty * ones(numel(target_dofs),1); end
                % ==============================================================
                % MASS MATRIX per LHS del sistema ode15s
                % ==============================================================
                M_state = blkdiag(speye(n_dofs), obj.M);
                
                % ==============================================================
                % SETUP SPECIFICO DEL MODELLO E JACOBIANO ANALITICO
                % ==============================================================
                mass_singular = 'no'; % Default: ode normale
                
                if strcmpi(args.ModelType, 'FOM') || strcmpi(args.ModelType, 'Rubin') || strcmpi(args.ModelType, 'MCB') || strcmpi(args.ModelType, 'MN')
                    % Tolleranze
                    
                    % Per i ROM basati su CMS, tolleranza assoluta più stringente
                    if strcmpi(args.ModelType, 'Rubin')
                    end
                    
                    % Se è MCB o MN (MacNeal), forziamo il solver in modalità DAE algebrica 
                    % poiché il blocco di massa all'interfaccia è nullo (M_bb = 0)
                    if strcmpi(args.ModelType, 'MCB') || strcmpi(args.ModelType, 'MN')
                        mass_singular = 'yes';
                    end
                    
                    % Assegnazione Handle (Nodi Standard poiché l'interfaccia è in testa per i ROM CMS)
                    der_handle = @(t_curr, y) state_space_standard(t_curr, y, obj.K, obj.C, target_dofs, gap_wall, k_penalty, F_handle);
                    jac_handle = @(t_curr, y) jacobian_standard(t_curr, y, obj.K, obj.C, target_dofs, gap_wall, k_penalty);
                    
                elseif strcmpi(args.ModelType, 'MC')
                    % Tolleranze
                    Pc = args.ProjectionMatrix;
                    if isempty(Pc)
                        error('For MC method, you must provide the ProjectionMatrix (Pc).');
                    end
                    
                    % OTTIMIZZAZIONE CRITICA per Milman-Chu
                    Pc_contact = Pc(target_dofs, :);

                    der_handle = @(t_curr, y) state_space_projected(t_curr, y, obj.K, obj.C, Pc_contact, gap_wall, k_penalty, F_handle);
                    jac_handle = @(t_curr, y) jacobian_projected(t_curr, y, obj.K, obj.C, Pc_contact, gap_wall, k_penalty);
                else
                    error('ModelType not recognized. Use ''FOM'', ''MC'', ''Rubin'', ''MCB'', or ''MN''.');
                end
                % ============ TOLLERANZE PESATE IN ENERGIA ============
                % AbsTol scalare non e' invariante rispetto alla base di riduzione:
                % applica criteri d'errore diversi a ROM diversi. L'energia meccanica
                % invece e' uno scalare fisico, quindi identica per tutte le basi.
                %   atol_q(i)  = epsE*sqrt(2*Eref/K_ii)
                %   atol_qd(i) = epsE*sqrt(2*Eref/M_ii)
                if isempty(args.Eref)
                    error('Parametro ''Eref'' mancante: calcolarlo una volta dal FOM.');
                end
                adaptive_reltol = args.RelTol;
                epsE = 1e-2 * adaptive_reltol;

                % rigidezza di contatto, caso peggiore (tutti i DOF attivi)
                dK_c = zeros(n_dofs, 1);
                if strcmpi(args.ModelType, 'MC')
                    dK_c = full(sum((Pc_contact.^2) .* k_penalty, 1)).';
                else
                    dK_c(target_dofs) = k_penalty;
                end

                dM = abs(full(diag(obj.M)));
                dM = max(dM, 1e-6 * median(dM(dM > 0)));          % guardia massless
                dK = abs(full(diag(obj.K))) + dK_c;
                dK = max(dK, (2*pi/tmax)^2 .* dM);                % floor sui modi lenti

                adaptive_abstol = epsE * [sqrt(2*args.Eref ./ dK); ...
                    sqrt(2*args.Eref ./ dM)];

                fprintf('  AbsTol energetica [%.2e .. %.2e] | RelTol %.1e\n', ...
                    min(adaptive_abstol), max(adaptive_abstol), adaptive_reltol);
                % ======================================================
                % ==============================================================
                % OPZIONI ODE15S CON JACOBIANO E FLAG SINGOLARITA'
                % ==============================================================
                options = odeset('RelTol', adaptive_reltol, 'AbsTol', adaptive_abstol, 'MaxStep', dt, ...
                    'Mass', M_state, 'MassSingular', mass_singular, 'Jacobian', jac_handle, 'Stats','on');

                % ==============================================================
                % ESECUZIONE INTEGRAZIONE
                % ==============================================================
                fprintf('Integrating %s model with analytical Jacobian (RelTol: %g)...\n', upper(args.ModelType), adaptive_reltol);
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
            else
                error('Linear simulation currently unsupported in this version.');
            end
            
            % ===============================================================
            % LOCAL FUNCTIONS: STATE-SPACE & ANALYTICAL JACOBIANS
            % ===============================================================
            
            % --- 1. FUNZIONI PER FOM, RUBIN, MCB E MN (Nodi Fisici Diretti) ---
            function f = state_space_standard(t_curr, y, K, C, contact_dofs, gap, k_pen, F_ext_handle)
                n = size(K, 1);
                q = y(1:n);
                qd = y(n+1:end);
                
                % Compenetrazione col segno (positiva o negativa)
                penetration = q(contact_dofs) - gap;
                
                % is_pen è VERO solo se il nodo supera la parete nella direzione corretta
                is_pen = (sign(gap) .* penetration) > 0;
                
                F_pen = zeros(n, 1);
                
                if any(is_pen)
                    % La forza generata restituisce il segno corretto per opporsi al moto
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
                    % La derivata della forza elastica K*(q-gap) è sempre positiva (+K)
                    K_pen = sparse(active_dofs, active_dofs, k_pen(is_pen), n, n);
                    K_eff = K + K_pen;
                else
                    K_eff = K;
                end
                
                Z = sparse(n, n);
                I = speye(n);
                J = [Z, I; -K_eff, -C];
            end

            % --- 2. FUNZIONI PER MC (Nodi Proiettati) ---
            function f = state_space_projected(t_curr, y, K, C, Pc_contact, gap, k_pen, F_ext_handle)
                n = size(K, 1);
                q = y(1:n);
                qd = y(n+1:end);
                
                q_phys_contact = Pc_contact * q;
                penetration = q_phys_contact - gap;
                
                is_pen = (sign(gap) .* penetration) > 0;
                
                F_pen_rom = zeros(n, 1);
                if any(is_pen)
                    F_pen_contact = k_pen(is_pen) .* penetration(is_pen);
                    Pc_active = Pc_contact(is_pen, :);
                    F_pen_rom = Pc_active' * F_pen_contact;
                end
                
                f = [qd; F_ext_handle(t_curr) - K*q - C*qd - F_pen_rom];
            end

            function J = jacobian_projected(~, y, K, C, Pc_contact, gap, k_pen)
                n = size(K, 1);
                q = y(1:n);
                
                q_phys_contact = Pc_contact * q;
                penetration = q_phys_contact - gap;
                
                is_pen = (sign(gap) .* penetration) > 0;
                
                if any(is_pen)
                    Pc_active = Pc_contact(is_pen, :);
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






