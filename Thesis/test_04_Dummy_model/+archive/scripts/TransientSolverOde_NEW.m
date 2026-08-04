classdef TransientSolverOde_NEW < handle
    properties
        M % Mass matrix
        K % Stiffness matrix
        C % Damping matrix
    end
    
    methods
        function obj = TransientSolverOde_NEW(M, K, varargin)
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
            addParameter(p, 'RelTol', []); % <--- NUOVO: Tolleranza Relativa
            addParameter(p, 'AbsTol', []); % <--- NUOVO: Tolleranza Assoluta
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
                
                % MASS MATRIX per LHS del sistema ode15s
                M_state = blkdiag(speye(n_dofs), obj.M);
                
                % ==============================================================
                % GESTIONE TOLLERANZE DA INPUT
                % ==============================================================
                adaptive_reltol = 1e-4; % Default base
                adaptive_abstol = 1e-5; % Default base
                if strcmpi(args.ModelType, 'Rubin') || strcmpi(args.ModelType, 'MCB') || strcmpi(args.ModelType, 'MN')
                    adaptive_abstol = 1e-8; 
                elseif strcmpi(args.ModelType, 'MC')
                    adaptive_abstol = 1e-10;
                end
                
                % Sovrascrittura se passate dall'utente (per Analisi Sensitività)
                if ~isempty(args.RelTol), adaptive_reltol = args.RelTol; end
                if ~isempty(args.AbsTol), adaptive_abstol = args.AbsTol; end
                
                % ==============================================================
                % SETUP SPECIFICO DEL MODELLO E JACOBIANO ANALITICO
                % ==============================================================
                mass_singular = 'no'; % Default: ode normale
                
                if strcmpi(args.ModelType, 'FOM') || strcmpi(args.ModelType, 'Rubin') || strcmpi(args.ModelType, 'MCB') || strcmpi(args.ModelType, 'MN')
                    
                    if strcmpi(args.ModelType, 'MCB') || strcmpi(args.ModelType, 'MN')
                        mass_singular = 'yes';
                    end
                    
                    der_handle = @(t_curr, y) state_space_standard(t_curr, y, obj.K, obj.C, target_dofs, gap_wall, k_penalty, F_handle);
                    jac_handle = @(t_curr, y) jacobian_standard(t_curr, y, obj.K, obj.C, target_dofs, gap_wall, k_penalty);
                    
                elseif strcmpi(args.ModelType, 'MC')
                    Pc = args.ProjectionMatrix;
                    if isempty(Pc), error('Missing ProjectionMatrix for MC.'); end
                    Pc_contact = Pc(target_dofs, :);
                    
                    der_handle = @(t_curr, y) state_space_projected(t_curr, y, obj.K, obj.C, Pc_contact, gap_wall, k_penalty, F_handle);
                    jac_handle = @(t_curr, y) jacobian_projected(t_curr, y, obj.K, obj.C, Pc_contact, gap_wall, k_penalty);
                else
                    error('ModelType not recognized.');
                end
                
                % ==============================================================
                % OPZIONI ODE15S
                % ==============================================================
                options = odeset('RelTol', adaptive_reltol, 'AbsTol', adaptive_abstol, 'MaxStep', dt, ...
                    'Mass', M_state, 'MassSingular', mass_singular, 'Jacobian', jac_handle);
                    
                % ==============================================================
                % ESECUZIONE INTEGRAZIONE
                % ==============================================================
                fprintf('Integrating %s [RelTol: %g, AbsTol: %g]...\n', upper(args.ModelType), adaptive_reltol, adaptive_abstol);
                tic;
                [t_out, y_out] = ode15s(der_handle, [0 tmax], y0, options);
                integration_time = toc;
                
                t = t_out';
                q_history = y_out(:, 1:n_dofs)';
                % rimosso il log del tempo qui per evitare spam nella console durante il loop
            else
                error('Linear simulation unsupported.');
            end
            
            % --- FUNZIONI LOCALI ---
            function f = state_space_standard(t_curr, y, K, C, contact_dofs, gap, k_pen, F_ext_handle)
                n = size(K, 1); q = y(1:n); qd = y(n+1:end);
                penetration = q(contact_dofs) - gap;
                is_pen = penetration > 0;
                F_pen = zeros(n, 1);
                if any(is_pen), F_pen(contact_dofs(is_pen)) = k_pen * penetration(is_pen); end
                f = [qd; F_ext_handle(t_curr) - K*q - C*qd - F_pen];
            end
            
            function J = jacobian_standard(~, y, K, C, contact_dofs, gap, k_pen)
                n = size(K, 1); q = y(1:n);
                penetration = q(contact_dofs) - gap;
                is_pen = penetration > 0;
                if any(is_pen)
                    active_dofs = contact_dofs(is_pen);
                    K_pen = sparse(active_dofs, active_dofs, k_pen, n, n);
                    K_eff = K + K_pen;
                else
                    K_eff = K;
                end
                J = [sparse(n, n), speye(n); -K_eff, -C];
            end

            function f = state_space_projected(t_curr, y, K, C, Pc_contact, gap, k_pen, F_ext_handle)
                n = size(K, 1); q = y(1:n); qd = y(n+1:end);
                q_phys_contact = Pc_contact * q;
                penetration = q_phys_contact - gap;
                is_pen = penetration > 0;
                F_pen_rom = zeros(n, 1);
                if any(is_pen)
                    Pc_active = Pc_contact(is_pen, :);
                    F_pen_rom = Pc_active' * (k_pen * penetration(is_pen));
                end
                f = [qd; F_ext_handle(t_curr) - K*q - C*qd - F_pen_rom];
            end

            function J = jacobian_projected(~, y, K, C, Pc_contact, gap, k_pen)
                n = size(K, 1); q = y(1:n);
                q_phys_contact = Pc_contact * q;
                penetration = q_phys_contact - gap;
                is_pen = penetration > 0;
                if any(is_pen)
                    Pc_active = Pc_contact(is_pen, :);
                    K_eff = K + Pc_active' * (k_pen * Pc_active);
                else
                    K_eff = K;
                end
                J = [sparse(n, n), speye(n); -K_eff, -C];
            end
        end
    end
end