classdef TransientSolver < handle
    properties
        M % Mass matrix
        K % Stiffness matrix
        C % Damping matrix
    end
    
    methods
        function obj = TransientSolver(M, K, varargin)
            obj.M = M;
            obj.K = K;
            if nargin > 2 && ~isempty(varargin{1})
                obj.C = varargin{1};
            else
                obj.C = sparse(size(K,1), size(K,2));
            end
        end
        
        function set_rayleigh_damping(obj, alpha, beta)
            if alpha == 0 && beta == 0
                obj.C = sparse(size(obj.K,1), size(obj.K,2));
            else
                obj.C = alpha * obj.M + beta * obj.K;
            end
        end
        
        function [t, q_history] = solve(obj, tmax, dt, q0, qd0, F_handle, varargin)
            p = inputParser;
            addParameter(p, 'ContactTargetDOF', []);
            addParameter(p, 'ContactGap', []);
            addParameter(p, 'ContactPenalty', []);
            addParameter(p, 'ProjectionMatrix', []); 
            addParameter(p, 'qdd0', []); 
            addParameter(p, 'ModelType', 'FOM'); % Scelta del modello: 'FOM', 'MC', 'Rubin'
            parse(p, varargin{:});
            args = p.Results;
            
            is_nonlinear = ~isempty(args.ContactTargetDOF) && ~isempty(args.ContactGap) && ~isempty(args.ContactPenalty);
            
            if isscalar(obj.C) && obj.C == 0
                obj.C = sparse(size(obj.K,1), size(obj.K,2));
            end
            
            % --- ASSEGNAZIONE INTELLIGENTE DI qdd0 ---
            if isempty(args.qdd0)
                qdd0_start = obj.M \ (F_handle(0) - obj.C * qd0 - obj.K * q0);
            else
                qdd0_start = args.qdd0;
            end
            
            if is_nonlinear
                target_dofs = args.ContactTargetDOF;
                gap_wall = args.ContactGap;
                k_penalty = args.ContactPenalty;
                
                % ==============================================================
                % BIFORCAZIONE DEI RESIDUI (Function Wrappers)
                % ==============================================================
                if strcmpi(args.ModelType, 'FOM')
                    fprintf('Integrazione modello FOM...\n');
                    res_handle = @(q, qd, qdd, t_curr) residual_fom(q, qd, qdd, t_curr, obj.M, obj.K, obj.C, target_dofs, gap_wall, k_penalty, F_handle);
                    
                elseif strcmpi(args.ModelType, 'MC')
                    fprintf('Integrazione modello ROM (Milman-Chu)...\n');
                    Pc = args.ProjectionMatrix;
                    if isempty(Pc)
                        error('Per il metodo MC devi fornire la ProjectionMatrix (Pc).');
                    end
                    res_handle = @(q, qd, qdd, t_curr) residual_mc(q, qd, qdd, t_curr, obj.M, obj.K, obj.C, Pc, target_dofs, gap_wall, k_penalty, F_handle);
                    
                elseif strcmpi(args.ModelType, 'Rubin')
                    fprintf('Integrazione modello ROM (Rubin Classico)...\n');
                    % In Rubin, target_dofs devono essere gli indici DEI NODI DI CONTATTO 
                    % ALL'INTERNO DEL VETTORE q RIDOTTO (es. i primi N gradi di libertà).
                    res_handle = @(q, qd, qdd, t_curr) residual_rubin(q, qd, qdd, t_curr, obj.M, obj.K, obj.C, target_dofs, gap_wall, k_penalty, F_handle);
                    
                else
                    error('ModelType non riconosciuto. Usa ''FOM'', ''MC'' o ''Rubin''.');
                end
                
                % --- IMPOSTAZIONI NEWMARK E INTEGRAZIONE ---
                tic;
                TI = ImplicitNewmark('timestep', dt, 'alpha', 0.1, 'linear', false, 'ATS', true, 'hmin', 1e-10);
                TI.tol = 1e-6;
                TI.MaxNRit = 15;
                
                TI.Integrate(q0, qd0, qdd0_start, tmax, res_handle);
                integration_time = toc;
                
                t = TI.Solution.time;
                q_history = TI.Solution.q;
                fprintf('Tempo di integrazione: %.2f s\n', integration_time);
            else
                error('Linear simulation currently unsupported in this ROM-ready version.');
            end
        end 
    end
end