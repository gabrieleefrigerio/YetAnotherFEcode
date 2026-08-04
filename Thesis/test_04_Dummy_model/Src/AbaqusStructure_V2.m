classdef AbaqusStructure_V2 < handle
    % ABAQUSSTRUCTURE_V2 Class to import an Abaqus mesh and create a YAFEC model.
    % Handles multiple contact interfaces (Top, Bottom, Left, Right).
    % Robust against missing node sets in the .inp file.
    properties
        % --- Input file ---
        filename = ''; 
        
        % --- Material and Thickness Parameters ---
        thickness = 10e-6;     % Out-of-plane thickness [m]
        elementType = 'QUAD8'; 
        E = 165e9;             % Young's modulus [Pa]
        rho = 2330;            % Density [kg/m^3]
        nu = 0.25;             % Poisson's ratio
        
        % --- YAFEC Objects and Output ---
        nodes
        elements
        bc_nodes
        
        % --- Contact Nodes (V2 Additions) ---
        contact_nodes_T = []; % Top
        contact_nodes_B = []; % Bottom
        contact_nodes_L = []; % Left
        contact_nodes_R = []; % Right
        
        MeshObj
        AssemblyObj
        K
        M
        C
        
        % --- Analysis Results ---
        frequencies
        mode_shapes
    end
    
    methods
        function obj = AbaqusStructure_V2(varargin)
            if nargin > 0 && mod(nargin,2) == 0
                for i = 1:2:nargin
                    obj.(varargin{i}) = varargin{i+1};
                end
            end
        end
        
        function build(obj)
            if isempty(obj.filename)
                error('Specify the filename before calling build()');
            end
            obj.import_mesh(obj.filename);
            obj.setup_yafec_model();
        end
        
        function import_mesh(obj, filename)
            meshinfo = abqmesh(filename);
            
            obj.nodes = meshinfo.nodes;
            obj.elements = meshinfo.elem{1};
            
            % Reset dei set (fondamentale se si richiama build() più volte)
            obj.bc_nodes = [];
            obj.contact_nodes_T = [];
            obj.contact_nodes_B = [];
            obj.contact_nodes_L = [];
            obj.contact_nodes_R = [];
            
            fid = fopen(filename, 'r');
            if fid == -1
                error('Unable to open file %s to read Sets.', filename);
            end
            
            current_set_name = '';
            while ~feof(fid)
                line = strtrim(fgetl(fid));
                if isempty(line) || startsWith(line, '**')
                    continue; 
                end
                
                if startsWith(line, '*Nset', 'IgnoreCase', true)
                    tokens = regexp(line, 'nset=([^,]+)', 'tokens', 'ignorecase');
                    if ~isempty(tokens)
                        current_set_name = strtrim(tokens{1}{1});
                    else
                        current_set_name = '';
                    end
                elseif startsWith(line, '*') && ~startsWith(line, '*Nset', 'IgnoreCase', true)
                    current_set_name = '';
                elseif ~isempty(current_set_name)
                    nums = sscanf(line, '%f,');
                    
                    % Smistamento nodi in base al nome del set
                    if strcmpi(current_set_name, 'BottomFixed') || strcmpi(current_set_name, 'TopFixed')
                        obj.bc_nodes = [obj.bc_nodes; nums];
                    elseif strcmpi(current_set_name, 'ContactInterface_T')
                        obj.contact_nodes_T = [obj.contact_nodes_T; nums];
                    elseif strcmpi(current_set_name, 'ContactInterface_B')
                        obj.contact_nodes_B = [obj.contact_nodes_B; nums];
                    elseif strcmpi(current_set_name, 'ContactInterface_L')
                        obj.contact_nodes_L = [obj.contact_nodes_L; nums];
                    elseif strcmpi(current_set_name, 'ContactInterface_R')
                        obj.contact_nodes_R = [obj.contact_nodes_R; nums];
                    end
                end
            end
            fclose(fid);
            
            % Rimozione di eventuali duplicati
            if ~isempty(obj.bc_nodes), obj.bc_nodes = unique(obj.bc_nodes); end
            if ~isempty(obj.contact_nodes_T), obj.contact_nodes_T = unique(obj.contact_nodes_T); end
            if ~isempty(obj.contact_nodes_B), obj.contact_nodes_B = unique(obj.contact_nodes_B); end
            if ~isempty(obj.contact_nodes_L), obj.contact_nodes_L = unique(obj.contact_nodes_L); end
            if ~isempty(obj.contact_nodes_R), obj.contact_nodes_R = unique(obj.contact_nodes_R); end
            
            % Output a schermo selettivo (stampa solo i set effettivamente trovati)
            fprintf('--- Mesh successfully imported from %s ---\n', filename);
            fprintf(' Total nodes: %d\n', size(obj.nodes, 1));
            if ~isempty(obj.bc_nodes), fprintf(' Fixed Nodes: %d\n', length(obj.bc_nodes)); end
            if ~isempty(obj.contact_nodes_T), fprintf(' Contact Nodes (Top): %d\n', length(obj.contact_nodes_T)); end
            if ~isempty(obj.contact_nodes_B), fprintf(' Contact Nodes (Bottom): %d\n', length(obj.contact_nodes_B)); end
            if ~isempty(obj.contact_nodes_L), fprintf(' Contact Nodes (Left): %d\n', length(obj.contact_nodes_L)); end
            if ~isempty(obj.contact_nodes_R), fprintf(' Contact Nodes (Right): %d\n', length(obj.contact_nodes_R)); end
        end
        
        function setup_yafec_model(obj)
            myMaterial = KirchoffMaterial();
            set(myMaterial, 'YOUNGS_MODULUS', obj.E, 'DENSITY', obj.rho, 'POISSONS_RATIO', obj.nu);
            myMaterial.PLANE_STRESS = true;
            
            switch upper(obj.elementType)
                case 'TRI3',  myConstructor = @()Tri3Element(obj.thickness, myMaterial);
                case 'TRI6',  myConstructor = @()Tri6Element(obj.thickness, myMaterial);
                case 'QUAD4', myConstructor = @()Quad4Element(obj.thickness, myMaterial);
                case 'QUAD8', myConstructor = @()Quad8Element(obj.thickness, myMaterial);
                otherwise, error('Element type not supported');
            end
            
            obj.MeshObj = Mesh(obj.nodes);
            obj.MeshObj.create_elements_table(obj.elements, myConstructor);
            
            if ~isempty(obj.bc_nodes)
                obj.MeshObj.set_essential_boundary_condition(obj.bc_nodes, 1:2, 0);
            end
            
            obj.AssemblyObj = Assembly(obj.MeshObj);
            obj.M = obj.AssemblyObj.mass_matrix();
            u0 = zeros(obj.MeshObj.nDOFs, 1);
            [obj.K, ~] = obj.AssemblyObj.tangent_stiffness_and_force(u0);
            
            obj.AssemblyObj.DATA.K = obj.K;
            obj.AssemblyObj.DATA.M = obj.M;
        end
        
        function plot_undeformed(obj, varargin)
            if isempty(obj.nodes) || isempty(obj.elements), error('Mesh not found.'); end
            switch upper(obj.elementType)
                case 'TRI3',  ind_elem_plot = 1:3;
                case 'TRI6',  ind_elem_plot = [1 4 2 5 3 6];
                case 'QUAD4', ind_elem_plot = 1:4;
                case 'QUAD8', ind_elem_plot = [1 5 2 6 3 7 4 8];
            end
            elementPlot = obj.elements(:, ind_elem_plot);
            
            if nargin > 1 && ~isempty(varargin{1}), colorField = varargin{1};
            else, colorField = zeros(size(obj.nodes, 1), 2); end
            
            figure('Name', 'Undeformed Structure', 'Color', 'w', 'Units', 'normalized', 'Position', [0.3 0.25 0.4 0.6]);
            PlotFieldonDeformedMesh(obj.nodes, elementPlot, colorField, 'factor', 1e-12);
            title('Undeformed Structure');
            colormap jet;
            colorbar;
            try
                clim([0, 1e-9]);
            catch
                caxis([0, 1e-9]);
            end
            set(findobj(gca, '-property', 'Marker'), 'Marker', 'none');
            axis equal; grid on;
        end
        
        function plot_mode(obj, mode_idx, scale_factor)
            if isempty(obj.frequencies) || isempty(obj.mode_shapes)
                error('No eigenmodes found. Run compute_eigenmodes() first.');
            end
            if mode_idx > length(obj.frequencies)
                error('Requested mode index exceeds the number of computed modes.');
            end
            
            switch upper(obj.elementType)
                case 'TRI3',  ind_elem_plot = 1:3;
                case 'TRI6',  ind_elem_plot = [1 4 2 5 3 6];
                case 'QUAD4', ind_elem_plot = 1:4;
                case 'QUAD8', ind_elem_plot = [1 5 2 6 3 7 4 8];
            end
            elementPlot = obj.elements(:, ind_elem_plot);
            
            if nargin < 3 || isempty(scale_factor)
                scale_factor = 500e-6 * 0.2;
            end
            
            v1 = reshape(obj.mode_shapes(:, mode_idx), 2, []).';
            figure('Name', ['Mode Shape ' num2str(mode_idx)], 'Color', 'w', ...
                'Units', 'normalized', 'Position', [0.3 0.25 0.4 0.6]);
            PlotMesh(obj.nodes, elementPlot, 0);
            hold on;
            PlotFieldonDeformedMesh(obj.nodes, elementPlot, v1, 'factor', scale_factor);
            colormap jet;
            colorbar;
            set(findobj(gca, '-property', 'Marker'), 'Marker', 'none');
            title(['\Phi_{' num2str(mode_idx) '} - Frequency = ' num2str(obj.frequencies(mode_idx), 4) ' Hz']);
            xlabel('X [m]');
            ylabel('Y [m]');
            axis equal;
            grid on;
        end
        
        function compute_eigenmodes(obj, n_modes)
            if isempty(obj.K) || isempty(obj.M)
                error('Model matrices not found. Run build() before calculating eigenmodes.');
            end
            
            Kc = obj.AssemblyObj.constrain_matrix(obj.K);
            Mc = obj.AssemblyObj.constrain_matrix(obj.M);
            [V0, om] = eigs(Kc, Mc, n_modes, 'SM');
            [obj.frequencies, ind] = sort(sqrt(diag(om)) / (2 * pi));
            V0 = V0(:, ind);
            
            for ii = 1:n_modes
                V0(:, ii) = V0(:, ii) / max(sqrt(sum(V0(:, ii).^2, 2)));
            end
            
            obj.mode_shapes = obj.AssemblyObj.unconstrain_vector(V0);
            fprintf('--- Modal Analysis Complete ---\n');
            for ii = 1:n_modes
                fprintf(' Mode %d: %.3f Hz\n', ii, obj.frequencies(ii));
            end
        end
        
        function dofs_constrained = get_specific_contact_dofs(obj, location, dir)
            % GET_SPECIFIC_CONTACT_DOFS Returns the constrained DOFs for a specific interface.
            % Inputs:
            %   - location: 'T' (Top), 'B' (Bottom), 'L' (Left), or 'R' (Right)
            %   - dir: 1 for X direction, 2 for Y direction
            
            if nargin < 3
                error('Please specify both location (''T'',''B'',''L'',''R'') and direction (1 or 2).');
            end
            
            switch upper(location)
                case 'T', target_nodes = obj.contact_nodes_T;
                case 'B', target_nodes = obj.contact_nodes_B;
                case 'L', target_nodes = obj.contact_nodes_L;
                case 'R', target_nodes = obj.contact_nodes_R;
                otherwise
                    error('Invalid location. Use ''T'', ''B'', ''L'', or ''R''.');
            end
            
            % Ritorna array vuoto se il nodeset non è presente (nessun errore/warning)
            if isempty(target_nodes)
                dofs_constrained = [];
                return;
            end
            
            dofs_global = (target_nodes - 1) * obj.MeshObj.nDOFPerNode + dir;
            dofs_constrained = obj.AssemblyObj.free2constrained_index(dofs_global);
            dofs_constrained = dofs_constrained(dofs_constrained > 0);
        end
        
        function F_c = create_constrained_force_vector(obj, target_node, dof_dir)
            dof_global = (target_node - 1) * obj.MeshObj.nDOFPerNode + dof_dir;
            F_full = zeros(obj.MeshObj.nDOFs, 1);
            F_full(dof_global) = 1;
            F_c = obj.AssemblyObj.constrain_vector(F_full);
        end
        
        function [C, alpha_ray, beta_ray] = compute_rayleigh_damping(obj, Q1, Q2)
            if isempty(obj.frequencies) || length(obj.frequencies) < 2
                obj.compute_eigenmodes(2);
            end
            
            w1 = obj.frequencies(1) * 2 * pi;
            w2 = obj.frequencies(2) * 2 * pi;
            
            zeta_1 = 1 / (2 * Q1);
            zeta_2 = 1 / (2 * Q2);
            
            alpha_ray = (2 * w1 * w2 * (zeta_1 * w2 - zeta_2 * w1)) / (w2^2 - w1^2);
            beta_ray  = (2 * (zeta_2 * w2 - zeta_1 * w1)) / (w2^2 - w1^2);
            
            if isempty(obj.K) || isempty(obj.M)
                error('K and M matrices not found. Make sure to run build() first.');
            end
            
            obj.C = alpha_ray * obj.M + beta_ray * obj.K;
            obj.AssemblyObj.DATA.C = obj.C;
            
            C = obj.C;
            
            fprintf('\n--- Rayleigh Damping Matrix ---\n');
            fprintf('Base frequencies: f1 = %.3f Hz, f2 = %.3f Hz\n', obj.frequencies(1), obj.frequencies(2));
            fprintf('Q1 = %g, Q2 = %g  ->  zeta_1 = %g, zeta_2 = %g\n', Q1, Q2, zeta_1, zeta_2);
            fprintf('alpha = %e\n', alpha_ray);
            fprintf('beta  = %e\n', beta_ray);
        end

        function plot_static_result(obj, U, scale_factor)
            % PLOT_STATIC_RESULT Visualizza la deformata di un'analisi statica
            % Input:
            %   - U: Vettore degli spostamenti (può essere vincolato o completo)
            %   - scale_factor: (Opzionale) Fattore di scala per enfatizzare la deformazione
            
            if nargin < 3 || isempty(scale_factor)
                scale_factor = 1; % Di default scala reale
            end
            
            % 1. Controllo dimensionale ed estensione del vettore ai GdL bloccati
            n_dofs_full = obj.MeshObj.nDOFs;
            
            if length(U) < n_dofs_full
                % Se il vettore è più corto del numero totale di GdL, è quello "ridotto"
                U_full = obj.AssemblyObj.unconstrain_vector(U);
            elseif length(U) == n_dofs_full
                % Se ha la stessa lunghezza, è già il vettore globale completo
                U_full = U;
            else
                error('Dimensione del vettore degli spostamenti non compatibile.');
            end
            
            % 2. Mappatura della tipologia di elemento per il plot
            switch upper(obj.elementType)
                case 'TRI3',  ind_elem_plot = 1:3;
                case 'TRI6',  ind_elem_plot = [1 4 2 5 3 6];
                case 'QUAD4', ind_elem_plot = 1:4;
                case 'QUAD8', ind_elem_plot = [1 5 2 6 3 7 4 8];
            end
            elementPlot = obj.elements(:, ind_elem_plot);
            
            % 3. Formattazione degli spostamenti per la funzione di YAFEC
            U_plot = reshape(U_full, 2, []).';
            
            % 4. Generazione della figura
            figure('Name', 'Static Analysis - Deformed Shape', 'Color', 'w', ...
                'Units', 'normalized', 'Position', [0.3 0.25 0.4 0.6]);
            
            % Plot della mesh indeformata in background
            PlotMesh(obj.nodes, elementPlot, 0);
            hold on;
            
            % Plot della mesh deformata
            PlotFieldonDeformedMesh(obj.nodes, elementPlot, U_plot, 'factor', scale_factor);
            
            colormap jet;
            colorbar;
            set(findobj(gca, '-property', 'Marker'), 'Marker', 'none');
            title(sprintf('Static Deformation (Scale Factor: %g)', scale_factor));
            xlabel('X [m]');
            ylabel('Y [m]');
            axis equal;
            grid on;
        end
    end
end