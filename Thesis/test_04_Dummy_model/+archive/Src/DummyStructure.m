classdef DummyStructure < handle
    % DUMMYSTRUCTURE Class to create and mesh the structure with YAFEC.
    % Models a suspension beam clamped at the ends with a
    % large lateral mass, assuming a 2D plane stress state with thickness.

    properties
        % --- Geometry [m] ---
        W_b = 3.0e-6;  % Beam width (suspension beam)
        L_b = 200.0e-6 * 2 + 500e-6;   % Total beam length
        Y_m = 200.0e-6;   % Y-coordinate where the mass starts
        H_m = 500e-6;   % Mass height
        W_m = 500e-6 - 3.0e-6;   % Mass width
        thickness = 10e-6; % Out-of-plane thickness

        % --- Discretization parameters (number of elements) ---
        nx_b = 2;    % Elem. along the beam width
        nx_m = 50;    % Elem. along the mass width
        ny_1 = 50;   % Elem. along Y for the lower portion of the beam
        ny_2 = 50;   % Elem. along Y for the central portion (and mass)
        ny_3 = 50;   % Elem. along Y for the upper portion of the beam

        elementType = 'QUAD8'; % Element type

        % --- Material (Polysilicon) ---
        E = 165e9;   % Young's Modulus [Pa]
        rho = 2330;  % Density [kg/m^3]
        nu = 0.25;   % Poisson's Ratio

        % --- YAFEC objects and output ---
        nodes
        elements
        bc_nodes
        MeshObj
        AssemblyObj
        K
        M
        C

        % --- Analysis Results ---
        frequencies  % Natural frequencies in Hz
        mode_shapes  % Normalized mode shapes
    end

    methods
        function obj = DummyStructure(varargin)
            % Constructor: allows overwriting properties
            if nargin > 0 && mod(nargin,2) == 0
                for i = 1:2:nargin
                    obj.(varargin{i}) = varargin{i+1};
                end
            end
        end

        function build(obj)
            % Main method that executes the entire procedure
            obj.generate_merged_mesh();
            obj.setup_yafec_model();
        end

        function generate_merged_mesh(obj)
            % Generates and merges the 4 rectangular portions

            % 1. Bottom Beam (0 <= x <= W_b,  0 <= y <= Y_m)
            [n1, e1, ~] = mesh_2Drectangle(obj.W_b, obj.Y_m, obj.nx_b, obj.ny_1, obj.elementType);

            % 2. Middle Beam (0 <= x <= W_b,  Y_m <= y <= Y_m + H_m)
            [n2, e2, ~] = mesh_2Drectangle(obj.W_b, obj.H_m, obj.nx_b, obj.ny_2, obj.elementType);
            n2(:,2) = n2(:,2) + obj.Y_m; % Translation along Y

            % 3. Top Beam (0 <= x <= W_b,  Y_m+H_m <= y <= L_b)
            L_top = obj.L_b - obj.Y_m - obj.H_m;
            [n3, e3, ~] = mesh_2Drectangle(obj.W_b, L_top, obj.nx_b, obj.ny_3, obj.elementType);
            n3(:,2) = n3(:,2) + obj.Y_m + obj.H_m;

            % 4. Big Mass (W_b <= x <= W_b+W_m,  Y_m <= y <= Y_m + H_m)
            [n4, e4, ~] = mesh_2Drectangle(obj.W_m, obj.H_m, obj.nx_m, obj.ny_2, obj.elementType);
            n4(:,1) = n4(:,1) + obj.W_b; % Translation along X
            n4(:,2) = n4(:,2) + obj.Y_m; % Translation along Y

            % --- Data Concatenation ---
            nodes_all = [n1; n2; n3; n4];

            % Update the connectivity index for the elements
            e2 = e2 + size(n1,1);
            e3 = e3 + size(n1,1) + size(n2,1);
            e4 = e4 + size(n1,1) + size(n2,1) + size(n3,1);
            elements_all = [e1; e2; e3; e4];

            % --- Removal of Duplicate Nodes ---
            % Coincident nodes on the interfaces are merged to create
            % a continuous mesh. We use a spatial tolerance.
            tol = 1e-9;
            nodes_round = round(nodes_all / tol) * tol;
            [~, ia, ic] = unique(nodes_round, 'rows', 'stable');

            obj.nodes = nodes_all(ia, :);       % Unique nodes (keeps original precision)
            obj.elements = ic(elements_all);    % Update connectivity

            % --- Identification of Boundary Conditions ---
            % Clamped at the bottom (y=0) and top (y=L_b) ends of the beam
            tol_bc = 1e-6;
            bottom_nodes = find(abs(obj.nodes(:,2) - 0) < tol_bc);
            top_nodes    = find(abs(obj.nodes(:,2) - obj.L_b) < tol_bc);

            obj.bc_nodes = [bottom_nodes; top_nodes];

            fprintf('--- Mesh successfully generated ---\n');
            fprintf(' Total nodes: %d\n', size(obj.nodes, 1));
            fprintf(' Total elements: %d\n', size(obj.elements, 1));
        end

        function setup_yafec_model(obj)
            % Configure material, mesh, and matrices in yafec
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
            obj.MeshObj.set_essential_boundary_condition(obj.bc_nodes, 1:2, 0);

            obj.AssemblyObj = Assembly(obj.MeshObj);
            obj.M = obj.AssemblyObj.mass_matrix();
            u0 = zeros(obj.MeshObj.nDOFs, 1);
            [obj.K, ~] = obj.AssemblyObj.tangent_stiffness_and_force(u0);

            obj.AssemblyObj.DATA.K = obj.K;
            obj.AssemblyObj.DATA.M = obj.M;
            % Rayleigh Damping
            alpha_ray = 30;
            beta_ray = 8.4e-9;
            obj.C = alpha_ray * obj.M + beta_ray * obj.K;
            obj.AssemblyObj.DATA.C = obj.C;
        end

        function compute_eigenmodes(obj, n_modes)
            % COMPUTE_EIGENMODES Calculates the first n_modes natural frequencies and mode shapes
            if isempty(obj.K) || isempty(obj.M)
                error('Model matrices not found. Run build() before calculating eigenmodes.');
            end

            % Constrain matrices using YAFEC built-in reduction
            Kc = obj.AssemblyObj.constrain_matrix(obj.K);
            Mc = obj.AssemblyObj.constrain_matrix(obj.M);

            % Solve the generalized eigenvalue problem for the smallest magnitudes ('SM')
            [V0, om] = eigs(Kc, Mc, n_modes, 'SM');

            % Sort frequencies (Hz) and corresponding mode shapes
            [obj.frequencies, ind] = sort(sqrt(diag(om)) / (2 * pi));
            V0 = V0(:, ind);

            % Normalize mode shapes based on maximum displacement magnitude
            for ii = 1:n_modes
                V0(:, ii) = V0(:, ii) / max(sqrt(sum(V0(:, ii).^2, 2)));
            end

            % Map back to full unconstrained system DOFs
            obj.mode_shapes = obj.AssemblyObj.unconstrain_vector(V0);

            fprintf('--- Modal Analysis Complete ---\n');
            for ii = 1:n_modes
                fprintf(' Mode %d: %.3f Hz\n', ii, obj.frequencies(ii));
            end
        end

        function dofs_constrained = get_right_wall_dofs(obj)
            % Trova i nodi sul bordo destro
            max_x = max(obj.nodes(:,1));
            tol = 1e-9;
            nodes_right = find(abs(obj.nodes(:,1) - max_x) < tol);

            % Calcola i GdL globali lungo X
            dofs_global = (nodes_right - 1) * obj.MeshObj.nDOFPerNode + 1;

            % Mappa sui GdL vincolati usando il metodo nativo di YAFEC
            dofs_constrained = obj.AssemblyObj.free2constrained_index(dofs_global);
            dofs_constrained = dofs_constrained(dofs_constrained > 0);
        end

        function F_c = create_constrained_force_vector(obj, target_node, dof_dir)
            % Crea un vettore forza concentrata unitaria e lo vincola
            dof_global = (target_node - 1) * obj.MeshObj.nDOFPerNode + dof_dir;
            F_full = zeros(obj.MeshObj.nDOFs, 1);
            F_full(dof_global) = 1;

            F_c = obj.AssemblyObj.constrain_vector(F_full);
        end

        % METHODS FOR PLOTS

        function plot_undeformed(obj, varargin)
            % PLOT_UNDEFORMED Visualizes the mesh of the undeformed structure.
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
            title(['Undeformed Structure']);
            colormap jet;

            cb = colorbar;

            try
                clim([0, 1e-9]);
            catch
                caxis([0, 1e-9]);
            end

            set(findobj(gca, '-property', 'Marker'), 'Marker', 'none');
            axis equal; grid on;
        end

        function plot_mode(obj, mode_idx, scale_factor)
            % PLOT_MODE Plots a specific deformed mode shape using color shading
            if isempty(obj.frequencies) || isempty(obj.mode_shapes)
                error('No eigenmodes found. Run compute_eigenmodes() first.');
            end
            if mode_idx > length(obj.frequencies)
                error('Requested mode index exceeds the number of computed modes.');
            end

            % Selection of nodes for the element perimeter plotting
            switch upper(obj.elementType)
                case 'TRI3',  ind_elem_plot = 1:3;
                case 'TRI6',  ind_elem_plot = [1 4 2 5 3 6];
                case 'QUAD4', ind_elem_plot = 1:4;
                case 'QUAD8', ind_elem_plot = [1 5 2 6 3 7 4 8];
            end
            elementPlot = obj.elements(:, ind_elem_plot);

            % Assign a default scale factor if empty (20% of total beam length)
            if nargin < 3 || isempty(scale_factor)
                scale_factor = obj.L_b * 0.2;
            end

            % Reshape displacement vector to match Nx2 node topology
            v1 = reshape(obj.mode_shapes(:, mode_idx), 2, []).';

            % Figure creation
            figure('Name', ['Mode Shape ' num2str(mode_idx)], 'Color', 'w', ...
                'Units', 'normalized', 'Position', [0.3 0.25 0.4 0.6]);

            % Plot original mesh outline as reference
            PlotMesh(obj.nodes, elementPlot, 0);
            hold on;

            % Plot the field on the deformed mesh
            PlotFieldonDeformedMesh(obj.nodes, elementPlot, v1, 'factor', scale_factor);

            colormap jet;
            colorbar;
            set(findobj(gca, '-property', 'Marker'), 'Marker', 'none');
            % Formatting
            title(['\Phi_{' num2str(mode_idx) '} - Frequency = ' num2str(obj.frequencies(mode_idx), 4) ' Hz']);
            xlabel('X [m]');
            ylabel('Y [m]');
            axis equal;
            grid on;
        end
function plot_node_transient(obj, t, y_full, node_idx, dof_dir, varargin)
            % PLOT_NODE_TRANSIENT Plots the time history of a specific node and DOF direction.
            % Accepts optional Name-Value pairs for custom titles and wall gap visualization.
            
            p = inputParser;
            addOptional(p, 'TitleStr', 'Transient Response', @ischar);
            addParameter(p, 'WallGap', [], @isnumeric); % Aggiunto per plottare il muro
            parse(p, varargin{:});
            title_base = p.Results.TitleStr;
            gap = p.Results.WallGap;
            
            % Map node index and direction (1=X, 2=Y) to the global DOF row index
            dof_global = (node_idx - 1) * obj.MeshObj.nDOFPerNode + dof_dir;
            
            if dof_global > size(y_full, 1)
                error('The requested node or DOF direction exceeds the size of the displacement matrix.');
            end
            
            dir_labels = {'X', 'Y'};
            
            figure('Name', title_base, 'Color', 'w', ...
                   'Units', 'normalized', 'Position', [0.2 0.5 0.5 0.35]);
            
            % Plot displacement in micrometers
            plot(t, y_full(dof_global, :) * 1e6, 'LineWidth', 2, 'Color', '#0072BD', 'DisplayName', 'Node Displacement');
            hold on;
            
            % Se è stato passato il parametro WallGap, disegna la linea rossa del muro
            if ~isempty(gap)
                yline(gap * 1e6, 'r--', 'LineWidth', 2, 'DisplayName', 'Wall');
                legend('Location', 'best');
            end
            
            title(sprintf('%s - Node %d, Dir %s', title_base, node_idx, dir_labels{dof_dir}));
            xlabel('Time [s]');
            ylabel('Displacement [\mum]');
            grid on;
        end
        
        function plot_contact_force(obj, t, y_full, contact_nodes, dof_dir, gap_wall, k_penalty, varargin)
            % PLOT_CONTACT_FORCE Computes and plots the total dynamic contact force.
            % It evaluates penetration for all provided contact nodes and sums the penalty forces.
            
            p = inputParser;
            addOptional(p, 'TitleStr', 'Total Contact Force', @ischar);
            parse(p, varargin{:});
            
            % Mappa tutti i nodi di contatto ai rispettivi GdL globali
            dofs_global = (contact_nodes - 1) * obj.MeshObj.nDOFPerNode + dof_dir;
            
            % Estrai la sottomatrice degli spostamenti solo per i nodi in contatto
            disps = y_full(dofs_global, :);
            
            % Calcola le compenetrazioni (azzera i valori negativi dove non c'è contatto)
            penetrations = max(0, disps - gap_wall);
            
            % Calcola la forza per ogni singolo nodo (Legge di Hooke della penalità)
            nodal_forces = k_penalty * penetrations;
            
            % Somma le forze lungo le colonne (somma su tutti i nodi per ogni istante di tempo)
            total_force = sum(nodal_forces, 1);
            
            figure('Name', p.Results.TitleStr, 'Color', 'w', ...
                   'Units', 'normalized', 'Position', [0.2 0.1 0.5 0.35]);
                   
            % Plotta la forza totale in micro-Newton
            plot(t, total_force * 1e6, 'LineWidth', 1.5, 'Color', '#D95319');
            
            title(p.Results.TitleStr);
            xlabel('Time [s]');
            ylabel('Contact Force [\muN]');
            grid on;
        end
    end
end