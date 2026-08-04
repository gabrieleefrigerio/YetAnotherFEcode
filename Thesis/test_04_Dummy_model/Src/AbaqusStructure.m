classdef AbaqusStructure < handle
    % ABAQUSSTRUCTURE Import an Abaqus mesh (.inp) and build the YAFEC model.
    %
    % Handles ANY number of contact interfaces. Node sets whose name starts
    % with 'ContactInterface' are discovered automatically and exposed through
    % a label:
    %
    %   *Nset, nset=ContactInterface      ->  label 'C'   (single interface)
    %   *Nset, nset=ContactInterface_T    ->  label 'T'
    %   *Nset, nset=ContactInterface_LEFT ->  label 'LEFT'
    %
    % The same code therefore runs both on the four-interface model
    % (DummyStructureAbaqus_V4.inp) and on a single-interface one: only the
    % list of discovered labels changes, not the logic of the main.
    %
    % The compact form '*Nset, ..., generate' is supported, where the data
    % line is a triple (first, last, increment) that must be expanded. Reading
    % those three numbers as three node IDs yields wrong contact sets.
    %
    % Contact API:
    %   obj.contact_labels            labels found in the file, sorted
    %   obj.get_contact_nodes(lbl)    node IDs of the interface
    %   obj.get_contact_dofs(lbl, d)  constrained DOFs, d = 1 (X) or 2 (Y)
    %   obj.describe_interfaces()     summary with the node coordinates
    %   obj.plot_contact_interfaces() visual check of the sets that were read

    properties
        % --- Input file ---
        filename = '';

        % --- Material and thickness ---
        thickness = 10e-6;     % Out-of-plane thickness [m]
        elementType = 'QUAD8';
        E = 165e9;             % Young's modulus [Pa]
        rho = 2330;            % Density [kg/m^3]
        nu = 0.25;             % Poisson's ratio

        % --- Mesh ---
        nodes
        elements
        bc_nodes

        % --- Contact interfaces (arbitrary number) ---
        contact_sets = struct();   % one field per label -> node vector
        contact_labels = {};       % sorted list of the labels found

        % --- YAFEC objects ---
        MeshObj
        AssemblyObj
        K
        M
        C

        % --- Analysis results ---
        frequencies
        mode_shapes
    end

    methods
        function obj = AbaqusStructure(varargin)
            if nargin > 0 && mod(nargin,2) == 0
                for i = 1:2:nargin
                    obj.(varargin{i}) = varargin{i+1};
                end
            end
        end

        function build(obj)
            if isempty(obj.filename)
                error('AbaqusStructure:NoFile', 'Set filename before calling build().');
            end
            obj.import_mesh(obj.filename);
            obj.setup_yafec_model();
        end

        % =================================================================
        function import_mesh(obj, filename)
            % IMPORT_MESH Read nodes, elements, boundary and contact node sets.
            meshinfo = abqmesh(filename);

            obj.nodes = meshinfo.nodes;
            obj.elements = meshinfo.elem{1};

            % Reset, needed when build() is called more than once
            obj.bc_nodes = [];
            obj.contact_sets = struct();
            obj.contact_labels = {};

            fid = fopen(filename, 'r');
            if fid == -1
                error('AbaqusStructure:CannotOpen', ...
                    'Cannot open file %s to read the node sets.', filename);
            end
            cleaner = onCleanup(@() fclose(fid));

            current_label = '';      % contact label of the current block
            current_is_bc = false;
            current_gen   = false;

            while ~feof(fid)
                line = strtrim(fgetl(fid));
                if isempty(line) || startsWith(line, '**')
                    continue;
                end

                if startsWith(line, '*Nset', 'IgnoreCase', true)
                    [current_label, current_is_bc, current_gen] = obj.parse_nset_header(line);

                elseif startsWith(line, '*')
                    % Any other keyword closes the current block
                    current_label = '';
                    current_is_bc = false;
                    current_gen   = false;

                elseif current_is_bc || ~isempty(current_label)
                    nums = sscanf(line, '%f,');
                    if current_gen
                        nums = obj.expand_generate(nums);
                    end
                    if current_is_bc
                        obj.bc_nodes = [obj.bc_nodes; nums];
                    else
                        obj.append_contact_nodes(current_label, nums);
                    end
                end
            end

            % Remove duplicates
            if ~isempty(obj.bc_nodes)
                obj.bc_nodes = unique(obj.bc_nodes);
            end
            for i = 1:numel(obj.contact_labels)
                lbl = obj.contact_labels{i};
                obj.contact_sets.(lbl) = unique(obj.contact_sets.(lbl));
            end
            obj.contact_labels = sort(obj.contact_labels);

            % --- Summary ---
            fprintf('--- Mesh imported from %s ---\n', filename);
            fprintf(' Total nodes: %d\n', size(obj.nodes, 1));
            if ~isempty(obj.bc_nodes)
                fprintf(' Fixed nodes: %d\n', numel(obj.bc_nodes));
            end
            if isempty(obj.contact_labels)
                fprintf(' No contact interface found.\n');
            else
                fprintf(' Contact interfaces: %d (%s)\n', ...
                    numel(obj.contact_labels), strjoin(obj.contact_labels, ', '));
                for i = 1:numel(obj.contact_labels)
                    lbl = obj.contact_labels{i};
                    fprintf('   %-6s -> %d nodes\n', lbl, numel(obj.contact_sets.(lbl)));
                end
            end
        end

        % =================================================================
        function setup_yafec_model(obj)
            % SETUP_YAFEC_MODEL Build the mesh, assembly and system matrices.
            myMaterial = KirchoffMaterial();
            set(myMaterial, 'YOUNGS_MODULUS', obj.E, 'DENSITY', obj.rho, 'POISSONS_RATIO', obj.nu);
            myMaterial.PLANE_STRESS = true;

            switch upper(obj.elementType)
                case 'TRI3',  myConstructor = @()Tri3Element(obj.thickness, myMaterial);
                case 'TRI6',  myConstructor = @()Tri6Element(obj.thickness, myMaterial);
                case 'QUAD4', myConstructor = @()Quad4Element(obj.thickness, myMaterial);
                case 'QUAD8', myConstructor = @()Quad8Element(obj.thickness, myMaterial);
                otherwise, error('AbaqusStructure:BadElement', ...
                        'Unsupported element type: %s', obj.elementType);
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

        % =================================================================
        %  CONTACT INTERFACES
        % =================================================================
        function n = n_interfaces(obj)
            % N_INTERFACES Number of contact interfaces found in the file.
            n = numel(obj.contact_labels);
        end

        function nodes_out = get_contact_nodes(obj, label)
            % GET_CONTACT_NODES Node IDs of interface 'label'.
            label = obj.check_label(label);
            nodes_out = obj.contact_sets.(label);
        end

        function dofs_constrained = get_contact_dofs(obj, label, dir)
            % GET_CONTACT_DOFS Constrained DOFs of interface 'label'.
            %   label : interface label ('T', 'R', 'C', ...)
            %   dir   : 1 for the X direction, 2 for Y
            %
            % Returns [] if the set exists but all of its DOFs are locked by
            % the essential boundary conditions.
            if nargin < 3
                error('AbaqusStructure:MissingDir', ...
                    'Both label and direction are required: get_contact_dofs(label, dir).');
            end
            if ~ismember(dir, [1 2])
                error('AbaqusStructure:BadDir', ...
                    'dir must be 1 (X) or 2 (Y), not %g.', dir);
            end

            target_nodes = obj.get_contact_nodes(label);
            if isempty(target_nodes)
                dofs_constrained = [];
                return;
            end

            dofs_global = (target_nodes - 1) * obj.MeshObj.nDOFPerNode + dir;
            dofs_constrained = obj.AssemblyObj.free2constrained_index(dofs_global);
            dofs_constrained = dofs_constrained(dofs_constrained > 0);
            dofs_constrained = dofs_constrained(:);
        end

        function describe_interfaces(obj)
            % DESCRIBE_INTERFACES Geometric summary of the interfaces found.
            % Meant to catch a misread node set immediately: the nodes of a
            % face must be aligned, i.e. have a nearly constant X or Y.
            if isempty(obj.contact_labels)
                fprintf('No contact interface.\n');
                return;
            end
            fprintf('\n--- Contact interfaces ---\n');
            fprintf('%-6s %6s   %-25s %-25s\n', 'Label', 'Nodes', 'X range [m]', 'Y range [m]');
            for i = 1:numel(obj.contact_labels)
                lbl = obj.contact_labels{i};
                n = obj.contact_sets.(lbl);
                if isempty(n), continue; end
                x = obj.nodes(n, 1);
                y = obj.nodes(n, 2);
                fprintf('%-6s %6d   [%9.3e %9.3e]  [%9.3e %9.3e]\n', ...
                    lbl, numel(n), min(x), max(x), min(y), max(y));
            end
            fprintf('\n');
        end

        % =================================================================
        %  ANALYSIS
        % =================================================================
        function compute_eigenmodes(obj, n_modes)
            % COMPUTE_EIGENMODES Free-interface normal modes of the model.
            if isempty(obj.K) || isempty(obj.M)
                error('AbaqusStructure:NoMatrices', ...
                    'Model matrices are missing. Run build() before compute_eigenmodes().');
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
            fprintf('--- Modal analysis complete (%d modes) ---\n', n_modes);
            n_show = min(n_modes, 10);
            for ii = 1:n_show
                fprintf(' Mode %d: %.3f Hz\n', ii, obj.frequencies(ii));
            end
            if n_modes > n_show
                fprintf(' ... Mode %d: %.3f Hz\n', n_modes, obj.frequencies(n_modes));
            end
        end

        function [C, alpha_ray, beta_ray] = compute_rayleigh_damping(obj, Q1, Q2)
            % COMPUTE_RAYLEIGH_DAMPING Rayleigh damping from two quality factors.
            % Also returns alpha and beta, which the massless ROMs use to build
            % an equivalent diagonal modal damping.
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
                error('AbaqusStructure:NoMatrices', ...
                    'K and M matrices are missing. Run build() first.');
            end

            obj.C = alpha_ray * obj.M + beta_ray * obj.K;
            obj.AssemblyObj.DATA.C = obj.C;

            C = obj.C;

            fprintf('\n--- Rayleigh damping ---\n');
            fprintf('Base frequencies: f1 = %.3f Hz, f2 = %.3f Hz\n', ...
                obj.frequencies(1), obj.frequencies(2));
            fprintf('Q1 = %g, Q2 = %g  ->  zeta_1 = %g, zeta_2 = %g\n', ...
                Q1, Q2, zeta_1, zeta_2);
            fprintf('alpha = %e\n', alpha_ray);
            fprintf('beta  = %e\n', beta_ray);
        end

        function F_c = create_constrained_force_vector(obj, target_node, dof_dir)
            % CREATE_CONSTRAINED_FORCE_VECTOR Unit load on one node and direction.
            dof_global = (target_node - 1) * obj.MeshObj.nDOFPerNode + dof_dir;
            F_full = zeros(obj.MeshObj.nDOFs, 1);
            F_full(dof_global) = 1;
            F_c = obj.AssemblyObj.constrain_vector(F_full);
        end

        % =================================================================
        %  PLOTS
        % =================================================================
        function plot_undeformed(obj, varargin)
            % PLOT_UNDEFORMED Undeformed mesh, optionally coloured by a field.
            if isempty(obj.nodes) || isempty(obj.elements)
                error('AbaqusStructure:NoMesh', 'Mesh not found.');
            end
            elementPlot = obj.elements(:, obj.plot_connectivity_index());

            if nargin > 1 && ~isempty(varargin{1})
                colorField = varargin{1};
            else
                colorField = zeros(size(obj.nodes, 1), 2);
            end

            figure('Name', 'Undeformed Structure', 'Color', 'w', 'Units', 'normalized', ...
                'Position', [0.3 0.25 0.4 0.6]);
            PlotFieldonDeformedMesh(obj.nodes, elementPlot, colorField, 'factor', 1e-12);
            title('Undeformed Structure');
            colormap jet;
            colorbar;
            try
                clim([0, 1e-9]);
            catch
                caxis([0, 1e-9]); %#ok<CAXIS>
            end
            set(findobj(gca, '-property', 'Marker'), 'Marker', 'none');
            axis equal; grid on;
        end

        function plot_contact_interfaces(obj)
            % PLOT_CONTACT_INTERFACES Mesh with the nodes of each interface marked.
            % Immediate visual check that the node sets were read correctly.
            if isempty(obj.contact_labels)
                error('AbaqusStructure:NoInterfaces', 'No contact interface to plot.');
            end
            elementPlot = obj.elements(:, obj.plot_connectivity_index());

            figure('Name', 'Contact Interfaces', 'Color', 'w', 'Units', 'normalized', ...
                'Position', [0.3 0.25 0.4 0.6]);
            PlotMesh(obj.nodes, elementPlot, 0);
            hold on;
            markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p'};
            h = gobjects(1, numel(obj.contact_labels));
            for i = 1:numel(obj.contact_labels)
                lbl = obj.contact_labels{i};
                n = obj.contact_sets.(lbl);
                if isempty(n), continue; end
                h(i) = plot(obj.nodes(n,1), obj.nodes(n,2), ...
                    markers{mod(i-1, numel(markers))+1}, ...
                    'MarkerSize', 8, 'LineWidth', 1.5, 'LineStyle', 'none', ...
                    'DisplayName', sprintf('%s (%d nodes)', lbl, numel(n)));
            end
            legend(h(isgraphics(h)), 'Location', 'bestoutside');
            title('Contact interfaces');
            xlabel('X [m]'); ylabel('Y [m]');
            axis equal; grid on;
        end

        function plot_mode(obj, mode_idx, scale_factor)
            % PLOT_MODE Mode shape mode_idx over the undeformed mesh.
            if isempty(obj.frequencies) || isempty(obj.mode_shapes)
                error('AbaqusStructure:NoModes', ...
                    'No mode available. Run compute_eigenmodes() first.');
            end
            if mode_idx > length(obj.frequencies)
                error('AbaqusStructure:BadModeIndex', ...
                    'Mode index beyond the number of computed modes.');
            end
            elementPlot = obj.elements(:, obj.plot_connectivity_index());

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
            title(['\Phi_{' num2str(mode_idx) '} - Frequency = ' ...
                num2str(obj.frequencies(mode_idx), 4) ' Hz']);
            xlabel('X [m]'); ylabel('Y [m]');
            axis equal; grid on;
        end

        function plot_static_result(obj, U, scale_factor)
            % PLOT_STATIC_RESULT Deformed shape of a static analysis.
            %   U            : displacements (constrained or full vector)
            %   scale_factor : scale factor (default 1, true scale)
            if nargin < 3 || isempty(scale_factor)
                scale_factor = 1;
            end

            n_dofs_full = obj.MeshObj.nDOFs;
            if length(U) < n_dofs_full
                U_full = obj.AssemblyObj.unconstrain_vector(U);
            elseif length(U) == n_dofs_full
                U_full = U;
            else
                error('AbaqusStructure:BadSize', ...
                    'Displacement vector size is not compatible with the model.');
            end

            elementPlot = obj.elements(:, obj.plot_connectivity_index());
            U_plot = reshape(U_full, 2, []).';

            figure('Name', 'Static Analysis - Deformed Shape', 'Color', 'w', ...
                'Units', 'normalized', 'Position', [0.3 0.25 0.4 0.6]);
            PlotMesh(obj.nodes, elementPlot, 0);
            hold on;
            PlotFieldonDeformedMesh(obj.nodes, elementPlot, U_plot, 'factor', scale_factor);
            colormap jet;
            colorbar;
            set(findobj(gca, '-property', 'Marker'), 'Marker', 'none');
            title(sprintf('Static Deformation (Scale Factor: %g)', scale_factor));
            xlabel('X [m]'); ylabel('Y [m]');
            axis equal; grid on;
        end
    end

    % =====================================================================
    methods (Access = private)
        function [label, is_bc, is_gen] = parse_nset_header(~, line)
            % PARSE_NSET_HEADER Extract from '*Nset, nset=NAME, ..., generate'
            % the set name, whether it is a boundary set, and whether it uses
            % the compact 'generate' form.
            label  = '';
            is_bc  = false;
            is_gen = ~isempty(regexpi(line, ',\s*generate\s*(,|$)', 'once'));

            tok = regexpi(line, 'nset\s*=\s*([^,]+)', 'tokens', 'once');
            if isempty(tok)
                return;
            end
            name = strtrim(tok{1});

            if strcmpi(name, 'BottomFixed') || strcmpi(name, 'TopFixed')
                is_bc = true;
                return;
            end

            % Any set starting with 'ContactInterface' is an interface.
            prefix = 'ContactInterface';
            if strncmpi(name, prefix, numel(prefix))
                suffix = name(numel(prefix)+1:end);
                suffix = regexprep(suffix, '^[_\-]', '');   % drop the separator
                if isempty(suffix)
                    label = 'C';        % single set, no suffix
                else
                    label = matlab.lang.makeValidName(suffix);
                end
            end
        end

        function nums = expand_generate(~, raw)
            % EXPAND_GENERATE Expand the (first, last, increment) triples of
            % the compact 'generate' form into explicit node IDs.
            raw = raw(:).';
            if isempty(raw) || mod(numel(raw), 3) ~= 0
                warning('AbaqusStructure:BadGenerate', ...
                    ['''generate'' line with %d values (expected a multiple of 3): ' ...
                     'read without expansion.'], numel(raw));
                nums = raw(:);
                return;
            end
            nums = [];
            for i = 1:3:numel(raw)
                first = raw(i);
                last  = raw(i+1);
                step  = raw(i+2);
                if step == 0, step = 1; end
                nums = [nums, first:step:last]; %#ok<AGROW>
            end
            nums = nums(:);
        end

        function append_contact_nodes(obj, label, nums)
            % APPEND_CONTACT_NODES Add nodes to a contact set, creating it if new.
            if ~isfield(obj.contact_sets, label)
                obj.contact_sets.(label) = [];
                obj.contact_labels{end+1} = label;
            end
            obj.contact_sets.(label) = [obj.contact_sets.(label); nums(:)];
        end

        function label = check_label(obj, label)
            % CHECK_LABEL Validate an interface label against the ones found.
            if ~ischar(label) && ~isstring(label)
                error('AbaqusStructure:BadLabel', ...
                    'The interface label must be a string.');
            end
            label = char(label);
            if ~isfield(obj.contact_sets, label)
                error('AbaqusStructure:NoSuchInterface', ...
                    'Interface ''%s'' is not present in the file. Available: %s', ...
                    label, strjoin(obj.contact_labels, ', '));
            end
        end

        function idx = plot_connectivity_index(obj)
            % PLOT_CONNECTIVITY_INDEX Node ordering used to draw each element type.
            switch upper(obj.elementType)
                case 'TRI3',  idx = 1:3;
                case 'TRI6',  idx = [1 4 2 5 3 6];
                case 'QUAD4', idx = 1:4;
                case 'QUAD8', idx = [1 5 2 6 3 7 4 8];
                otherwise, error('AbaqusStructure:BadElement', ...
                        'Unsupported element type: %s', obj.elementType);
            end
        end
    end
end
