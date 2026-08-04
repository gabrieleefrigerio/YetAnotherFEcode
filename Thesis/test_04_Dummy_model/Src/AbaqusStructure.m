classdef AbaqusStructure < handle
    % ABAQUSSTRUCTURE Importa una mesh Abaqus (.inp) e costruisce il modello YAFEC.
    %
    % Gestisce un numero QUALSIASI di interfacce di contatto. I node set del
    % file .inp il cui nome inizia per 'ContactInterface' vengono scoperti
    % automaticamente e resi accessibili tramite un'etichetta:
    %
    %   *Nset, nset=ContactInterface      ->  etichetta 'C'   (interfaccia unica)
    %   *Nset, nset=ContactInterface_T    ->  etichetta 'T'
    %   *Nset, nset=ContactInterface_LEFT ->  etichetta 'LEFT'
    %
    % Lo stesso codice gira quindi sia sul modello a una faccia di contatto
    % (DummyStructureAbaqus.inp) sia su quello a quattro (..._V4.inp): cambia
    % solo la lista di etichette trovate, non la logica del main.
    %
    % Supporta la forma compatta '*Nset, ..., generate', dove la riga dati e'
    % una terna (primo, ultimo, incremento) e va espansa. Leggerne i tre numeri
    % come tre ID di nodo produce set di contatto sbagliati.
    %
    % API contatto:
    %   obj.contact_labels           etichette trovate nel file, ordinate
    %   obj.get_contact_nodes(lbl)   ID dei nodi dell'interfaccia
    %   obj.get_contact_dofs(lbl, d) GdL vincolati, d = 1 (X) oppure 2 (Y)
    %   obj.describe_interfaces()    riepilogo con le coordinate dei nodi
    %   obj.plot_contact_interfaces() controllo visivo dei set letti

    properties
        % --- File di input ---
        filename = '';

        % --- Materiale e spessore ---
        thickness = 10e-6;     % Spessore fuori piano [m]
        elementType = 'QUAD8';
        E = 165e9;             % Modulo di Young [Pa]
        rho = 2330;            % Densita' [kg/m^3]
        nu = 0.25;             % Coefficiente di Poisson

        % --- Mesh ---
        nodes
        elements
        bc_nodes

        % --- Interfacce di contatto (numero arbitrario) ---
        contact_sets = struct();   % un campo per etichetta -> vettore di nodi
        contact_labels = {};       % elenco ordinato delle etichette trovate

        % --- Oggetti YAFEC ---
        MeshObj
        AssemblyObj
        K
        M
        C

        % --- Risultati dell'analisi ---
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
                error('Specificare filename prima di chiamare build().');
            end
            obj.import_mesh(obj.filename);
            obj.setup_yafec_model();
        end

        % =================================================================
        function import_mesh(obj, filename)
            meshinfo = abqmesh(filename);

            obj.nodes = meshinfo.nodes;
            obj.elements = meshinfo.elem{1};

            % Reset (necessario se build() viene richiamato piu' volte)
            obj.bc_nodes = [];
            obj.contact_sets = struct();
            obj.contact_labels = {};

            fid = fopen(filename, 'r');
            if fid == -1
                error('Impossibile aprire il file %s per la lettura dei Set.', filename);
            end
            cleaner = onCleanup(@() fclose(fid));

            current_label = '';      % etichetta contatto del blocco corrente
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
                    % Qualunque altra keyword chiude il blocco corrente
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

            % Rimozione duplicati
            if ~isempty(obj.bc_nodes)
                obj.bc_nodes = unique(obj.bc_nodes);
            end
            for i = 1:numel(obj.contact_labels)
                lbl = obj.contact_labels{i};
                obj.contact_sets.(lbl) = unique(obj.contact_sets.(lbl));
            end
            obj.contact_labels = sort(obj.contact_labels);

            % --- Riepilogo ---
            fprintf('--- Mesh importata da %s ---\n', filename);
            fprintf(' Nodi totali: %d\n', size(obj.nodes, 1));
            if ~isempty(obj.bc_nodes)
                fprintf(' Nodi vincolati: %d\n', numel(obj.bc_nodes));
            end
            if isempty(obj.contact_labels)
                fprintf(' Nessuna interfaccia di contatto trovata.\n');
            else
                fprintf(' Interfacce di contatto: %d (%s)\n', ...
                    numel(obj.contact_labels), strjoin(obj.contact_labels, ', '));
                for i = 1:numel(obj.contact_labels)
                    lbl = obj.contact_labels{i};
                    fprintf('   %-6s -> %d nodi\n', lbl, numel(obj.contact_sets.(lbl)));
                end
            end
        end

        % =================================================================
        function setup_yafec_model(obj)
            myMaterial = KirchoffMaterial();
            set(myMaterial, 'YOUNGS_MODULUS', obj.E, 'DENSITY', obj.rho, 'POISSONS_RATIO', obj.nu);
            myMaterial.PLANE_STRESS = true;

            switch upper(obj.elementType)
                case 'TRI3',  myConstructor = @()Tri3Element(obj.thickness, myMaterial);
                case 'TRI6',  myConstructor = @()Tri6Element(obj.thickness, myMaterial);
                case 'QUAD4', myConstructor = @()Quad4Element(obj.thickness, myMaterial);
                case 'QUAD8', myConstructor = @()Quad8Element(obj.thickness, myMaterial);
                otherwise, error('Tipo di elemento non supportato: %s', obj.elementType);
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
        %  INTERFACCE DI CONTATTO
        % =================================================================
        function n = n_interfaces(obj)
            n = numel(obj.contact_labels);
        end

        function nodes_out = get_contact_nodes(obj, label)
            % GET_CONTACT_NODES ID dei nodi dell'interfaccia 'label'.
            label = obj.check_label(label);
            nodes_out = obj.contact_sets.(label);
        end

        function dofs_constrained = get_contact_dofs(obj, label, dir)
            % GET_CONTACT_DOFS GdL vincolati dell'interfaccia 'label'.
            %   label : etichetta dell'interfaccia ('T', 'R', 'C', ...)
            %   dir   : 1 per la direzione X, 2 per la Y
            %
            % Restituisce [] se il set esiste ma tutti i suoi GdL sono bloccati
            % dalle condizioni al contorno.
            if nargin < 3
                error('Servono etichetta e direzione: get_contact_dofs(label, dir).');
            end
            if ~ismember(dir, [1 2])
                error('dir deve valere 1 (X) oppure 2 (Y), non %g.', dir);
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
            % DESCRIBE_INTERFACES Riepilogo geometrico delle interfacce trovate.
            % Serve ad accorgersi subito di un node set letto male: i nodi di una
            % faccia devono essere allineati, cioe' avere X oppure Y quasi costante.
            if isempty(obj.contact_labels)
                fprintf('Nessuna interfaccia di contatto.\n');
                return;
            end
            fprintf('\n--- Interfacce di contatto ---\n');
            fprintf('%-6s %6s   %-25s %-25s\n', 'Label', 'Nodi', 'X range [m]', 'Y range [m]');
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
        %  ANALISI
        % =================================================================
        function compute_eigenmodes(obj, n_modes)
            if isempty(obj.K) || isempty(obj.M)
                error('Matrici del modello assenti. Eseguire build() prima di compute_eigenmodes().');
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
            fprintf('--- Analisi modale completata (%d modi) ---\n', n_modes);
            n_show = min(n_modes, 10);
            for ii = 1:n_show
                fprintf(' Modo %d: %.3f Hz\n', ii, obj.frequencies(ii));
            end
            if n_modes > n_show
                fprintf(' ... Modo %d: %.3f Hz\n', n_modes, obj.frequencies(n_modes));
            end
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
                error('Matrici K e M assenti. Eseguire build() prima.');
            end

            obj.C = alpha_ray * obj.M + beta_ray * obj.K;
            obj.AssemblyObj.DATA.C = obj.C;

            C = obj.C;

            fprintf('\n--- Smorzamento di Rayleigh ---\n');
            fprintf('Frequenze di base: f1 = %.3f Hz, f2 = %.3f Hz\n', obj.frequencies(1), obj.frequencies(2));
            fprintf('Q1 = %g, Q2 = %g  ->  zeta_1 = %g, zeta_2 = %g\n', Q1, Q2, zeta_1, zeta_2);
            fprintf('alpha = %e\n', alpha_ray);
            fprintf('beta  = %e\n', beta_ray);
        end

        function F_c = create_constrained_force_vector(obj, target_node, dof_dir)
            dof_global = (target_node - 1) * obj.MeshObj.nDOFPerNode + dof_dir;
            F_full = zeros(obj.MeshObj.nDOFs, 1);
            F_full(dof_global) = 1;
            F_c = obj.AssemblyObj.constrain_vector(F_full);
        end

        % =================================================================
        %  PLOT
        % =================================================================
        function plot_undeformed(obj, varargin)
            if isempty(obj.nodes) || isempty(obj.elements), error('Mesh non trovata.'); end
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
            % PLOT_CONTACT_INTERFACES Mesh con i nodi di ogni interfaccia evidenziati.
            % Controllo visivo immediato che i node set siano stati letti bene.
            if isempty(obj.contact_labels)
                error('Nessuna interfaccia di contatto da disegnare.');
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
                    'DisplayName', sprintf('%s (%d nodi)', lbl, numel(n)));
            end
            legend(h(isgraphics(h)), 'Location', 'bestoutside');
            title('Interfacce di contatto');
            xlabel('X [m]'); ylabel('Y [m]');
            axis equal; grid on;
        end

        function plot_mode(obj, mode_idx, scale_factor)
            if isempty(obj.frequencies) || isempty(obj.mode_shapes)
                error('Nessun modo disponibile. Eseguire compute_eigenmodes() prima.');
            end
            if mode_idx > length(obj.frequencies)
                error('Indice di modo oltre il numero di modi calcolati.');
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
            title(['\Phi_{' num2str(mode_idx) '} - Frequency = ' num2str(obj.frequencies(mode_idx), 4) ' Hz']);
            xlabel('X [m]'); ylabel('Y [m]');
            axis equal; grid on;
        end

        function plot_static_result(obj, U, scale_factor)
            % PLOT_STATIC_RESULT Deformata di un'analisi statica.
            %   U            : spostamenti (vettore vincolato o completo)
            %   scale_factor : fattore di scala (default 1, scala reale)
            if nargin < 3 || isempty(scale_factor)
                scale_factor = 1;
            end

            n_dofs_full = obj.MeshObj.nDOFs;
            if length(U) < n_dofs_full
                U_full = obj.AssemblyObj.unconstrain_vector(U);
            elseif length(U) == n_dofs_full
                U_full = U;
            else
                error('Dimensione del vettore degli spostamenti non compatibile.');
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
            % Estrae dalla riga '*Nset, nset=NOME, ..., generate' il nome del set,
            % se e' un set di vincolo e se usa la forma compatta 'generate'.
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

            % Qualsiasi set che inizia per 'ContactInterface' e' un'interfaccia.
            prefix = 'ContactInterface';
            if strncmpi(name, prefix, numel(prefix))
                suffix = name(numel(prefix)+1:end);
                suffix = regexprep(suffix, '^[_\-]', '');   % via il separatore
                if isempty(suffix)
                    label = 'C';        % set unico, senza suffisso
                else
                    label = matlab.lang.makeValidName(suffix);
                end
            end
        end

        function nums = expand_generate(~, raw)
            % Espande le terne (primo, ultimo, incremento) della forma 'generate'.
            raw = raw(:).';
            if isempty(raw) || mod(numel(raw), 3) ~= 0
                warning('AbaqusStructure:BadGenerate', ...
                    ['Riga ''generate'' con %d valori (attesi multipli di 3): ' ...
                     'letta senza espansione.'], numel(raw));
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
            if ~isfield(obj.contact_sets, label)
                obj.contact_sets.(label) = [];
                obj.contact_labels{end+1} = label;
            end
            obj.contact_sets.(label) = [obj.contact_sets.(label); nums(:)];
        end

        function label = check_label(obj, label)
            if ~ischar(label) && ~isstring(label)
                error('L''etichetta dell''interfaccia deve essere una stringa.');
            end
            label = char(label);
            if ~isfield(obj.contact_sets, label)
                error('AbaqusStructure:NoSuchInterface', ...
                    'Interfaccia ''%s'' non presente nel file. Disponibili: %s', ...
                    label, strjoin(obj.contact_labels, ', '));
            end
        end

        function idx = plot_connectivity_index(obj)
            switch upper(obj.elementType)
                case 'TRI3',  idx = 1:3;
                case 'TRI6',  idx = [1 4 2 5 3 6];
                case 'QUAD4', idx = 1:4;
                case 'QUAD8', idx = [1 5 2 6 3 7 4 8];
                otherwise, error('Tipo di elemento non supportato: %s', obj.elementType);
            end
        end
    end
end
