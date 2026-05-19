classdef MicroCantilever < handle
    % MICROCANTILEVER Class to simulate a MEMS beam in YAFEC
    
    properties
        % Geometry [m]
        L           % Beam length
        W           % In-plane width (Y axis in YAFEC)
        thickness   % Out-of-plane thickness
        
        % Material Properties (Silicon)
        E      % Young's modulus [Pa]
        rho     % Density [kg/m^3]
        nu      % Poisson's ratio
        
        % Mesh Parameters
        nx        % Number of elements along X
        ny         % Number of elements along Y
        elementType % Most accurate element type for bending
        
        % Raw Mesh Data for Plotting
        nodes
        elements
        
        % YAFEC Objects
        MeshObj
        AssemblyObj
        K               % Stiffness Matrix (Full Order)
        M               % Mass Matrix (Full Order)
        Kc              % Constrained stiffness matrix (without blocked DOFs)
        Mc              % Constrained mass matrix
    end
    
    methods
        % --- CONSTRUCTOR ---
        function obj = MicroCantilever(cfg)
            % Initializes parameters from the configuration struct
            if nargin > 0
                % 1. Geometry
                obj.L = cfg.geom.L;
                obj.W = cfg.geom.W;
                obj.thickness = cfg.geom.thickness;
                
                % 2. Material
                obj.E = cfg.material.E;
                obj.rho = cfg.material.rho;
                obj.nu = cfg.material.nu;
                
                % 3. Mesh
                obj.nx = cfg.mesh.nx;
                obj.ny = cfg.mesh.ny;
                obj.elementType = cfg.mesh.elementType;
            end
        end
        
        % --- METHOD 1: BUILD YAFEC MODEL ---
        function buildModel(obj)
            % 1. Material Definition
            mat = KirchoffMaterial();
            set(mat, 'YOUNGS_MODULUS', obj.E, 'DENSITY', obj.rho, 'POISSONS_RATIO', obj.nu);
            mat.PLANE_STRESS = true; 
            
            % 2. Element and Mesh
            elemConstructor = @()Quad8Element(obj.thickness, mat);
            [obj.nodes, obj.elements, nset] = mesh_2Drectangle(obj.L, obj.W, obj.nx, obj.ny, obj.elementType);
            
            obj.MeshObj = Mesh(obj.nodes);
            obj.MeshObj.create_elements_table(obj.elements, elemConstructor);
            
            % 3. Boundary conditions (Clamped on the left edge)
            % Assuming nset{1} contains the nodes on the left edge (X=0)
            obj.MeshObj.set_essential_boundary_condition(nset{1}, 1:2, 0);
            
            % 4. Assembly
            obj.AssemblyObj = Assembly(obj.MeshObj);
            obj.M = obj.AssemblyObj.mass_matrix();
            
            u0 = zeros(obj.MeshObj.nDOFs, 1);
            [obj.K, ~] = obj.AssemblyObj.tangent_stiffness_and_force(u0);
            
            % 5. Constrained matrices for subsequent calculations
            obj.Kc = obj.AssemblyObj.constrain_matrix(obj.K);
            obj.Mc = obj.AssemblyObj.constrain_matrix(obj.M);
            
            fprintf('Model successfully built. Active degrees of freedom: %d\n', size(obj.Kc,1));
        end
        
        % --- METHOD 2: MODAL ANALYSIS (FULL ORDER) ---
        function [freqs, modes] = runFullOrderModal(obj, numModes)
            % Solves the eigenvalue problem
            [V_c, om] = eigs(obj.Kc, obj.Mc, numModes, 'SM');
            
            % Sorts and computes frequencies in Hz
            [freqs, ind] = sort(sqrt(diag(om)) / (2*pi));
            V_c = V_c(:, ind);
            
            % Normalizes and reconstructs the full unconstrained vectors
            modes = zeros(obj.MeshObj.nDOFs, numModes);
            for ii = 1:numModes
                V_c(:,ii) = V_c(:,ii) / max(sqrt(sum(V_c(:,ii).^2, 2)));
                modes(:,ii) = obj.AssemblyObj.unconstrain_vector(V_c(:,ii));
            end
            
            % Prints results
            fprintf('\n--- Natural Frequencies (Full Order) ---\n');
            for i = 1:numModes
                fprintf('Mode %d: %.2f Hz\n', i, freqs(i));
            end
        end
        
        % --- METHOD 3: PLOT MODE SHAPE ---
        function plotMode(obj, modeVector, titleStr)
            % Extracts nodes and elements directly from the object
            nodesPlot = obj.nodes;
            elementPlot = obj.elements(:, 1:4);
            
            figure;
            PlotMesh(nodesPlot, elementPlot, 0);
            
            % Reshapes the vector for plotting
            v_plot = reshape(modeVector, 2, []).';
            PlotFieldonDeformedMesh(nodesPlot, elementPlot, v_plot, 'factor', obj.L*0.2);
            title(titleStr);
        end
    end
end