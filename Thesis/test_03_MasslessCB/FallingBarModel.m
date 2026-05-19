classdef FallingBarModel < handle
    % FALLINGBARMODEL Generates the mesh and the M, K matrices for the benchmark 
    % of the falling 1D bar (Monjaraz Tec et al. Section 7.1).
    
    properties
        Assembly        % YAFEC Assembly object
        Mesh            % YAFEC Mesh object
        M_full          % Full mass matrix
        K_full          % Full stiffness matrix
        nodes           % Node coordinates
        nl_dof          % Global index of the contact DOF (longitudinal)
        nl_node         % Contact node ID
        elements        % Element connectivity
    end
    
    methods
        function obj = FallingBarModel()
            % 1. DATA FROM THE PAPER (Section 7.1)
            E       = 900;           % Young's modulus
            nu      = 0.3;           % Poisson's ratio (irrelevant in pure 1D)
            rho     = 1.0;           % Density
            L       = 10.0;          % Bar length
            % Fix width and height to 1 to have Area = 1 and inertia I = 1/12
            % This way the axial matrix EA/L coincides with E/L
            width   = 1.0;        
            h       = 1.0;      
            
            % 2. MATERIAL AND ELEMENTS SETUP (YAFEC)
            myMaterial = KirchoffMaterial();
            set(myMaterial, 'YOUNGS_MODULUS', E, 'DENSITY', rho, 'POISSONS_RATIO', nu, 'DAMPING_MODULUS', 0);
            
            myElementConstructor = @()BeamElement(width, h, myMaterial);
            
            % 3. 1D MESH CREATION (YAFEC)
            % The paper requires \Delta x = 10^-2. On L = 10, this is 1000 elements.
            dx = 1e-2;
            nx = round(L / dx); 
            x_coords = linspace(0, L, nx+1)';
            nodes_coord = [x_coords, zeros(nx+1, 1)]; 
            
            elements_array = [(1:nx)', (2:nx+1)'];
            
            obj.nodes = nodes_coord;
            obj.elements = elements_array;
            
            obj.Mesh = Mesh(nodes_coord);
            obj.Mesh.create_elements_table(elements_array, myElementConstructor);
            
            % 4. BOUNDARY CONDITIONS (Purely 1D kinematics)
            % The bar is in free flight axially (free-free).
            % We block DOFs 2 (Uy) and 3 (RotZ) for all nodes to
            % force strictly longitudinal wave propagation.
            for i = 1:(nx+1)
                obj.Mesh.set_essential_boundary_condition(i, [2, 3], 0);
            end
            
            % 5. MATRICES ASSEMBLY
            u0 = zeros(obj.Mesh.nDOFs, 1);
            obj.Assembly = Assembly(obj.Mesh);
            obj.M_full = obj.Assembly.mass_matrix();
            [obj.K_full, ~] = obj.Assembly.tangent_stiffness_and_force(u0);
            
            % 6. NONLINEAR DOF IDENTIFICATION (Contact)
            % The impact occurs at the bottom end x = 0 (node 1).
            obj.nl_node = 1;
            
            % The contact degree of freedom is the axial displacement (Ux)
            node_dofs = obj.Mesh.get_DOF_from_location(nodes_coord(obj.nl_node, :));
            obj.nl_dof = node_dofs(1); 
        end
        
        function display_linear_modes(obj)
            % Display natural frequencies (Longitudinal free-free modes)
            Kc = obj.Assembly.constrain_matrix(obj.K_full);
            Mc = obj.Assembly.constrain_matrix(obj.M_full);
            
            % Note: being free-free, we will have a zero-frequency mode (rigid body).
            % We add a very small shift to avoid numerical warnings on singular matrices
            shift = 1e-6;
            [~, om2] = eigs(Kc + shift*Mc, Mc, 7, 'smallestabs');
            
            % Remove the effect of the shift
            eigvals = diag(om2) - shift;
            eigvals(eigvals < 0) = 0; % Correction for numerical zero
            
            freqs = sort(sqrt(eigvals) / (2*pi));
            
            fprintf('\n--- BAR NATURAL FREQUENCIES (LONGITUDINAL) ---\n');
            for i = 1:7
                fprintf('Mode %d: %.2f Hz\n', i, freqs(i));
            end
            fprintf('----------------------------------------------\n\n');
        end
    end
end