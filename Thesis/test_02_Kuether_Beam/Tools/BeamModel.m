classdef BeamModel < handle
    % BEAMMODEL Generates the mesh and the M, K matrices of the linear beam 
    % using 1D BeamElements from YAFEC.
    
    properties
        Assembly        % YAFEC Assembly object
        Mesh            % YAFEC Mesh object
        M_full          % Full mass matrix
        K_full          % Full stiffness matrix
        nodes           % Node coordinates
        nl_dof          % Global index of the transverse degree of freedom at the center
        nl_node         % ID of the center node
        elements        % Element connectivity
    end
    
    methods
        function obj = BeamModel()
            % 1. DATA FROM THE PAPER 
            E       = 204.8e9;       % Young's modulus [Pa]
            nu      = 0.28;          % Poisson's ratio 
            rho     = 7866;          % Density [kg/m^3]
            
            L       = 0.2286;        % Length [m]
            h       = 0.000787;      % Thickness (height in the 2D plane) [m]
            width   = 0.0127;        % Width (out-of-plane) [m]
            
            % 2. MATERIAL AND ELEMENT SETUP (YAFEC)
            myMaterial = KirchoffMaterial();
            % We also set DAMPING_MODULUS to 0 to avoid initialization 
            % errors inside the YAFEC BeamElement
            set(myMaterial, 'YOUNGS_MODULUS', E, 'DENSITY', rho, 'POISSONS_RATIO', nu, 'DAMPING_MODULUS', 0);
            
            % The BeamElement takes: (width_b, height_h, Material) as input
            myElementConstructor = @()BeamElement(width, h, myMaterial);
            
            % 3. 1D MESH CREATION (YAFEC)
            % Create a line of nodes along the X axis
            nx = 40; % 40 elements along the beam
            x_coords = linspace(0, L, nx+1)';
            nodes_coord = [x_coords, zeros(nx+1, 1)]; % All nodes lie at Y = 0
            
            % Element topology (connection between node i and node i+1)
            elements_array = [(1:nx)', (2:nx+1)'];
            
            obj.nodes = nodes_coord;
            obj.elements = elements_array;
            
            obj.Mesh = Mesh(nodes_coord);
            obj.Mesh.create_elements_table(elements_array, myElementConstructor);
            
            % 4. BOUNDARY CONDITIONS (Simply Supported for Beam elements)
            % Note: There are 3 DOFs per node (Ux=1, Uy=2, RotZ=3).
            node_left  = 1;
            node_right = nx + 1;
            
            % Left support (pinned): Blocks Ux and Uy, leaves rotation free (DOF 3)
            obj.Mesh.set_essential_boundary_condition(node_left, [1, 2], 0);  
            % Right support (pinned): Blocks Ux and Uy, leaves rotation free (DOF 3)
            obj.Mesh.set_essential_boundary_condition(node_right, [1, 2], 0);  
            
            % 5. MATRIX ASSEMBLY
            u0 = zeros(obj.Mesh.nDOFs, 1);
            obj.Assembly = Assembly(obj.Mesh);
            obj.M_full = obj.Assembly.mass_matrix();
            [obj.K_full, ~] = obj.Assembly.tangent_stiffness_and_force(u0);
            
            % 6. IDENTIFICATION OF THE NONLINEAR DOF
            % Find the node exactly in the middle (L/2)
            obj.nl_node = find(abs(nodes_coord(:,1) - L/2) < 1e-6);
            
            if isempty(obj.nl_node)
                error('Impact node not found at the center of the beam!');
            end
            
            % Get the degrees of freedom associated with the node and take the second one (Y direction)
            node_dofs = obj.Mesh.get_DOF_from_location(nodes_coord(obj.nl_node, :));
            obj.nl_dof = node_dofs(2);
        end
        
        function display_linear_modes(obj)
            % Calculates and prints the natural frequencies to the command window
            Kc = obj.Assembly.constrain_matrix(obj.K_full);
            Mc = obj.Assembly.constrain_matrix(obj.M_full);
            
            % Use 'smallestabs' for compatibility with newer versions of eigs
            [~, om2] = eigs(Kc, Mc, 7, 'smallestabs');
            freqs = sort(sqrt(diag(om2)) / (2*pi));
            
            fprintf('\n--- BEAM NATURAL FREQUENCIES (1D BEAM ELEMENTS) ---\n');
            for i = 1:7
                fprintf('Mode %d: %.2f Hz\n', i, freqs(i));
            end
            fprintf('---------------------------------------------------\n\n');
        end
        
        function plot_mode(obj, mode_idx)
            % PLOT_MODE Calculates and displays the bending mode shape
            
            num_calc = max(10, mode_idx);
            Kc = obj.Assembly.constrain_matrix(obj.K_full);
            Mc = obj.Assembly.constrain_matrix(obj.M_full);
            
            [Phi_c, om2] = eigs(Kc, Mc, num_calc, 'smallestabs');
            [freqs, sort_idx] = sort(sqrt(diag(om2)) / (2*pi));
            Phi_c = Phi_c(:, sort_idx);
            
            f0 = freqs(mode_idx);
            phi_c_mod = Phi_c(:, mode_idx);
            phi_full = obj.Assembly.unconstrain_vector(phi_c_mod);
            
            % Normalize with respect to the maximum transverse displacement
            phi_full = phi_full / max(abs(phi_full));
            
            % Extract only the Y coordinates (the DOFs in the Y direction are 2, 5, 8, 11...)
            n_nodes = size(obj.nodes, 1);
            y_disp = phi_full(2:3:end); % Takes the 2nd DOF of each block of 3
            
            % Visual amplification factor
            beam_length = max(obj.nodes(:,1)) - min(obj.nodes(:,1));
            S = 1.2; 
            y_plot = y_disp * S;
            
            % Plot Creation
            figure('Name', sprintf('Bending Mode %d (Beam Element)', mode_idx), 'Color', 'w');
            
            % Undeformed beam (dashed black line)
            plot(obj.nodes(:,1), zeros(n_nodes, 1), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Undeformed');
            hold on;
            
            % Deformed beam (blue line with red scatter points on nodes)
            plot(obj.nodes(:,1), y_plot, 'b-', 'LineWidth', 2, 'DisplayName', sprintf('Mode %d', mode_idx));
            scatter(obj.nodes(:,1), y_plot, 30, 'r', 'filled', 'HandleVisibility', 'off');
            
            title(sprintf('\\Phi_{%d} - Frequency = %.2f Hz', mode_idx, f0), 'FontSize', 14);
            xlabel('X [m]');
            ylabel('Mode Shape (Visual Scale)');
            % --- SET MARGINS FOR BETTER VISUALIZATION ---
            % X-axis margin: 5% of beam length before and after the supports
            x_margin = 0.05 * beam_length;
            xlim([min(obj.nodes(:,1)) - x_margin, max(obj.nodes(:,1)) + x_margin]);
            
            % Y-axis margin: 20% extra space above and below the max peak
            y_margin = max(abs(y_plot)) + 5*max(abs(y_plot)); 
            ylim([-(y_margin), (y_margin)]); 
            
            legend('Location', 'best');
            grid on; box on;
        end
    end
end