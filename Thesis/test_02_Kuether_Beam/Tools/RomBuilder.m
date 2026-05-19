classdef RomBuilder < handle
    % ROMBUILDER Builds the reduction basis (ROM) with Milman-Chu modes

    properties
        Assembly    % Reference to the YAFEC Assembly
        P           % Transformation matrix P = [Phi, Psi]
        numModes    % Number of linear modes to include
        nl_dof      % Index of the degree of freedom where the non-linearity acts
        spring_k    % Stiffness of the contact spring
        include_MC  % Boolean flag to include MC vector in P
    end

    methods
        function obj = RomBuilder(yafec_assembly, num_linear_modes, nonlinear_dof, spring_stiffness, include_MC)
            obj.Assembly = yafec_assembly;
            obj.numModes = num_linear_modes;
            obj.nl_dof = nonlinear_dof;

            % If the spring stiffness is not provided, default to open gap (0)
            if nargin < 4
                obj.spring_k = 0;
            else
                obj.spring_k = spring_stiffness;
            end

            % Se il flag MC non è fornito, di default includiamo il vettore
            if nargin < 5
                obj.include_MC = true;
            else
                obj.include_MC = include_MC;
            end
        end

        function build(obj)
            % 1. Retrieve full matrices from YAFEC
            M_full = obj.Assembly.mass_matrix();
            K_full = obj.Assembly.stiffness_matrix();

            % 2. APPLICATION OF CONSTRAINTS
            Mc = obj.Assembly.constrain_matrix(M_full);
            Kc = obj.Assembly.constrain_matrix(K_full);

            % 3. Calculation of CONSTRAINED linear modes (Phi_c)
            num_to_extract = obj.numModes;
            [Phi_all, D] = eigs(Kc, Mc, num_to_extract, 'smallestabs');

            % Sort by increasing frequency
            [~, sort_idx] = sort(diag(D));
            Phi_all = Phi_all(:, sort_idx);

            % --- ODD MODES FILTER ---
            % We select indices 1, 3, 5... up to numModes
            odd_indices = 1:2:(obj.numModes);
            Phi_c = Phi_all(:, odd_indices);

            % Normalization with respect to the constrained mass
            for i = 1:size(Phi_c,2)
                Phi_c(:,i) = Phi_c(:,i) / sqrt(Phi_c(:,i)' * Mc * Phi_c(:,i));
            end

            if obj.include_MC
                % 4. Calculation of the discontinuous Milman-Chu vector (Psi_c)
                n_dofs_full = obj.Assembly.Mesh.nDOFs;
                F_unit = zeros(n_dofs_full, 1);
                F_unit(obj.nl_dof) = 1;

                F_unit_c = obj.Assembly.constrain_vector(F_unit);
                % Create a full stiffness matrix just for the contact spring
                K_spring_full = sparse(n_dofs_full, n_dofs_full);
                K_spring_full(obj.nl_dof, obj.nl_dof) = obj.spring_k;

                % Constrain the spring matrix to match the active DOFs in YAFEC
                K_spring_c = obj.Assembly.constrain_matrix(K_spring_full);

                % Add the spring stiffness to the beam stiffness
                Kc_closed = Kc + K_spring_c;
                % Calculate the static response with the closed gap
                Psi_c = Kc_closed \ F_unit_c;

                % Base arricchita
                V_initial_c = [Phi_c, Psi_c];
            else
                % Base standard (solo modi lineari)
                V_initial_c = Phi_c;
            end

            % 5. Orthonormalization (Gram-Schmidt)
            P_c = obj.gram_schmidt(V_initial_c, Mc);

            % 6. EXPANSION TO THE FULL SPACE
            obj.P = obj.Assembly.unconstrain_vector(P_c);
        end

        function P_ortho = gram_schmidt(obj, V, M)
            n_vecs = size(V, 2);
            P_ortho = zeros(size(V));
            for i = 1:n_vecs
                v_curr = V(:,i);
                for j = 1:i-1
                    proj_coeff = (P_ortho(:,j)' * M * v_curr) / (P_ortho(:,j)' * M * P_ortho(:,j));
                    v_curr = v_curr - proj_coeff * P_ortho(:,j);
                end
                P_ortho(:,i) = v_curr / sqrt(v_curr' * M * v_curr);
            end
        end

        function plot_milman_chu_comparison(obj)
            original_k = obj.spring_k;
            n_dofs_full = obj.Assembly.Mesh.nDOFs;

            % Retrieve mass and stiffness matrices from YAFEC
            M_full = obj.Assembly.mass_matrix();
            K_full = obj.Assembly.stiffness_matrix();

            Mc = obj.Assembly.constrain_matrix(M_full);
            Kc = obj.Assembly.constrain_matrix(K_full);

            F_unit = zeros(n_dofs_full, 1);
            F_unit(obj.nl_dof) = 1;
            F_unit_c = obj.Assembly.constrain_vector(F_unit);

            % 1. CALCULATION OF THE FIRST LINEAR MODE (Phi_1)
            [Phi_all, D] = eigs(Kc, Mc, 1, 'smallestabs');
            % Expand the mode shape vector to the full space
            Phi_full = obj.Assembly.unconstrain_vector(Phi_all(:, 1));

            % 2. CALCULATION WITH CONTACT (k active)
            K_spring_full_on = sparse(n_dofs_full, n_dofs_full);
            K_spring_full_on(obj.nl_dof, obj.nl_dof) = original_k;
            K_spring_c_on = obj.Assembly.constrain_matrix(K_spring_full_on);

            Kc_closed_on = Kc + K_spring_c_on;
            Psi_c_contact = Kc_closed_on \ F_unit_c;

            % Expansion to the full space
            Psi_full_contact = obj.Assembly.unconstrain_vector(Psi_c_contact);

            % 3. CALCULATION WITHOUT CONTACT (k = 0)
            Psi_c_no_contact = Kc \ F_unit_c;

            % Expansion to the full space
            Psi_full_no_contact = obj.Assembly.unconstrain_vector(Psi_c_no_contact);

            % 4. RETRIEVE COORDINATES AND VERTICAL DISPLACEMENTS (Uy)
            nodes_coord = obj.Assembly.Mesh.nodes;
            x_coords = nodes_coord(:, 1);

            % Assuming 3 degrees of freedom per node (Ux=1, Uy=2, RotZ=3)
            uy_indices = 2:3:n_dofs_full;

            uy_contact = Psi_full_contact(uy_indices);
            uy_no_contact = Psi_full_no_contact(uy_indices);
            uy_mode1 = Phi_full(uy_indices);

            % --- Scale the mode shape to match the peak of the "No Contact" displacement ---
            % This makes the shape comparison visually clear and meaningful
            peak_no_contact = max(abs(uy_no_contact));
            [~, max_idx] = max(abs(uy_mode1));
            scale_factor = (peak_no_contact / uy_mode1(max_idx)) * sign(uy_no_contact(max_idx));
            uy_mode1_scaled = uy_mode1 * scale_factor;

            %% 5. PLOT WITH DUAL Y-AXIS (yyaxis)
            figure('Name', 'Milman-Chu & Linear Mode Comparison', 'Color', 'w');
            hold on;
            grid on;

            % --- Left axis (Blue): Displacement without contact and 1st Linear Mode ---
            yyaxis left
            plot(x_coords, uy_no_contact, '-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
            plot(x_coords, uy_mode1_scaled, '--g', 'LineWidth', 2); % 1st Mode in dashed Green
            ylabel('Displacement Without Contact & Mode 1 [m] (Left Axis)', 'Color', 'k');
            ax = gca;
            ax.YColor = 'b';

            % --- Right axis (Red): Displacement with contact (highly rigid) ---
            yyaxis right
            plot(x_coords, uy_contact, '-^r', 'LineWidth', 2, 'MarkerFaceColor', 'r');
            ylabel('Displacement With Contact [m] (Right Axis)', 'Color', 'r');
            ax.YColor = 'r';

            title('Comparison of Milman-Chu Vectors (\Psi_c) and 1^{st} Linear Mode', 'FontSize', 12);
            xlabel('Beam X Coordinate [m]');

            % Identify the coordinate of the nonlinear node
            nl_node_idx = find(abs(nodes_coord(:,1) - max(x_coords)/2) < 1e-6, 1);
            if ~isempty(nl_node_idx)
                x_contact = nodes_coord(nl_node_idx, 1);
                xline(x_contact, '--k', 'Contact (L/2)', 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.5);
            end

            legend({'Without Contact (k = 0)', '1^{st} Linear Mode (Scaled)', ...
                sprintf('With Contact (k = %.1f N/m)', original_k)}, ...
                'Location', 'southoutside', 'Orientation', 'horizontal');

            hold off;
        end
        function display_rom_frequencies(obj)
            % DISPLAY_ROM_FREQUENCIES 
            if isempty(obj.P)
                error('La matrice di proiezione P è vuota. Esegui build() prima di chiamare questo metodo.');
            end

            M_full = obj.Assembly.mass_matrix();
            K_full = obj.Assembly.stiffness_matrix();

            M_r = obj.P' * M_full * obj.P;
            K_r = obj.P' * K_full * obj.P;

            % Simmetrizzazione numerica per prevenire instabilità in eig()
            % M_r = (M_r + M_r') / 2;
            % K_r = (K_r + K_r') / 2;

            % Compute eigenfrequenies
            [~, D] = eig(full(K_r), full(M_r));
            omega2 = diag(D);

            % Pulizia di eventuali artefatti numerici macchina
            % omega2(omega2 < 0) = 0;

            freqs_Hz = sort(sqrt(omega2) / (2*pi));

            % 4. Stampa formattata su Command Window
            fprintf('\n=======================================================\n');
            fprintf('        ROM NATURAL FREQUENCIES    \n');
            fprintf('=======================================================\n');
            for i = 1:length(freqs_Hz)
                if obj.include_MC && i == length(freqs_Hz)
                    % L'ultima frequenza è quella artificiale generata dal vettore ortogonalizzato
                    fprintf(' Modo %-2d (Alta Freq. da MC) : %10.2f Hz\n', i, freqs_Hz(i));
                else
                    fprintf(' Modo %-2d                     : %10.2f Hz\n', i, freqs_Hz(i));
                end
            end
            fprintf('=======================================================\n\n');
        end
    end
end