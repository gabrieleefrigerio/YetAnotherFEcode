%% --- 3. INIZIALIZZAZIONE STRUTTURA BASE ---
fprintf('Setting up the model...\n');
DummyStruct = AbaqusStructure();
DummyStruct.filename = 'DummyStructureAbaqus_V4.inp'; 
DummyStruct.elementType = 'TRI3';           
DummyStruct.build();
max_phi = max(array_linModes);
fprintf('Estrazione modale di %d modi...\n', max_phi);
DummyStruct.compute_eigenmodes(max_phi);
% Extracting Structure's data
Mc = DummyStruct.AssemblyObj.constrain_matrix(DummyStruct.M);
Kc = DummyStruct.AssemblyObj.constrain_matrix(DummyStruct.K);
n_dofs_fom = size(Mc, 1);
k_base = max(diag(Kc)); % Rigidezza base da moltiplicare nel ciclo
% Extracting contact DoFs
contact_dofs = DummyStruct.get_contact_dofs(1); 
contact_nodes = DummyStruct.contact_nodes;
coord_Y_contact_nodes = DummyStruct.nodes(contact_nodes, 2); 
% GdL Global for X and Y
contact_nodes_global_X = (contact_nodes - 1) * DummyStruct.MeshObj.nDOFPerNode + 1;
contact_nodes_global_Y = (contact_nodes - 1) * DummyStruct.MeshObj.nDOFPerNode + 2;
% Forcing
n_dof_per_node = DummyStruct.MeshObj.nDOFPerNode;
dir_vector = zeros(n_dofs_fom, 1);
dir_vector(1:n_dof_per_node:n_dofs_fom) = 1; 
F_spatial_fom = Mc * dir_vector;
F_fom_handle = @(t) F_spatial_fom * shock_amp * sin(pi * t / t_shock) * (t <= t_shock);
% Initial Conditions
q0 = zeros(n_dofs_fom, 1);
qd0 = zeros(n_dofs_fom, 1);


%% --- RISOLUZIONE PROBLEMA STATICO PURO CON SPOSTAMENTO IMPOSTO ---

% 1. Inizializzazione: matrice K vincolata e vettore forze nullo
K_static = Kc; 
F_static = zeros(n_dofs_fom, 1); 

% 2. Definizione dei GdL da controllare e dello spostamento
dofs_imposti = contact_nodes_global_X; % Usiamo i GdL Y dei nodi di contatto
valore_spostamento = -15e-6;            % Inserisci qui lo spostamento desiderato

% 3. Metodo della Penale
% Usiamo una rigidezza fittizia molto grande per "forzare" il nodo
alfa_penale = k_base * 1e8; 

% Aggiungiamo la penale sulla diagonale di K
K_static = K_static + sparse(dofs_imposti, dofs_imposti, alfa_penale, n_dofs_fom, n_dofs_fom);

% Applichiamo la forza fittizia corrispondente
F_static(dofs_imposti) = alfa_penale * valore_spostamento;

% 4. Risoluzione del Sistema
fprintf('Calcolo della soluzione statica...\n');
u_static = K_static \ F_static;

% 5. Calcolo delle Reazioni Vincolari (Opzionale ma utile)
% Moltiplicando la matrice di rigidezza originale per gli spostamenti calcolati,
% otteniamo le forze necessarie per mantenere quello spostamento.
F_reazioni = Kc * u_static;
Forza_totale_necessaria = sum(F_reazioni(dofs_imposti));

fprintf('Spostamento imposto di %g completato.\n', valore_spostamento);
fprintf('Forza di reazione totale sui nodi spostati: %g\n', Forza_totale_necessaria);
% Ad esempio, ingrandendo la deformazione di 10 volte per renderla più visibile
DummyStruct.plot_static_result(u_static,1);

% 
% DummyStruct.plot_mode(1);
% 
% DummyStruct.plot_mode(2);