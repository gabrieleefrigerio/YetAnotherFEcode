clear; close all; clc;

%% CARICAMENTO DATI FEM
[filename, pathname] = uigetfile('*.mat', 'Select the FEM mesh file');
if isequal(filename,0)
    error('Nessun file selezionato.');
end

data     = load(fullfile(pathname, filename));
nodes    = data.nodes;
elements = double(data.elements);    % converte da int64 a double

fprintf('File caricato: %s\n', filename);
fprintf('Numero di nodi: %d\n', size(nodes,1));
fprintf('Numero di elementi: %d\n', size(elements,1));

%% DEFINIZIONE MATERIALE
E   = 168e3;
rho = 2.33e-15;
nu  = 0.23;
thickness = 30;

myMaterial = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);

%% COSTRUZIONE DELLA MESH
myMesh = Mesh(nodes);

% Usa elemento wedge a 15 nodi (già incluso in YAFEC)
myElementConstructor = @() Wed15Element(myMaterial);

% Crea la tabella degli elementi
myMesh.create_elements_table(elements, myElementConstructor);

%% PLOT MESH (wireframe tipo YAFEC)
figure;
h = PlotMesh(nodes, elements);
legend('Mesh');
axis equal; axis on;
xlabel('X [μm]'); ylabel('Y [μm]'); zlabel('Z [μm]');
title('3D FEM Mesh');
view(25,35);
% per top view:
% view(0,90);

%% SELEZIONE NODI DA VINCOLARE
tol = 1e-4; % tolleranza numerica per confrontare le coordinate z
fixedNodes = find( abs(nodes(:,3) - 52.3) < tol | abs(nodes(:,3) + 2.18) < tol );
fprintf('Numero di nodi vincolati: %d\n', numel(fixedNodes));

%% APPLICAZIONE BOUNDARY CONDITIONS
% Vincolo totale: blocca spostamenti in x, y, z
myMesh.set_essential_boundary_condition(fixedNodes, 1:3, 0);

%% VISUALIZZAZIONE
figure;
PlotMesh(nodes, elements);
hold on;
plot3(nodes(fixedNodes,1), nodes(fixedNodes,2), nodes(fixedNodes,3), ...
      'ro', 'MarkerSize', 5, 'DisplayName', 'Fixed nodes');
legend('Mesh','Fixed nodes');
title('Nodi vincolati');
axis equal; axis on;
view(25,35);

%% ASSEMBLY
myAssembly = Assembly(myMesh);

M  = myAssembly.mass_matrix();
u0 = zeros(myMesh.nDOFs,1);
[K,~] = myAssembly.tangent_stiffness_and_force(u0);

C = sparse(size(M,1), size(M,2));

% Applica i vincoli
Kc = myAssembly.constrain_matrix(K);
Mc = myAssembly.constrain_matrix(M);
Cc = myAssembly.constrain_matrix(C);

% Salva (opzionale)
myAssembly.DATA.Kc = Kc;
myAssembly.DATA.Mc = Mc;
myAssembly.DATA.Cc = Cc;

%% MODI DI VIBRARE
n_VMs = 50;                 % numero di modi
opts  = struct('issym',true,'isreal',true,'tol',1e-8,'maxit',1e3);
[Vc,om] = eigs(Kc, Mc, n_VMs, 'SM', opts);   % Kc * v = om * Mc * v

f0 = sqrt(diag(om))/(2*pi);                  % frequenze [Hz]
[ f0, ind ] = sort(f0);
Vc = Vc(:,ind);

% Normalizza i modi
for ii = 1:n_VMs
    Vc(:,ii) = Vc(:,ii) / max( sqrt(sum(reshape(Vc(:,ii),3,[]).^2,1)) );
end

% Rimuovi i vincoli (spazio pieno)
V0 = myAssembly.unconstrain_vector(Vc);

%% PLOT DEI PRIMI n MODI
mod = 1;                    % scegli il modo
for i=1:mod
    figure('units','normalized','position',[.2 .1 .6 .8]);
    % PlotMesh(nodes, elements, 0); serve a plottare la struttura indeformata
    % campo spostamenti 3D (Nx3)
    v1 = reshape(V0(:,i), 3, []).';
    PlotFieldonDeformedMesh(nodes, elements, v1, 'factor', 10);

    title(['\Phi_' num2str(i) '  f = ' num2str(f0(i),3) ' Hz']);
    view(25,35); axis equal; axis on;
end

% %% SELEZIONE GRAFICA DEI NODI D'IMPATTO
% figure;
% PlotMesh(nodes, elements);
% axis equal; axis on;
% view(25,35);
% title({'Ruota e zooma la vista','Poi premi INVIO per entrare in modalità selezione'});
% rotate3d on;
%
% disp('Ruota e zooma la figura per trovare la zona desiderata.');
% disp('Quando sei pronto a selezionare, premi INVIO.');
% pause;   % attende che premi Invio
%
% rotate3d off;        % blocca la vista attuale per selezionare
% title('Clicca sui nodi da selezionare (premi INVIO per terminare)');
% disp('Ora clicca sui nodi da selezionare. Premi INVIO per terminare.');
%
% [x_sel, y_sel] = ginput();   % clic sui nodi, poi Invio
%
% % Trova i nodi più vicini ai clic
% impactNodes = [];
% for i = 1:length(x_sel)
%     [~, idx] = min((nodes(:,1)-x_sel(i)).^2 + (nodes(:,2)-y_sel(i)).^2);
%     impactNodes = [impactNodes; idx];
% end
% impactNodes = unique(impactNodes);
%
% fprintf('Nodi selezionati per impatto: %d\n', numel(impactNodes));
%
% hold on;
% plot3(nodes(impactNodes,1), nodes(impactNodes,2), nodes(impactNodes,3), ...
%       'ro', 'MarkerSize',6, 'LineWidth',1.5, 'DisplayName','Impact nodes');
% legend('Mesh','Impact nodes');
% title('Nodi d''impatto selezionati');

%% SELEZIONE NODI DI IMPATTO PER COORDINATE

tol = 1e-3;   % tolleranza numerica sulle coordinate

% Esempi di range:
x_min = -202.875;   x_max = -172.875;    % seleziona una fascia in x
y_min = -115;       y_max = -105;        % seleziona una fascia in y
z_min = 30.3;       z_max = 30.3;        % seleziona una fascia in z

impactNodes_z = find( ...
    nodes(:,1) >= x_min - tol & nodes(:,1) <= x_max + tol & ...
    nodes(:,2) >= y_min - tol & nodes(:,2) <= y_max + tol & ...
    nodes(:,3) >= z_min - tol & nodes(:,3) <= z_max + tol );

fprintf('Nodi selezionati per impatto: %d\n', numel(impactNodes_z));
ImpactDOF   = myMesh.get_DOF_from_nodeIDs(impactNodes_z);
ImpactDOF_z = ImpactDOF(:,3);
idx_ImpactDOF_constr = map_full_to_constrained(myAssembly, ImpactDOF_z);

% plot nodi selezionati
figure;
PlotMesh(nodes, elements);
hold on;
plot3(nodes(impactNodes_z,1), nodes(impactNodes_z,2), nodes(impactNodes_z,3), ...
      'bo', 'MarkerSize',3, 'LineWidth',1.5, 'DisplayName','Impact nodes');
legend('Mesh','Impact nodes');
title('Nodi che possono impattare');
axis equal; axis on; view(25,35);

%% VARIABILI DEL MURO
wall_distance_z = 2;  % [μm] distanza del muro dai nodi che impattano
k_tilde = 1000;       % contact stiffness

%% INITIAL CONDITIONS: deform the structure along a certain mode.
idx_mode = 1;
scale    = 10;
q0   = Vc(:,idx_mode) * scale;
qd0  = zeros(size(q0, 1),1);    % gia vincolati
qdd0 = zeros(size(q0, 1),1);

%% INTEGRATION TIME
T    = 1/f0(idx_mode);   % period of the choosen mode
h    = T/500;            % time step choosen to have 1000 points/samples in a period
tmax = 0.8*T;            % total simulation time is three periods

%% FULL-ORDER DYNAMICS
tic
residual = @(q,qd,qdd,t) residual_linear_impact(q,qd,qdd,t, myAssembly, wall_distance_z, k_tilde, ImpactDOF_z);  % nonlinear
TI_NL_alpha = GeneralizedAlpha('timestep',h,'rho_inf',0.7,'linear',false);
TI_NL_alpha.Integrate(q0,qd0,qdd0,tmax,residual);
time_full_order = toc;

% PLOT
TI_NL_alpha.Solution.u = myAssembly.unconstrain_vector(TI_NL_alpha.Solution.q);
Z_displacement = TI_NL_alpha.Solution.u(ImpactDOF_z,:) + nodes(impactNodes_z,3);
figure; hold on; grid on
plot(TI_NL_alpha.Solution.time, Z_displacement, 'LineWidth', 1)
yline(wall_distance_z + z_min, 'k', 'Wall')
xlabel('Time [s]'); ylabel('Z position [μm]');
title('Full-order dynamics with generalized-α')
legend(compose('node %d', impactNodes_z), 'Location', 'bestoutside')

%% ANIMAZIONE DEL MOTO FULL-ORDER
figure('units','normalized','position',[.2 .1 .6 .8]);
axis equal; axis on; grid on;
xlabel('X [μm]'); ylabel('Y [μm]'); zlabel('Z [μm]');
title('Full-order dynamic response');
view(25,35);
hold on;

% Parametri animazione
factor = 1;        % fattore di amplificazione dello spostamento
skip   = 5;        % mostra 1 frame ogni 10 step di integrazione
nFrames = length(1:skip:length(TI_NL_alpha.Solution.time));

% Posizione iniziale
u_full = myAssembly.unconstrain_vector(TI_NL_alpha.Solution.q(:,1));
u_disp = reshape(u_full, 3, []).';
hPlot  = PlotFieldonDeformedMesh(nodes, elements, factor*u_disp, 'factor', 1);
drawnow;

% Aggiornamento frame-by-frame
for k = 1:skip:length(TI_NL_alpha.Solution.time)
    cla;   % pulisce gli oggetti grafici
    qk     = TI_NL_alpha.Solution.q(:,k);
    u_full = myAssembly.unconstrain_vector(qk);
    u_disp = reshape(u_full, 3, []).';
    PlotFieldonDeformedMesh(nodes, elements, factor*u_disp, 'factor', 1);

    % Disegna il muro di contatto
    hold on
    z_wall = nodes(impactNodes_z(1),3) + wall_distance_z;
    fill3([min(nodes(:,1)) max(nodes(:,1)) max(nodes(:,1)) min(nodes(:,1))], ...
          [min(nodes(:,2)) min(nodes(:,2)) max(nodes(:,2)) max(nodes(:,2))], ...
          [z_wall z_wall z_wall z_wall], ...
          [0.8 0.8 0.8], 'FaceAlpha',0.5, 'EdgeColor','none');

    title(sprintf('Time = %.2e s', TI_NL_alpha.Solution.time(k)));
    axis equal; axis tight; view(35, 25);
    drawnow;
end

%% --- NOTE SU FULL E CONSTRAINED DOFs ---------------------------------
%
%  In YAFEC ogni modello esiste in due spazi:
%
%   1) FULL SPACE:
%        - Include TUTTI i DOF del modello, anche quelli vincolati.
%        - Dimensione: n_full = myMesh.nDOFs
%        - Vettori in questo spazio (es. u_full) contengono anche
%          i DOF bloccati a zero dalle condizioni al contorno.
%
%   2) CONSTRAINED SPACE:
%        - Include SOLO i DOF liberi (cioè non vincolati).
%        - Dimensione: n_c = size(Kc,1)
%        - Vettori in questo spazio (es. q) rappresentano
%          solo i gradi di libertà effettivamente incogniti.
%
%  La corrispondenza tra i due spazi è data da:
%       myAssembly.Mesh.EBC.unconstrainedDOFs
%  che contiene gli INDICI FULL dei DOF liberi.
%
%  Le due funzioni chiave per il passaggio di spazio sono:
%       q      = myAssembly.constrain_vector(u_full);     % FULL -> CONSTRAINED
%       u_full = myAssembly.unconstrain_vector(q);        % CONSTRAINED -> FULL
%
%  ---------------------------------------------------------------------
%  DOF DI INTERFACCIA E INTERNI
%
%  Gli indici dei DOF di interfaccia (boundary) e interni (internal)
%  cambiano a seconda dello spazio in cui stiamo lavorando:
%
%   - boundaryDOFs_full : indici nel vettore FULL (u_full)
%   - boundaryDOFs      : indici nel vettore CONSTRAINED (q)
%
%  Il mapping tra i due si ottiene con:
%       boundaryDOFs      = map_full_to_constrained(myAssembly, boundaryDOFs_full);
%       boundaryDOFs_full = myAssembly.Mesh.EBC.unconstrainedDOFs(boundaryDOFs);
%
%  Nel metodo di Rubin si lavora SEMPRE nello spazio CONSTRAINED,
%  ma per plottare o visualizzare i risultati bisogna tornare allo spazio FULL.
%
%  ---------------------------------------------------------------------

%% SCELGO DOFs INTERNI E DI INTERFACCIA

% % DOF di contatto nel FULL
% boundaryDOFs_full = ImpactDOF_z;
%
% % Mappa FULL → CONSTRAINED (rimuove i DOF vincolati a terra)
% boundaryDOFs = map_full_to_constrained(myAssembly, boundaryDOFs_full);
% boundaryDOFs = unique(boundaryDOFs(:), 'stable');
%
% % DOF interni = tutti i DOF liberi esclusi i boundary
% n_c   = size(myAssembly.DATA.Kc,1);
% all_c = (1:n_c).';
% InternalDOFs = setdiff(all_c, boundaryDOFs, 'stable');
%
% % Ordine utile [boundary; internal]
% ord_c = [boundaryDOFs; InternalDOFs];

% io farei cosi':
% Boundary DOFs
boundaryDOFs = ImpactDOF_z;
idx_boundaryDOFs_full   = boundaryDOFs;
idx_boundaryDOFs_constr = map_full_to_constrained(myAssembly, idx_boundaryDOFs_full);

% Internal DOFs
idx_internalDOFs_full   = setdiff(1:myMesh.nDOFs, idx_boundaryDOFs_full, 'stable');
idx_internalDOFs_constr = map_full_to_constrained(myAssembly, idx_internalDOFs_full);

% Il vettore delle coordinate fisiche/reali e constrained sara' diviso in [boundary; internal]
ord_constr = [idx_boundaryDOFs_constr; idx_internalDOFs_constr];

%% STEP 2 — Modi free-interface Φ̄ (da Vc) e loro separazione b/i

% % 1) prendi i primi n_modes_free modi dalla soluzione vincolata
% n_modes_free = 50;                  % scegli tu (20–100 tipico)
% Phi_bar_c    = Vc(:, 1:n_modes_free);   % [n_c x n_modes_free], spazio CONSTRAINED
%
% % 2) porta i modi nello spazio FULL per poterli poi separare sui DOF desiderati
% Phi_bar_full = myAssembly.unconstrain_vector(Phi_bar_c);   % [n_full x n_modes_free]
%
% % 3) ricostruisci un selettore FULL con ordine [ub; ui]
% %    ub, ui sono nel CONSTRAINED → servono i FULL id dei DOF liberi
% u_full_free = myAssembly.Mesh.EBC.unconstrainedDOFs(:);    % FULL ids dei DOF liberi
%
% boundaryDOFs_full = u_full_free(boundaryDOFs);   % FULL ids di ub
% InternalDOFs_full = u_full_free(InternalDOFs);   % FULL ids di ui
%
% % 4) separa i blocchi Φ̄_b e Φ̄_i secondo l'ordine standard Rubin [ub; ui]
% Phi_bar_b = Phi_bar_full(boundaryDOFs_full, :);   % [n_b x n_modes_free]
% Phi_bar_i = Phi_bar_full(InternalDOFs_full, :);   % [n_i x n_modes_free]

% io farei cosi':
n_modes_free = 50;   % scelgo il numero id modes che voglio nella base ridotta

% divido le modes nei temrini associati ai DOFs interni e boundary

Phi_bar_full_internal = V0(idx_internalDOFs_full, 1:n_modes_free);
Phi_bar_full_boundary = V0(idx_boundaryDOFs_full, 1:n_modes_free);

Phi_bar_constr_internal = Vc(idx_internalDOFs_constr, 1:n_modes_free);
Phi_bar_constr_boundary = Vc(idx_boundaryDOFs_constr, 1:n_modes_free);

%% STEP 3 — RAMs Ψ_r (constrained)

% b   = boundaryDOFs(:);          % DOF boundary (constrained)
% i   = InternalDOFs(:);          % DOF interni (constrained)
% n_b = numel(b);
%
% Kc = myAssembly.DATA.Kc;
%
% % Partizione
% K_bb = Kc(b,b);
% K_bi = Kc(b,i);
% K_ib = Kc(i,b);
% K_ii = Kc(i,i);
%
% % RAMs: risolvi K_ii * Psi_i = -K_ib
% Psi_i = -(K_ii \ K_ib);         % [n_i x n_b]
%
% % Assembla Ψ_r nello spazio CONSTRAINED, ordine [ub; ui]
% Psi_r_constr      = sparse(size(Kc,1), n_b);
% Psi_r_constr(b,:) = +speye(n_b);
% Psi_r_constr(i,:) = Psi_i;
%
% % (opzionale) porta in FULL per plotting
% Psi_r_full = myAssembly.unconstrain_vector(Psi_r_constr);
%
% % (opzionale) plot dei primi RAM
% n_plot = min(3,n_b);
% for k = 1:n_plot
%     v = reshape(Psi_r_full(:,k), 3, []).';      % 3D: [Nx3]
%     figure; PlotFieldonDeformedMesh(nodes, elements, v, 'factor', 10);
%     title(sprintf('Residual Attachment Mode %d', k));
% end
%
% % Check rapido (facoltativo): ||K*[Ψ_r] - [0; -K_ib]|| sulla partizione
% res_check = norm(K_ii*Psi_i + K_ib, 'fro')/max(1,norm(K_ib,'fro'));
% fprintf('RAM residual (norm. Fro): %.2e\n', res_check);

% io farei cosi':
% decido di calcolare le RAMs constrained
n_b = numel(idx_boundaryDOFs_constr);

% partiziono la matrice di rigidezza constrained
K_bb = Kc(idx_boundaryDOFs_constr, idx_boundaryDOFs_constr);
K_bi = Kc(idx_boundaryDOFs_constr, idx_internalDOFs_constr);
K_ib = Kc(idx_internalDOFs_constr, idx_boundaryDOFs_constr);
K_ii = Kc(idx_internalDOFs_constr, idx_internalDOFs_constr);

% Schur complement (solo se serve per controllo)
S = K_bb - K_bi * (K_ii \ K_ib);

% Inizializzo il vettore delle residual attachment modes
Psi_r_boundary_constr = zeros(n_b, n_b);
Psi_r_internal_constr = zeros(numel(idx_internalDOFs_constr), n_b);

tic
for k = 1:n_b
    e_k    = zeros(n_b,1);
    e_k(k) = -1;                          % unit load su ogni DOF di interfaccia
    psi_b  = S \ e_k;                     % risposta ai DOF di interfaccia
    psi_i  = -(K_ii \ (K_ib * psi_b));    % risposta ai DOF interni
    Psi_r_boundary_constr(:,k) = psi_b;
    Psi_r_internal_constr(:,k) = psi_i;
end

% assemblaggio finale
Psi_r_constr = zeros(size(Kc,1), n_b);
Psi_r_constr(idx_internalDOFs_constr,:) = Psi_r_internal_constr;
Psi_r_constr(idx_boundaryDOFs_constr,:) = Psi_r_boundary_constr;
time_RAMs_computation = toc;

% Riporta nello spazio full per fare il plot
Psi_r_full          = myAssembly.unconstrain_vector(Psi_r_constr);
Psi_r_internal_full = Psi_r_full(idx_internalDOFs_constr,:);
Psi_r_boundary_full = Psi_r_full(idx_boundaryDOFs_constr,:);

% PLOT
n_plot = min(4, size(Psi_r_full,2));   % plotta i primi 4 al massimo

for mod = 1:n_plot
    figure('units','normalized','position',[.2 .1 .6 .8])
    v1 = reshape(Psi_r_full(:,mod),3,[]).';
    PlotFieldonDeformedMesh(nodes,elements,v1,'factor',100);
    title(['Residual Attachment Mode ' num2str(mod)])
    drawnow
end

%% STEP 4 — Costruzione e riduzione Rubin (uso dei dati già calcolati)

% % Ordine coerente [ub; ui]
% ord_c = [boundaryDOFs; InternalDOFs];
%
% % Base Rubin già definita
% T_R_c = [ speye(numel(boundaryDOFs))   sparse(numel(boundaryDOFs), n_modes_free);
%           Psi_i                        Phi_bar_i ];
%
% % Matrici ridotte
% % Mc = myAssembly.DATA.Mc;  Kc = myAssembly.DATA.Kc;  Cc = myAssembly.DATA.Cc;
% Mr = T_R_c.' * Mc(ord_c,ord_c) * T_R_c;
% Kr = T_R_c.' * Kc(ord_c,ord_c) * T_R_c;
% Cr = T_R_c.' * Cc(ord_c,ord_c) * T_R_c;
%
% % ---- Stampa riepilogo ----
% fprintf('\nProiezione completata:\n');
% fprintf('  Mr: %d x %d\n', size(Mr,1), size(Mr,2));
% fprintf('  Kr: %d x %d\n', size(Kr,1), size(Kr,2));
% fprintf('  Cr: %d x %d\n', size(Cr,1), size(Cr,2));
%
% Dim_base_reduced = n_modes_free + n_b;
% fprintf('Numero totale DOF ridotti: %d ( %d modi interni free + %d residual attachment modes )\n', ...
%          Dim_base_reduced, n_modes_free, n_b);

% io farei cosi':
% la cosa importante nel costruire la T e' la coerenza. Ovvero per i
% termini che formeranno la T devo usarli tutti full o tutti constrained.
% Ad esempio io la faccio constrained

T_rubin_constr = [     speye(n_b)                sparse(n_b, n_modes_free);
                   Psi_r_internal_constr         Phi_bar_constr_internal    ];

% Matrici ridotte
Mr = T_rubin_constr.' * Mc(ord_constr,ord_constr) * T_rubin_constr;
Kr = T_rubin_constr.' * Kc(ord_constr,ord_constr) * T_rubin_constr;
Cr = T_rubin_constr.' * Cc(ord_constr,ord_constr) * T_rubin_constr;

% ---- Stampa riepilogo ----
fprintf('\nProiezione completata:\n');
fprintf('  Mr: %d x %d\n', size(Mr,1), size(Mr,2));
fprintf('  Kr: %d x %d\n', size(Kr,1), size(Kr,2));
fprintf('  Cr: %d x %d\n', size(Cr,1), size(Cr,2));

Dim_base_reduced = n_modes_free + n_b;
fprintf('Numero totale DOF ridotti: %d ( %d modi interni free + %d residual attachment modes )\n', ...
         Dim_base_reduced, n_modes_free, n_b);

%% REDUCED-ORDER DYNAMICS  (Rubin)

% Condizioni iniziali nel ridotto
q0_ord  = q0(ord_constr);
qd0_ord = qd0(ord_constr);
qr0     = T_rubin_constr \ q0_ord;
qrd0    = T_rubin_constr \ qd0_ord;
qrdd0   = zeros(size(qr0));

tic
% Integrazione con Generalized-α
residual_RO = @(qr,qrd,qrdd,t) residual_rubin_contact( ...
    qr,qrd,qrdd,t, Mr, Cr, Kr, wall_distance_z, k_tilde, n_b);

TI_RED = GeneralizedAlpha('timestep',h,'rho_inf',0.7,'linear',false);
TI_RED.Integrate(qr0, qrd0, qrdd0, tmax, residual_RO);
time_reduced_order = toc;

% Ricostruzione in spazio vincolato e poi completo
q_c_red = T_rubin_constr * TI_RED.Solution.q;          % [ub; ui]
Qc      = zeros(size(Kc,1), size(q_c_red,2));          % spazio constrained globale
Qc(ord_constr,:) = q_c_red;
U_red   = myAssembly.unconstrain_vector(Qc);           % spazio full

% Spostamento dei nodi d'impatto
Z_red = U_red(ImpactDOF_z,:) + nodes(impactNodes_z,3);

% --- Confronto full-order vs reduced-order ---
figure; hold on; grid on
plot(TI_NL_alpha.Solution.time, Z_displacement, 'LineWidth', 1.2, 'DisplayName', 'Full order')
plot(TI_RED.Solution.time,      Z_red,    '--', 'LineWidth', 1.2, 'DisplayName', 'Reduced (Rubin)')
yline(wall_distance_z + z_min, 'k', 'Wall')
xlabel('Time [s]'); ylabel('Z position [μm]');
title('Full-order vs Reduced-order (Rubin)')
legend('Location','bestoutside')


%% ==================== FUNZIONI ====================

function [r, drdqdd, drdqd, drdq, c0] = residual_linear_impact(q,qd,qdd, ...
    t, myAssembly, wall_distance_z, k_tilde, ImpactDOF_z)
% NOTE: q,qd,qdd sono già VINCOLATI.
% Mc,Cc,Kc sono già vincolate. Fext_full è in DOF completi.

% ricostruisci spostamenti non vincolati
% q contiene gli spostamenti dei nodi non vincolati,
% ed e' quello che ci serve per l'integrazione
% u invece contiene la posizione di tutti i nodi,
% inclusi anche quelli bloccati a terra

u_full      = myAssembly.unconstrain_vector(q);
penetration = u_full(ImpactDOF_z) - wall_distance_z;

F_contact = zeros(myAssembly.Mesh.nDOFs, 1);

% applica forza elastica solo dove c'è penetrazione
hit = penetration > 0;
if any(hit)
    F_contact(ImpactDOF_z(hit)) = + k_tilde * penetration(hit);   % ho dovuto cambiare il segno qui ma non ho capito perche
end

% forza di contatto vincolata
F_contact_c = myAssembly.constrain_vector(F_contact);

% Residual is computed according to the formula above:
F_inertial = myAssembly.DATA.Mc * qdd;
F_damping  = myAssembly.DATA.Cc * qd;
F_elastic  = myAssembly.DATA.Kc * q;
r = F_inertial + F_damping + F_elastic - F_contact_c ;
drdqdd = myAssembly.DATA.Mc;
drdqd  = myAssembly.DATA.Cc;
drdq   = myAssembly.DATA.Kc;
c0 = norm(F_inertial) + norm(F_damping) + norm(F_elastic) + norm(F_contact_c);
end

function idx_c = map_full_to_constrained(myAssembly, idx_full)
    % Converte indici FULL in CONSTRAINED (solo DOF liberi)
    unconstrained_full = myAssembly.Mesh.EBC.unconstrainedDOFs(:);
    [tf, loc] = ismember(idx_full(:), unconstrained_full);
    idx_c = loc(tf);   % scarta i DOF vincolati
end

function [r, drdqdd, drdqd, drdq, c0] = residual_rubin_contact( ...
    qr, qrd, qrdd, t, Mr, Cr, Kr, wall_distance_z, k_tilde, n_b)

    % --- DOF di interfaccia (primi n_b del ridotto) ---
    u_b = qr(1:n_b);

    % --- Penetrazione rispetto al muro ---
    % Tutti i nodi di contatto hanno stesso z_min di riferimento
    penetration = u_b - wall_distance_z;

    % --- Forza di contatto (solo dove c'è penetrazione) ---
    hit = penetration > 0;
    F_r = zeros(size(qr));
    if any(hit)
        F_r(hit) = + k_tilde * penetration(hit);   % stesso segno del full-order
    end

    % --- Residuo ridotto ---
    r       = Mr*qrdd + Cr*qrd + Kr*qr - F_r;
    drdqdd  = Mr;
    drdqd   = Cr;
    drdq    = Kr;
    c0      = norm(Mr*qrdd) + norm(Cr*qrd) + norm(Kr*qr) + norm(F_r);
end