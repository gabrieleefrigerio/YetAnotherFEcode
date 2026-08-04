%% =====================================================================
%  SCRIPT: plot_results.m
%  Obiettivo: Caricare e plottare i risultati delle simulazioni salvate
%  =====================================================================
clearvars -except DummyStruct Mc Kc Cc; % Puliamo il workspace mantenendo le matrici pesanti se esistono
clc;

fprintf('=== POST-PROCESSING RISULTATI ===\n');

% --- 1. Selezione Interattiva della Cartella ---
% Apre un popup per farti scegliere la cartella dentro "results"
starting_folder = fullfile(pwd, 'results');
if ~exist(starting_folder, 'dir')
    error('La cartella "results" non esiste. Nessun salvataggio trovato.');
end

selected_dir = uigetdir(starting_folder, 'Seleziona la cartella con i risultati da plottare');

if selected_dir == 0
    fprintf('Operazione annullata dall''utente.\n');
    return;
end
fprintf('Cartella selezionata: %s\n', selected_dir);

% --- 2. Caricamento dei File ---
% Inizializziamo dei flag per sapere cosa abbiamo trovato
has_fom_newmark = false;
has_fom_ode = false;
has_rom_rubin = false;

if exist(fullfile(selected_dir, 'FOM_Newmark.mat'), 'file')
    load(fullfile(selected_dir, 'FOM_Newmark.mat'));
    has_fom_newmark = true;
    fprintf('Caricato: FOM_Newmark.mat\n');
end

if exist(fullfile(selected_dir, 'FOM_ODE15s.mat'), 'file')
    load(fullfile(selected_dir, 'FOM_ODE15s.mat'));
    has_fom_ode = true;
    fprintf('Caricato: FOM_ODE15s.mat\n');
end

if exist(fullfile(selected_dir, 'ROM_Rubin.mat'), 'file')
    load(fullfile(selected_dir, 'ROM_Rubin.mat'));
    has_rom_rubin = true;
    fprintf('Caricato: ROM_Rubin.mat\n');
end

if ~(has_fom_newmark || has_fom_ode || has_rom_rubin)
    error('Nessun file .mat valido trovato nella cartella selezionata.');
end

% --- 3. Recupero Parametri Strutturali (Failsafe) ---
% Se hai appena aperto MATLAB, DummyStruct non esiste. Dobbiamo rigenerarla
% per sapere quali sono i nodi da plottare.
if ~exist('DummyStruct', 'var')
    fprintf('DummyStruct non trovata nel workspace. Rigenerazione geometria...\n');
    DummyStruct = DummyStructure('elementType', 'QUAD4');
    DummyStruct.build();
    Kc = DummyStruct.AssemblyObj.constrain_matrix(DummyStruct.K);
end

% Parametri fisici (Se avevi cambiato questi valori nel main, aggiornali anche qui!)
gap_wall = 1.5e-6; 
k_contact = max(diag(Kc)) * 10; 

% Identificazione dei nodi per i plot
nodo_spia = length(DummyStruct.nodes);
nodi_parete = find(abs(DummyStruct.nodes(:,1) - max(DummyStruct.nodes(:,1))) < 1e-9);

dof_spia_global = (nodo_spia - 1) * DummyStruct.MeshObj.nDOFPerNode + 1;
dofs_parete_global = (nodi_parete - 1) * DummyStruct.MeshObj.nDOFPerNode + 1;

% --- 4. Generazione Grafico Cinematica (Spostamento) ---
figure('Name', ['Cinematica - ' nome_cartella(selected_dir)], 'Color', 'w', 'Units', 'normalized', 'Position', [0.1 0.4 0.4 0.4]);
hold on;

if has_fom_newmark
    plot(t_fom_newmark, y_fom_newmark(dof_spia_global, :)*1e6, 'b', 'LineWidth', 2, 'DisplayName', 'FOM (Newmark)');
end
if has_fom_ode
    plot(t_fom_ode, y_fom_ode(dof_spia_global, :)*1e6, 'g-.', 'LineWidth', 1.5, 'DisplayName', 'FOM (ode15s)');
end
if has_rom_rubin
    plot(t_rubin, y_rom_rubin(dof_spia_global, :)*1e6, 'r--', 'LineWidth', 1.5, 'DisplayName', 'ROM (Rubin)');
end

yline(gap_wall*1e6, 'k:', 'LineWidth', 2, 'DisplayName', 'Muro (Gap)');
title(sprintf('Transitorio di Impatto - Nodo #%d', nodo_spia));
xlabel('Tempo [s]'); ylabel('Spostamento X [\mum]');
legend('Location', 'best'); grid on;
hold off;

% --- 5. Generazione Grafico Forza di Contatto ---
figure('Name', ['Forze di Contatto - ' nome_cartella(selected_dir)], 'Color', 'w', 'Units', 'normalized', 'Position', [0.5 0.4 0.4 0.4]);
hold on;

if has_fom_newmark
    pen_fom = max(0, y_fom_newmark(dofs_parete_global, :) - gap_wall);
    forza_fom = sum(k_contact * pen_fom, 1);
    plot(t_fom_newmark, forza_fom*1e6, 'b', 'LineWidth', 2, 'DisplayName', 'FOM (Newmark)');
end

if has_fom_ode
    pen_fom_ode = max(0, y_fom_ode(dofs_parete_global, :) - gap_wall);
    forza_fom_ode = sum(k_contact * pen_fom_ode, 1);
    plot(t_fom_ode, forza_fom_ode*1e6, 'g-.', 'LineWidth', 1.5, 'DisplayName', 'FOM (ode15s)');
end

if has_rom_rubin
    pen_rom = max(0, y_rom_rubin(dofs_parete_global, :) - gap_wall);
    forza_rom = sum(k_contact * pen_rom, 1);
    plot(t_rubin, forza_rom*1e6, 'r--', 'LineWidth', 1.5, 'DisplayName', 'ROM (Rubin)');
end

title('Forza di Contatto Totale');
xlabel('Tempo [s]'); ylabel('Forza [\muN]');
legend('Location', 'best'); grid on;
hold off;

fprintf('--- Plot completati con successo! ---\n');

% =========================================================================
% FUNZIONE HELPER LOCALE
% =========================================================================
function nome = nome_cartella(percorso_completo)
    % Estrae solo il nome finale della cartella (es. '2026-06-02_10-15-30')
    [~, nome, ~] = fileparts(percorso_completo);
end