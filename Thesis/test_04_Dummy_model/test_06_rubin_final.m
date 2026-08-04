%% =====================================================================
%  TEST FINALE: PARAMETRIZZAZIONE DELLA BASE E COSTO ONLINE
% ---------------------------------------------------------------------
%  Cinque configurazioni a SOTTOSPAZIO COSTANTE (Rubin) piu' il riferimento
%  Milman-Chu, che genera lo stesso sottospazio in coordinate diverse.
%
%    A   Rubin grezza   + ramo standard    -> baseline
%    B   Rubin grezza   + ramo proiettato  -> costo del ramo (B vs A)
%    C   Rubin scalata  + ramo proiettato  -> effetto scaling (C vs B)
%    C'  Rubin scalata  + ramo standard    -> configurazione di produzione
%    D   Milman-Chu + Gram-Schmidt         -> riferimento M-ortonormale
%
%  Confronti controllati:
%    B vs A   : solo il ramo cambia            -> costo O(b*r) vs O(b)
%    C vs B   : solo lo scaling cambia         -> effetto condizionamento
%    C' vs C  : solo il ramo cambia, e la AbsTol energetica e' IDENTICA
%               perche' nel ramo proiettato dK_c = k*sum(Pc_c.^2) e le prime
%               b colonne di Pc_c valgono diag(1/d_j), cioe' esattamente il
%               k_run = k/d_j^2 passato al ramo standard.
%    C' vs D  : Rubin vs Milman-Chu, entrambi ben condizionati
%
%  INVARIANTE: il GRE deve essere identico (entro il rumore di roundoff)
%  in A, B, C, C'. Sono cambi base invertibili dello stesso sottospazio.
%
%  Richiede: RomRubinTest.m, TransientSolverOde_V2_test.m, RomMC.m
% =====================================================================
clear; close all; clc;

%% --- 1. INPUTS ---
array_linModes = 200;          % es. [10 50 100 150 200] per lo sweep completo
array_QFactor  = 1000;
array_k_mult   = 10;

use_fixed_abstol = false;      % true -> AbsTol scalare identica a tutte le cfg

impulse_g         = 1e5;
impulse_amp       = impulse_g * 9.81;
impulse_angle_deg = 0;
impulse_sign      = 1;
t_shock = 10e-7;
dt      = 0.5e-8;
tmax    = (1e-3)/2;

gap_T =  5e-6;    gap_B = -1.5e-6;
gap_L = -5e-6;    gap_R =  1.5e-6;

bench_RelTol = 1e-9;

%% --- 2. DIRECTORY ---
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM');
save_dir  = fullfile('results_V2', sprintf('RubinFinal_%s', timestamp));
if ~exist(save_dir, 'dir'), mkdir(save_dir); end
fprintf('Directory creata: %s\n\n', save_dir);

%% --- 3. STRUTTURA BASE ---
fprintf('Setting up the model...\n');
DummyStruct = AbaqusStructure_V2();
DummyStruct.filename    = 'DummyStructureAbaqus_V4.inp';
DummyStruct.elementType = 'TRI3';
DummyStruct.build();

max_phi = max(array_linModes);
fprintf('Estrazione modale di %d modi...\n', max_phi);
DummyStruct.compute_eigenmodes(max_phi);

Mc = DummyStruct.AssemblyObj.constrain_matrix(DummyStruct.M);
Kc = DummyStruct.AssemblyObj.constrain_matrix(DummyStruct.K);
n_dofs_fom = size(Mc, 1);
k_base = max(diag(Kc));

%% --- 3.1 INTERFACCE ---
dofs_T = DummyStruct.get_specific_contact_dofs('T', 2);
dofs_B = DummyStruct.get_specific_contact_dofs('B', 2);
dofs_L = DummyStruct.get_specific_contact_dofs('L', 1);
dofs_R = DummyStruct.get_specific_contact_dofs('R', 1);

% L'ordine fissato qui e' l'unica verita': viene ereditato da contactDofs,
% dalle colonne di Psi_res, dalle coordinate ridotte e da scaleD.
contact_dofs = [dofs_T; dofs_B; dofs_L; dofs_R];
gaps_array   = [repmat(gap_T, length(dofs_T), 1); ...
                repmat(gap_B, length(dofs_B), 1); ...
                repmat(gap_L, length(dofs_L), 1); ...
                repmat(gap_R, length(dofs_R), 1)];
n_bnd = length(contact_dofs);

assert(numel(unique(contact_dofs)) == n_bnd, 'contact_dofs contiene duplicati.');
assert(numel(gaps_array) == n_bnd, 'gaps_array disallineato da contact_dofs.');

Interfaces = struct();
labels     = {'T','B','L','R'};
nodes_list = {DummyStruct.contact_nodes_T, DummyStruct.contact_nodes_B, ...
              DummyStruct.contact_nodes_L, DummyStruct.contact_nodes_R};
for i = 1:numel(labels)
    l = labels{i};  n = nodes_list{i};
    Interfaces.(l).nodes = n;
    if ~isempty(n)
        Interfaces.(l).global_X = (n-1) * DummyStruct.MeshObj.nDOFPerNode + 1;
        Interfaces.(l).global_Y = (n-1) * DummyStruct.MeshObj.nDOFPerNode + 2;
        Interfaces.(l).coord_X  = DummyStruct.nodes(n, 1);
        Interfaces.(l).coord_Y  = DummyStruct.nodes(n, 2);
    end
end

%% --- 3.2 FORZANTE E CONDIZIONI INIZIALI ---
n_dof_per_node = DummyStruct.MeshObj.nDOFPerNode;
impulse_dir    = impulse_sign * [cosd(impulse_angle_deg); sind(impulse_angle_deg)];

dir_vector = zeros(n_dofs_fom, 1);
dir_vector(1:n_dof_per_node:n_dofs_fom) = impulse_dir(1);
dir_vector(2:n_dof_per_node:n_dofs_fom) = impulse_dir(2);

F_spatial_fom = Mc * dir_vector;
F_fom_handle  = @(t) F_spatial_fom * impulse_amp * sin(pi*t/t_shock) * (t <= t_shock);

q0  = zeros(n_dofs_fom, 1);
qd0 = zeros(n_dofs_fom, 1);

v_max = impulse_amp * 2 * t_shock / pi;
m_eff = dir_vector' * Mc * dir_vector;
Eref  = 0.5 * m_eff * v_max^2;          % scalare fisico, indipendente dalla base

t_common = 0 : 10*dt : tmax;

fprintf('Impulso: %.1e g @ %.1f deg (sign %+d)\n', impulse_g, impulse_angle_deg, impulse_sign);
fprintf('Eref = %.4e J   (v_max = %.4f m/s)\n', Eref, v_max);
fprintf('Interfaccia: %d DOF | rapporto b/(2r) atteso a phi=%d: %.1f%%\n\n', ...
        n_bnd, max_phi, 100*n_bnd/(2*(max_phi+n_bnd)));

%% --- 4. CONFIGURAZIONI ---
cfgs = struct( ...
  'name',   {'A_raw_std','B_raw_proj','C_scaled_proj','Cp_scaled_std','D_mc_gs'}, ...
  'method', {'rubin',    'rubin',     'rubin',        'rubin',         'mc'    }, ...
  'scaling',{false,      false,       true,           true,            false   }, ...
  'branch', {'std',      'proj',      'proj',         'std',           'proj'  }, ...
  'enable', {true,       true,        true,           true,            true    });

summary = table();

%% --- 5. RUN ---
for Q = array_QFactor
    DummyStruct.compute_rayleigh_damping(Q, Q);

    for k_mult = array_k_mult
        k_contact = k_base * k_mult;

        epsE = 1e-2 * bench_RelTol;
        abstol_fixed = epsE * min(sqrt(2*Eref/k_contact), sqrt(2*Eref/m_eff));

        for phi = array_linModes
            fprintf('\n===== Phi %d | Q %d | Kmult %g =====\n', phi, Q, k_mult);

            % ---------- build Rubin (una volta sola) ----------
            tic_off = tic;
            rom_r = RomRubinTest(DummyStruct, phi, contact_dofs, false);
            rom_r.build();
            off_rubin = toc(tic_off);
            fprintf(' Offline Rubin: %.2f s\n', off_rubin);

            s_if = rom_r.diag_info.sigma_iface;
            fprintf(' sigma interfaccia (norm.): '); fprintf('%.2e ', s_if(1:min(end,10)));
            fprintf('...\n rango effettivo @1e-6: %d / %d | cond(T1_bb) = %.2e\n', ...
                    sum(s_if > 1e-6), numel(s_if), rom_r.diag_info.cond_T1bb);

            % ---------- build MC (se richiesto) ----------
            rom_mc = []; Pc_mc = []; off_mc = NaN;
            if any([cfgs(strcmp({cfgs.method},'mc')).enable])
                tic_off = tic;
                rom_mc = RomMC(DummyStruct, phi, contact_dofs, k_contact, 1);
                rom_mc.build();
                Pc_mc = zeros(n_dofs_fom, size(rom_mc.P, 2));
                for i = 1:size(Pc_mc,2)
                    Pc_mc(:,i) = DummyStruct.AssemblyObj.constrain_vector(rom_mc.P(:,i));
                end
                off_mc = toc(tic_off);
                Mr_chk = full(Pc_mc' * Mc * Pc_mc);
                fprintf(' Offline MC: %.2f s | r = %d | cond(Mr) = %.3e\n', ...
                        off_mc, size(Pc_mc,2), cond(Mr_chk));
            end

            % ---------- loop configurazioni ----------
            for c = 1:numel(cfgs)
                cfg = cfgs(c);
                if ~cfg.enable, continue; end
                fprintf('\n>>> %s  (metodo %s, scaling %d, ramo %s)\n', ...
                        cfg.name, cfg.method, cfg.scaling, cfg.branch);

                % --- base, matrici ridotte, proiezione ---
                switch cfg.method
                    case 'rubin'
                        rom_r.set_scaling(cfg.scaling);
                        [Mr, Kr, Cr] = rom_r.get_reduced_matrices();
                        Pc_r     = rom_r.Pc;
                        d_if     = rom_r.scaleD(1:n_bnd);
                        offline_time = off_rubin;
                        cond_Mr  = rom_r.diag_info.cond_Mr;
                        cond_Kr  = rom_r.diag_info.cond_Kr;
                    case 'mc'
                        Mr = rom_mc.P' * DummyStruct.M * rom_mc.P;
                        Kr = rom_mc.P' * DummyStruct.K * rom_mc.P;
                        Cr = rom_mc.P' * DummyStruct.C * rom_mc.P;
                        Mr = (Mr+Mr')/2; Kr = (Kr+Kr')/2; Cr = (Cr+Cr')/2;
                        Pc_r     = Pc_mc;
                        d_if     = ones(n_bnd,1);
                        offline_time = off_mc;
                        cond_Mr  = cond(full(Mr));
                        cond_Kr  = cond(full(Kr));
                end

                % --- CI proiettate: forma valida anche quando Mr ~= I ---
                q0_r  = Mr \ (Pc_r' * (Mc * q0));
                qd0_r = Mr \ (Pc_r' * (Mc * qd0));
                F_r_handle = @(t) Pc_r' * F_fom_handle(t);

                % --- ramo del solutore e semantica di ContactTargetDOF ---
                %   'std'  -> COORDINATE RIDOTTE (1:n_bnd). Con scaling attivo
                %             q_fis = y/d, quindi gap e k vanno riscalati:
                %                gap'_i = d_i*gap_i     k'_i = k/d_i^2
                %             I segni si conservano perche' d_i > 0.
                %   'proj' -> DOF FISICI constrained (contact_dofs). La
                %             proiezione assorbe i fattori 1/d automaticamente.
                switch lower(cfg.branch)
                    case 'std'
                        assert(strcmp(cfg.method,'rubin'), ...
                            'Il ramo standard richiede un''interfaccia diagonale (solo Rubin).');
                        model_type  = 'Rubin';
                        target_dofs = (1:n_bnd)';
                        gaps_run    = d_if .* gaps_array;
                        k_run       = k_contact ./ d_if.^2;
                        extra_args  = {};
                    case 'proj'
                        model_type  = 'MC';
                        target_dofs = contact_dofs;
                        gaps_run    = gaps_array;
                        k_run       = k_contact * ones(n_bnd,1);
                        extra_args  = {'ProjectionMatrix', Pc_r};
                end

                if use_fixed_abstol
                    extra_args = [extra_args, {'AbsTolOverride', abstol_fixed}]; %#ok<AGROW>
                end

                % --- integrazione ---
                tic;
                solverROM = TransientSolverOde_V2_test(Mr, Kr, Cr);
                [t_rom, q_rom, ode_stats] = solverROM.solve(tmax, dt, q0_r, qd0_r, F_r_handle, ...
                    'ContactTargetDOF', target_dofs, ...
                    'ContactGap',       gaps_run, ...
                    'ContactPenalty',   k_run, ...
                    'ModelType',        model_type, ...
                    'Eref',             Eref, ...
                    'RelTol',           bench_RelTol, ...
                    'OutputTimes',      t_common, ...
                    extra_args{:});
                cpu_time = toc;

                % --- ricostruzione fisica (valida per qualunque base) ---
                q_c        = Pc_r * q_rom;
                y_rom_full = DummyStruct.AssemblyObj.unconstrain_vector(q_c);

                y_contact = struct();
                for loc = labels
                    l = loc{1};
                    if ~isempty(Interfaces.(l).nodes)
                        y_contact.(l).X = y_rom_full(Interfaces.(l).global_X, :);
                        y_contact.(l).Y = y_rom_full(Interfaces.(l).global_Y, :);
                    end
                end

                % --- amplificazione delle coordinate (cancellazione) ---
                amp = max(vecnorm(q_rom)) / ...
                      max(max(sqrt(abs(sum(q_c .* (Mc*q_c), 1)))), realmin);

                nfail_pct = 100 * ode_stats.nfailed / ode_stats.nsteps;
                lu_per_jac = ode_stats.ndecomps / max(ode_stats.npds, 1);

                fprintf('  CPU %.2f s | cond(Mr) %.3e | falliti %.2f%% | LU/Jac %.1f | amp %.2e\n', ...
                        cpu_time, cond_Mr, nfail_pct, lu_per_jac, amp);

                diag_info = struct('cfg_name', cfg.name, 'method', cfg.method, ...
                    'branch', cfg.branch, 'scaling', cfg.scaling, ...
                    'r', size(Mr,1), 'cond_Mr', cond_Mr, 'cond_Kr', cond_Kr, ...
                    'amp_coords', amp, 'cpu_time', cpu_time, ...
                    'offline_time', offline_time, 'stats', ode_stats);

                summary = [summary; table( ...
                    string(cfg.name), phi, size(Mr,1), cond_Mr, cond_Kr, cpu_time, ...
                    ode_stats.nsteps, nfail_pct, ode_stats.ndecomps, ...
                    ode_stats.npds, lu_per_jac, amp, ...
                    'VariableNames', {'cfg','phi','r','condMr','condKr','cpu', ...
                                      'nsteps','failedPct','nLU','nJac','LUperJac','amp'})]; %#ok<AGROW>

                file_name = sprintf('ROM_%s_Phi%03d_Q%04d_K%g.mat', cfg.name, phi, Q, k_mult);
                save(fullfile(save_dir, file_name), 't_rom', 'y_contact', ...
                     'Interfaces', 'cpu_time', 'offline_time', 'diag_info');
            end
        end
    end
end

%% --- 6. RIEPILOGO ---
fprintf('\n================ RIEPILOGO ================\n');
disp(summary);
writetable(summary, fullfile(save_dir, 'summary.csv'));
save(fullfile(save_dir, 'summary.mat'), 'summary');

fprintf('\nConfronti:\n');
fprintf('  B vs A  : costo del ramo proiettato, O(b*r) vs O(b)\n');
fprintf('  C vs B  : effetto dello scaling sul condizionamento\n');
fprintf('  C'' vs C : costo del ramo, a AbsTol identica\n');
fprintf('  C'' vs D : Rubin scalata vs Milman-Chu, entrambi ben condizionati\n');
fprintf('Il GRE deve essere identico in A, B, C, C''.\n');
