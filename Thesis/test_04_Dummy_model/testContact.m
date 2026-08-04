%% =====================================================================
%  TEST: si puo' portare G' dentro la parentesi di Macaulay?
%  Da eseguire DOPO la costruzione di Pc_mc (o Pc_mt) nel main.
%
%  Ipotesi di lavoro:
%     Q_c(q) = k_p * G' * <-(g0 + G*q)>          (forma corretta)
%     Q_c(q) =? k_p *      <-(G'*g0 + (G'*G)*q)> (forma con G' dentro)
%
%  Convenzione uniforme: si riscrive il gap in modo che g < 0 <=> penetrazione,
%  assorbendo il segno del muro (gaps_array puo' essere positivo o negativo).
% =====================================================================

%% --- 0. Costruzione di G e g0 nella convenzione uniforme ---
% Il test usa la base MC; per MT sostituire Pc_mc -> Pc_mt.
Pc_test = Pc_mc;

G_raw = Pc_test(contact_dofs, :);      % n_c x m
s     = sign(gaps_array);              % direzione del muro per ciascun nodo
G     = -s .* G_raw;                   % n_c x m, convenzione g<0 = penetrazione
g0    = abs(gaps_array);               % n_c x 1, distanze iniziali (tutte > 0)

nc = size(G, 1);
m  = size(G, 2);

kp = k_contact;

fprintf('\n==============================================\n');
fprintf('  TEST SCAMBIO G'' <-> MACAULAY\n');
fprintf('==============================================\n');
fprintf('n_c = %d nodi candidati | m = %d coordinate ridotte\n', nc, m);

%% --- 1. Verifica della convenzione sulla configurazione iniziale ---
% A q = 0 la struttura e' a riposo: tutti i gap devono risultare aperti (g > 0).
q_rest = zeros(m, 1);
g_rest = g0 + G * q_rest;

fprintf('\n--- 1. Configurazione di riposo (q = 0) ---\n');
fprintf('min(g) = %.4e  [deve essere > 0: nessun contatto a riposo]\n', min(g_rest));
fprintf('nodi con g < 0 : %d  [atteso 0]\n', sum(g_rest < 0));

%% --- 2. Confronto delle due forme sullo stato di riposo ---
% Forma corretta: taglio prima (per nodo), somma dopo (via G').
Q_ok  = kp * (G.' * max(0, -(g0 + G*q_rest)));

% Forma con G' dentro: somma prima, taglio dopo.
Q_bad = kp * max(0, -(G.'*g0 + (G.'*G)*q_rest));

fprintf('\n--- 2. Forza di contatto a riposo ---\n');
fprintf('|Q| forma corretta      = %.6e   [deve essere 0]\n', norm(Q_ok));
fprintf('|Q| forma con G'' dentro = %.6e   [se > 0: forza spuria]\n', norm(Q_bad));

if norm(Q_ok) > 0
    warning('La forma corretta non e'' nulla a riposo: controllare la convenzione dei segni.');
end

%% --- 3. Confronto su stati con contatto attivo ---
% Si eccita ciascun modo con ampiezza crescente finche' qualche nodo penetra,
% cosi' da confrontare le due forme in regime di contatto reale.
fprintf('\n--- 3. Stati con contatto attivo ---\n');
fprintf('%6s %8s %14s %14s %12s\n', 'modo', 'n_att', '|Q| corretta', '|Q| G'' dentro', 'err. rel.');

rng(0);
n_trial   = 8;
err_rel   = nan(n_trial, 1);
n_act_vec = nan(n_trial, 1);

for k = 1:n_trial
    % direzione casuale nello spazio ridotto
    dq = randn(m, 1);
    dq = dq / norm(dq);

    % scalatura tale da portare almeno un nodo in penetrazione
    Gd = G * dq;
    neg = Gd < -eps;
    if ~any(neg),  continue,  end
    alpha = 1.2 * min(g0(neg) ./ (-Gd(neg)));   % 20% oltre il primo contatto
    q_try = alpha * dq;

    g_try = g0 + G * q_try;
    act   = g_try < 0;

    Q_ok_k  = kp * (G.' * max(0, -g_try));
    Q_bad_k = kp * max(0, -(G.'*g0 + (G.'*G)*q_try));

    n_act_vec(k) = sum(act);
    err_rel(k)   = norm(Q_bad_k - Q_ok_k) / max(norm(Q_ok_k), eps);

    fprintf('%6d %8d %14.4e %14.4e %12.3e\n', ...
        k, n_act_vec(k), norm(Q_ok_k), norm(Q_bad_k), err_rel(k));
end

fprintf('\nerrore relativo mediano: %.3e\n', median(err_rel(~isnan(err_rel))));

%% --- 4. Condizione algebrica di validita' dello scambio ---
% Lo scambio G'<.> = <G'.> e' esatto se e solo se ogni colonna di G ha al piu'
% un elemento non nullo (e non negativo): in tal caso la somma su i contiene
% un solo termine e taglio e somma commutano.
tol     = 1e-10 * max(abs(G(:)));
nnz_col = sum(abs(G) > tol, 1);

fprintf('\n--- 4. Struttura di G ---\n');
fprintf('nonzeri per colonna: min %d | mediana %d | max %d\n', ...
    min(nnz_col), round(median(nnz_col)), max(nnz_col));
fprintf('densita'' di G: %.1f%%\n', 100 * nnz(abs(G) > tol) / numel(G));

if max(nnz_col) <= 1
    fprintf('=> G e'' una matrice di selezione: lo scambio e'' LECITO.\n');
    fprintf('   (ma in tal caso G*q e'' un''indicizzazione, non una proiezione)\n');
else
    fprintf('=> G ha colonne con piu'' nonzeri: lo scambio NON e'' lecito.\n');
end

%% --- 5. Modalizzazione di g0 ---
% Esiste q0 tale che G*q0 = -g0 ? In tal caso g = G*(q - q0) e g0 sparisce.
q0_mod = -G \ g0;
res_g0 = norm(g0 + G*q0_mod) / norm(g0);

fprintf('\n--- 5. Modalizzazione di g0 ---\n');
fprintf('residuo relativo ||g0 + G*q0|| / ||g0|| = %.3e\n', res_g0);
if res_g0 < 1e-10
    fprintf('=> g0 e'' esattamente rappresentabile: g = G*(q - q0).\n');
else
    fprintf('=> g0 NON e'' nel range di G: la modalizzazione sarebbe approssimata.\n');
end