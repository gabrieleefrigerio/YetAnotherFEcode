%% Model Parameters
m = 1;       % Massa [kg]
k = 1000;    % Rigidezza [N/m]
xsi = 0;
x0 = 0.1;    % Spostamento iniziale [m]
v0 = 0;      % Velocità iniziale [m/s]
%% Reference solution:
c = 2*xsi*sqrt(k*m);    
% Definizione del sistema dX/dt = A*X
A = [0, 1; 
    -k/m, -c/m];

% ode45
sys_ode = @(t, X) A * X;

t_span = [0 0.5]; 
X0 = [x0; v0];
[t_ode, X_ode] = ode45(sys_ode, t_span, X0);
x_ref = X_ode(:, 1); % Estraiamo lo spostamento

%% YaFEc solution:
% ode45
% Creazione delle matrici di sistema per yaFEc
M_mat = m;
K_mat = k;
C_mat = c;
F_ext = 0; % Nessuna forzante, risposta libera

%% YaFEc solution: Implicit Newmark
% 1. Condizioni iniziali complete
q0 = x0;
qd0 = v0;
% Calcolo dell'accelerazione iniziale dall'equilibrio a t=0: M*a0 + C*v0 + K*x0 = F_ext
qdd0 = (F_ext - C_mat*qd0 - K_mat*q0) / M_mat;

% 2. Definizione del time step (h)
% Calcoliamo il periodo naturale del sistema per avere una buona risoluzione
omega_n = sqrt(k/m);
T_n = 2*pi/omega_n;
h = T_n / 50; % 50 passi di integrazione per ogni periodo di oscillazione

% 3. Definizione della funzione Residuo
% Invece di usare un oggetto Assembly, scriviamo direttamente l'equazione di equilibrio
residual_lin = @(q, qd, qdd, t) my_linear_residual(q, qd, qdd, t, M_mat, C_mat, K_mat, F_ext);
% 4. Inizializzazione dell'integratore YaFEc
% Nota: impostiamo 'alpha' a 0. In molti schemi derivati da Newmark (come HHT), 
% alpha introduce smorzamento numerico. Poiché stai testando un sistema c=0 
% e vuoi confrontarlo con ode45, alpha=0 garantisce nessuna dissipazione artificiale.
TI = ImplicitNewmark('timestep', h, 'alpha', 0, 'linear', true);

% 5. Integrazione nel tempo
tmax = t_span(2); % Allineiamo il tempo finale a quello di ode45 (0.5 s)
TI.Integrate(q0, qd0, qdd0, tmax, residual_lin);

% 6. Estrazione dei risultati
t_yafec = TI.Solution.time;
x_yafec = TI.Solution.q; % Vettore degli spostamenti calcolati da YaFEc

%% Confronto visivo (Plotting)
figure('Name', 'Confronto ODE45 vs YaFEc (Newmark)', 'Color', 'w');

% Plot ODE45 (linea solida nera)
plot(t_ode, x_ref, 'k-', 'LineWidth', 1.5, 'DisplayName', 'ode45');
hold on;

% Plot YaFEc Newmark (linea tratteggiata rossa)
plot(t_yafec, x_yafec, 'r--', 'LineWidth', 2, 'DisplayName', 'YaFEc (Implicit Newmark)');

% Formattazione grafico
grid on;
xlabel('Tempo [s]', 'Interpreter', 'latex');
ylabel('Spostamento $x(t)$ [m]', 'Interpreter', 'latex');
title('Risposta Libera Sistema 1-DOF: ode45 vs Implicit Newmark', 'Interpreter', 'latex');
legend('show', 'Location', 'best', 'Interpreter', 'latex');
axis tight;

% figure;
% plot(t_ode, x_ref, 'k-', 'LineWidth', 1.5); hold on;
% plot(t_yafec, x_yafec, 'r--', 'LineWidth', 1.5);
% xlabel('Tempo [s]');
% ylabel('Spostamento [m]');
% legend('ode45 (Esatta)', 'yaFEc (Implicita)');
% title('Validazione Integratore Temporale 1-DOF');
% grid on;

function [r, drdqdd, drdqd, drdq] = my_linear_residual(q, qd, qdd, t, M, C, K, F)
    % 1. Calcolo del residuo (Forze d'inerzia + smorzamento + elastiche - Forze esterne)
    r = M*qdd + C*qd + K*q - F;
    
    % 2. Calcolo dei Jacobiani (Derivate del residuo)
    drdqdd = M;  % Derivata rispetto a q_ddot
    drdqd  = C;  % Derivata rispetto a q_dot
    drdq   = K;  % Derivata rispetto a q
end