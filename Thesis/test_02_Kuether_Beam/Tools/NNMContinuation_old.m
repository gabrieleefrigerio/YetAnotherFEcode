classdef NNMContinuation_old < handle
    properties
        ReducedAssembly 
        spring          
        P               
        SystemNLvib     
        
        X_out           
        energies        
        frequencies
        rms_displacements % <--- NUOVA PROPRIETA' AGGIUNTA
        modal_amplitudes % <--- NUOVA PROPRIETA'
    end
    
    methods
        function obj = NNMContinuation(reduced_assembly, spring_obj, P_matrix)
            obj.ReducedAssembly = reduced_assembly;
            obj.spring = spring_obj;
            obj.P = P_matrix;
            
            % Iniettiamo il wrapper intelligente per la direzione della molla
            obj.ReducedAssembly.DATA.fnl_CUSTOM = @(q) obj.nl_force_wrapper(q);
            
            n_dofs = size(P_matrix, 2);
            F_ext_null = zeros(n_dofs, 1);
            obj.SystemNLvib = FE_system(obj.ReducedAssembly, F_ext_null, 'custom');
        end
        
        function [Kt, f_nl] = nl_force_wrapper(obj, q)
    u_full = obj.P * q;
    u_nl = u_full(obj.spring.dof_idx);
    
    f_full = zeros(size(u_full, 1), 1);
    K_full = sparse(size(u_full, 1), size(u_full, 1));
    
    % --- REGOLARIZZAZIONE DEL CONTATTO ---
    % g è la compenetrazione
    g = u_nl - obj.spring.a;
    
    % epsilon controlla la "dolcezza" della transizione. 
    % Deve essere qualche ordine di grandezza più piccolo del gap 'a'
eps_smooth = obj.spring.a * 0.10; 
    if eps_smooth == 0, eps_smooth = 1e-6; end % safety fallback
    
    % Forza regolarizzata: se g<<0 è quasi 0; se g>>0 è k*g.
    % È una curva iperbolica che smussa l'angolo dell'impatto.
    f_val = (obj.spring.k / 2) * (g + sqrt(g^2 + eps_smooth^2));
    
    % Rigidezza tangente: derivata analitica esatta della formula sopra.
    % Passa da 0 a k in modo continuo.
    k_val = (obj.spring.k / 2) * (1 + g / sqrt(g^2 + eps_smooth^2));
    
    % Assegniamo i valori (niente più blocco IF, la formula vale sempre)
    f_full(obj.spring.dof_idx) = f_val;
    K_full(obj.spring.dof_idx, obj.spring.dof_idx) = k_val;
    
    f_nl = obj.P' * f_full;
    Kt = obj.P' * K_full * obj.P;
end
        
function solve(obj, mode_idx, log10a_start, log10a_end)
    M_r = obj.ReducedAssembly.DATA.M;
    K_r = obj.ReducedAssembly.DATA.K;
    [Phi_red, Om2_red] = eig(K_r, M_r);
    [om_lin_red, sort_idx] = sort(sqrt(diag(Om2_red)));
    Phi_red = Phi_red(:, sort_idx);
    
    om_start = om_lin_red(mode_idx);
    phi_start = Phi_red(:, mode_idx);
    
    % --- NUOVO: Orientiamo il modo iniziale verso la molla fisica (+a) ---
    phi_phys = obj.P * phi_start;
    if phi_phys(obj.spring.dof_idx) < 0
        phi_start = -phi_start; % Ribaltiamo il modo se punta dalla parte sbagliata
    end
    
            
            % Normalizziamo in modo standard
            [~, inorm] = max(abs(phi_start)); 
            phi_start = phi_start / phi_start(inorm); 
            
            m = size(M_r, 1);
            nnorm = setdiff(1:2*m, [inorm, inorm+m]);
            ys0 = [phi_start; zeros(m, 1)]; 
            a0 = ys0(inorm);
            x0 = [ys0(nnorm)/a0; om_start; 0];
            
            % --- Parametri dello Shooting Method ---
            Ntd = 500; % Alta risoluzione temporale per l'impatto severo
            Np = 1;                    
            qscl = obj.spring.a;       
            fscl = mean(diag(K_r)) * qscl; 
            
            % Passo della continuazione (piccolo per seguire meglio la curva)
            ds = 0.00001;                 
            
            dscale = [ones(length(x0)-2, 1); om_start; 1e0; 1e0];
            Sopt = struct('Dscale', dscale, 'dynamicDscale', 1, 'stepmax', 10000, 'dsmin', ds/100000, 'dsmax', 100*ds);
            
            fprintf('Avvio Calcolo NNM %d (Shooting Method)...\n', mode_idx);
            
            obj.X_out = solve_and_continue(x0, ...
                @(X) shooting_residual(X, obj.SystemNLvib, Ntd, Np, 'NMA', qscl, fscl, inorm), ...
                log10a_start, log10a_end, ds, Sopt);
                
            fprintf('Calcolo completato!\n');
            obj.compute_energies(inorm);
        end
        
        function compute_energies(obj, inorm)
            m = size(obj.P, 2);
            om_sh = obj.X_out(end-2, :);
            log10a_sh = obj.X_out(end, :);
            a_sh = 10.^log10a_sh;
            
            nnorm = setdiff(1:2*m, [inorm, inorm+m]);
            
            obj.frequencies = om_sh / (2*pi);
            obj.energies = zeros(1, size(obj.X_out, 2));
            obj.modal_amplitudes = zeros(m, size(obj.X_out, 2));
            % Inizializziamo l'array dell'RMS
            obj.rms_displacements = zeros(1, size(obj.X_out, 2)); 
            
            M_r = obj.ReducedAssembly.DATA.M;
            K_r = obj.ReducedAssembly.DATA.K;
            
            dir_mult = sign(obj.P(obj.spring.dof_idx, 1));
            if dir_mult == 0, dir_mult = 1; end
            
            for i = 1:size(obj.X_out, 2)
                ai = a_sh(i);
                ys_nnorm_i = obj.X_out(1:end-3, i) * ai;
                
                ys = zeros(2*m, 1);
                ys(inorm) = ai;
                ys(nnorm) = ys_nnorm_i;
                
                qi = ys(1:m);
                vi = ys(m+1:end);
                
                % Salviamo il valore assoluto della partecipazione di ogni modo
                obj.modal_amplitudes(:, i) = abs(qi); % <--- SALVA QUI
                % Ricostruzione dello spostamento fisico dell'intera trave
                x_full = obj.P * qi;
                
                % --- CALCOLO RMS SPAZIALE ---
                % Calcola la radice del valore quadratico medio (RMS) di tutti i GDL
                obj.rms_displacements(i) = sqrt(mean(x_full.^2)); 
                
                % (Nota: Se preferisci plottare solo l'RMS temporale del nodo di contatto, 
                % de-commenta la riga seguente e commenta quella sopra)
                % obj.rms_displacements(i) = abs(x_full(obj.spring.dof_idx)) / sqrt(2);
                
                % Calcolo energia (Mantenuto per archivio)
                E_lin = 0.5 * vi' * M_r * vi + 0.5 * qi' * K_r * qi;
                
                u_nl = x_full(obj.spring.dof_idx) * dir_mult; 
                
                E_nl = 0;
                if u_nl > obj.spring.a
                    E_nl = 0.5 * obj.spring.k * (u_nl - obj.spring.a)^2;
                end
                
                obj.energies(i) = E_lin + E_nl;
            end
        end
        
        function plot_backbone(obj, style_str, display_name)
            if nargin < 2, style_str = 'b-'; end
            if nargin < 3, display_name = 'NNM Shooting'; end
            
            % Usiamo ancora semilogx perché gli spostamenti coprono
            % diversi ordini di grandezza (da micrometri a millimetri)
            semilogx(obj.energies, obj.frequencies, style_str, 'LineWidth', 2, 'DisplayName', display_name);
            grid on;
            xlabel('RMS Spostamento [m]');
            ylabel('Frequenza [Hz]');
        end
    end
end