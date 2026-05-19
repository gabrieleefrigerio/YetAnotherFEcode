function [r, drdqdd, drdqd, drdq, c0] = residual_fom_contact(u_c, ud_c, udd_c, t, F_ext_c, M_c, K_c, Assembly, spring)
    % Smorzamento nullo
    C_c = sparse(size(M_c, 1), size(M_c, 2));
    
    % --- 1. Forze interne lineari (TRAVE PURAMENTE LINEARE) ---
    % Sostituiamo il calcolo non-lineare di YaFEc con la semplice legge di Hooke
    F_elastic_beam_c = K_c * u_c;
    
    % --- 2. Calcolo forza e rigidezza non lineare (MOLLA DI CONTATTO) ---
    u_full = Assembly.unconstrain_vector(u_c);
    u_nl = u_full(spring.dof_idx);
    
    f_full = zeros(size(u_full, 1), 1);
    K_spring_full = sparse(size(u_full, 1), size(u_full, 1));
    
    g = u_nl - spring.a;
    if g > 0
        f_val = spring.k * g;
        k_val = spring.k;
    else
        f_val = 0;
        k_val = 0;
    end
    
    f_full(spring.dof_idx) = f_val;
    K_spring_full(spring.dof_idx, spring.dof_idx) = k_val;
    
    f_nl_spring_c = Assembly.constrain_vector(f_full);
    K_spring_c = Assembly.constrain_matrix(K_spring_full);
    
    % --- 3. Forza esterna ---
    F_external = F_ext_c(t);
    
    % --- 4. Calcolo Vettore Residuo ---
    F_inertial = M_c * udd_c;
    F_damping  = C_c * ud_c;
    F_elastic  = F_elastic_beam_c + f_nl_spring_c; 
    
    r = F_inertial + F_damping + F_elastic - F_external;
    
    % --- 5. Jacobiani per Newton-Raphson ---
    drdqdd = M_c;
    drdqd  = C_c;
    drdq   = K_c + K_spring_c; % La rigidezza è K costante + la molla
    
    c0 = norm(F_inertial) + norm(F_damping) + norm(F_elastic) + norm(F_external) + eps;
end