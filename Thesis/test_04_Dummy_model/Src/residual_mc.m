function [r, drdqdd, drdqd, drdq, c0] = residual_mc(q, qd, qdd, t_curr, M, K, C, Pc, target_dofs_phys, gap_wall, k_penalty, F_handle)
    F_structural = K * q;
    F_damping = C * qd;
    F_inertia = M * qdd;
    
    n_dofs_rom = length(q);
    n_phys = size(Pc, 1);
    
    F_penalty_phys = zeros(n_phys, 1);
    
    % Proiezione nello spazio fisico per "vedere" il muro
    q_phys = Pc * q;
    current_disps = q_phys(target_dofs_phys); 
    is_penetrating = current_disps > gap_wall;
    
    if any(is_penetrating)
        active_dofs = target_dofs_phys(is_penetrating);
        penetrations = current_disps(is_penetrating) - gap_wall;
        
        F_penalty_phys(active_dofs) = k_penalty * penetrations;
        K_penalty_phys = sparse(active_dofs, active_dofs, k_penalty, n_phys, n_phys);
        
        % Proiezione inversa nel ROM
        F_penalty_rom = Pc' * F_penalty_phys;
        K_penalty_rom = Pc' * K_penalty_phys * Pc;
    else
        F_penalty_rom = zeros(n_dofs_rom, 1);
        K_penalty_rom = sparse(n_dofs_rom, n_dofs_rom);
    end
    
    F_ext_current = F_handle(t_curr); 
    r = F_inertia + F_damping + F_structural + F_penalty_rom - F_ext_current;
    
    drdqdd = M;
    drdqd  = C;
    drdq   = K + K_penalty_rom; 
    
    F_elastic_tot = F_structural + F_penalty_rom;
    c0 = norm(F_inertia) + norm(F_damping) + norm(F_elastic_tot) + norm(F_ext_current);
    if c0 == 0; c0 = 1; end
end