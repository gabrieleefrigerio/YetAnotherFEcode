function [r, drdqdd, drdqd, drdq, c0] = residual_fom(q, qd, qdd, t_curr, M, K, C, target_dofs, gap_wall, k_penalty, F_handle)
    F_structural = K * q;
    F_damping = C * qd;
    F_inertia = M * qdd;
    
    n_dofs = length(q);
    F_penalty = zeros(n_dofs, 1);
    
    % Mappatura diretta: q contiene gli spostamenti fisici
    current_disps = q(target_dofs); 
    is_penetrating = current_disps > gap_wall;
    
    if any(is_penetrating)
        active_dofs = target_dofs(is_penetrating);
        penetrations = current_disps(is_penetrating) - gap_wall;
        
        F_penalty(active_dofs) = k_penalty * penetrations;
        K_penalty_mat = sparse(active_dofs, active_dofs, k_penalty, n_dofs, n_dofs);
    else
        K_penalty_mat = sparse(n_dofs, n_dofs);
    end
    
    F_ext_current = F_handle(t_curr); 
    r = F_inertia + F_damping + F_structural + F_penalty - F_ext_current;
    
    drdqdd = M;
    drdqd  = C;
    drdq   = K + K_penalty_mat; 
    
    F_elastic_tot = F_structural + F_penalty;
    c0 = norm(F_inertia) + norm(F_damping) + norm(F_elastic_tot) + norm(F_ext_current);
    if c0 == 0; c0 = 1; end
end