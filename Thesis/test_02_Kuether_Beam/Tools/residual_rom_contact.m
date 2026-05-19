function [r, drdqdd, drdqd, drdq, c0] = residual_rom_contact(q, qd, qdd, t, ReducedAssembly, spring, P, F_ext)
    % RESIDUAL_ROM_CONTACT Evaluates the residual and Jacobians for the ROM with impact
    
    M = ReducedAssembly.DATA.M;
    K = ReducedAssembly.DATA.K;
    C = ReducedAssembly.DATA.C;
    
    % External and non-linear forces
    F_external = F_ext(t);
    f_nl_red = spring.evaluate_reduced(q, P);
    K_nl_red = spring.jacobian_reduced(q, P);
    
    % Physical force decomposition (standard YAFEC approach)
    F_inertial = M * qdd;
    F_damping  = C * qd;
    F_elastic  = K * q + f_nl_red; % Linear + Non-linear (Spring)
    
    % 1. Calculate Residual vector r
    r = F_inertial + F_damping + F_elastic - F_external;
    
    % 2. Calculate directional derivatives (Jacobians)
    drdqdd = M;
    drdqd  = C;
    drdq   = K + K_nl_red;
    
    % 3. Calculate normalization factor c0
    % Add 'eps' (MATLAB's machine epsilon) to prevent division by zero 
    % when the beam is stationary.
    c0 = norm(F_inertial) + norm(F_damping) + norm(F_elastic) + norm(F_external) + eps;
end