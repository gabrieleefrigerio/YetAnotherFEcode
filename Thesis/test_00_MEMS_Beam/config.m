function cfg = config()
    % --- Parametri di Base ---
    cfg.material.E = 169e9;
    cfg.material.rho = 2330;
    cfg.material.nu = 0.22;

    cfg.geom.W = 15e-6;
    cfg.geom.thickness = 20e-6;
    
    cfg.mesh.nx = 20;
    cfg.mesh.ny = 4;
    cfg.mesh.elementType = 'QUAD8'; % Most accurate element type for bending

    cfg.run_type = 'single_run'; % 'single_run' o 'sweep_length'
    
    % single_run', usa questo:
    cfg.geom.L = 2e-3; 
    
    % sweep, usa questi limiti:
    cfg.sweep.L_min = 1e-3;
    cfg.sweep.L_max = 3e-3;
    cfg.sweep.steps = 10;
end