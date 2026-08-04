% =========================================================================
% SCRIPT FOR MODEL VALIDATION (DUMMYSTRUCTURE)
% =========================================================================

DummyStruct = AbaqusStructure();
DummyStruct.filename = 'DummyStructureAbaqus.inp';
DummyStruct.elementType = 'TRI3';
DummyStruct.build();
DummyStruct.compute_eigenmodes(100);
n_modes = length(DummyStruct.frequencies);

Kc = DummyStruct.AssemblyObj.constrain_matrix(DummyStruct.K);
Mc = DummyStruct.AssemblyObj.constrain_matrix(DummyStruct.M);

% CONDITION NUMBER
cond_K = condest(Kc);
cond_M = condest(Mc);

fprintf('Condition Number for Kc (Stiffness): %e\n', cond_K);
fprintf('Condition Number for Mc (Mass):      %e\n', cond_M);

% RESIDUALS (EIG PROBLEM)
fprintf('\n--- 2. Residui dei Problemi agli Autovalori ---\n');
fprintf('            || Kc*phi - lambda*Mc*phi ||            \n\n');

for ii = 1:n_modes
    omega = DummyStruct.frequencies(ii) * 2 * pi;
    lambda = omega^2;
    
    % get mode
    phi_full = DummyStruct.mode_shapes(:, ii);
    phi_c = DummyStruct.AssemblyObj.constrain_vector(phi_full);
    
    % residual
    residual_vec = Kc * phi_c - lambda * (Mc * phi_c);
    res_norm = norm(residual_vec);
    
    % Term for normalization
    stiffness_term_norm = norm(Kc * phi_c);
    
    if stiffness_term_norm > 0
        rel_res_norm = res_norm / stiffness_term_norm;
    else
        rel_res_norm = NaN;
    end
    
    fprintf('Mode %3d (%10.3f Hz) | Residiual: %e\n', ...
            ii, DummyStruct.frequencies(ii), rel_res_norm);
           
end

