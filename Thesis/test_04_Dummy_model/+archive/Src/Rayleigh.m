% ==========================================================
% STIMA DEI COEFFICIENTI DI RAYLEIGH DAI FATTORI DI QUALITA'
% ==========================================================
Q1 = 1000; 
Q2 = 1000; 

zeta_1 = 1 / (2 * Q1);
zeta_2 = 1 / (2 * Q2);


[~, D_fom] = eigs(Kc, Mc, 2, 'smallestabs');
omega = sort(sqrt(diag(D_fom)));
w1 = omega(1);
w2 = omega(2);

alpha_rayleigh = (2 * w1 * w2 * (zeta_1 * w2 - zeta_2 * w1)) / (w2^2 - w1^2);
beta_rayleigh = (2 * (zeta_2 * w2 - zeta_1 * w1)) / (w2^2 - w1^2);

fprintf('\n--- Rayleigh Coefficients ---\n');
fprintf('alpha = %e\n', alpha_rayleigh);
fprintf('beta  = %e\n', beta_rayleigh);
