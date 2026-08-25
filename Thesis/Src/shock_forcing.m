function shock = shock_forcing(Struct, cfg)
%SHOCK_FORCING Half-sine acceleration shock, as a separable forcing.
%
%   shock = SHOCK_FORCING(Struct, cfg)
%
% The forcing is a fixed spatial vector times a scalar function of time.
% Keeping the two factors apart lets every ROM project the spatial part ONCE,
% instead of redoing a dense (r x n_dofs) product at each ODE evaluation.
%
% cfg fields used: impulse_g, impulse_dir, t_shock, g_value (acceleration of
% gravity in the units of the mesh), dt, tmax, output_stride.
%
% OUTPUT
%   shock.F_spatial  [n_free x 1] spatial part, M * direction
%   shock.profile    @(t) scalar time profile, already scaled by the amplitude
%   shock.handle     @(t) F_spatial * profile(t)
%   shock.Eref       reference energy for the energy-weighted AbsTol
%   shock.t_out      output grid shared by every model, so that the
%                    post-processing compares histories without interpolating
%   shock.q0, .qd0   initial conditions (rest)
%
% See also TRANSIENTSOLVERODE.

    Mc     = Struct.AssemblyObj.constrain_matrix(Struct.M);
    n_dofs = size(Mc, 1);
    nD     = Struct.MeshObj.nDOFPerNode;

    if numel(cfg.impulse_dir) ~= nD
        error('shock_forcing:BadDirection', ...
            'impulse_dir has %d components but the model has %d DOFs per node.', ...
            numel(cfg.impulse_dir), nD);
    end

    dir_unit = cfg.impulse_dir(:) / norm(cfg.impulse_dir);
    dir_vec  = zeros(n_dofs, 1);
    for d = 1:nD
        dir_vec(d:nD:n_dofs) = dir_unit(d);
    end

    amp = cfg.impulse_g * cfg.g_value;

    shock            = struct();
    shock.direction  = dir_unit;
    shock.amplitude  = amp;
    shock.F_spatial  = Mc * dir_vec;
    shock.profile    = @(t) amp * sin(pi*t/cfg.t_shock) * (t <= cfg.t_shock);
    shock.handle     = @(t) shock.F_spatial * shock.profile(t);
    shock.q0         = zeros(n_dofs, 1);
    shock.qd0        = zeros(n_dofs, 1);
    shock.t_out      = 0 : cfg.output_stride*cfg.dt : cfg.tmax;

    % Reference energy: the shock imparts a velocity to the effective mass
    % along the impulse direction, and half its kinetic energy is the scalar
    % the tolerances are weighted on.
    v_max      = amp * 2 * cfg.t_shock / pi;
    m_eff      = full(dir_vec' * Mc * dir_vec);
    shock.Eref = 0.5 * m_eff * v_max^2;

    fprintf('Impulse: %.1e g along [%s]\n', cfg.impulse_g, num2str(dir_unit', ' %.3g'));
    fprintf('Eref = %.4e   (v_max = %.4g)\n', shock.Eref, v_max);
    fprintf('FOM: %d DOFs, %d output instants\n', n_dofs, numel(shock.t_out));
end
