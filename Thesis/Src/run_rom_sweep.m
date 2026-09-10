function run_rom_sweep(Struct, contact, shock, cfg, run_dir)
%RUN_ROM_SWEEP Reduced-order runs over every requested combination.
%
%   RUN_ROM_SWEEP(Struct, contact, shock, cfg, run_dir)
%
% Sweeps quality factor, penalty multiplier, number of linear modes, method
% and number of CC modes, and writes one .mat per combination.
%
% The contact needs no special case per method. Whatever the basis, the
% operator in the solved coordinates is
%
%       N_run = N * Pc          (and N * Pc * T_CC after interface reduction)
%
% with Pc mapping the reduced coordinates to the physical DOFs. Gaps and
% penalty therefore stay physical, and a scaled basis such as Rubin's needs no
% rescaling of its own: its scaling already sits inside Pc.
%
% See also BUILD_ROM, INTERFACE_REDUCTION, TRANSIENTSOLVERODE,
% TRANSIENTSOLVERMASSLESS.

    rom_list = {'MT', 'MC', 'CB', 'Rubin', 'MCB', 'MN'};
    rom_list = rom_list(cellfun(@(m) cfg.run.(m), rom_list) == 1);
    if isempty(rom_list), return; end

    % Methods whose basis keeps the interface as physical coordinates at the
    % head of the reduced vector, and which can therefore take a secondary
    % interface reduction. The massless ones are excluded because their M_bb
    % is zero, so the CC eigenproblem has no mass matrix to work with.
    ir_eligible = {'CB', 'Rubin'};

    fprintf('\n=========================================\n');
    fprintf('                 ROM\n');
    fprintf('=========================================\n');

    max_phi = max(cfg.array_linModes);
    labels  = fieldnames(contact.Interfaces)';

    % The Guyan interface pencil does not depend on the method or on the
    % sweep, so it is computed once here rather than inside the loops.
    Kbb_guyan = []; Mbb_guyan = [];
    if cfg.interface_reduction.enabled && strcmpi(cfg.interface_reduction.basis, 'guyan')
        fprintf('Computing the Guyan interface pencil (%d constraint modes)...\n', ...
            contact.n_bnd);
        tic_g = tic;
        [Kbb_guyan, Mbb_guyan] = guyan_interface_pencil(Struct, contact.dofs);
        fprintf('  done in %.1f s\n', toc(tic_g));
    end

    [Q_pairs, f_anchor] = damping_spec(cfg);
    for i_Q = 1:size(Q_pairs, 2)
        Q = Q_pairs(1, i_Q);            % also the tag carried by the file name

        % compute_rayleigh_damping does two things: it updates Struct.C, used
        % by the penalty ROMs, and it returns alpha and beta, from which the
        % massless ROMs build an equivalent diagonal modal damping. Using the
        % same pair keeps the damping identical across all methods.
        [~, alpha_ray, beta_ray] = Struct.compute_rayleigh_damping( ...
            Q_pairs(1, i_Q), Q_pairs(2, i_Q), f_anchor);
        rayleigh = struct('alpha', alpha_ray, 'beta', beta_ray);

        % Rayleigh fixes alpha and beta from two anchors and extrapolates the
        % rest, so the damping at the highest retained mode is a consequence,
        % never a choice. Worth printing before reading any convergence plot,
        % because it affects every method equally.
        w_hi = 2*pi * Struct.frequencies(max_phi);
        z_hi = 0.5*(alpha_ray/w_hi + beta_ray*w_hi);
        fprintf('  zeta(mode %d, %.4g Hz) = %.4f | requested zeta at the anchors = %.4f, %.4f\n', ...
            max_phi, Struct.frequencies(max_phi), z_hi, ...
            1/(2*Q_pairs(1, i_Q)), 1/(2*Q_pairs(2, i_Q)));
        if z_hi > 1
            warning('ROMSWEEP:Overdamped', ...
                'Rayleigh makes the high modes OVERDAMPED (zeta_%d = %.2f).', ...
                max_phi, z_hi);
        end

        for k_mult = cfg.array_k_mult
            k_contact = contact.k_base * k_mult;

            for phi = cfg.array_linModes
                for im = 1:numel(rom_list)
                    model = rom_list{im};

                    % The massless models use exact set-valued contact, so the
                    % penalty does not apply to them and the k_mult loop would
                    % be degenerate. They run at the first value only.
                    is_massless = any(strcmp(model, {'MCB', 'MN'}));
                    if is_massless && k_mult ~= cfg.array_k_mult(1), continue; end

                    fprintf('\n--- %s | Phi %d | Q %d | k_mult %g ---\n', ...
                        model, phi, Q, k_mult);

                    tic_offline = tic;
                    [rom, Pc, Mr, Kr, Cr] = build_rom(Struct, model, phi, ...
                        contact, k_contact, rayleigh);
                    offline_base = toc(tic_offline);

                    % Spatial part of the forcing projected once: the time
                    % dependence is a scalar, so there is no need to redo this
                    % product at every ODE function evaluation.
                    F_r_full = Pc' * shock.F_spatial;

                    % n_cc = 0 is the baseline without interface reduction,
                    % always run so the comparison term is present.
                    if cfg.interface_reduction.enabled && any(strcmp(model, ir_eligible))
                        cc_list = [0, cfg.array_ccModes];
                    else
                        cc_list = 0;
                    end

                    for n_cc = cc_list
                        ir = apply_interface_reduction(Mr, Kr, Cr, F_r_full, ...
                            rom, model, n_cc, cfg, contact, Pc, ...
                            Kbb_guyan, Mbb_guyan, offline_base);

                        fprintf('  [%s] reduced size %d\n', ir.tag, size(ir.K, 1));

                        z0       = zeros(size(ir.K, 1), 1);
                        F_handle = @(t) ir.F * shock.profile(t);

                        tic;
                        if is_massless
                            [t, q_rom, lambda, info] = integrate_massless( ...
                                ir, rom, contact, Pc, z0, F_handle, cfg);
                        else
                            [t, q_rom] = integrate_penalty(ir, rom, contact, ...
                                Pc, n_cc, z0, F_handle, shock, cfg, k_contact);
                            lambda = []; info = [];
                        end
                        cpu_time = toc;

                        % Back through the interface transformation first and
                        % only then through Pc: applying T (a few hundred rows)
                        % to the time history is far cheaper than forming Pc*T.
                        % T comes from the reduction rather than being rebuilt
                        % as blkdiag(Phi_CC, I) here, because the refined basis
                        % is NOT block diagonal - it couples the interface to
                        % the substructure modes.
                        if n_cc == 0
                            q_full = q_rom;
                        else
                            q_full = ir.info.T * q_rom;
                        end
                        y_contact = extract_contact_response(Struct, ...
                            contact.Interfaces, labels, Pc * q_full);

                        save_rom_result(run_dir, contact, model, ir, phi, ...
                            n_cc, Q, k_mult, t, y_contact, cpu_time, ...
                            lambda, info);
                    end
                end
            end
        end
    end
end

% =====================================================================
function ir = apply_interface_reduction(Mr, Kr, Cr, F_r_full, rom, model, ...
                                        n_cc, cfg, contact, Pc, ...
                                        Kbb_guyan, Mbb_guyan, t_off)
%APPLY_INTERFACE_REDUCTION Secondary modal reduction of the interface block.
% Returns the untouched matrices when n_cc is zero, which is the baseline.

    ir = struct('M', Mr, 'K', Kr, 'C', Cr, 'F', F_r_full, 'Phi_CC', [], ...
                'info', [], 'mode', 'none', 'tag', model, 'offline', t_off);
    if n_cc == 0, return; end

    tic_ir = tic;

    % The Guyan pencil is known in PHYSICAL interface coordinates, while
    % interface_reduction works in the ROM's own interface coordinates. With
    % u_phys = A*u_rom (A being the contact rows of the boundary block of Pc,
    % the identity for CB and the scaling for Rubin), the congruence A'*(.)*A
    % moves the pencil across and leaves the eigenvalues unchanged, so Phi_CC
    % comes back directly in ROM coordinates.
    if strcmpi(cfg.interface_reduction.basis, 'guyan')
        A = Pc(contact.dofs, 1:rom.n_bnd);
        pencil = struct('K', A' * Kbb_guyan * A, 'M', A' * Mbb_guyan * A);
    else
        pencil = [];
    end

    do_static = isfield(cfg.interface_reduction, 'static_correction') && ...
                cfg.interface_reduction.static_correction;

    % Refinement of the CC modes (Ahn et al. 2022). Independent of the contact
    % compliance above: that one changes the contact law, this one changes the
    % basis. They compose, and either can be used alone.
    do_refine = isfield(cfg.interface_reduction, 'refine') && ...
                cfg.interface_reduction.refine;

    % Optional: give every contact face the same number of CC modes, instead of
    % letting the frequency pooling decide the split. n_cc is distributed as
    % evenly as possible over the faces (exactly equal when it is a multiple of
    % their number), so "k modes per interface" is array_ccModes = k*n_faces.
    cc_alloc = [];
    equal_pf = isfield(cfg.interface_reduction, 'equal_per_face') && ...
               cfg.interface_reduction.equal_per_face;
    if equal_pf && strcmpi(cfg.interface_reduction.mode, 'per_interface')
        nf = numel(contact.blocks);
        base = floor(n_cc / nf);
        rem  = n_cc - base*nf;                 % spread the remainder one-each
        cc_alloc = base * ones(1, nf);
        cc_alloc(1:rem) = cc_alloc(1:rem) + 1; % faces are ordered as in contact.blocks
    end

    [ir.M, ir.K, ir.C, ir.Phi_CC, ir.info] = interface_reduction( ...
        Mr, Kr, Cr, rom.n_bnd, n_cc, cfg.interface_reduction.mode, ...
        contact.blocks, pencil, do_static, cc_alloc, do_refine);

    ir.mode = ir.info.mode;
    switch ir.mode
        case 'global',        ir.tag = [model 'IRG'];
        case 'per_interface', ir.tag = [model 'IRP'];
        otherwise,            ir.tag = [model 'IR'];
    end
    % A correction changes the model, so it must change the file name too:
    % without this a corrected and an uncorrected run write the SAME file and
    % the second silently overwrites the first. The suffixes are letters only,
    % which is what the post-processing method parser accepts, so each variant
    % appears as its own method next to the plain one. They can both be on:
    %   SC  series contact compliance (corrects the contact law)
    %   RF  refined CC modes (corrects the basis)
    if ~isempty(ir.info)
        if isfield(ir.info, 'R_res') && ~isempty(ir.info.R_res)
            ir.tag = [ir.tag 'SC'];
        end
        if isfield(ir.info, 'refined') && ir.info.refined
            ir.tag = [ir.tag 'RF'];
        end
    end

    % T' applied to the projected forcing. Taken from the reduction, not rebuilt
    % from Phi_CC: the refined transformation carries an interface/substructure
    % coupling block that blkdiag(Phi_CC, I) does not have.
    ir.F = ir.info.T' * F_r_full;
    ir.offline = t_off + toc(tic_ir);
end

% =====================================================================
function [t, q] = integrate_penalty(ir, rom, contact, Pc, n_cc, z0, ...
                                    F_handle, shock, cfg, k_contact)
%INTEGRATE_PENALTY Penalty contact with ode15s, one path for every method.

    N_run = contact.N * Pc;

    % Operator mapping the ROM's INTERFACE coordinates onto the constraints,
    % captured before the CC reduction is applied to it. It is the right map for
    % the compliance below: R_res comes out of the interface pencil, so it lives
    % in the ROM's interface coordinates, which for a scaled basis such as
    % Rubin's are NOT the physical ones. Using the physical contact operator
    % here silently produced a compliance twelve orders of magnitude too small.
    Nb_rom = N_run(:, 1:rom.n_bnd);

    if n_cc > 0
        % Right-multiplication by the interface transformation. Taken from the
        % reduction rather than assumed block diagonal: after refinement it is
        % not, so N_run(:,1:n_bnd)*Phi_CC would drop the coupling block.
        N_run = N_run * ir.info.T;
    end

    % Compliance of the interface content removed by the truncation, mapped onto
    % the contact constraints: one row and column per constraint.
    Scomp = [];
    if n_cc > 0 && isfield(ir.info, 'R_res') && ~isempty(ir.info.R_res)
        Scomp = full(Nb_rom * ir.info.R_res * Nb_rom');
        Scomp = (Scomp + Scomp') / 2;
    end

    solver = TransientSolverOde(ir.M, ir.K, ir.C);
    [t, q] = solver.solve(cfg.tmax, cfg.dt, z0, z0, F_handle, ...
        'ContactCompliance', Scomp, ...
        'ContactOperator', N_run, ...
        'ContactGap',      contact.gaps, ...
        'ContactPenalty',  k_contact, ...
        'Label',           ir.tag, ...
        'Eref',            shock.Eref, ...
        'RelTol',          cfg.RelTol, ...
        'OutputTimes',     shock.t_out);
end

% =====================================================================
function [t, q, lambda, info] = integrate_massless(ir, rom, contact, Pc, ...
                                                   z0, F_handle, cfg)
%INTEGRATE_MASSLESS Exact set-valued contact (Monjaraz-Tec et al. 2022).
% Solver convention: g = g0 + W'*q_b, contact when g <= 0. The penetration of
% the operator is the opposite of that gap, hence W = -N_b' and g0 = gaps.

    n_bnd = rom.n_bnd;
    if n_bnd ~= numel(contact.gaps)
        error('ROMSWEEP:BndMismatch', ...
            'n_bnd = %d but there are %d contact constraints.', ...
            n_bnd, numel(contact.gaps));
    end

    % This scheme solves the contact on the static boundary partition alone,
    % so the operator must not reach the modal coordinates.
    N_rom = contact.N * Pc;
    reach = norm(full(N_rom(:, n_bnd+1:end)), 'fro');
    if reach > 1e-8 * norm(full(N_rom(:, 1:n_bnd)), 'fro')
        error('ROMSWEEP:ContactOnModal', ...
            ['The contact operator reaches the modal block (||.|| = %.3e). ' ...
             'The massless formulation requires the contact to act on the ' ...
             'boundary partition only.'], reach);
    end
    W = -N_rom(:, 1:n_bnd)';

    solver = TransientSolverMassless(ir.M, ir.K, ir.C, n_bnd, W, contact.gaps);
    [t, q, lambda, info] = solver.solve(cfg.tmax, cfg.dt, z0, z0, F_handle);

    % The massless solver integrates at fixed dt and knows nothing about
    % OutputTimes, so its output is brought onto the common grid. Since that
    % grid has a step of output_stride*dt, its instants are an EXACT subset of
    % the solver grid: sampling introduces no interpolation.
    idx    = 1 : cfg.output_stride : numel(t);
    t      = t(idx);
    q      = q(:, idx);
    lambda = lambda(:, idx);
end

% =====================================================================
function save_rom_result(run_dir, contact, model, ir, phi, n_cc, Q, k_mult, ...
                         t, y_contact, cpu_time, lambda, info)
%SAVE_ROM_RESULT One file per combination, named after it.
% The variable names are the contract with the post-processing.

    Interfaces   = contact.Interfaces;
    model_tag    = ir.tag;
    ir_mode      = ir.mode;
    ir_info      = ir.info;
    offline_time = ir.offline;
    n_modes      = phi;

    if n_cc == 0
        fname = sprintf('ROM_%s_Phi%03d_Q%g_K%g.mat', ir.tag, phi, Q, k_mult);
    else
        fname = sprintf('ROM_%s_Phi%03d_CC%03d_Q%g_K%g.mat', ...
            ir.tag, phi, n_cc, Q, k_mult);
    end

    save(fullfile(run_dir, fname), ...
        't', 'y_contact', 'Interfaces', 'cpu_time', 'offline_time', ...
        'model', 'model_tag', 'n_modes', 'Q', 'k_mult', 'lambda', 'info', ...
        'n_cc', 'ir_mode', 'ir_info');
end
