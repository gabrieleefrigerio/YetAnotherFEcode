function run_fom(Struct, contact, shock, cfg, run_dir)
%RUN_FOM Reference full-order runs, one per quality factor and penalty.
%
%   RUN_FOM(Struct, contact, shock, cfg, run_dir)
%
% Also runs contact_activity on each result and stores it alongside: knowing
% whether the contact is flat or on a corner is what decides whether a
% reduction study is worth starting, so it should never be an afterthought.
%
% See also TRANSIENTSOLVERODE, CONTACT_ACTIVITY, EXTRACT_CONTACT_RESPONSE.

    fprintf('\n=========================================\n');
    fprintf('                 FOM\n');
    fprintf('=========================================\n');

    Mc = Struct.AssemblyObj.constrain_matrix(Struct.M);
    Kc = Struct.AssemblyObj.constrain_matrix(Struct.K);
    labels = fieldnames(contact.Interfaces)';

    for Q = cfg.array_QFactor
        Struct.compute_rayleigh_damping(Q, Q);
        Cc = Struct.AssemblyObj.constrain_matrix(Struct.C);

        for k_mult = cfg.array_k_mult
            k_contact = contact.k_base * k_mult;
            fprintf('\n[FOM] Q = %d | k_mult = %g\n', Q, k_mult);

            tic;
            solver = TransientSolverOde(Mc, Kc, Cc);
            [t, q] = solver.solve(cfg.tmax, cfg.dt, shock.q0, shock.qd0, shock.handle, ...
                'ContactOperator', contact.N, ...
                'ContactGap',      contact.gaps, ...
                'ContactPenalty',  k_contact, ...
                'Label',           'FOM', ...
                'Eref',            shock.Eref, ...
                'RelTol',          cfg.RelTolFOM, ...
                'OutputTimes',     shock.t_out);
            cpu_time = toc;

            y_contact = extract_contact_response(Struct, contact.Interfaces, labels, q);

            fprintf('\n--- Contact activity ---\n');
            activity = contact_activity(y_contact, contact.Interfaces, labels, t); %#ok<NASGU>

            Interfaces   = contact.Interfaces; %#ok<NASGU>
            model        = 'FOM';              %#ok<NASGU>
            n_modes      = size(Mc, 1);        %#ok<NASGU>
            offline_time = 0;                  %#ok<NASGU>
            save(fullfile(run_dir, sprintf('FOM_Q%04d_K%g.mat', Q, k_mult)), ...
                't', 'y_contact', 'Interfaces', 'cpu_time', 'offline_time', ...
                'model', 'n_modes', 'Q', 'k_mult', 'activity');
        end
    end
end
