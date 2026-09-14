function res = contact_resolution(y_contact, Interfaces, labels, t, h)
%CONTACT_RESOLUTION Is the time step small enough for the impacts in this run?
%
%   res = CONTACT_RESOLUTION(y_contact, Interfaces, labels, t, h)
%
% Answers, from ONE finished run and with no reference to compare against,
% the question a fixed-step integrator cannot answer by itself: whether h
% resolves the contact events or steps over them.
%
% A Newton iteration that converges says the algebraic problem of that step
% was solved. It says NOTHING about h being small enough, because the scheme
% will happily converge onto a step that jumped across an entire impact. The
% quantity that fails first is not the trajectory but the IMPACT COUNT: on
% this model the events last about 60 ns, and once fewer than roughly twelve
% steps fall inside one, adjacent impacts start merging into one.
%
% INPUT
%   y_contact   the per-interface response struct
%   Interfaces  the metadata struct from run_config
%   labels      interface names to examine
%   t           output time grid [s]
%   h           the integration step [s] (NOT the output spacing)
%
% OUTPUT
%   res.n_events        events detected, per interface and total
%   res.dur_median      median event duration [s]
%   res.steps_median    median integration steps per event
%   res.steps_min       the worst resolved event
%   res.thin            events resolved by fewer than 12 steps
%   res.verdict         'adequate' | 'marginal' | 'inadequate'
%
% The output grid must also resolve the events, otherwise the durations
% measured here are quantised by the sampling rather than by the physics; the
% function says so when that is the case.
%
% See also CONTACT_ACTIVITY, TRANSIENTSOLVERNEWMARK.

    THIN = 12;                     % steps per event below which impacts merge
    dt_out = t(2) - t(1);

    dur = []; per_face = zeros(1, numel(labels));
    for i = 1:numel(labels)
        g  = Interfaces.(labels{i}).gap_nodes(1);
        p  = max(y_contact.(labels{i}).normal, [], 1) - g;
        act = p > 0;
        d = diff([false act false]);
        i0 = find(d == 1);  i1 = find(d == -1) - 1;
        per_face(i) = numel(i0);
        dur = [dur, (i1 - i0 + 1) * dt_out]; %#ok<AGROW>
    end

    res.n_events     = sum(per_face);
    res.per_face     = per_face;
    res.dur_median   = median(dur);
    res.dur_min      = min(dur);
    res.steps_median = res.dur_median / h;
    res.steps_min    = res.dur_min / h;
    res.thin         = nnz(dur / h < THIN);
    res.samples_median = res.dur_median / dt_out;

    if isempty(dur)
        res.verdict = 'no contact';
    elseif res.steps_min >= THIN
        res.verdict = 'adequate';
    elseif res.steps_median >= THIN
        res.verdict = 'marginal';
    else
        res.verdict = 'inadequate';
    end

    fprintf('\n--- contact resolution ---\n');
    fprintf('  events %d over %d interfaces | step h = %.3g ns\n', ...
        res.n_events, numel(labels), h*1e9);
    fprintf('  duration: median %.1f ns, shortest %.1f ns\n', ...
        res.dur_median*1e9, res.dur_min*1e9);
    fprintf('  integration steps per event: median %.1f, worst %.1f\n', ...
        res.steps_median, res.steps_min);
    fprintf('  events under %d steps: %d of %d\n', THIN, res.thin, res.n_events);
    fprintf('  VERDICT: %s\n', upper(res.verdict));
    if res.samples_median < 4
        fprintf(['  [!] the OUTPUT grid samples an event %.1f times: the durations\n' ...
                 '      above are quantised by the sampling, not by the physics.\n' ...
                 '      Lower cfg.output_stride before trusting them.\n'], ...
                 res.samples_median);
    end
end
