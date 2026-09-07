function trk = rom_tracking_time(t, err_pct, thresholds, varargin)
%ROM_TRACKING_TIME For how long does a ROM stay within a given error?
%
%   trk = ROM_TRACKING_TIME(t, err_pct, thresholds)
%   trk = ROM_TRACKING_TIME(..., 'HoldFraction', 0.01, 'Label', 'CB Phi200')
%
% A global relative error compresses a whole time history into one number,
% and on a vibro-impact problem that number depends on a window chosen by
% hand: too long and it saturates on the chaotic decorrelation, too short and
% it only measures the first impact. This function asks the engineering
% question instead - up to when can this ROM be trusted - which is monotone,
% does not saturate, and needs no window.
%
% Two answers are returned, because they are not the same question:
%
%   T        the trustworthy horizon: the first instant at which the error
%            crosses the threshold AND STAYS above it. The hold requirement
%            is what keeps a single spike at a zero crossing from cutting the
%            horizon in half.
%
%   frac     the fraction of the whole window spent below the threshold. A
%            ROM that leaves the band during one impact and comes back is
%            penalised by T but not by frac, so a large gap between the two
%            means the error is episodic rather than accumulating.
%
% CENSORING. When the error never crosses the threshold, T is the length of
% the run and censored is true. That value is a lower bound, not the horizon,
% and must be reported as "> T", never averaged in with the uncensored ones.
%
% CHOOSING THE THRESHOLDS. In the exponential regime the error grows as
% exp(lambda*t), so the horizon moves only logarithmically with the
% threshold: one decade of threshold buys ln(10)/lambda of extra time. On the
% 3D model lambda = 9.6e4 1/s was measured between two FOM runs at different
% tolerances, i.e. one decade is worth 24 us. The threshold therefore does
% not change the RANKING of the methods, but it does change the numbers a
% lot, which is why the default is a set of three rather than a single value:
%
%   0.1 %   demanding, separates the good ROMs from one another
%   1   %   engineering grade, the headline number
%   10  %   still qualitatively right, separates the broken ones
%
% The lower end is bounded by the integration floor: on the 3D model two FOMs
% at RelTol 1e-7 and 1e-8 differ by 5e-7 % at 1 us and 6e-6 % at 10 us, so a
% threshold below about 1e-4 % would be measuring the integrator, not the ROM.
%
% INPUT
%   t           [1 x n_time] time grid
%   err_pct     [1 x n_time] error in per cent, ALREADY normalised. Pass the
%               same curve that gets plotted (gre_t in the post-processing,
%               normalised by the peak of the reference) so that the number
%               and the figure can never disagree.
%   thresholds  [1 x k] per cent, default [0.1 1 10]
%
% OUTPUT
%   trk.thresholds  [1 x k]
%   trk.T           [1 x k] trustworthy horizon [s]
%   trk.censored    [1 x k] logical, true when the threshold was never crossed
%   trk.frac        [1 x k] fraction of the window spent below the threshold
%
% See also ROM_DIVERGENCE, CONTACT_ACTIVITY.

    p = inputParser;
    addParameter(p, 'HoldFraction', 0.01);
    addParameter(p, 'Label', '');
    addParameter(p, 'Verbose', false);
    parse(p, varargin{:});
    opt = p.Results;

    if nargin < 3 || isempty(thresholds)
        thresholds = [0.1 1 10];
    end

    t   = t(:)';
    err = err_pct(:)';
    if numel(t) ~= numel(err)
        error('rom_tracking_time:SizeMismatch', ...
            't has %d samples but the error curve has %d.', numel(t), numel(err));
    end

    n      = numel(t);
    n_hold = min(max(1, round(opt.HoldFraction * n)), n);
    k      = numel(thresholds);

    trk.thresholds = thresholds(:)';
    trk.T          = nan(1, k);
    trk.censored   = false(1, k);
    trk.frac       = nan(1, k);

    for i = 1:k
        above = err > thresholds(i);

        % Forward-looking count, so a crossing counts only once the error has
        % stayed above the threshold for the whole hold window.
        sustained = movsum(double(above), [0 n_hold-1]) >= n_hold;
        j = find(sustained, 1);

        if isempty(j)
            trk.T(i)        = t(end);
            trk.censored(i) = true;
        else
            trk.T(i) = t(j);
        end
        trk.frac(i) = nnz(~above) / n;
    end

    if opt.Verbose
        fprintf('--- tracking time: %s ---\n', opt.Label);
        for i = 1:k
            if trk.censored(i), pre = '>'; else, pre = ' '; end
            fprintf('  within %6.3g %%  %s%8.2f us   (below for %5.1f%% of the run)\n', ...
                trk.thresholds(i), pre, 1e6*trk.T(i), 100*trk.frac(i));
        end
    end
end
