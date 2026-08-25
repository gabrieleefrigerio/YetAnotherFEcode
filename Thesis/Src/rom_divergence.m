function div = rom_divergence(t, y_fom, y_rom, varargin)
%ROM_DIVERGENCE Is the error of this ROM still meaningful, or has it saturated?
%
%   div = ROM_DIVERGENCE(t, y_fom, y_rom)
%   div = ROM_DIVERGENCE(..., 'Plot', true, 'Label', 'CB n_cc=8')
%
% Vibro-impact is chaotic: two trajectories that differ by any amount, however
% small, separate exponentially until they decorrelate, after which their
% distance stops growing and settles at a level that depends only on the size
% of the response, not on the quality of the ROM. A global relative error
% computed over a window longer than that saturation time therefore measures
% nothing about the ROM. On the 2D dummy model this was visible directly:
% every method and every shock direction returned the same 48.1 to 48.5 per
% cent at n_cc = 8, which is the signature of a saturated metric, not of
% seven models that happen to be equally good.
%
% This function measures where that boundary is, so that an error can be
% quoted over a window where it still discriminates.
%
%   d(t) = || y_rom(t) - y_fom(t) ||   over the interface
%
%   d_sat      saturation level, the median of d over the last third
%   T_sat      first instant at which d reaches half of d_sat. Read it as
%              genuine decorrelation only when the ROM does saturate within
%              the run: for a well converged model the error never really
%              plateaus, and T_sat then marks the last large divergence
%              event instead, which is why several good models can come back
%              with almost the same value.
%   window     suggested window for a discriminating error: t <= T_sat
%   saturated  true if the run is longer than T_sat, i.e. an error computed
%              over the whole run is already contaminated
%   lambda     exponential separation rate, fitted on the growth phase, and
%              ONLY when that phase spans at least one decade of d. Otherwise
%              NaN: a model that saturates almost immediately has no
%              exponential regime to fit, and fitting one anyway produces a
%              number that changes with the window chosen, as three different
%              window choices did here, one of them even reversing the
%              ordering of the methods.
%   horizon   1/lambda, the predictability time, NaN when lambda is
%
% INPUT
%   t       [1 x n_time]
%   y_fom   [n x n_time] reference response
%   y_rom   [n x n_time] reduced response on the same grid
%
% OUTPUT  div, plus the printed summary
%
% See also CONTACT_ACTIVITY.

    p = inputParser;
    addParameter(p, 'Plot', false);
    addParameter(p, 'Label', 'ROM');
    parse(p, varargin{:});
    opt = p.Results;

    t = t(:)';
    if ~isequal(size(y_fom), size(y_rom))
        error('rom_divergence:SizeMismatch', ...
            'FOM [%s] and ROM [%s] do not have the same size.', ...
            num2str(size(y_fom)), num2str(size(y_rom)));
    end

    d     = sqrt(sum((y_rom - y_fom).^2, 1));
    scale = sqrt(sum(y_fom.^2, 1));

    n     = numel(t);
    tail  = max(1, round(2*n/3)) : n;
    d_sat = median(d(tail));

    j = find(d >= 0.5*d_sat, 1);
    if isempty(j), j = n; end
    T_sat = t(j);

    div.d           = d;
    div.d_sat       = d_sat;
    div.T_sat       = T_sat;
    div.window      = [t(1) T_sat];
    div.saturated   = t(end) > T_sat;
    div.rel_full    = norm(y_rom - y_fom, 'fro') / norm(y_fom, 'fro');
    in_win          = t <= T_sat;
    div.rel_window  = norm(y_rom(:, in_win) - y_fom(:, in_win), 'fro') / ...
                      norm(y_fom(:, in_win), 'fro');
    div.sat_level   = d_sat / max(median(scale(tail)), realmin);

    % Exponential rate, only where there is a decade to fit on.
    grow = d > 0 & t <= T_sat;
    lambda = NaN;
    if nnz(grow) > 5
        dg = d(grow);
        if max(dg)/min(dg) >= 10          % at least one decade
            cf = polyfit(t(grow), log(dg), 1);
            lambda = cf(1);
        end
    end
    div.lambda  = lambda;
    div.horizon = 1/lambda;

    fprintf('--- %s ---\n', opt.Label);
    fprintf('  saturation level     %.3e  (%.1f%% of the response)\n', ...
        d_sat, 100*div.sat_level);
    fprintf('  T_sat                %.2f us\n', 1e6*T_sat);
    if div.saturated
        verdict = 'SATURATED, the global error does not discriminate';
    else
        verdict = 'still within the deterministic regime';
    end
    fprintf('  simulated window     %.2f us -> %s\n', 1e6*t(end), verdict);
    if ~isnan(lambda)
        fprintf('  lambda %.3e 1/s | predictability horizon %.1f us\n', ...
            lambda, 1e6*div.horizon);
    else
        fprintf('  lambda not estimable: no exponential stretch spanning a decade\n');
    end
    fprintf('  relative error: whole window %.4f%% | up to T_sat %.4f%%\n', ...
        100*div.rel_full, 100*div.rel_window);

    if opt.Plot
        figure('Color', 'w', 'Position', [100 100 720 420]);
        semilogy(1e6*t, max(d, realmin), 'LineWidth', 1.3); hold on; grid on
        yline(d_sat, '--', 'saturation', 'LineWidth', 1.1);
        xline(1e6*T_sat, ':', 'T_{sat}', 'LineWidth', 1.1);
        xlabel('t [\mus]'); ylabel('||y_{ROM} - y_{FOM}||');
        title(sprintf('%s  -  divergence', opt.Label));
    end
end
