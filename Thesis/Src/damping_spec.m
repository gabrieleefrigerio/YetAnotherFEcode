function [Q_pairs, f_anchor] = damping_spec(cfg)
%DAMPING_SPEC The damping cases to run, read from a run configuration.
%
%   [Q_pairs, f_anchor] = DAMPING_SPEC(cfg)
%
% Turns the damping fields of a configuration into the arguments
% FeStructure.compute_rayleigh_damping expects, so every main and every sweep
% reads them the same way.
%
% cfg.Q_freq   [fa fb]  anchor FREQUENCIES [Hz]. Omitted -> the first two
%              natural frequencies, which is what every earlier run used.
%
% The anchors are frequencies and never mode indices: Rayleigh gives
% zeta_n = alpha/(2*w_n) + beta*w_n/2, in which the index has cancelled, so an
% index would only be a way of naming a frequency through the spectrum.
%
% HOW THE QUALITY FACTORS ARE READ is decided by whether Q_freq was given, and
% by nothing else:
%
%   Q_freq PRESENT   cfg.array_QFactor holds EXACTLY TWO values, the quality
%                    factor at the first anchor and the one at the second. Row
%                    or column makes no difference. No sweep in this form: to
%                    compare two dampings, run the main twice.
%
%   Q_freq ABSENT    cfg.array_QFactor is a sweep: one quality factor per case,
%                    applied at both default anchors. This is the legacy form
%                    and it keeps producing exactly what it used to.
%
% The orientation is deliberately ignored. An earlier version distinguished a
% [2 x N] column, one case with two anchors, from a [1 x N] row, N cases
% sharing one value: one semicolon apart, both valid, so a typo silently
% swapped the physics for a sweep.
%
% Rayleigh fixes alpha and beta from two anchors and extrapolates everywhere
% else, so on a shock problem the anchors matter more than the quality factors
% themselves: pinned on modes 1 and 2 of the 3D model they sit ten per cent
% apart and leave the whole contact band, 100 kHz to 5 MHz, unspecified.
%
% Example, the parametrisation used at TDK:
%   cfg.Q_freq        = [30e3, 5e6];
%   cfg.array_QFactor = [300, 5000];      % Q = 300 at 30 kHz, 5000 at 5 MHz
%
% OUTPUT
%   Q_pairs   [2 x N] one column per case: Q at the first anchor, Q at the
%             second. N = 1 whenever Q_freq is given.
%   f_anchor  [1 x 2] anchor frequencies [Hz], or [] for the first two modes.
%
% See also FESTRUCTURE/COMPUTE_RAYLEIGH_DAMPING.

    if ~isfield(cfg, 'array_QFactor') || isempty(cfg.array_QFactor)
        error('DAMPING_SPEC:NoQ', 'cfg.array_QFactor is missing or empty.');
    end
    QF = cfg.array_QFactor;
    if any(QF(:) <= 0) || any(~isfinite(QF(:)))
        error('DAMPING_SPEC:BadQ', 'Every quality factor must be positive and finite.');
    end

    if ~isfield(cfg, 'Q_freq') || isempty(cfg.Q_freq)
        % Legacy sweep on the first two natural frequencies.
        if ~isvector(QF)
            error('DAMPING_SPEC:BadShape', ...
                ['Without cfg.Q_freq, cfg.array_QFactor is a sweep and must be a ' ...
                 'vector, got [%d x %d].'], size(QF,1), size(QF,2));
        end
        QF = QF(:)';
        Q_pairs  = [QF; QF];
        f_anchor = [];

        % The file names carry the quality factor of each case, so two cases
        % sharing it would silently overwrite each other's results.
        if numel(unique(QF)) < numel(QF)
            error('DAMPING_SPEC:AmbiguousTag', ...
                ['Two damping cases share the same quality factor (%s): their ' ...
                 'result files would collide.'], mat2str(QF));
        end
        return
    end

    f_anchor = cfg.Q_freq(:)';
    if numel(f_anchor) ~= 2
        error('DAMPING_SPEC:BadFreq', ...
            'cfg.Q_freq must hold exactly two frequencies [Hz], got %d.', ...
            numel(f_anchor));
    end
    if numel(QF) ~= 2
        error('DAMPING_SPEC:QCountMismatch', ...
            ['With cfg.Q_freq set, cfg.array_QFactor must hold exactly two ' ...
             'quality factors, one per anchor, got %d.\n' ...
             'To compare several dampings, run the main once per case.'], numel(QF));
    end
    Q_pairs = QF(:);                     % order follows Q_freq
end
