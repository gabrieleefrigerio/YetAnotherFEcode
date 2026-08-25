function act = contact_activity(y_contact, Interfaces, labels, t, varargin)
%CONTACT_ACTIVITY What the contact actually does, before any reduction study.
%
%   act = CONTACT_ACTIVITY(y_contact, Interfaces, labels, t)
%   act = CONTACT_ACTIVITY(..., 'Plot', true, 'Coord', 1)
%
% Answers, per interface, the questions that decide whether reducing the
% interface can work at all:
%
%   how many nodes ever touch, and how many at the same time
%   when the first contact happens
%   what fraction of the time the interface is loaded (duty cycle)
%   whether the contact is FLAT or on a CORNER
%
% The last one is the important one. On the 2D dummy model the same mesh, the
% same walls and the same shock amplitude gave 81 active nodes out of 81 with
% a shock at 0 degrees, and 4 out of 81 at 45 degrees, and the interface
% reduction degraded 11800 times more in the second case. A flat contact loads
% the face with a smooth, low-wavenumber profile that a handful of CC modes
% represent well; a corner contact concentrates the force on a few nodes,
% which is a high-wavenumber profile, and truncating the interface basis makes
% the face infinitely stiff in exactly the directions the load needs. Run this
% before any n_cc sweep: if the contact is on a corner, the sweep will only
% measure how fast the model degrades.
%
% The flat/corner indicator is the fraction of the face loaded at the same
% time, max(simultaneously active) / n_nodes: at or above one half the face
% lands essentially parallel to the wall, well below it touches on an edge.
% A supporting figure, spread_ratio, is the spread of the penetration across
% the whole face at the deepest instant, normalized by the maximum: order 1
% for a face landing parallel, orders of magnitude more when the far end
% stays wide open.
%
% INPUT
%   y_contact   as returned by extract_contact_response, needing the field
%               .normal (older runs without it fall back on the axis given
%               by Interfaces.<label>.dir)
%   Interfaces  interface metadata from the main
%   labels      cell array of the interfaces to analyse
%   t           time vector [1 x n_time], in seconds
%
% OPTIONS
%   'Plot'   false by default. Draws the activity raster, one row per node
%            ordered along the face, plus the count of simultaneously active
%            nodes.
%   'Coord'  coordinate used to order the nodes along the face in the raster.
%            Default: the one with the largest spread.
%
% OUTPUT  act.<label> with fields
%   n_nodes, n_ever, n_max_simultaneous, active_fraction, first_contact,
%   duty_cycle, max_penetration, spread_ratio,
%   regime ('flat' | 'corner' | 'none')
%
% See also EXTRACT_CONTACT_RESPONSE, CONTACT_OPERATOR, INTERFACE_REDUCTION.

    p = inputParser;
    addParameter(p, 'Plot', false);
    addParameter(p, 'Coord', []);
    parse(p, varargin{:});
    opt = p.Results;

    t = t(:)';
    act = struct();
    n_if = numel(labels);

    if opt.Plot
        fig = figure('Color', 'w', 'Position', [80 60 1100 240*ceil((n_if+1)/2)]);
        tiledlayout(fig, ceil((n_if+1)/2), 2, 'TileSpacing', 'compact', ...
            'Padding', 'compact');
    end
    total = zeros(size(t));

    for i = 1:n_if
        lbl = labels{i};
        I   = Interfaces.(lbl);
        yc  = y_contact.(lbl);

        if isfield(yc, 'normal')
            un = yc.normal;
        else
            ax = {'X', 'Y', 'Z'};
            un = sign(I.gap) * yc.(ax{I.dir});
        end

        % Gap as a positive distance, whatever convention the run stored.
        g = abs(I.gap);
        if isfield(I, 'gap_nodes') && ~isempty(I.gap_nodes)
            g = I.gap_nodes(:);
        end

        pen  = un - g;                    % > 0 where the node has crossed
        hit  = pen > 0;
        anyt = any(hit, 1);
        total = total + sum(hit, 1);

        a.n_nodes            = size(un, 1);
        a.n_ever             = nnz(any(hit, 2));
        a.n_max_simultaneous = max(sum(hit, 1));
        a.duty_cycle         = mean(anyt);
        a.max_penetration    = max(max(pen(:)), 0);

        % Flat or corner. The discriminator is how much of the face carries
        % the load at the same time, because that is what sets the wavenumber
        % content of the contact force profile, and hence whether a few CC
        % modes can represent it. Measuring the tilt of the face instead would
        % be misleading: the spread of the displacement over the face is
        % dominated by the rigid translation of the whole body, which says
        % nothing about the contact.
        a.active_fraction = a.n_max_simultaneous / a.n_nodes;

        % Supporting figure: spread of the penetration across the WHOLE face,
        % normalized by its maximum, at the instant of deepest penetration. A
        % face landing parallel gives order 1, a face touching on an edge
        % gives orders of magnitude more, the far end being wide open.
        if any(anyt)
            [~, j_peak] = max(max(pen, [], 1));
            pk = pen(:, j_peak);
            a.spread_ratio  = (max(pk) - min(pk)) / max(max(pk), realmin);
            a.first_contact = t(find(anyt, 1));
            if a.active_fraction >= 0.5
                a.regime = 'flat';
            else
                a.regime = 'corner';
            end
        else
            a.spread_ratio  = NaN;
            a.first_contact = NaN;
            a.regime        = 'none';
        end
        act.(lbl) = a;

        fprintf(['%-6s %4d nodes | touching %4d | max at once %4d (%5.1f%%) | ' ...
                 'first %8.3f us | active %5.2f%% | spread %8.3g -> %s\n'], ...
            lbl, a.n_nodes, a.n_ever, a.n_max_simultaneous, ...
            100*a.active_fraction, 1e6*a.first_contact, 100*a.duty_cycle, ...
            a.spread_ratio, upper(a.regime));

        if opt.Plot
            s = pick_coord(I, opt.Coord);
            [ss, ord] = sort(s);
            hitp = hit(ord, :);
            nexttile; hold on; grid on; box on
            [ni, ti] = find(hitp);
            if ~isempty(ni)
                scatter(1e6*t(ti), ss(ni), 6, 'filled', ...
                    'MarkerFaceColor', [0.85 0.15 0.15]);
            end
            if numel(ss) > 1, ylim([min(ss) max(ss)]); end
            xlim(1e6*[t(1) t(end)]);
            % The coordinates are left as they are: the mesh may be in metres
            % (2D case) or in micrometres (TDK mesh). Time instead is in
            % seconds in both, so the microsecond scale always applies.
            xlabel('t [\mus]'); ylabel('position along the face');
            title(sprintf('%s  -  %d/%d nodes  -  %s', lbl, a.n_ever, ...
                a.n_nodes, upper(a.regime)));
        end
    end

    if opt.Plot
        nexttile([1 2]); hold on; grid on; box on
        plot(1e6*t, total, 'k-', 'LineWidth', 1.4);
        xlabel('t [\mus]'); ylabel('nodes in contact');
        title(sprintf('Simultaneously active nodes (max %d)', max(total)));
        xlim(1e6*[t(1) t(end)]);
    end
end

% =====================================================================
function s = pick_coord(I, forced)
%PICK_COORD Coordinate along which the nodes of a face are ordered.
% By default the one with the largest spread, which for a flat face is a
% direction lying in the face rather than the one normal to it.
    C = [];
    for f = {'coord_X', 'coord_Y', 'coord_Z'}
        if isfield(I, f{1}), C = [C, I.(f{1})(:)]; end %#ok<AGROW>
    end
    if isempty(C)
        s = (1:numel(I.nodes))';
        return
    end
    if ~isempty(forced)
        s = C(:, forced);
    else
        [~, d] = max(max(C, [], 1) - min(C, [], 1));
        s = C(:, d);
    end
end
