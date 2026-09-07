function plot_contact_nodes(Struct, Interfaces, labels, y_contact)
%PLOT_CONTACT_NODES The contact nodes on the mesh, marked by what they did.
%
%   PLOT_CONTACT_NODES(Struct, Interfaces, labels)
%   PLOT_CONTACT_NODES(Struct, Interfaces, labels, y_contact)
%
% Without y_contact it simply shows where the contact nodes are, one marker
% per interface. With it, each node is coloured by its maximum penetration and
% the ones that never touched are left hollow, which makes the flat/corner
% verdict immediate: a face landing parallel comes out uniformly coloured, a
% face touching on an edge shows a single lit line with the rest hollow.
%
% This is the spatial counterpart of the raster that contact_activity draws
% with 'Plot', true, which shows the same events against time instead.
%
% INPUT
%   Struct      FeStructure, already built
%   Interfaces  interface metadata, needing .nodes and .gap_nodes (or .gap)
%   labels      cell array of the interfaces to draw
%   y_contact   optional, from extract_contact_response
%
% See also CONTACT_ACTIVITY, EXTRACT_CONTACT_RESPONSE, FESTRUCTURE/PLOT_NODE_SETS.

    if nargin < 4, y_contact = []; end
    if ischar(labels) || isstring(labels), labels = cellstr(labels); end

    figure('Name', 'Contact nodes', 'Color', 'w');
    PlotMesh(Struct.nodes, Struct.plot_connectivity(), 0);
    hold on;

    is3d = Struct.n_dim == 3;
    mk   = {'o', 's', 'd', '^', 'v', '>', '<', 'p'};
    h    = gobjects(1, 2*numel(labels));
    k    = 0;
    pen_all = [];

    % First pass: the peak penetration of every node, so that one colour scale
    % serves every interface and the faces stay comparable.
    peak = cell(1, numel(labels));
    for i = 1:numel(labels)
        lbl = labels{i};
        if isempty(y_contact) || ~isfield(y_contact, lbl)
            peak{i} = [];
            continue
        end
        g = gap_of(Interfaces.(lbl));
        peak{i} = max(y_contact.(lbl).normal - g, [], 2);
        pen_all = [pen_all; peak{i}(peak{i} > 0)]; %#ok<AGROW>
    end
    has_pen = ~isempty(pen_all);

    for i = 1:numel(labels)
        lbl = labels{i};
        ids = Interfaces.(lbl).nodes(:);
        c   = Struct.nodes(ids, :);
        m   = mk{mod(i-1, numel(mk)) + 1};

        if isempty(peak{i})
            k = k + 1;
            h(k) = draw(c, is3d, m, [], sprintf('%s (%d)', lbl, numel(ids)));
            continue
        end

        touched = peak{i} > 0;

        % Never touched: hollow, so they read as present but idle.
        if any(~touched)
            k = k + 1;
            h(k) = draw(c(~touched,:), is3d, m, [], ...
                sprintf('%s: never (%d)', lbl, nnz(~touched)));
            set(h(k), 'MarkerFaceColor', 'none', 'MarkerEdgeColor', [.55 .58 .62]);
        end

        % Touched: filled, coloured by peak penetration.
        if any(touched)
            k = k + 1;
            h(k) = draw(c(touched,:), is3d, m, peak{i}(touched), ...
                sprintf('%s: touched (%d)', lbl, nnz(touched)));
        end
    end

    if has_pen
        cb = colorbar;
        cb.Label.String = 'peak penetration';
        clim([0 max(pen_all)]);
    end

    legend(h(isgraphics(h)), 'Location', 'bestoutside');
    title('Contact nodes'); axis equal; grid on;
    if is3d, view(3); end
    xlabel('x'); ylabel('y'); if is3d, zlabel('z'); end
end

% =====================================================================
function h = draw(c, is3d, marker, cdata, name)
%DRAW One group of nodes, scattered or plotted depending on the colouring.
    if isempty(cdata)
        if is3d
            h = plot3(c(:,1), c(:,2), c(:,3), marker);
        else
            h = plot(c(:,1), c(:,2), marker);
        end
        set(h, 'MarkerSize', 7, 'LineWidth', 1.3, 'LineStyle', 'none', ...
            'DisplayName', name);
    else
        if is3d
            h = scatter3(c(:,1), c(:,2), c(:,3), 55, cdata, 'filled', marker);
        else
            h = scatter(c(:,1), c(:,2), 55, cdata, 'filled', marker);
        end
        set(h, 'MarkerEdgeColor', 'k', 'DisplayName', name);
    end
end

% =====================================================================
function g = gap_of(I)
%GAP_OF Per-node gaps, falling back on the scalar when a run stored only that.
    if isfield(I, 'gap_nodes') && ~isempty(I.gap_nodes)
        g = I.gap_nodes(:);
    else
        g = abs(I.gap);
    end
end
