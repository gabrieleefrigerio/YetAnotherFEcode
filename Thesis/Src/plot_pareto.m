function fig = plot_pareto(T, gre_desc, gre_floor_pct, ttl)
%PLOT_PARETO Accuracy against online cost, with every point identifiable.
%
%   fig = PLOT_PARETO(T, gre_desc, gre_floor_pct, ttl)
%
% T is the summary table (method, phi, n_cc, gre_ref, cpu).
%
% The earlier version drew one line per method through the points sorted by
% CPU time. That is misleading twice over: sorting by cost is not a direction
% the sweep moves along, so the line zigzags and suggests a trend that does not
% exist, and nothing on the figure says which (phi, n_cc) a point is.
%
% Here a line joins the points that differ by ONE swept quantity, so following
% it means something:
%   with interface reduction   one line per (method, phi), n_cc increasing
%   without it                 one line per method, phi increasing
%
% Each point carries its n_cc as a label and its marker grows with it, each
% line is tagged with its phi at the expensive end, and the data tips hold the
% full identification: in the .fig, clicking a point says exactly which run it
% is.
%
% The vertical range follows the DATA. Forcing the integration floor into view
% costs nine empty decades when the errors sit between 1 and 300 per cent,
% which is most of the figure spent on nothing; when the floor falls outside
% the data it is stated in the subtitle instead of drawn.

    if nargin < 4 || isempty(ttl), ttl = 'Accuracy-cost trade-off (online)'; end

    methods_u = unique(T.method, 'stable');
    cols = lines(max(numel(methods_u), 3));
    mkr  = {'o','s','^','d','v','>','<','p'};

    fig = figure('Name','Accuracy-cost Pareto','Color','w','Position',[100 100 980 620]);
    hold on; grid on; box on
    set(gca, 'XScale','log', 'YScale','log');

    h_leg = gobjects(1, numel(methods_u));

    for m = 1:numel(methods_u)
        is_m = strcmp(T.method, methods_u{m});
        c    = cols(m, :);
        mk   = mkr{mod(m-1, numel(mkr)) + 1};

        % --- reduced interface: one line per phi, walking along n_cc ---
        phis = unique(T.phi(is_m & T.n_cc > 0));
        for p = 1:numel(phis)
            s = is_m & T.n_cc > 0 & T.phi == phis(p);
            if ~any(s), continue; end
            [nc, ord] = sort(T.n_cc(s));
            cpu = T.cpu(s);  cpu = cpu(ord);
            gre = T.gre_ref(s); gre = gre(ord);
            h = plot(cpu, gre, '-', 'Color', c, 'LineWidth', 1.4, ...
                'HandleVisibility', 'off');
            add_tips(h, methods_u{m}, phis(p), nc);
            markers(cpu, gre, nc, c, mk);
            labels(cpu, gre, nc, [0.3 0.3 0.3], 7, 'normal');
            % phi once per line, at its cheap end: the expensive end is where the
            % n_cc = n_bnd point sits on top of the un-reduced run of the same
            % method, so a tag there overprints the other one.
            [~, ie] = min(cpu);
            tag(cpu(ie), gre(ie), ['\phi=' num2str(phis(p))], c);
        end

        % --- full interface (no reduction): one line, walking along phi ---
        s = is_m & T.n_cc == 0;
        if any(s)
            [ph, ord] = sort(T.phi(s));
            cpu = T.cpu(s);  cpu = cpu(ord);
            gre = T.gre_ref(s); gre = gre(ord);
            h = plot(cpu, gre, '-', 'Color', c, 'LineWidth', 2.4, ...
                'HandleVisibility', 'off');
            add_tips(h, methods_u{m}, ph, zeros(size(ph)));
            plot(cpu, gre, mk, 'Color', c, 'MarkerFaceColor', c, ...
                'MarkerSize', 9, 'HandleVisibility', 'off');
            % phi is tagged only if this method has no reduced line to carry it:
            % at n_cc = n_bnd the reduced run IS this run, the two markers sit on
            % top of each other, and tagging both just overprints the label.
            if isempty(phis)
                for i = 1:numel(cpu)
                    tag(cpu(i), gre(i), ['\phi=' num2str(ph(i))], c);
                end
            end
        end

        h_leg(m) = plot(NaN, NaN, [mk '-'], 'Color', c, 'MarkerFaceColor', c, ...
            'LineWidth', 1.6, 'MarkerSize', 7, 'DisplayName', methods_u{m});
    end

    % --- vertical range on the data, floor drawn only if it lands inside ---
    g = T.gre_ref(isfinite(T.gre_ref) & T.gre_ref > 0);
    lo = 10^floor(log10(min(g))); hi = 10^ceil(log10(max(g)));
    sub = 'line = one sweep in n_{cc} at fixed \phi   |   marker size and label = n_{cc}   |   filled = full interface';
    if gre_floor_pct >= lo
        yline(gre_floor_pct, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Integration floor');
        legend([h_leg, findobj(gca,'Type','ConstantLine')'], 'Location','southwest');
    else
        sub = [sub newline '\rm\fontsize{8}integration floor at ' ...
               num2str(gre_floor_pct, '%.1e') ' %, below the range shown'];
        legend(h_leg, 'Location','southwest');
    end
    ylim([lo hi]);

    xlabel('Online CPU time [s]');
    ylabel(gre_desc);
    title({ttl, ['\rm\fontsize{9}' sub]});
    % The interactive toolbar is exported into the PNG otherwise.
    ax_h = gca;
    try, ax_h.Toolbar.Visible = 'off'; catch; end
end

% =====================================================================
function markers(cpu, gre, nc, c, mk)
% Marker area grows with n_cc, so the eye reads the interface size directly.
    sz = 5 + 7 * (log2(max(nc, 1)) / log2(256));
    for i = 1:numel(cpu)
        plot(cpu(i), gre(i), mk, 'Color', c, 'MarkerFaceColor', 'w', ...
            'MarkerSize', sz(i), 'LineWidth', 1.3, 'HandleVisibility', 'off');
    end
end

function labels(cpu, gre, nc, col, fs, wt)
% Offset ABOVE the marker, not centred on it: a white box centred on the point
% hides the very marker whose size is meant to encode n_cc. Drawn on top of
% everything by being called after all the lines of a series are in place.
    for i = 1:numel(cpu)
        text(cpu(i), gre(i)*1.14, num2str(nc(i)), 'FontSize', fs, 'Color', col, ...
            'FontWeight', wt, 'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom', 'Margin', 0.5, ...
            'BackgroundColor','w', 'Clipping','on');
    end
end

function tag(x, y, str, c)
    text(x, y, ['  ' str], 'Color', c, 'FontSize', 8, 'FontWeight','bold', ...
        'VerticalAlignment','middle', 'Margin', 0.5, ...
        'BackgroundColor','w', 'Clipping','on');
end

function add_tips(h, method, phi, ncc)
% Full identification on click. Version-sensitive, hence the guard: without it
% the figure is still correct, only less interrogable.
    try
        n = numel(h.XData);
        if isscalar(phi), phi = repmat(phi, n, 1); end
        h.DataTipTemplate.DataTipRows = [ ...
            dataTipTextRow('method', repmat({method}, n, 1)), ...
            dataTipTextRow('phi',    phi(:)), ...
            dataTipTextRow('n_cc',   ncc(:)), ...
            dataTipTextRow('cpu [s]',  h.XData(:)), ...
            dataTipTextRow('GRE [%]',  h.YData(:))];
    catch
    end
end
