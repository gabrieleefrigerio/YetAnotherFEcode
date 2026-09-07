function fig = plot_shapes(Struct, V, varargin)
%PLOT_SHAPES Draw shape vectors on the mesh, 2D or 3D.
%
%   fig = PLOT_SHAPES(Struct, V)
%   fig = PLOT_SHAPES(Struct, V, 'Labels', {...}, 'Scale', s, 'Cols', n)
%
% V holds one shape per COLUMN, on ALL the degrees of freedom of the model
% (the unconstrained numbering, i.e. n_nodes*nDOFPerNode rows). That is the
% numbering of Struct.mode_shapes and of the P property every ROM class
% exposes, so the same call draws any of them:
%
%   eigenmodes            V = Struct.mode_shapes
%   any ROM basis         V = rom.P                (CB, Rubin, MC, MCB, MN)
%   constraint modes      V = rom.P(:, 1:rom.n_bnd)      for a CB ROM
%   interface (CC) modes  V = rom.P(:, 1:rom.n_bnd) * Phi_CC
%
% The last one is the useful way to look at an interface reduction: Phi_CC
% lives on the contact DOFs alone, and multiplying by the constraint modes
% extends it to the whole body, which is the field that CC mode actually
% represents.
%
% OPTIONS
%   'Labels'  cellstr, one per column. Default: 'shape k'.
%   'Scale'   amplification. Default: each shape is first normalised to unit
%             largest nodal displacement, then drawn at 20% of the model size,
%             so the picture is independent of the units and of how the vector
%             happened to be scaled.
%   'Cols'    tiles per row. Default: ceil(sqrt(k)).
%   'Title'   figure title.
%
% See also FESTRUCTURE/PLOT_MODE, PLOT_CONTACT_NODES.

    p = inputParser;
    addParameter(p, 'Labels', {});
    addParameter(p, 'Scale',  []);
    addParameter(p, 'Cols',   []);
    addParameter(p, 'Title',  'Shapes');
    addParameter(p, 'MaxStrain', 0.3);
    parse(p, varargin{:});
    opt = p.Results;

    nD  = Struct.MeshObj.nDOFPerNode;
    nDo = size(Struct.nodes, 1) * nD;
    if size(V, 1) ~= nDo
        error('PlotShapes:BadSize', ...
            ['V must have one row per DOF of the FULL model (%d = %d nodes x %d), ' ...
             'got %d.\nA basis expressed on the FREE DOFs has to be expanded ' ...
             'first: use rom.P instead of rom.Pc, or\n' ...
             '  V = Struct.AssemblyObj.unconstrain_vector(V_free).'], ...
            nDo, size(Struct.nodes,1), nD, size(V,1));
    end

    k = size(V, 2);
    if isempty(opt.Cols), opt.Cols = ceil(sqrt(k)); end
    if isempty(opt.Labels)
        opt.Labels = arrayfun(@(i) sprintf('shape %d', i), 1:k, 'Uni', 0);
    end
    if numel(opt.Labels) ~= k
        error('PlotShapes:BadLabels', '%d labels for %d shapes.', numel(opt.Labels), k);
    end

    span = max(max(Struct.nodes, [], 1) - min(Struct.nodes, [], 1));
    if isempty(opt.Scale), opt.Scale = 0.2 * span; end

    rows = ceil(k / opt.Cols);
    fig = figure('Name', opt.Title, 'Color', 'w', ...
                 'Position', [60 60 min(360*opt.Cols, 1600) min(320*rows, 950)]);
    tl = tiledlayout(fig, rows, opt.Cols, 'TileSpacing', 'compact', 'Padding', 'compact');
    conn = Struct.plot_connectivity();
    ax = gobjects(1, k);
    failed = false(1, k);

    % Node pairs from every element's own node list, computed ONCE (it only
    % depends on the mesh). This is the general safety check: a linear "plot
    % at scale s" moves node p to node(p,:) + s*v(p,:), and that stops being a
    % believable small picture once s*|v(p,:)-v(q,:)| becomes comparable to
    % the UNDEFORMED distance between two nodes p, q of the same element - the
    % element folds onto itself and the renderer's triangulation breaks.
    %
    % This subsumes a pure rigid rotation as a special case (a rotating body
    % keeps every true distance exactly fixed, so linearising a large angle is
    % what manufactures the strain, concentrated far from the axis) but it ALSO
    % catches shapes that are not rotations at all: Rubin's own basis columns
    % are ATTACHMENT modes, i.e. flexibility under a unit FORCE, not a unit
    % DISPLACEMENT, so their local scaling can be wildly uneven even with no
    % rotation involved - measured at strain ~20 (2000%) on every one of six
    % modes, which is exactly the failure a rotation-only check cannot see.
    E  = Struct.elements;
    np = size(E, 2);
    pr = zeros(0, 2);
    for a = 1:np-1
        for b = a+1:np
            pr = [pr; E(:,a), E(:,b)]; %#ok<AGROW>
        end
    end
    pr  = unique(sort(pr, 2), 'rows');
    pr  = pr(pr(:,1) ~= pr(:,2), :);
    dn  = Struct.nodes(pr(:,1), :) - Struct.nodes(pr(:,2), :);
    Ln  = sqrt(sum(dn.^2, 2));
    keep = Ln > 1e-12 * span;
    pr = pr(keep, :); Ln = Ln(keep);

    for i = 1:k
        ax(i) = nexttile(tl);
        v = reshape(V(:, i), nD, []).';

        % Normalise on the largest NODAL displacement, not on the largest DOF:
        % a vector scaled by its largest single component draws a lopsided
        % picture whenever one direction dominates.
        amp = max(sqrt(sum(v.^2, 2)));
        if amp > 0, v = v / amp; end

        % A single linear scale cannot be safe for every shape at once (see
        % above), so it is shrunk for THIS TILE ONLY until the worst-case
        % relative displacement between two nodes of the same element stays
        % under a threshold well short of self-intersection.
        s_eff = opt.Scale;
        dv = v(pr(:,1), :) - v(pr(:,2), :);
        strain_max = max(sqrt(sum(dv.^2, 2)) ./ Ln);
        if s_eff * strain_max > opt.MaxStrain
            s_eff = opt.MaxStrain / max(strain_max, eps);
        end

        try
            PlotFieldonDeformedMesh(Struct.nodes, conn, v, 'factor', s_eff);
            set(findobj(ax(i), '-property', 'Marker'), 'Marker', 'none');
        catch ME
            failed(i) = true;
            cla(ax(i));
            text(0.5, 0.5, sprintf('render failed\n(%s)', ME.identifier), ...
                'Units', 'normalized', 'HorizontalAlignment', 'center', ...
                'Color', [0.6 0 0], 'FontSize', 8);
        end

        ttl = opt.Labels{i};
        if s_eff < opt.Scale * 0.999
            ttl = sprintf('%s  (%.0f%% scale)', ttl, 100*s_eff/opt.Scale);
        end
        title(ttl, 'Interpreter', 'none', 'FontSize', 9);
        axis equal off;
        if Struct.n_dim == 3, view(3); end
    end

    if any(failed)
        warning('PlotShapes:RenderFailed', ...
            'Could not render shape(s) %s even at the reduced scale.', ...
            mat2str(find(failed)));
    end

    % One colour bar for the whole figure. The plotting helper adds its own to
    % every tile, and six of them take more room than the shapes while saying
    % nothing: each shape is normalised to unit displacement, so the scale is
    % the same everywhere and only the RELATIVE amplitude carries meaning.
    % Wrapped in try/catch: a tile that failed to render can leave the figure's
    % graphics scene in a state where touching another object's properties
    % (colorbar.Label included) throws, and that must not take the whole plot
    % down when the shapes themselves are already on screen.
    try
        delete(findobj(fig, 'Type', 'colorbar'));
        colormap(fig, jet);
        for i = 1:k
            if isgraphics(ax(i)), clim(ax(i), [0 1]); end
        end
        cb = colorbar(ax(end));
        cb.Layout.Tile = 'east';
        cb.Label.String = 'relative nodal displacement';
    catch
    end

    % Two lines, the second empty: with 'axis equal' the 3D tiles do not all
    % have the same height, so the tallest one pushes its title up into the
    % layout title. The blank line is the reliable way to keep them apart.
    title(tl, {opt.Title, ''}, 'FontWeight', 'bold');
end

