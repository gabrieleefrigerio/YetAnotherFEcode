function animate_contact_3d(Struct, Interfaces, labels, y_contact, t, varargin)
%ANIMATE_CONTACT_3D The stoppers moving against their walls, in 3D.
%
%   ANIMATE_CONTACT_3D(Struct, Interfaces, labels, y_contact, t)
%   ANIMATE_CONTACT_3D(..., 'Amplify', 20, 'GifFile', 'contact.gif')
%
% Draws the mesh once as a static grey skin for context, then animates the
% contact nodes with their displacement amplified, lighting up the ones that
% are touching. Each wall is drawn as a translucent plane at its own gap,
% amplified by the same factor so that node and wall stay consistent.
%
% Only the contact nodes move. That is not a simplification for its own sake:
% the result files store the response at those nodes alone, 232 of 7476 here.
% To animate the whole deformed mesh the run has to save the full field, which
% run_fom does when cfg.save_full_field is set; this function then takes it
% through the 'Displacement' option.
%
% OPTIONS
%   'Amplify'      displacement scale factor (default 20). Displacements are
%                  microns on a structure hundreds of microns wide, so without
%                  amplification nothing is visible.
%   'Stride'       use one frame every N (default 1)
%   'GifFile'      path to write an animated GIF; omitted, it just plays
%   'DelayTime'    seconds per GIF frame (default 0.08)
%   'View'         [az el] (default [-35 20])
%   'ShowMesh'     draw the static mesh skin (default true)
%   'Displacement' [n_dofs x n_time] full field, to deform the whole mesh
%
% See also CONTACT_ACTIVITY, PLOT_CONTACT_NODES, EXTRACT_CONTACT_RESPONSE.

    p = inputParser;
    addParameter(p, 'Amplify', 20);
    addParameter(p, 'Stride', 1);
    addParameter(p, 'GifFile', '');
    addParameter(p, 'DelayTime', 0.08);
    addParameter(p, 'View', [-35 20]);
    addParameter(p, 'ShowMesh', false);
    addParameter(p, 'Displacement', []);
    parse(p, varargin{:});
    opt = p.Results;

    if ischar(labels) || isstring(labels), labels = cellstr(labels); end
    t = t(:)';
    frames = 1 : opt.Stride : numel(t);
    amp = opt.Amplify;

    % ---- full field: keep the drawn frames, on all the mesh DOFs -------
    % run_fom stores the response on the FREE DOFs, so a constrained model
    % gives fewer rows than the mesh has, and deforming the mesh with it
    % straight away would silently mismatch. Expanding costs a copy, so only
    % the frames actually drawn are kept.
    if ~isempty(opt.Displacement)
        n_full = numel(Struct.nodes);
        opt.Displacement = opt.Displacement(:, frames);
        if size(opt.Displacement, 1) ~= n_full
            opt.Displacement = Struct.AssemblyObj.unconstrain_vector(opt.Displacement);
        end
        if size(opt.Displacement, 1) ~= n_full
            error('ANIM3D:BadField', ...
                ['The displacement field has %d rows: the mesh has %d nodes ' ...
                 'x %d DOFs = %d.'], size(opt.Displacement, 1), ...
                size(Struct.nodes, 1), Struct.MeshObj.nDOFPerNode, n_full);
        end
    end

    % ---- per-interface geometry, gaps and penetration -----------------
    n_if = numel(labels);
    X0 = cell(1,n_if); U = cell(1,n_if); PEN = cell(1,n_if);
    nrm = cell(1,n_if); gap = cell(1,n_if);
    for i = 1:n_if
        I  = Interfaces.(labels{i});
        yc = y_contact.(labels{i});
        X0{i} = Struct.nodes(I.nodes(:), :);

        if isfield(I,'gap_nodes') && ~isempty(I.gap_nodes)
            gap{i} = I.gap_nodes(:);
        else
            gap{i} = abs(I.gap) * ones(size(X0{i},1),1);
        end
        if isfield(I,'normal') && ~isempty(I.normal)
            nv = I.normal(:)'; nv = nv / norm(nv);
        else
            nv = [0 0 1];
        end
        nrm{i} = nv;

        % Full nodal displacement if stored, otherwise the normal component
        % put back along the normal. Either way the contact motion is right.
        if isfield(yc, 'U') && ~isempty(yc.U)
            U{i} = yc.U;
        else
            un = yc.normal;
            U{i} = reshape(reshape(un, [], 1) * nv, size(un,1), 3, size(un,2));
        end
        PEN{i} = yc.normal - gap{i};
    end

    % ---- figure and static context ------------------------------------
    fig = figure('Name','Contact animation','Color','w','Position',[80 60 1000 760]);
    ax = axes(fig); hold(ax,'on');
    if opt.ShowMesh
        try
            PlotMesh(Struct.nodes, Struct.plot_connectivity(), 0);
            set(findobj(ax,'Type','patch'), 'FaceAlpha', 0.06, 'EdgeAlpha', 0.12);
            set(findobj(ax,'-property','Marker'), 'Marker', 'none');
        catch
            warning('ANIM3D:NoMesh', 'Mesh skin could not be drawn; continuing without it.');
        end
    end

    % Deformed mesh, when the full field is available
    hMesh = [];
    if ~isempty(opt.Displacement)
        conn = Struct.plot_connectivity();
        hMesh = patch('Faces', conn(:,1:min(3,size(conn,2))), ...
                      'Vertices', Struct.nodes, 'FaceColor', [.7 .78 .88], ...
                      'FaceAlpha', 0.25, 'EdgeColor', [.5 .55 .62], 'EdgeAlpha', .25);
    end

    % ---- walls, one translucent plane per interface -------------------
    for i = 1:n_if
        drawWall(ax, X0{i}, nrm{i}, amp*mean(gap{i}));
    end

    % ---- moving node clouds -------------------------------------------
    hFree = gobjects(1,n_if); hHit = gobjects(1,n_if);
    for i = 1:n_if
        hFree(i) = plot3(ax, nan, nan, nan, 'o', 'MarkerSize', 4, ...
            'MarkerFaceColor', [.45 .55 .70], 'MarkerEdgeColor','none', 'LineStyle','none');
        hHit(i)  = plot3(ax, nan, nan, nan, 'o', 'MarkerSize', 9, ...
            'MarkerFaceColor', [.85 .15 .12], 'MarkerEdgeColor','k', 'LineStyle','none');
    end

    axis(ax,'equal'); grid(ax,'on'); box(ax,'on'); view(ax, opt.View);
    xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'z');
    pad = 0.06 * max(range(Struct.nodes,1));
    xlim(ax, [min(Struct.nodes(:,1))-pad max(Struct.nodes(:,1))+pad]);
    ylim(ax, [min(Struct.nodes(:,2))-pad max(Struct.nodes(:,2))+pad]);
    zlim(ax, [min(Struct.nodes(:,3))-pad-amp*3 max(Struct.nodes(:,3))+pad+amp*3]);
    ttl = title(ax, '', 'FontWeight','bold');

    % ---- frames --------------------------------------------------------
    for k = 1:numel(frames)
        j = frames(k);
        n_act = 0;
        for i = 1:n_if
            xd = X0{i} + amp * U{i}(:,:,j);
            a  = PEN{i}(:,j) > 0;
            n_act = n_act + nnz(a);
            set(hFree(i), 'XData', xd(~a,1), 'YData', xd(~a,2), 'ZData', xd(~a,3));
            set(hHit(i),  'XData', xd( a,1), 'YData', xd( a,2), 'ZData', xd( a,3));
        end
        if ~isempty(hMesh)
            % Column k, not j: the field was subsampled to the drawn frames.
            d = reshape(opt.Displacement(:,k), size(Struct.nodes,2), []).';
            set(hMesh, 'Vertices', Struct.nodes + amp*d);
        end
        set(ttl, 'String', sprintf('t = %7.3f \\mus   -   nodes in contact: %3d / %d   (x%g)', ...
            1e6*t(j), n_act, sum(cellfun(@(x) size(x,1), X0)), amp));
        drawnow limitrate

        if ~isempty(opt.GifFile)
            fr = getframe(fig);
            [A, map] = rgb2ind(frame2im(fr), 256);
            if k == 1
                imwrite(A, map, opt.GifFile, 'gif', 'LoopCount', Inf, ...
                    'DelayTime', opt.DelayTime);
            else
                imwrite(A, map, opt.GifFile, 'gif', 'WriteMode', 'append', ...
                    'DelayTime', opt.DelayTime);
            end
        end
    end

    if ~isempty(opt.GifFile)
        fprintf('Animation written to %s (%d frames)\n', opt.GifFile, numel(frames));
    end
end

% =====================================================================
function drawWall(ax, X0, nv, offset)
%DRAWWALL A translucent rectangle at the wall position of one interface.
% Spanned by the two in-plane directions of the patch, offset along the
% normal by the (amplified) gap.
    c = mean(X0, 1);
    d = X0 - c;
    [~, ~, V] = svd(d, 'econ');       % principal directions of the patch
    e1 = V(:,1)'; e2 = V(:,2)';
    s1 = 1.15 * max(abs(d*e1'));  s2 = 1.15 * max(abs(d*e2'));
    if s1 == 0, s1 = 1; end
    if s2 == 0, s2 = 1; end
    p0 = c + offset*nv;
    P  = [p0 - s1*e1 - s2*e2
          p0 + s1*e1 - s2*e2
          p0 + s1*e1 + s2*e2
          p0 - s1*e1 + s2*e2];
    patch(ax, 'Vertices', P, 'Faces', [1 2 3 4], ...
        'FaceColor', [.25 .25 .28], 'FaceAlpha', 0.28, 'EdgeColor', [.3 .3 .3]);
end
