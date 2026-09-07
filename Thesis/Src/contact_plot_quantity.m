function [v, wall, ylab] = contact_plot_quantity(yc, nrm, gap, mode, d_ref)
%CONTACT_PLOT_QUANTITY The quantity drawn in a per-stopper panel, and its wall.
%
%   [v, wall, ylab] = CONTACT_PLOT_QUANTITY(yc, nrm, gap, mode, d_ref)
%
% INPUT
%   yc     one interface of a y_contact struct (fields normal, U, ...)
%   nrm    wall normal of that interface
%   gap    positive gap of that interface
%   mode   'physical' | 'clearance' | 'normal'
%   d_ref  reference direction shared by every panel, normally the shock one
%
% OUTPUT
%   v      [n_nodes x n_t] the quantity to draw
%   wall   scalar ordinate of the wall in that quantity
%   ylab   axis label
%
% 'physical' projects the nodal displacement on ONE reference direction,
% shared by every panel, and puts the wall at its signed position along it:
% a wall whose normal agrees with the reference sits at +gap, one opposing it
% at -gap. That is what makes a mass driven along +Z rise in all eight panels
% while only the upper walls are approached. Both alternatives were tried and
% both mislead: 'normal' sends the bottom panels downwards and 'clearance'
% sends the top ones downwards, in each case for a sign convention rather than
% for anything the structure does.
%
% 'physical' needs the normal to be parallel to the reference direction, which
% holds for stoppers facing along the drive axis. An oblique wall has no single
% line in that coordinate, so it falls back to the normal component and says so
% through its label.

    un = yc.normal;

    switch lower(mode)
        case 'clearance'
            v = gap - un;  wall = 0;  ylab = 'clearance';

        case 'physical'
            algn = dot(nrm(:)', d_ref(:)');
            if abs(abs(algn) - 1) < 1e-6 && isfield(yc, 'U') && ~isempty(yc.U)
                U = yc.U;                                  % [n_nodes x n_dim x n_t]
                v = squeeze(sum(U .* reshape(d_ref, 1, [], 1), 2));
                if size(U,1) == 1, v = reshape(v, 1, []); end
                wall = sign(algn) * gap;
                ylab = sprintf('u along [%g %g %g]', d_ref);
            else
                v = un;  wall = gap;  ylab = 'u_n (oblique)';
            end

        otherwise
            v = un;  wall = gap;  ylab = 'u_n';
    end
end
