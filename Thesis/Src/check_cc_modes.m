function check_cc_modes(cfg, contact)
%CHECK_CC_MODES Validate the requested CC mode counts against the interface.
%
%   CHECK_CC_MODES(cfg, contact)
%
% Catches a misconfigured sweep before any time is spent integrating: asking
% for more CC modes than there are interface DOFs is meaningless, and asking
% for exactly n_bnd is the useful control point, since it applies the change
% of basis without truncating anything and must therefore reproduce the
% untruncated baseline.

    if ~cfg.interface_reduction.enabled, return; end

    bad = cfg.array_ccModes(cfg.array_ccModes < 1 | ...
                            cfg.array_ccModes > contact.n_bnd);
    if ~isempty(bad)
        error('check_cc_modes:OutOfRange', ...
            ['cfg.array_ccModes contains values out of range: %s. ' ...
             'They must lie between 1 and n_bnd = %d.'], ...
            mat2str(bad), contact.n_bnd);
    end

    fprintf('Interface reduction: %s / %s basis | CC modes %s (n_bnd = %d)\n', ...
        cfg.interface_reduction.mode, cfg.interface_reduction.basis, ...
        mat2str(cfg.array_ccModes), contact.n_bnd);

    if ~any(cfg.array_ccModes == contact.n_bnd)
        fprintf(['  [note] n_cc = %d is not in the list. It is the control ' ...
                 'point that must\n         reproduce the baseline, and it ' ...
                 'is worth including.\n'], contact.n_bnd);
    end
end
