function [rom, Pc, Mr, Kr, Cr] = build_rom(Struct, model, phi, contact, k_contact, rayleigh)
%BUILD_ROM Build one reduction basis and its reduced matrices.
%
%   [rom, Pc, Mr, Kr, Cr] = BUILD_ROM(Struct, model, phi, contact, k_contact, rayleigh)
%
%   model  'MT' | 'MC' | 'CB' | 'Rubin' | 'MCB' | 'MN'
%   phi    number of retained linear modes
%
% Pc always maps the reduced coordinates to the CONSTRAINED physical DOFs.
% MT and MC expose their basis on the global (unconstrained) DOFs, so it is
% brought back here; the CMS classes already expose a constrained Pc. Having
% one convention is what lets the caller write the contact operator as N*Pc
% for every method without a special case.
%
% See also ROMMC, ROMCB, ROMRUBIN, ROMMCB, ROMMN.

    switch model
        case 'MT'                                   % modal truncation
            rom = RomMC(Struct, phi, contact.dofs, k_contact, 0);
            rom.build();
        case 'MC'                                   % Milman-Chu
            rom = RomMC(Struct, phi, contact.dofs, k_contact, 1);
            rom.build();
        case 'CB'
            rom = RomCB(Struct, phi, contact.dofs);
            rom.build();
        case 'Rubin'
            rom = RomRubin(Struct, phi, contact.dofs);
            rom.build();
        case 'MCB'
            rom = RomMCB(Struct, phi, contact.dofs);
            rom.build(rayleigh);
            rom.check();
        case 'MN'
            rom = RomMN(Struct, phi, contact.dofs);
            rom.build(rayleigh);
            rom.check();
        otherwise
            error('build_rom:BadModel', 'Unknown reduction method ''%s''.', model);
    end

    if any(strcmp(model, {'MT', 'MC'}))
        n_dofs = size(Struct.AssemblyObj.constrain_matrix(Struct.K), 1);
        Pc = zeros(n_dofs, size(rom.P, 2));
        for ic = 1:size(Pc, 2)
            Pc(:, ic) = Struct.AssemblyObj.constrain_vector(rom.P(:, ic));
        end
    else
        Pc = rom.Pc;
    end

    [Mr, Kr, Cr] = rom.get_reduced_matrices();
end
