clear; close all; clc;

% 1. Carico tutte le configurazioni
cfg = config();

myMems = MicroCantilever(cfg);
myMems.buildModel();
[f_fom, Phi_fom] = myMems.runFullOrderModal(3);


myMems.plotMode(Phi_fom(:,1), sprintf('Modo 1 - Full Order: %.2f Hz', f_fom(1)));