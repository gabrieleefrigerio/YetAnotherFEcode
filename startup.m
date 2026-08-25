clc

[filepath, name, ext] = fileparts(mfilename('fullpath'));

addpath(genpath(fullfile(filepath, 'src')));
addpath(genpath(fullfile(filepath, 'external')));
addpath(genpath(fullfile(filepath, 'examples', 'Meshes')));

% Thesis library shared by the test cases. Plain addpath, not genpath: the
% folder is flat, and genpath would also pick up any results/ subfolder.
% The mains resolve it themselves from mfilename, but that only works when the
% whole file is run - not when a selection is evaluated - and the
% post-processing scripts do not resolve it at all, so it belongs here.
addpath(fullfile(filepath, 'Thesis', 'Src'));

disp('              _____ _____     ')
disp('  _   _  __ _|  ___| ____|___ ')
disp(' | | | |/ _` | |_  |  _| / __|')
disp(' | |_| | (_| |  _| | |__| (__ ')
disp('  \__, |\__,_|_|   |_____\___|')
disp('  |___/       YetAnotherFEcode')
fprintf('\n\n')