% Temporal smoothing (TS) trials for the vector sensor
% Summer 2020
% Will Howard, {wwhoward}@vt.edu
% wireless @ VT
clc;clear;close all

addpath('configs')
config_name = '2d_ml_linpol';

% Multipath Scenarios
par.K = 2;                  % Keep K <= 4: This is the total incident paths on the sensor. This is overwritten if par.forcePath is not empty. 
par.forceMulti = false;     % If true, enables path states of [1], [1 2], etc
par.forcePath = [1 1];    % set par.forcePath = [] for random path configuration. Overrides par.K, so be careful using this
par.type = '2d';            % Problem dimensionality in {'1d', '2d'}. Effects MUSIC-type searches as well as DoA ground truth assignments. 
par.minSep = pi/10;         % Minimum separation between azimuth AoA (in radians)
par.aziRange = [pi/20, 39*pi/20]; % Specify range that azimuths must be in
par.eleRange = [pi/20, 19*pi/20]; % Same for elevations
par.polType = 'rnd_lin';        % {'rnd_lin', 'rnd'}
par.forcePol = [];          % Not implemented yet, don't use

% Signal settings
par.signal_length = 2^13;   % Long enough to fit blocks*(snapshot+interboock)
par.mod = 'QPSK';           % Modulation scheme in {'QPSK'}, more added if needed
par.interblock = [5 5];     % Delay between TS blocks. Samples randomly in this range [a b]. 
par.pathdelay = [0 0];      % Delay between different recieved paths. If this is a range [a,b], samples randomly. 
par.blocks = 4;             % Number of blocks for TS averaging. (par.K+1 is a good default)
par.blockSweep = 1:10;      % Specific to runtype 'blocksweep': varies number of blocks each run
par.snapshot = 2^10;         % Snapshot window for signals
par.SNR = 10;               % Static SNR for runtype 'single', overwritten for sweeps
par.snrSweep = 20:2:30;       % Variable SNR for runtypes other than 'single'

% Estimation & Statistics
par.res = 2^5 / pi;         % Resolution for MUSIC-type estimators
par.recombination = 'rnd';  % {'rnd','max'}
par.scrub = true;           % Additional variable for ml generation - removes RMSE outliers - NOTE: If true, then X (output from ml_gen) takes on the "clean" values
par.sampling = 1;           % WORK IN PROGRESS DON'T CHANGE - Calculate statistics for this percent of trials (includes DoA estimation & statistics) - lower value speeds up processing. Only implemented in "ml_gen"                   

% Simulation
par.Trials = 10000;          % During a sweep, how many trials for each parameter
par.runtype = 'ml_gen';     % {'single','snr_sweep','block_sweep','ml_gen'} - names are fairly self-explainitory
par.accelerate = 1;         % See readme

% Save configuration
par.saveFlag = 1;           % Should anything be saved? NOTE: Overridden to 0 for runtype 'single'
par.saveLight = 0;          % Should everything be saved or just statistics? 
par.saveName = "2path_highSNR_100k_1d";  % Appends this to the name for ease of identification (edited later to be unique)
par.savePlots = 1;          % Save plots as fig, eps, png

[par.DOA_function, par.DOA_function_inverse, par.phi, par.tht] = getManifold(par);


save_name="configs/"+config_name;
if isfile(save_name+".mat")
    display('Warning: Previous file overwritten')
end
save(save_name, "par");
display("Saved as " + config_name);