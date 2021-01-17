% Temporal smoothing (TS) trials for the vector sensor
% Summer 2020
% Will Howard, {wwhoward}@vt.edu
% wireless @ VT
clc;clear;close all

addpath('configs')
config_name = 'associate_evaluate';

% Multipath Scenarios
par.K = 3;                  % Keep K <= 4: This is the total incident paths on the sensor. This is overwritten if par.forcePath is not empty. 
par.forceMulti = false;     % If true, enables path states of [1], [1 2], etc
par.forcePath = [1 1 1];    % set par.forcePath = [] for random path configuration. Overrides par.K, so be careful using this
par.powerDecay = 1;         % Fraction of power lost on each multipath
par.decayType = 'none';      % {'exp', 'rnd', 'none'}. Set to 'none' to ignore. 'exp' generates powers in an exponential sequence: p1 = powerDecay^0, p2=powerDecay^1, ... while 'rand' generates random losses in [par.powerDecay, 1] for path=2,.. i.e. p1 = 1, p2=powerDecay*rand, ..

% Signals
par.signal_length = 2^15;   % Long enough to fit blocks*(snapshot+interboock)
par.mod = 'QPSK';           % Modulation scheme in {'QPSK'}, more added if needed
par.interblock = [0 20];     % Delay between TS blocks. Samples randomly in this range [a b]. 
par.pathdelay = [0 0];      % Delay between different recieved paths. If this is a range [a,b], samples randomly. 
par.blocks = 5;             % Number of blocks for TS averaging. (par.K+1 is a good default)
par.blockSweep = 1:10;      % Specific to runtype 'blocksweep': varies number of blocks each run
par.snapshot = 2^10;         % Snapshot window for signals
par.SNR = 30;               % Static SNR for runtype 'single', overwritten for sweeps
par.snrSweep = -20:2:40;       % Variable SNR for runtypes other than 'single'
par.type = '2d';            % Problem dimensionality in {'1d', '2d'}. Effects MUSIC-type searches as well as DoA ground truth assignments. 
par.minSep = pi/10;         % Minimum separation between azimuth AoA (in radians)
par.aziRange = [pi/20, 39*pi/20]; % Specify range that signal azimuths must be in
par.eleRange = [pi/20, 19*pi/20]; % Same for signal elevations
par.polType = 'rnd';        % {'rnd_lin', 'rnd'}
par.forcePol = [];          % Force all signals to this polarization
par.manifold = 'sim';       % {'sim', '8in'} Which manifold to use for data simulation. Sim is normal, 8in is ASI

% Estimation & Statistics
par.aziEstRange = [0,2*pi]; % Bounds for estimating (should only effect music-type estimators)
par.eleEstRange = [0,pi];   % Bounds for estimating (should only effect music-type estimators)
par.res = 2^8 / pi;         % Resolution for MUSIC-type estimators
par.recombination = 'rnd';  % {'rnd','max'}
par.scrub = true;           % Additional variable for ml generation - removes RMSE outliers - NOTE: If true, then X (output from ml_gen) takes on the "clean" values
par.sampling = 1;           % WORK IN PROGRESS DON'T CHANGE - Calculate statistics for this percent of trials (includes DoA estimation & statistics) - lower value speeds up processing. Only implemented in "ml_gen"                   

% Simulation
par.Trials = 100;          % During a sweep, how many trials for each parameter
par.runtype = 'realtime';     % {'single','snr_sweep','block_sweep','ml_gen'} - names are fairly self-explainitory
par.estimator = ["music"];   % Vector of strings, can be {'music', 'dr-music','fusion','ml'}
par.smooth = ["ts", "nts"];  %Whether to do temporal smoothing (ts) or not. Can also do both (ts,nts)
par.accelerate = 1;         % See readme

% Fusion
par.fusion_interval = 1;   % How many degrees per fusion interval (in degrees for human intepretability)
par.fusion_res = 10;       % How many points to test in each interval

% Paths to model directories (dummies for now)
par.interval_model_path = 0;
par.order_model_path = 0;
par.doa_model_path = 0; 

% Save configuration
par.saveFlag = 1;           % Should anything be saved? NOTE: Overridden to 0 for runtype 'single'
par.saveLight = 0;          % Should everything be saved or just statistics? 
par.saveName = "for_jet";  % Appends this to the name for ease of identification (edited later to be unique)
par.savePlots = 1;          % Save plots as fig, eps, png


%% 
%---------------------\\-------------------------
% Ensure that the specified parameters are compatible
%---------------------\\------------------------- 
if ~isempty(par.forcePath) % Check that the number of paths matches the forced path configuration, if it exists. 
    par.K = length(par.forcePath);
end
if any(strcmp(par.estimator, "music")) % Check that the polarization is defaulted if MUSIC is being used (might change later, if MUSIC polarization estimation is added)
    warning("For now, MUSIC only searches over azimuth and elevation. Setting polarization to a default...")
    par.forcePol = [pi/4, 0];
end
if strcmp(par.manifold, '8in') && any(strcmp(par.estimator, 'drmusic')) % Make sure to use MUSIC for the 8in manifold, since we have the 4 parameter manifold but not the 2 parameter
    warning('ASI manifold not compatible with DR-MUSIC estimator. Switching to MUSIC. ')
    par.estimator = "music";
end
if strcmp(par.manifold, '8in') % Fix signal parameters to be within the bounds of the 8in manifold
    if par.eleEstRange(1) < pi/4; par.eleEstRange(1) = pi/4; end
    if par.eleEstRange(2) > 10*pi/18; par.eleEstRange(2) = 10*pi/18; end
    if par.eleRange(1) < pi/4; par.eleRange(1) = pi/4; end
    if par.eleRange(2) > 10*pi/18; par.eleRange(2) = 10*pi/18; end
end
if mod(180,par.fusion_interval)~=0 % Calculate an appropriate fusion interval, that is close to the requested value. Also convert to radians. 
    par.fusion_interval = 180/floor(180/par.fusion_interval);
end
par.fusion_interval = par.fusion_interval/180 * pi; % convert to radians for compatability and speed

par.fusion_res = par.fusion_res/par.fusion_interval;

if exist('par2', 'var') % If pre-saved parameters were loaded, inject them here. 
sameFields = string(intersect(fieldnames(par2), fieldnames(par)));
    for f = 1:length(sameFields)
        par.(sameFields(f)) = par2.(sameFields(f));    
    end
end

%%
%---------------------\\-------------------------
% Calculate the appropriate manifolds
%---------------------\\-------------------------

%if any(strcmp(par.estimator,'drmusic')) || any(strcmp(par.estimator,'music')) % Get the necessary manifolds (single for just MUSIC, interval-based for fusion
[par.AOA_function, par.AOA_function_H, par.phi, par.tht, par.DOA_function, par.DOA_function_H] = get_manifold(par);


save_name="configs/"+config_name;
if isfile(save_name+".mat")
    display('Warning: Previous file overwritten')
end
save(save_name, "par");
display("Saved as " + config_name);