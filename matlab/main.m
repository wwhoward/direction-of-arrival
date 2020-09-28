% Temporal smoothing (TS) trials for the vector sensor
% Summer 2020
% Last updated: 9/22/2020
% Will Howard, {wwhoward}@vt.edu
% wireless @ VT

% Notation:
% A multipath state of [1 1] indicates one signal on two paths
% Similarly, [1 2] indicates two signals on one path each and [1 1 2 3] indicates one signal on two paths + two more signals on one path each

clc,clear,close all
addpath('data', 'Helper', 'functions', 'runtype', 'configs', 'Helper/mtimesx', 'Helper/Eig3Folder'); 
par2 = load('sept_2_hires_2d'); % Load saved parameters. Note that if the saved file doesn't have par.x, it will be replaced by whatever par.x is below. 
par2 = par2.par;


% Set parameters! 
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
par.powerDecay = 0.5;       % Set to 1 to ignore.  Fraction of power lost on each multipath (i.e. first path has powerDecay^0, 2nd powerDecay^1, ...

% Signal settings
par.signal_length = 2^13;   % Long enough to fit blocks*(snapshot+interboock)
par.mod = 'QPSK';           % Modulation scheme in {'QPSK'}, more added if needed
par.interblock = [5 5];     % Delay between TS blocks. Samples randomly in this range [a b]. 
par.pathdelay = [0 0];      % Delay between different recieved paths. If this is a range [a,b], samples randomly. 
par.blocks = 4;             % Number of blocks for TS averaging. (par.K+1 is a good default)
par.blockSweep = 1:10;      % Specific to runtype 'blocksweep': varies number of blocks each run
par.snapshot = 2^10;         % Snapshot window for signals
par.SNR = 10;               % Static SNR for runtype 'single', overwritten for sweeps
par.snrSweep = 28:2:30;       % Variable SNR for runtypes other than 'single'

% Estimation & Statistics
par.res = 2^6 / pi;         % Resolution for MUSIC-type estimators
par.recombination = 'rnd';  % {'rnd','max'}
par.scrub = true;           % Additional variable for ml generation - removes RMSE outliers - NOTE: If true, then X (output from ml_gen) takes on the "clean" values
par.sampling = 1;           % WORK IN PROGRESS DON'T CHANGE - Calculate statistics for this percent of trials (includes DoA estimation & statistics) - lower value speeds up processing. Only implemented in "ml_gen"                   

% Simulation
par.Trials = 100000;          % During a sweep, how many trials for each parameter
par.runtype = 'single';     % {'single','snr_sweep','block_sweep','ml_gen'} - names are fairly self-explainitory
par.accelerate = 1;         % See readme

% Save configuration
par.saveFlag = 1;           % Should anything be saved? NOTE: Overridden to 0 for runtype 'single'
par.saveLight = 0;          % Should everything be saved or just statistics? 
par.saveName = "2path_highSNR_100k_1d";  % Appends this to the name for ease of identification (edited later to be unique)
par.savePlots = 1;          % Save plots as fig, eps, png

[par.DOA_function, par.DOA_function_inverse, par.phi, par.tht] = getManifold(par);

if exist('par2', 'var')
sameFields = string(intersect(fieldnames(par2), fieldnames(par)));
for f = 1:length(sameFields)
    par.(sameFields(f)) = par2.(sameFields(f));    
end
end
%---------------------\\-------------------------

switch par.runtype
    case 'single'
        [paths, signal, est_ts, est_nts, stats_ts, stats_nts] = run_single(par);
        par.saveFlag=0;
    case 'ml_gen'
        [X, Y, RMSE, stats] = run_ml_gen(par);
    case 'snr_sweep'
        [paths, signal, est, stats] =  run_snr_sweep(par);
    case 'block_sweep'
        run_block_sweep(par);
end

if par.saveFlag
    i=0;
    savename = "data/"+date+"_"+par.saveName + "_" + string(i);
    while isfile(savename+".mat")
        i=i+1;
        savename = "data/"+date+"_"+par.saveName+"_"+string(i);
    end
    if par.saveLight
        save(savename, "stats")
    else
        save(savename)
    end
    
    clc
    display('Data saved as '+savename)
    
    if par.savePlots
        mkdir(savename)
        FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
        for f = 1:length(FigList)
            savefig(FigList(f), fullfile(savename, [get(FigList(f),'Name'), '.fig'])); %savename+string(get(FigList(f),'Name'))+'.fig') 
            saveas(FigList(f), fullfile(savename, [get(FigList(f),'Name'), '.png'])); %savename+string(get(FigList(f),'Name'))+'.png') 
        end
    end
end