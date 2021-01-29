% Temporal smoothing (TS) trials for the vector sensor
% Summer 2020
% Will Howard, {wwhoward}@vt.edu
% wireless @ VT

% Specify parameters (struct par) for DoA code
% Each parameter is explained below


clc;clear;close all

addpath('configs')
config_name = 'test';

%%
%---------------------\\-------------------------
% Major Parameters: 
par.K = 3;                  % Keep K <= 4: This is the total incident paths on the sensor. This is overwritten if par.forcePath is not empty. 
par.forcePath = [1 2 3];    % set par.forcePath = [] for random path configuration. Overrides par.K, so be careful using this
par.type = '1d';            % Problem dimensionality in {'1d', '2d'}. Effects MUSIC-type searches as well as DoA ground truth assignments. 
par.Trials = 100;           % During a sweep, how many trials for each parameter
par.runtype = 'single';  % {'single','snr_sweep','block_sweep','ml_gen'} - names are fairly self-explainitory
par.estimator = ["drmusic"];  % Vector of strings, can be {'music', 'drmusic','fusion','ml'}
par.manifold = 'sim';       % {'sim', '8in'} Which manifold to use for data simulation. Sim is normal, 8in is ASI


% Minor Parameters: (less likely to need modification)
%%
% Multipath Scenarios
par.forceMulti = false;     % If true, enables path states of [1], [1 2], etc
par.powerDecay = 1;         % Fraction of power lost on each multipath (1 = no decay)
par.decayType = 'none';     % {'exp', 'rnd', 'none'}. Set to 'none' to ignore. 'exp' generates powers in an exponential sequence: p1 = powerDecay^0, p2=powerDecay^1, ... while 'rand' generates random losses in [par.powerDecay, 1] for path=2,.. i.e. p1 = 1, p2=powerDecay*rand, ..

% Signals
par.signal_length = 2^15;   % Long enough to fit blocks*(snapshot+interboock)
par.mod = 'QPSK';           % Modulation scheme in {'QPSK'}, more added if needed
par.interblock = [0 20];    % Delay between TS blocks. Samples randomly in this range [a b]. 
par.pathdelay = [0 0];      % Delay between different recieved paths. If this is a range [a,b], samples randomly. 
par.blocks = par.K+1;       % Number of blocks for TS averaging. (par.K+1 is a good default)
par.blockSweep = 1:10;      % Specific to runtype 'blocksweep': varies number of blocks each run
par.snapshot = 2^10;        % Snapshot window for signals
par.SNR = 30;               % Static SNR for runtype 'single', overwritten for sweeps
par.snrSweep = -20:2:40;    % Variable SNR for runtypes other than 'single'
par.minSep = pi/10;         % Minimum separation between azimuth AoA (in radians)
par.aziRange = [pi/20, 39*pi/20]; % Specify range that signal azimuths must be in
par.eleRange = [pi/20, 19*pi/20]; % Same for signal elevations
par.polType = 'rnd';        % {'rnd_lin', 'rnd'}
par.forcePol = [];          % Force all signals to this polarization

% Estimation & Statistics
par.aziEstRange = [0,2*pi]; % Bounds for estimating (should only effect music-type estimators)
par.eleEstRange = [0,pi];   % Bounds for estimating (should only effect music-type estimators)
par.res = 2^8 / pi;         % Resolution for MUSIC-type estimators
par.recombination = 'rnd';  % {'rnd','max'}
par.scrub = true;           % Additional variable for ml generation - removes RMSE outliers - NOTE: If true, then X (output from ml_gen) takes on the "clean" values
par.sampling = 1;           % WORK IN PROGRESS DON'T CHANGE - Calculate statistics for this percent of trials (includes DoA estimation & statistics) - lower value speeds up processing. Only implemented in "ml_gen"                   

% Simulation
par.smooth = ["ts", "nts"]; %Whether to do temporal smoothing (ts) or not. Can also do both (ts,nts)
par.accelerate = 1;         % See readme

% Fusion
par.fusion.interval = 2;    % How many degrees per fusion interval (in degrees for human intepretability)
par.fusion.res = 5;         % How many points to test in each interval
par.fusion.model_order_type = "two_step"; % 1 > one step, 2 > two step
par.fusion.doa_type = "classifier"; % classifier or regression

% Paths to model directories (dummies for now)
par.fusion.interval_model_path = 0;
par.fusion.order_model_path = 0;
par.fusion.doa_model_path = 0; 
par.fusion.two_step_path = "models\2_step_classifier_trained_model";
par.fusion.one_step_path = "models\18_classifier_trained_model";

% Save configuration
par.saveFlag = 1;           % Should anything be saved? NOTE: Overridden to 0 for runtype 'single'
par.saveLight = 0;          % Should everything be saved or just statistics? 
par.saveName = "test";   % Appends this to the name for ease of identification (edited later to be unique)
par.savePlots = 1;          % Save plots as fig, eps, png


%% 
%---------------------\\-------------------------
% Ensure that the specified parameters are compatible
%---------------------\\------------------------- 
%%
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
if mod(180,par.fusion.interval)~=0 % Calculate an appropriate fusion interval, that is close to the requested value. Also convert to radians. 
    par.fusion.interval = 180/floor(180/par.fusion.interval);
end
par.fusion.interval = par.fusion.interval/180 * pi; % convert to radians for compatability and speed

par.fusion.res = par.fusion.res/par.fusion.interval;

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
%%
%if any(strcmp(par.estimator,'drmusic')) || any(strcmp(par.estimator,'music')) % Get the necessary manifolds (single for just MUSIC, interval-based for fusion
[par.AOA_function, par.AOA_function_H, par.phi, par.tht, par.DOA_function, par.DOA_function_H] = get_manifold(par);
%end
% if 1==0%any(strcmp(par.estimator,'fusion')) %change back to just fusion
% [par.intervals.DOA_function, par.intervals.DOA_function_H, par.intervals.AOA, par.intervals.AOA_H, par.intervals.phi, par.intervals.tht] = get_manifold_on_intervals(par);
% end
%%
%---------------------\\-------------------------
% If estimating using ML models, load those models
%---------------------\\-------------------------
%%
% Check if ML being used, load appropriate models to par
if any(strcmp(par.estimator, "fusion"))
    S = warning('off','nnet_cnn:internal:cnn:analyzer:NetworkAnalyzer:NetworkHasWarnings');
    switch par.fusion.model_order_type
        case "one_step"
            modelfile = par.fusion.one_step_path + "\18classifier.h5";
            weightfile = par.fusion.one_step_path + "\18classifier.json";
            par.fusion.model_order.classifier           = importKerasNetwork(weightfile, 'WeightFile', modelfile, 'OutputLayerType', 'Classification');
        case "two_step"
            par.fusion.model_order.step_one             = importKerasNetwork(par.fusion.two_step_path+"\1_STEP\model_TF.json", 'WeightFile', par.fusion.two_step_path+"\1_STEP\model_TF.h5", 'OutputLayerType', 'Classification');
            par.fusion.model_order.step_two.two_sig     = importKerasNetwork(par.fusion.two_step_path+"\2_STEP_NSIG_2\model_TF.json", 'WeightFile', par.fusion.two_step_path+"\2_STEP_NSIG_2\model_TF.h5", 'OutputLayerType', 'Classification');
            par.fusion.model_order.step_two.three_sig   = importKerasNetwork(par.fusion.two_step_path+"\2_STEP_NSIG_3\model_TF.json", 'WeightFile', par.fusion.two_step_path+"\2_STEP_NSIG_3\model_TF.h5", 'OutputLayerType', 'Classification');
            par.fusion.model_order.step_two.four_sig    = importKerasNetwork(par.fusion.two_step_path+"\2_STEP_NSIG_4\model_TF.json", 'WeightFile', par.fusion.two_step_path+"\2_STEP_NSIG_4\model_TF.h5", 'OutputLayerType', 'Classification');
            par.fusion.model_order.step_two.five_sig    = importKerasNetwork(par.fusion.two_step_path+"\2_STEP_NSIG_5\model_TF.json", 'WeightFile', par.fusion.two_step_path+"\2_STEP_NSIG_5\model_TF.h5", 'OutputLayerType', 'Classification');
    end
    disp("Ignore those warnings for now");
end