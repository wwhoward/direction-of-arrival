# direction-of-arrival
Project code for DoA estimation

Will W Howard
wwhoward@vt.edu


MATLAB

main.m contains four "runtypes" : 'single', 'snr_sweep', 'block_sweep', 'ml_gen'

'single' 	- Runs a single trial with given parameters (see following section), and presents plots of 
		temporal smoothing vs. no smoothing

'snr_sweep' 	- Runs for length(par.snrSweep)*par.Trials, then presents plots of temporally smoothed DR-MUSIC RMSE 
		averaged for each SNR value. 

'block_sweep' 	- Runs for length(par.blockSweep)*par.Trials, then presents plots of temporally smoothed DR-MUSIC RMSE
		averaged for each SNR value

'ml_gen' 	- Easy way to generate a par.Trials-sized dataset for training ML models, as well as providing RMSE 
		values to validate the dataset. Optional settings can configure whether any trials are removed from the 
		resulting dataset. 

Parameters

All of these parameters are contained in the struture par. 

% Signal settings
par.signal_length = 2^13;   % Long enough to fit blocks*(snapshot+interboock)
par.mod = 'QPSK';           % Modulation scheme in {'QPSK'}, more added if needed
par.interblock = [5 5];     % Delay between TS blocks. Samples randomly in this range [a b]. 
par.pathdelay = [0 0];      % Delay between different recieved paths. If this is a range [a,b], samples randomly. 
par.blocks = 3;             % Number of blocks for TS averaging. (par.K+1 is a good default)
par.blockSweep = 1:10;      % Specific to runtype 'blocksweep': varies number of blocks each run
par.snapshot = 2^9;         % Snapshot window for signals
par.SNR = 10;               % Static SNR for runtype 'single', overwritten for sweeps
par.snrSweep = 10:30;       % Variable SNR for runtypes other than 'single'

% Estimation & Statistics
par.res = 2^6 / pi;         % Resolution for MUSIC-type estimators
par.recombination = 'rnd';  % {'rnd','max'}
par.scrub = true;           % Additional variable for ml generation - removes RMSE outliers - NOTE: If true, then X (output from ml_gen) takes on the "clean" values
par.sampling = 1;           % WORK IN PROGRESS DON'T CHANGE - Calculate statistics for this percent of trials (includes DoA estimation & statistics) - lower value speeds up processing. Only implemented in "ml_gen"                   

% Simulation
par.Trials = 1000;          % During a sweep, how many trials for each parameter
par.runtype = 'single';     % {'single','snr_sweep','block_sweep','ml_gen'} - names are fairly self-explainitory

% Save configuration
par.saveFlag = 1;           % Should anything be saved? NOTE: Overridden to 0 for runtype 'single'
par.saveLight = 0;          % Should everything be saved or just statistics? 
par.saveName = "ml_clean";  % Appends this to the name for ease of identification (edited later to be unique)
par.savePlots = 0;          % Save plots as fig, eps, png


How Temporal Smoothing works

Temporal smoothing is implemented with the goal of removing multipath-induced correlation from recieved signals. 

Temporal smoothing relies on several "blocks" taken from a recieved signal, with the time displacement between the beginning of each
being greater than the channel coherence time. The covariance matrix for each block is calculated, then averaged together to form the 
signal covariance. 

We simulate the coherence time of the channel by multiplying each block by a different random complex gain. 


Multipath Simulation

We use the following notation to represent multipath scenarios: Each source is assigned a number, starting with 1, then the number of paths for that source is represented by multiplicity of it's number. So, for one signal on two paths, we have [1 1]. 

To simulate multipath reception of a signal, we use the delay-then-sum technique. Path delays are specified as a parameter (above), and can be randomized in a range. 