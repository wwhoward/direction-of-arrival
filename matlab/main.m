% Temporal smoothing (TS) trials for the vector sensor
% Summer 2020
% Last updated: 9/11/2020
% Will Howard, {wwhoward}@vt.edu
% wireless @ VT

% Notation:
% A multipath state of [1 1] indicates one signal on two paths
% Similarly, [1 2] indicates two signals on one path each and [1 1 2 3] indicates one signal on two paths + two more signals on one path each



clc,clear,close all
addpath('data', 'Helper'); 

%% Set parameters! 
% Multipath Scenarios
par.K = 3;                  % Keep K <= 4: This is the total incident paths on the sensor. This is overwritten if par.forcePath is not empty. 
par.forceMulti = false;     % If true, enables path states of [1], [1 2], etc
par.forcePath = [1 1];    % set par.forcePath = [] for random path configuration. Overrides par.K, so be careful using this
par.type = '1d';            % Problem dimensionality in {'1d', '2d'}. Effects MUSIC-type searches as well as DoA ground truth assignments. 
par.minSep = pi/10;         % Minimum separation between azimuth AoA (in radians)

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






%% 

switch par.runtype
    case 'single'
        
    par.saveFlag = 0; % Nothing to save for single runs, so set saveFlag to 0
    
    if ~isempty(par.forcePath)
        par.K = length(par.forcePath);
    end
    
    paths = assign_paths(par); % Assign directions of arrival and multipath state


    signal = transmitter(paths, par); % Assign waveforms to each signal


    signal = reciever(signal, paths, par); % Recieve waveforms with the specified multipath state


    signal = temporal_smooth(signal, par); % Calculate the temporally smoothed covariance matrix


    est_ts = est_drmusic(signal.ts, par); % Estimate the DoA for the TS case
    est_nts = est_drmusic(signal.R, par); % Estimate the DoA without smoothing

    stats_ts = calc_stats(est_ts, signal, par, paths);
    stats_nts = calc_stats(est_nts, signal, par, paths);

    % Plot the results! 
    
    str1 = 'Temporal Smoothing - '+string(par.K)+' Signals - DR-MUSIC';
    str2 = 'No Smoothing - '+string(par.K)+' Signals - DR-MUSIC';
    plott(est_ts, stats_ts, paths, str1)
    plott(est_nts, stats_nts, paths, str2)

    display(stats_ts.rmse_deg)
    display(stats_ts.rmsce_deg)
    display(stats_nts.rmse_deg)
    display(stats_nts.rmsce_deg)

%% Now try a complicated run where we average the results from 500 trials of path=[1,1,2]

    case 'snr_sweep'
    for snr=1:length(par.snrSweep)
        par.SNR = par.snrSweep(snr);
        for trial = 1:par.Trials
            paths = assign_paths(par);
            signal = transmitter(paths, par);
            signal = reciever(signal, paths, par);
            signal = temporal_smooth(signal, par);

            est_ts = est_drmusic(signal.ts, par);
            est_nts = est_drmusic(signal.R, par);

            stats.ts(trial, snr) = calc_stats(est_ts, signal, par, paths);
            stats.nts(trial, snr) = calc_stats(est_nts, signal, par, paths);

            RMSCE_ts(trial, snr) = stats.ts(trial, snr).rmsce_deg;
            RMSCE_nts(trial, snr) = stats.nts(trial, snr).rmsce_deg;

            RMSE_ts(trial, snr) = stats.ts(trial, snr).rmse_deg;
            RMSE_nts(trial, snr) = stats.nts(trial, snr).rmse_deg;
            
            norecom_RMSE_ts(trial, snr) = stats.ts(trial,snr).norecom_rmse_deg;
            norecom_RMSE_nts(trial,snr) = stats.nts(trial,snr).norecom_rmse_deg;
            
            recom_percent_ts(trial,snr) = stats.ts(trial,snr).recom_percent;
            recom_percent_nts(trial,snr) = stats.nts(trial,snr).recom_percent;
    %         if trial == 1 && snr == 1
    %             plott(est_ts, stats.ts, paths, 'Temporal Smoothing - Four Signals - DR-MUSIC')
    %         end
        end
        display(par.SNR)
    end

    % plott(est_ts, stats.ts, paths, 'Temporal Smoothing - Four Signals - DR-MUSIC')
    Recom_percent = mean(recom_percent_ts);
    % Plot those things
    rmse_ts_ave = mean(RMSE_ts);
    rmse_nts_ave = mean(RMSE_nts);
    norecom_rmse_ts_ave = mean(norecom_RMSE_ts);
    norecom_rmse_nts_ave = mean(norecom_RMSE_nts);

    rmsce_ts_ave = mean(RMSCE_ts);
    rmsce_nts_ave = mean(RMSCE_nts);

    figure('Name','fig_rmsce'); 
    semilogy(par.snrSweep, rmsce_ts_ave); hold on; grid on;
    semilogy(par.snrSweep, rmsce_nts_ave);
    title('Temporal Smoothing vs. No Smoothing')
    legend('Temporal Smoothing, 4 Blocks', 'No Smoothing')
    ylim([10^-2 10^3])
    xlabel('SNR')
    ylabel('Corrected RMSE')
    annotation('textbox', [0.2, 0.2, 0.2, 0.2], 'string', string(par.forcePath))

    figure('Name','fig_rmse'); 
    semilogy(par.snrSweep, rmse_ts_ave); hold on; grid on;
    semilogy(par.snrSweep, rmse_nts_ave);
    title('Temporal Smoothing vs. No Smoothing')
    legend('Temporal Smoothing, 4 Blocks', 'No Smoothing')
    ylim([10^-2 10^3])
    xlabel('SNR')
    ylabel('RMSE')
    annotation('textbox', [0.2, 0.2, 0.1, 0.1], 'string', 'One signal, four paths')
    
    figure('Name','fig_norecom_rmse'); 
    semilogy(par.snrSweep, norecom_rmse_ts_ave); hold on; grid on;
    semilogy(par.snrSweep, norecom_rmse_nts_ave);
    title('Temporal Smoothing vs. No Smoothing : No Recombinations')
    legend('Temporal Smoothing, 4 Blocks', 'No Smoothing')
    ylim([10^-2 10^3])
    xlabel('SNR')
    ylabel('RMSE')
    annotation('textbox', [0.2, 0.2, 0.1, 0.1], 'string', 'One signal, two paths')

%% Let's try sweeping the number of blocks
    case 'block_sweep'
    for b=1:length(par.blockSweep)
        par.blocks = par.blockSweep(b);
        for trial = 1:par.Trials
            paths = assign_paths(par);
            signal = transmitter(paths, par);
            signal = reciever(signal, paths, par);
            signal = temporal_smooth(signal, par);

            est_ts = est_drmusic(signal.ts, par);
            est_nts = est_drmusic(signal.R, par);

            stats.ts(trial, b) = calc_stats(est_ts, signal, par, paths);
            stats.nts(trial, b) = calc_stats(est_nts, signal, par, paths);

            RMSCE_ts(trial, b) = stats.ts(trial, b).rmsce_deg;
            RMSCE_nts(trial, b) = stats.nts(trial, b).rmsce_deg;

            RMSE_ts(trial, b) = stats.ts(trial, b).rmse_deg;
            RMSE_nts(trial, b) = stats.nts(trial, b).rmse_deg;
    %         if trial == 1 && snr == 1
    %             plott(est_ts, stats.ts, paths, 'Temporal Smoothing - Four Signals - DR-MUSIC')
    %         end
        end
        display(par.blocks)
    end

    rmse_ts_ave = mean(RMSE_ts);
    rmse_nts_ave = mean(RMSE_nts);

    rmsce_ts_ave = mean(RMSCE_ts);
    rmsce_nts_ave = mean(RMSCE_nts);

    figure('Name','fig_rmsce'); 
    semilogy(rmsce_ts_ave); hold on; grid on;
    %xticks(1:6)
    %semilogy(par.blockSweep, rmsce_nts_ave);
    title('Temporal Smoothing - Three Signals - DR-MUSIC')
    ylim([10^-1 10^3])
    xlabel('Number of Blocks')
    ylabel('RMSE')
    
    
%% Generate ML dataset
    case 'ml_gen'
    clc    
    X = zeros(par.Trials, 42);
    if isempty(par.forcePath)
        Y = zeros(par.Trials, par.K);
    else
        Y = zeros(par.Trials, length(par.forcePath));
    end
    
    sample=1;
    for trial = 1:par.Trials
        
        % Display Progress
        if mod(trial, round(0.1*par.Trials))==0
            clc
            fprintf("%i%% complete", 100*trial/par.Trials)
        end
        par.SNR = (par.snrSweep(end)-par.snrSweep(1)).*rand+par.snrSweep(1);

        paths = assign_paths(par);
        signal = transmitter(paths, par);
        signal = reciever(signal, paths, par);
        signal = temporal_smooth(signal, par);
        
        S2 = triu(signal.ts);
        S2 = nonzeros(S2);
        S = [real(S2);imag(S2)];
        
        X(trial,:) = S;
        Y(trial,:) = paths.AoA(:,2);
        
        if rand<=par.sampling
            est_ts = est_drmusic(signal.ts, par);
            est_nts = est_drmusic(signal.R, par);

            stats.ts(sample) = calc_stats(est_ts, signal, par, paths);
            stats.nts(sample) = calc_stats(est_nts, signal, par, paths);

            RMSCE_ts(sample) = stats.ts(sample).rmsce_deg;
            RMSCE_nts(sample) = stats.nts(sample).rmsce_deg;

            RMSE_ts(sample) = stats.ts(sample).rmse_deg;
            RMSE_nts(sample) = stats.nts(sample).rmse_deg;

            norecom_RMSE_ts(sample) = stats.ts(sample).norecom_rmse_deg;
            norecom_RMSE_nts(sample) = stats.nts(sample).norecom_rmse_deg;

            recom_percent_ts(sample) = stats.ts(sample).recom_percent;
            recom_percent_nts(sample) = stats.nts(sample).recom_percent;
            
            sample=sample+1;
        end
    end
    
    if par.scrub==1 % Filter outliers away
        dirty_X = X;
        dirty_Y = Y;
        
        clean_idx = find(RMSE_ts<5);
        X   = X(clean_idx, :);
        Y   = Y(clean_idx, :);
    end
end



%% Save shtuff
if par.saveFlag
    i=0;
    savename = "data/"+date+par.saveName + "_" + string(i);
    while isfile(savename+".mat")
        i=i+1;
        savename = "data/"+date+par.saveName+"_"+string(i);
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
        end
    end
end

%% Functions (local for now, might seperate this large file into several in the future)
function [] = plott(x, stats, paths, t)
figure
plt = plot(x.phi, x.spectrum);
hold on; 

for p=1:length(paths.AoA(:,2)); l_plt(p)=xline(paths.AoA(p,2), '--r'); end
for e=1:min(length(paths.AoA(:,2)), length(x.peaks)); p_plt(e)=plot(x.peaks(e), x.peakvals(e), 'ro', 'MarkerSize', 10); end
title(t)
legend([plt, l_plt(1), p_plt(1)],'Spatial Spectrum', 'True Angle', 'Estimated Angle');
set(gca,'XTick',0:pi/2:2*pi)
set(gca,'XTickLabel',{'0','pi/2','pi','3*pi/2','2*pi'})
xlabel('Azimuth')
end

function [] = plot_rmse(rmse, snrSweep)
    % Work item

end

function stats = calc_stats(est, signal, par, paths)
switch par.type
    case '1d'
        u = paths.AoA(:,2)';
        uh = est.peaks;

        for i=1:length(u)
            if ~isempty(uh)
                [stats.err(i), closestIndex] = min(abs(u(i) - uh.'));
                closestValue(i) = uh(closestIndex);
                uh(closestIndex) = [];
                stats.corrected_err(i) = pi/180 * abs(CalculateAngleDifference(180/pi * u(i), 180/pi * closestValue(i), 'azi')); % this function requires units of degrees
                stats.recom_flag(i) = 0;
            else % if there are fewer than K estimates, we set the error of the extra ones to max (pi)
                if par.recombination == 'max'
                    stats.err(i) = pi;
                    stats.corrected_err(i) = pi;
                    closestValue = [];
                    stats.recom_flag(i) = 1;
                elseif par.recombination == 'rnd'
                    stats.err(i) = pi*rand;
                    stats.corrected_err(i) = stats.err(i);
                    closestValue = [];
                    stats.recom_flag(i) = 1;
                end
            end
            
        end
        
        
        stats.u = u;
        stats.uh = closestValue;
        
        stats.mse = sum(stats.err.^2)/length(stats.err);
        stats.rmse = sqrt(stats.mse);
        
        stats.err_deg = stats.err/pi * 180;
        stats.mse_deg = sum(stats.err_deg.^2)/length(stats.err_deg);
        stats.rmse_deg = sqrt(stats.mse_deg);
        
        
        stats.msce = sum(stats.corrected_err.^2)/length(stats.corrected_err);
        stats.rmsce = sqrt(stats.msce);
        
        stats.corrected_err_deg = stats.corrected_err/pi * 180;
        stats.msce_deg = sum(stats.corrected_err_deg.^2)/length(stats.corrected_err_deg);
        stats.rmsce_deg = sqrt(stats.msce_deg);
        
        
        % Now calculate these for only the non-recom cases
        stats.recom_percent = sum(stats.recom_flag)/length(stats.recom_flag);
        stats.norecom_mse = sum(stats.err(stats.recom_flag==0).^2)/length(stats.err(stats.recom_flag==0)); 
        stats.norecom_rmse = sqrt(stats.norecom_mse);
        
        stats.norecom_err_deg = stats.err(stats.recom_flag==0)/pi * 180;
        stats.norecom_mse_deg = sum(stats.norecom_err_deg.^2)/length(stats.norecom_err_deg);
        stats.norecom_rmse_deg = sqrt(stats.norecom_mse_deg);
        
        stats.SNR = par.SNR;
end
end

function est = est_drmusic(R, par)
% Input:
%   R: 6x6 Covariance Matrix

switch par.type
    case '1d'
        est.phi = 0:1/par.res:2*pi;
        for ph = 1:size(est.phi,2)
            [est.DOA_function(:,:,ph),~,~] = VectorSensor([pi/2,est.phi(ph)],[pi/4,0]);
        end
        
        [eigvect, eigval]=eig(R);
        [~,idx] = sort(diag(eigval));
        NoiseSpace = eigvect(:,idx(1:length(idx)-par.K));
        
        for i=1:size(est.DOA_function, 3)
            est.spectrum(i) = 1/min(real(eig((est.DOA_function(:,:,i)'*NoiseSpace*NoiseSpace'*est.DOA_function(:,:,i))))); % This should be real since B should be hermitian, but due to computational error, there will be some small imaginary part. Take abs() or real() to correct
        end
        
        [peakvals, peak_idx] = findpeaks(est.spectrum);
        [est.peakvals, idx] = sort(peakvals, 'descend');
        peak_idx = peak_idx(idx);
        est.peaks = est.phi(peak_idx);        
end
end

function signal = temporal_smooth(signal, par)

R_ts = zeros(6,6);
r = zeros(par.blocks, 6,6);
for b=1:par.blocks
    r(b,:,:) = 1/par.snapshot * squeeze(signal.rx(b,:,:))*squeeze(signal.rx(b,:,:))';
    R_ts(:,:) = 1/par.blocks * squeeze(r(b,:,:)) + R_ts(:,:);
end
signal.ts = R_ts;
end

function signal = reciever(signal, paths, par)
% Take modulated message, simulate multipath propegation, add AWGN

tx = signal.tx;
rx = zeros(par.blocks, 6, par.snapshot);
for b=1:par.blocks
    for k=1:paths.sources
        for p=1:paths.multi(k)
            if p ~= 1
                h = sqrt(0.5)*(randn+1j*randn);
            else
                h = 1;
            end
            a = paths.signal_vector(k).path(p,:)';
            del = randi(par.pathdelay);
            block = randi(par.interblock);
            sig = tx(k, 1+(b-1)*block+(p-1)*del : par.snapshot + (b-1)*block + (p-1)*del);
            rx(b, :,:) = h*a*sig + squeeze(rx(b,:,:));
%             display('got here')
            % rx(b,6,par.snapshot)=sqrt(0.5)*(randn+1j*randn)*paths.signal_vector(k).path(p,:)'*tx(k, b*par.interblock+(k-1)*par.pathdelay:par.snapshot-1+b*par.interblock+(k-1)*par.pathdelay)+rx(b,:,:);
        end
    end
end
signal.rx = awgn(rx, par.SNR);

signal.R = 1/par.snapshot * squeeze(signal.rx(1,:,:))*squeeze(signal.rx(1,:,:))';

end

function signal = transmitter(paths, par)
% Create messages and modulate them
% Number of bits transmitted per frame is set to be 1000. For QPSK
% modulation, this corresponds to 500 symbols per frame.
if nargin < 2
    par.mod = 'QPSK';
end

switch par.mod
    case 'QPSK'
        modu = comm.QPSKModulator( ...
            'BitInput',    true, ...
            'PhaseOffset', pi/4);
end

for k=1:paths.sources
    msg(k, :) = randi([0, 1], par.signal_length, 1);
    tx(k, :) = modu(msg(k,:)');
end


signal.tx = tx;
end

function paths = assign_paths(par)
% Creates path object based on inputs:
% par.K: total number of paths
% par.type: '1d or '2d'

% outputs path object with parameters
% paths.sources : total number of sources
% paths.multi : vector showing how many paths per source
% paths.AoA : angles of arrival, random on the uniform sphere
% paths.signal_vector(k).path(p) : signal vector for the p'th arrival of the k'th signal
K = par.K;

if par.forceMulti && isempty(par.forcePath)
    if K==1
        options=[1];
        r=1;
    elseif K==2
        options=[1,1];
        r=randi(1);
    elseif K==3
        options=[1,1,2;1,1,1];
        r=randi(2);
    elseif K==4
        options=[1,1,2,3;1,1,2,2;1,1,1,2;1,1,1,1];
        r=randi(4);
    elseif K==5
        options=[1,1,2,3,4;1,1,1,2,3;1,1,1,1,2;1,1,1,1,1;1,1,2,2,3;1,1,1,2,2];
        r=randi(6);
    end
elseif isempty(par.forcePath)
    if K==1
        options=[1];
        r=1;
    elseif K==2
        options=[1,2;1,1];
        r=randi(K);
    elseif K==3
        options=[1,2,3;1,1,2;1,1,1];
        r=randi(K);
    elseif K==4
        options=[1,2,3,4;1,1,2,3;1,1,2,2;1,1,1,2;1,1,1,1];
        r=randi(5);
    elseif K==5
        options=[1,2,3,4,5;1,1,2,3,4;1,1,1,2,3;1,1,1,1,2;1,1,1,1,1;1,1,2,2,3;1,1,1,2,2];
        r=randi(7);
    end
else
    options = par.forcePath;
    r=1;
end
Q=options(r,:);

paths.sources = max(Q);

for s = 1:paths.sources
    paths.multi(s) = sum(Q(:)==s);
end

% Assign azimuth and elevation (or if type = '1d', just azimuth)
switch par.type
    case '2d' % Needs minSep still
        paths.AoA = [];
        for k=1:paths.sources
            for p=1:paths.multi(k)
                [azi,ele] = RandUniformSphere([1,1]);
                paths.AoA(:, :) = [paths.AoA [azi, ele]'];
                
            end
        end
    case '1d'
        paths.AoA = [];
        Azi = assignAzi(par);
        minSepFlag=0;
        if min(diff(sort(Azi)))<par.minSep
            minSepFlag=1;
        end
        while ~isempty(par.minSep) && minSepFlag==1
            Azi = assignAzi(par);
            if min(diff(sort(Azi)))>=par.minSep
                minSepFlag=0;
            end
        end
        for k=1:paths.sources
            for p = 1:paths.multi(k)
                azi = Azi(1);
                ele = pi/2 + 0.01*rand;
                paths.AoA = [paths.AoA,; [ele, azi]];
                [~,~,paths.signal_vector(k).path(p,:)] = VectorSensor([ele,azi], [pi/2*rand, 2*pi*rand-pi]);
                Azi(1)=[];
            end
            
                
        end
        
end

end

function azi = assignAzi(par)
    if isempty(par.forcePath)
        k = par.K;
    else
        k = length(par.forcePath);
    end
    for i=1:k
        azi(i) = rand*2*pi;
    end
end