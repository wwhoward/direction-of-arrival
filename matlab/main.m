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
par.forcePath = [1 1];      % set par.forcePath = [] for random path configuration. Overrides par.K, so be careful using this
par.type = '1d';            % Problem dimensionality in {'1d', '2d'}. Effects MUSIC-type searches as well as DoA ground truth assignments. 

% Signal settings
par.signal_length = 2^13;   % Long enough to fit blocks*(snapshot+interboock)
par.mod = 'QPSK';           % Modulation scheme in {'QPSK'}, more added if needed
par.interblock = [5 10];    % Delay between TS blocks. Samples randomly in this range [a b]. 
par.pathdelay = [1 3];      % Delay between different recieved paths. If this is a range [a,b], samples randomly. 
par.blocks = 3;             % Number of blocks for TS averaging. (par.K+1 is a good default)
par.blockSweep = 1:10;      % Specific to runtype 'blocksweep': varies number of blocks each run
par.snapshot = 2^9;         % Snapshot window for signals
par.SNR = 20;               % Static SNR for runtype 'single'
par.snrSweep = 10:30;       % Variable SNR for runtypes other than 'single'

% Estimation & Statistics
par.res = 2^12 / pi;        % Resolution for MUSIC-type estimators
par.recombination = 'rnd';  % {'rnd','max'}

% Simulation
par.Trials = 10000;         % During a sweep, how many trials for each parameter
par.runtype = 'single';     % {'single','snr_sweep','block_sweep','ml_gen'} - names are fairly self-explainitory

% Save configuration
par.saveFlag = 0;           % Should anything be saved? NOTE: Overridden to 0 for runtype 'single'
par.saveLight = 1;          % Should everything be saved or just statistics? 
par.saveName = "ml_test";   % Appends this to the name for seperability (edited later to be unique)
par.savePlots = 0;          % Save plots as fig, eps, png






%% 

switch par.runtype
    case 'single'
        
    par.saveFlag = 0; % Nothing to save for single runs, so set saveFlag to 0
    paths = assign_paths(par); % Assign directions of arrival and multipath state


    signal = transmitter(paths, par); % Assign waveforms to each signal


    signal = reciever(signal, paths, par); % Recieve waveforms with the specified multipath state


    signal = temporal_smooth(signal, par); % Calculate the temporally smoothed covariance matrix


    est_ts = est_drmusic(signal.ts, par); % Estimate the DoA for the TS case
    est_nts = est_drmusic(signal.R, par); % Estimate the DoA without smoothing

    stats_ts = calc_stats(est_ts, signal, par, paths);
    stats_nts = calc_stats(est_nts, signal, par, paths);

    % Plot the results! 
    plott(est_ts, stats_ts, paths, 'Temporal Smoothing - Three Signals - DR-MUSIC')
    plott(est_nts, stats_nts, paths, 'No Smoothing - Three Signals - DR-MUSIC')

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
    X = zeros(par.Trials, 42);
    if isempty(par.forcePath)
        Y = zeros(par.Trials, par.K);
    else
        Y = zeros(par.Trials, length(par.forcePath));
    end
    
    for trial = 1:par.Trials
        par.SNR = (par.snrSweep(end)-par.snrSweep(1)).*rand+par.snrSweep(1);

        paths = assign_paths(par);
        signal = transmitter(paths, par);
        signal = reciever(signal, paths, par);
        signal = temporal_smooth(signal, par);

        est_ts = est_drmusic(signal.ts, par);
        est_nts = est_drmusic(signal.R, par);

        stats.ts(trial) = calc_stats(est_ts, signal, par, paths);
        stats.nts(trial) = calc_stats(est_nts, signal, par, paths);
        
        S2 = triu(signal.ts);
        S2 = nonzeros(S2);
        S = [real(S2);imag(S2)];
        
        X(trial,:) = S;
        Y(trial,:) = paths.AoA(:,2);
        
        RMSCE_ts(trial) = stats.ts(trial).rmsce_deg;
        RMSCE_nts(trial) = stats.nts(trial).rmsce_deg;

        RMSE_ts(trial) = stats.ts(trial).rmse_deg;
        RMSE_nts(trial) = stats.nts(trial).rmse_deg;

        norecom_RMSE_ts(trial) = stats.ts(trial).norecom_rmse_deg;
        norecom_RMSE_nts(trial) = stats.nts(trial).norecom_rmse_deg;

        recom_percent_ts(trial) = stats.ts(trial).recom_percent;
        recom_percent_nts(trial) = stats.nts(trial).recom_percent;
        
        
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
    if par.savePlots
        mkdir(savename)
        FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
        for f = 1:length(FigList)
            savefig(FigList(f), fullfile(savename, [get(FigList(f),'Name'), '.fig'])); %savename+string(get(FigList(f),'Name'))+'.fig') 
        end
    end
end
%% Functions that make it do the work
function [] = plott(x, stats, paths, t)
figure
plot(x.phi, x.spectrum)
hold on; 

for p=1:length(paths.AoA(:,2)); xline(paths.AoA(p,2), '--r'); end
for e=1:length(x.peaks); plot(x.peaks(e), x.peakvals(e), 'ro', 'MarkerSize', 10); end
title(t)
set(gca,'XTick',0:pi/2:2*pi)
set(gca,'XTickLabel',{'0','pi/2','pi','3*pi/2','2*pi'})
xlabel('Azimuth')
end

function [] = plot_rmse(rmse, snrSweep)
    

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
        
        [est.peakvals, peak_idx] = findpeaks(est.spectrum);
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
    case '2d'
        paths.AoA = [];
        for k=1:paths.sources
            for p=1:paths.multi(k)
                [azi,ele] = RandUniformSphere([1,1]);
                paths.AoA(:, :) = [paths.AoA [azi, ele]'];
                
            end
        end
    case '1d'
        paths.AoA = [];
        for k=1:paths.sources
            for p = 1:paths.multi(k)
                azi = rand*2*pi;
                ele = pi/2 + 0.01*rand;
                paths.AoA = [paths.AoA,; [ele, azi]];
                [~,~,paths.signal_vector(k).path(p,:)] = VectorSensor([ele,azi], [pi/2*rand, 2*pi*rand-pi]);
            end
        end
end

end