% Temporal smoothing (TS) trials for the vector sensor
% Summer 2020
% Will Howard, {wwhoward}@vt.edu
% wireless @ VT

function [paths, signal, est_ts, est_nts, stats_ts, stats_nts] = run_single(par)
tic

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
plott(est_ts, paths, par, str1)
plott(est_nts, paths, par, str2)

%elapsedTime = toc
switch par.type
    case '2d'
        display("TS Root mean square angular error: " + string(stats_ts.rmsae_deg) + " degrees (with " + string(sum(stats_ts.recom_flag)) + " recombinations)")
        display("NTS Root mean square angular error: " + string(stats_nts.rmsae_deg) + " degrees (with " + string(sum(stats_nts.recom_flag)) + " recombinations)")
    case '1d'
        display("TS Root mean square error: " + string(stats_ts.rmse_deg) + " degrees (with " + string(sum(stats_ts.recom_flag)) + " recombinations)")
        display("NTS Root mean square error: " + string(stats_nts.rmse_deg) + " degrees (with " + string(sum(stats_nts.recom_flag)) + " recombinations)")

end
%     display(stats_ts.rmse_deg)
%     display(stats_ts.rmsce_deg)
%     display(stats_nts.rmse_deg)
%     display(stats_nts.rmsce_deg)

end