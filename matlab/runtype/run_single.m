% Temporal smoothing (TS) trials for the vector sensor
% Summer 2020
% Will Howard, {wwhoward}@vt.edu
% wireless @ VT

function [paths, signal, est, stats] = run_single(par)
tic

if ~isempty(par.forcePath)
    par.K = length(par.forcePath);
end

paths = assign_paths(par); % Assign directions of arrival and multipath state

signal = transmitter(paths, par); % Assign waveforms to each signal

signal = reciever(signal, paths, par); % Recieve waveforms with the specified multipath state

signal = temporal_smooth(signal, par); % Calculate the temporally smoothed covariance matrix

est = estimate(signal, par, paths);

for estimator = par.estimator
    for smooth = par.smooth
        stats.(estimator).(smooth) = calc_stats(est.(estimator).(smooth), signal, par, paths);
    end
end

if strcmp(par.type, '2d') && any(strcmp(par.estimator, 'music')) % can't do assignment for 1d, can't do MUSIC for rndpol. Could write dr-music on intervals, but only doing assignment for ASI so no purpose. 
    for estimator = par.estimator
        [est.associate.weights.(estimator), est.associate.idx.(estimator), est.associate.assignment.(estimator)] = toy_associate(par, paths, est, signal, stats.(estimator).ts);
        %[weights2, idx2] = toy_associate(par, paths, est, signal, stats.music.ts);
    end
end

plott(est, paths, par, stats)

for estimator = par.estimator
    if estimator == "music"
        str_e = "MUSIC   ";
    elseif estimator == "drmusic"
        str_e = "DR-MUSIC";
    elseif estimator == "fusion"
        str_e = "Fusion  ";
    end
    for smooth = par.smooth
        if smooth == "ts"
            str_s = "Temporal Smoothing";
        elseif smooth == "nts"
            str_s = "No Smoothing      ";
        end
        switch par.type
            case '2d'
                fprintf('%s \t %s \t mean square angular error: %f degrees \t (with %i recombinations) \n', ...
                    str_s, str_e, ...
                    stats.(estimator).(smooth).rmsae_deg, sum(stats.(estimator).(smooth).recom_flag));
                
                %disp(est.(estimator).(smooth).smoothing + " " + est.(estimator).(smooth).label + " " + "mean square angular error: " + ...
                    %string(stats.(estimator).(smooth).rmsae_deg) + " degrees (with " + string(sum(stats.(estimator).(smooth).recom_flag)) + ...
                    %" recombinations)")
            case '1d'
                disp(est.(estimator).(smooth).smoothing + " " + est.(estimator).(smooth).label + " " + "mean square error: " + ...
                    string(stats.(estimator).(smooth).rmse_deg) + " degrees (with " + string(sum(stats.(estimator).(smooth).recom_flag)) + ...
                    " recombinations)")
        end
    end
end

% switch par.type
%     case '2d'
%         display("TS Root mean square angular error: " + string(stats_ts.rmsae_deg) + " degrees (with " + string(sum(stats_ts.recom_flag)) + " recombinations)")
%         display("NTS Root mean square angular error: " + string(stats_nts.rmsae_deg) + " degrees (with " + string(sum(stats_nts.recom_flag)) + " recombinations)")
%     case '1d'
%         display("TS Root mean square error: " + string(stats_ts.rmse_deg) + " degrees (with " + string(sum(stats_ts.recom_flag)) + " recombinations)")
%         display("NTS Root mean square error: " + string(stats_nts.rmse_deg) + " degrees (with " + string(sum(stats_nts.recom_flag)) + " recombinations)")
% 
% end
%     display(stats_ts.rmse_deg)
%     display(stats_ts.rmsce_deg)
%     display(stats_nts.rmse_deg)
%     display(stats_nts.rmsce_deg)

end