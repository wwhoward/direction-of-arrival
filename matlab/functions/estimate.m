% Temporal smoothing (TS) trials for the vector sensor
% Fall 2020
% Will Howard, {wwhoward}@vt.edu
% wireless @ VT

function est = estimate(signal, par, paths) % Ultimately, don't send over paths
%ESTIMATE Summary of this function goes here
%   Detailed explanation goes here


for estimator = par.estimator
    for smooth = par.smooth
        if smooth == "ts"
            R = signal.ts;
            str = "Temporal Smoothing";
        elseif smooth == "nts"
            R = signal.R;
            str = "No Smoothing";
        end
        switch estimator
            case 'drmusic'
                est.(estimator).(smooth) = est_drmusic(R, par);                
                est.(estimator).(smooth).label = "DR-MUSIC";
            case 'music'
                est.(estimator).(smooth) = est_music(R, par);
                est.(estimator).(smooth).label = "MUSIC";
            case 'fusion'
                est.(estimator).(smooth) = est_fusion(R, par, paths, signal); % Ultimately, don't send over paths
                est.(estimator).(smooth).label = "FUSION";
        end        
        est.(estimator).(smooth).smoothing = str;
    end
end
end

