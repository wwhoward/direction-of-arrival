function [X, Y, RMSE, stats] = run_ml_gen(par)
clc
X = zeros(par.Trials, 42);
X_nts = X;
if isempty(par.forcePath)
    Y = zeros(par.Trials, par.K);
else
    Y = zeros(par.Trials, length(par.forcePath));
end

sample=1;
switch par.type
    case '1d'
        RMSCE_ts = zeros(1, par.Trials);
        
        RMSE_ts = zeros(1, par.Trials);
        
        norecom_RMSE_ts = zeros(1, par.Trials);
        
        recom_percent_ts = zeros(1, par.Trials);
    case '2d'
        RMSAE_ts = zeros(1, par.Trials);
        
        norecom_RMSAE_ts = zeros(1, par.Trials);
end
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
    
    S2_nts = triu(signal.R);
    S2_nts = nonzeros(S2_nts);
    S_nts = [real(S2_nts); imag(S2_nts)];
    
    X(trial,:) = S;
    Y(trial,:) = paths.AoA(:,2);
    
    X_nts(trial,:) = S_nts;
    
    if rand<=par.sampling
        est_ts = est_drmusic(signal.ts, par);
        
        stats.ts(sample) = calc_stats(est_ts, signal, par, paths);
        if par.type == '1d'            
            RMSCE_ts(sample) = stats.ts(sample).rmsce_deg;
            
            RMSE_ts(sample) = stats.ts(sample).rmse_deg;
            
            norecom_RMSE_ts(sample) = stats.ts(sample).norecom_rmse_deg;
            
            recom_percent_ts(sample) = stats.ts(sample).recom_percent;
            
        elseif par.type == '2d'
            RMSAE_ts(sample) = stats.ts(sample).rmsae_deg;
            
            norecom_RMSAE_ts(sample) = stats.ts(sample).norecom_rmsae_deg;
        end
        sample=sample+1;
    end
end
switch par.type
    case '1d'
        RMSE = RMSE_ts;
    case '2d'
        RMSE = RMSAE_ts;
end
if par.scrub==1 % Filter outliers away
    dirty_X = X;
    dirty_Y = Y;
    
    clean_idx = find(RMSE<5);
    X   = X(clean_idx, :);
    Y   = sort(Y(clean_idx, :),2);
    
    X_nts_clean = X_nts(clean_idx, :);
end

end

