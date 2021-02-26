function [X, Y, Y_order, RMSE, stats] = run_ml_gen(par)
clc
X = zeros(par.Trials, 42);
X_nts = X;

sample=1;
switch par.type
    case '1d'
        RMSCE_ts = zeros(1, par.Trials);
        
        RMSE_ts = zeros(1, par.Trials);
        
        norecom_RMSE_ts = zeros(1, par.Trials);
        
        recom_percent_ts = zeros(1, par.Trials);
        if isempty(par.forcePath)
            Y = zeros(par.Trials, max(par.K_range));
        else
            Y = zeros(par.Trials, length(par.forcePath));
        end
    case '2d'
        RMSAE_ts = zeros(1, par.Trials);
        
        norecom_RMSAE_ts = zeros(1, par.Trials);
        if isempty(par.forcePath)
            Y = zeros(par.Trials, par.K*2);
        else
            Y = zeros(par.Trials, par.K*2);
        end
end
for trial = 1:par.Trials
    
    % Display Progress
    if mod(trial, round(0.01*par.Trials))==0
        clc
        fprintf("%i%% complete", 100*trial/par.Trials)
    end
    par.SNR = (par.snrSweep(end)-par.snrSweep(1)).*rand+par.snrSweep(1);
    if ~isempty(par.K_range)
        par.K   = par.K_range(trial);
    end
    
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
    switch par.type 
        case '1d'
            Y(trial,:) = [paths.AoA(:,2), zeros(1, size(Y,2)-par.K)];
        case '2d'
            Y(trial,:) = [sort(paths.AoA(:,1));sort(paths.AoA(:,2)); zeros(size(Y,2) - par.K*2, 1)]';
    end
    
    
    X_nts(trial,:) = S_nts;
    r = rand;
    if r<=par.sampling
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
    else
        stats = [];
    end
end
if par.sampling  % TODO: bad! fix this.
    switch par.type
        case '1d'
            RMSE = RMSE_ts;
        case '2d'
            RMSE = RMSAE_ts;
    end
end
if par.scrub==1 % Filter outliers away
    dirty_X = X;
    dirty_Y = Y;
    
    clean_idx = find(RMSE<5);
    X   = X(clean_idx, :);
    if par.type == '1d'
        Y   = [pi/2 * ones(size(Y(clean_idx, :))), sort(Y(clean_idx, :), 2)];
    elseif par.type == '2d'
        Y   = [sort(Y(clean_idx,1:par.K),2), sort(Y(clean_idx,par.K+1:end),2)];
    end
    
    X_nts_clean = X_nts(clean_idx, :);
else
    if par.type == '1d'
        Y   = [pi/2 * ones(size(Y(:, :))), sort(Y(:, :), 2)];
    elseif par.type == '2d'
        %Y   = [sort(Y(:,1:par.K),2), sort(Y(:,par.K+1:end),2)];
        
    end
end
%Y_doa = Y;
Y_order = par.K_range';

if ~exist('RMSE', 'var')
    RMSE = [];
end
end

