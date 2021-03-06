function [] = run_block_sweep(par)
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
end

