% Temporal smoothing (TS) trials for the vector sensor
% Summer 2020
% Will Howard, {wwhoward}@vt.edu
% wireless @ VT

function [paths, signal, est, stats] = run_snr_sweep(par)
if par.type == '2d'
    for snr=1:length(par.snrSweep)
        par.SNR = par.snrSweep(snr);
        for trial = 1:par.Trials
            paths = assign_paths(par);
            signal = transmitter(paths, par);
            signal = reciever(signal, paths, par);
            signal = temporal_smooth(signal, par);
            
            est.ts(trial, snr) = est_drmusic(signal.ts, par);
            est.nts(trial, snr) = est_drmusic(signal.R, par);
            
            stats.ts(trial, snr) = calc_stats(est.ts(trial, snr), signal, par, paths);
            stats.nts(trial, snr) = calc_stats(est.nts(trial, snr), signal, par, paths);
            
            azi_RMSE_ts(trial, snr) = stats.ts(trial, snr).azi_rmse_deg;
            ele_RMSE_ts(trial, snr) = stats.ts(trial, snr).ele_rmse_deg;
            RMSAE_ts(trial, snr) = stats.ts(trial, snr).rmsae_deg;
            azi_RMSE_nts(trial, snr) = stats.nts(trial, snr).azi_rmse_deg;
            ele_RMSE_nts(trial, snr) = stats.nts(trial, snr).ele_rmse_deg;
            RMSAE_nts(trial, snr) = stats.nts(trial, snr).rmsae_deg;
            
            azi_norecom_RMSE_ts(trial, snr) = stats.ts(trial,snr).norecom_azi_rmse_deg;
            ele_norecom_RMSE_ts(trial, snr) = stats.ts(trial,snr).norecom_ele_rmse_deg;
            norecom_RMSAE_ts(trial, snr) = stats.ts(trial, snr).norecom_rmsae_deg;
            azi_norecom_RMSE_nts(trial,snr) = stats.nts(trial,snr).norecom_azi_rmse_deg;
            ele_norecom_RMSE_nts(trial,snr) = stats.nts(trial,snr).norecom_ele_rmse_deg;
            norecom_RMSAE_nts(trial, snr) = stats.nts(trial, snr).norecom_rmsae_deg;
            
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
    azi_rmse_ts_ave = mean(azi_RMSE_ts);
    ele_rmse_ts_ave = mean(ele_RMSE_ts);
    azi_rmse_nts_ave = mean(azi_RMSE_nts);
    ele_rmse_nts_ave = mean(ele_RMSE_nts);
    azi_norecom_rmse_ts_ave = mean(azi_norecom_RMSE_ts);
    ele_norecom_rmse_ts_ave = mean(ele_norecom_RMSE_ts);
    azi_norecom_rmse_nts_ave = mean(azi_norecom_RMSE_nts);
    ele_norecom_rmse_nts_ave = mean(ele_norecom_RMSE_nts);
    rmsae_ts_ave = mean(RMSAE_ts);
    rmsae_nts_ave = mean(RMSAE_nts);
    rmsae_norecom_ts_ave = mean(norecom_RMSAE_ts);
    rmsae_norecom_nts_ave = mean(norecom_RMSAE_nts);
    
    figure('Name','angular_rmse');
    semilogy(par.snrSweep, rmsae_ts_ave); hold on; grid on;
    semilogy(par.snrSweep, rmsae_nts_ave);
    title('Temporal Smoothing: Angular RMSE')
    legend("Temporal Smoothing, "+par.blocks+" Blocks", 'No Smoothing')
    ylim([10^-2 10^3])
    xlabel('SNR')
    ylabel('RMSE, Degrees')
    
    figure('Name','norecom_angular_rmse');
    semilogy(par.snrSweep, rmsae_norecom_ts_ave); hold on; grid on;
    semilogy(par.snrSweep, rmsae_norecom_nts_ave);
    title('Temporal Smoothing: Angular RMSE without recombinations')
    legend("Temporal Smoothing, "+par.blocks+" Blocks", 'No Smoothing')
    ylim([10^-2 10^3])
    xlabel('SNR')
    ylabel('RMSE, Degrees')
    
    figure('Name','norecom_azi_rmse');
    semilogy(par.snrSweep, azi_norecom_rmse_ts_ave); hold on; grid on;
    semilogy(par.snrSweep, azi_norecom_rmse_nts_ave);
    title('Temporal Smoothing: Azimuth RMSE without recombinations')
    legend("Temporal Smoothing, "+par.blocks+" Blocks", 'No Smoothing')
    ylim([10^-2 10^3])
    xlabel('SNR')
    ylabel('RMSE, Degrees')
    
    figure('Name','norecom_ele_rmse');
    semilogy(par.snrSweep, ele_norecom_rmse_ts_ave); hold on; grid on;
    semilogy(par.snrSweep, ele_norecom_rmse_nts_ave);
    title('Temporal Smoothing: Elevation RMSE without recombinations')
    legend("Temporal Smoothing, "+par.blocks+" Blocks", 'No Smoothing')
    ylim([10^-2 10^3])
    xlabel('SNR')
    ylabel('RMSE, Degrees')
    
    figure('Name','azi_fig_rmse');
    semilogy(par.snrSweep, azi_rmse_ts_ave); hold on; grid on;
    semilogy(par.snrSweep, azi_rmse_nts_ave);
    title('Temporal Smoothing: Azimuth RMSE')
    legend("Temporal Smoothing, "+par.blocks+" Blocks", 'No Smoothing')
    ylim([10^-2 10^3])
    xlabel('SNR')
    ylabel('RMSE, Degrees')
    
    figure('Name','ele_fig_rmse');
    semilogy(par.snrSweep, ele_rmse_ts_ave); hold on; grid on;
    semilogy(par.snrSweep, ele_rmse_nts_ave);
    title('Temporal Smoothing: Elevation RMSE')
    legend("Temporal Smoothing, "+par.blocks+" Blocks", 'No Smoothing')
    ylim([10^-2 10^3])
    xlabel('SNR')
    ylabel('RMSE, Degrees')
    
elseif par.type == '1d'
    for snr=1:length(par.snrSweep)
        par.SNR = par.snrSweep(snr);
        for trial = 1:par.Trials
            paths = assign_paths(par);
            signal = transmitter(paths, par);
            signal = reciever(signal, paths, par);
            signal = temporal_smooth(signal, par);
            
            est.ts(trial, snr) = est_drmusic(signal.ts, par);
            est.nts(trial, snr) = est_drmusic(signal.R, par);
            
            stats.ts(trial, snr) = calc_stats(est.ts(trial, snr), signal, par, paths);
            stats.nts(trial, snr) = calc_stats(est.nts(trial, snr), signal, par, paths);
            
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
    legend("Temporal Smoothing, "+par.blocks+" Blocks", 'No Smoothing')
    ylim([10^-2 10^3])
    xlabel('SNR')
    ylabel('Corrected RMSE')
    
    figure('Name','fig_rmse');
    semilogy(par.snrSweep, rmse_ts_ave); hold on; grid on;
    semilogy(par.snrSweep, rmse_nts_ave);
    title('Temporal Smoothing vs. No Smoothing')
    legend("Temporal Smoothing, "+par.blocks+" Blocks", 'No Smoothing')
    ylim([10^-2 10^3])
    xlabel('SNR')
    ylabel('RMSE')
    
    figure('Name','fig_norecom_rmse');
    semilogy(par.snrSweep, norecom_rmse_ts_ave); hold on; grid on;
    semilogy(par.snrSweep, norecom_rmse_nts_ave);
    title('Temporal Smoothing vs. No Smoothing : No Recombinations')
    legend("Temporal Smoothing, "+par.blocks+" Blocks", 'No Smoothing')
    ylim([10^-2 10^3])
    xlabel('SNR')
    ylabel('RMSE')
end

end

