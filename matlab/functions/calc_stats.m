% Temporal smoothing (TS) trials for the vector sensor
% Summer 2020
% Will Howard, {wwhoward}@vt.edu
% wireless @ VT

function stats = calc_stats(est, signal, par, paths)
switch par.type
    case '1d' % need to separate into azi err, ele err
        u = paths.AoA(:,2)';
        uh = est.peaks_azi;
        
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
                    closestValue(i) = pi*rand;
                    stats.recom_flag(i) = 1;
                elseif par.recombination == 'rnd'
                    stats.err(i) = pi*rand;
                    stats.corrected_err(i) = stats.err(i);
                    closestValue(i) = pi*rand;
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
        
    case '2d'
        azi_u = paths.AoA(:,2)';
        ele_u = paths.AoA(:,1)';
        %azi_uh = zeros(size(azi_u));
        %ele_uh = zeros(size(azi_u));
        azi_uh = est.peaks_azi;
        ele_uh = est.peaks_ele;
        stats.azi_uh = azi_uh;
        stats.ele_uh = ele_uh;
        stats.peak_val = est.peak_val;
        
        stats.recom_flag = zeros(size(azi_u));
        if  isempty(azi_uh) % Sloppy, don't assign any angles if no detections, just set error to rand and recom to 1. 
            azi_uh = rand(size(azi_u))*pi;
            ele_uh = rand(size(ele_u))*pi/2;
            stats.recom_flag = ones(size(azi_u));
        end
        
        % Assign estimates to truth
        
        dist_matrix = abs(acos(sin(ele_u') * sin(ele_uh) + cos(ele_u)' * cos(ele_uh).* cos(min(abs(azi_u'-azi_uh), abs(azi_u'-azi_uh - pi)))));
        
        peak_idx = zeros(size(azi_u));
        azi_uh_ordered = zeros(size(azi_u));
        ele_uh_ordered = zeros(size(ele_u));
        stats.recom_flag = zeros(size(azi_u));
        err = zeros(size(azi_u));
        unused_tru = 1:length(azi_u);
        unused_est = 1:length(azi_uh);
        new_dist = dist_matrix;
        
        % My bipartite 1:1 matching algo, bad compared to Hungarian (munkres) so replaced, but keep for the memories
        % todo Only use max peaks for matching
        for p = 1:min(length(azi_uh), length(azi_u))
            [tru_idx, est_idx] = find(dist_matrix == min(new_dist(:)));
            tru_idx = min(intersect(tru_idx, unused_tru));
            est_idx = min(intersect(est_idx, unused_est));
            peak_idx(tru_idx) = est_idx;
            err(tru_idx) = dist_matrix(tru_idx, est_idx);
            azi_uh_ordered(tru_idx) = azi_uh(est_idx);
            ele_uh_ordered(tru_idx) = ele_uh(est_idx);
            unused_tru = setdiff(unused_tru, tru_idx);
            unused_est = setdiff(unused_est, est_idx);
            new_dist = dist_matrix(unused_tru, unused_est);
        end
        while ~isempty(find(err==0, 1))
            bad_idx = find(err==0);
            azi_uh_ordered(bad_idx) = rand(size(bad_idx))*pi;
            ele_uh_ordered(bad_idx) = rand(size(bad_idx))*pi/2;
            err(bad_idx) = rand(size(bad_idx))*pi;
            stats.recom_flag(bad_idx) = 1;
        end
        
        %Do the same as above, but using the top K peaks & the Hungarian algo
        
%         [~, peak_order] = sort(stats.peak_val);
%         dist_matrix_topK = abs(acos(sin(ele_u(peak_order(1:min(length(ele_u),par.K)))') * sin(ele_uh(peak_order(1:min(length(ele_uh),par.K)))) + cos(ele_u(peak_order(1:min(length(ele_u),par.K))))' * cos(ele_uh(peak_order(1:min(length(ele_uh),par.K)))).* cos(min(abs(azi_u(peak_order(1:min(length(ele_uh),par.K)))'-azi_uh(peak_order(1:min(length(ele_uh),par.K)))), abs(azi_u(peak_order(1:min(length(ele_uh),par.K)))'-azi_uh(peak_order(1:min(length(ele_uh),par.K))) - pi)))));
%         dist_matrix_topK = abs(acos(sin(ele_uh(peak_order(1:min(length(ele_uh),par.K)))
%         assignment = munkres(dist_matrix);
        
        stats.peak_idx = peak_idx;
        
        stats.angular_err = err;
        stats.azi_err = min(abs(azi_u - azi_uh_ordered), abs(azi_u - azi_uh_ordered - pi));
        stats.ele_err = min(abs(ele_u - ele_uh_ordered), abs(ele_u - ele_uh_ordered - pi/2));
        
        
        stats.u = [azi_u; ele_u];
        stats.uh = [azi_uh_ordered; ele_uh_ordered];
        
        stats.azi_mse = sum(stats.azi_err.^2)/length(stats.azi_err);
        stats.ele_mse = sum(stats.ele_err.^2)/length(stats.ele_err);
        stats.azi_rmse = sqrt(stats.azi_mse);
        stats.ele_rmse = sqrt(stats.ele_mse);
        stats.msae = sum(stats.angular_err.^2)/length(stats.angular_err);
        stats.rmsae = sqrt(stats.msae);
        
        stats.azi_err_deg = stats.azi_err/pi * 180;
        stats.ele_err_deg = stats.ele_err/pi * 180;
        stats.azi_mse_deg = sum(stats.azi_err_deg.^2)/length(stats.azi_err_deg);
        stats.ele_mse_deg = sum(stats.ele_err_deg.^2)/length(stats.ele_err_deg);
        stats.azi_rmse_deg = sqrt(stats.azi_mse_deg);
        stats.ele_rmse_deg = sqrt(stats.ele_mse_deg);
        stats.angular_err_deg = stats.angular_err/pi * 180;
        stats.msae_deg = sum(stats.angular_err_deg.^2)/length(stats.angular_err_deg);
        stats.rmsae_deg = sqrt(stats.msae_deg);
        
        
        % Now calculate these for only the non-recom cases
        stats.recom_percent = sum(stats.recom_flag)/length(stats.recom_flag);
        stats.norecom_azi_mse = sum(stats.azi_err(stats.recom_flag==0).^2)/length(stats.azi_err(stats.recom_flag==0));
        stats.norecom_ele_mse = sum(stats.ele_err(stats.recom_flag==0).^2)/length(stats.ele_err(stats.recom_flag==0));
        stats.norecom_azi_rmse = sqrt(stats.norecom_azi_mse);
        stats.norecom_ele_rmse = sqrt(stats.norecom_ele_mse);
        stats.norecom_msae = sum(stats.angular_err(stats.recom_flag==0).^2)/length(stats.angular_err(stats.recom_flag==0));
        stats.norecom_rmsae = sqrt(stats.norecom_msae);
        
        stats.norecom_azi_err_deg = stats.azi_err(stats.recom_flag==0)/pi * 180;
        stats.norecom_ele_err_deg = stats.ele_err(stats.recom_flag==0)/pi * 180;
        stats.norecom_azi_mse_deg = sum(stats.norecom_azi_err_deg.^2)/length(stats.norecom_azi_err_deg);
        stats.norecom_ele_mse_deg = sum(stats.norecom_ele_err_deg.^2)/length(stats.norecom_ele_err_deg);
        stats.norecom_azi_rmse_deg = sqrt(stats.norecom_azi_mse_deg);
        stats.norecom_ele_rmse_deg = sqrt(stats.norecom_ele_mse_deg);
        stats.norecom_angular_err_deg = stats.angular_err(stats.recom_flag==0)/pi * 180;
        stats.norecom_msae_deg = sum(stats.norecom_angular_err_deg.^2)/length(stats.norecom_angular_err_deg);
        stats.norecom_rmsae_deg = sqrt(stats.norecom_msae_deg);
        
        stats.SNR = par.SNR;
        
end
end
