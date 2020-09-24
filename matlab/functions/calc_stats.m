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
        if isempty(est.peaks_azi)
            azi_uh = rand(size(azi_u))*pi;
            ele_uh = rand(size(azi_u))*pi;
        else
            azi_uh = est.peaks_azi;
            ele_uh = est.peaks_ele;
        end
        azi_uh_ordered = zeros(size(azi_u));
        ele_uh_ordered = zeros(size(ele_u));
        
        % Assign estimates to truth
        dist_matrix = abs(acos(sin(ele_u') * sin(ele_uh) + cos(ele_u') * cos(ele_uh).* cos(azi_u'-azi_uh)));
        
        [err, peak_idx] = min(dist_matrix, [], 2);
        peak_idx = reshape(peak_idx, size(azi_u));
        azi_uh_ordered = azi_uh(peak_idx);
        ele_uh_ordered = ele_uh(peak_idx);
        stats.recom_flag = zeros(size(azi_u));
        
        if length(peak_idx) ~= length(unique(peak_idx))
            [~, w] = unique(peak_idx, 'stable');
            dup = max(setdiff(1:numel(peak_idx), w));
            dup_val = peak_idx(dup);
            dup_idx = find(peak_idx==peak_idx(dup));
            good_idx = find(err==min(err(dup_idx)));
            bad_idx = setdiff(dup_idx, good_idx);
            if length(azi_u) == length(azi_uh)
                peak_idx(bad_idx) = setdiff(1:numel(peak_idx), unique(peak_idx));
                err(bad_idx) = diag(dist_matrix(bad_idx, peak_idx(bad_idx)));
                azi_uh_ordered = azi_uh(peak_idx);
                ele_uh_ordered = ele_uh(peak_idx);
            else
                peak_idx(dup_idx)=0;
                peak_idx(good_idx) = dup_val;
                err(bad_idx) = rand*pi;
                azi_uh_ordered(bad_idx) = rand(1, length(bad_idx))*pi;
                ele_uh_ordered(bad_idx) = rand(1, length(bad_idx))*pi/2;
                stats.recom_flag(1,[bad_idx]) = 1;
            end
        end
        
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
