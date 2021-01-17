% Temporal smoothing (TS) trials for the vector sensor
% Summer 2020
% Will Howard, {wwhoward}@vt.edu
% wireless @ VT

function est = temp_est_fusion(R, par, paths)
% Input:
%   R: 6x6 Covariance Matrix

% if par.resType == 'measured'
%     ele_crb = 10^(-par.SNR/10) / (2*par.blocks*par.snapshot);
%     res = 
% end






detection_indicator = [0 0 1 0 0 1 0 0];
detection_index = find(detection_indicator==1);
num_intervals = length(detection_indicator); % For consistancy the real number is this minus one
points_in_interval = 2^4;

for ind = 1:length(detection_index)
    detection = detection_index(ind);
    phi=linspace(360*(detection-1)/num_intervals,360*detection/num_intervals,points_in_interval);
    phi=phi(1:end-1);
end
 



switch par.type
    case '1d'
        phi = par.phi;
        DOA_function = par.DOA_function;
        DOA_function_inverse = par.DOA_function_inverse;
        
        [eigvect, eigval]=eig(R);
        [~,idx] = sort(diag(eigval));
        NoiseSpace = eigvect(:,idx(1:length(idx)-par.K));
        Noise = NoiseSpace*NoiseSpace';
        
        if par.accelerate
            est.spectrum = 1./min(real(eig2(mtimesx(mtimesx(DOA_function_inverse, Noise), DOA_function))));
        else
            for i=1:size(DOA_function, 3)
                est.spectrum(i) = 1/min(real(eig((DOA_function_inverse(:,:,i)*Noise*DOA_function(:,:,i))))); % This should be real since B should be hermitian, but due to computational error, there will be some small imaginary part. Take abs() or real() to correct
            end
        end
            
        [peakvals, peak_idx] = findpeaks(est.spectrum);
        [est.peak_val, idx] = sort(peakvals, 'descend');
        peak_idx = peak_idx(idx);
        est.peaks_azi = phi(peak_idx);    % azimuth indices of peaks
        
    case '2d'
        phi = par.phi;
        tht = par.tht;
        DOA_function = par.DOA_function;
        DOA_function_inverse = par.DOA_function_inverse;
        
        [eigvect, eigval]=eig(R);
        [~, idx] = sort(diag(eigval));
        NoiseSpace = eigvect(:, idx(1:length(idx)-par.K));
        Noise = (NoiseSpace*NoiseSpace');
        
        if par.accelerate
            est.spectrum = reshape(1./min(real(eig2(mtimesx(mtimesx(DOA_function_inverse, Noise), DOA_function)))), size(DOA_function,3), size(DOA_function,4));
        else
            est.spectrum = zeros(size(DOA_function,3), size(DOA_function,4));
            for ph=1:size(DOA_function, 3)
                for th=1:size(DOA_function,4)
                    est.spectrum(ph,th) = 1/min(real(eig((DOA_function_inverse(:,:,ph,th)*Noise*DOA_function(:,:,ph,th)))));
                end
            end
        end
        
        peaks_idx = FastPeakFind(est.spectrum);
        if isempty(peaks_idx)
            [Xmax, Imax] = extrema2(est.spectrum);
            [x_idx, y_idx] = ind2sub(size(est.spectrum),Imax(1:min(10, length(Imax))));
            peaks_idx = [y_idx';x_idx'];
            peaks_idx = peaks_idx(:)';
        end
        azi_peaks = peaks_idx(2:2:end);
        ele_peaks = peaks_idx(1:2:end);
        [est.peak_val, sorted_idx] = sort(diag(est.spectrum(azi_peaks, ele_peaks)), 'descend');
        est.peaks_azi = phi(azi_peaks(sorted_idx));
        est.peaks_ele = tht(ele_peaks(sorted_idx));
        
        
        
end
end

