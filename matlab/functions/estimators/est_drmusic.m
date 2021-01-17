% Temporal smoothing (TS) trials for the vector sensor
% Summer 2020
% Will Howard, {wwhoward}@vt.edu
% wireless @ VT

function est = est_drmusic(R, par, azi_interval, ele_interval)
% Input:
%   R: 6x6 Covariance Matrix

% if par.resType == 'measured'
%     ele_crb = 10^(-par.SNR/10) / (2*par.blocks*par.snapshot);
%     res = 
% end

if nargin == 2
    azi_interval = [0,2*pi];
    ele_interval = [0,pi];
end


switch par.type
    case '1d'
        phi = par.phi;
        DOA_function = par.DOA_function;
        DOA_function_H = par.DOA_function_H;
        
        [eigvect, eigval]=eig(R);
        [~,idx] = sort(diag(eigval));
        NoiseSpace = eigvect(:,idx(1:length(idx)-par.K));
        Noise = NoiseSpace*NoiseSpace';
        
        if par.accelerate
            est.spectrum = 1./min(real(eig2(mtimesx(mtimesx(DOA_function_H, Noise), DOA_function))));
        else
            for i=1:size(DOA_function, 3)
                est.spectrum(i) = 1/min(real(eig((DOA_function_H(:,:,i)*Noise*DOA_function(:,:,i))))); % This should be real since B should be hermitian, but due to computational error, there will be some small imaginary part. Take abs() or real() to correct
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
        DOA_function_H = par.DOA_function_H;
        
        [eigvect, eigval]=eig(R);
        [~, idx] = sort(diag(eigval));
        NoiseSpace = eigvect(:, idx(1:length(idx)-par.K));
        Noise = (NoiseSpace*NoiseSpace');
        
        if par.accelerate
            est.spectrum = reshape(1./min(real(eig2(mtimesx(mtimesx(DOA_function_H, Noise), DOA_function)))), size(DOA_function,3), size(DOA_function,4));
        else
            est.spectrum = zeros(size(DOA_function,3), size(DOA_function,4));
            for ph=1:size(DOA_function, 3)
                for th=1:size(DOA_function,4)
                    est.spectrum(ph,th) = 1/min(real(eig((DOA_function_H(:,:,ph,th)*Noise*DOA_function(:,:,ph,th)))));
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

