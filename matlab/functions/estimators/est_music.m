% Temporal smoothing (TS) trials for the vector sensor
% Summer 2020
% Will Howard, {wwhoward}@vt.edu
% wireless @ VT

function est = est_music(R, par, azi_interval, ele_interval, manifold, manifold_inverse)
% Input:
%   R: 6x6 Covariance Matrix

% if par.resType == 'measured'
%     ele_crb = 10^(-par.SNR/10) / (2*par.blocks*par.snapshot);
%     res =
% end

if nargin == 2
    azi_interval = [0, 2*pi];
    ele_interval = [0,pi];
    manifold = par.AOA_function;
    manifold_inverse = par.AOA_function_H;
end


switch par.type
    case '1d'
        phi = par.phi;
        AOA_function = manifold;
        AOA_function_H = manifold_inverse;
        
        [eigvect, eigval]=eig(R);
        [~,idx] = sort(diag(eigval));
        NoiseSpace = eigvect(:,idx(1:length(idx)-par.K));
        Noise = NoiseSpace*NoiseSpace';
        
        if par.accelerate && 0==1
            est.spectrum = squeeze(1./abs(mtimesx(mtimesx(AOA_function_H, Noise), reshape(AOA_function,6,1,size(AOA_function,2)))))';
        else
            for i=1:size(AOA_function, 2)
                % est.spectrum(i) = squeeze(1/abs(AOA_function_H(:,:,i)*Noise*AOA_function(:,i))); % This should be real since B should be hermitian, but due to computational error, there will be some small imaginary part. Take abs() or real() to correct
                est.spectrum(i) = 1/abs(AOA_function_H(:,:,i)*Noise*AOA_function(:,i)); % This should be real since B should be hermitian, but due to computational error, there will be some small imaginary part. Take abs() or real() to correct
            end
        end
        
        [peakvals, peak_idx] = findpeaks(est.spectrum);
        [est.peak_val, idx] = sort(peakvals, 'descend');
        peak_idx = peak_idx(idx);
        est.peaks_azi = phi(peak_idx);    % azimuth indices of peaks
        
    case '2d'
        phi = par.phi(par.phi >= azi_interval(1) & par.phi <= azi_interval(2));
        tht = par.tht(par.tht >= ele_interval(1) & par.tht <=ele_interval(2));
        AOA_function = manifold;%(:,par.phi >= azi_interval(1) & par.phi <= azi_interval(2), par.tht >= ele_interval(1) & par.tht <=ele_interval(2));
        AOA_function_H = manifold_inverse;%(:,:,par.phi >= azi_interval(1) & par.phi <= azi_interval(2), par.tht >= ele_interval(1) & par.tht <=ele_interval(2));
        
        [eigvect, eigval]=eig(R);
        [~, idx] = sort(diag(eigval));
        NoiseSpace = eigvect(:, idx(1:length(idx)-par.K));
        Noise = (NoiseSpace*NoiseSpace');
        
        if par.accelerate
            est.spectrum = squeeze(1./abs(mtimesx(mtimesx(AOA_function_H, Noise), reshape(AOA_function,6,1,size(AOA_function,2), size(AOA_function,3)))));
        else % This step needs updating to MUSIC rather than DR-MUSIC
            est.spectrum = zeros(size(AOA_function,3), size(AOA_function,4));
            for ph=1:size(DOA_function, 3)
                for th=1:size(DOA_function,4)
                    est.spectrum(ph,th) = 1/min(real(eig((AOA_function_H(:,:,ph,th)*Noise*AOA_function(:,:,ph,th)))));
                end
            end
        end
        if nargin == 2
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
end

