% Temporal smoothing (TS) trials for the vector sensor
% Summer 2020
% Will Howard, {wwhoward}@vt.edu
% wireless @ VT

function [AOA, AOA_H, phi, tht, DOA_function, DOA_function_H] = get_manifold(par, res, range_azi, range_ele)

if nargin == 1
    range_azi = par.aziEstRange;
    range_ele = par.eleEstRange;
    res = par.res;
end

switch par.manifold
    case 'sim'
        switch par.type
            case '1d'
                phi = range_azi(1):1/res:range_azi(2);
                tht = pi/2;

                DOA_function = zeros(6, 2, size(phi,2));
                DOA_function_H = zeros(2, 6, size(phi,2));
                AOA = zeros(6, size(phi,2));
                AOA_H = zeros(1, 6, size(phi,2));
                for ph = 1:size(phi,2)
                    [DOA_function(:,:,ph),~,AOA(:,ph)] = VectorSensor([pi/2,phi(ph)],[pi/4,0]);
                    DOA_function_H(:,:,ph) = DOA_function(:,:,ph)';
                    AOA_H(:,:,ph) = AOA(:,ph)';
                end
            case '2d'
                phi = range_azi(1):1/res:range_azi(2);
                tht = range_ele(1):1/res:range_ele(2);

                DOA_function = zeros(6, 2, size(phi,2), size(tht,2));
                DOA_function_H = zeros(2, 6, size(phi,2), size(tht,2));
                AOA = zeros(6, size(phi,2), size(tht,2));
                AOA_H = zeros(1, 6, size(phi,2), size(tht,2));
                for ph = 1:size(phi,2)
                    for th = 1:size(tht,2)
                        [DOA_function(:,:,ph,th),~,AOA(:,ph,th)] = VectorSensor([tht(th), phi(ph)], [pi/4, 0]);
                        DOA_function_H(:,:,ph,th) = DOA_function(:,:,ph,th)';
                        AOA_H(1, :,ph,th) = AOA(:,ph,th)';
                    end
                end
        end
    case '8in'
        load('Helper/8inSim/mat/v_pol_max.mat', 'arr');
        %manifold = permute(arr, [2,1,3]);
        manifold = arr;
        DOA_function = 0;
        DOA_function_H = 0;
        switch par.type
            case '1d'
                phi = range_azi(1):1/res:range_azi(2);
                tht = pi/2;
                
                AOA = zeros(6, size(phi,2));
                for i=1:6
                    AOA(i,:) = spline(0:pi/180:(2*pi-pi/180), squeeze(manifold(10, :, i)), phi);
                end
                AOA_H = conj(reshape(AOA, 1, 6, size(phi,2)));
            case '2d'
                phi = range_azi(1):1/res:range_azi(2);
                tht = range_ele(1):1/res:range_ele(2);
                
                AOA = zeros(6, size(phi,2), size(tht,2));
                for i=1:6
                    AOA(i,:,:) = interp2(0:pi/180:(2*pi-pi/180), pi/4:5*pi/180:10*pi/18, squeeze(manifold(:,:,i)), phi,tht', 'spline')';
                end
                AOA_H = conj(reshape(AOA, 1, 6, size(phi,2), size(tht,2))); %Might not be the right way to do this
        end
%         tmpAOA = AOA-min(AOA(:));
%         AOA = tmpAOA./max(tmpAOA(:));
%         
%         tmpAOA_H = AOA_H-min(AOA_H(:));
%         AOA_H = tmpAOA_H./max(tmpAOA_H(:));
end
end

