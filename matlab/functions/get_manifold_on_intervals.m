% Temporal smoothing (TS) trials for the vector sensor
% Summer 2020
% Will Howard, {wwhoward}@vt.edu
% wireless @ VT

% 'sim' 1d & 2d done, need to do 8in next. 

function [DOA_function, DOA_function_H, AOA, AOA_H, phi, tht] = get_manifold_on_intervals(par)

ele_endpoints = [par.eleEstRange(1):par.fusion_interval:(par.eleEstRange(2)-par.fusion_interval) ; par.eleEstRange(1)+par.fusion_interval:par.fusion_interval:par.eleEstRange(2)]';
azi_endpoints = [par.aziEstRange(1):par.fusion_interval:(par.aziEstRange(2)-par.fusion_interval) ; par.aziEstRange(1)+par.fusion_interval:par.fusion_interval:par.aziEstRange(2)]';


switch par.manifold
    case 'sim'
        switch par.type
            case '1d'
                phi = zeros(length(azi_endpoints), par.fusion_res);
                
                DOA_function = zeros(6, 2, par.fusion_res, length(azi_endpoints));
                DOA_function_H = zeros(2, 6, size(phi,2), length(azi_endpoints));
                AOA = zeros(6, size(phi,2), length(azi_endpoints));
                AOA_H = zeros(1, 6, size(phi,2), length(azi_endpoints));
                for m=1:length(azi_endpoints)
                    phi(m,:) = azi_endpoints(m,1):1/par.fusion_res:azi_endpoints(m,2);
                    tht = pi/2;
                    
%                     tmp_DOA_function = zeros(6, 2, size(phi,2));
%                     tmp_DOA_function_H = zeros(2, 6, size(phi,2));
%                     tmp_AOA = zeros(6, size(phi,2));
%                     tmp_AOA_H = zeros(1, 6, size(phi,2));
                    
                    for ph = 1:size(phi,2)
                        [DOA_function(:,:,ph,m),~,AOA(:,ph,m)] = VectorSensor([pi/2,phi(m,ph)],[pi/4,0]);
                        DOA_function_H(:,:,ph,m) = squeeze(DOA_function(:,:,ph,m))';
                        AOA_H(:,:,ph,m) = squeeze(AOA(:,ph,m))';
                    end
                end
            case '2d'
                phi = zeros(length(azi_endpoints), par.fusion_res);
                tht = zeros(length(ele_endpoints), par.fusion_res);
                
                DOA_function = zeros(6, 2, size(phi,2), size(tht,2), length(azi_endpoints), length(ele_endpoints));
                DOA_function_H = zeros(2, 6, size(phi,2), size(tht,2), length(azi_endpoints), length(ele_endpoints));
                AOA = zeros(6, size(phi,2), size(tht,2), length(azi_endpoints), length(ele_endpoints));
                AOA_H = zeros(1, 6, size(phi,2), size(tht,2), length(azi_endpoints), length(ele_endpoints));
                
                for m=1:length(azi_endpoints)
                    for n=1:length(ele_endpoints)
                        phi(m,:) = azi_endpoints(m,1):1/par.fusion_res:azi_endpoints(m,2);
                        tht(n,:) = ele_endpoints(n,1):1/par.fusion_res:ele_endpoints(n,2);
                        
                        
                        for ph = 1:size(phi,2)
                            for th = 1:size(tht,2)
                                [DOA_function(:,:,ph,th,m,n),~,AOA(:,ph,th,m,n)] = VectorSensor([tht(n,th), phi(m,ph)], [pi/4, 0]);
                                DOA_function_H(:,:,ph,th,m,n) = squeeze(DOA_function(:,:,ph,th,m,n))';
                                AOA_H(1, :,ph,th) = squeeze(AOA(:,ph,th,m,n))';
                            end
                        end
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
                phi = par.aziEstRange(1):1/par.fusion_res:par.aziEstRange(2);
                tht = pi/2;
                
                AOA = zeros(6, size(phi,2));
                for i=1:6
                    AOA(i,:) = spline(0:pi/180:(2*pi-pi/180), squeeze(manifold(10, :, i)), phi);
                end
                AOA_H = conj(reshape(AOA, 1, 6, size(phi,2)));
            case '2d'
                phi = par.aziEstRange(1):1/par.fusion_res:par.aziEstRange(2);
                tht = par.eleEstRange(1):1/par.fusion_res:par.eleEstRange(2);
                
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

