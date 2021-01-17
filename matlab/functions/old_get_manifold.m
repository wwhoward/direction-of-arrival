% Temporal smoothing (TS) trials for the vector sensor
% Summer 2020
% Will Howard, {wwhoward}@vt.edu
% wireless @ VT

function [DOA_function, DOA_function_H, phi, tht] = getManifold(par)

ele_endpoints = [0:par.fusion_interval:(pi-par.fusion_interval) ; par.fusion_interval:par.fusion_interval:pi]';
azi_endpoints = [0:par.fusion_interval:((2*pi)-par.fusion_interval) ; par.fusion_interval:par.fusion_interval:(2*pi)]';

switch par.type
    case '1d'
        phi = 0:1/par.res:2*pi;
        tht = pi/2;
        
        for i=1:2*pi/par.fusion_interval
            fusion_phi(i,:) = azi_endpoints(i,1):1/par.fusion_res:azi_endpoints(i,2);
        end
        fusion_ele = 
        
        DOA_function = zeros(6, 2, size(phi,2));
        DOA_function_intervals = zeros(par.fusion_interval, 6, 2, size(fusion_phi,2)/par.fusion_interval);
        
        DOA_function_H = zeros(2, 6, size(phi,2));
        DOA_function_interval_H = zeros(par.fusion_interval, 6, 2, size(fusion_phi,2)/par.fusion_interval);
        
        for ph = 1:size(phi,2)
            [DOA_function(:,:,ph),~,~] = VectorSensor([pi/2,phi(ph)],[pi/4,0]);
            DOA_function_H(:,:,ph) = DOA_function(:,:,ph)';
        end
    case '2d'
        phi = 0:1/par.res:2*pi;
        tht = 0:1/par.res:pi;
        
        DOA_function = zeros(6, 2, size(phi,2), size(tht,2));
        DOA_function_H = zeros(2, 6, size(phi,2), size(tht,2));
        for ph = 1:size(phi,2)
            for th = 1:size(tht,2)
                [DOA_function(:,:,ph,th),~,~] = VectorSensor([tht(th), phi(ph)], [pi/4, 0]);
                DOA_function_H(:,:,ph,th) = DOA_function(:,:,ph,th)';
            end
        end
end
end

