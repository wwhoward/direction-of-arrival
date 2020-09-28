% Temporal smoothing (TS) trials for the vector sensor
% Summer 2020
% Will Howard, {wwhoward}@vt.edu
% wireless @ VT

function [DOA_function, DOA_function_inverse, phi, tht] = getManifold(par)

switch par.type
    case '1d'
        phi = 0:1/par.res:2*pi;
        tht = pi/2;
        
        DOA_function = zeros(6, 2, size(phi,2));
        DOA_function_inverse = zeros(2, 6, size(phi,2));
        for ph = 1:size(phi,2)
            [DOA_function(:,:,ph),~,~] = VectorSensor([pi/2,phi(ph)],[pi/4,0]);
            DOA_function_inverse(:,:,ph) = DOA_function(:,:,ph)';
        end
    case '2d'
        phi = 0:1/par.res:2*pi;
        tht = 0:1/par.res:pi;
        
        DOA_function = zeros(6, 2, size(phi,2), size(tht,2));
        DOA_function_inverse = zeros(2, 6, size(phi,2), size(tht,2));
        for ph = 1:size(phi,2)
            for th = 1:size(tht,2)
                [DOA_function(:,:,ph,th),~,~] = VectorSensor([tht(th), phi(ph)], [pi/4, 0]);
                DOA_function_inverse(:,:,ph,th) = DOA_function(:,:,ph,th)';
            end
        end
end
end

