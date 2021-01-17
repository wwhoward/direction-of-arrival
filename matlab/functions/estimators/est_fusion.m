function est = est_fusion(R, par, paths, signal)
% 1) Estimate model order
% 2) Estimate indicator vector
% 3) Estimate MUSIC in indicated intervals

%%
%-----------------//-----------------
% Model Order Estimation goes here
k_hat = par.K;
%-----------------//-----------------

%%
%-----------------//-----------------
% Indicator Vector Estimation goes here, requires estimate of model order
% todo: Replace all of this with model based estimation!!!
est_prelim = est_drmusic(R, par);
stats = calc_stats(est_prelim, signal, par, paths);

ele_endpoints = [0:par.fusion_interval:(pi-par.fusion_interval) ; par.fusion_interval:par.fusion_interval:pi]';
azi_endpoints = [0:par.fusion_interval:((2*pi)-par.fusion_interval) ; par.fusion_interval:par.fusion_interval:(2*pi)]';

ele_ind_est = sum(((stats.uh(2,:)>ele_endpoints(:,1)) & (stats.uh(2,:) < ele_endpoints(:,2))),2);
azi_ind_est = sum(((stats.uh(1,:)>azi_endpoints(:,1)) & (stats.uh(1,:) < azi_endpoints(:,2))),2);
%-----------------//-----------------

%%
%-----------------//-----------------
% MUSIC on intervals & association

%-----------------//-----------------
switch par.type
    case '1d'
        azi_ind_est = paths.azi_intind;
        
        azi_index = find(azi_ind_est == 1);
        
        azi_interval = [par.fusion_interval*(azi_index-1), par.fusion_interval*azi_index];
    case '2d'
        
        azi_index = find(azi_ind_est == 1);
        ele_index = find(ele_ind_est == 1);
        
        azi_interval = [par.fusion_interval*(azi_index-1), par.fusion_interval*azi_index];
        ele_interval = [par.fusion_interval*(ele_index-1), par.fusion_interval*ele_index];
end

azi_index = find(azi_ind_est == 1);
ele_index = find(ele_ind_est == 1);
azi_n = length(azi_index);
ele_n = length(ele_index);

azi_interval = [par.fusion_interval*(azi_index-1), par.fusion_interval*azi_index];
ele_interval = [par.fusion_interval*(ele_index-1), par.fusion_interval*ele_index];

weights = zeros(azi_n, ele_n);
idx = zeros(2,azi_n, ele_n);

for azi=1:azi_n
    for ele=1:ele_n
        [AOA, AOA_H, phi, tht] = get_manifold(par, par.fusion_res, azi_interval(azi,:), ele_interval(ele,:));
        spectrum = est_music(signal.ts, par, azi_interval(azi,:), ele_interval(ele,:), AOA, AOA_H);
        weights(azi,ele) = max(spectrum.spectrum(:));
        [x,y] = find(spectrum.spectrum==weights(azi,ele));
        idx(:,azi, ele) = [phi(x),tht(y)];
    end
end

% Use hungarian to find the matching

assignment = munkres(-1*weights); %negative to maximize cost
est = idx(:,assignment);
%-----------------//-----------------
% 

end

