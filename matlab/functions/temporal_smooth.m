% Temporal smoothing (TS) trials for the vector sensor
% Summer 2020
% Will Howard, {wwhoward}@vt.edu
% wireless @ VT

function signal = temporal_smooth(signal, par)

R_ts = zeros(6,6);
r = zeros(par.blocks, 6,6);
for b=1:par.blocks
    r(b,:,:) = 1/par.snapshot * squeeze(signal.rx(b,:,:))*squeeze(signal.rx(b,:,:))'; % Calculate covariance for this block, weighted by the length of the window
    R_ts(:,:) = 1/par.blocks * squeeze(r(b,:,:)) + R_ts(:,:); % TS covariance is the weighted sum of block covariances
end
signal.ts = R_ts;
end

