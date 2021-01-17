% Temporal smoothing (TS) trials for the vector sensor
% Summer 2020
% Will Howard, {wwhoward}@vt.edu
% wireless @ VT

function [azi_, ele_, pol1_, pol2_] = CalcCramerRao(t1, t2, t3, t4, noisePower, totalWindow)
% Output is in units of rad^2
t2 = t2 - pi/2; % This is because while in our simulation, elevation is in [0,pi], the CRLB requires elevation in [-pi/2, pi/2].

d = 2*(cos(t2)^2 + 2*sin(t4)^2*sin(t2)^2)*sin(t3)^2*cos(t3)^2 + (4*cos(t3)^2 -1)*sin(t2)^2*sin(t4)^2;

azi_ = noisePower / totalWindow * (sin(t3)^2 * cos(t3)^2) / d;

ele_ = noisePower / (2*totalWindow);

pol1_ = noisePower / (4*totalWindow) * (2*(sin(t2)^2+1)*sin(t3)^2*cos(t3)^2 + (4*cos(t3)^2 - 1)*sin(t2)^2*sin(t4)^2)/d;

pol2_ = noisePower / (2*totalWindow) * (cos(t2)^2 + 8*sin(t2)^2 * cos(t3)^2 * sin(t4)^2) / d;

end