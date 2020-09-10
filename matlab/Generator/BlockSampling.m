% Generate a number of rx signals with multipath present. Specify coherence
% time, number of paths
% 
% Approximates coherence time by seperating into blocks of fixed path gains
% wwhoward, wireless@vt

clc,clear,close all

function [rx] = GenMultipath(n_paths, t_coh, n_samples, path_delay)

if nargin == 3
    path_delay = zeros(1, n_paths);
elseif nargin == 4
    if length(path_delay) ~= n_paths
        error('path_delay isn''t long enough')
    end
end

y = DigitalSignalGenerator(Ns,fc,Tsym,a,Lpulse,par)




end