% Temporal smoothing (TS) trials for the vector sensor
% Summer 2020
% Will Howard, {wwhoward}@vt.edu
% wireless @ VT

function signal = transmitter(paths, par)
% Create messages and modulate them
% Number of bits transmitted per frame is set to be 1000. For QPSK
% modulation, this corresponds to 500 symbols per frame.
if nargin < 2
    par.mod = 'QPSK';
end

switch par.mod
    case 'QPSK'
        modu = comm.QPSKModulator( ...
            'BitInput',    true, ...
            'PhaseOffset', pi/4);
end

for k=1:paths.sources
    msg(k, :) = randi([0, 1], par.signal_length, 1);
    tx(k, :) = modu(msg(k,:)');
end


signal.tx = tx;
end


