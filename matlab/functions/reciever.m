% Temporal smoothing (TS) trials for the vector sensor
% Summer 2020
% Will Howard, {wwhoward}@vt.edu
% wireless @ VT

function signal = reciever(signal, paths, par)
% Take modulated message, simulate multipath propegation, add AWGN

tx = signal.tx;
rx = zeros(par.blocks, 6, par.snapshot);
for b=1:par.blocks
    for k=1:paths.sources
        for p=1:paths.multi(k)
            if paths.multi ~= 1
                h = sqrt(0.5)*(randn+1j*randn); % Random complex path gain if multipath is present
            else
                h = 1; % Unit path gain if multipath is not present
            end
            a = paths.signal_vector(k).path(p,:)'; % Array manifold for desired AoA of this source&path
            del = randi(par.pathdelay); % Path delay, random in specified range
            block = randi(par.interblock); % Block delay, random in specified range
            sig = tx(k, 1+(b-1)*block+(p-1)*del : par.snapshot + (b-1)*block + (p-1)*del); % Window the source appropriately, given desired delays
            rx(b, :,:) = par.powerDecay^p * h*a*sig + squeeze(rx(b,:,:)); % Sum cumulative signal with current source, multiplying by array manifold & path gain
        end
    end
end
signal.rx = awgn(rx, par.SNR, 'measured'); % Add AWGN at specified SNR

signal.R = 1/par.snapshot * squeeze(signal.rx(1,:,:))*squeeze(signal.rx(1,:,:))'; % Calculate covariance matrix for the first block

end
