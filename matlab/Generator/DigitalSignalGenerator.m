% Subject:  Digital Signal Generator Function
% Date:     May 23, 2020
% Name:     Daniel Tait

% Improved the original "SignalGenerator" function to include different pulse shaping and modulation.

function y = DigitalSignalGenerator(Ns,fc,Tsym,a,Lpulse,par)

if nargin < 6
    signalType = 'sc'; %{'mc','sc'}
    pilotEnable = 0;
    pilotSpacing = 4;
    fd = 0;
    pulseType = 'rc';
else
    signalType = par.signalType;
    pilotEnable = par.pilotEnable;
    if pilotEnable
        pilotSpacing = par.pilotSpacing;
    end
    fd = par.fd;
    pulseType = par.pulseType; %{'rc','rrc','sqr'}
end

% Create Pulse Shape
%a = 1;
%Lpulse = 16;
[P,delay] = PulseShape(Tsym,a,Lpulse,pulseType);
Nsp = Tsym;

% Random Symbols
Nbuf = 16;
Nsym = floor(Ns/Nsp);
Nsymtot = Nsym + Nbuf;
% k = log2(16);
% modScheme = 'qam';
if isempty(par.M)
    k = log2(4);
else
    k = log2(par.M);
end
if isempty(par.modScheme)
    modScheme = 'psk';
else
    modScheme = par.modScheme;
end

switch signalType
    case 'sc'
        ys = sqrt(Tsym)*DigitalModulation(Nsymtot,k,modScheme);
        
        % Pilots
        if pilotEnable
            ys(1:pilotSpacing:end) = ys(1);
        end
        
        % Pulse Shaping
        if rem(Nsp,2)
            Nsep = Nsp-2;
        else
            Nsep = Nsp-1;
        end
        y = [ys.',zeros(Nsymtot,Nsep)];
        y = reshape(y.',1,[]);
        y = conv(P,y);
        y = y(delay:end-floor(Nsep/2)-floor(Nsp/2));

        % Bandpass modulation
        Ny = length(y);
        n = 0:Ny-1;
        y = y.*exp(1j*2*pi*fc*n);
        offset = randi([1,Nsp]);
        y = y(offset:offset+Ns-1);
        phase = exp(1j*2*pi*rand);
        %y = phase*y;
        
    case 'mc'
        yss = sqrt(Tsym)*DigitalModulation(2*Nsymtot,k,modScheme);
        yss = reshape(yss,2,[]);
        if pilotEnable
            %yss(2,1:par.pilotSpacing:end) = yss(1,1:par.pilotSpacing:end);
            yss(1,1:pilotSpacing:end) = yss(1);
            yss(2,1:pilotSpacing:end) = yss(1);
        end
        for i=1:2
            ys = yss(i,:);
            % Pulse Shaping
            if rem(Nsp,2)
                Nsep = Nsp-2;
            else
                Nsep = Nsp-1;
            end
            z = [ys.',zeros(Nsymtot,Nsep)];
            z = reshape(z.',1,[]);
            z = conv(P,z);
            z = z(delay:end-floor(Nsep/2)-floor(Nsp/2));

            % Bandpass modulation
            if i == 1
                offset = randi([1,Nsp]);
                Ny = length(z);
                n = 0:Ny-1;
                carrier = exp(1j*2*pi*(fc-fd/2)*n);
            end
            if i == 2
                carrier = exp(1j*2*pi*(fc+fd/2)*n);
            end
            %zz(i,:) = z;
            z = z.*carrier;
            z = z(offset:offset+Ns-1);
            phase = exp(1j*2*pi*rand);
            %y = phase*y;
            
            if i == 1
                y = zeros(1,length(z));
            end
            y = y + z/sqrt(2);
        end
end











