% Subject:  Cyclostationarity Signal Generator
% Date:     May 23, 2020
% Name:     Daniel Tait

% Function that creates some cyclostationary signal.
clear all; close all; clc; format compact;

%% Input

Ns = 10*1024;
fc = 0;
Tsym = 64;

%% Function

% Create Pulse Shape
a = 0;
Lpulse = 10;
pulseType = 'rc';
[P,delay] = PulseShape(Tsym,a,Lpulse,pulseType);
Nsp = Tsym;

% Random Symbols
Nbuf = 5;
Nsym = floor(Ns/Nsp);
Nsymtot = Nsym + Nbuf;
k = 2;
modScheme = 'psk';
ys = sqrt(Tsym)*DigitalModulation(Nsymtot,k,modScheme);

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

%% Plot Setup

pltPosition = [0.1 0.1 0.35 0.3; 
                0.1 0.5 0.35 0.3;
                0.5 0.1 0.35 0.3
                0.5 0.5 0.35 0.3];
pltx = 1;
pltPosition(:,1) = pltPosition(:,1) + pltx;
%[NumFigures,~] = size(pltPosition);
NumFigures = 3;

for i=1:NumFigures
    hf = figure(i);
    clf;
    set(hf, 'units', 'normalized');
    set(hf, 'paperpositionmode', 'auto');
    set(hf, 'position', pltPosition(i,:));
end

%% Output

ny = 0:length(y)-1;

figure(1)
plot(ny,abs(y))

figure(2)
plot(ny,real(y))
hold on
plot(ny,imag(y))

figure(3)
Y = fftshift(fft(y));
plot(abs(Y))

energy = y*y'