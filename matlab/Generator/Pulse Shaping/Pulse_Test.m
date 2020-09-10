% Subject:  Pulse Shaping Test
% Date:     May 23, 2020
% Name:     Daniel Tait

% Testing the pulse-shaping functions.
clear all; close all; clc; format compact;

%% Input

Tsym = 1;
L = 128;
Nsym = 2;
alphas = [0, 0.5, 1];
fs = L/Tsym;
lineColors = ['b','r','g','k','c'];
i=1;
legendString = cell(1,4);

for a=alphas

    %[PTx,t,delay] = PulseRootRaisedCosine(a,L,Nsym);
    [PTx,t] = PulseSquare(L);
    %[PTx,t] = PulseRaisedCosine(a,L,Nsym);
    PRx = PTx;
    %y = conv(PTx,PRx,'same');
    y = PTx;
    
    subplot(1,2,1);
    t = Tsym*t;
    plot(t,y,lineColors(i));
    hold on;

    subplot(1,2,2);
    Y = fftshift(fft(y));
    plot(abs(Y),lineColors(i));
    hold on;
    
    i = i+1;
    
end



