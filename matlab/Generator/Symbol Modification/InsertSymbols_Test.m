% Subject:  Insert Pilot Function
% Date:     June 6, 2020
% Name:     Daniel Tait

% Inserts pilot symbols into a matrix of symbols

%%
clear all; close all; clc; format compact;

% Generate symbols
Nsym = 32;
k = log2(16);
modScheme = 'qam';
s = DigitalModulation(Nsym,k,modScheme);
M = 2^k;
pilotIdx = randi([1,M]);
p = GetSymbol(pilotIdx,k,modScheme);

% OFDM
Nsc = 4;
s = reshape(s,Nsc,[]);

% comb type pilots => adjacent pilots in time
% block type pilots => adjacent pilots in frequency
% null symbol => zero

% Parameters
signatureEnable = 0;
signatureIndex = [1:Nsc];

preambleEnable = 0;
preambleLength = 3;
preambleIndex = [3,4,5];

pilotEnable = 1;
pilotType = 'block';
pilotSpacing = 4;

% Cyclostationary signature
if signatureEnable
    signatureIndex(signatureIndex > Nsc) = [];
    for i = signatureIndex
        s(i,:) = s(signatureIndex(1),:);
    end
end

% Preamble
if preambleEnable
    p = DigitalModulation(preambleLength,k,modScheme);
    preambleIndex(preambleIndex > Nsc) = [];
    s(:,1:preambleLength) = 0;
    for i = preambleIndex
        s(i,1:preambleLength) = p;
    end
end

% Pilots
if pilotEnable
    [Nsc,Nsymsc] = size(s);
    if preambleEnable
        tmp = preambleLength + 1;
    else
        tmp = 1;
    end
    switch pilotType
        case 'block'
            p = DigitalModulation(Nsc,k,modScheme).';
            for i = tmp:pilotSpacing:Nsymsc
                s(:,i) = p;
            end
        case 'comb'
            p = DigitalModulation(Nsymsc,k,modScheme);
            for i = tmp:pilotSpacing:Nsc
                s(i,:) = p;
            end
        otherwise
            error("Invalid pilot type.");
    end
end

s









