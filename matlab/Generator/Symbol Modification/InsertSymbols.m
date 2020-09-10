function y = InsertSymbols(s,k,modScheme,par)

[Nsc,N] = size(s);

M = 2^k;
% pilotIdx = randi([1,M]);
% p = GetSymbol(pilotIdx,k,modScheme);

% % Parameters
% signatureEnable = 0;
% signatureIndex = [1,2,Nsc-1,Nsc];
% 
% preambleEnable = 0;
% preambleLength = 64;
% preambleIndex = 1:Nsc;
% 
% pilotEnable = 1;
% pilotType = 'block';
% % pilotType = 'comb';
% pilotSpacing = 1;

% Parameters
signatureEnable = par.signatureEnable;
signatureIndex = par.signatureIndex;
signatureIndexShift = par.signatureIndexShift;
signatureIndexShiftValue = par.signatureIndexShiftValue;

preambleEnable = par.preambleEnable;
preambleLength = par.preambleLength;
preambleIndex = par.preambleIndex;

pilotEnable = par.pilotEnable;
pilotType = par.pilotType;
pilotSpacing = par.pilotSpacing;

% Cyclostationary signature
if signatureEnable
    signatureIndex(signatureIndex > Nsc) = [];
    for i = signatureIndex
        s(i,:) = s(signatureIndex(1),:);
    end
    if ~isempty(signatureIndexShift)
        for i = signatureIndexShift
            s(i,signatureIndexShiftValue:end) = s(signatureIndex(1),1:end-signatureIndexShiftValue+1);
        end
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
            if par.pilotRepeatSymbol
                p(2:end) = p(1);
            end
            for i = tmp:pilotSpacing:Nsymsc
                s(:,i) = p;
            end
        case 'comb'
            p = DigitalModulation(Nsymsc,k,modScheme);
            if par.pilotRepeatSymbol
                p(2:end) = p(1);
            end
            if isempty(par.pilotCarrierIndex)
                for i = 1:pilotSpacing:Nsc
                    s(i,:) = p;
                end
            else
                for i = par.pilotCarrierIndex
                    s(i,:) = p;
                end
            end
        otherwise
            error("Invalid pilot type.");
    end
end

y = s;
