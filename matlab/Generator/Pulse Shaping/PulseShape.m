function [P,filtDelay] = PulseShape(Ns,a,Nspan,pulseType)

if nargin <= 3
    pulseType = 'rrc';
end

% Input
L = floor(Ns/2);    %oversampling factor
Nsym = 2*Nspan;     %filter span

switch pulseType
    case 'sqr'
        [P,t,filtDelay] = PulseSquare(L);
    case 'rc'
        [P,t,filtDelay] = PulseRaisedCosine(a,L,Nsym);
    case 'rrc'
        [P,t,filtDelay] = PulseRootRaisedCosine(a,L,Nsym);
    otherwise
        error("Invalid pulse type.");
end