function [P,t,filtDelay] = PulseSquare(L)

t = -1:1/L:1;
P = ones(1,length(t));
P(1) = 0;
P = P/sqrt(P*P');
filtDelay = (length(P)-1)/2;