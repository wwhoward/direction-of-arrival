function [P,t,filtDelay] = PulseRaisedCosine(a,L,Nsym)

Tsym = 1;
t = -(Nsym/2):1/L:(Nsym/2);
A = sin(pi*t/Tsym)./(pi*t/Tsym);
B = cos(pi*a*t/Tsym);
P = A.*B./(1-(2*a*t/Tsym).^2);
P(ceil(length(P)/2)) = 1;
temp = (a/2)*sin(pi/(2*a));
P(t==Tsym/(2*a)) = temp;
P(t==-Tsym/(2*a)) = temp;
P = P/sqrt(P*P');
filtDelay = (length(P)-1)/2;