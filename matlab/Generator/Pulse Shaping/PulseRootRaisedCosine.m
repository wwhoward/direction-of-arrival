function [P,t,filtDelay] = PulseRootRaisedCosine(a,L,Nsym)

Tsym = 1;
t = -(Nsym/2):1/L:(Nsym/2);
num = sin(pi*t*(1-a)/Tsym)+...
    ((4*a*t/Tsym).*cos(pi*t*(1+a)/Tsym));
den = pi*t.*(1-(4*a*t/Tsym).^2)/Tsym;
P = 1/sqrt(Tsym)*num./den;
P(ceil(length(P)/2)) = 1/sqrt(Tsym)*((1-a)+4*a/pi);
temp = (a/sqrt(2*Tsym))*((1+2/pi)*sin(pi/(4*a))...
    +(1-2/pi)*cos(pi/(4*a)));
P(t==Tsym/(4*a)) = temp;
P(t==-Tsym/(4*a)) = temp;
P = P/sqrt(P*P');
filtDelay = (length(P)-1)/2;