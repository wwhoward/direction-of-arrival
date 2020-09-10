function [gtht,gpol, g] = VectorSensor(aoa,pol)

% aoa - [elevation,azimuth]
%   elevation on [0,pi)
%   azimuth on [0,2pi)
% pol - [gamma,eta]
%   gamma (aux pol phase angle) on [0,pi/2)
%   eta (pol phase difference) on [-pi,pi)
%   eta = 0deg => linear pol
%   gam = 45deg => circular pol, eta = +/-90deg => LCP/RCP


gtht = [cos(aoa(1))*cos(aoa(2)),    -sin(aoa(2));
        cos(aoa(1))*sin(aoa(2)),    cos(aoa(2));
        -sin(aoa(1)),               0;
        -sin(aoa(2)),               -cos(aoa(1))*cos(aoa(2));
        cos(aoa(2)),                -cos(aoa(1))*sin(aoa(2));
        0,                          sin(aoa(1))];
    
gpol = [sin(pol(1))*exp(1j*pol(2));
        cos(pol(1))];
    
g = gtht * gpol;
