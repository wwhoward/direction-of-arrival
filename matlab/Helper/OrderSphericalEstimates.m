function [phiDeg_ord,thtDeg_ord] = OrderSphericalEstimates(phiDeg,thtDeg,phiDeg_h,thtDeg_h)

N = length(phiDeg_h);

phi = deg2rad(phiDeg-180);
tht = deg2rad(thtDeg-90);
phi_h = deg2rad(phiDeg_h-180);
tht_h = deg2rad(thtDeg_h-90);
[x,y,z] = sph2cart(phi,tht,ones(1,N));
[xh,yh,zh] = sph2cart(phi_h,tht_h,ones(1,N));

u = [x;y;z];
uh = [xh;yh;zh];
uhh = OrderDistance(u,uh);

[phi_ord,tht_ord] = cart2sph(uhh(1,:),uhh(2,:),uhh(3,:));

phiDeg_ord = rad2deg(phi_ord)+180;
thtDeg_ord = rad2deg(tht_ord)+90;