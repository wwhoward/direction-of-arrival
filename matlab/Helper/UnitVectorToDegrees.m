function [phiDeg,thtDeg] = UnitVectorToDegrees(uh)

[~,K] = size(uh);
phiDeg = zeros(1,K);
thtDeg = zeros(1,K);
for k=1:K
    uhk = uh(:,k);
    [u_azi,u_ele,~] = cart2sph(uhk(1),uhk(2),uhk(3));
    phiDeg(k) = rad2deg(u_azi);
    if phiDeg(k) < 0
        phiDeg(k) = phiDeg(k) + 360;
    end
    thtDeg(k) = -rad2deg(u_ele) + 90;
end

end