function [u_h_ord,B2] = OrderDistance(u,u_h)

Nsig = size(u,2);
Au = repmat(u,[1,1,Nsig]);
Au_h = repmat(u_h,[1,1,Nsig]);
for n=1:Nsig-1
    Au_h(:,:,n+1) = circshift(Au_h(:,:,n+1),[0,n,0]);
end

B = permute(acos(dot(Au,Au_h,1)),[2,3,1]);
B2 = B;
colIdx = repmat(1:Nsig,[1,2]);
u_h_ord = zeros([3,Nsig]);
IN = eye(Nsig);
for n=1:Nsig
    [~,Ic] = min(B(:));
    [Ir,Ic] = ind2sub(size(B),Ic);
    idx = colIdx(Nsig+1+Ir-Ic);
    B = B + 100*circshift(IN,idx-1);
    B(Ir,:) = 100;
    u_h_ord(:,Ir) = u_h(:,idx);
end