function [uh,fh,R] = EstimatorESPRIT(y1,y2,Nd,K)

Z1 = y1;
Z2 = y2;

ph = [];

N = max(size(Z1));
R = [Z1*Z1', Z1*Z2';...
     Z2*Z1', Z2*Z2']/N;

% Signal Subspace Matrix
[E,Ee] = eig(R);
E1 = E(1:6,:);
E2 = E(7:12,:);
E1 = E1(:,12-K+1:end);
E2 = E2(:,12-K+1:end);

% Total Least Squares
[U,S,V] = svd([E1,E2]'*[E1,E2]);
V12 = U(1:K,K+1:end);
V22 = U(K+1:end,K+1:end);
Psi = -V12*(V22\eye(length(V22)));

% Least Squares
% Psi = inv(E1'*E1)*(E1'*E2);

[Tinv,Phi] = eig(Psi);
fh = angle(diag(Phi))/Nd/2/pi;
Ah = 0.5*(E1*Tinv + E2*Tinv*inv(Phi));
% Ah = 0.5*(E1*inv(T) + E2*inv(T*Phi'));

uh = zeros(3,K);
for k=1:K
    ek = Ah(1:3,k).';
    hk = Ah(4:6,k).';
    uh(:,k) = normalize(real(cross(ek,conj(hk))),'norm');
end
