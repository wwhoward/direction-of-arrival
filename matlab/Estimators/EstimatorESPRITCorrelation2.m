function [uh,fh,ph] = EstimatorESPRITCorrelation2(R,K)

% Signal Subspace Matrix
[E,Ee] = eig(R);
% E1 = E(1:6,:);
% E2 = E(7:12,:);
% E1 = E1(:,12-K+1:end);
% E2 = E2(:,12-K+1:end);

% Total Least Squares
[U,S,V] = svd([E]'*[E]);
V12 = U(1:K,K+1:end);
V22 = U(K+1:end,K+1:end);
PSI = -V12*(V22\eye(length(V22)));

[T,PHI] = eig(PSI);
Ah = 0.5*(E1*T + E2*T*conj(PHI));

uh = zeros(3,K);
for k=1:K
    ek = Ah(1:3,k).';
    hk = Ah(4:6,k).';
    uh(:,k) = normalize(real(cross(ek,conj(hk))),'norm');
end
