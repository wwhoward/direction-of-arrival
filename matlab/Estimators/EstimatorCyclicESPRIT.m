function uh = EstimatorCyclicESPRIT(Rc1,Rc2,K)

% Function
[E1,lambda1] = eig(Rc1);
[~,idx] = sort(abs(diag(lambda1)),'descend');
E1 = E1(:,idx(1:K));
% lambda1
%E1 = E1(:,6-K+1:end);

[E2,lambda2] = eig(Rc2);
[~,idx] = sort(abs(diag(lambda2)),'descend');
E2 = E2(:,idx(1:K));
% lambda2
%E2 = E2(:,6-K+1:end);

% Total Least Squares
[U,S,V] = svd([E1,E2]'*[E1,E2]);
V12 = U(1:K,K+1:end);
V22 = U(K+1:end,K+1:end);
PSI = -V12*(V22\eye(length(V22)));

[T,PHI] = eig(PSI);
%fh = angle(diag(PHI))/Nd/2/pi;
Ah = 0.5*(E1*T + E2*T*conj(PHI));

uh = zeros(3,K);
for k=1:K
    ek = Ah(1:3,k).';
    hk = Ah(4:6,k).';
    uh(:,k) = normalize(real(cross(ek,conj(hk))),'norm');
end

% Bad result
% R1 = Rcy(:,:,2,2);
% R2 = Rcy(:,:,3,3) + Rcy(:,:,1,3);

% Good Result
% R1 = Rcy(:,:,3,1);
% R2 = Rcy(:,:,3,3);

% Good Result
% R1 = Rcy(:,:,1,1)+Rcy(:,:,2,1)+Rcy(:,:,3,1);
% R2 = Rcy(:,:,5,3)+Rcy(:,:,6,3)+Rcy(:,:,7,3);

% R1 = Rcy(:,:,1,1)+Rcy(:,:,3,1); % was using this before the change
% R2 = Rcy(:,:,5,3)+Rcy(:,:,7,3);

% Good Result
% R1 = Rcy(:,:,3,2);
% R2 = Rcy(:,:,3,3);