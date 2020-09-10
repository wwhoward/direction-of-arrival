function [uh,fh] = EstimatorESPRITCorrelationMatrix(R,delta,K,par)

bool = 0;
M = length(R)/2;
if nargin == 4
    bool = 1;
end

% Signal Subspace Matrix
[E,Ee] = eig(R);
[~,idx] = sort(abs(diag(Ee)),'descend');
E = E(:,idx);
E1 = E(1:M,:);
E2 = E(M+1:2*M,:);
Es1 = E1(:,1:K);
Es2 = E2(:,1:K);

% Total Least Squares
[U,S,V] = svd([Es1,Es2]'*[Es1,Es2]);
V12 = U(1:K,K+1:end);
V22 = U(K+1:end,K+1:end);
Psi = -V12*(V22\eye(length(V22)));

[T,Phi] = eig(Psi);
fh = angle(diag(Phi))/delta/2/pi;
Ah = (Es1*T + Es2*T*conj(Phi))/2;
% Ah = (Es1*T + Es2*T*conj(Phi))/2;
% Ah = Es1*T;
% Ah = Ah(1:6,:);
% if bool
% %     for k=1:K
% %         Es1o(:,k) = mean(reshape(Es1(:,k),M,[]),2);
% %         Es2o(:,k) = mean(reshape(Es2(:,k),M,[]),2);
% %     end
% %     Ah = (Es1o*T + Es2o*T*conj(Phi));
% 
%     for k=1:K
%         Aho(:,k) = mean(reshape(Ah(:,k),M,[]),2);
%     end
%     Ah = Aho;
% end
Ah = ( (abs(Ah)).^(1/2) ).*exp(1j*angle(Ah));

uh = zeros(3,K);
for k=1:K
    ek = Ah(1:3,k).';
    hk = Ah(4:6,k).';
    uh(:,k) = normalize(real(cross(ek,conj(hk))),'norm');
end
