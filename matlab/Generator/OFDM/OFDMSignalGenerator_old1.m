function y = OFDMSignalGeneratorOrthogonality(Ns,fc,tempPar)

% First OFDM signal Generator script, before changes were made
% Should have a difference in the frequency spacing

if nargin <= 2
    tempPar = 1;
end

K = 8;     % Upsampling factor

% Digital Modulation
Nsym = floor(Ns/K);
k = log2(4);
modScheme = 'psk';
s = DigitalModulation(Nsym,k,modScheme);

% OFDM
Nsc = 64;

% DFT Matrix
Nc = Nsc;
nr = 1:1/K:Nc;
Nr = length(nr);
nr = repmat(nr.',[1,Nc]);
nc = 1:Nc;
nc = repmat(nc,[Nr,1]);
Dn = exp(-1j*2*pi*(nr-1).*(nc-1)/(Nsc-1)).'/sqrt(Nsc);

% Cyclic Prefix
Tcp = 0.2;
In = eye(Nr);
Ncp = floor(Tcp*Nr);
Incp = In(:,end-Ncp+1:end);
Tcp = [Incp.';In.'];

% Zero Padding
% Tzp = 0.1;
% In = eye(Nr);
% Ncp = floor(Tzp*Nr);
% Incp = In(:,end-Ncp+1:end);
% Tcp = [0*Incp.';In.'];

s = reshape(s,Nsc,[]);

% Cyclostationary Signature
%s(end,:) = s(1,:);
%s(end-1,:) = s(2,:);

% if tempPar == 1
%     s(end,:) = s(1,:);
%     s(2,:) = s(1,:);
%     s(end-1,:) = s(1,:);
% end

if tempPar == 1
    for k=1:Nsc
        s(k,:) = s(1,:);
    end
end


% s(Nsc/2+1,:) = s(Nsc/2,:);
% s(Nsc/2-1,:) = s(Nsc/2,:);
% s(Nsc/2+2,:) = s(Nsc/2,:);

% s(Nsc/2-2,:) = s(1,:);
% s(Nsc/2+1,:) = s(1,:);
% s(Nsc/2-1,:) = s(1,:);
% s(Nsc/2+2,:) = s(1,:);

xt = Tcp*Dn'*s;
%xt = Dn'*s;
y = reshape(xt,1,[]);
n = 0:length(y)-1;
y = y.*exp(1j*2*pi*(fc-1/K/2)*n);

y = y(1:Ns);
