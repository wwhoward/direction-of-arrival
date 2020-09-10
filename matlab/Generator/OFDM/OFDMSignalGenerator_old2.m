function y = OFDMSignalGenerator(Ns,fc,cyclicParameter)

if nargin < 3
    cyclicParameter = 1;
end

% OFDM
Nsc = 16;
W = 0.5;
Tsym = 64;

% Digital Modulation
Nsym = Nsc*Ns/Tsym;
k = log2(16);
modScheme = 'qam';
s = DigitalModulation(Nsym,k,modScheme);

% Cyclic Prefix and Zero Padding
guardLength = 0.0;
Tg = floor(guardLength*Tsym);
Tsc = Tsym - Tg;
Nsymsc = Nsym/Nsc;

% Cyclostationary Signature
s = reshape(s,Nsc,[]);

% if cyclicParameter == 1
%     for k=1:Nsc
%         s(k,:) = s(1,:);
%     end
% %     s(Nsc,:) = s(1,:);
% end

% Method 1
y = repmat(reshape(s.',1,[]),[Tsc,1]);
y = reshape(y,[],Nsc).';

% Method 2
% y = zeros(Nsc,Nsymsc*Tsc);
% idx = 1:Tsc:Nsymsc*Tsc;
% y(:,idx) = s;

% Subcarrier Modulation
n = 0:Tsc-1;
fsc = 0:1/Tsc:(Nsc-1)/Tsc;
Dn = exp(1j*2*pi*fsc'*n);
Dn = repmat(Dn,[1,Nsymsc]);
y = y.*Dn;

% Prefix Insertion
guardType = 'cp';
xt = InsertCyclicPrefix(y,Tg,Tsym,guardType);

% Carrier Frequency
y = sum(xt,1);
n = 0:length(y)-1;
y = y.*exp(1j*2*pi*(fc-Nsc/2/Tsc)*n);

y = sqrt(Ns)*y/sqrt(y*y');
