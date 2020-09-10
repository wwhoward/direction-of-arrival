% Subject:  OFDM Function 2
% Date:     June 5, 2020
% Name:     Daniel Tait

% Creates an OFDM signal.

%%
% Creating another OFDM function that works better, so there can be a window.
clear all; %close all; clc; format compact;

% Generate symbols
Nsym = 32*1024;
k = log2(16);
modScheme = 'qam';
s = DigitalModulation(Nsym,k,modScheme);

% OFDM parameters
Nsc = 16;
Tsym = 128;
b = 0.1;                % Window parameter
guardType = 'cp';       % Guard {'cp','zp'}
guardLength = 0.25;
Tg = floor(guardLength*Tsym);
Tsc = Tsym - Tg;
Nsymsc = Nsym/Nsc;

% Symbols
s = reshape(s,Nsc,[]);
so = s;
so(:,2:2:end) = 0;
tmp = circshift(so,1,2);
tmp(:,1) = 0;
so = so + tmp;
se = s;
se(:,1:2:end) = 0;
tmp = circshift(se,-1,2);
tmp(:,end) = 0;
se = se + tmp;

% Reshaping for carrier modulation
ye = repmat(reshape(se.',1,[]),[Tsym,1]);
ye = reshape(ye,[],Nsc).';
yo = repmat(reshape(so.',1,[]),[Tsym,1]);
yo = reshape(yo,[],Nsc).';

% Centering symbols over the windows
yec = circshift(ye,Tsym/2,2);
yo(:,1:Tsym/2) = repmat(yo(:,end),[1,Tsym/2]);
yoc = circshift(yo,-Tsym/2,2);

% Creating window
w = RaisedCosineWindow(2*Tsym,b);
w = w/sqrt(w*w');
wo = circshift(w,-Tsym/2);
wo = repmat(wo,[1,round(Nsymsc/2)]);
we = circshift(w,Tsym/2);
we = repmat(we,[1,round(Nsymsc/2)]);
wo(end-Tsym/2:end) = 0;
we(1:Tsym/2) = 0;
wo = wo(1:length(yo));
we = we(1:length(ye));

% Carrier modulation
nsc = 0:Tsc-1;
nright = (nsc(end)+1):(nsc(end)+Tsym/2);
if Tg > 0
    ng = (Tsc-Tg):Tsc-1;
    nleft = (ng(1)-Tsym/2):(ng(1)-1);
    n = [nleft,ng,nsc,nright];
else
    ng = [];
    nleft = (nsc(1)-Tsym/2):(nsc(1)-1);
    n = [nleft,nsc,nright];
end
fo = 0.1;
fsc = 0:1/Tsc:(Nsc-1)/Tsc;
fsc = fsc + fo;

% Guard insertion
switch guardType
    case 'cp'
        D1 = exp(1j*2*pi*fsc'*n);
    case 'zp'
        D1 = [zeros(Nsc,length([nleft,ng])),exp(1j*2*pi*fsc'*[nsc,nright])];
    otherwise
        error("Invalid guard type.");
end
D = repmat(D1,[1,round(Nsymsc/2)]);
D = D(:,1:length(yo));
De = circshift(D,Tsym/2,2);
Do = circshift(D,-Tsym/2,2);

% Combine even and odd symbols
xe = ye.*De.*we;
xo = yo.*Do.*wo;
xt = xo + xe;

% Scaling the subcarriers
p = ones(Nsc,1);
p(4:5) = 2;
xt = xt.*sqrt(p);

% Summing subcarriers and shifting
fc = 0;
y = sum(xt,1);
t = 0:length(y)-1;
y = y.*exp(1j*2*pi*(fc-Nsc/2/Tsc-fo)*t);

% figure(3)
% plot(we,'-')
% hold on
% plot(abs(yec(1,:)),'-')
% title('even window')
% 
% figure(4)
% plot(wo,'-')
% hold on
% plot(abs(yoc(1,:)),'-')
% title('odd window')

figure(1)
for i=1:Nsc
    %xt = Tcp*Dn'*so;
    z = xt(i,:);
    L = 1024;
    [Sz,ff] = cpsd(z,z,L,round(L/4),[],1,'centered');
    plot(ff,Sz);
    hold on
end

%figure(2)
L = 1024;
[Sx,ff] = cpsd(y,y,L,round(L/4),[],1,'centered');
semilogy(ff,Sx);
%ylim([1e-5,1])
hold on

% figure
% plot(angle(D1(3,:)))
% figure
% plot(real(D1(3,:)))
% hold on
% plot(imag(D1(3,:)))
% 
% figure
% plot(real(xe(3,:)))

%%
% clear all; close all; clc; format compact;
% 
% Ns = 8*1024;
% fc = 0;
% y = OFDMSignalGenerator(Ns,fc);
% 
% figure(1)
% L = 1024;
% [Sx,ff] = cpsd(y,y,L,round(L/4),[],1,'centered');
% semilogy(ff,Sx);
% hold on