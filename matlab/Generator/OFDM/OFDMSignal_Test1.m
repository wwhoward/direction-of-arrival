% Subject:  OFDM Function
% Date:     May 24, 2020
% Name:     Daniel Tait

% Creates an OFDM signal.

%%
% Testing frequency orthogonality.
clear all; close all; clc; format compact;

N = 1024;
T = 100;
f1 = 1/T;
f2 = 2/T;
n = 1:N-1;
x1 = exp(1j*2*pi*rand)*exp(1j*2*pi*f1*n);
x2 = exp(1j*2*pi*rand)*exp(1j*2*pi*f2*n);

temp = x1(1:T).*conj(x2(1:T));
E = sum(temp)/T


%%
% OFDM in the comms book that I own.
clear all; close all; clc; format compact;

% Digital Modulation
Nsym = 10*1024;
k = log2(16);
modScheme = 'qam';
s = DigitalModulation(Nsym,k,modScheme);

% OFDM Parameters
Nsc = 64;
Ncp = 16;

xt = ifft(s,Nsc);
xt = [xt(end-Ncp+1:end),xt];

Y = fftshift(fft(y));

plot(abs(Y))

%%
% OFDM test.
clear all; close all; clc; format compact;

% Digital Modulation
Nsym = 128*1024;
k = log2(16);
modScheme = 'qam';
s = DigitalModulation(Nsym,k,modScheme);

% OFDM
k = 8;
Nsc = 64;
Nf = k*Nsc;
Ncp = k*16;
In = eye(Nf);
Incp = In(:,end-Ncp+1:end);
Tcp = [Incp.';In.'];                        %add cyclic prefix
Dn = CreateDFTMatrix([Nf,Nf],Nf)/sqrt(Nsc);

so = [s;zeros(k-1,Nsym)];
so = reshape(so,Nf,[]);
%so = reshape([s,zeros(Nf-Nsc),Nf,[]);
% so = reshape([so'; zeros(size(so'))],[],1);
% so = reshape(so,Nf,[]);
xt = Tcp*Dn'*so;
y = reshape(xt,1,[]);

Y = fftshift(fft(y));
plot(abs(Y))
% L = 256; 
% [Sx,ff] = cpsd(x,x,L,round(L/4),[],1,'centered');
% plot(n,x1,n,x2);

%%
clear all; close all; clc; format compact;
Nr = 8;
Nc = 4;
k = 4;
Nf = k*Nr;
Ns = Nr*Nc;
x = randn(1,Ns);
x = [x;zeros(k-1,Ns)]
y = reshape(x,Nf,[])

%y = permute(reshape(x,Nf,[]),[2,1])
imagesc(abs(y))
figure
imagesc(abs(x))

%%
% Creating OFDM function.
clear all; close all; clc; format compact;

% Digital Modulation
Nsym = 32*1024;
k = log2(16);
modScheme = 'qam';
s = DigitalModulation(Nsym,k,modScheme);

% OFDM
Nsc = 16;

% DFT Matrix
upsamp = 8;
Nc = Nsc;
nr = 1:1/upsamp:Nc;
Nr = length(nr);
nr = repmat(nr.',[1,Nc]);
nc = 1:Nc;
nc = repmat(nc,[Nr,1]);
df = 1.18;
%Dn = exp(-1j*2*pi*(nr-1).*(nc-1)/(Nsc-1)).'/sqrt(Nsc);
Dn = exp(-1j*2*pi*(nr-1).*(nc-1)/Nsc).'/sqrt(Nsc);

% Cyclic Prefix
Tcp = 0.0;
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
zo = zeros(size(s));

figure(1)
for i=1:Nsc

    %xt = Tcp*Dn'*so;
    z = zo;
    z(i,:) = s(i,:);
    xt = Tcp*Dn'*z;
    y = reshape(xt,1,[]);

    L = upsamp*512;
    [Sx,ff] = cpsd(y,y,L,round(L/4),[],1,'centered');
    plot(ff,Sx);
    hold on

end

xt = Tcp*Dn'*s;
y = reshape(xt,1,[]);
n = 0:length(y)-1;
y = y.*exp(-1j*2*pi*n/upsamp/2);

figure(2)
L = 8*1024;
[Sx,ff] = cpsd(y,y,L,round(L/4),[],1,'centered');
semilogy(ff,Sx);
hold on

%%
% Creating another OFDM function that works better, so that the frequencies are orthogonal.
clear all; close all; clc; format compact;

% Digital Modulation
Nsym = 32*1024;
k = log2(16);
modScheme = 'qam';
s = DigitalModulation(Nsym,k,modScheme);

% OFDM
Nsc = 16;
W = 0.5;
Tsym = 64;

% Cyclic Prefix and Zero Padding
guardLength = 0.0;
Tg = floor(guardLength*Tsym);
Tsc = Tsym - Tg;
Nsymsc = Nsym/Nsc;

% Symbols
s = reshape(s,Nsc,[]);

% Method 1
y = repmat(reshape(s.',1,[]),[Tsc,1]);
y = reshape(y,[],Nsc).';

% Method 2
% y = zeros(Nsc,Nsymsc*Tsc);
% idx = 1:Tsc:Nsymsc*Tsc;
% y(:,idx) = s;

% 
n = 0:Tsc-1;
fsc = 0:1/Tsym:(Nsc-1)/Tsym;
Dn = exp(1j*2*pi*fsc'*n);
Dn = repmat(Dn,[1,Nsymsc]);
y = y.*Dn;

% Prefix Insertion
guardType = 'cp';
xt = InsertCyclicPrefix(y,Tg,Tsym,guardType);

figure(1)
for i=1:Nsc
    %xt = Tcp*Dn'*so;
    z = y(i,:);
    L = 1024;
    [Sz,ff] = cpsd(z,z,L,round(L/4),[],1,'centered');
    plot(ff,Sz);
    hold on
end

fc = 0;
y = sum(xt,1);
n = 0:length(y)-1;
y = y.*exp(1j*2*pi*(fc-Nsc/2/Tsc)*n);

figure(2)
L = 512;
[Sx,ff] = cpsd(y,y,L,round(L/4),[],1,'centered');
semilogy(ff,Sx);
hold on


%%
clear all; close all; clc; format compact;

Ns = 8*1024;
fc = 0;
y = OFDMSignalGenerator(Ns,fc);

figure(1)
L = 1024;
[Sx,ff] = cpsd(y,y,L,round(L/4),[],1,'centered');
semilogy(ff,Sx);
hold on