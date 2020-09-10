% Subject:  OFDM Function 2
% Date:     June 5, 2020
% Name:     Daniel Tait

% Creates an OFDM signal.

%%
% Creating another OFDM function that works better, so there can be a window.
clear all; close all; clc; format compact;

% Digital Modulation
Nsym = 32 + 4;
k = log2(16);
modScheme = 'qam';
s = DigitalModulation(Nsym,k,modScheme);

% OFDM
Nsc = 4;
W = 0.5;
Tsym = 128;

% Cyclic Prefix and Zero Padding
guardLength = 0.25;
Tg = floor(guardLength*Tsym);
Tsc = Tsym - Tg;
Nsymsc = Nsym/Nsc;
nsc = 0:Tsc-1;
ng = (Tsc-Tg):Tsc-1;
n = [ng,nsc];

% Symbols
s = reshape(s,Nsc,[])
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

% Guard Insertion
guardType = 'cp';
fsc = 0:1/Tsc:(Nsc-1)/Tsc;
switch guardType
    case 'cp'
        % Reshaping for carrier modulation
        ye = repmat(reshape(se.',1,[]),[Tsym,1]);
        ye = reshape(ye,[],Nsc).';
        yo = repmat(reshape(so.',1,[]),[Tsym,1]);
        yo = reshape(yo,[],Nsc).';
        
        % Window
        b = 0.4;
        w = RaisedCosineWindow(2*Tsym,b);
        w = w/sqrt(w*w');
        wo = circshift(w,-Tsym/2);
        wo = repmat(wo,[1,round(Nsymsc/2)]);
        we = circshift(w,Tsym/2);
        we = repmat(we,[1,round(Nsymsc/2)]);
        wo(end-Tsym/2:end) = 0;
        we(1:Tsym/2) = 0;
        wo = wo(1:length(xe));
        we = we(1:length(xe));
        
        % Centering symbols over the windows
        yec = circshift(ye,Tsym/2,2);
        yo(1:Tsym) = yo(end);
        yoc = circshift(yo,-Tsym/2,2);
        
        % Carrier Modulation
        D = exp(1j*2*pi*fsc'*n);
        D = repmat(D,[1,Nsymsc]);
        xe = ye.*D;
        xo = yo.*D;

%         % Reshaping for carrier modulation
%         y = repmat(reshape(s.',1,[]),[Tsym,1]);       % Previous way
%         y = reshape(y,[],Nsc).';
%         
%         % Carrier Modulation
%         D = exp(1j*2*pi*fsc'*n);
%         D = repmat(D,[1,Nsymsc]);
%         y = y.*D;
%         xt = y;
    case 'zp'
        % Reshaping for carrier modulation
        y = repmat(reshape(s.',1,[]),[Tsc,1]);
        y = reshape(y,[],Nsc).';
        
        % Carrier Modulation
        D = exp(1j*2*pi*fsc'*nsc);
        D = repmat(D,[1,Nsymsc]);
        y = y.*D;
        
        % Zero Padding Insertion
        xt = InsertCyclicPrefix(y,Tg,Tsym,guardType);
    otherwise
        error("Invalid guard type.");
end

figure(1)
plot(we,'-')
hold on
plot(abs(yec(1,:)),'-')

figure(3)
plot(wo,'-')
hold on
plot(abs(yoc(1,:)),'-')

% figure(1)
% for i=1:Nsc
%     %xt = Tcp*Dn'*so;
%     z = y(i,:);
%     L = 1024;
%     [Sz,ff] = cpsd(z,z,L,round(L/4),[],1,'centered');
%     plot(ff,Sz);
%     hold on
% end

fc = 0;
y = sum(xt,1);
n = 0:length(y)-1;
y = y.*exp(1j*2*pi*(fc-Nsc/2/Tsc)*n);

figure(2)
L = 1024;
[Sx,ff] = cpsd(y,y,L,round(L/4),[],1,'centered');
semilogy(ff,Sx);
hold on

hold off

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