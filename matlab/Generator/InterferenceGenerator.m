function [y,f] = InterferenceGenerator(N,fbounds)

f = (0:(N-1))/N-0.5;
H = zeros(1,N);

fl = fbounds(1);
fh = fbounds(2);
fc = abs(fh-fl)/2;
fcenter = fl+fc;
b = 0.25;
T = 1/2/fc;
for i=1:N
    fo = abs(f(i));
    bd1 = (1-b)/2/T;
    bd2 = (1+b)/2/T;
    if fo <= bd1
        H(i) = 1;
    else
        if fo <= bd2
            H(i) = 0.5*(1+cos(pi*T/b*(fo-bd1)));
        end
    end
end

n = 0:N-1;
x = (randn(1,N)+1j*randn(1,N))/sqrt(2);
%x = randn(1,N);

%h = fftshift(ifft(H));
%z = sigConvolution(x,h);
z = ifftshift(ifft(fftshift(fft(x)).*H));
Z = fftshift(fft(z));

y = z.*exp(1j*2*pi*(fcenter+0.5)*n);

y = y/sqrt(y*y')*sqrt(N);



