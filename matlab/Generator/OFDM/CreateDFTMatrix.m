function D = CreateDFTMatrix(v,N)

Nr = v(1);
Nc = v(2);
nr = 1:Nr;
nr = repmat(nr.',[1,Nc]);
nc = 1:Nc;
nc = repmat(nc,[Nr,1]);
D = exp(-1j*2*pi*(nr-1).*(nc-1)/N);

end

