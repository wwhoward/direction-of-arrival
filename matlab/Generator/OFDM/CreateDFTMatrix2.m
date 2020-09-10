function D = CreateDFTMatrix2(Nsc,upsamp,N)

Nr = upsamp*Nsc;
Nc = Nsc;

nr = 0:1/upsamp:Nc-1;
nr = repmat(nr.',[1,Nc]);
nc = 0:Nc-1;
nc = repmat(nc,[Nc,1]);

D = exp(-1j*2*pi*(nr-1).*(nc-1)/N);

end

