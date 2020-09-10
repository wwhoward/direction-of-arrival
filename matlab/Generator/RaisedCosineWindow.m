function H = RaisedCosineWindow(N,b)

%N = 2*T;
f = linspace(-0.5,0.5,N);
H = zeros(1,N);

fl = -0.25;
fh = 0.25;
fc = abs(fh-fl)/2;
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

end

