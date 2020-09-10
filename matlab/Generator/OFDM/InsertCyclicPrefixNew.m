function y = InsertCyclicPrefixNew(x,Tg,Tsym,type)

Tsc = Tsym - Tg;
[Nsc,~] = size(x);

y = reshape(x.',Tsc,[]).';
switch type
    case 'cp'
        y = [y(:,(end-Tg+1):end),y];
    case 'zp'
        y = [zeros(length(y),Tg),y];
    otherwise
        error('Invalid guard type.');
end
y = reshape(y.',1,[]).';
y = reshape(y,[],Nsc).';