function s = GetSymbol(idx,k,modscheme)

% modscheme = {'psk','qam'}

if rem(k,2) ~= 0 && strcmp(modscheme,'qam')
    error("QAM modulation my be 4, 16, 64, etc. M=2^k.");
end
M = 2^k;
if idx > M || idx <= 0
    error("idx must be one of the Mth symbols.");
end

switch modscheme
    case 'psk'
        n = (0:M-1);
        y = exp(1j*2*pi*n/M);
    case 'qam'
        sqM = sqrt(M);
        y = (0:sqM-1);
        y = 2*(y-sqM/2+0.5);
        y = repmat(y,sqM,1);
        y = y - 1j*y.';
        Esym = sum(sum(y.*conj(y)))/M;
        y = y/sqrt(Esym);
    otherwise
        error("Invalid modulation scheme.");
end
s = y(idx);

end