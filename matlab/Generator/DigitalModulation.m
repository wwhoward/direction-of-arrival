function y = DigitalModulation(Nsym,k,modScheme)

% modscheme = {'psk','qam'}

if rem(k,2) ~= 0 && strcmp(modScheme,'qam')
    error("QAM modulation my be 4, 16, 64, etc. M=2^k.");
end

M = 2^k;
sym = randi([0,M-1],1,Nsym);

% modulation
switch modScheme
    case 'psk'
        y = exp(1j*2*pi*sym/M);
        % no need to normalize for average symbol energy
    case 'qam'
        y = zeros(Nsym,1);
        sqM = sqrt(M);
        for n=1:Nsym
            f = floor((sym(n)/sqM));
            y(n) = f;
            if rem(f,2) == 0
                y(n) = y(n) + rem(sym(n),sqM)*1j;
            else
                y(n) = y(n) + (sqM-1)*1j - rem(sym(n),sqM)*1j;
            end
        end
        y_off = (sqrt(M)/2-0.5)*(1+1j);
        y = 2*(y - y_off);
        % normalizing for average symbol energy
        yqam = (0:sqM-1);
        yqam = 2*(yqam-sqM/2+0.5);
        yqam = repmat(yqam,sqM,1);
        yqam = yqam - 1j*yqam.';
        Esym = sum(sum(yqam.*conj(yqam)))/M;
        y = y/sqrt(Esym);
        y = y.';
    otherwise
        error("Invalid modulation scheme.");
end