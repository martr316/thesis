function f = freq_func(r,freq_Hz,modes,w_m,Y,n_m,X_m,X_m_d,Y_fft,NFFT)

epsilon_u_piezo = zeros(length(modes),length(freq_Hz));

%iterator = 1;

r_star = zeros(1,3);

for n= 1:modes
    r_star(n) = r(n)./4;
end

Q_m = 1./(w_m(1:3).*n_m(1:3).*r_star);

% for w = freq_Hz;
%     epsilon_u_piezo(iterator) = (1/Y)*sum(((n_m(1:modes).*X_m(1:modes).*X_m_d(1:modes)).*(1-((w*2*pi)./w_m(1:modes)).^2+1i.*((w*2*pi)./w_m(1:modes)).*(1./Q_m(1:modes))))./...
%         (((1-((w*2*pi)./w_m(1:modes)).^2).^2)+(((w*2*pi)./w_m(1:modes)).*(1./Q_m(1:modes))).^2));
%     iterator = iterator + 1;
% end

for n = 1:modes;
    epsilon_u_piezo(n,:) = (1/Y)*(((n_m(n)*X_m(n)*X_m_d(n))*(1-((freq_Hz*2*pi)/w_m(n)).^2+1i*((freq_Hz*2*pi)/w_m(n))*(1/Q_m(n))))./...
        (((1-((freq_Hz*2*pi)/w_m(n)).^2).^2)+(((freq_Hz*2*pi)/w_m(n))*(1/Q_m(n))).^2));
end

f = sum(((10.^((11 + 20*log10(abs(sum(epsilon_u_piezo,1))))./20))' - (2*abs(Y_fft(1903:NFFT/2+1)))).^2);

end