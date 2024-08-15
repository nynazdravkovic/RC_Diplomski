function [P,fi] = snagaRCDrugi(i, A, M1_full, num_avg,dataset, t, B, fwhm,finj,PSve,fazaSve,tspan,K)
%t je jedan vremenski korak
input = A*M1_full(1,num_avg,:)*dataset(i); %ovo mi je visina impulsa za inpute dataseta(i)
w0 = 2*pi*finj;

o = fwhm/sqrt(8*log(2));
mi = 1/B*(1:length(input))+1/B*0.5;
p = 0.01*1e-3;

for k = 1:length(mi)
    p = p + exp(-(t-mi(k)).^2/o^2)*input(k)*1e-3;
end


if i ==1
    P = p; 
    fi = 0;
else
%fja za resavanje difjna uzima vrednosti t van tspana pa mora da se
%interpolira snaga i faza
    Pint = interp1(tspan,PSve(i-1,:),t, 'spline'); 
    % phaset = interp1(tspan,fazaSve(i-1,:),t, 'linear'); 
    % Ein = sqrt(Pint)*exp(1j*(phaset+t*w0))*sqrt(K);
    % E = sqrt(p).*exp(1j*w0*t) + Ein; 
    % P = abs(E).^2;
    % fi = angle(E);
    P = Pint*K+p;
end

end