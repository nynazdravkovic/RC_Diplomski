function [P,fi] = snagaRCDrugi(input, t, B, fwhm,finj,Pin,phaseIn,tspan,K)
%t je jedan vremenski korak

o = fwhm/sqrt(8*log(2));
mi = 1/B*(1:length(input))+1/B*0.5;
% P = 0;
w = 2*pi*finj;
% w=0;
% Pint = interp1(tspan,Pin,t);
% phaset = interp1(tspan,phaseIn,t);
% % phaset = 0;
% Ein = sqrt(Pint)*exp(1j*(phaset+t*w))*K;
% % Ein = (Pint);
Ein = 0;
p = 0;
for i = 1:length(mi)
    p = p + exp(-(t-mi(i)).^2/o^2)*input(i)*1e-3;
end
E = Ein + sqrt(p)*exp(1j*w*t);
% E = Ein + p;
P = abs(E)^2;
% P = E;
% fi = phaset;
fi = angle(E);
end