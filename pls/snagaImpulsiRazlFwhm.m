function P = snagaImpulsiRazlFwhm(t,FWHM)
dP = 0.07e-3; %visina prvog impulsa
% P0 = 0.02e-3; %base snaga
P0 = 0;
o = FWHM./sqrt(8*log(2))*1e-9;
s = size(o);
endmi = 1e-9*s(1);
mi = 10e-9:1e-9:endmi;

if t < 9.5e-9 
    P = 0;
else
    P = P0;
    for i = 1:length(mi)
        % P = P + exp(-1/2*(t-mi(i)).^2/o(i)^2)*Pin(i);
        P = P + exp(-1/2*(t-mi(i)).^2/o(i)^2)*dP*i;
    end

end 
end