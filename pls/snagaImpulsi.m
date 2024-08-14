function P = snagaImpulsi(t)
dP = 0.02e-3; %visina prvog impulsa
% P0 = 0.02e-3; %base snaga
P0 = 0;
FWHM = 0.1e-9;
o = FWHM/sqrt(8*log(2));
mi = 10e-9:1e-9:150e-9;

if t < 9e-9 
    P = 0;
else
    P = P0;
    for i = 1:length(mi)
        P = P + exp(-1/2*(t-mi(i)).^2/o^2)*dP*i;
    end

end 
end