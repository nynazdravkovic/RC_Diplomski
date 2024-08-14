function P = snagaImpulsiDole(t)
dP = 0.007e-3; %visine impulsa
P0 = 0.5e-3; %base snaga
FWHM = 0.5e-9;
o = FWHM/sqrt(8*log(2));
% mi1 = 20e-9;
% mi2 = 25e-9;
% mi3 = 30e-9;
% mi4 = 35e-9;
% mi5 = 40e-9;
% mi6 = 45e-9;
% mi7 = 50e-9;
mi = 20e-9:5e-9:145e-9;

% k=1;

% P = P0 + exp(-1/2*(t-mi(1)).^2/o^2)*dP + exp(-1/2*(t-mi(2)).^2/o^2)*2*dP+ exp(-1/2*(t-mi(3)).^2/o^2)*3*dP+ exp(-1/2*(t-mi(4)).^2/o^2)*4*dP+ exp(-1/2*(t-mi(5)).^2/o^2)*5*dP+ exp(-1/2*(t-mi(6)).^2/o^2)*6*dP+ exp(-1/2*(t-mi(7)).^2/o^2)*7*dP;
% P = P0 + exp(-1/2*(t-mi1).^2/o^2)*dP;
% disp('sta')
if t < 10e-9 
    P = P0;
else
    P = 0.2e-3;
    for i = 1:length(mi)
        P = P - exp(-1/2*(t-mi(i)).^2/o^2)*dP*i;
    end
end
end