close all
clear all 
clc

%fundamentalne const
h = 6.62607015e-34;
c0 = 299792458;
q = 1.60217663e-19;
fsr = 0.200328e12; %intermodaln

%const
Gamaqw = 1.589e-2;
neff = 3.9; %3.19
ntr = 1.5e18;

%gubici
K0 = 20; %20.676 menjam
K1 = 0; 
nrw =3.59;
R = ((neff-1)/(neff+1))^2;

%dimenzije
%L = 235e-4;
L = c0/(2*fsr*neff)*1e2;

W = 1.5e-4;
Nqw = 4; % 3
Hqw = 8.5e-9*1e2;
Hbin = 20e-9*1e2;
Hbout = 30e-9*1e2;
Hsch = 60e-9*1e2;

%racunanje gth
Gama = Nqw*Gamaqw;
vg = c0/neff*1e2;

alphaI = K0 + Gama * K1 * ntr;
alphaM = 1/L * log(1/R);
alphaTot = (alphaM + alphaI)/Gama; 
gth = alphaTot;
disp(gth) 
taup = 1/Gama/vg/gth;
eta0 = 0.5*alphaM/(alphaM + alphaI);

%racunanje potrebnih vr
Hw = Nqw*Hqw;
Hb = 2*(Hsch+Hbout)+Hw+(Nqw-1)*Hbin;
Vtot = W*Hb*L;
Vw = W*Hw*L;
Vqw = W*Hqw*L;

%konstante za brzinske
etainj = 0.8; %0.8

Ab = 6.7e7;
Bb = 0.7e-10;
Cb = 3e-28;

Aw = 1.1e8;
Cw = 0.207e-27;% 1.1e-27menjala

eps = 2.85e-17;
betaTilda = 3.74e-24;

taubw = 5.6e-12;
tauwb = 37.4e-12; %

betasp = Gamaqw*vg*betaTilda/Vqw;


%otvaranje datoteka
load("C:\Users\ninaz\Desktop\diplomski\brzinskeJne\RsptotVSn3QW.mat");
load("C:\Users\ninaz\Desktop\diplomski\brzinskeJne\nkoncrange3QW.mat");
load("C:\Users\ninaz\Desktop\diplomski\brzinskeJne\gain3QW.mat");
load("C:\Users\ninaz\Desktop\diplomski\brzinskeJne\photonenergy3QW.mat");
eng = photonenergy;
Rspsve = Expression1;
nsve = nkoncrange;

[gainNkonc,eng0,nth] = gOdE(gainN,eng, nsve, gth);
centralnoLambda = c0/(q*eng0/h)*1e9;
save("konstante.mat")

%pozivanje ODE

tspan = linspace(0,10e-9,1000);
y0 = [0,ntr,0];
I = [0.18,2.38,3.93,5.48,6.23,6.98,7.63,8.18,8.63,9.08,9.48,9.93,10.48,11.03,12.03,12.58,13.28,13.83,14.23,14.83,15.38,16.08,16.68,17.23,18.13,19.03];
PIzmereno = [2.95E-03,2.69E-02,6.96E-02,1.57E-01,2.30E-01,3.44E-01,5.10E-01,7.50E-01,1.11E+00,1.82E+00,3.47E+00,1.27E+01,5.33E+01,9.92E+01,1.91E+02,2.35E+02,2.85E+02,3.31E+02,3.63E+02,4.13E+02,4.57E+02,5.31E+02,5.75E+02,6.20E+02,6.97E+02,7.63E+02];
S = zeros(1,length(I));

for i=1:length(I)
    [S(i),k, j] = resiBrzinskeG(Rspsve, nsve,gainNkonc,I(i), tspan, y0,Ab,Aw,Bb,Cb,Cw, Gama, Vtot,Vw, taup,betasp,eps,etainj,q,taubw,tauwb,vg);
end

P = eta0 * h * c0 * 1e2 / centralnoLambda * 1e7 * S * Vw / taup / Gama;

koef = (P(26)-P(20))*1e3/(I(26)-I(20));

%plotovanje
figure()
plot(I,P*1e3, LineWidth=1.2)
hold on
scatter(I,PIzmereno/1e3,20,"filled")
% plot(I,PIzmereno/1e3,'.',LineWidth=4);
grid on
grid minor
legend("model","izmereno")
title('IP kka');
xlabel('I [mA]');
ylabel('P[mW]');
% 
% %mod krive:
I = [20, 30, 40];
figure()
for i=1:length(I)
    [S,nw,nb] = resiBrzinskeG(Rspsve, nsve,gainNkonc,I(i), tspan, y0,Ab,Aw,Bb,Cb,Cw, Gama, Vtot,Vw, taup,betasp,eps,etainj,q,taubw,tauwb,vg);
    [f,mod] = modKrive(nb, nw, S, nsve, gainNkonc, Rspsve, Vtot, Ab, Aw, Bb, Cb, Cw, taubw, tauwb, Vw, Vtot, vg, eps, Gama, betasp, etainj, q, taup);
    plot(f, mod, LineWidth=1.2)
    hold on
end

a = [-35.55555516, -29.7777788, -32.77777745, -35.55555516, -38.55555593];
a = a - a(1);
f1 = [0,4990859198,6471663931,7294332996,8555758475];

b = [-35.55555516, -29.7777788, -32.44444497, -35.55555516, -38.55555593];
b = b - b(1);
f0 = [0, 7157221485, 10009141325, 10639853019, 12806215306];

scatter(f1,a,"filled");
hold on
scatter(f0,b,"filled");
grid on
grid minor
legend('I = 20mA', 'I = 30mA','I = 40mA')
title('Modulaciona kriva');
ylabel('|dS/dI| [dB]');
xlabel('f [Hz]');
