clear all 
close all
clc

load("konstante.mat");


df = -150e9/2/pi; %detuning
m = 3; % mod u koji se injektuje, moze da bude poz i neg broj

alpha = 4;
I = 1.2*9.64234e-3;

% tau = 100e-9;

[gainNkonc,eng0,nth] = gOdE(gainN,eng, nsve,gth);

f0 = q*eng0/h; 

finj = f0 + fsr * m + df;

lambdainj = c0/finj*1e9;

%m je mod u koji injektujemo a fsr je intermodalni prostor
brMod = 70; %sa jedne strane ako je centralni mod indeksa 0
ind = (-brMod:brMod);
modovi = f0 + ind * fsr; %modovi od -brMod do brMod
modoviEng = modovi * h / q;

% interpolacija gaina:
fitovaniGain = zeros(length(ind), 5); %matrica koeficijenata za svaki mod
a = 799;
% a = length(nsve);
for j = (1:length(ind))
    g = zeros(size(gainN(1:799,1))); %ovde upisujem g(n) za odredjeno lambda
    for k = 1:length(nsve(1:a))
        g(k) = interp1(eng, gainN(k,:), modoviEng(j));
    end 
    f = fit(nsve(1:a)', g,'poly4');
    fitovaniGain(j,1) = f.p1;
    fitovaniGain(j,2) = f.p2;
    fitovaniGain(j,3) = f.p3;
    fitovaniGain(j,4) = f.p4;
    fitovaniGain(j,5) = f.p5;
end
grid on 
grid minor

kc = c0*1e2/2/neff/L * log(1+(1-R)/sqrt(R));

tspan = linspace(0,150e-9,150000);

poc = ones(1,brMod*2+1);
y0 = [1e16,1e16,1e8*poc,0];

[t,y] = ode45(@(t,y) odeInjLocking(t,y, Rspsve, nsve, I,Ab,Aw,Bb,Cb,Cw, Gama, Vtot, Vw, betasp, eps, etainj, q, taubw, tauwb, vg, alpha, m, kc, df, alphaI,  ind, brMod, fitovaniGain, snagaImpulsiDole(t), 1, eta0, h, lambdainj, taup, c0), tspan, y0);

Pin = zeros(size(t));

for i = 1:length(Pin)
    Pin(i) = snagaImpulsiDole(t(i));
end

Pout = eta0 * h * f0 * y(:,brMod+1+m+2) * Vw / taup / Gama;

figure
yyaxis left
plot(t*1e9, Pin*1e3, 'LineWidth',1.1)
ylabel('Pin[mW]')
yyaxis right
plot(t*1e9,Pout*1e3, 'LineWidth',1.1)
grid on 
grid minor
ylabel('Pout[mW]')
xlabel('t [ns]')
title('Mod u koji injektujemo')
legend('Pin','Pout')
