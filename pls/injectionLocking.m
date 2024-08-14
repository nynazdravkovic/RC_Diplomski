close all
clear all 
clc
%ovaj kod resava diferencijalne jednacine vise modova pri injektovanjoj
%snazi koja je prvo 0, onda stepenasto raste pa se stepenasto spusta.
%Rezultat je histerezis. 

load("konstante.mat"); %konstante se dobijaju ranovanjem koda pavlove konst, tu se deklarisu sve konst, crta se mod i pi kriva


df = -150e9/2/pi; %detuning
m = 3; % mod u koji se injektuje, moze da bude poz i neg broj

alpha = 4;
Ith = 9.64234e-3;
I = 1.2*Ith;

tau = 100e-9;
%racunanje energije centralnog moda preko gth, vracanje g(n, eng = eng0) i
%racunanje nth:
[gainNkonc,eng0,nth] = gOdE(gainN,eng, nsve,gth); 
%% interpolacija gaina za potrebne modove

f0 = q*eng0/h; %frekv centralnog moda
finj = f0 + fsr * m + df; %frekv injekcije 

lambdainj = c0/finj*1e9; 

brMod = 20; %broj modova sa jedne strane ako je centralni mod indeksa 0
ind = (-brMod:brMod); %indeksi svih modova
modovi = f0 + ind * fsr; %frekvencije svih modova od -brMod do brMod
modoviEng = modovi * h / q; %energije svih modova

% interpolacija gaina:

fitovaniGain = zeros(length(ind), 5); %matrica koeficijenata za svaki mod
figure
a = 799; %fituje se do konc n = 4e18, a 799 je njen indeks u nsve
% a = length(nsve);
for j = (1:length(ind))
    g = zeros(size(gainN(1:a,1))); %ovde upisujem g(n) za odredjeno lambda
    for k = 1:length(nsve(1:a))
        g(k) = interp1(eng, gainN(k,:), modoviEng(j));
    end 

    % f = fit(nsve(1:a)', g,'poly9');
    f = fit(nsve(1:a)', g,'poly4');
    fitovaniGain(j,1) = f.p1;
    fitovaniGain(j,2) = f.p2;
    fitovaniGain(j,3) = f.p3;
    fitovaniGain(j,4) = f.p4;
    fitovaniGain(j,5) = f.p5;
    % fitovaniGain(j,6) = f.p6;
    % fitovaniGain(j,7) = f.p7;
    % fitovaniGain(j,8) = f.p8;
    % fitovaniGain(j,9) = f.p9;
    % fitovaniGain(j,10) = f.p10;
    if j == brMod+1
        hold on
        plot(nsve(1:a), g, '--')
        hold on 
        plot(nsve(1:a), fitovaniGain(j,1)*nsve(1:a).^4 + fitovaniGain(j,2)*nsve(1:a).^3 +fitovaniGain(j,3)*nsve(1:a).^2 +fitovaniGain(j,4)*nsve(1:a) +fitovaniGain(j,5));
    end
end


kc = c0*1e2/2/neff/L * log(1+(1-R)/sqrt(R));

%% stacionarni rezim
%palim prvi laser
tspan1 = linspace(0,tau,1000);

poc = ones(1,brMod*2+1);
y01 = [1e16,1e16,1e8*poc,0];

%sustinski nebitni brojevi ovde, samo su potrebni kao arg odeInjLocking
brStep = 101;
dt = 20e-9;


[t,y] = ode45(@(t,y) odeInjLocking(t,y, Rspsve, nsve, I,Ab,Aw,Bb,Cb,Cw, Gama, Vtot, Vw, betasp, eps, etainj, q, taubw, tauwb, vg, alpha, m, kc, df, alphaI,  ind, brMod, fitovaniGain,snaga(t, brStep, dt), 0, eta0, h, lambdainj, taup, c0), tspan1, y01);
%uzimam resenja stacionarnog rezima kao pocetne uslove za injekciju
%% histerezisi

tspan = linspace(0,dt*brStep,1000*brStep);

y0 = y(end,:);

[t,y] = ode45(@(t,y) odeInjLocking(t,y, Rspsve, nsve, I,Ab,Aw,Bb,Cb,Cw, Gama, Vtot, Vw, betasp, eps, etainj, q, taubw, tauwb, vg, alpha, m, kc, df, alphaI,  ind, brMod, fitovaniGain,snaga(t, brStep, dt), 1, eta0, h, lambdainj, taup, c0), tspan, y0);

figure
plot(t,y(:,brMod+1+m+2))
grid on 
grid minor
ylabel('Sout [W]')
xlabel('t [s]')
title('Mod u koji injektujemo')


%crtanje histerezisa
Pin = zeros(size(tspan));
for i = 1:length(Pin)
    Pin(i) = snaga(t(i), brStep, dt);
end

Pin = Pin(1000:1000:end);
Sout = y(1000:1000:end, brMod+3+m);
Pout = eta0 * h * f0 * Sout * Vw / taup / Gama;
figure
plot(Pin*1000, Pout*1e3)
grid on 
grid minor
xlabel('Pin [mW]')
ylabel('Pout [mW]')
title('Histerezis')
