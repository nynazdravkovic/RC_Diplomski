close all 
clear all 
clc

load("konstante.mat");

[gainNkonc,eng0,nth] = gOdE(gainN,eng, nsve,gth);

f0 = q*eng0/h; 


%m je mod u koji injektujemo a fsr je intermodalni prostor
brMod = 10; %sa jedne strane ako je centralni mod indeksa 0
ind = (-brMod:brMod);
modovi = f0 + ind * fsr; %modovi od -brMod do brMod
modoviEng = modovi * h / q; 

% interpolacija gaina:
% fitovaniGain = zeros(length(ind), 10); %matrica koeficijenata za svaki mod
fitovaniGain = zeros(length(ind), 5); %matrica koeficijenata za svaki mod
% figure
a = 799;
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
I = [12e-3, 20e-3, 25e-3];
tau = 20e-9;
%naredne konstante su nebitne ali su potrebni argumenti funkcije koja
%resava dif jne koju koristim
alpha = 4;
m = 1;
df = 10e-9;
f0 = q*eng0/h; 
finj = f0 + fsr * m + df;

lambdainj = c0/finj*1e9;tspan1 = linspace(0,tau,1000);

poc = ones(1,brMod*2+1);
y01 = [1e16,1e16,1e8*poc,0];
% spektri = zeros(length(I), length(ind));
figure
for i = 1:length(I)
    [t,y] = ode45(@(t,y) odeInjLocking(t,y, Rspsve, nsve, I(i),Ab,Aw,Bb,Cb,Cw, Gama, Vtot, Vw, betasp, eps, etainj, q, taubw, tauwb, vg, alpha, m, kc, df, alphaI,  ind, brMod, fitovaniGain,snaga(t, 101, 20e-9), 0, eta0, h, lambdainj, taup, c0), tspan1, y01);
    
    % figure
    % plot(t,y(:,brMod+1))
    % hold on
    % plot(t,y(:,brMod+1+1))
    % hold on
    % plot(t,y(:,brMod+1+2))
    % hold on
    % plot(t,y(:,brMod+1+3))
    % hold on
    % plot(t,y(:,brMod+1+4))
    hold on
    P = eta0 * h * f0 * y(end,3:length(ind)+2) * Vw / taup / Gama;
    scatter(modovi,10*log10(P*1e3),"filled")
end

xlabel('frekvencija');
ylabel('Snaga [dBm]')
title('Spektar')
legend('I = 12 mA','I = 20 mA','I = 25 mA')
grid on
grid minor
