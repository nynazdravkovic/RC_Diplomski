% clear all 
%close all
clc
clear t 
load("konstante.mat");
% eps = 0;

df = -180e9/2/pi; %detuning
m = -9; % mod u koji se injektuje, moze da bude poz i neg broj
alpha = 4;
I = 1.2*9.64234e-3;
brMod = 40; %sa jedne strane ako je centralni mod indeksa 0

% tau = 100e-9;

[gainNkonc,eng0,nth] = gOdE(gainN,eng, nsve,gth);

f0 = q*eng0/h; 

finj = f0 + fsr * m + df;

lambdainj = c0/finj*1e9;
%%
%m je mod u koji injektujemo a fsr je intermodalni prostor
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


kc = c0*1e2/2/neff/L * log(1+(1-R)/sqrt(R));
save('C:\Users\ninaz\Desktop\diplomski\brzinskeJne\pls\RC\konstanteZaRC.mat')
tspan = linspace(0,150e-9,150000);

poc = ones(1,brMod*2+1);

y0 = [1e16,1e16,1e8*poc,0];
%%
[t,y] = ode45(@(t,y) odeInjLocking(t,y, Rspsve, nsve, I,Ab,Aw,Bb,Cb,Cw, Gama, Vtot, Vw, betasp, eps, etainj, q, taubw, tauwb, vg, alpha, m, kc, df, alphaI,  ind, brMod, fitovaniGain, snagaImpulsi(t), 1, eta0, h, lambdainj, taup, c0), tspan, y0);


Pin = zeros(size(t));

for i = 1:length(Pin)
    Pin(i) = snagaImpulsi(t(i));
end

Pout = eta0 * h * f0 * y(:,brMod+1+m+2) * Vw / taup / Gama;
faza = y(:,end);
save('C:\Users\ninaz\Desktop\diplomski\brzinskeJne\prviProlazPout.mat','Pout','faza','c0','t')
%%
[yuPrvi,xuPrvi,fwhmuPrvi,proms] = findpeaks(abs(Pin(8000:147500))*1e3, t(8000:147500)*1e9, 'MinPeakDistance', 0.8);
[yiPrvi,xiPrvi,fwhmiPrvi,proms] = findpeaks(abs(Pout(8000:147000))*1e3, t(8000:147000)*1e9, 'MinPeakDistance', 0.8);
% save('dt180mneg9fwhm100dP005Pbase0.mat','Pin','t','Pout','yuPrvi','xuPrvi','yiPrvi','xiPrvi', "fwhmiPrvi")
%%
figure
yyaxis left
plot(t*1e9, Pin*1e3, 'LineWidth',1.1)
ylabel('Pin[mW]')
hold on
scatter(xuPrvi,yuPrvi,'x')
hold on
yyaxis right
plot(t*1e9,Pout*1e3, 'LineWidth',1.1)
hold on
scatter(xiPrvi,yiPrvi,'x')
grid on 
grid minor
ylabel('Pout[mW]')
xlabel('t [ns]')
title('Mod u koji injektujemo')
legend('Pin','pik Pin','Pout','pik Pout')

%% trazenje pikova

figure
plot(yuPrvi(1:138), yiPrvi(1:138), 'LineWidth',1.1)
% hold on 
% plot(yu25, yi25, 'LineWidth',1.1)
grid on 
grid minor
ylabel('Pout [mW]')
xlabel('Pi [mW]')
title('Zavisnost visine izlaznih od visine ulaznih pikova')
% legend('drugi prolaz nakon 80ps','25ps')(8000:147500)(8000:147500)
% 
%% drugi prolaz
load("C:\Users\ninaz\Desktop\diplomski\brzinskeJne\pls\nakonProlaska.mat")
Poutprviprolaz = POut(1:131073);
Poutprviprolaz(1) = Poutprviprolaz(2);
% Poutprviprolaz = Pout;
tspan = linspace(0,150e-9,150000);
tspan = tspan(1:131073);
[t,y] = ode45(@(t,y) odeInjLockingDrugiProlaz(t,y, Rspsve, nsve, I,Ab,Aw,Bb,Cb,Cw, Gama, Vtot, Vw, betasp, eps, etainj, q, taubw, tauwb, vg, alpha, m, kc, df, alphaI,  ind, brMod, fitovaniGain, Poutprviprolaz, tspan, 1, eta0, h, lambdainj, taup, c0), tspan, y0);
PoutDrugi = eta0 * h * f0 * y(:,brMod+1+m+2) * Vw / taup / Gama;
%%
[yuDrugiAT,xuDrugiAT,fwhmuDrugiAT,proms] = findpeaks(abs(Poutprviprolaz(8000:147500))*1e3, t(8000:147500)*1e9, 'MinPeakDistance', 0.8);
[yiDrugiAT,xiDrugiAT,fwhmiDrugiAT,proms] = findpeaks(abs(PoutDrugi(8000:147500))*1e3, t(8000:147500)*1e9, 'MinPeakDistance', 0.8);
% save('DRUGIdt180mneg9fwhm100dP005Pbase0.mat','PoutDrugi','t','Poutprviprolaz','yuDrugiAT','xuDrugiAT','yiDrugiAT','xiDrugiAT', "fwhmiDrugiAT")
%%
figure
yyaxis left
plot(t*1e9, Poutprviprolaz*1e3, 'LineWidth',1.1)
ylabel('Pin[mW]')
hold on
scatter(xuDrugiAT,yuDrugiAT,'x')
hold on
yyaxis right
plot(t*1e9,PoutDrugi*1e3, 'LineWidth',1.1)
hold on
scatter(xiDrugiAT,yiDrugiAT,'x')
grid on 
grid minor
ylabel('Pout[mW]')
xlabel('t [ns]')
title('Mod u koji injektujemo')
legend('Pin','pik Pin','Pout','pik Pout')
%% trazenje pikova

figure
plot(yuDrugiAT(1:138), yiDrugiAT(1:138), 'LineWidth',1.1)
grid on 
grid minor
ylabel('Pout [mW]')
xlabel('Pi [mW]')
title('Zavisnost visine izlaznih od visine ulaznih pikova')
% legend('drugi prolaz nakon 80ps','25ps')
% 

%%
% %% poredjenje sa razl fwhm
% poc = ones(1,brMod*2+1);
% 
% y0 = [1e16,1e16,1e8*poc,0];
% 
% [t,y] = ode45(@(t,y) odeInjLocking(t,y, Rspsve, nsve, I,Ab,Aw,Bb,Cb,Cw, Gama, Vtot, Vw, betasp, eps, etainj, q, taubw, tauwb, vg, alpha, m, kc, df, alphaI,  ind, brMod, fitovaniGain, snagaImpulsi(t), 1, eta0, h, lambdainj, taup, c0), tspan, y0);
% 
% 
% Pin = zeros(size(t));
% 
% for i = 1:length(Pin)
%     Pin(i) = snagaImpulsiRazlFwhm(t(i),fwhmiPrvi, yiPrvi);
% end
% 
% PoutF30 = eta0 * h * f0 * y(:,brMod+1+m+2) * Vw / taup / Gama;
% %%
% 
% [yuPrviF30,xuPrviF30,fwhmuPrviF30,proms] = findpeaks(abs(Pin(8000:147500))*1e3, t(8000:147500)*1e9, 'MinPeakDistance', 0.8);
% [yiPrviF30,xiPrviF30,fwhmiPrviF30,proms] = findpeaks(abs(PoutF30(8000:147000))*1e3, t(8000:147000)*1e9, 'MinPeakDistance', 0.8);
% save('F30dt100m4fwhm100dP002Pbase0.mat','Pin','t','PoutF','yuPrviF','xuPrviF','yiPrviF','xiPrviF', "fwhmiPrviF")
% figure
% yyaxis left
% plot(t*1e9, Pin*1e3, 'LineWidth',1.1)
% ylabel('Pin[mW]')
% hold on
% scatter(xuPrviF30,yuPrviF30,'x')
% hold on
% yyaxis right
% plot(t*1e9,PoutF30*1e3, 'LineWidth',1.1)
% hold on
% scatter(xiPrviF30,yiPrviF30,'x')
% grid on 
% grid minor
% ylabel('Pout[mW]')
% xlabel('t [ns]')
% title('Mod u koji injektujemo')
% legend('Pin','pik Pin','Pout','pik Pout')
% 
% %% trazenje pikova
% 
% figure
% plot(yuPrviF30(1:128), yiPrviF30(1:128), 'LineWidth',1.1)
% % hold on 
% % plot(yu25, yi25, 'LineWidth',1.1)
% grid on 
% grid minor
% ylabel('Pout [mW]')
% xlabel('Pi [mW]')
% title('Zavisnost visine izlaznih od visine ulaznih pikova')
% % legend('drugi prolaz nakon 80ps','25ps')
%%
poc = ones(1,brMod*2+1);

y0 = [1e16,1e16,1e8*poc,0];

[t,y] = ode45(@(t,y) odeInjLocking(t,y, Rspsve, nsve, I,Ab,Aw,Bb,Cb,Cw, Gama, Vtot, Vw, betasp, eps, etainj, q, taubw, tauwb, vg, alpha, m, kc, df, alphaI,  ind, brMod, fitovaniGain, snagaImpulsiRazlFwhm(t,fwhmiPrvi), 1, eta0, h, lambdainj, taup, c0), tspan, y0);


Pin = zeros(size(t));

for i = 1:length(Pin)
    Pin(i) = snagaImpulsiRazlFwhm(t(i),fwhmiPrvi);
end

PoutF = eta0 * h * f0 * y(:,brMod+1+m+2) * Vw / taup / Gama;
%%
[yuPrviF,xuPrviF,fwhmuPrviF,proms] = findpeaks(abs(Pin(8000:147500))*1e3, t(8000:147500)*1e9, 'MinPeakDistance', 0.8);
[yiPrviF,xiPrviF,fwhmiPrviF,proms] = findpeaks(abs(PoutF(8000:147000))*1e3, t(8000:147000)*1e9, 'MinPeakDistance', 0.8);
figure
yyaxis left
plot(t*1e9, Pin*1e3, 'LineWidth',1.1)
ylabel('Pin[mW]')
hold on
scatter(xuPrviF,yuPrviF,'x')
hold on
yyaxis right
plot(t*1e9,PoutF*1e3, 'LineWidth',1.1)
hold on
scatter(xiPrviF,yiPrviF,'x')
grid on 
grid minor
ylabel('Pout[mW]')
xlabel('t [ns]')
title('Mod u koji injektujemo')
legend('Pin','pik Pin','Pout','pik Pout')


%% plotovanje svih s kriva zajedno 
figure
% plot(yuDrugi(1:120),yiDrugi(1:120))
% hold on
plot(yuPrvi(1:138),yiPrvi(1:138), 'LineWidth',1.1)
hold on
plot(yuDrugiAT(1:138),yiDrugiAT(1:138), 'LineWidth',1.1)
hold on 
plot(yuPrviF(1:129), yiPrviF(1:129), 'LineWidth',1.1)

% % hold on
% % plot(yu300,yi300)
% % hold on
% % plot(yu400,yi400)
% % hold on 
% % plot(yu500,yi500)
% legend('Drugi prolaz','Prvi prolaz', 'Drugi prolaz kada je prvi atenuiran 2 puta')
legend('Prvi prolaz', 'Drugi prolaz',  'Privi prolaz sa sirinama impulsa nakon prvog prolaza')

grid on
grid minor
ylabel('Pout [mW]')
xlabel('Pi [mW]')
title('Zavisnost visine izlaznih od visine ulaznih pikova')
