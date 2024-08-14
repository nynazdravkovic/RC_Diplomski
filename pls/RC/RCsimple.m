% treba da:
% otvori fajl koji mi je mladen poslao 
% Uzme prvu vrednost iz dataseta i masku 
clc
load('C:\Users\ninaz\Desktop\diplomski\brzinskeJne\RC_numeric_data\masks_30nodes_100averg.mat')
load('C:\Users\ninaz\Desktop\diplomski\brzinskeJne\RC_numeric_data\SantaFe_laser_dataset_norm.mat')
load('C:\Users\ninaz\Desktop\diplomski\brzinskeJne\RC_numeric_data\weights_A_0_8_K_0_2_nodes_30.mat')
%konst za laser
load('konstanteZaRC.mat')
%konst za rc
A = 0.8;
K = 0.2;
NMSE = 0.0979;
Nv = 30;
J_input = A*M1_full(1,1,:)*SantaFe_laser_dataset_norm(1);
% napravi impulse koji imaju sirinu od 100ps a visinu koja odg J_input
B = 1/1e-9;
fwhm = 100e-12; 
tend = 1/B*length(J_input);
tspan = linspace(0,tend, length(J_input)*1000);
%Pin
Pin = zeros(1,length(tspan));
for i = 1:length(tspan)
    Pin(i) = snagaRC(J_input,tspan(i),B,fwhm);
end
figure
plot(tspan,Pin)
title('Ulazni impulsi')
grid on 
grid minor
%%
% putim impulse kroz nelin elem 
poc = ones(1,brMod*2+1);
y0 = [1e16,1e16,1e8*poc,0];

[t,y] = ode45(@(t,y) odeInjLocking(t,y, Rspsve, nsve, I,Ab,Aw,Bb,Cb,Cw, Gama, Vtot, Vw, betasp, eps, etainj, q, taubw, tauwb, vg, alpha, m, kc, df, alphaI,  ind, brMod, fitovaniGain, snagaRC(J_input,t,B,fwhm), 1, eta0, h, lambdainj, taup, c0), tspan, y0);

P1 = eta0 * h * f0 * y(:,brMod+1+m+2) * Vw / taup / Gama;

faza1 = y(:,end);

%%
% [yuPrvi,xuPrvi,fwhmuPrvi,proms] = findpeaks(abs(Pin(8000:147500))*1e3, t(8000:147500)*1e9, 'MinPeakDistance', 0.8);
% [yiPrvi,xiPrvi,fwhmiPrvi,proms] = findpeaks(abs(Pout(8000:147000))*1e3, t(8000:147000)*1e9, 'MinPeakDistance', 0.8);
% save('dt180mneg9fwhm100dP005Pbase0.mat','Pin','t','Pout','yuPrvi','xuPrvi','yiPrvi','xiPrvi', "fwhmiPrvi")
%%
figure
yyaxis left
plot(t*1e9, Pin*1e3, 'LineWidth',1.1)
ylabel('Pin[mW]')
hold on
% scatter(xuPrvi,yuPrvi,'x')
% hold on
yyaxis right
plot(t*1e9,Pout*1e3, 'LineWidth',1.1)
% hold on
% scatter(xiPrvi,yiPrvi,'x')
grid on 
grid minor
ylabel('Pout[mW]')
xlabel('t [ns]')
title('Mod u koji injektujemo')
legend('Pin','pik Pin','Pout','pik Pout')
%% sabiram polja
J_input1 = A*M1_full(2,1,:)*SantaFe_laser_dataset_norm(2);

Pin = zeros(1,length(tspan));
fiin = zeros(1,length(tspan));
for i = 1:length(tspan)
    [Pin(i),fiin(i)] = snagaRCDrugi(J_input,tspan(i),B,fwhm,finj,Pout,faza,tspan,K);
end
% figure
% plot(tspan,faza)
% hold on
% plot(tspan,fiin,"--")
% title('Faza signala')

%%
figure
plot(tspan,Pin,'--')
title('Impulsi nakon sabiranja')
hold on
plot(tspan,p1)
% hold on
% plot(tspan,sabSnage)
legend( 'faza(t) = 0','faza(t) != 0')
% title('Ulazni impulsi')
xlabel('vreme (s)')
ylabel('snaga (W)')
grid on 
grid minor

% poc uslovi su mi kraj prethodnoh stanja 
% y0 = y(end,:);
% 
% 
% [t,y] = ode45(@(t,y) odeInjLocking(t,y, Rspsve, nsve, I,Ab,Aw,Bb,Cb,Cw, Gama, Vtot, Vw, betasp, eps, etainj, q, taubw, tauwb, vg, alpha, m, kc, df, alphaI,  ind, brMod, fitovaniGain, snagaRC(J_input,t,B,fwhm), 1, eta0, h, lambdainj, taup, c0), tspan, y0);
% 
% Pout = eta0 * h * f0 * y(:,brMod+1+m+2) * Vw / taup / Gama;
% faza = y(:,end);
% 

% pusim kroz vlakno
% napravim maska * koef * data (2)
% saberem polja sqrt(Pin) * exp(1i*w*t) + sqrt(Pout) * exp(1i*w*t+faza)
%%
Ein = zeros(1,length(tspan));
for i = 1:length(tspan)
    Ein(i) = sqrt(Pout(i))*exp(1j*(faza(i)));
end

figure
plot(t,angle(Ein))
title('Zavisnost faze od vremena svedena na opseg od -Pi, do Pi')
grid on
grid minor
ylabel('faza')
xlabel('vreme (s)')