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
B = 1/5e-8;
fwhm = 100e-12; 

num_buff = SantaFe_laser_dataset_norm(3102:3201);
Output_test = SantaFe_laser_dataset_norm(3202:4201);

tend = 1/B*length(M1_full(1,1,:));
tspan = linspace(0,tend, length(J_input)*500); %100 tacaka po bitslotu
dataset = num_buff;
PSve = zeros(length(dataset), length(tspan));
fazaSve = zeros(length(dataset), length(tspan));

%% stacionarni rezim
tspan1 = linspace(0,10e-9,1000);

poc = ones(1,brMod*2+1);
y01 = [1e16,1e16,1e8*poc,0];

[t,y] = ode45(@(t,y) odeInjLocking(t,y, Rspsve, nsve, I,Ab,Aw,Bb,Cb,Cw, Gama, Vtot, Vw, betasp, eps, etainj, q, taubw, tauwb, vg, alpha, m, kc, df, alphaI,  ind, brMod, fitovaniGain, 0, 0, eta0, h, lambdainj, taup, c0), tspan1, y01);
%%
for i = 1:length(dataset)
    % napravi impulse koji imaju sirinu fwhm, bitske brzine B, 
    % a visinu koja odg J_input
    % i sabere ih sa prethodnim impulsima expet ako je i = 0
    disp(i)
    PIn =  ulazniImpulsi(i, K, tspan, M1_full, B, fwhm, dataset, A, 1, PSve, fazaSve, finj);
    tekst = "Ulazni impulsi u vremenu 100 ps, iteracija " + num2str(i);
    plotuj(tspan,PIn,tekst,'vreme [s]', 'snaga [W]');
    %pustim impulse kroz nelin elem 
    % if i == 1
    %     poc = ones(1,brMod*2+1);
    %     y0 = [1e16,1e16,1e8*poc,0];
    % else
    %     y0 = y(end,:);
    % end
    y0 = y(end,:);
    [t,y] = ode45(@(t,y) odeInjLocking(t,y, Rspsve, nsve, I,Ab,Aw,Bb,Cb,Cw, Gama, Vtot, Vw, betasp, eps, etainj, q, taubw, tauwb, vg, alpha, m, kc, df, alphaI,  ind, brMod, fitovaniGain, snagaRCDrugi(i,A, M1_full, 1,dataset, t, B, fwhm,finj,PSve,fazaSve,tspan,K), 1, eta0, h, lambdainj, taup, c0), tspan, y0);
    PSve(i,:) = eta0 * h * f0 * y(:,brMod+1+m+2) * Vw / taup / Gama;
    tekst = "Izlazni impulsi u vremenu 100 ps, iteracija " + num2str(i);
    plotuj(tspan,PSve(i,:),tekst,'vreme [s]', 'snaga [W]');
    fazaSve(i,:) = y(:,end);
end
save('B5e8Psve100ps30dbm.mat','PSve')