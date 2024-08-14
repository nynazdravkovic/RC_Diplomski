close all 
clc

%impuls odredjene sirine
% lambda0 = 1550e-9; %centralni mod u um 
% lambda0 = lambdainj*1e-9;
finj = 1.916004195344080e+14;
c0 = 299792458;
f0 = finj;
w0 = 2*pi*f0;
lambda0 = c0/f0;
% 
% P0 = 1e-3; %snaga u mW
% o0 = 1e-8;
% w1 = 2*pi*c0/1.5e-6;
% w2 = 2*pi*c0/1.6e-6;
% w = linspace(w2,w1,1000);
% t1 = linspace(0,0.5e-9,10000);
% o = 100e-12;
% 
% P1 = 1e-3*exp(-(t1-0.25e-9).^2/o^2);
% P1 = fotoni; 
% faza1 = zeros(length(t1),1);
%parametri vlakna
D = -17e-12/1e-9/1e3; %ostavljeni su m
%alpha = 0.2/log10(exp(1))/10;
alpha=0;
L = 1000; %u m
nr = 1.5; %indeks prelamanja u vlaknu 
hbar = 1.054571817e-34;


% load('C:\Users\ninaz\Desktop\diplomski\brzinskeJne\pls\prviProlazPout.mat');
% load('prviProlazPout.mat');
%izvlaci mi Pout, fazu, c0 i t
% t1 = t; 
EIn = sqrt(P1).*exp(1i*faza1);
% EIn = sqrt(P1);
PIn = P1; %Pin u vlakno je Pout iz akt fje
t1 = t; 
%% interpolacija tacaka za odgovarajuci Fs
% delta_f = 9.5367e6;
delta_t = 1.1444e-13;
n = 18;

% n = 15;
% delta_t = (2.5e-8-t1(1))/2^n;

Ls = 2^n;
delta_f = 1/Ls/delta_t;

Fs = Ls*delta_f;

% delta_t = 1/Ls/delta_f;

tNovo = (1:Ls)*delta_t+t1(1);

fazaNovo = interp1(t1,faza1,tNovo,'linear');
PInNovo = interp1(t1,P1,tNovo,'linear');
figure
plot(tNovo,PInNovo);
hold on
plot(t1,PIn,'--')
legend('interp','pravo')
figure
plot(tNovo,fazaNovo);
hold on
plot(t1,faza1,'--')


EInInterp = sqrt(PInNovo).*exp(1j*fazaNovo);

%% fft i mnozenje sa transfer fjom

F = fft(EInInterp);
f = (-Ls/2:Ls/2-1)*Fs/Ls;
% 
% figure
% plot(f/1e9,abs(fftshift(F)), "LineWidth",1.2);
% xlabel('f [GHz]')
% ylabel('|F(f)|')
% title('Frekvencijski domen signala')
% xlim([-100,200])
% grid on

FSh = fftshift(F);

w = f*2*pi;
beta1 = nr/c0; % 1/vg
% beta1 = 0;
beta2 = -lambda0^2*D/2/pi/c0;
beta0 = 0;
beta = beta0 + 0.5*beta2.*((w).^2) +beta1.*(w);
H = exp(-(alpha/2+1i.*beta).*L);
H(1)=exp(-(alpha/2)*L);
%H = exp(1i*D*lambda0^2*L.*w.^2/4/pi/c0 - alpha/2*L);
Y1 = H.*FSh;
figure
plot(f/1e9,10*log10(abs(Y1)), "LineWidth",1.3);
hold on
plot(f/1e9,10*log10(abs(FSh)), "--","LineWidth",1.4);
xlabel('f [GHz]')
ylabel('|E*H| [dB]')
% xlim([-300,300])
legend("nakon vlakna", "pre vlakna")
title('Frekvencijski domen signala nakon prolaska kroz vlakno')
grid on

% figure
% plot(f/1e9,angle(Y1), "LineWidth",1.3);
% hold on
% plot(f/1e9,angle(FSh), "--","LineWidth",1.4);
% xlabel('f [GHz]')
% ylabel('phase(E*H)')
% xlim([-300,300])
% legend("nakon vlakna", "pre vlakna")
% title('Faza signala nakon prolaska kroz vlakno')
% grid on


%% inveryan fft
K = fftshift(Y1);
X = ifft(K);
% pomeraj = 9.64734e-9;
% pomeraj = 6.5189e-9;
pomeraj = 0.96498e-8;
% pomeraj = 0.8e-10;
pomeri = ceil(Fs*pomeraj);
X1 = circshift(X,pomeri);
% X1(tNovo<=1e-8) = EInInterp(tNovo<=1e-8);
figure
plot(tNovo,abs(X1).^2, "LineWidth",1.2)
% hold on
% plot(tNovo,normDisp, "LineWidth",1.2)
hold on
plot(tNovo,abs(EInInterp).^2,"--", "LineWidth",1.2)
xlabel("Vreme [s]")
ylabel("|Ein|^2")
title("Snaga signala pre i nakon vlakna")
grid on
grid minor
xlim([1.35e-8,1.64e-8])
legend('nakon vlakna', 'na ulazu u vlakno')
save('nakonProlaska.mat','X','tNovo')
%% 
figure
plot(tNovo,angle(X1), "LineWidth",1.2)
hold on
plot(t1,angle(EIn),"--", "LineWidth",1.2)
xlabel("Vreme [s]")
ylabel("angle(Ein)")
title("Faza signala pre i nakon vlakna")
grid on
grid minor
xlim([1.35e-8,1.64e-8])
legend('nakon vlakna', 'na ulazu u vlakno')
save('nakonProlaska.mat','X','tNovo')
% 
%%
[pks,locs,widths] =findpeaks(abs(X1).^2,tNovo);
widths
o2 = o*sqrt(1+(beta2*L/o^2)^2)/sqrt(2);
fwhm2 = o2*sqrt(8*log(2))
%% smanjujem broj tacaka 
%interpoliram na vreme koje imam u injectionImpulsi. Ovaj nacin
%interpolacije je jedino radio. 

fazaOut = interp1(tNovo,angle(X1),t1,'spline');
POut = interp1(tNovo,abs(X1),t1,'cubic');
POut = POut.^2;

% Y3 = interp1(tNovo,abs(X1),t1,'linear').*exp(1i*interp1(tNovo,angle(X1),t1,'linear'));
Y3 = sqrt(POut).*exp(1i*fazaOut);

% Y4 = circshift(Y3,length(t1)-13800);
% Y3 = interp1(tNovo,X,t1);

% Y3 = downsample(X,8);
% tNovo1 = downsample(tNovo,8);
figure
plot(t1,POut, "LineWidth",1.2)
hold on
plot(t1,PIn,"--", "LineWidth",1.2)
xlabel("Vreme [s]")
ylabel("|Ein|^2")
title("Snaga signala pre i nakon vlakna downsamlovan")
grid on
grid minor
% xlim([1.2e-8,1.3e-8])
legend('nakon vlakna', 'na ulazu u vlakno')
% save('nakonProlaska.mat','Y3','POut','fazaOut')

figure
plot(t1,angle(Y3), "LineWidth",1.2)
hold on
plot(t1,angle(EIn),"--", "LineWidth",1.2)
xlabel("Vreme [s]")
ylabel("angle(Ein)")
title("Faza signala pre i nakon vlakna")
grid on
grid minor
% xlim([0.99e-8,1.01e-8])
legend('nakon vlakna', 'na ulazu u vlakno')
% save('nakonProlaska.mat','X','tNovo')
