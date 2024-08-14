close all 
clc


tspan = [0 25e-9];
y0 = [0,0,0,0];
tau = 5e-9;
f=6e9;
Ith = 0.015;

[t,y] = ode45(@(t,y) odefun2(t,y,tau,f,Ith), tspan, y0);

figure(1)
plot(t,y(:,3))
grid on
grid minor
title('Koncentracija fotona sa modulisanom strujom');
xlabel('Vreme');

%vreme kada modulisemo struju
indeksi = find(t>(tau+1/f));
t1 = transpose(t(t>(tau+1/f)));
fotoni = transpose(y(indeksi(1):end,3));

figure()
plot(t1,fotoni)
grid on
grid minor
title('Koncentracija fotona sa modulisanom strujom');
xlabel('Vreme');

faza = transpose(y(indeksi(1):end,4));
figure()
plot(t1,faza)
grid on
grid minor
title('Faza sa modulisanom strujom');
xlabel('Vreme');

% figure(2)
% plot(t,y(:,1))
% grid on
% grid minor
% title('Koncentracije nosilaca u barijeri');
% xlabel('Vreme');
% 
% figure(3)
% plot(t,y(:,3))
% grid on
% grid minor
% title('Evolucija faze');
% xlabel('Vreme');
