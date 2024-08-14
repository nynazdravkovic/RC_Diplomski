function dydt = odeInjLocking2(t,y, Rspsve, n, I,Ab,Aw,Bb,Cb,Cw, Gama, Vtot, Vw, betasp, eps, etainj, q, taubw, tauwb, vg, alpha, m, kc, Sinj, df, alphaI,  ind, brMod, fitovaniGain)

%m je mod u koji injektujemo a fsr je intermodalni prostor
% 
% brMod = 20; %sa jedne strane ako je centralni mod indeksa 0
% i = (-brMod:brMod);
% % frekvencije = q*energy/h; %sve frekvencije koje imam u fajlu 
% modovi = f0 + i*fsr; %modovi od -20 do 20 
k = isreal(y(2));

if k ~= 1
    display(k)
end

dydt=zeros(length(ind)+3,1); 
y(2) = abs(y(2));
% alphaI = K0 + Gama * K1 * abs(y(2));
Rsp = interp1(n, Rspsve, y(2));
taub = (Ab+Bb*y(1)+Cb*y(1)^2)^(-1);
tauw = (Aw+Rsp/y(2)+Cw*y(2)^2)^(-1);
g = zeros(size(ind));

suma = 0;
suma2 = sum(y(3:length(ind)+2));

for j = (1:length(ind))
    g(j) = fitovaniGain(j,1)*y(2)^4 + fitovaniGain(j,2)*y(2)^3 + fitovaniGain(j,3)*y(2)^2 + fitovaniGain(j,4)*y(2)^1 + fitovaniGain(j,5);
    suma = suma + g(j)*y(2+j)/(1+eps*suma2); %ovo je za sumu u dnwdt
end

dydt(1) = etainj*I/q/Vtot - y(1)/taub -y(1)/taubw + y(2)*Vw/tauwb/Vtot;

dydt(2) = y(1)*Vtot/taubw/Vw - y(2)/tauw - y(2)/tauwb - vg*suma;
%dodati sumu svuda i p(t)
for j = (1:length(ind))
    if j~=brMod+1+m
        dydt(2+j) = Gama*vg*g(j)*y(2+j)/(1+eps*suma2) + betasp*Gama*Rsp - vg*alphaI*y(2+j);
    else
        dydt(2+j) = Gama*vg*g(j)*y(2+j)/(1+eps*suma2) + betasp*Gama*Rsp - vg*alphaI*y(2+j)+ 2*kc*sqrt(Sinj*y(2+j))*cos(y(3+length(ind)));
    end 
end

dydt(3+length(ind)) = Gama*vg*g(brMod+1+m)*alpha/(1+eps*suma2)/2 - alpha/2*alphaI*vg - 1*pi*df - kc*sqrt(Sinj/y(brMod+1+m))*sin(y(3+length(ind)));

end