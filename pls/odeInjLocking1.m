function dydt = odeInjLocking1(t,y, Rspsve, n ,I,Ab,Aw,Bb,Cb,Cw, Gama, Vtot, Vw, betasp, eps, etainj, q, taubw, tauwb, vg, alpha, m, alphaI, ind, brMod,fitovaniGain)

dydt=zeros(length(ind)+3,1);

% alphaI = K0 + Gama * K1 * abs(y(2));
Rsp = interp1(n, Rspsve, abs(y(2)));
taub = (Ab+Bb*y(1)+Cb*y(1)^2)^(-1);
tauw = (Aw+Rsp/y(2)+Cw*y(2)^2)^(-1);
g = zeros(size(ind));
suma = 0;

for j = (1:length(ind))
    g(j) = fitovaniGain(j,1)*y(2)^4 + fitovaniGain(j,2)*y(2)^3 + fitovaniGain(j,3)*y(2)^2 + fitovaniGain(j,4)*y(2)^1 + fitovaniGain(j,5);
    suma = suma + g(j)*y(2+j)/(1+eps*y(2+j)); %ovo je za sumu u dnwdt
end


dydt(1) = etainj*I/q/Vtot - y(1)/taub -y(1)/taubw + y(2)*Vw/tauwb/Vtot;

dydt(2) = y(1)*Vtot/taubw/Vw - y(2)/tauw - y(2)/tauwb - vg*suma;

for j = (1:length(ind))
    dydt(2+j) = Gama*vg*g(j)*y(2+j)/(1+eps*(y(2+j))) + betasp*Gama*Rsp - vg*alphaI*y(2+j);
end

dydt(3+length(ind)) = Gama*vg*g(brMod+1+m)*alpha/(1+eps*(y(brMod+1+m+2)))/2 - alpha/2*alphaI*vg;
