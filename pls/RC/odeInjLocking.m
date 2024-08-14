function dydt = odeInjLocking(t,y, Rspsve, n, I,Ab,Aw,Bb,Cb,Cw, Gama, Vtot, Vw, betasp, eps, etainj, q, taubw, tauwb, vg, alpha, m, kc, df, alphaI,  ind, brMod, fitovaniGain,P, injektuj, eta0, h, lambdainj, taup, c0)

if(injektuj ==1)
    Sinj = P/(eta0 * h * c0 * 1e2 / lambdainj * 1e7 * Vw / taup / Gama); %m/s pretvaram u cm/s i nm u ncm
else
    Sinj = 0; 
    df = 0;
end
%gledam da li imam injekciju u ovom slucaju ili ne

k = isreal(y(2));
if k ~= 1
    display(k)
end

dydt=zeros(length(ind)+3,1); 
Rsp = interp1(n, Rspsve, abs(y(2)));

taub = (Ab+Bb*y(1)+Cb*y(1)^2)^(-1);
tauw = (Aw+Rsp/y(2)+Cw*y(2)^2)^(-1);

g = zeros(size(ind));

suma = 0;
suma2 = sum(y(3:(length(ind)+2)));

for j = (1:length(ind))
    % g(j) = fitovaniGain(j,1)*y(2)^6 + fitovaniGain(j,2)*y(2)^5 + fitovaniGain(j,3)*y(2)^4 + fitovaniGain(j,4)*y(2)^3 + fitovaniGain(j,5)*y(2)^2 + fitovaniGain(j,6)*y(2) + fitovaniGain(j,7);
    g(j) = fitovaniGain(j,1)*y(2)^4 + fitovaniGain(j,2)*y(2)^3 + fitovaniGain(j,3)*y(2)^2 + fitovaniGain(j,4)*y(2) + fitovaniGain(j,5);
    % g(j) = fitovaniGain(j,1)*y(2)^9 + fitovaniGain(j,2)*y(2)^8 + fitovaniGain(j,3)*y(2)^7 + fitovaniGain(j,4)*y(2)^6 + fitovaniGain(j,5)*y(2)^5 + fitovaniGain(j,6)*y(2)^4 + fitovaniGain(j,7)*y(2)^3 +fitovaniGain(j,8)*y(2)^2 +fitovaniGain(j,9)*y(2) +fitovaniGain(j,10);

    suma = suma + g(j)*y(2+j)/(1+eps*suma2); %ovo je za sumu u dnwdt
end

dydt(1) = etainj*I/q/Vtot - y(1)/taub -y(1)/taubw + y(2)*Vw/tauwb/Vtot;

dydt(2) = y(1)*Vtot/taubw/Vw - y(2)/tauw - y(2)/tauwb - vg*suma;

for j = (1:length(ind))
    if j~=brMod+1+m
        dydt(2+j) = Gama*vg*g(j)*y(2+j)/(1+eps*suma2) + betasp*Gama*Rsp - 1/taup*y(2+j);
    elseif j==brMod+1+m
        % dydt(2+j) = Gama*vg*g(j)*y(2+j)/(1+eps*suma2) + betasp*Gama*Rsp - 1/taup*y(2+j)+ 2*kc*sqrt(Sinj*y(2+j))*cos(y(3+length(ind)))*ceil(heaviside(t-tau));
        dydt(2+j) = Gama*vg*g(j)*y(2+j)/(1+eps*suma2) + betasp*Gama*Rsp - 1/taup*y(2+j)+ 2*kc*sqrt(Sinj*y(2+j))*cos(y(3+length(ind)));

    end 
end

% dydt(3+length(ind)) = Gama*vg*g(brMod+1+m)*alpha/(1+eps*suma2)/2 - alpha/2*1/taup - 2*pi*df - kc*sqrt(Sinj/y(brMod+1+m+2))*sin(y(3+length(ind)))*ceil(heaviside(t-tau));
dydt(3+length(ind)) = Gama*vg*g(brMod+1+m)*alpha/(1+eps*suma2)/2 - alpha/2*1/taup - 2*pi*df - kc*sqrt(Sinj/y(brMod+1+m+2))*sin(y(3+length(ind)));

end