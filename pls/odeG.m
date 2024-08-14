function dydt = odeG(t,y, Rspsve, n,gainNkonc,I,Ab,Aw,Bb,Cb,Cw, Gama, Vtot,Vw, taup ,betasp,eps,etainj,q,taubw,tauwb,vg)
% function dydt = odeG(t,y, Rspsve, n,gainNkonc,I,Vqw, Vw, Vtot, Gamaqw, Gama, vg, K0, K1)

dydt=zeros(3,1);
%load("konstante.mat", "Ab","Aw","Bb","Cb","Cw", "Gama", "K0","K1","Vtot","Vw", "alphaI","betasp","eps","etainj","q","taubw","tauwb","vg");

%ako hocemo da sve bude fiksirano sem Rst komentarisi do njega:
gth = interp1(n, gainNkonc, y(2));
% disp(gth)

% Rsp = interp1(n, Rspsve, y(2),"linear","extrap");
Rsp = interp1(n, Rspsve, y(2));
taub = (Ab+Bb*y(1)+Cb*y(1)^2)^(-1);
tauw = (Aw+Rsp/y(2)+Cw*y(2)^2)^(-1);

dydt(1) = etainj*I/q/Vtot - y(1)/taub -y(1)/taubw + y(2)*Vw/tauwb/Vtot;
dydt(2) = y(1)*Vtot/taubw/Vw - y(2)/tauw - y(2)/tauwb - vg*gth*y(3)/(1+eps*y(3));
dydt(3) = Gama*vg*gth*y(3)/(1+eps*y(3)) + betasp*Gama*Rsp - 1/taup*y(3);


end