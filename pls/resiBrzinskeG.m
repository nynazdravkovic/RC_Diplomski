function [Ssr, nwsr, nbsr] = resiBrzinskeG(Rspsve, nsve,gainNkonc,I, tspan, y0, Ab,Aw,Bb,Cb,Cw, Gama,Vtot,Vw, taup,betasp,eps,etainj,q,taubw,tauwb,vg)
%resiBrzinske resava brzinske jednačine iz fajla odeG i vraća finalne
%koncentracije u stacionarnom stanju 
[t,y] = ode45(@(t,y) odeG(t,y,Rspsve, nsve,gainNkonc,I*1e-3,Ab,Aw,Bb,Cb,Cw, Gama, Vtot,Vw, taup,betasp,eps,etainj,q,taubw,tauwb,vg), tspan, y0);
s = y(:,3);
nw = y(:,2);
nb = y(:,1);
Ssr = s(end);
nwsr = nw(end);
nbsr = nb(end);
end