function dydt = odefun2(t,y,tau, f, Ith)
taub=2.5e-9;
tauw=1.5e-9;
taubw=14.5e-12;
tauwb=70e-12;

omega=1.5e-15;
taup=2.3e-12;
eps = 2e-17;
gama=0.08;
n0=4.2e18;
Rsp=1.2e12;
alpha=3.5;
q=1.6e-19;
vg=8.5e9;

niinj=0.85;
Vtot=1.03e-10;
Vw=1.68e-11;

dydt=zeros(4,1);
dydt(1) = niinj/q/Vtot*4*(Ith + heaviside(t-tau)*Ith*sin(2*pi*f*t))-y(1)/taub-y(1)/taubw+y(2)*Vw/tauwb/Vtot;
dydt(2) = y(1)*Vtot/taubw/Vw-y(2)/tauw-y(2)/tauwb-vg*omega*(y(2)-n0)*y(3)/(1+eps*y(3));
dydt(3) = gama*vg*omega*(y(2)-n0)*y(3)/(1+eps*y(3))-y(3)/taup+gama*Rsp/Vtot;
dydt(4) = 1/2*alpha*(gama*vg*omega*(y(2)-n0)/(1+eps*y(3))-1/taup);

end