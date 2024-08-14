function [f, mod] = modKrive(nb, nw, S, nsve, gainNkonc, Rspsve, Vtot, Ab, Aw, Bb, Cb, Cw, taubw, tauwb, Vw, Vb, vg, eps, Gama, betasp, etainj, q, taup)
%funkcija koja mi daje 20*log10(abs(dS/dI)) a u argumentima su joj stac
%vrednosti koje su dobijene kao resenje brzinskih jednacina koriscenjem
%neke struje
g = interp1(nsve, gainNkonc, nw);
dRspdn = diff(Rspsve)./diff(nsve.');
posl = dRspdn(end);
dRspdn(end+1)=posl;
dRspdnsr = interp1(nsve, dRspdn, nw);

dgdn = diff(gainNkonc)./diff(nsve.');
posl = dgdn(end);
dgdn(end+1)=posl;

dgdnsr = interp1(nsve, dgdn, nw);
f = linspace(5e6, 20e9, 1000);
ksi = 1i*2*pi*f;
mod = zeros(1,length(ksi));
for i = 1:length(ksi)
    hi11 = -Ab-2*nb*Bb-3*nb^2*Cb-1/taubw - ksi(i);
    hi12 = Vw/tauwb/Vb;
    hi13 = 0;

    hi21 = Vb/taubw/Vw;
    hi22 = -Aw-dRspdnsr-3*Cw*nw^2-1/tauwb-vg*S/(1+eps*S)*dgdnsr - ksi(i);
    hi23 = -vg*g*(1-S*eps/(1+S*eps))/(1+S*eps);

    hi31 = 0;
    hi32 = Gama*vg*dgdnsr*S/(1+eps*S)+Gama*betasp*dRspdnsr;
    hi33 = Gama*vg*g*(1/(1+eps*S)-S*eps/(1+eps*S)^2) - 1/taup - ksi(i);

    hi = [[hi11, hi12, hi13]; [hi21, hi22, hi23];[hi31, hi32, hi33]];
    res = [-etainj/q/Vtot;0;0];
    %x mi je matruca [dnb/dI, dnw/dI, dS/dI]
    a = linsolve(hi,res);
    mod(i) = 20*log10(abs(a(3)));
end
mod = mod - mod(1);
end