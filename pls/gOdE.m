function [gainNkonc,centralnaEng,nth] = gOdE(gainN,photonenergy, nkoncrange,gth)

%za svaku koncentraciju nalazim maksimalan gain
maksGain = max(gainN,[],2);
indeksTh=find(maksGain>gth,1); %ovo mi je indeks nth zapravo

gainNEng = gainN(indeksTh,:); %ovo je g(eng, n=nth)
% Lambda = 1239.8./photonenergy; %prebacujem energiju u lambda

indeksEng = find(gainNEng == max(gainNEng));
centralnaEng = photonenergy(indeksEng); %nasla sam centralnu talasnu duzinu
nth = interp1(gainN(:,indeksEng), nkoncrange, gth);
gainNkonc = gainN(:,indeksEng); %g(lambda=lambdaCentralno, n)
%gth, neff, L
end