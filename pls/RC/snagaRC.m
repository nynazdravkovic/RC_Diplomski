function P = snagaRC(input, t, B, fwhm)
%t je jedan vremenski korak

o = fwhm/sqrt(8*log(2));
mi = 1/B*(1:length(input))+1/B*0.5;
P = 0;
for i = 1:length(mi)
    P = P + exp(-(t-mi(i)).^2/o^2)*input(i)*1e-3;
end

end