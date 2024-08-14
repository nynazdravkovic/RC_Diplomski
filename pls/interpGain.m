function [gain1] = interpGain(gain,n, enegry, hni)
%vraca gain od koncentracije za trazeno hni

sz = size(gain(:,100));

gain1 = zeros(sz);
    for i = 1:length(n)
        gain1(i) =  interp1(enegry, gain(i,:), hni);
    end  
end