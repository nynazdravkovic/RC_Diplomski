function Pout = snagaImpulsiDrugiProlaz(P, tspan, t)
if t < 5e-9 
    Pout = 0;
else
    ind = find(tspan==t);    
    Pout = P(ind);
    disp(['t: ', num2str(t), ', ind: ', num2str(ind), ', Pout: ', num2str(Pout)]);

end     

end