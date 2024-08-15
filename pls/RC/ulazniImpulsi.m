function P = ulazniImpulsi(i, K, tspan,M1_full, B, fwhm, dataset, A, num_avg, PSve, fazaSve, finj)
% Ova funkcija napravi gausove impulse u vremenu tspan za odgovarajuci podatak u datasetu
% (SantaFe_laser_dataset_norm(i)) pomnozeno sa maskom, odn rezultat su
% gausovi impulsi sirine fwhm, brzine B, u vremenu tspan koje krece od 0!!
% num_avg mi definise koju masku uzimam od onih 100 koje sluze za srednju
% vrednost
w0 = 2*pi*finj;

J_input = A*M1_full(1,num_avg,:)*dataset(i); %ovo mi je visina impulsa za inpute celog dataseta
Pin = zeros(1,length(tspan));
for k = 1:length(tspan)
    Pin(k) = snagaRC(J_input,tspan(k),B,fwhm);
end
%sabiram impulse sa onima koji su izasli iz nelinearnosti 
disp(i)
if i == 1
    E = sqrt(Pin).*exp(1j*w0*tspan); 
    P = abs(E).^2; 
    disp('usao')
else
    Ein = sqrt(PSve(i-1,:)).*exp(1j*(fazaSve(i-1,:)+tspan*w0))*sqrt(K);
    E = sqrt(Pin).*exp(1j*w0*tspan) + Ein; 
    P = abs(E).^2;
    % fi = angle(E);
    % E = sqrt(Pin).*exp(1j*w0*tspan) + K*sqrt(PSve(i-1,:)).*exp(1j*w0*tspan+fazaSve(i-1,:)); 
    % P = abs(E).^2; 
end

end

