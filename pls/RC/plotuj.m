function a = plotuj(x,y,naziv,xosa,yosa)
f = figure('visible','off');
plot(x,y,'LineWidth',1.1)
title(naziv)
xlabel(xosa)
ylabel(yosa)
grid on 
grid minor
a = 'plot';
saveas(f,naziv,'jpg');

end