function P = snaga(t, brStep, dt)
% p1 = -20; %pocetni stepenik
% p2 = -10; %najvisi stepenik
% dp = -0.2;
% dP = 10^(dp/10)/1e3;

Pstart = 0.0001e-3;
Pend = 0.5e-3;
% Pstart = 0.1e-3;
% Pend = 10e-3;
top_step = (brStep + 1) / 2;
t_step = ceil(t / dt);

dP = (Pend-Pstart)/((brStep + 1) / 2);

if (t_step <= top_step)
    P = Pstart + t_step * dP;
else
    P = Pstart + dP * top_step - (t_step - top_step) * dP;
end

end