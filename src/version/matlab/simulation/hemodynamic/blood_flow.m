function res = blood_flow(t, cycle_time)
    t = mod(t,cycle_time);
    S = 0.070;
    res = (1.4/(S.*sqrt(pi))).*exp(-((t - 0.36).^2)/(2.*S^2)) + 5;
end


