function dQ  =  dQ(t, cycle_time, Q_dt)
    
    dt = Q_dt;
    t = mod(t,cycle_time);
    dQ = (Q(t + dt, cycle_time) - Q(t, cycle_time)) / dt;
    
end