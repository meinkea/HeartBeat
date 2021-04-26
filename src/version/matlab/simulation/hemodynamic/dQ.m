function dQ  =  dQ(t, cycle_time, Q_dt)
    
    dt = Q_dt;
    t = mod(t,cycle_time);
    dQ = (blood_flow(t + dt, cycle_time) - blood_flow(t, cycle_time)) / dt;
    
end