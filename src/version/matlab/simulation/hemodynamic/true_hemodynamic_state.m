function [true_pressure_state, true_time_stamp] = true_hemodynamic_state(hemo_pressure, x0_init, t_sim_start, dt_sim, t_sim_end, vargz, cycles_to_skip)
  cycle_time = vargz(4);
  
  t_for_sim = (t_sim_start):dt_sim:(t_sim_end +cycles_to_skip*cycle_time);
  
  options = odeset('RelTol',1e-5,'Stats','on','OutputFcn',@odeplot);
  
  % ode23s is an implicit solver to solve the stiff equation, and it being used to replace the normal ode45
  [t_sim, p_sim] = ode23s(hemo_pressure, t_for_sim, x0_init, options, vargz);
  
  true_pressure_state = p_sim(1+ (cycle_time*cycles_to_skip)/dt_sim : end);
  
  t_sim = t_sim(1+ (cycle_time*cycles_to_skip)/dt_sim : end);
  
  true_time_stamp = t_sim - cycle_time*cycles_to_skip + dt_sim;
  
  for i = 1:length(true_pressure_state)
    true_pressure_state(i,2) = vargz(1);
    true_pressure_state(i,3) = vargz(2);
    true_pressure_state(i,4) = vargz(3);
  end

end