function [tm_true_stamp, true_measurement] = true_hemodynamic_measurement(true_state_to_measurement, true_states, tm_stamp)
  
    
  true_measurement = true_state_to_measurement(true_states, 0.0);
  
  tm_true_stamp  = tm_stamp;
  
end