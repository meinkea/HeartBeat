function measurement = hemodynamic_state_to_measurement(state, Vargz_h)
  
  H = [1 0 0 0];
  
  measurement = H * state';
  
end