function res = hemodynamic__true_state_to_meas(x, vargz_h)
  
  H = [1 0 0 0];
  %H = [[1 0 0 0];...
  %     [0 1 0 0];...
  %     [0 0 1 0];...
  %     [0 0 0 1]];
  
  res = H * x';
  
end