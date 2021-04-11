function res = linear_with_parameter_meas_to_state(y, vargz)

  H = [[1 0 0];...
       [0 1 0]];
     
  res = H*y';

end
