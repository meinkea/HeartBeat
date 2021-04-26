function [res, res_types] = linear_with_two_parameter__meas_to_state(y, y_types, vargz)
  
  H = [[1 0 0 0];...
       [0 1 0 0];...
       [0 0 1 0]];
  
  res = H*y';
  res_types = {y_types{[1,2,3]}};
  
end
