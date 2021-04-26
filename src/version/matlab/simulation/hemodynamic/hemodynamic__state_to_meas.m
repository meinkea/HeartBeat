function [res, res_types] = hemodynamic__state_to_meas(x, x_types, vargz_h)
  
  H = [1 0 0 0];
%   H = [[1 0 0 0];...
%        [0 0 0 1]];
  
  res = H * x';
%   res_types = {x_types{[1, 4]}};
  res_types = {x_types{[1]}};
  
end