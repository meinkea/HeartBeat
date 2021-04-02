% v 1.0.2
function Sigmas_h = h_func(x,  Vargz_h)
%  H = [[1 0];...
%       [0 1]];
  
  H = [[1 0 0 0];...
       [0 0 1 0]];
  
  Sigmas_h = H *(x');
  Sigmas_h = Sigmas_h';
  
  
end