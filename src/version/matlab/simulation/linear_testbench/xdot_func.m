% v1.0.3
function ret = func_f(t, y, vargz)
  %ret = [y(2)*vargz(1); 0.0];  
  ret = [y(2)*vargz(1); 0.0; y(4)*vargz(1); 0.0];
  %ret = [y(2)*vargz(1); 0.0; 0.0];
end