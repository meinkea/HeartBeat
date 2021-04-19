%------------------------------------------------------------------
function p = pressure_function(t,p_init, vargz)
  %==============================================================
  % Solves differential equation in the interval t with initial
  % guess y0
  %==============================================================
  
  R1 = vargz(1);
  R2 = vargz(2);
  C  = vargz(3);
  cycle_time = vargz(4);
  Q_dt = vargz(5);
  
  p = (1/C * ( blood_flow(t,cycle_time) + C*dQ(t,cycle_time,Q_dt) - (p_init / R2) + (blood_flow(t,cycle_time)*(R1/R2))));
  
end
%------------------------------------------------------------------