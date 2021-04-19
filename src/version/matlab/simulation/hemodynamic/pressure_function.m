%------------------------------------------------------------------
function dP = pressure_function(t,P, vargz)
  %==============================================================
  % Solves differential equation in the interval t with initial
  % guess y0
  %==============================================================
  
  R1 = vargz(1);
  R2 = vargz(2);
  C  = vargz(3);
  cycle_time = vargz(4);
  Q_dt = vargz(5);
  
  dP = (1/C * ( (1+ (R1/R2))*blood_flow(t,cycle_time) + C*R1*dQ(t,cycle_time,Q_dt) - (P / R2))));
  
end
%------------------------------------------------------------------