%------------------------------------------------------------------
function p = pressure_function(t,p_init, vargz)
  %==============================================================
  % Solves differential equation in the interval t with initial
  % guess y0
  %==============================================================
  
  Rc = vargz(1);
  Rd = vargz(2);
  C  = vargz(3);
  cycle_time = vargz(4);
  Q_dt = vargz(5);
  
  p = (1/C * ( Q(t,cycle_time) + C*dQ(t,cycle_time,Q_dt) - (p_init / Rd) + (Q(t,cycle_time)*(Rc/Rd))));
  
end
%------------------------------------------------------------------