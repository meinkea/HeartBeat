function res = linear_with_parameter_sys(t, y, vargz)

  res(1) = -2*y(1) +2*y(2) +vargz(1);
  res(2) = -4*y(1) +2*y(2);
  
  res = res';

end