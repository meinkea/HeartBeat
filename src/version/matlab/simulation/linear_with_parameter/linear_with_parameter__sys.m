function res = linear_with_parameter__sys(t, y, vargz)

  a = vargz(1);
  
  res(1) = -2*y(1) +2*y(2) +a;
  res(2) = -4*y(1) +2*y(2);
  
  res = res';

end