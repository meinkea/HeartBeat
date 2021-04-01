% v 1.0.2
function  Sigmas_h  = state_to_measurements(func, Sigmas, Vargz_h)
  Sigmas_h = func(Sigmas,  Vargz_h);
end