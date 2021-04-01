% v 1.0.2
function  Sigmas_h  = state_to_measurement(func, Sigmas, Vargz_h)
  Sigmas_h = func(Sigmas,  Vargz_h);
end