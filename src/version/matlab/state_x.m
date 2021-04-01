% v 1.0.2
function x = state_x(mean,K,y)
  x = mean' + K*y';
end