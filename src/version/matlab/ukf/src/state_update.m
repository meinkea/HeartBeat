% v 1.0.2
function x = state_update(mean,K,y)
  x = mean' + K*y';
end