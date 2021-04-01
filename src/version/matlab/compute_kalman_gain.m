% v 1.0.2
function K = compute_kalman_gain(P_bar_xh, P_bar_hh)
  P_bar_hh_inv = inv(P_bar_hh);
  K = P_bar_xh * P_bar_hh_inv;
end