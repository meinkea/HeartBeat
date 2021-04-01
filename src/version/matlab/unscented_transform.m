% v 1.0.2
function [x_bar, P_bar] = unscentedTransform(Wm, Wc, Sigmas, P_design)
  
  x_bar = Wm * Sigmas;
  
  P_bar = cross_covariance(Wc, x_bar, Sigmas, x_bar, Sigmas);
  
  P_bar = P_bar + P_design;
  
end