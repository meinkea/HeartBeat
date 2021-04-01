% v 1.0.2
function [y_true,zP,t] = Noisy_Measurements(func, Vargz, x, t_start, dt, t_end, StaD)
%==============================================================
% Adds Noise to the process model to simulate experimental
% measurments
%==============================================================
  options = odeset('RelTol',1e-5,'Stats','on','OutputFcn',@odeplot);
  
  % Time
  t = t_start:dt:t_end;
  
  % True Measurements
  [t, y_true] = ode45(func, t, x, options, Vargz);
  
  % Pressure Measuments
  zP = y_true + sqrt(StaD)*randn(size(y_true));
end