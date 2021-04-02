% v 1.0.2
function sim_meas = noisy_measurements(true_meas, StaD)
%==============================================================
% Adds Noise to the true measurements to simulate experimental
% measurments
%==============================================================
  
  % Pressure Measuments
  sim_meas = true_meas + sqrt(StaD)*randn(size(true_meas));
end