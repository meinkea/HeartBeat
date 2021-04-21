function [res, res_probability_distribution_type] = hemodynamic__get_measurement(t, x, sol, vargz)
  
  tm = vargz{1}; % time series
  zm = vargz{2}; % measurement series
  
  index = find(tm ==t); % find index for measurment corresponding to time stamp
  
  func = vargz{3}; % state function
  
  
  res = [zm(:, index)'];
  
  res_probability_distribution_type = {'lognormal'};
  
  
%   if length(sol(1,:)) <= 5
%     %pseudo_meas_0 = vargz{4};
%     %pseudo_meas_1 = vargz{5};
%     pseudo_meas_2 = vargz{6};
%   else
%     t0 = tm(index-1);
%     
%     P = sol(1,end);
%     R1 = sol(2,end);
%     R2 = sol(3,end);
%     C = sol(4,end);
%     
%     P_dot = pressure_function(t0, P, [R1 R2 C 0.8 0.00001]);
%     
%     Q = blood_flow(t0, 0.8);
%     Q_dot = dQ(t0, 0.8, 0.00001);
%     
%     pseudo_meas_2 =  ( (1+R1/R2)*Q - (P/R2) ) / (P_dot - C*R1*Q_dot );
%     
%     if pseudo_meas_2 < 0.0000001
%       pseudo_meas_2 = 0.0000001;
%     end
%   end
%   
%   res = [zm(:, index)' pseudo_meas_2];
%   
%   res_probability_distribution_type = {'lognormal' 'lognormal'};
  
end