function res = linear_with_parameter__get_measurement(t, x, sol, vargz)

  tm = vargz{1}; % time series
  zm = vargz{2}; % measurement series

  index = find(tm ==t); % find index for measurment corresponding to time stamp
  
  if length(sol(1,:)) <= 1
    pseudo_meas = vargz{3};
  else
    pseudo_meas = ( (sol(1,end)-sol(1,end-1)) / (tm(index)-tm(index-1)) )  - (-2*sol(1,end) +2*sol(2,end));
  end
  
  res = [zm(:, index)' pseudo_meas];
  
end