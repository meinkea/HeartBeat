function res = hemodynamic__get_measurement(t, x, sol, vargz)
  
  tm = vargz{1}; % time series
  zm = vargz{2}; % measurement series
  
  index = find(tm ==t); % find index for measurment corresponding to time stamp
  
  func = vargz{3}; % state function
  
  if length(sol(1,:)) <= 1
    pseudo_meas_0 = vargz{4};
    %pseudo_meas_1 = vargz{5};
    %pseudo_meas_2 = vargz{6};
  else
    pseudo_meas_0 =  sol(4,end) * (func(t, sol(1,end), [sol(2,end) sol(3,end) sol(4,end) vargz{7} vargz{8}]))  /  ( (sol(1,end)-sol(1,end-1)) / (tm(index)-tm(index-1)) ) ;
  end
  
  res = [zm(:, index)' pseudo_meas_0];
  
end