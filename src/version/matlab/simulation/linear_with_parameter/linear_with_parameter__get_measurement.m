function res = linear_with_parameter__get_measurement(t, x, vargz)

  tm = vargz(1, :); % time series
  zm = vargz([2 3], :); % measurement series

  index = find(tm ==t); % find index for measurment corresponding to time stamp
  
  res = [zm(:, index)' 1.00];
  
end