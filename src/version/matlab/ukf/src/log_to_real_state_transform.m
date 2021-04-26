function x_real = log_to_real_state_transform(x_log, x_probability_distribution_types)
  
  lognormal_index = find(strcmp(x_probability_distribution_types, 'lognormal') ==1);
  
  x_real = x_log;
  
  for I = 1:length(x_real(:,1))
    x_real(I, lognormal_index) = exp(x_log(I, lognormal_index));
  end
  
end