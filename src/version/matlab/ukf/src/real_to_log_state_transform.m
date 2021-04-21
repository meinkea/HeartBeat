function x_log = real_to_log_state_transform(x_real, x_probability_distribution_types)
  
  lognormal_index = find(strcmp(x_probability_distribution_types, 'lognormal') ==1);
  
  x_log = x_real;
  
  ln = @log;  % why is natural log not ln???
  for I = 1:length(x_real(:,1))
    x_log(I, lognormal_index) = ln(x_real(I, lognormal_index));
  end
  
end