% v 1.0.2
function [P_bar_ab] = cross_covariance(Wc, x_a, sigmas_a,  x_b, sigmas_b)
  na = length(x_a);
  nb = length(x_b);
  
  NA = length(sigmas_a(:,1));
  NB = length(sigmas_b(:,1));

  if NA ~= NB
    ME = MException('ERR:: sigmas_a and sigmas_b Dimensions do not match!');
    throw(ME)
  end
  
  
  P_bar_ab = zeros(na,nb);
  
  parfor i = 1:NA;
    da = sigmas_a(i,:) - x_a;
    db = sigmas_b(i,:) - x_b;
    
    cross_covariance_at_point = (da' * db);
    
    P_bar_ab = P_bar_ab + Wc(i) * cross_covariance_at_point;
  end
  
end