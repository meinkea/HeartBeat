% v 1.0.2
function sigmaPoints = compute_sigma_points(n, lambda, x, P)
  U = chol((n + lambda)*P);
  %U = U'; % U was wrong transpose version
    
  sigmaPoints = zeros(2*n+1, n);
    
  sigmaPoints(1,:) = x';
  for i = 1:n
    sigmaPoints(i+1,:) = x' + U(i,:);
    sigmaPoints(n + i + 1,:) = x' - U(i,:);
  end
end