% v 1.0.2
function  Sigmas  = stateTransition(func, Vargz, sigmaPoints, t_0, dt, t_1)
  options = odeset('RelTol',1e-5,'Stats','off','OutputFcn',@odeplot);
  prediction_time_points = t_0:dt:t_1;
  
  Sigmas = zeros(size(sigmaPoints));
  parfor i = 1:size(sigmaPoints,1)
    [~, res] = ode45(func, prediction_time_points, sigmaPoints(i,:)', options, Vargz);
    Sigmas(i,:) = res(end,:)';
  end
  %Sigmas = sigmaPoints;
end