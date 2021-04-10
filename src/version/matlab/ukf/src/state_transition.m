% v 1.0.2
function  Sigmas  = state_transition(func, Vargz, sigmaPoints, t_0, dt, t_1)
  options = odeset('RelTol',1e-5,'Stats','off','OutputFcn',@odeplot);
  prediction_time_points = t_0:dt:t_1;
  
  Sigmas = zeros(size(sigmaPoints));
  parfor i = 1:size(sigmaPoints,1)
    % For 'pesudo states' for the parameter estimates
    
    [~, res] = ode23s(func, prediction_time_points, sigmaPoints(i,1), options, [sigmaPoints(i,2) sigmaPoints(i,3) sigmaPoints(i,4) Vargz(4) Vargz(5)]);
    Sigmas(i,:) = [res(end) sigmaPoints(i,2) sigmaPoints(i,3) sigmaPoints(i,4)];
  end
  %Sigmas = sigmaPoints;
end