% v 1.0.2
function  Sigmas  = state_transition(func, Vargz, sigmaPoints, state_types, t_0, dt, t_1)
  options = odeset('RelTol',1e-5,'Stats','off','OutputFcn',@odeplot);
  
  if dt == 0
    prediction_time_points = [t_0 t_1];
  else
    prediction_time_points = t_0:dt:t_1;
  end
  
  % Get indexes
  state_index = find(strcmp(state_types, 'state') ==1);              % gets index of states
  parameter_index = find(strcmp(state_types, 'parameter') ==1);      % gets index of parameters
  constraint_index = find(strcmp(state_types, 'constraint') ==1);    % gets index of constraints
  
  
  Sigmas = zeros(size(sigmaPoints));
  parfor i = 1:size(sigmaPoints,1)
    % For 'pesudo states' for the parameter estimates

    [~, res] = ode23s(func, prediction_time_points, sigmaPoints(i,state_index)', options, [sigmaPoints(i,parameter_index) Vargz]);
    Sigmas(i,:) = [res(end, 1:length(res(state_index,end))) sigmaPoints(i,parameter_index)];
    % The code "[res(end, 1:length(res(state_index,end)))" accomidates pulling multiple state variable forcast from ODE solver
    %   'end' gets the final results (the row) of states, for the forcast
    %   '1:length(res(state_index,end)))' pulls the row out (1 to num of forcasted state variable
  end
  %Sigmas = sigmaPoints;
end