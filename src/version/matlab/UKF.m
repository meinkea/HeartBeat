% v 1.5.0

clear all
close all
%==========================================================================
%=============================SET Constants================================
%==========================================================================

addpath(genpath('./ukf/')); % import from sub directories (installs the toolkit)

%addpath('./simulation/hemodynamic/');
addpath('./simulation/linear_with_parameter/');
%addpath('./user_design/linear_testbench/');
addpath('./user_design/hemodynamic/');

%==============================Plots=======================================

Test_plots = 1; %flag
profiling_flag = 0;

%%
%============================== SIMULATION ==============================
% ------------------------- Simulation Settings -------------------------
% These are
%x0_init = 8.777e+04; % True intial start
% c = 0.000750061683; % unit conversion 1 dyn/cm^2=0.0007750061 mmHg
% x0_init = 80/c;
% 
% meas_variance = 0.5*10000000;
% 
% Cylces = 3;
% 
% % Measurement Time
% t_sim_start = 0.0;         % start time
% t_sim_end = Cylces * 0.8;  % end time
% dt_sim = 1/500;            % time between measurements



% -------------------------------- Sim --------------------------------



% Rc_true = 1600;  % Silly
% Rd_true = 13000; %
% C_true = 0.000025; %
% 
% cycle_time = 0.8;
% Q_dt = 1e-4;
% 
% cycles_to_skip = 4;
% 
% % additional arguments for hemo_pressure (alias of pressure_function)
% true_hemodynamic_state_argz = [Rc_true, Rd_true, C_true, cycle_time, Q_dt];
% 
% % Generate true states
% [true_pressure_state, true_time_stamp] = true_hemodynamic_state(@pressure_function,  x0_init,  t_sim_start, dt_sim, t_sim_end,  true_hemodynamic_state_argz, cycles_to_skip);
% 
% % Generate true measurements.  Returns also the new time series that accounts
% [tm_true_stamp, true_measurement] = true_hemodynamic_measurement(@hemodynamic_state_to_measurement, true_pressure_state, true_time_stamp);
% 
% -- Sim measurements
% zm = noisy_measurements(true_measurement, meas_variance);
% 
% -- Sim runtime
% tm = tm_true_stamp;


% Measurement Time

t_sim_start = 0.0;         % start time
t_sim_end = 4.0;  % end time
dt_sim = 0.01;            % time between measurements

tm = t_sim_start:dt_sim:t_sim_end;
x_init = [1 2]

true_meas(1,:) = (3*sin(2*tm)+3*cos(2*tm)-1)/2;
true_meas(2,:) = 3*cos(2*tm)-1;

zm = noisy_measurements(true_meas, 0.0);


%%
% ======================= SYSTEM DESIGN from USER =======================
%Pref  = 80/c;      %mmHg
%Rcref = 1600;
%Rdref = 13000;
%Cref  = 0.000025;
% 
% 
% % -- State Vector
% x_user = { 80/c  'state'    ;...
%            Rcref 'parameter';...
%            Rdref 'parameter';...
%            Cref  'parameter'};
% 
% % -- State Estimator Variance and Covariance Matrix
% P_user = [[ 8.0  0.0  0.0  0.0];...
%           [ 0.0  8.0  0.0  0.0];...
%           [ 0.0  0.0  8.0  0.0];...
%           [ 0.0  0.0  0.0  8.0]];
%   
% % -- Process model
% process_model_function_user = @pressure_function;
% 
% % UKF Process model prediction Time
% t_last_update = t_sim_start; % time at last update
% t_prior = 0.0;               % time at prior, this is updated in loop
% dtp = 0;                     % time per step in sigma predecition projection
% 
% 
% dt_user = dtp;
% 
% % -- Process Model Variance and Covariance Matrix
% a_ = 0.001;
% f_ = 0.01;
% Q_user = ([[  a_,   0.0,   0.0,   0.0];...
%            [ 0.0, f_*a_,   0.0,   0.0];...
%            [ 0.0,   0.0, f_*a_,   0.0];...
%            [ 0.0,   0.0,   0.0,  f_*a_]]);
% 
% % State Function Model for State Estimator
% state_to_measurment_function = @hemodynamic_state_to_measurement;
% 
% % State Function Model parameters (vargz for state_transition())
% state_transition_Vargz = [cycle_time Q_dt];
% 
% % -- Measurement Variance and Covariance Matrix
% R_user = [0.005]; % For X meas only
% 
% 
% % --- Sigma Points tuning parameters ---
% 
% alpha =  0.1;
% beta  =  2.0;
% kappa = -1.0; % prolly should be set to n-3, so -1


% -- State Vector
x_user = {1.0 'state';...
          2.0 'state';...
          1.0 'parameter'};
        
% -- State Estimator Variance and Covariance Matrix
P_user = [[ 0.001  0.0     0.0        ];...
          [ 0.0    0.001  0.0        ];...
          [ 0.0    0.0     0.00000001]];...

% -- Process model
process_model_function_user = @linear_with_parameter_sys;

% UKF Process model prediction Time
t_last_update = tm(1); % time at last update
t_prior = 0.0;               % time at prior, this is updated in loop
dtp = 0;                     % time per step in sigma predecition projection

% -- Process Model Variance and Covariance Matrix
Q_user = [[ 1.0  0.0   0.0   ];...
          [ 0.0   1.0  0.0   ];...
          [ 0.0   0.0   0.001 ]];...

% State Function Model parameters (vargz for state_transition())
state_transition_Vargz = [0.0];

% -- Measurement Variance and Covariance Matrix
R_user = [[0.01 0.00];...
          [0.00 0.01]]*0.00001;


% --- Sigma Points tuning parameters ---
alpha =  0.1;
beta  =  2.0;
kappa = -1.0; % prolly should be set to n-3, so -1

%%
%==============Making Measuments and comparing to True funtion=============

%setting seed
rng('default')


if Test_plots == 1
    figure
    %plot(x_True(1),x_True(3),'o')
    hold on
    %plot(y_True, x_True, 'k', 'LineWidth', 3)
    plot(tm, zm, 'k', 'LineWidth', 3)
    grid on
    xlabel('Time (s)')
    ylabel('Pressure (log world)')
    shg
end


%==========================Initialize Matrices=============================
x = [x_user{:,1}]';
x_types = {x_user{:,2}};
P = P_user;
Q = Q_user;
R = R_user;

process_model_function = process_model_function_user;


% number of state vectors
n = length( x );

% number of measured observations
nmo = length(R);

[sigmaPoints,weights,lambda,P_bar_xx,P_bar_hh] = initialize(n,alpha,kappa,nmo); % this should be broken up
sol = x;



% These are the weigths used in the uncented transform
% They are static and never change
[Wc,Wm] =  compute_weights(n,lambda,alpha,beta);



% for compa


for i = 1:(length(tm)-1)
    if profiling_flag == 1
      tic
    end
    % --- PREDICT ---
    % Get time for prior ( which is the time of the next measurement )
    t_prior = tm(i+1); % this needs to be decoupled from sim
    
    
    sigmaPoints = compute_sigma_points(n, lambda, x, P); % Checked
    if profiling_flag == 1
      disp('sigmaPoints')
      toc
    end

    if profiling_flag == 1
      tic
    end
    % Project sigma points using ODE solver
    Sigmas_f = state_transition(process_model_function, state_transition_Vargz, sigmaPoints, x_types, t_last_update, dtp, t_prior); % Checked
    if profiling_flag == 1
      disp('state_transition')
      toc
    end
    
    if profiling_flag == 1
      tic
    end
    % _bar denotes that it came from unscentedTransform
    [ x_bar_xx, P_bar_xx ] =  unscented_transform(Wm,Wc, Sigmas_f, Q); % Checked
    if profiling_flag == 1
      disp('unscented_transform')
      toc
    end
    
    if profiling_flag == 1
      tic
    end
    % Recompute NEW sigmas at the new predicted state
    Sigmas_x = compute_sigma_points(n,lambda,x_bar_xx',P_bar_xx);
    if profiling_flag == 1
      disp('compute_sigma_points')
      toc
    end
    % --- UPDATE ---
    
    if profiling_flag == 1
      tic
    end
    % _h denotes measurement version of prediction
    % >>> !!! THIS WAS A BUG THAT WAS NEVER FIXED !!! <<<
    %Sigmas_h = state_to_measurment(hemodynamic_state_to_measurement, Sigmas_x, 0);
    %Sigmas_h = hemodynamic_state_to_measurement(Sigmas_x,  0.0);
    Sigmas_h = linear_with_parameter_meas_to_state(Sigmas_x, 0.0);
    
    
    Sigmas_h = Sigmas_h';
    if profiling_flag == 1
      disp('hemodynamic_state_to_measurement')
      toc
    end
    
    if profiling_flag == 1
      tic
    end
    % computeMeasurementMeanAndCovariance_Pz is just a redundant copy of unscentedTransform
    [x_bar_hh, P_bar_hh] = unscented_transform(Wm,Wc, Sigmas_h, R);
    
    % Computes the covariance between the projected states and the measurments
    P_bar_xh = cross_covariance(Wm, x_bar_xx, Sigmas_x, x_bar_hh, Sigmas_h);
    
    % This is where the magic happens ;)
    K = compute_kalman_gain(P_bar_xh, P_bar_hh);
    
    % difference between the actual measurement and where the projected
    % state thinks what the measurement should have been
    y = residual_y(zm(i), x_bar_hh);
    
    % Uses the Kalman gain to adjust the system's covariance matrix
    %   This tells you the accuracy of the states
    P = coveriance_update(K,P_bar_xx,P_bar_hh);
    
    % Uses the Kalman gain to adjust the system's state estimate
    %   This will be someplace between the uncented transform's project
    %   state and the measurement state
    x = state_update(x_bar_xx,K,y); % wrong Sigmas where going in here!!!
    x(3) = 1.0;
    
    sol = [sol, x];
    
    % update the last update timestamp
    t_last_update = t_prior; % this needs to be decoupled from sim
    if profiling_flag == 1
      disp('last calcs')
      toc
    end
end



% some plotting
state_index = find(strcmp(x_types, 'state') ==1);              % gets index of states
parameter_index = find(strcmp(x_types, 'parameter') ==1);      % gets index of parameters
constraint_index = find(strcmp(x_types, 'constraint') ==1);    % gets index of constraints


hold off
subplot(1,3,1)
plot(tm, zm(state_index,:) ,'o','LineWidth',1)
hold on
plot(tm, sol(state_index,:) ,'-','LineWidth',2)

subplot(1,3,2)
plot(tm, sol(parameter_index,:) ,'-','LineWidth',2)

subplot(1,3,3)
[res_sys_state, res_sys_time] = true_hemodynamic_state(process_model_function_user,  x0_init,  t_sim_start, dt_sim, 10 * 0.8,  [sol(2,end), sol(3,end), sol(4,end), cycle_time, Q_dt], 0)
plot(res_sys_time, res_sys_state(:,1)','-','LineWidth',2)