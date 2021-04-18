% v 1.5.0

clear all
close all
%==========================================================================
%========================= UKF Import Toolkit =============================
%==========================================================================

addpath(genpath('./ukf/')); % import from sub directories (installs the toolkit)

% -- Toolkit Flags

periodic_plotting = 1;
profiling_flag = 0;
dynamic_uncertainty = 0;

%%
%==========================================================================
%========================= User  System Setup =============================
%==========================================================================

% ===================== Import User System Directory =====================

addpath('./simulation/linear_with_parameter/');

% ======================== Measurement Simulation ========================

% -- Measurement Time
t_sim_start =  0.00;  % start time
t_sim_end   = 14.00;  % end time
dt_sim      =  0.01;  % time between measurements

tm = t_sim_start:dt_sim:t_sim_end;

% -- True Inital State
x_init = [1 2];

true_meas(1,:) = (3*sin(2*tm)+3*cos(2*tm)-1)/2;
true_meas(2,:) =  3*cos(2*tm)-1;

% -- Noisy Measurements
zm = noisy_measurements(true_meas, 0.005);


% ============================ USER UKF DESIGN ============================

% STATE DESIGN

% -- State Vector
x_user = {  zm(1) 'state';...
            zm(2) 'state';...
           -300.0 'parameter'};
        
% -- State Estimator Variance and Covariance Matrix
P_user = [[ 10000.0      0.0      0.0 ];...
          [     0.0  10000.0      0.0 ];...
          [     0.0      0.0  20000.0 ]];...


% PROCESS DESIGN

% -- Process Model
process_model_function__user = @linear_with_parameter__process_model;


% -- Process Model Forcast Time
t_last_update = tm(1); % time at last update
t_prior = 0.0;         % time at prior, this is updated in loop
dtp = 0;               % time per step in sigma predecition projection

% -- Process Model Variance and Covariance Matrix
Q_user = [[ 0.0001   0.0      0.0  ];...
          [ 0.0      0.0001   0.0  ];...
          [ 0.0      0.0      0.001]];

% -- State Function Model parameters (vargz for state_transition())
state_transition__vargz = [0];


% MEASUREMENT DESIGN

% -- Measurement function
get_measurement__user = @linear_with_parameter__get_measurement;

get_measurement__vargz = {tm zm x_user{3,1}};

% -- Measurement Variance and Covariance Matrix
R_user = [[0.1  0.00 0.00];...
          [0.00 0.1  0.00];...
          [0.00 0.00 0.01]];

% -- State to Measurement function
meas_to_state__user = @linear_with_parameter__meas_to_state;
meas_to_state__vargz = [0];


% UKF PARAMETER DESIGN

% -- Sigma Points tuning parameters --
alpha =  0.1;
beta  =  2.0;
kappa =  0.0;




%%


%==========================================================================
%================== >>> = UKF PARAMETER ESTIMATION = <<< ==================
%==========================================================================


% ======================== Initialize Matrices ===========================
x = [x_user{:,1}]';
x_types = {x_user{:,2}};
P = P_user;
Q = Q_user;
R = R_user;

process_model_function = process_model_function__user;


% number of state vectors
n = length( x );

% number of measured observations
nmo = length(R);

[sigmaPoints,weights,lambda,P_bar_xx,P_bar_hh] = initialize(n,alpha,kappa,nmo); % this should be broken up
sol = x;



% These are the weigths used in the uncented transform
% They are static and never change
[Wc,Wm] =  compute_weights(n,lambda,alpha,beta);

param_mes = [];
sol_K = [];


tic_step__periodic_plotting = 0;
tic_step__dynamic_uncertainty = 0;


disp('  ');
disp('    $$$$$$      $$$$$ ');
disp('  $$$$$$$$$$  $$$$$$$$$ ');
disp(' $$$$$$$$$$$$$$$$$$$$$$$$ ');
disp('$$$    /\          /\  $$$ ');
disp('$$$_  /  \__/\__  /  \_$$$ ');
disp('$$$ \/          \/     $$$ ');
disp(' $$$$$$$$$$$$$$$$$$$$$$$$ ');
disp('  $$$$$$$$$$$$$$$$$$$$$$ ');
disp('    $$$$$$$$$$$$$$$$$$ ');
disp('      $$$$$$$$$$$$$$$ ');
disp('        $$$$$$$$$$$$ ');
disp('          $$$$$$$$ ');
disp('           $$$$$$ ');
disp('            $$$$ ');
disp('             $$ ');
disp('                ');
disp('            $ ');
disp('           $$$ ');
disp('          $$$$$ ');
disp('         $$$$$$$ ');
disp('         $$$$$$$ ');
disp('          $$$$$ ');
disp('                 $ ');
disp('                $$$ ');
disp('               $$$$$ ');
disp('              $$$$$$$ ');
disp('              $$$$$$$ ');
disp('               $$$$$ ');
disp('  ');



% ============================= Main Loop ================================
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
    % The state_transition function can take and state equation function making the code general
    %
    % process_model_function - handle for the state equation function
    %                          must have the form: myStateEquFunction(t, y, argvz)
    % sigmaPoints - this it the intial sigma points
    Sigmas_f = state_transition(process_model_function, state_transition__vargz, sigmaPoints, x_types, t_last_update, dtp, t_prior); % Checked
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
    Sigmas_h = meas_to_state__user(Sigmas_x, meas_to_state__vargz);
    
    
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
    
    sol_K = [sol_K, K];
    
    % Get measurement for user defined measurement function
    meas = get_measurement__user(t_prior, x_bar_hh, sol, get_measurement__vargz);
    
    % difference between the actual measurement and where the projected
    % state thinks what the measurement should have been
    y = residual_y(meas, x_bar_hh);
    
    % Uses the Kalman gain to adjust the system's covariance matrix
    %   This tells you the accuracy of the states
    P = coveriance_update(K,P_bar_xx,P_bar_hh);
    
    % Uses the Kalman gain to adjust the system's state estimate
    %   This will be someplace between the uncented transform's project
    %   state and the measurement state
    x = state_update(x_bar_xx,K,y); % wrong Sigmas where going in here!!!
    
    
    sol = [sol, x];
    
    % update the last update timestamp
    t_last_update = t_prior; % this needs to be decoupled from sim
    if profiling_flag == 1
      disp('last calcs')
      toc
    end
    
    param_mes = [param_mes, meas(3)];
    
    
    % This is additional user defined functionality
    
    % periodic plotting
    if (periodic_plotting == 1) && (tic_step__periodic_plotting < tm(i))
      
      % some plotting
      state_index = find(strcmp(x_types, 'state') ==1);              % gets index of states
      parameter_index = find(strcmp(x_types, 'parameter') ==1);      % gets index of parameters
      constraint_index = find(strcmp(x_types, 'constraint') ==1);    % gets index of constraints
      
      hold off
      subplot(1,2,1)
      plot(tm(1:length(zm(state_index,:))), zm(state_index,:) ,'o','LineWidth',1)
      hold on
      plot(tm(1:length(sol(state_index,:))), sol(state_index,:) ,'-','LineWidth',2)
      
      subplot(1,2,2)
      plot(tm(1:length(param_mes([1],:))), param_mes([1],:) ,'o','LineWidth',1)
      
      plot(tm(1:length(sol(parameter_index,:))), sol(parameter_index,:) ,'-','LineWidth',2)
      
      tic_step__periodic_plotting = tic_step__periodic_plotting + 1;
      
      pause(1)
    end
    
    if (dynamic_uncertainty == 1)
      if (tic_step__dynamic_uncertainty < tm(i))
        tic_step__dynamic_uncertainty = tic_step__dynamic_uncertainty + 1;
      end
      
      if (tic_step__dynamic_uncertainty == 3)
        Q = [[ 0.0001    0.0       0.0       ];...
             [ 0.0       0.0001    0.0       ];...
             [ 0.0       0.0       0.000001  ]];
      elseif (tic_step__dynamic_uncertainty == 6)
        Q = [[ 0.0001    0.0       0.0        ];...
             [ 0.0       0.0001    0.0        ];...
             [ 0.0       0.0       0.0000001  ]];
      elseif (tic_step__dynamic_uncertainty == 9)
        Q = [[ 0.00001   0.0       0.0        ];...
             [ 0.0       0.00001   0.0        ];...
             [ 0.0       0.0       0.0000001  ]];
      elseif (tic_step__dynamic_uncertainty == 9)
        Q = [[ 0.000001  0.0       0.0        ];...
             [ 0.0       0.000001  0.0        ];...
             [ 0.0       0.0       0.0000001  ]];
      elseif (tic_step__dynamic_uncertainty == 12)
        Q = [[ 0.000001  0.0       0.0        ];...
             [ 0.0       0.000001  0.0        ];...
             [ 0.0       0.0       0.00000001 ]];
      end
    end
end


% some plotting
state_index = find(strcmp(x_types, 'state') ==1);              % gets index of states
parameter_index = find(strcmp(x_types, 'parameter') ==1);      % gets index of parameters
constraint_index = find(strcmp(x_types, 'constraint') ==1);    % gets index of constraints


hold off
subplot(1,2,1)
plot(tm(1:length(zm(state_index,:))), zm(state_index,:) ,'o','LineWidth',1)
hold on
plot(tm(1:length(sol(state_index,:))), sol(state_index,:) ,'-','LineWidth',2)

subplot(1,2,2)
plot(tm(1:length(sol(parameter_index,:))), sol(parameter_index,:) ,'-','LineWidth',2)

%subplot(1,3,3)
%[res_sys_state, res_sys_time] = true_hemodynamic_state(process_model_function_user,  x0_init,  t_sim_start, dt_sim, t_sim_end,  , 0)
%plot(res_sys_time, res_sys_state(:,1)','-','LineWidth',2)