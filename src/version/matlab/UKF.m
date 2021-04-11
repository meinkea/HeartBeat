% v 1.5.0

clear all
close all
%==========================================================================
%=============================SET Constants================================
%==========================================================================

addpath(genpath('./ukf/')); % import from sub directories (installs the toolkit)

addpath('./simulation/hemodynamic/');
%addpath('./user_design/linear_testbench/');
addpath('./user_design/hemodynamic/');

%==============================Plots=======================================

Test_plots = 1; %flag

%%
%============================== SIMULATION ==============================
% ------------------------- Simulation Settings -------------------------
% These are
%x0_init = 8.777e+04; % True intial start
c = 0.000750061683; % unit conversion 1 dyn/cm^2=0.0007750061 mmHg
x0_init = 80/c;

meas_variance = 0.5;%*10000000;

Cylces = 3;

% Measurement Time
t_sim_start = 0.0;         % start time
t_sim_end = Cylces * 0.8;  % end time
dt_sim = 1/500;            % time between measurements

% UKF prediction Time
t_last_update = t_sim_start; % time at last update
t_prior = 0.0;               % time at prior, this is updated in loop
dtp = 0;           % time per step in sigma predecition projection

% -------------------------------- Sim --------------------------------



Rc_true = 1600;  % Silly
Rd_true = 13000; %
C_true = 0.000025; %

cycle_time = 0.8;
Q_dt = 1e-4;

cycles_to_skip = 4;

% additional arguments for hemo_pressure (alias of pressure_function)
true_hemodynamic_state_argz = [Rc_true, Rd_true, C_true, cycle_time, Q_dt];

% Generate true states
[true_pressure_state, true_time_stamp] = true_hemodynamic_state(@pressure_function,  x0_init,  t_sim_start, dt_sim, t_sim_end,  true_hemodynamic_state_argz, cycles_to_skip);

% Generate true measurements.  Returns also the new time series that accounts
[tm_true_stamp, true_measurement] = true_hemodynamic_measurement(@hemodynamic_state_to_measurement, true_pressure_state, true_time_stamp);

zm = noisy_measurements(true_measurement, meas_variance);

tm = tm_true_stamp;


%%
% ======================= SYSTEM DESIGN from USER =======================
Pref  = 80/c;      %mmHg
Rcref = 1600;
Rdref = 13000;
Cref  = 0.000025;


% -- State Vector
x_user = { 80/c  'state'    ;...
           Rcref 'parameter';...
           Rdref 'parameter';...
           Cref  'parameter'};

% State Estimator variance
P_user = [[ 8.0  0.0  0.0  0.0];...
          [ 0.0  8.0  0.0  0.0];...
          [ 0.0  0.0  8.0  0.0];...
          [ 0.0  0.0  0.0  8.0]];
  
% Process Model Variance and Covariance Matrix
process_model_function_user = @pressure_function;

dt_user = dtp;

a_ = 0.0001;
f_ = 0.01;
Q_user = ([[  a_,   0.0,   0.0,   0.0];...
           [ 0.0, f_*a_,   0.0,   0.0];...
           [ 0.0,   0.0, f_*a_,   0.0];...
           [ 0.0,   0.0,   0.0,  f_*a_]]);

% State Function Model for State Estimator
state_to_measurment_function = @hemodynamic_state_to_measurement;

% State Function Model parameters (vargz for state_transition())
state_transition_Vargz = [cycle_time Q_dt];


R_user = [0.005]; % For X meas only


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
    tic
    % --- PREDICT ---
    % Get time for prior ( which is the time of the next measurement )
    t_prior = tm(i+1); % this needs to be decoupled from sim
    
    
    sigmaPoints = compute_sigma_points(n, lambda, x, P); % Checked
    disp('sigmaPoints')
    toc

    tic
    % Project sigma points using ODE solver
    Sigmas_f = state_transition(process_model_function, state_transition_Vargz, sigmaPoints, x_types, t_last_update, dtp, t_prior); % Checked
    disp('state_transition')
    toc
    
    tic
    % _bar denotes that it came from unscentedTransform
    [ x_bar_xx, P_bar_xx ] =  unscented_transform(Wm,Wc, Sigmas_f, Q); % Checked
    %disp('unscented_transform')
    %toc
    
    tic
    % Recompute NEW sigmas at the new predicted state
    Sigmas_x = compute_sigma_points(n,lambda,x_bar_xx',P_bar_xx);
    %disp('compute_sigma_points')
    %toc
    % --- UPDATE ---
    
    tic
    % _h denotes measurement version of prediction
    %Sigmas_h = state_to_measurment(hemodynamic_state_to_measurement, Sigmas_x, 0);
    Sigmas_h = hemodynamic_state_to_measurement(Sigmas_x,  0.0);
    Sigmas_h = Sigmas_h';
    disp('hemodynamic_state_to_measurement')
    toc
    
    tic
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
    
    
    sol = [sol, x];
    
    % update the last update timestamp
    t_last_update = t_prior; % this needs to be decoupled from sim
    disp('last calcs')
    toc
end

toc

% some plotting
hold off
subplot(2,3,1)
plot(tm, zm(1,:) ,'o','LineWidth',1)
hold on
plot(tm, sol(1,:) ,'-','LineWidth',2)

subplot(2,3,2)
plot(tm, sol(2,:) ,'-','LineWidth',2)

subplot(2,3,3)
plot(tm, sol(3,:) ,'-','LineWidth',2)

subplot(2,3,4)
plot(tm, sol(4,:) ,'-','LineWidth',2)

[res_sys_time, res_sys_state] = true_hemodynamic_state(@pressure_function,  x0_init,  t_sim_start, dt_sim, 10 * 0.8,  [sol(2,end), sol(3,end), sol(4,end), cycle_time, Q_dt], 0)
subplot(2,3,5)
plot(res_sys_time, res_sys_state ,'-','LineWidth',2)