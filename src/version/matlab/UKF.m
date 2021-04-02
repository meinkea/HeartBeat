% v 1.1.0

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

Test_plots = 0; %flag

%%
%============================== SIMULATION ==============================
% ------------------------- Simulation Settings -------------------------
% These are
%x0_init = 8.777e+04; % True intial start
c = 0.000750061683; % unit conversion 1 dyn/cm^2=0.0007750061 mmHg
x0_init = 80/c;

meas_variance = 0.5;

% Measurement Time
t_sim_start = 0.0;           % start time
t_sim_end = 3 * 0.8;           % end time
dt_sim = 1/500;                % time between measurements

% UKF prediction Time
t_last_update = t_sim_start; % time at last update
t_prior = 0.0;            % time at prior, this is updated in loop
dtp = dt_sim*0.01;                % time per step in sigma predecition projection

% -------------------------------- Sim --------------------------------



Rc_true = 1600;  % Silly
Rd_true = 13000; %
C_true = 2.5e-5; %

cycle_time = 0.8;
Q_dt = 1e-4;

cycles_to_skip = 4;

true_hemodynamic_state_argz = [Rc_true, Rd_true, C_true, cycle_time, Q_dt]; % additional arguments for hemo_pressure (alias of pressure_function)
[true_pressure_state, true_time_stamp] = true_hemodynamic_state(@pressure_function,  x0_init,  t_sim_start, dt_sim, t_sim_end,  true_hemodynamic_state_argz, cycles_to_skip);


[tm_true_stamp, true_measurement] = true_hemodynamic_measurement(@hemodynamic_state_to_measurement, true_pressure_state, true_time_stamp);
6
zm = noisy_measurements(true_measurement, meas_variance);

tm = tm_true_stamp;


%%
% ======================= SYSTEM DESIGN from USER =======================
Pref  = 50;      %mmHg
Cref  = (2e-5);    
Rcref = 1500;
Rdref = 15000;


% -- State Vector
x = [[              80/c];...
     [   log2(500/Rcref)];...
     [log2(0.00065/Cref)];...
     [ log2(25000/Rdref)]];

% State Estimator variance
P = [[ 8.0  0.0  0.0  0.0];...
     [ 0.0  8.0  0.0  0.0];...
     [ 0.0  0.0  8.0  0.0];...
     [ 0.0  0.0  0.0  8.0]];
   
% Process Model Variance and Covariance Matrix
process_model_function = @pressure_function;

dt = dtp;
Q = ([[(dt^7)/100, (dt^6)/60, (dt^5)/20, (dt^4)/8];...
      [(dt^6)/60,  (dt^5)/20, (dt^4)/8,  (dt^3)/6];...
      [(dt^5)/20,  (dt^4)*8,  (dt^3)/3,  (dt^2)/2];...
      [(dt^4)/8,   (dt^3)/6,  (dt^2)/2,  (dt^1)]]);

% State Function Model for State Estimator
state_to_measurment_function = @hemodynamic_state_to_measurement;

% State Function Model parameters (vargz for state_transition())
state_transition_Vargz = [0.0 0.0 0.0 cycle_time Q_dt];


R = [0.005]; % For X meas only

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
    plot(x_True(1),x_True(3),'o')
    hold on
    %plot(y_True, x_True, 'k', 'LineWidth', 3)
    plot(zm(:,1),zm(:,2), 'k', 'LineWidth', 3)
    grid on
    xlabel('Time (s)')
    ylabel('Pressure (log world)')
end


%==========================Initialize Matrices=============================

% number of state vectors
n = length(x);

% number of measured observations
nmo = length(R);

[sigmaPoints,weights,lambda,P_bar_xx,P_bar_hh] = initialize(n,alpha,kappa,nmo); % this should be broken up
sol =x;



% These are the weigths used in the uncented transform
% They are static and never change
[Wc,Wm] =  compute_weights(n,lambda,alpha,beta);

for i = 1:(length(tm)-1)
    % --- PREDICT ---
    % Get time for prior ( which is the time of the next measurement )
    t_prior = tm(i+1); % this needs to be decoupled from sim
    
    sigmaPoints = compute_sigma_points(n, lambda, x, P); % Checked
    
    
    % Project sigma points using ODE solver
    Sigmas_f = state_transition(process_model_function, state_transition_Vargz, sigmaPoints, t_last_update, dtp, t_prior); % Checked
    
    % _bar denotes that it came from unscentedTransform
    [ x_bar_xx, P_bar_xx ] =  unscented_transform(Wm,Wc, Sigmas_f, Q); % Checked
    
    % Recompute NEW sigmas at the new predicted state
    Sigmas_x = compute_sigma_points(n,lambda,x_bar_xx',P_bar_xx);
    
    % --- UPDATE ---
    
    % _h denotes measurement version of prediction
    %Sigmas_h = state_to_measurment(hemodynamic_state_to_measurement, Sigmas_x, 0);
    Sigmas_h = hemodynamic_state_to_measurement(Sigmas_x,  0.0);
    Sigmas_h = Sigmas_h';
    
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
    
    
    sol = [sol,x];
    
    % update the last update timestamp
    t_last_update = t_prior; % this needs to be decoupled from sim
end


% some plotting
hold off
plot(zm_sol(:,1), zm_sol(:,3) ,'o','LineWidth',3)
hold on
plot(sol(1,:)', sol(3,:) ,'-','LineWidth',4)



