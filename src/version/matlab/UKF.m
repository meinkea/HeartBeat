% v 1.0.4

clear all
close all
%==========================================================================
%=============================SET Constants================================
%==========================================================================



%==============================Plots=======================================

Test_plots = 0; %flag

%%
%============================Initialize State==============================

% State Vector
x = [[ 30.0];...
     [-10.0];...
     [ 50.0];...
     [-10.0]];

% State Estimator variance
P = [[32.0 15.0  0.0  0.0];...
     [15.0 40.0  0.0  0.0];...
     [ 0.0  0.0 32.0 15.0];...
     [ 0.0  0.0 15.0 40.0]];

% Process Model Variance and Covariance Matrix
process_model_function = @xdot_func;

Q = [[0.02  0.01  0.00  0.00];...
     [0.01  0.02  0.00  0.00];...
     [0.00  0.00  0.02  0.01];...
     [0.00  0.00  0.01  0.02]];


% Measurment Variance
state_to_measurment_function = @h_func;

R = [[0.3^2 0.0  ];...
     [0.0   0.3^2]]; % For X meas only


n = length(x);

%number of measured Observations
nmo = 2;

% --- Sigma Points tuning parameters ---

alpha = 0.3;
beta = 2.0;
kappa = -1; % prolly should be set to n-3, so -1

  

%%
%===============Initial Guess and Parameters for DE solver=================
% Measurement Time

tm_start = 0.0; % start time
tm_end = 100.0; % end time
dtm = 1.0;      % time between measurements


%R = 0.005;
% Noise from process model
%Qp = 0;
%Q = [[0.02  0.01];...
%     [0.01  0.02]];
Q = [[0.02  0.01  0.00  0.00];...
     [0.01  0.02  0.00  0.00];...
     [0.00  0.00  0.02  0.01];...
     [0.00  0.00  0.01  0.02]];
%==============Making Measuments and comparing to True funtion=============

%setting seed
rng('default')

x0_init = [0.0 1.0, 0.0, 1.0]'; % True intial start
meas_variance = 0.5;

[x_True, zm_sol, tm] = noisy_measurements(process_model_function, [dtm], x0_init, tm_start, dtm, tm_end, meas_variance);

zm = zm_sol(:,[1 3]); % For position meas only

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

[sigmaPoints,weights,lambda,P_bar_xx,P_bar_hh] = initialize(n,alpha,kappa,nmo); % this should be broken up
sol =x;

% UKF prediction Time
t_last_update = tm(1); % time at last update
t_prior = tm(1);       % time at prior, this is updated in loop
dtp = 0.1;             % time per step in sigma predecition projection


[Wc,Wm] =  compute_weights(n,lambda,alpha,beta); % this never changes, moved before loop

for i = 1:(length(tm)-1)
    % --- PREDICT ---
    % Get time for prior ( which is the time of the next measurement )
    
    t_prior = tm(i+1);
    %P
    sigmaPoints = compute_sigma_points(n, lambda, x, P); % Checked
    % dt = 0.002
    
    % Project sigma points using ODE solver
    Sigmas_f = state_transition(process_model_function, [dtm], sigmaPoints, t_last_update, dtp, t_prior); % Checked
    
    % _bar denotes that it came from unscentedTransform
    [ x_bar_xx, P_bar_xx ] =  unscented_transform(Wm,Wc, Sigmas_f, Q); % Checked
    
    % Recompute NEW sigmas at the new predicted state
    Sigmas_x = compute_sigma_points(n,lambda,x_bar_xx',P_bar_xx);
    
    % --- UPDATE ---
    
    % _h denotes measurement version of prediction
    Sigmas_h = state_to_measurment(state_to_measurment_function, Sigmas_x, 0); 
    
    % computeMeasurementMeanAndCovariance_Pz is just a redundant copy of unscentedTransform
    [x_bar_hh, P_bar_hh] = unscented_transform(Wm,Wc, Sigmas_h, R);
    
    P_bar_xh = cross_covariance(Wm, x_bar_xx, Sigmas_x, x_bar_hh, Sigmas_h);
    
    K = compute_kalman_gain(P_bar_xh, P_bar_hh);
    
    y = residual_y(zm(i,:), x_bar_hh);
    
    P = coveriance_update(K,P_bar_xx,P_bar_hh);
    
    x = state_update(x_bar_xx,K,y); % wrong Sigmas where going in here!!!
    
    
    sol = [sol,x];
    
    % update the last update timestamp
    t_last_update = t_prior;
end

hold off
plot(zm_sol(:,1), zm_sol(:,3) ,'o','LineWidth',3)
hold on
plot(sol(1,:)', sol(3,:) ,'-','LineWidth',4)


%plot(sol(2:end), sol(4,:) ,'*','LineWidth',4)
%ylim([-15 125])

