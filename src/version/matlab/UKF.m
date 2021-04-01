% v 1.0.2

clear all
close all
%==========================================================================
%=============================SET Constants================================
%==========================================================================



%==============================Plots=======================================

Test_plots = 1; %flag

%============================Initialize State==============================
% x = [log2(80/Pref)  log2((500*c)/Rcref)  log2((0.00065/c)/Cref) log2((25000*c)/Rdref)...
%      log2((500*c)/Rcref)+log2((25000*c)/Rdref)];

%x = [[ 30.0];...
%     [-10.0]];

%P = [[32.0 15.0];...
%     [15.0 40.0]];
   
   
x = [[ 30.0];...
     [-10.0];...
     [ 50.0];...
     [-10.0]];

P = [[32.0 15.0  0.0  0.0];...
     [15.0 40.0  0.0  0.0];...
     [ 0.0  0.0 32.0 15.0];...
     [ 0.0  0.0 15.0 40.0]];

n = length(x);

%number of measured Observations
nmo = 2;

%=========================Sigma Points tuning parameters===================

alpha = 0.3;
beta = 2.0;
kappa = -1; % prolly should be set to n-3, so -1

%=========================Uncertainty and noise parameters=================
% variance in noise from sensors
%R = [[0.3^2 0.0  ];...
%     [0.0   0.3^2]];
   
R = [[0.3^2 0.0  ];...
     [0.0   0.3^2]]; % For X meas only


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

%[x_True, zm, tm] = Noisy_Measurements(@Pressure_Funtion, [dtm], x0_init, tm_start, dtm, tm_end, meas_variance);

[x_True, zm_sol, tm] = Noisy_Measurements(@Pressure_Funtion, [dtm], x0_init, tm_start, dtm, tm_end, meas_variance);
%[y_True, zm, tm] = Noisy_Measurements(@Pressure_Funtion, [dtm], x0_init, tm_start, dtm, tm_end, meas_variance);

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


[Wc,Wm] =  computeWeights(n,lambda,alpha,beta); % this never changes, moved before loop

for i = 1:(length(tm)-1)
    % --- PREDICT ---
    % Get time for prior ( which is the time of the next measurement )
    
    t_prior = tm(i+1);
    %P
    sigmaPoints = computeSigmaPoints(n, lambda, x, P); % Checked
    % dt = 0.002
    
    % Project sigma points using ODE solver
    Sigmas_f = stateTransition(@Pressure_Funtion, [dtm], sigmaPoints, t_last_update, dtp, t_prior); % Checked
    
    % _bar denotes that it came from unscentedTransform
    [ x_bar_xx, P_bar_xx ] =  unscentedTransform(Wm,Wc, Sigmas_f, Q); % Checked
    [ x_bar_xx, P_bar_xx ] =  unscented_transform(Wm,Wc, Sigmas_f, Q); % Checked
    
    % Recompute NEW sigmas at the new predicted state
    Sigmas_x = computeSigmaPoints(n,lambda,x_bar_xx',P_bar_xx);
    
    % --- UPDATE ---
    
    % _h denotes measurement version of prediction
    Sigmas_h = computeMeasurementSigmas_Z(@h_func, Sigmas_x, 0); 
    
    % computeMeasurementMeanAndCovariance_Pz is just a redundant copy of unscentedTransform
    [x_bar_hh, P_bar_hh] = unscentedTransform(Wm,Wc, Sigmas_h, R);
    
    P_bar_xh = cross_covariance(Wm, x_bar_xx, Sigmas_x, x_bar_hh, Sigmas_h);
    
    K = compute_Kalman_gain(P_bar_xh, P_bar_hh);
    
    y = Residual_y(zm(i,:), x_bar_hh);
    
    P = CoverianceUpdate(K,P_bar_xx,P_bar_hh);
    
    x = state_x(x_bar_xx,K,y); % wrong Sigmas where going in here!!!
    
    
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

