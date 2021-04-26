clear all
close all

c = 0.000750061683; % unit conversion 1 dyn/cm^2=0.0007750061 mmHg
x0_init = 80/c;

meas_variance = 0.5;%*10000000;

Cylces = 20;

% Measurement Time
t_sim_start = 0.0;         % start time
t_sim_end = Cylces * 0.8;  % end time
dt_sim = 1/500;            % time between measurements

% UKF prediction Time
t_last_update = t_sim_start; % time at last update
t_prior = 0.0;               % time at prior, this is updated in loop
dtp = dt_sim*0.01;           % time per step in sigma predecition projection

% -------------------------------- Sim --------------------------------



Rc_true = 1600;  % Silly
Rd_true = 13000; %
C_true = 0.000025; %

cycle_time = 0.8;
Q_dt = 1e-5;

cycles_to_skip = 0;

% additional arguments for hemo_pressure (alias of pressure_function)
true_hemodynamic_state_argz = [Rc_true, Rd_true, C_true, cycle_time, Q_dt];

starz = 1;
casez = 1;

parfor I = starz:casez
  disp(I)
  
  %if C_true-0.000001*I == (0.000000 || -0.000001)
  %  continue
  %end
  
  true_hemodynamic_state_argz = [Rc_true, Rd_true, C_true, cycle_time, Q_dt];
  [true_pressure_state, true_time_stamp] = true_hemodynamic_state(@pressure_function,  x0_init,  t_sim_start, dt_sim, t_sim_end,  true_hemodynamic_state_argz, cycles_to_skip);
  hold on
  
  sols(I,:) = true_pressure_state(:,1)
  
end

hold off
close all

for I = starz:casez
  hold on
  plot(sols(I,:))
end

for I = starz:casez
  disp(I)
  sols(I,1) - sols(I,end)
end

