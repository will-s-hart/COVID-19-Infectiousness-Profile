clear all; close all; clc;
rng(1)

addpath('../../Data')
addpath('../../Functions/General')
addpath('../../Functions/Mechanistic')

% Load data
load('../../Data/input_data.mat','data_struct_observed','k_inc','gamma','f_inc','k_I')

t_s1_mat = cell2mat((data_struct_observed.t_s1_data)')';
t_s2_mat = cell2mat((data_struct_observed.t_s2_data)')';
t_s1 = t_s1_mat(:,end/2);
t_s2 = t_s2_mat(:,end/2);

data_struct = rmfield(data_struct_observed,{'t_s1_data','t_s2_data'});
data_struct.t_s1_data = t_s1; data_struct.t_s2_data = t_s2;

b_cond_form = @(t_tost,t_inc1,theta)b_cond_form_mechanistic(t_tost,t_inc1,get_params_constinf(theta,k_inc,gamma,k_I));
ll_form = @(theta)sum(log_likelihood_indiv_mechanistic(f_inc,@(t_tost,t_inc1)b_cond_form(t_tost,t_inc1,theta),data_struct));

theta_init = [1,-0.8];

options = optimoptions(@fminunc,'OptimalityTolerance',1e-10,'StepTolerance',1e-10);
theta_best = fminunc(@(theta)-ll_form(theta),theta_init,options);

ll_best = ll_form(theta_best);

% Save results
save('../../Results/param_fit_ml_constinf.mat','ll_best')

rmpath('../../Data')
rmpath('../../Functions/General')
rmpath('../../Functions/Mechanistic')