clear all; close all; clc;
rng(1)

addpath('../../Data')
addpath('../../Functions/General')
addpath('../../Functions/Ferretti')

% Load data
load('../../Data/input_data.mat','data_struct_observed','f_inc','inc_shape','inc_scale')

m_inc = inc_shape*inc_scale;

t_s1_mat = cell2mat((data_struct_observed.t_s1_data)')';
t_s2_mat = cell2mat((data_struct_observed.t_s2_data)')';
t_s1 = t_s1_mat(:,end/2);
t_s2 = t_s2_mat(:,end/2);

data_struct = rmfield(data_struct_observed,{'t_s1_data','t_s2_data'});
data_struct.t_s1_data = t_s1; data_struct.t_s2_data = t_s2;

b_cond_form = @(t_tost,t_inc1,theta)b_cond_form_ferretti(t_tost,t_inc1,get_params_ferretti(theta,m_inc));
ll_form = @(theta)sum(log_likelihood_indiv_ferretti(f_inc,@(t_tost,t_inc1)b_cond_form(t_tost,t_inc1,theta),data_struct));

theta_init = [-5,0.5,2.5];

options = optimoptions(@fminunc,'OptimalityTolerance',1e-10,'StepTolerance',1e-10);
theta_best = fminunc(@(theta)-ll_form(theta),theta_init,options);

ll_best = ll_form(theta_best);

% Save results
save('../../Results/param_fit_ml_ferretti.mat','ll_best')

rmpath('../../Data')
rmpath('../../Functions/General')
rmpath('../../Functions/Ferretti')