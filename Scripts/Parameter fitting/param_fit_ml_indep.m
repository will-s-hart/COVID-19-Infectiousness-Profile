clear all; close all; clc;
rng(1)

addpath('../../Data')
addpath('../../Functions/General')
addpath('../../Functions/Indep')

% Load data
load('../../Data/input_data.mat','data_struct_observed','f_inc')

t_s1_mat = cell2mat((data_struct_observed.t_s1_data)')';
t_s2_mat = cell2mat((data_struct_observed.t_s2_data)')';
t_s1 = t_s1_mat(:,end/2);
t_s2 = t_s2_mat(:,end/2);

data_struct = rmfield(data_struct_observed,{'t_s1_data','t_s2_data'});
data_struct.t_s1_data = t_s1; data_struct.t_s2_data = t_s2;

f_gen_form = @(t_gen,theta)gampdf(t_gen,exp(2*theta(1))/exp(theta(2)),exp(theta(2))/exp(theta(1)));
ll_form = @(theta)sum(log_likelihood_indiv_indep(f_inc,@(t_gen)f_gen_form(t_gen,theta),data_struct));

theta_init = [1.7,1.7];

options = optimoptions(@fminunc,'OptimalityTolerance',1e-10,'StepTolerance',1e-10);
theta_best = fminunc(@(theta)-ll_form(theta),theta_init,options);

ll_best = ll_form(theta_best);

% Save results
save('../../Results/param_fit_ml_indep.mat','ll_best')

rmpath('../../Data')
rmpath('../../Functions/General')
rmpath('../../Functions/Indep')