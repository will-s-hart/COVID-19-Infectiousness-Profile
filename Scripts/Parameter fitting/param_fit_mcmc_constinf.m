% Fit the parameters of the constant infectiousness model using data
% augmentation MCMC.

% The vector of unknown model parameters, theta  = [log(k_E),log(mu)], and
% the symptom onset times for each source and recipient, are updated in
% alternating steps of the model fitting procedure.

clear all; close all; clc;
rng(1)

addpath('../../Data')
addpath('../../Functions/General')
addpath('../../Functions/Mechanistic')

% Load data
load('../../Data/input_data.mat','data_struct_observed','k_inc','gamma','f_inc','k_I')

% Number of steps in chain
no_steps = 2500000;

% Function handle giving infectiousness at time since onset
% t_tost, conditional on incubation period of the source t_inc1, for
% parameters theta
b_cond_form = @(t_tost,t_inc1,theta)b_cond_form_mechanistic(t_tost,t_inc1,get_params_constinf(theta,k_inc,gamma,k_I));

% Function handle giving the log-likelihood for parameters theta and
% augmented data data_struct_augmented
ll_indiv_form = @(theta,data_struct_augmented)log_likelihood_indiv_mechanistic(f_inc,@(t_tost,t_inc1)b_cond_form(t_tost,t_inc1,theta),data_struct_augmented);

% Standard deviations of proposal distributions for model parameters
sd_prop_lkE = 0.225;
sd_prop_lmu = 1*sd_prop_lkE;
theta_prop_sd = [sd_prop_lkE,sd_prop_lmu];

% Initial values of model parameters
theta_init = [1,-0.75];

% Run MCMC
tic
[theta_mat,ll_vec,acceptance_rate] = fit_params_mcmc(ll_indiv_form,data_struct_observed,no_steps,theta_prop_sd,theta_init);
toc

% Save results
save('../../Results/param_fit_mcmc_constinf.mat','theta_mat','ll_vec','acceptance_rate')

rmpath('../../Data')
rmpath('../../Functions/General')
rmpath('../../Functions/Mechanistic')