% Fit the parameters of the Ferretti model using data augmentation MCMC.

% The vector of unknown model parameters, theta  =
% [mu,log(sigma),log(alpha)], and the symptom onset times for each source
% and recipient, are updated in alternating steps of the model fitting
% procedure.

clear all; close all; clc;
rng(2)

addpath('../../Data')
addpath('../../Functions/General')
addpath('../../Functions/Ferretti')

% Load data
load('../../Data/input_data.mat','data_struct_observed','f_inc','inc_shape','inc_scale')

m_inc = inc_shape*inc_scale;

% Number of steps in chain
no_steps = 2500000;

% Function handle giving infectiousness at time since onset
% t_tost, conditional on incubation period of the source t_inc1, for
% parameters theta
b_cond_form = @(t_tost,t_inc1,theta)b_cond_form_ferretti(t_tost,t_inc1,get_params_ferretti(theta,m_inc));

% Function handle giving the log-likelihood for parameters theta and
% augmented data data_struct_augmented
ll_indiv_form = @(theta,data_struct_augmented)log_likelihood_indiv_ferretti(f_inc,@(t_tost,t_inc1)b_cond_form(t_tost,t_inc1,theta),data_struct_augmented);

% Standard deviations of proposal distributions for model parameters
sd_prop_mu = 0.5;
sd_prop_lsigma = 0.05*sd_prop_mu;
sd_prop_lalpha = 0.5*sd_prop_mu;
theta_prop_sd = [sd_prop_mu,sd_prop_lsigma,sd_prop_lalpha];

% Initial values of model parameters
theta_init = [-2.5,0.5,2];

% Run MCMC
tic
[theta_mat,ll_vec,acceptance_rate] = fit_params_mcmc(ll_indiv_form,data_struct_observed,no_steps,theta_prop_sd,theta_init);
toc

% Save results
save('../../Results/param_fit_mcmc_ferretti.mat','theta_mat','ll_vec','acceptance_rate')

rmpath('../../Data')
rmpath('../../Functions/General')
rmpath('../../Functions/Ferretti')