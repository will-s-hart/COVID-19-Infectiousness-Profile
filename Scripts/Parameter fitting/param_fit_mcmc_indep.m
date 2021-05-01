% Fit the parameters of the independent transmission and symptoms model
% using data augmentation MCMC

% The vector of unknown model parameters, theta  = [log(m_gen),log(V_gen)],
% where m_gen and V_gen are the mean and variance of the generation time
% distribution, respectively, and the symptom onset times for each source
% and recipient, are updated in alternating steps of the model fitting
% procedure.

clear all; close all; clc;
rng(1)

addpath('../../Data')
addpath('../../Functions/General')
addpath('../../Functions/Indep')

% Load data
load('../../Data/input_data.mat','data_struct_observed','f_inc')

% Number of steps in chain
no_steps = 2500000;

% Function handele giving the generation time density at time since
% infection t_gen, conditional on incubation period of the source t_inc1,
% for parameters theta
f_gen_form = @(t_gen,theta)gampdf(t_gen,exp(2*theta(1))/exp(theta(2)),exp(theta(2))/exp(theta(1)));

% Function handle giving the log-likelihood for parameters theta and
% augmented data data_struct_augmented
ll_indiv_form = @(theta,data_struct_augmented)log_likelihood_indiv_indep(f_inc,@(t_gen)f_gen_form(t_gen,theta),data_struct_augmented);

% Standard deviations of proposal distributions for model parameters
sd_prop_lmean = 0.1;
sd_prop_lvar = 5*sd_prop_lmean;
theta_prop_sd = [sd_prop_lmean,sd_prop_lvar];

% Initial values of model parameters
theta_init = [1.75,1.5];

% Run MCMC
tic
[theta_mat,ll_vec,acceptance_rate] = fit_params_mcmc(ll_indiv_form,data_struct_observed,no_steps,theta_prop_sd,theta_init);
toc

% Save results
save('../../Results/param_fit_mcmc_indep.mat','theta_mat','ll_vec','acceptance_rate')

rmpath('../../Data')
rmpath('../../Functions/General')
rmpath('../../Functions/Indep')