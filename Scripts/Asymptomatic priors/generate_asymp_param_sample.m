% Obtain samples from the assumed prior distributions for the proportion
% and relative infectiousness of asymptomatic hosts.

clear all; close all; clc;

addpath('../../Results')

% Get required sample length
load('../../Results/mcmc_posterior_indep','prob_presymp_post')
sample_size = length(prob_presymp_post);

% Sample from assumed distribution for the proportion of asymptomatic cases
load('../../Results/asymp_proportion_params','a_best','b_best')
p_A_sample = betarnd(a_best,b_best,sample_size,1);

% Sample from assumed distribution for the relative infectiousness of
% asymptomatic hosts
load('../../Results/asymp_infectiousness_params','mu_best','sigma_best')
x_A_sample = lognrnd(mu_best,sigma_best,sample_size,1);

% Save samples
save('../../Results/asymp_param_samples','p_A_sample','x_A_sample');

rmpath('../../Results')