% Import posterior parameter distributions for the variable infectiousness
% model obtained using data augmentation MCMC, and calculate the posterior
% distribution of the proportion of presymptomatic transmissions.

clear all; close all; clc;

addpath('../../Data')
addpath('../../Results')
addpath('../../Functions/Mechanistic')

% Load known model parameters
load('../../Data/input_data.mat','k_inc','gamma','k_I')

% Load output of MCMC fitting procedure
load('../../Results/param_fit_mcmc_varinf.mat','theta_mat','ll_vec','acceptance_rate')

% Calculate posterior distributions of individual model parameters,
% assuming uniform priors
k_E_post = exp(theta_mat(:,1));
k_P_post = k_inc - k_E_post;
mu_post = exp(theta_mat(:,2));
alpha_post = exp(theta_mat(:,3));

% Calculate posterior distribution for the proportion of presymptomatic
% transmissions
prob_presymp_post = (alpha_post.*k_P_post/(k_inc*gamma))./((alpha_post.*k_P_post/(k_inc*gamma))+1./mu_post);

% Remove burn-in and thin posteriors
no_steps = length(k_E_post);
indices_kept = (500001):100:no_steps;

theta_post = theta_mat(indices_kept,:);
k_E_post = k_E_post(indices_kept);
k_P_post = k_P_post(indices_kept);
mu_post = mu_post(indices_kept);
alpha_post = alpha_post(indices_kept);
prob_presymp_post = prob_presymp_post(indices_kept);

% Median and 95% confidence bounds for the proportion of presymptomatic
% transmissions
quantile(prob_presymp_post,[0.025,0.5,0.975])

% Find posterior mean parameters
theta_point = mean(theta_mat(indices_kept,:));
k_E_point = exp(theta_point(1));
k_P_point = k_inc - k_E_point;
mu_point = exp(theta_point(2));
alpha_point = exp(theta_point(3));
% k_E_point = mean(k_E_post);
% k_P_point = k_inc - k_E_point;
% mu_point = mean(mu_post);
% alpha_point = mean(alpha_post);
% theta_point = log([k_E_point,mu_point,alpha_point]);

params_point = get_params_varinf(theta_point,k_inc,gamma,k_I);
prob_presymp_point = (alpha_point.*k_P_point/(k_inc*gamma))./((alpha_point.*k_P_point/(k_inc*gamma))+1./mu_point);

% Save results
save('../../Results/mcmc_posterior_varinf','theta_point','params_point','prob_presymp_point','prob_presymp_post','theta_post')

rmpath('../../Data')
rmpath('../../Results')
rmpath('../../Functions/Mechanistic')