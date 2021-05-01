% Import posterior parameter distributions for the independent transmission
% and symptoms model obtained using data augmentation MCMC, and calculate
% the posterior distribution of the proportion of presymptomatic
% transmissions.

clear all; close all; clc;

addpath('../../Data')
addpath('../../Results')
addpath('../../Functions/Indep')

% Load known model parameters
load('../../Data/input_data.mat','inc_shape','inc_scale')

% Load output of MCMC fitting procedure
load('../../Results/param_fit_mcmc_indep.mat','theta_mat','ll_vec','acceptance_rate')

% Calculate posterior distributions of individual model parameters,
% assuming uniform priors
gen_mean_post = exp(theta_mat(:,1));
gen_var_post = exp(theta_mat(:,2));
gen_sd_post = sqrt(gen_var_post);
gen_shape_post = (gen_mean_post.^2)./gen_var_post;
gen_scale_post = gen_var_post./gen_mean_post;

% Remove burn-in and thin posteriors
no_steps = length(gen_mean_post);
indices_kept = (500001):100:no_steps;

gen_mean_post = gen_mean_post(indices_kept);
gen_var_post = gen_var_post(indices_kept);
gen_shape_post = gen_shape_post(indices_kept);
gen_scale_post = gen_scale_post(indices_kept);

% Calculate posterior distribution for the proportion of presymptomatic
% transmissions
prob_presymp_post = get_presymp_trans_probs_indep(gen_shape_post,gen_scale_post,inc_shape,inc_scale);

% Median and 95% confidence bounds for the proportion of presymptomatic
% transmissions
quantile(prob_presymp_post,[0.025,0.5,0.975])

% Find posterior mean parameters
theta_point = mean(theta_mat(indices_kept,:));
gen_mean_point = exp(theta_point(1));
gen_var_point = exp(theta_point(2));
gen_sd_point = sqrt(gen_var_point);
gen_shape_point = (gen_mean_point.^2)./gen_var_point;
gen_scale_point = gen_var_point./gen_mean_point;
params_point = [gen_shape_point,gen_scale_point];
prob_presymp_point = get_presymp_trans_probs_indep(params_point(1),params_point(2),inc_shape,inc_scale);

% Save results
save('../../Results/mcmc_posterior_indep','params_point','prob_presymp_point','prob_presymp_post','gen_shape_post','gen_scale_post')

rmpath('../../Data')
rmpath('../../Results')
rmpath('../../Functions/Indep')