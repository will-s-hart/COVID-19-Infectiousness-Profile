% Calculate distributions of the generation time, time from onset of
% symptoms to transmission (TOST) and serial interval for the constant
% infectiousness model, using the parameter values for which the likelihood
% is maximised.

% This script requires the Chebfun package to run (freely available at
% https://www.chebfun.org/download/).

clear all; close all; clc;
splitting on

addpath('../../Data')
addpath('../../Results')
addpath('../../Functions/Mechanistic')

% Load best-fitting parameter values
load('../../Results/mcmc_posterior_constinf.mat','params_point')
params_constinf = params_point;

% Generation time
f_gen_constinf = get_gen_dist_mechanistic(params_constinf);

% TOST
t_range = [-500,-50,-25,0,25,50,500];
f_tost_constinf = chebfun(@(t)f_tost_form_mechanistic(t,params_constinf),t_range);

% Serial interval
f_serial_constinf = get_serial_dist_mechanistic(params_constinf);

% Save results
save('../../Results/gen_tost_serial_constinf.mat','f_gen_constinf','f_tost_constinf','f_serial_constinf')

rmpath('../../Data')
rmpath('../../Results')
rmpath('../../Functions/Mechanistic')