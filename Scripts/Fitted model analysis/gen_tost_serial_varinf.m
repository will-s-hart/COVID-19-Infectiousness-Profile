% Calculate distributions of the generation time, time from onset of
% symptoms to transmission (TOST) and serial interval for the variable
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
load('../../Results/mcmc_posterior_varinf.mat','params_point')
params_varinf = params_point;

% Generation time
f_gen_varinf = get_gen_dist_mechanistic(params_varinf);

% TOST
t_range = [-500,-50,-25,0,25,50,500];
f_tost_varinf = chebfun(@(t)f_tost_form_mechanistic(t,params_varinf),t_range);

% Serial interval
f_serial_varinf = get_serial_dist_mechanistic(params_varinf);

% Save results
save('../../Results/gen_tost_serial_varinf.mat','f_gen_varinf','f_tost_varinf','f_serial_varinf')

rmpath('../../Data')
rmpath('../../Results')
rmpath('../../Functions/Mechanistic')