% Calculate distributions of the generation time, time from onset of
% symptoms to transmission (TOST) and serial interval for the independent
% transmission and symptoms model, using the parameter values for which the
% likelihood is maximised.

% This script requires the Chebfun package to run (freely available at
% https://www.chebfun.org/download/).

clear all; close all; clc;
splitting on

addpath('../../Data')
addpath('../../Results')

% Load incubation period distribution, and the distribution of minus the
% incubation period
load('../../Data/input_data.mat','f_inc_cheb','f_m_inc_cheb')

% Load best-fitting parameter values
load('../../Results/mcmc_posterior_indep.mat','params_point')
params_indep = params_point;

% Generation time
t_range = [-500,-50,-25,0,25,50,500];
f_gen_indep = @(t,params)gampdf(t,params_indep(1),params_indep(2));
f_gen_indep = chebfun(f_gen_indep,t_range);

% TOST
f_tost_indep = conv(f_gen_indep,f_m_inc_cheb,'same','old');

% Serial interval
f_serial_indep = conv(f_inc_cheb,f_tost_indep,'same','old');

% Save results
save('../../Results/gen_tost_serial_indep.mat','f_gen_indep','f_tost_indep','f_serial_indep')

rmpath('../../Data')
rmpath('../../Results')