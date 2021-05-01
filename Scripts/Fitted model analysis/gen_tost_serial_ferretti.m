% Calculate distributions of the generation time, time from onset of
% symptoms to transmission (TOST) and serial interval for the Ferretti
% model, using the parameter values for which the likelihood is maximised.

% This script requires the Chebfun package to run (freely available at
% https://www.chebfun.org/download/).

clear all; close all; clc;
splitting on

addpath('../../Data')
addpath('../../Results')
addpath('../../Functions/Ferretti')

% Load incubation period distribution
load('../../Data/input_data.mat','f_inc','f_inc_cheb')

% Load best-fitting parameter values
load('../../Results/mcmc_posterior_ferretti.mat','params_point')
params_ferretti = params_point;

% Function handle for infectiousness at time since onset t_tost,
% conditional on incubation period t_inc1, under the best-fitting
% parameters
b_cond_ferretti = @(t_tost,t_inc1)b_cond_form_ferretti(t_tost,t_inc1,params_ferretti);

% Generation time
f_gen_ferretti = get_gen_dist_ferretti(f_inc,b_cond_ferretti);

% TOST
f_tost_ferretti1 = get_tost_dist_ferretti(f_inc,b_cond_ferretti);
f_tost_ferretti = chebfun(f_tost_ferretti1,[-50,-25,0,25,50],100); %calculate Chebyshev interpolant at fewer gridpoints to ensure efficient calculatio of serial interval distribution

% Serial interval
f_serial_ferretti = conv(f_inc_cheb,f_tost_ferretti,'same','old');

% Save results
save('../../Results/gen_tost_serial_ferretti.mat','f_gen_ferretti','f_tost_ferretti','f_serial_ferretti')

rmpath('../../Data')
rmpath('../../Results')
rmpath('../../Functions/Ferretti')