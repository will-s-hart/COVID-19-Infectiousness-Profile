% Import and format raw transmission pair data, and set up incubation
% period distribution and assumed model parameter values.

% This script requires the Chebfun package to run (freely available at
% https://www.chebfun.org/download/).

clear all; close all; clc;
splitting on

% Import transmission pair data
T = readtable('transmission_pair_data.xlsx');

% Middle of day of symptom onset for sources and recipients
t_s1_mid = T.t_s1;
t_s2_mid = T.t_s2;

% Empirical serial intervals
serial_data = t_s2_mid-t_s1_mid;

% Start and end of day of symptom onset for sources and recipients
t_s1L = t_s1_mid-0.5;
t_s1R = t_s1_mid+0.5;
t_s2L = t_s2_mid-0.5;
t_s2R = t_s2_mid+0.5;

% Left and right bounds for infection times
t_i1L = T.t_i1L-0.5; t_i1L(isnan(t_i1L))=-inf;
t_i1R = T.t_i1R+0.5; t_i1R(isnan(t_i1R))=t_s1R(isnan(t_i1R));
t_i2L = T.t_i2L-0.5; t_i2L(isnan(t_i2L))=t_i1L(isnan(t_i2L));
t_i2R = T.t_i2R+0.5; t_i2R(isnan(t_i2R))=t_s2R(isnan(t_i2R));

% Remove impossible cases where recipient infected before source
t_i1R(t_i2R<t_i1R) = t_i2R(t_i2R<t_i1R);

% Remove impossible cases where infection occurs after symptom onset
t_i1R(t_s1R<t_i1R) = t_s1R(t_s1R<t_i1R);
t_i2R(t_s2R<t_i2R) = t_s2R(t_s2R<t_i2R);

% Set maximum allowed incubation period
t_inc_max = 30;

% Remove possible incubation periods longer than the allowed maximum
t_i1L = max(t_i1L,t_s1L-t_inc_max);
t_i2L = max(t_i2L,t_s2L-t_inc_max);

if min(t_i1R-t_i1L)<0
        error('No plausible range of values for infection time of at least one source')
end
if min(t_i2R-t_i2L)<0
        error('No plausible range of values for infection time of at least one recipient')
end

% Time step used for numerical integration
dt = 0.125;

% Create cell arrays containing vectors of possible (discretised) infection
% and onset times for each host
t_s1_data = cell(length(t_s1_mid),1);
t_s2_data = cell(length(t_s1_mid),1);
t_i1_data = cell(length(t_s1_mid),1);
t_i2_data = cell(length(t_s1_mid),1);

for i = 1:length(t_s1_data)
    t_s1_data{i} = ((t_s1L(i)+dt):dt:t_s1R(i))';
    t_s2_data{i} = ((t_s2L(i)+dt):dt:t_s2R(i))';
    t_i1_data{i} = ((t_i1L(i)+(dt/2)):dt:(t_i1R(i)-(dt/2)))';
    t_i2_data{i} = ((t_i2L(i)+(dt/2)):dt:(t_i2R(i)-(dt/2)))';
end

% Create structure array containing observed data
data_struct_observed.t_s1_data = t_s1_data;
data_struct_observed.t_s2_data = t_s2_data;
data_struct_observed.t_i1_data = t_i1_data;
data_struct_observed.t_i2_data = t_i2_data;
data_struct_observed.dt = dt;

% Store minimum and maximum possible values for epidemiological time
% intervals
data_struct_observed.t_inc_min = min([t_s1L-t_i1R;t_s2L-t_i2R]); %minimum may be less than 0 for ease of interpolation when evaluating likelihood
data_struct_observed.t_inc_max = max([t_s1R-t_i1L;t_s2R-t_i2L]); %not quite equal to max allowed inc period set above due to uncertainty in onset times
data_struct_observed.t_tost_min = min(t_i2L-t_s1R);
data_struct_observed.t_tost_max = max(t_i2R-t_s1L);
data_struct_observed.t_gen_min = min(t_i2L-t_i1R); %minimum may be less than 0 for ease of interpolation when evaluating likelihood
data_struct_observed.t_gen_max = max(t_i2R-t_i1L);

% Assumed parameters of the incubation period distribution
inc_shape_lauer = 5.807;
inc_scale_lauer = 0.948;
inc_shape = inc_shape_lauer;
inc_scale = inc_scale_lauer;

% Function handle for the incubation period distribution
f_inc = @(t_inc)gampdf(t_inc,inc_shape,inc_scale);

% Chebfun for the incubation period distributions
t_range = [-500,-50,-25,0,25,50,500];
f_inc_cheb = chebfun(f_inc,t_range);

% Chebfun for the distribution of minus the incubation period
f_m_inc_cheb = chebfun(@(t)f_inc(-t),t_range);

% Assumed parameter values for the mechanistic model
k_inc = inc_shape;
gamma = 1/(k_inc*inc_scale);
k_I = 1;

% Save data to .mat file
save('input_data.mat','data_struct_observed','serial_data','k_inc','gamma','f_inc','f_inc_cheb','f_m_inc_cheb','inc_shape','inc_scale','k_I')