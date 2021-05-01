% Calculate both the proportion of contacts identified for different
% contact elicitation windows, and the reduction in transmission through
% isolation for different times between either exposure or symptom onset
% and isolation.

% This script requires the Chebfun package to run (freely available at
% https://www.chebfun.org/download/).

clear all; close all; clc;
splitting on

addpath('../../Functions/Indep')
addpath('../../Functions/Mechanistic')
addpath('../../Results')

load('../../Results/gen_tost_serial_indep.mat','f_gen_indep','f_tost_indep','f_serial_indep')
load('../../Results/gen_tost_serial_constinf.mat','f_gen_constinf','f_tost_constinf','f_serial_constinf')
load('../../Results/gen_tost_serial_varinf.mat','f_gen_varinf','f_tost_varinf','f_serial_varinf')
load('../../Results/gen_tost_serial_ferretti.mat','f_gen_ferretti','f_tost_ferretti','f_serial_ferretti')

% Proportion of infectious presymptomatic contacts found when
% contacts of a symptomatic host are traced

F_tost_indep = cumsum(f_tost_indep);
F_tost_constinf = cumsum(f_tost_constinf);
F_tost_varinf = cumsum(f_tost_varinf);
F_tost_ferretti = cumsum(f_tost_ferretti);

t_traced_from = chebfun('t',[0,7]);
contacts_traced_indep = (F_tost_indep(0)-F_tost_indep(-t_traced_from))/(F_tost_indep(0));
contacts_traced_constinf = (F_tost_constinf(0)-F_tost_constinf(-t_traced_from))/(F_tost_constinf(0));
contacts_traced_varinf = (F_tost_varinf(0)-F_tost_varinf(-t_traced_from))/(F_tost_varinf(0));
contacts_traced_ferretti = (F_tost_ferretti(0)-F_tost_ferretti(-t_traced_from))/(F_tost_ferretti(0));

% Transmissions prevented when contacts isolated at different
% times since exposure

F_gen_indep = cumsum(f_gen_indep);
F_gen_constinf = cumsum(f_gen_constinf);
F_gen_varinf = cumsum(f_gen_varinf);
F_gen_ferretti = cumsum(f_gen_ferretti);

t_isolated_con = chebfun('t',[0,10]);
transm_prev_con_indep = 1-F_gen_indep(t_isolated_con);
transm_prev_con_constinf = 1-F_gen_constinf(t_isolated_con);
transm_prev_con_varinf = 1-F_gen_varinf(t_isolated_con);
transm_prev_con_ferretti = 1-F_gen_ferretti(t_isolated_con);

% Transmissions prevented when symptomatic hosts isolated at different
% times since infection

t_isolated_sym = chebfun('t',[0,7]);
transm_prev_sym_indep = 1-F_tost_indep(t_isolated_sym);
transm_prev_sym_constinf = 1-F_tost_constinf(t_isolated_sym);
transm_prev_sym_varinf = 1-F_tost_varinf(t_isolated_sym);
transm_prev_sym_ferretti = 1-F_tost_ferretti(t_isolated_sym);

% Obtain credible intervals for estimates at specific time points

load('../../Data/input_data.mat','k_inc','gamma','k_I','inc_shape','inc_scale')
load('../../Results/mcmc_posterior_indep','prob_presymp_post','gen_shape_post','gen_scale_post')
prob_presymp_post_indep = prob_presymp_post;
load('../../Results/mcmc_posterior_varinf','prob_presymp_post','theta_post')
prob_presymp_post_varinf = prob_presymp_post;

trans_prev_sym_1_indep_post = 1-F_tost_indep_vec(1,gen_shape_post,gen_scale_post,inc_shape,inc_scale);
trans_prev_sym_1_indep_CI = quantile(trans_prev_sym_1_indep_post,[0.025,0.975]);

trans_prev_sym_1_varinf_post = 1-F_tost_varinf_vec(1,theta_post,k_inc,gamma,k_I);
trans_prev_sym_1_varinf_CI = quantile(trans_prev_sym_1_varinf_post,[0.025,0.975]);

F_tost_m2_indep_post = F_tost_indep_vec(-2,gen_shape_post,gen_scale_post,inc_shape,inc_scale);
con_trac_indep_2_post = (prob_presymp_post_indep-F_tost_m2_indep_post)./prob_presymp_post_indep;
con_trac_indep_2_CI = quantile(con_trac_indep_2_post,[0.025,0.975]);

F_tost_m2_varinf_post = F_tost_varinf_vec(-2,theta_post,k_inc,gamma,k_I);
con_trac_varinf_2_post = (prob_presymp_post_varinf-F_tost_m2_varinf_post)./prob_presymp_post_varinf;
con_trac_varinf_2_CI = quantile(con_trac_varinf_2_post,[0.025,0.975]);

F_tost_m4_varinf_post = F_tost_varinf_vec(-4,theta_post,k_inc,gamma,k_I);
con_trac_varinf_4_post = (prob_presymp_post_varinf-F_tost_m4_varinf_post)./prob_presymp_post_varinf;
con_trac_varinf_4_CI = quantile(con_trac_varinf_4_post,[0.025,0.975]);

save('../../Results/contact_tracing','contacts_traced_indep','contacts_traced_constinf','contacts_traced_varinf','contacts_traced_ferretti','transm_prev_con_indep','transm_prev_con_constinf','transm_prev_con_varinf','transm_prev_con_ferretti')
save('../../Results/contact_tracing','transm_prev_sym_indep','transm_prev_sym_constinf','transm_prev_sym_varinf','transm_prev_sym_ferretti','-append')

rmpath('../../Functions/Indep')
rmpath('../../Functions/Mechanistic')
rmpath('../../Results')