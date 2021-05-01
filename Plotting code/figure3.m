% Produce the panels in Fig. 3 of our manuscript.

% This script requires the Violinplot-Matlab package to produce the violin
% plots (freely available at https://github.com/bastibe/Violinplot-Matlab),
% and requires the export_fig package (freely available at
% https://github.com/altmany/export_fig) to export the figure panels to pdf
% files.

clear all; close all; clc;

addpath('../Results')

load('../Results/mcmc_posterior_indep','prob_presymp_post')
prob_presymp_post_indep = prob_presymp_post;
load('../Results/mcmc_posterior_constinf','prob_presymp_post')
prob_presymp_post_constinf = prob_presymp_post;
load('../Results/mcmc_posterior_varinf','prob_presymp_post')
prob_presymp_post_varinf = prob_presymp_post;
load('../Results/mcmc_posterior_ferretti','prob_presymp_post')
prob_presymp_post_ferretti = prob_presymp_post;

% Colorscheme

c1 = [0, 0.4470, 0.7410]; c2 = [0.8500, 0.3250, 0.0980];
c3 = [0.9290, 0.6940, 0.1250]; c4 = [0.4940, 0.1840, 0.5560];

% Comparison of posterior proportions of presymptomatic transmission

for k = 1:2
figsetup(k)
end

figure(1); hold on;
v = violinplot(100*[prob_presymp_post_varinf,prob_presymp_post_constinf,prob_presymp_post_ferretti,prob_presymp_post_indep],{'','','',''},'ShowData',false,'ViolinAlpha',1);
v(1).ViolinColor = c1; v(2).ViolinColor = c2; v(3).ViolinColor = c3; v(4).ViolinColor = c4; 
v1 = v(1).ViolinPlot; v2 = v(2).ViolinPlot; v3 = v(3).ViolinPlot; v4 = v(4).ViolinPlot;
ylim([30,100])
ylabel({'Transmissions before';'symptom onset (%)'})
l = legend([v1,v2,v3,v4],{'Variable infectiousness','Constant infectiousness','Ferretti',['Independent transmission' newline 'and symptoms']});
l.FontSize = 14;

% Range of possible total proportions of non-symptomatic transmissions

load('../Results/asymp_param_samples','p_A_sample','x_A_sample');

prob_not_symp_sample_indep = (p_A_sample.*x_A_sample+(1-p_A_sample).*prob_presymp_post_indep)./(p_A_sample.*x_A_sample+(1-p_A_sample));
prob_not_symp_sample_constinf = (p_A_sample.*x_A_sample+(1-p_A_sample).*prob_presymp_post_constinf)./(p_A_sample.*x_A_sample+(1-p_A_sample));
prob_not_symp_sample_varinf = (p_A_sample.*x_A_sample+(1-p_A_sample).*prob_presymp_post_varinf)./(p_A_sample.*x_A_sample+(1-p_A_sample));
prob_not_symp_sample_ferretti = (p_A_sample.*x_A_sample+(1-p_A_sample).*prob_presymp_post_ferretti)./(p_A_sample.*x_A_sample+(1-p_A_sample));

figure(2); hold on;
v = violinplot(100*[prob_not_symp_sample_varinf,prob_not_symp_sample_constinf,prob_not_symp_sample_ferretti,prob_not_symp_sample_indep],{'','','',''},'ShowData',false,'ViolinAlpha',1);
v(1).ViolinColor = c1; v(2).ViolinColor = c2; v(3).ViolinColor = c3; v(4).ViolinColor = c4;
ylim([30,100])
ylabel({'Total non-symptomatic';'transmissions (%)'})

for k = 1:2
figsetup(k)
end

figure(1); export_fig Figures/fig3a.pdf -nocrop -painters -transparent
figure(2); export_fig Figures/fig3b.pdf -nocrop -painters -transparent

rmpath('../Results')