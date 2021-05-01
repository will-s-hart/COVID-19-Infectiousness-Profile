% Produce the panels in Fig. S2 of our manuscript.

% This script requires the Chebfun package to run (freely available at
% https://www.chebfun.org/download/), and requires the export_fig package
% (freely available at https://github.com/altmany/export_fig) to export the
% figure panels to pdf files.

clear all; close all; clc;
splitting on

addpath('../Results')

load('../Results/asymp_proportion_params','a_best','b_best')
load('../Results/asymp_infectiousness_params','mu_best','sigma_best')

p_A_range = [0,1];
p_A_prior = chebfun(@(t)betapdf(t,a_best,b_best),p_A_range);

x_A_range = [0,1.5];
x_A_prior = chebfun(@(t)lognpdf(t,mu_best,sigma_best),x_A_range);

% Colorscheme

c1 = [0, 0.4470, 0.7410]; c2 = [0.8500, 0.3250, 0.0980];
c3 = [0.9290, 0.6940, 0.1250]; c4 = [0.4940, 0.1840, 0.5560];

% Plot distributions

for k = 1:2
figsetup(k)
end

figure(1); hold on;
plot(p_A_prior,'color',c1,'linewidth',3)
xticks(0:0.25:1)
ylim([0,max(ylim)])
xlabel('Proportion of asymptomatic cases')
ylabel('Density')

figure(2); hold on;
plot(x_A_prior,'color',c1,'linewidth',3)
ylim([0,max(ylim)])
xlabel({'Relative infectiousness of';'asymptomatic hosts'})
ylabel('Density')

for k = 1:2
figsetup(k)
end

figure(1); export_fig Figures/figS2a.pdf -nocrop -transparent
figure(2); export_fig Figures/figS2b.pdf -nocrop -transparent

rmpath('../Results')