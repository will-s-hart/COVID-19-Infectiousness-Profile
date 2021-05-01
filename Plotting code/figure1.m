% Produce the panels in Fig. 1B of our manuscript.

% This script requires the Chebfun package to run (freely available at
% https://www.chebfun.org/download/), and requires the export_fig package
% (freely available at https://github.com/altmany/export_fig) to export the
% figure panels to pdf files.

clear all; close all; clc;

addpath('../Results')

rng(3);

load('../Results/mcmc_posterior_varinf','ll_best','params_point')
params = params_point;

load('../Results/gen_tost_serial_indep.mat','f_gen_indep')

% Colorscheme

c1 = [0, 0.4470, 0.7410];
c6 = [0.3010, 0.7450, 0.9330];

% Plot figure panels

for k_inc = 1:2
figsetup(k_inc)
end

figure(1); hold on;
plot(f_gen_indep,'color',c1,'linewidth',3)
xlim([0,15])
ylim([0,max(ylim)])
xticks([]);
yticks([]);
xlabel('Time since infection')
ylabel('Infectiousness')

gamma = params(1); mu = params(2); k_inc = params(3); k_E = params(4); k_P = k_inc-k_E;
k_I = params(5); alpha = params(6);

y_E_ex = gamrnd(k_E,1/(k_inc*gamma));
y_P_ex = gamrnd(k_P,1/(k_inc*gamma));
y_I_ex = gamrnd(k_I,1/(k_I*mu));

figure(2); hold on;
stairs([0,y_E_ex,y_E_ex+y_P_ex,y_E_ex+y_P_ex+y_I_ex,12],[0,alpha,1,0,0],'color',c1,'linewidth',3)
stairs([0,y_E_ex,y_E_ex+y_P_ex,y_E_ex+y_P_ex+y_I_ex,12],[0,((alpha*y_P_ex+y_I_ex)/(y_P_ex+y_I_ex))*[1,1],0,0],'--','color',c6,'linewidth',3)
xlim([0,12])
ylim([0,2.5])
xticks([]);
yticks([]);
xlabel('Time since infection')
ylabel('Infectiousness')

for k_inc = 1:2
figsetup(k_inc)
end

set(gca, 'Layer', 'bottom')

figure(1); export_fig Figures/fig1a.pdf -transparent
figure(2); export_fig Figures/fig1b.pdf -transparent

rmpath('../Results')