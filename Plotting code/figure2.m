% Produce the panels in Fig. 2 of our manuscript.

% This script requires the Chebfun package to run (freely available at
% https://www.chebfun.org/download/), and requires the export_fig package
% (freely available at https://github.com/altmany/export_fig) to export the
% figure panels to pdf files.

clear all; close all; clc;
splitting on

addpath('../Data')
addpath('../Results')

load('../Data/input_data.mat','serial_data')

load('../Results/gen_tost_serial_indep.mat','f_gen_indep','f_tost_indep','f_serial_indep')
load('../Results/gen_tost_serial_constinf.mat','f_gen_constinf','f_tost_constinf','f_serial_constinf')
load('../Results/gen_tost_serial_varinf.mat','f_gen_varinf','f_tost_varinf','f_serial_varinf')
load('../Results/gen_tost_serial_ferretti.mat','f_gen_ferretti','f_tost_ferretti','f_serial_ferretti')

load('../Results/param_fit_ml_indep','ll_best')
ll_indep = ll_best;
load('../Results/param_fit_ml_constinf','ll_best')
ll_constinf = ll_best;
load('../Results/param_fit_ml_varinf','ll_best')
ll_varinf = ll_best;
load('../Results/param_fit_ml_ferretti','ll_best')
ll_ferretti = ll_best;

% AIC

AIC_indep = 2*(2-ll_indep);
AIC_constinf = 2*(2-ll_constinf);
AIC_varinf = 2*(3-ll_varinf);
AIC_ferretti = 2*(3-ll_ferretti);

AIC_best = min([AIC_indep,AIC_constinf,AIC_varinf,AIC_ferretti]);

dAIC_indep = AIC_indep - AIC_best;
dAIC_constinf = AIC_constinf - AIC_best;
dAIC_varinf = AIC_varinf - AIC_best;
dAIC_ferretti = AIC_ferretti - AIC_best;

% Colorscheme

c1 = [0, 0.4470, 0.7410]; c2 = [0.8500, 0.3250, 0.0980];
c3 = [0.9290, 0.6940, 0.1250]; c4 = [0.4940, 0.1840, 0.5560];

% Plot distributions

for k = 1:3
figsetup(k)
end

figure(1); hold on;
p2=plot(f_gen_constinf,'color',c2,'linewidth',3);
p1=plot(f_gen_varinf,'color',c1,'linewidth',3);
p3=plot(f_gen_ferretti,'-.','color',c3,'linewidth',2.5);
p4=plot(f_gen_indep,'-.','color',c4,'linewidth',2.5);
xlim([0,15])
ylim([0,max(ylim)])
xlabel('Generation time (days)')
ylabel('Density')
l = legend([p1,p2,p3,p4],{'Variable infectiousness','Constant infectiousness','Ferretti',['Independent transmission' newline 'and symptoms']});
l.FontSize = 14;
l.Position = [0.4510    0.6690    0.4070    0.1770];

figure(2); hold on;
plot(f_tost_constinf,'color',c2,'linewidth',3)
plot(f_tost_varinf,'color',c1,'linewidth',3)
plot(f_tost_ferretti,'-.','color',c3,'linewidth',2.5)
plot(f_tost_indep,'-.','color',c4,'linewidth',2.5)
xlim([-10,10]);
ylim([0,0.3])
xlabel({'Time from onset of symptoms';' to transmission (days)'})
ylabel('Density')

figure(3); hold on;
histogram(serial_data,'normalization','pdf','facecolor',0.8*[1,1,1],'edgecolor',0.8*[1,1,1])
plot(f_serial_constinf,'color',c2,'linewidth',3)
plot(f_serial_varinf,'color',c1,'linewidth',3)
plot(f_serial_ferretti,'-.','color',c3,'linewidth',2.5)
plot(f_serial_indep,'-.','color',c4,'linewidth',2.5)
xlim([-5,20]);
xlabel('Serial interval (days)')
ylabel('Density')

for k = 1:3
figsetup(k)
end

figure(1); export_fig Figures/fig2a.pdf -nocrop -transparent
figure(2); export_fig Figures/fig2b.pdf -nocrop -transparent
figure(3); export_fig Figures/fig2c.pdf -nocrop -painters -transparent

rmpath('../Data')
rmpath('../Results')