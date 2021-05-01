% Produce the panels in Fig. S1 of our manuscript.

% This script requires the Chebfun package to run (freely available at
% https://www.chebfun.org/download/), and requires the export_fig package
% (freely available at https://github.com/altmany/export_fig) to export the
% figure panels to pdf files.

clear all; close all; clc;
splitting on

addpath('../Data')
addpath('../Results')
addpath('../Functions/General')

load('../Data/input_data.mat','serial_data')

load('../Results/gen_tost_serial_indep.mat','f_serial_indep')
load('../Results/gen_tost_serial_constinf.mat','f_serial_constinf')
load('../Results/gen_tost_serial_varinf.mat','f_serial_varinf')
load('../Results/gen_tost_serial_ferretti.mat','f_serial_ferretti')

% Discretise serial interval distributions

t_discr = (-5):20;
serial_dist_discr_indep = discretise_serial_distrib(f_serial_indep,t_discr);
serial_dist_discr_constinf = discretise_serial_distrib(f_serial_constinf,t_discr);
serial_dist_discr_varinf = discretise_serial_distrib(f_serial_varinf,t_discr);
serial_dist_discr_ferretti = discretise_serial_distrib(f_serial_ferretti,t_discr);

% Colorscheme

c1 = [0, 0.4470, 0.7410]; c2 = [0.8500, 0.3250, 0.0980];
c3 = [0.9290, 0.6940, 0.1250]; c4 = [0.4940, 0.1840, 0.5560];

% Plot distributions

for k = 1:4
figsetup(k)
end

figure(1); hold on;
h = histogram(serial_data,'normalization','pdf','facecolor',0.8*[1,1,1],'edgecolor',0.8*[1,1,1]);
p1 = plot(serial_dist_discr_varinf.t,serial_dist_discr_varinf.p,'*','color',c1,'linewidth',3,'markersize',15);
xlim([-5,20]);
xlabel('Serial interval (days)')
ylabel('Probability')
l = legend([h,p1],'Data','Variable infectiousness');
l.FontSize = 14;
l.Position = [0.4700    0.7230    0.3700    0.0750];

figure(2); hold on;
h = histogram(serial_data,'normalization','pdf','facecolor',0.8*[1,1,1],'edgecolor',0.8*[1,1,1]);
p2 = plot(serial_dist_discr_constinf.t,serial_dist_discr_constinf.p,'*','color',c2,'linewidth',3,'markersize',15);
xlim([-5,20]);
xlabel('Serial interval (days)')
ylabel('Probability')
l = legend([h,p2],'Data','Constant infectiousness');
l.FontSize = 14;
l.Position = [0.4560    0.7230    0.3840    0.0750];

figure(3); hold on;
h = histogram(serial_data,'normalization','pdf','facecolor',0.8*[1,1,1],'edgecolor',0.8*[1,1,1]);
p3 = plot(serial_dist_discr_ferretti.t,serial_dist_discr_ferretti.p,'*','color',c3,'linewidth',3,'markersize',15);
xlim([-5,20]);
xlabel('Serial interval (days)')
ylabel('Probability')
l = legend([h,p3],'Data','Ferretti');
l.FontSize = 14;
l.Position = [0.6670    0.7230    0.1710    0.0750];

figure(4); hold on;
h = histogram(serial_data,'normalization','pdf','facecolor',0.8*[1,1,1],'edgecolor',0.8*[1,1,1]);
p4 = plot(serial_dist_discr_indep.t,serial_dist_discr_indep.p,'*','color',c4,'linewidth',3,'markersize',15);
xlim([-5,20]);
xlabel('Serial interval (days)')
ylabel('Probability')
l = legend([h,p4],'Data',['Independent transmission' newline 'and symptoms']);
l.FontSize = 14;
l.Position = [0.4310    0.6890    0.4070    0.1090];

for k = 1:4
figsetup(k)
end

figure(1); export_fig Figures/figS1a.pdf -nocrop -painters -transparent
figure(2); export_fig Figures/figS1b.pdf -nocrop -painters -transparent
figure(3); export_fig Figures/figS1c.pdf -nocrop -painters -transparent
figure(4); export_fig Figures/figS1d.pdf -nocrop -painters -transparent

rmpath('../Data')
rmpath('../Results')
rmpath('../Functions/General')