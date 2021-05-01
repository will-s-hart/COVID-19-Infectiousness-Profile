% Produce the panels in Fig. 4 of our manuscript.

% This script requires the Chebfun package to run (freely available at
% https://www.chebfun.org/download/), and requires the export_fig package
% (freely available at https://github.com/altmany/export_fig) to export the
% figure panels to pdf files.

clear all; close all; clc;

addpath('../Results')

load('../Results/contact_tracing','contacts_traced_indep','contacts_traced_constinf','contacts_traced_varinf','contacts_traced_ferretti','transm_prev_con_indep','transm_prev_con_constinf','transm_prev_con_varinf','transm_prev_con_ferretti')
load('../Results/contact_tracing','transm_prev_sym_indep','transm_prev_sym_constinf','transm_prev_sym_varinf','transm_prev_sym_ferretti')

% Colorscheme

c1 = [0, 0.4470, 0.7410]; c2 = [0.8500, 0.3250, 0.0980];
c3 = [0.9290, 0.6940, 0.1250]; c4 = [0.4940, 0.1840, 0.5560];

% Plot figures

for k = 1:3
figsetup(k)
end

% Isolation of symptomatic hosts

figure(1); hold on;
p2=plot(100*transm_prev_sym_constinf,'color',c2,'linewidth',3);
p1=plot(100*transm_prev_sym_varinf,'color',c1,'linewidth',3);
p3=plot(100*transm_prev_sym_ferretti,'-.','color',c3,'linewidth',2.5);
p4=plot(100*transm_prev_sym_indep,'-.','color',c4,'linewidth',2.5);
xticks(0:7);
ylim([0,50])
xlabel({'Time from symptom onset to isolation';'(days)'})
ylabel('Transmissions prevented (%)')
l = legend([p1,p2,p3,p4],{'Variable infectiousness','Constant infectiousness','Ferretti',['Independent transmission' newline 'and symptoms']});
l.FontSize = 14;

% Time to go back when contact tracing

figure(2); hold on;
p2=plot(100*contacts_traced_constinf,'color',c2,'linewidth',3);
p1=plot(100*contacts_traced_varinf,'color',c1,'linewidth',3);
p3=plot(100*contacts_traced_ferretti,'-.','color',c3,'linewidth',2.5);
p4=plot(100*contacts_traced_indep,'-.','color',c4,'linewidth',2.5);
xticks(0:7);
ylim([0,100])
xlabel({'Contact elicitation window before';'symptom onset (days)'})
ylabel({'Presymptomatic infectious contacts';'found (%)'})

% Isolation of infected contacts

figure(3); hold on;
p2=plot(100*transm_prev_con_constinf,'color',c2,'linewidth',3);
p1=plot(100*transm_prev_con_varinf,'color',c1,'linewidth',3);
p3=plot(100*transm_prev_con_ferretti,'-.','color',c3,'linewidth',2.5);
p4=plot(100*transm_prev_con_indep,'-.','color',c4,'linewidth',2.5);
xticks(0:2:10);
ylim([0,100])
xlabel('Time from exposure to isolation (days)')
ylabel('Onward transmissions prevented (%)')

for k = 1:3
figsetup(k)
end

figure(1); export_fig Figures/fig4a.pdf -nocrop -transparent
figure(2); export_fig Figures/fig4b.pdf -nocrop -transparent
figure(3); export_fig Figures/fig4c.pdf -nocrop -transparent

rmpath('../Results')