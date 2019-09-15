% Fig 10C: barplot of Ave S1 synaptic variable before and during stimulus,
% under conditons of precise prediction (PE=0) or wrong predicton (PE>0). 

clear all 
close all
warning off 

global base T dt pre N

base = 2000; % baseline time length (at least 2000)
simul_t = 2500; % simulation time after baseline
T = simul_t+base; % Total simulation time in ms
dt = 0.1;    % Timestep for Euler's Method
N = T/dt; 
t = linspace(0,T,N+1);
noise=1;
pre = (base-200)/dt:(1300+base)/dt;% time to present when plotting
plot_dynamics=0;
plot_comp_rs=0;

itr=15;
itrmin=1.8;
itrmax=2.2;
rep=1;
[~,recell_wro,~] = trials('naive', 'prediction_wrong', 'no feedback', ...
    {'beforeS1','duringS1'}, plot_dynamics, plot_comp_rs, itr, itrmin, itrmax, rep, noise);
[~,recell_pre,~] = trials('naive', 'prediction_precise','no feedback' ,...
    {'beforeS1','duringS1'}, plot_dynamics,plot_comp_rs, itr, itrmin, itrmax, rep, noise);
wro_before = recell_wro{1}(:,1); wro_dur = recell_wro{2}(:,1);
pre_before = recell_pre{1}(:,1); pre_dur = recell_pre{2}(:,1);
bar_plot([wro_before wro_dur],[pre_before pre_dur],0.08,[0.40,0.43,0.45,0.49],...
    {'****','****','****','****'},{'w/o Prediction','w/ Prediction'},...
    'Ave. S1 synaptic activ. (a.u.)',{[0 0.6 0.6;0 0.6 0.6],[0.8 0.0 0.4;0.8 0.0 0.4]},'prediction')
set(gca,'fontsize',15)
axis([0.75 1.75 0 0.6])
% p1=anova1([beforearr_e s1duringarr_e]);
% p2=anova1([beforearr_p s1duringarr_p]);
% p3=anova1([beforearr_e beforearr_p]);
% p4=anova1([s1duringarr_e s1duringarr_p]);