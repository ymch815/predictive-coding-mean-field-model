% Fig 8G: barplot for comparison of Ave. ACC synaptic variable of subpopulations
% under high stimulus

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
itrmin=4.8;
itrmax=5.2;
rep=1;
[~,recell_n,~] = trials('naive', 'evoked', 'no feedback', {'peakACC'}, plot_dynamics,...
    plot_comp_rs, itr, itrmin, itrmax, rep, noise);
[~,recell_c,~] = trials('chronic', 'evoked', 'no feedback', {'peakACC'}, plot_dynamics,...
    plot_comp_rs, itr, itrmin, itrmax, rep, noise);
peakarr_n=recell_n{1}; peakarr_c=recell_c{1};
bar_plot(peakarr_n,peakarr_c,0.08,[0.76,0.55],{' **','****'},...
    {'w/ S1 (N)','w/ S1 (C)','w/o S1 (N)','w/o S1 (C)'},'ACC synaptic activ. (a.u.)',{[0 0.4470 0.7410;...
    0.9290 0.6940 0.1250],[0.8500 0.3250 0.0980;0.4940 0.1840 0.5560]},'chronic')
set(gca,'fontsize',15)
axis([0.8 1.7 0 0.8])
% p1=anova1([peakarr_n(:,2) peakarr_c(:,2)]);
% p2=anova1([peakarr_n(:,3) peakarr_c(:,3)]);