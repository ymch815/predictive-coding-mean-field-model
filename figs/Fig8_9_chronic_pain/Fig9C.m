% Fig 9C: fraction sustain time with different z(0)

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

itr=100;
itrmin=0.6;
itrmax=3.6;
rep=1;

% load data OR calculate data

% load data directly
data=load('chronic_sustain.mat');
sustainarr_n=data.sustainarr_n; sustainarr_c=data.sustainarr_c; xarr=data.xarr;

% calculate data
% [~,recell_n,~] = trials('naive', 'spontaneous', 'no feedback', {'sustain'}, plot_dynamics,...
%     plot_comp_rs, itr, itrmin, itrmax, rep, noise);
% [xarr,recell_c,~] = trials('chronic', 'spontaneous', 'no feedback', {'sustain'}, plot_dynamics,...
%     plot_comp_rs, itr, itrmin, itrmax, rep, noise);
% sustainarr_n=recell_n{1}; sustainarr_c=recell_c{1};

% plot
grpsize=10;
pop_plot(xarr,{sustainarr_n,sustainarr_c},grpsize,'total')
myplot('Top-down expectation z(0)', 'Ave. ACC synaptic activ. (a.u.)', 20, {'Total (Naive)','Total (Chronic)'},...
    20,'Northeast',[xarr(1)-0.1 xarr(end)+0.1 min(sustainarr_n(:,1))-0.1 max(sustainarr_c(:,1))+0.1])
% % % save('chronic_sustain.mat','xarr','sustainarr_n','sustainarr_c')