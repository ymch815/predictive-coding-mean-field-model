% Fig D&E: Ave. ACC synaptic variable under different z(0) in placebo
% effect

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
itrmin=-3;
itrmax=3;
rep=1;

% load data OR calculate data

% load data directly
data=load('chronic_placebo.mat');
placeboarr_n=data.placeboarr_n; placeboarr_c=data.placeboarr_c;
xarr=data.xarr;


% calculate data
% [~,recell_n,~] = trials('naive', 'placebo', 'no feedback', {'peakACC'}, plot_dynamics,...
%     plot_comp_rs, itr, itrmin, itrmax, rep, noise);
% [xarr,recell_c,~] = trials('chronic', 'placebo', 'no feedback', {'peakACC'}, plot_dynamics,...
%     plot_comp_rs, itr, itrmin, itrmax, rep, noise);
% placeboarr_n=recell_n{1}; placeboarr_c=recell_c{1};

% plot
grpsize=10;
pop_plot(xarr,{placeboarr_n,placeboarr_c},grpsize,'total')
myplot('Top-down expectation z(0)', 'Ave. ACC synaptic activ. (a.u.)', 20, {'Total (Naive)','Total (Chronic)'},...
    20,'Northeast',[xarr(1)-0.1 xarr(end)+0.1 min(placeboarr_n(:,1)) max(placeboarr_c(:,1))+0.1])
pop_plot(xarr,{placeboarr_n,placeboarr_c},grpsize,'subp')
myplot('Top-down expectation z(0)', 'Ave. ACC synaptic activ. (a.u.)', 20, {'w/ S1 (Naive)',...
    'w/ S1 (Chronic)','w/o S1 (Naive)','w/o S1 (Chronic)'},...
    20,'Northwest',[xarr(1)-0.1 xarr(end)+0.1 min(placeboarr_n(:,3)) max(placeboarr_c(:,2))])
% % % save('chronic_placebo.mat','xarr','placeboarr_n','placeboarr_c')