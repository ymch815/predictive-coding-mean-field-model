% Fig S4: data of evoked, spontaneous and placebo conditions, 
% E-H, for spontaneous condition
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
data=load('supp_spon.mat');
recell_spon_n_peak=data.recell_spon_n_peak; 
recell_spon_n_max=data.recell_spon_n_max;
recell_spon_c_peak=data.recell_spon_c_peak; 
recell_spon_c_max=data.recell_spon_c_max;
xarr_spon=data.xarr_spon;

% calculate data
% [~,recell_spon_n,~] = trials('naive', 'spontaneous', 'no feedback', ...
%     {'peakACC','maximum'}, plot_dynamics, plot_comp_rs, itr, itrmin, itrmax, rep, noise);
% [xarr_spon,recell_spon_c,~] = trials('chronic', 'spontaneous', 'no feedback', ...
%     {'peakACC','maximum'}, plot_dynamics, plot_comp_rs, itr, itrmin, itrmax, rep, noise);
% recell_spon_n_peak = recell_spon_n{1}; recell_spon_n_max=recell_spon_n{2};
% recell_spon_c_peak = recell_spon_c{1}; recell_spon_c_max=recell_spon_c{2};

% plot: ave ACC
grpsize=10;
pop_plot(xarr_spon,{recell_spon_n_peak,recell_spon_c_peak},grpsize,'total')
myplot('Top-down expectation z(0)', 'Ave. ACC synaptic activ. (a.u.)', 20, {'Total (Naive)','Total (Chronic)'},...
    20,'Northwest',[xarr_spon(1)-0.1 xarr_spon(end)+0.1 min(recell_spon_n_peak(:,1)) max(recell_spon_c_peak(:,1))+0.1])
pop_plot(xarr_spon,{recell_spon_n_peak,recell_spon_c_peak},grpsize,'subp')
myplot('Top-down expectation z(0)', 'Ave. ACC synaptic activ. (a.u.)', 20, {'w/ S1 (Naive)',...
    'w/ S1 (Chronic)','w/o S1 (Naive)','w/o S1 (Chronic)'},...
    20,'Northwest',[xarr_spon(1)-0.1 xarr_spon(end)+0.1 min(recell_spon_n_peak(:,3)) max(recell_spon_c_peak(:,2))+0.1])
% plot: maximum
pop_plot(xarr_spon,{recell_spon_n_max,recell_spon_c_max},grpsize,'total')
myplot('Top-down expectation z(0)', 'Max ampl. (a.u.) of ACC syn. act', 20, {'Total (Naive)','Total (Chronic)'},...
    20,'Northwest',[xarr_spon(1)-0.1 xarr_spon(end)+0.1 min(recell_spon_n_max(:,1)) max(recell_spon_c_max(:,1))+0.1])
pop_plot(xarr_spon,{recell_spon_n_max,recell_spon_c_max},grpsize,'subp')
myplot('Top-down expectation z(0)', 'Max ampl. (a.u.) of ACC syn. act', 20, {'w/ S1 (Naive)',...
    'w/ S1 (Chronic)','w/o S1 (Naive)','w/o S1 (Chronic)'},...
20,'Northwest',[xarr_spon(1)-0.1 xarr_spon(end)+0.1 min(recell_spon_n_max(:,3)) max(recell_spon_c_max(:,2))+0.1])
% % % save('supp_spon','xarr_spon','recell_spon_n_peak','recell_spon_n_max',...
% % %     'recell_spon_c_peak','recell_spon_c_max')