% Fig S4: data of evoked, spontaneous and placebo conditions, 
% I-L, for placebo condition
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
data=load('supp_plb.mat');
recell_plb_n_max=data.recell_plb_n_max; 
recell_plb_n_lat=data.recell_plb_n_lat;
recell_plb_c_max=data.recell_plb_c_max; 
recell_plb_c_lat=data.recell_plb_c_lat;
xarr_plb=data.xarr_plb;

% calculate data
% [~,recell_plb_n,~] = trials('naive', 'placebo', 'no feedback', ...
%     {'maximum','latency'}, plot_dynamics, plot_comp_rs, itr, itrmin, itrmax, rep, noise);
% [xarr_plb,recell_plb_c,~] = trials('chronic', 'placebo', 'no feedback', ...
%     {'maximum','latency'}, plot_dynamics, plot_comp_rs, itr, itrmin, itrmax, rep, noise);
% recell_plb_n_max = recell_plb_n{1}; recell_plb_n_lat=recell_plb_n{2};
% recell_plb_c_max = recell_plb_c{1}; recell_plb_c_lat=recell_plb_c{2};

% plot: maximum
grpsize=10;
pop_plot(xarr_plb,{recell_plb_n_max,recell_plb_c_max},grpsize,'total')
myplot('Top-down expectation z(0)', 'Max ampl. (a.u.) of ACC syn. act.', 20, {'Total (Naive)','Total (Chronic)'},...
    20,'Northwest',[xarr_plb(1)-0.1 xarr_plb(end)+0.1 min(recell_plb_n_max(:,1))-0.1 max(recell_plb_c_max(:,1))+0.1])
pop_plot(xarr_plb,{recell_plb_n_max,recell_plb_c_max},grpsize,'subp')
myplot('Top-down expectation z(0)', 'Max ampl. (a.u.) of ACC syn. act.', 20, {'w/ S1 (Naive)',...
    'w/ S1 (Chronic)','w/o S1 (Naive)','w/o S1 (Chronic)'},...
    20,'Northwest',[xarr_plb(1)-0.1 xarr_plb(end)+0.1 min(recell_plb_n_max(:,3))-0.1 max(recell_plb_c_max(:,2))+0.1])
% plot: latency
pop_plot(xarr_plb,{recell_plb_n_lat,recell_plb_c_lat},grpsize,'total')
myplot('Top-down expectation z(0)', 'Latency from stim. onset to max ACC syn. act (ms)', 20, {'Total (Naive)','Total (Chronic)'},...
    20,'Northeast',[xarr_plb(1)-0.1 xarr_plb(end)+0.1 min(recell_plb_c_lat(:,1))-20 max(recell_plb_c_lat(:,1))+20])
pop_plot(xarr_plb,{recell_plb_n_lat,recell_plb_c_lat},grpsize,'subp')
myplot('Top-down expectation z(0)', 'Latency from stim. onset to max ACC syn. act (ms)', 20, {'w/ S1 (Naive)',...
    'w/ S1 (Chronic)','w/o S1 (Naive)','w/o S1 (Chronic)'},...
20,'Northeast',[xarr_plb(1)-0.1 xarr_plb(end)+0.1 min(recell_plb_c_lat(:,3))-20 max(recell_plb_c_lat(:,3))+20])
% save('supp_plb','xarr_plb','recell_plb_n_max','recell_plb_n_lat',...
%     'recell_plb_c_max','recell_plb_c_lat')