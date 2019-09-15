% Fig S4: data of evoked, spontaneous and placebo conditions, 
% A-D, for evoked condition
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
itrmin=1.3;
itrmax=5;
rep=1;

% load data OR calculate data

% load data directly
data=load('supp_evo.mat');
recell_evo_n_max=data.recell_evo_n_max; 
recell_evo_n_lat=data.recell_evo_n_lat;
recell_evo_c_max=data.recell_evo_c_max; 
recell_evo_c_lat=data.recell_evo_c_lat;
xarr_evo=data.xarr_evo;

% calculate data
% [~,recell_evo_n,~] = trials('naive', 'evoked', 'no feedback', ...
%     {'maximum','latency'}, plot_dynamics, plot_comp_rs, itr, itrmin, itrmax, rep, noise);
% [xarr_evo,recell_evo_c,~] = trials('chronic', 'evoked', 'no feedback', ...
%     {'maximum','latency'}, plot_dynamics, plot_comp_rs, itr, itrmin, itrmax, rep, noise);
% recell_evo_n_max = recell_evo_n{1}; recell_evo_n_lat=recell_evo_n{2};
% recell_evo_c_max = recell_evo_c{1}; recell_evo_c_lat=recell_evo_c{2};


% plot: maximum
grpsize=10;
pop_plot(xarr_evo,{recell_evo_n_max,recell_evo_c_max},grpsize,'total')
myplot('Stimulus amplitude', 'Max ampl. (a.u.) of ACC syn. act.', 20, {'Total (Naive)','Total (Chronic)'},...
    20,'Northwest',[xarr_evo(1)-0.1 xarr_evo(end)+0.1 min(recell_evo_n_max(:,1)) max(recell_evo_c_max(:,1))+0.1])
pop_plot(xarr_evo,{recell_evo_n_max,recell_evo_c_max},grpsize,'subp')
myplot('Stimulus amplitude', 'Max ampl. (a.u.) of ACC syn. act.', 20, {'w/ S1 (Naive)',...
    'w/ S1 (Chronic)','w/o S1 (Naive)','w/o S1 (Chronic)'},...
    20,'Northwest',[xarr_evo(1)-0.1 xarr_evo(end)+0.1 min(recell_evo_n_max(:,3)) max(recell_evo_c_max(:,2))+0.2])
% plot: latency
pop_plot(xarr_evo,{recell_evo_n_lat,recell_evo_c_lat},grpsize,'total')
myplot('Stimulus amplitude', {'Latency from stim. onset','to max ACC syn. act (ms)'}, 20, {'Total (Naive)','Total (Chronic)'},...
    20,'Northeast',[xarr_evo(1)-0.1 xarr_evo(end)+0.1 min(recell_evo_n_lat(:,1)) max(recell_evo_c_lat(:,1))+0.1])
pop_plot(xarr_evo,{recell_evo_n_lat,recell_evo_c_lat},grpsize,'subp')
myplot('Stimulus amplitude', {'Latency from stim. onset','to max ACC syn. act (ms)'}, 20, {'w/ S1 (Naive)',...
    'w/ S1 (Chronic)','w/o S1 (Naive)','w/o S1 (Chronic)'},...
    20,'Northeast',[xarr_evo(1)-0.1 xarr_evo(end)+0.1 min(recell_evo_n_lat(:,3)) max(recell_evo_c_lat(:,2))])
% % % save('supp_evo','xarr_evo','recell_evo_n_max','recell_evo_n_lat',...
% % %     'recell_evo_c_max','recell_evo_c_lat')