% Fig 10A & 10B: Ave S1 synaptic variable before / during stimulus under
% different expectation

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
itrmin=1;
itrmax=5;
rep=1;

% load data OR calculate data

% load data directly
beforedata=load('beforeS1.mat');
duringdata=load('duringS1.mat');
xarr = beforedata.xarr; 
beforeS1 = beforedata.beforeS1; duringS1_2 = duringdata.duringS1_2; 
duringS1_25 = duringdata.duringS1_25; duringS1_3 = duringdata.duringS1_3;

% calculate data
% [~,recell_2,~] = trials('naive', 'prediction_const=2', 'no feedback', ...
%     {'beforeS1','duringS1'}, plot_dynamics, plot_comp_rs, itr, itrmin, itrmax, rep, noise);
% [~,recell_25,~] = trials('naive', 'prediction_const=2.5', 'no feedback', ...
%     {'duringS1'}, plot_dynamics, plot_comp_rs, itr, itrmin, itrmax, rep, noise);
% [xarr,recell_3,~] = trials('naive', 'prediction_const=3', 'no feedback', ...
%     {'duringS1'}, plot_dynamics, plot_comp_rs, itr, itrmin, itrmax, rep, noise);
% beforeS1 = recell_2{1}; duringS1_2 = recell_2{2}; 
% duringS1_25 = recell_25{1}; duringS1_3 = recell_3{1}; 

% plot
grpsize=10;
pop_plot(xarr,{beforeS1},grpsize,'other')
myplot('Top-down expectation z(0)', 'Ave. S1 synaptic activ. (a.u.)', 20, {'before onset'},...
    20,'Northeast',[xarr(1)-0.1 xarr(end)+0.1 min(beforeS1(:,1))-0.02 max(beforeS1(:,1))+0.02])
pop_plot(xarr,{duringS1_2,duringS1_25,duringS1_3},grpsize,'other')
myplot('Top-down expectation z(0)', 'Ave. S1 synaptic activ. (a.u.)', 20, {'x=2.0','x=2.5','x=3.0'},...
    20,'Northwest',[xarr(1)-0.1 xarr(end)+0.1 min(duringS1_2(:,1))-0.02 max(duringS1_2(:,1))+0.02])
% % % save('beforeS1.mat','xarr','beforeS1')
% % % save('duringS1.mat','duringS1_2','duringS1_25','duringS1_3')