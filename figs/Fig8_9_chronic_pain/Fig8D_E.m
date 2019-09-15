% Fig D&E: Ave. ACC synaptic variable with different stimulus, in naive and
% chronic conditions
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
data=load('chronic_peak.mat');
peakarr_n=data.peakarr_n; peakarr_c=data.peakarr_c; xarr=data.xarr;

% calculate data
% [~,recell_n,~] = trials('naive', 'evoked', 'no feedback', {'peakACC'}, plot_dynamics,...
%     plot_comp_rs, itr, itrmin, itrmax, rep, noise);
% [xarr,recell_c,~] = trials('chronic', 'evoked', 'no feedback', {'peakACC'}, plot_dynamics,...
%     plot_comp_rs, itr, itrmin, itrmax, rep, noise);
% peakarr_n=recell_n{1}; peakarr_c=recell_c{1};

% plot

grpsize=10;
pop_plot(xarr,{peakarr_n,peakarr_c},grpsize,'total')
myplot('Stimulus Amplitude', 'Ave. ACC synaptic activ. (a.u.)', 20, {'Total (Naive)','Total (Chronic)'},...
    20,'Northeast',[xarr(1)-0.1 xarr(end)+0.1 min(peakarr_n(:,1)) max(peakarr_c(:,1))+0.1])
pop_plot(xarr,{peakarr_n,peakarr_c},grpsize,'subp')
myplot('Stimulus Amplitude', 'Ave. ACC synaptic activ. (a.u.)', 20, {'w/ S1 (Naive)',...
    'w/ S1 (Chronic)','w/o S1 (Naive)','w/o S1 (Chronic)'},...
    20,'Northwest',[xarr(1)-0.1 xarr(end)+0.1 min(peakarr_n(:,3)) max(peakarr_c(:,2))])
% % % save('chronic_peak.mat','xarr','peakarr_n','peakarr_c')