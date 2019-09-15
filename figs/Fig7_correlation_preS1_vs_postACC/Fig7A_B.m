% Fig 7 A and B: dynamics of evoked pain 
clear all 
close all
warning off 

global base T dt pre N

base = 2000; % baseline time length (at least 2000)
simul_t = 2500; % simulation time after baseline
T = simul_t+base; % Total simulation time in ms
dt = 0.1;    % Timestep for Euler's Method
N = T/dt; 
noise=1;
pre = (base+200)/dt:(1700+base)/dt;% time to present when plotting
plot_dynamics=1;
plot_comp_rs=0;
itr=1;
itrmin=2;
itrmax=2;
rep=1;
[~,~,~] = trials('naive', 'evoked', 'no feedback', {'preS1'}, plot_dynamics,...
    plot_comp_rs, itr, itrmin, itrmax, rep, noise);