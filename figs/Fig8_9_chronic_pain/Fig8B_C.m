% Fig 8B&C: the dynamics of evoked pain, in naive and chronic conditions
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
plot_dynamics=0; % plot the dynamics or not
plot_comp_rs=0; % plot the correlation between r and s or not 

for xvalue=[3,4]
itr=1;
itrmin=xvalue;
itrmax=xvalue;
rep=1;
[~,~,dynamics_n] = trials('naive', 'evoked', 'no feedback', {}, plot_dynamics,...
    plot_comp_rs, itr, itrmin, itrmax, rep, noise);
[~,~,dynamics_c] = trials('chronic', 'evoked', 'no feedback', {}, plot_dynamics,...
    plot_comp_rs, itr, itrmin, itrmax, rep, noise);
figure; hold on
plot(t(pre),dynamics_n(1,pre),'b-','LineWidth',2)
plot(t(pre),dynamics_c(1,pre),'r-','LineWidth',2)
myplot('Time(ms)', 'ACC synaptic activ. (a.u.)', 20, {'Naive','Chronic'}, 20,'Northeast',...
    [pre(1)*dt pre(end)*dt min(min(dynamics_n(1,pre)),min(dynamics_c(1,pre)))-0.05...
    max(max(dynamics_n(1,pre)),max(dynamics_c(1,pre)))+0.05])
end