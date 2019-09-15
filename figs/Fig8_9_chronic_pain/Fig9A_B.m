% Fi9 A&B: dynamics of spontaneous pain trial under naive and chronic
% conditons, with different expectation level z(0)

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

% Fig A and B: comparison of spon
for zvalue=[2,3]
itr=1;
itrmin=zvalue;
itrmax=zvalue;
rep=1;
[~,~,dynamics_n] = trials('naive', 'spontaneous', 'no feedback', {}, plot_dynamics,...
    plot_comp_rs, itr, itrmin, itrmax, rep, noise);
[~,~,dynamics_c] = trials('chronic', 'spontaneous', 'no feedback', {}, plot_dynamics,...
    plot_comp_rs, itr, itrmin, itrmax, rep, noise);
figure; hold on
plot(t(pre),dynamics_n(1,pre),'b-','LineWidth',2)
plot(t(pre),dynamics_c(1,pre),'r-','LineWidth',2)
plot(linspace(t(pre(1)),t(pre(end)),100),ones(1,100)*0.91,'k--')
myplot('Time(ms)', 'ACC synaptic activ. (a.u.)', 20, {'Naive','Chronic'}, 20,'Northeast',...
    [pre(1)*dt pre(end)*dt min(min(dynamics_n(1,pre)),min(dynamics_c(1,pre)))-0.05...
    max(max(dynamics_n(1,pre)),max(dynamics_c(1,pre)))+0.05])
end
