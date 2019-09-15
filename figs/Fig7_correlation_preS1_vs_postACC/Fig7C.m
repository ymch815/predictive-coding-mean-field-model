% Fig 7C: correlation between preS1 and post ACC (evoked)
clear all 
close all
warning off 

global base T dt pre N

base = 2000; % baseline time length (at least 2000)
simul_t = 2500; % simulation time after baseline
T = simul_t+base; % Total simulation time in ms
dt = 0.1;    % Timestep for Euler's Method
N = T/dt; 
pre = (base+200)/dt:(1700+base)/dt;% time to present when plotting
noise=0;
plot_dynamics=0;
plot_comp_rs=0;
itr=50;
itrmin=1.3;
itrmax=3;
rep=1;

% load data OR calculate data

% load directly
data = load('cor_evoked.mat');
xarr = data.xarr;
preS1 = data.preS1; 
postACC = data.postACC; 

% calculate data
% [xarr,recell,~] = trials('naive', 'evoked', 'no feedback', {'preS1','postACC'}, plot_dynamics,...
%     plot_comp_rs, itr, itrmin, itrmax, rep, noise);
% preS1 = recell{1}(:,1); postACC = recell{2}(:,1);

% plot
figure; hold on; box on; grid off;
scatter(preS1,postACC,40,xarr,'filled')
xlabel('Ave. s (pre-S1)'); ylabel('Ave. s (post-ACC)')
set(gca,'fontsize',20);set(gca,'linewidth',2)
axis([min(preS1)-0.002 max(preS1)+0.002 min(postACC)-0.002 max(postACC)+0.002])
h = colorbar;
set(get(h,'label'),'string','x(0)','fontsize',20); 
[rho,pval]=corr(preS1,postACC);
fprintf("the correlation is: %4.2f \n p value is: %4.4f \n",rho,pval)
% % % save('cor_evoked.mat','xarr','preS1','postACC')