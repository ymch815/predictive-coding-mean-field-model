% Fig 8H: barplot of Ave. ACC synaptic variable during BASELINE
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

itr=15;
itrmin=2;
itrmax=2;
rep=1;
[~,recell_n,~] = trials('naive', 'evoked', 'no feedback', {'baselineACC'}, plot_dynamics,...
    plot_comp_rs, itr, itrmin, itrmax, rep, noise);
[~,recell_c,~] = trials('chronic', 'evoked', 'no feedback', {'baselineACC'}, plot_dynamics,...
    plot_comp_rs, itr, itrmin, itrmax, rep, noise);
basearr_n=recell_n{1}; basearr_c=recell_c{1};
meanbasen=mean(basearr_n,1);
relerrn=std(basearr_n,1)/sqrt(size(basearr_n,1));
meanbasec=mean(basearr_c,1);
relerrc=std(basearr_c,1)/sqrt(size(basearr_c,1));
[hbar,~]=barwitherr([relerrn(1) relerrc(1)],[1 1.5],...
    [meanbasen(1) meanbasec(1)]);
hold on;
space=0.075;
plot(linspace(1-space,1+space,20),ones(20,1)*0.35,'k-','linewidth',2)
text(1-space/2,0.355,' ****','fontsize',20);
set(hbar(1),'FaceColor','flat','CData',[0 0.4470 0.7410;0.8500 0.3250 0.0980],...
    'EdgeColor','none','BarWidth',0.5);
set(gca,'xtick',[1,1.5],'TickLength',[0 0]);
set(gca,'xticklabel',{'Naive','Chronic'},'fontsize',20)
ylabel('ACC synaptic activ. (a.u.)','fontsize',20)
axis([0.5 2. 0 meanbasec(1)+relerrc(1)+0.02])
box on;grid off; set(gca,'linewidth',2)
% p1=anova1([basearr_n(:,1) basearr_c(:,1)]);