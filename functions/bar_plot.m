% This function generate barplot for our arrays
function bar_plot(arr_n,arr_c,space,position,textname,xlabelname,ylabelname,color,type)
meann=mean(arr_n,1);
relerrn=std(arr_n,1)/sqrt(size(arr_n,1));
meanc=mean(arr_c,1);
relerrc=std(arr_c,1)/sqrt(size(arr_c,1));
if strcmp('chronic',type)
[hbar,~]=barwitherr([relerrn(2) relerrc(2);relerrn(3) relerrc(3)],[1 1.5],...
    [meann(2) meanc(2); meann(3) meanc(3)]);
hold on;
plot(linspace(1-space,1+space,20),ones(20,1)*position(1),'k-','linewidth',2)
plot(linspace(1.5-space,1.5+space,20),ones(20,1)*position(2),'k-','linewidth',2)
text(1-space/2,position(1)+0.02,textname{1},'fontsize',20);
text(1.5-space/2,position(2)+0.02,textname{2},'fontsize',20);
set(gca,'xtick',[1-space 1+space,1.5-space 1.5+space],'TickLength',[0 0]);
set(gca,'xticklabel',xlabelname)
elseif strcmp('prediction',type)
[hbar,~]=barwitherr([relerrn(1) relerrc(1);relerrn(2) relerrc(2)],[1 1.5],...
    [meann(1) meanc(1); meann(2) meanc(2)]);
hold on;
plot(linspace(1-space,1+space,20),ones(20,1)*position(1),'k-','linewidth',2)
text(1-space/2,position(1)+0.02,textname{1},'fontsize',20);
plot(linspace(1.5-space,1.5+space,20),ones(20,1)*position(2),'k-','linewidth',2)
text(1.5-space/2,position(2)+0.02,textname{2},'fontsize',20);
plot(linspace(1-space,1.5-space,20),ones(20,1)*position(3),'k-','linewidth',2)
text(1.2,position(3)+0.02,textname{3},'fontsize',20);
plot(linspace(1+space,1.5+space,20),ones(20,1)*position(4),'k-','linewidth',2)
text(1.3,position(4)+0.02,textname{4},'fontsize',20);
set(gca,'xtick',[1 1.5],'TickLength',[0 0]);
set(gca,'xticklabel',xlabelname)
legend({'before onset','after onset'})
end
set(hbar(1),'FaceColor','flat','EdgeColor','none','CData',color{1});
set(hbar(2),'FaceColor','flat','EdgeColor','none','CData',color{2});
ylabel(ylabelname);
box on; grid off; 
set(gca,'fontsize',20);set(gca,'linewidth',2)
end