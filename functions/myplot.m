% this is the function for adjusting plotting formats

function myplot(xlabelname, ylabelname, fontsize, legendname, lgfontsize,lglocation,xylim)
box on;
grid off; 
xlabel(xlabelname);
ylabel(ylabelname);
legend(legendname,'fontsize',20);
set(legend,'FontSize',lgfontsize,'location',lglocation);
set(gca,'fontsize',fontsize);set(gca,'linewidth',2)
axis(xylim)
end