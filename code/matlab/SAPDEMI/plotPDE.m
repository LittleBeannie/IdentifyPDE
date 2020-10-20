function [] = plotPDE(xData,tNum, U, mark)
figure
for i = floor(linspace(1,tNum,5))
    if i == 1
        plot(xData,U(i,:),'Color',[0,0,0],'LineWidth',1.5)
    else
        plot(xData,U(i,:),'Color',[0 0.4470 0.7410,i/tNum],'LineWidth',1.5)
    end
    hold on 
end
ylim([-2,2.5])
xlabel('$x$','Interpreter','Latex')
ylabel(mark,'Interpreter','Latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',30)
