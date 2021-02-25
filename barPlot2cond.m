function h = barPlot2cond(metric,condLab,titlename,yLabName)

hold on;
h = bar(mean(metric),'linestyle','none');
errorbar(mean(metric),std(metric),'linestyle','none','color','k');
% boxplot(metric,'colors','k','symbol','o')
set(gca,'xticklabel',condLab,'xtick',1:2); box off; set(gca,'tickdir','out');
xlim([0.5 2.5]);
title(titlename),ylabel(yLabName); 

