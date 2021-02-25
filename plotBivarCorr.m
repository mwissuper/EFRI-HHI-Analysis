function [h,p,r,test] = plotBivarCorr(v1,v2,subj_array,subj_array_force,colors,temp,xname,yname)

for i = 1:length(subj_array_force) % must plot one at a time to get diff color point
    subj = subj_array_force(i);
    indC = find(subj_array == subj,1,'first'); % Use same color as kinem data
    plot(v1(i),v2(i),'.','markersize',14,'color',colors(indC,:));
end

[p,r,test] = corr2vars(v1(:),v2); % Checks if data is normal and does appropriate test
if p < 0.05
    c = polyfit(v1(:),v2,1);
    plot(v1(:),polyval(c,v1(:)),'k--');
    titlename = sprintf('%s \np = %.2f, rho = %.2f',temp,p,r); 
else
    titlename = sprintf('%s \np = %.2f',temp,p); 
end
title(titlename),box off;
xlabel(xname);ylabel(yname);
h = gcf;