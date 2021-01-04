% Test zscore and standardization of regression coefficients from 
% https://www.mathworks.com/matlabcentral/answers/369305-how-to-standardize-unstandardized-beta-coefficients

% create IVs x1 and x2 and DV y as unstandardized variables 
x1=rand(100,1)*30;
x2=rand(100,1)*2;
y=rand(100,1)*10;
% create column of ones
X=[ones(size(x1)) x1 x2];
% perform regression on these
beta1=regress(y,X)
% standardize variables and perform regression again
X=[ones(size(x1)) zscore(x1) zscore(x2)]; % attach ones
beta2=regress(y,X)


beta1(2)*std(x1)
beta2(2)

