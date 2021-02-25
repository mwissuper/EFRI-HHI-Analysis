function [lagMax,rmax] = getXcorr(x,y,sample_rate)

xmr = x - nanmean(x);
ymr = y - nanmean(y);

[r, lags] = xcorr(xmr,ymr,sample_rate,'normalized'); % Constrain window for xcorr to +/- 1s
ind = find(abs(r) == max(abs(r)));
lagMax = lags(ind)/sample_rate;
rmax = r(ind);