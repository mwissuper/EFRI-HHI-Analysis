function [perAposBpos,perAnegBpos,perAposBneg,perAnegBneg,meanABpos,meanABneg] = signProdPer(A,B)

indPosA = find(A >= 0);
indNegA = find(A < 0);
indAposBpos = find(B(indPosA,1) >= 0);
indAnegBpos = find(B(indNegA,1) >= 0);
indAposBneg = find(B(indPosA,1) < 0);
indAnegBneg = find(B(indNegA,1) < 0);
% Calculate percent of data with sign convention
perAposBpos = length(indAposBpos)/length(A);
perAnegBpos = length(indAnegBpos)/length(A);
perAposBneg = length(indAposBneg)/length(A);
perAnegBneg = length(indAnegBneg)/length(A);
% Calculate mean pos and neg power
P = A.*B;
indPpos = find(P>=0);
indPneg = find(P<0);
meanABpos = nanmean(P(indPpos));
meanABneg = nanmean(P(indPneg));