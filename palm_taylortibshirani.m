function T = palm_taylortibshirani(G,df2,plm,c)
% Computes the combined statistic. It is calculated
% the fastest equivalent way as the original statistic.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Aug/2013
% http://brainder.org

P = 1-palm_gcdf(G,plm.tmp.rC(c),df2);
[~,tmp] = sort(P);
[~,prank] = sort(tmp);
T = sum(1-P.*(plm.nY+1)./prank)/plm.nY;