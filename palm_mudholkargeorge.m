function T = palm_mudholkargeorge(G,df2,plm,c)
% Computes the combined statistic. It is calculated
% the fastest equivalent way as the original statistic.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Aug/2013
% http://brainder.org

P = palm_gcdf(G,plm.tmp.rC(c),df2);
mhcte = sqrt(3*(5*plm.nY+4)/plm.nY/(5*plm.nY+2))/pi;
T = mhcte*sum(log(P./(1-P)),1);