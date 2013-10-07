function T = palm_pearsondavid(G,df2,plm,c)
% Computes the combined statistic. It is calculated
% the fastest equivalent way as the original statistic.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Aug/2013
% http://brainder.org

P = palm_gcdf(G,plm.tmp.rC(c),df2);
T = -2*min(sum(log(1-P),1),sum(log(P),1));
