function T = palm_wilkinson(G,df2,plm,c)
% Computes the combined statistic. It is calculated
% the fastest equivalent way as the original statistic.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Aug/2013
% http://brainder.org

T = sum((1-palm_gcdf(G,plm.tmp.rC(c),df2)) <= plm.npcparm);