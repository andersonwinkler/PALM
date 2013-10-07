function T = palm_fristonnichols(G,df2,plm,c)
% Computes the combined statistic. It is calculated
% the fastest equivalent way as the original statistic.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Aug/2013
% http://brainder.org

% This will compute the stat in terms of min(1-P)
T = min(palm_gcdf(G,plm.tmp.rC(c),df2),[],1);
