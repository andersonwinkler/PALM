function T = palm_tippett(G,df2,plm,c)
% Computes the combined statistic. It is calculated
% the fastest equivalent way as the original statistic.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Aug/2013
% http://brainder.org

% Compute in terms of max(1-P), as each stat may have different df2.
T = max(palm_gcdf(G,plm.tmp.rC(c),df2),[],1);