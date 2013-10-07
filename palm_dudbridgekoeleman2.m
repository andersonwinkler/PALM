function T = palm_dudbridgekoeleman2(G,df2,plm,c)
% Computes the combined statistic. It is calculated
% the fastest equivalent way as the original statistic.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Aug/2013
% http://brainder.org

df2 = bsxfun(@times,ones(size(G)),df2);
P = -log10(1-palm_gcdf(G,plm.tmp.rC(c),df2));
[~,tmp] = sort(G,1,'descend');
[~,tmp] = sort(tmp);
P(tmp > plm.npcparm) = 0;
P(P < -log10(plm.npcparm2)) = 0;
T = sum(P,1);
