function T = palm_darlingtonhayes(G,df2,plm,c)
% Computes the combined statistic. It is calculated
% the fastest equivalent way as the original statistic.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Aug/2013
% http://brainder.org

df2     = bsxfun(@times,ones(size(G)),df2);
[~,tmp] = sort(G,1,'descend');
[~,tmp] = sort(tmp);
idx     = tmp <= plm.npcparm;
G       = reshape(G(idx),  [plm.npcparm size(G,2)]);
df2     = reshape(df2(idx),[plm.npcparm size(df2,2)]);
P       = palm_gcdf(G,plm.tmp.rC(c),df2);
Z       = norminv(P);
T       = mean(Z,1);