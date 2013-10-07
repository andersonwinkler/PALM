function P = palm_wilkinsonp(T,plm)
% Computes the parametric p-value for the combined statistic.
% It is calculated assuming that the statistic was computed 
% usign the function for the respective method.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Aug/2013
% http://brainder.org

% Work with logs
lfac = zeros(plm.nY+1,1);
for f = 1:plm.nY,
    lfac(f+1) = log(f) + lfac(f);
end
lalpha  = log(plm.npcparm);
l1alpha = log(1-plm.npcparm);

% Binomial expansion (truncated)
P = zeros(size(T));
for k = 1:plm.nY,
    lp1 = lfac(plm.nY+1) - lfac(k+1) - lfac(plm.nY-k+1);
    lp2 = k*lalpha;
    lp3 = (plm.nY-k)*l1alpha;
    P = P + (k>=T).*exp(lp1+lp2+lp3);
end
P = 1-P;