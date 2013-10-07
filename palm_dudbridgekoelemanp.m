function P = palm_dudbridgekoelemanp(T,plm)
% Computes the parametric p-value for the combined statistic.
% It is calculated assuming that the statistic was computed 
% usign the function for the respective method.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Aug/2013
% http://brainder.org

lT = -T;
lfac = zeros(plm.nY+1,1);
for f = 1:plm.nY,
    lfac(f+1) = log10(f) + lfac(f);
end

P = zeros(size(lT));
lp1 = lfac(plm.nY+1) - lfac(plm.npcparm+2) - lfac(plm.nY-plm.npcparm) + log10(plm.npcparm+2); % all logs
for v = 1:numel(lT);
    P(v) = quad(@(t)dkint(t,lp1,lT(v),plm.nY,plm.npcparm,lfac(1:plm.npcparm)),eps,1);
end
P = 1-P;

function q = dkint(t,lp1,lT,K,r,lfac)
lp2 = (K-r-1).*log10(1-t);
ltr = r.*log10(t);
L1  = real(lp1 + lp2 + ltr);
s1  = (lT > ltr).*10.^(L1);
j   = (1:r)';
lp3 = lT + (j-1)*log10(r*log10(t)-lT) ...
    - repmat(lfac(j),[1 numel(t)]);
L2  = real(lp1 + repmat(lp2,[r 1]) + lp3);
s2  = (lT <= ltr).*sum(10.^(L2));
q   = s1 + s2;