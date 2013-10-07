function P = palm_dudbridgekoeleman2p(T,plm)
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
lT = -T;
lfac = zeros(plm.nY+1,1);
for f = 1:plm.nY,
    lfac(f+1) = log10(f) + lfac(f);
end

P = zeros(1,size(T,2));
for k = 1:plm.npcparm,
    kk = (plm.nY-k)*log(1-plm.npcparm2);
    if isnan(kk), kk = 0; end
    p1 = exp(lfac(plm.nY+1) - lfac(k+1) - lfac(plm.nY-k+1) + kk);
    p2 = awtk(lT,plm.npcparm2,k,lfac(1:k));
    P = P + p1.*p2;
end
if k < plm.nY,
    lp1 = lfac(plm.nY+1) - lfac(plm.npcparm+2) - lfac(plm.nY-plm.npcparm) + log(plm.npcparm+2); % all logs
    for v = 1:numel(lT);
        P(v) = P(v) + quad(@(t)dkint(t,lp1,lT(v),plm.nY,plm.npcparm,lfac(1:plm.npcparm)),eps,plm.npcparm2);
    end
end
P = 1-P;

function q = dkint(t,lp1,lT,K,r,lfac)
lp2 = (K-r-1).*log(1-t);
ltr = r.*log(t);
L1  = real(lp1 + lp2 + ltr);
s1  = (lT > ltr).*exp(L1);
j   = (1:r)';
lp3 = lT + (j-1)*log(r*log(t)-lT) ...
    - repmat(lfac(j),[1 numel(t)]);
L2  = real(lp1 + repmat(lp2,[r 1]) + lp3);
s2  = (lT <= ltr).*sum(exp(L2));
q   = s1 + s2;

function A = awtk(lw,t,k,lfac)
ltk = k.*log(t);
tk = real(exp(ltk));
s = (1:k)';
L = bsxfun(@plus,lw,...
    bsxfun(@minus,(s-1)*log(k*log(t)-lw),lfac(s)));
S = sum(real(exp(L)),1);
A = (lw <= ltk).*S + (lw > ltk).*tk;
