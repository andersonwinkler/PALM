function P = palm_zaykinp(T,plm)
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
lalpha   = log10(plm.npcparm);
l1alpha  = log10(1-plm.npcparm);

P = zeros(size(lT));
for k = 1:plm.nY, % outer sum
    
    % First and second factors (logs)
    lp1 = lfac(plm.nY+1) - lfac(k+1) - lfac(plm.nY-k+1);
    lp2 = (plm.nY-k)*l1alpha;
    
    % Compute or not the inner sum depending on the
    % indicator function in Eqn.1
    Tsmall = lT <= k*lalpha;
    Tlarge = ~ Tsmall;
    
    % Third factor
    p3 = 0;
    lnum = log10(k*lalpha - lT(Tsmall));
    for j = 1:k, % inner sum
        p3 = p3 + 10.^(lT(Tsmall) + (j-1).*lnum - lfac(j));
    end
    lp3small = log10(p3);
    lp3large = k*lalpha;
    
    % Increment the probability
    P(Tsmall) = P(Tsmall) + 10.^(lp1 + lp2 + lp3small);
    P(Tlarge) = P(Tlarge) + 10.^(lp1 + lp2 + lp3large);
end
P = 1-P;