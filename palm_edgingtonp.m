function P = palm_edgingtonp(T,plm)
% Computes the parametric p-value for the combined statistic.
% It is calculated assuming that the statistic was computed 
% usign the function for the respective method.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Aug/2013
% http://brainder.org

lfac = zeros(plm.nY+1,1);
for f = 1:plm.nY,
    lfac(f+1) = log(f) + lfac(f);
end

T = 10.^(-T);
fT   = floor(T);
mxfT = max(fT(:));
P = zeros(size(T));
for j = 0:mxfT,
    p1  = (-1)^j;
    lp2 = - lfac(j+1) - lfac(plm.nY-j+1);
    lp3 = plm.nY*log(T-j);
    P = P + (j<=fT).*p1.*exp(lp2+lp3);
end
P = 1-P;