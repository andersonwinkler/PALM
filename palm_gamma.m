function pvals = palm_gamma(G,mu,sigsq,gamm1,rev)
% Return the p-values for a Gamma distribution, parameterised by
% its first three moments.
%
% pvals = palm_gamma(G,mu,s2,gamm1,rev)
% 
% - G     : Statistics for which p-values are to be computed.
% - mu    : Distribution mean.
% - sigsq : Distribution standard deviation.
% - gamm1 : Distribution skewness.
% - rev   : Use if lower values of the statistic are evidence in
%           favour of the alternative.
% - pvals : p-values.
% 
% References:
% * Mielke PW, Berry KJ, Brier GW. Application of Multi-Response
%   Permutation Procedures for Examining Seasonal Changes in
%   Monthly Mean Sea-Level Pressure Patterns. Mon Weather Rev.
%   1981;109(1):120-126.
% * Minas C, Montana G. Distance-based analysis of variance:
%   Approximate inference. Stat Anal Data Min. 2014;7(6):450-470.
% 
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% May/2015
% http://brainder.org

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% PALM -- Permutation Analysis of Linear Models
% Copyright (C) 2015 Anderson M. Winkler
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Note that there are no argument checking for speed, but
% sizes of all inputs need to be the same, or the moments need to
% be all scalars.

if gamm1 == 0,
    
    % If not skewed, use a normal approximation.
    G     = (G - mu)./sigsq.^.5;
    pvals = erfc(G/sqrt(2))/2;
    
else
    
    % For the negative skewness, flip the sign of G and mu. It's still
    % possible that G will be negative even with the skewness flipped to
    % the positive side. These cases are treated later.
    sg    = sign(gamm1);
    G     = sg.*G;
    mu    = sg.*mu;
    gamm1 = abs(gamm1);
    
    % Standardise G, so that all becomes a function of the skewness.
    G     = (G - mu)./sigsq.^.5;
    
    % Gamma distribution parameters (Minas & Montana, 2014).
    kpar  = 4/gamm1.^2;
    tpar  = gamm1/2;
    cpar  = -2/gamm1;
     
    % Actual p-value. If there are negatives here, the probability can
    % have an imaginary part, and the real is above unity. Taking the
    % min deals with the issue.
    if rev,
        tail = 'lower';
    else
        tail = 'upper';
    end
    pvals = min(1,gammainc((G-cpar)./tpar,kpar,tail));
end
