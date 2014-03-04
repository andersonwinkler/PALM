function pvals = palm_gpval(G,df1,df2)
% Wrapper to convert a pivotal statistic computed with
% 'pivotal.m' (or simplifications) to a parametric p-value.
%
% Note that this doesn't try to deal with precision
% issues when the p-value is close to 1. Use
% instead palm_gcdf to get these values right, i.e., the
% p-vals that *don't usually matter* are then close to zero.
%
% Usage:
% pvals = palm_gpval(G,df1,df2)
%
% Inputs:
% G        : G or Z statistic.
% df1, df2 : Degrees of freedom (non infinite).
%            df1 must be a scalar.
%            For z, use df1 = 0.
%            For Chi2, use df1 = -1, and df2 as the df.
%
% Outputs:
% pvals    : Parametric p-values based on a
%            t, F, z or Chi2 distribution.
%
% _____________________________________
% Anderson Winkler
% FMRIB / University of Oxford
% Aug/2013
% http://brainder.org

% Make sure the sizes match
df2 = bsxfun(@times,ones(size(G)),df2);

if df1 > 1,
    
    % G or F, via conversion to Beta
    B = (df1.*G./df2)./(1+df1.*G./df2);
    pvals = betainc(1-B,df2/2,df1/2);
    
elseif df1 == 1,
    
    % Student's t, Aspin-Welch v
    pvals = tcdf(-G,df2);
    
elseif df1 == 0,
    
    % Normal distribution
    pvals = normcdf(-G);
    
elseif df1 < 0,
    
    % Chi^2, via upper Gamma incomplete for precision and speed
    pvals = gammainc(G/2,df2/2,'upper');
end
