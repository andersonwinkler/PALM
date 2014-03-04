function gcdf = palm_gcdf(G,df1,df2)
% Wrapper to convert a pivotal statistic computed with
% 'pivotal.m' (or simplifications) to a parametric p-value.
% The output is 1-p, i.e. the CDF.
% 
% Note that this doesn't try to deal with precision
% issues when the value of the cdf is close to 1. Use
% instead palm_gpval to get these values right, i.e., the
% p-vals *that usually matter* are then close to zero.
% 
% Usage:
% cdf = palm_gcdf(G,df1,df2)
% 
% Inputs:
% G        : G or Z statistic.
% df1, df2 : Degrees of freedom (non infinite).
%            df1 must be a scalar
%            For z, use df1 = 0.
%            For Chi2, use df1 = -1, and df2 as the df.
% 
% Outputs:
% cdf      : Parametric cdf (1-p), based on a
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
    
    % G or F
    B = (df1.*G./df2)./(1+df1.*G./df2);
    gcdf = betainc(B,df1/2,df2/2);

elseif df1 == 1,
    
    % Student's t, Aspin's v
    gcdf = tcdf(G,df2);
    
elseif df1 == 0,
    
    % Normal distribution
    gcdf = normcdf(G);
    
elseif df1 < 0,
    
    % Chi^2, via lower Gamma incomplete for precision and speed
    gcdf = gammainc(G/2,df2/2);
    
end
