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

% Note that for speed, there's no argument checking,
% and some lines are repeated inside the conditions.

if df1 > 1,
    
    % G or F
    df2 = bsxfun(@times,ones(size(G)),df2);
    B = (df1.*G./df2)./(1+df1.*G./df2);
    gcdf = betainc(B,df1/2,df2/2);

elseif df1 == 1,
    
    % Student's t, Aspin's v
    df2 = bsxfun(@times,ones(size(G)),df2);
    ic = df2 == 1;
    in = df2 > 1e7;
    ig = ~(ic|in);
    gcdf = zeros(size(G));
    if any(ig(:)),
        gcdf(ig) = betainc(1./(1+G(ig).^2./df2(ig)),df2(ig)/2,.5)/2;
    end
    ig = G > 0 & ig;
    gcdf(ig) = 1 - gcdf(ig);
    if any(ic(:)),
        gcdf(ic) = .5 + atan(G(ic))/pi;
    end
    if any(in(:)),
        gcdf(ic) = palm_gcdf(G(in),0);
    end

elseif df1 == 0,
    
    % Normal distribution
    gcdf = erfc(-G/sqrt(2))/2;
    
elseif df1 < 0,
    
    % Chi^2, via lower Gamma incomplete for precision and speed
    df2 = bsxfun(@times,ones(size(G)),df2);
    gcdf = gammainc(G/2,df2/2);
    
end
