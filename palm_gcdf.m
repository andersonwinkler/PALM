function pval = palm_gcdf(G,df1,df2)
% Wrapper to convert a pivotal statistic computed with
% 'pivotal.m' (or simplifications) to a parametric p-value.
% 
% Usage:
% pval = gstat2pval(G,df1,df2)
% 
% Inputs:
% G        : Pivotal statistic.
% df1, df2 : Degrees of freedom.
% 
% Outputs:
% pval     : Parametric p-value (1-p), based on a
%            t or F distribution.
% 
% _____________________________________
% Anderson Winkler and Tom Nichols
% FMRIB / University of Oxford
% Aug/2013
% http://brainder.org

df2 = bsxfun(@times,ones(size(G)),df2);
if df1 == 1,
    pval = tcdf(G,df2);
else
    pval = fcdf(G,df1,df2);
end
