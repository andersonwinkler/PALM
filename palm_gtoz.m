function Z = palm_gtoz(G,df1,df2)
% Convert a G-statistic (or any of its particular cases)
% to a z-statistic (normally distributed).
% 
% Usage:
% Z = gstat2pval(G,df1,df2)
% 
% Inputs:
% G        : G statistic.
% df1, df2 : Degrees of freedom (non-infinite).
% 
% Outputs:
% Z        : Z-score
% 
% _____________________________________
% Anderson Winkler
% FMRIB / University of Oxford
% Jan/2014
% http://brainder.org

df2 = bsxfun(@times,ones(size(G)),df2);
Z = zeros(size(G));
if df1 == 1,
    
    % Deal with precision issues working on each
    % tail separately
    idx = G > 0;
    Z( idx) = -norminv(tcdf(-G( idx),df2( idx)));
    Z(~idx) =  norminv(tcdf( G(~idx),df2(~idx)));
    
else
    
    % G-vals above the upper half are treated as
    % "upper tail"; otherwise, "lower tail".
    thr = finv(.5,df1,df2);
    idx = G > thr;
    
    % Convert G-distributed variables to Beta-distributed
    % variables with parameters a=df1/2 and b=df2/2
    B = (df1.*G./df2)./(1+df1.*G./df2);
    a = df1/2;
    b = df2/2;
    
    % Convert to Z through a Beta incomplete function
    Z( idx) = -norminv(betainc(1-B( idx),b(idx),a));
    Z(~idx) =  norminv(betainc(  B(~idx),a,b(~idx)));
end
