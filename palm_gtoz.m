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
% If df2 = NaN and df1 = 1, G is treated as Pearson's r.
% If df2 = NaN and df1 > 1, G is treated as R^2.
% 
% _____________________________________
% Anderson Winkler
% FMRIB / University of Oxford
% Jan/2014
% http://brainder.org

% If df2 isn't supplied, use a Fisher r-to-z trasnformation
if isnan(df2),
    
    % If rank(C) > 1, i.e., df1 > 1, this is R^2.
    if df1 > 1,
        G = 2*G - 1;        
    end
    
    % Apply Fisher's r-to-z transform
    Z = atanh(G);
    
else
    siz = size(G);
    Z   = zeros(siz);
    df2 = bsxfun(@times,ones(siz),df2);
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
end