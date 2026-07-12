function X = deflation_regression(X, s, dimred)

% Performs deflation by subtracting the estimated source contribution to the observations as:
% 
%    X' = X - h*s
% 
% The source direction h is estimated via the least squares solution to the linear regression
% problem:
%
%    h_opt = arg min_h ||X - h*s||^2 = X*s'/(s*s').
%
%
% SYNTAX: Xn = deflation_regression(X, s, dimred);
%
%
% OUTPUT:
%         Xn     : observed data after subtraction (one signal per row, one sample per column).
%
%
% INPUTS:
%         X      : observed data (one signal per row, one sample per column)
%
%         s      : estimated source (row vector with one sample per column)
%
%         dimred : perform dimensionality reduction if parameter different from zero.
%
%
% Please, report any bugs, comments or suggestions to <a href = "mailto:zarzoso@i3s.unice.fr">zarzoso(a)i3s.unice.fr</a>.
%
%
% HISTORY:
% 
%    <Please add modification date here>: - <please add modification details here>
%
% -- 2014/11/21: Version 3 release ----------------------------------------------------------------
%
% -- 2010/02/16: Version 2 release ----------------------------------------------------------------
%
% -- 2008/03/31: Version 1 release ----------------------------------------------------------------
%
%    2008/03/26: - added this help
%
%    2008/03/18: - if required, perform dimensionality reduction via QR decomposition
%   
%    2008/03/13: - created by Vicente Zarzoso
%                  (I3S Laboratory, University of Nice Sophia Antipolis, CNRS, France).
 


s2 = s*s';          % extracted source power times sample size

if abs(s2) > eps    % don't perform subtraction if estimated component is null

h = X*s'/s2;        % source direction estimated via least squares 


if dimred

% with dimensionality reduction (old version)  *************************
    n = length(h);
    Q = [h, eye(n, n-1)];          
    [Q, R] = qr(Q);                 
    Q = Q(:,2:n);   % orthonormal basis of orhogonal subspace of h
    X = Q'*X;       % remaining contribution with dimensionality reduction   

else
    
% without dimensionality reduction
X = X - h*s;        

% if dimred   %%% an alternative version?  *** TO BE TESTED ***
%     [n, T] = size(X);
%     [V, S, U] = svd(X', 0); %'economy' SVD
%     hU = abs(h'*U)
%     diagS = diag(S)
%     X = sqrt(T)*V(:, 1:n-1)';
%     pause
% end % if dimred

end % if dimred

end % if s2 > eps
