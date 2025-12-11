function A = palm_reducedf(Y,M,k)
% Given Y and X of a general linear model:
% 
%       Y = M*b + e
% 
% produces a matrix A such that the model:
% 
%    A'*Y = A'*M*b + A'*e
% 
% has the following properties:
% 
% [1] Same least squares solutions b:
%         b = M\Y = (A'*M)\(A'*Y)
% [2] Same sums of squares of the residuals:
%        (e'*e) = (A'*e)'*(A'*e)
% [3] Same sums of squares of data:
%        (Y'*Y) = (A'*Y)'*(A'*Y)
% [4] Same sum of squares and of products of design:
%        (M'*M) = (A'*M)'*(A'*M)
% [5] Degrees of freedom k as specified,
%        1 <= k <= N-rank(M)
% 
% If k = 0, [1] and [4] also hold.
% 
% Usage:
% A = palm_reducedf(Y,M,k)
% 
% Inputs:
% Y : Dependent data
% M : Design matrix
% k : Degrees of freedom of the resulting model
%
% Outputs:
% A : Transformation matrix
% 
% _____________________________________
% Anderson M. Winkler
% UTRGV
% Dec/2025
% http://brainder.org

[Qm,~,~] = svd(M,'econ');
Mn       = null(Qm');
e        = Y - Qm*(Qm\Y);
% For typical designs, we want e/norm(e) so that e enters in
% the 2nd svd with same norm as the columns of Mn. However, if
% the residuals are zero (within the eps), then norming
% becomes a problem (for example, ssq(A'*Y) > ssq(Y)).
if e'*e > eps(10)
    e    = e/norm(e);
else
    e    = [];
end
[Qe,~,~] = svds([Mn e],k); % see note below
A        = [Qm Qe];