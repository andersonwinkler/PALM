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
Mn       = null(M');
e        = Y - M*(M\Y);
[Qe,~,~] = svds([Mn e/norm(e)],k);
A        = [Qm Qe];