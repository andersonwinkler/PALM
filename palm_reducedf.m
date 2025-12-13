function A = palm_reducedf(Y,M,k,po)
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
% Y  : Dependent data
% M  : Design matrix
% k  : Degrees of freedom of the resulting model
% po : Ensure Y is the positive octant? Default is false.
%
% Outputs:
% A : Transformation matrix
%
% _____________________________________
% Anderson M. Winkler
% UTRGV
% Dec/2025
% http://brainder.org

% Parse arguments
narginchk(3,4);
if nargin == 3
    po = false; % positive octant for Y?
end

% Space spanned by the model
[Qm,~,~] = svd(M,'econ');

% Space not spanned by the model
% For typical designs, we want the normed residuals included in the
% null space of the design, and with equal weight, so we want to
% use e/norm(e) in the second SVD. However, if the residuals are zero
% (within the eps), then norming becomes a problem (for example, it
% causes ssq(A'*Y) > ssq(Y)). So here we only normalize and include 
% the residuals if they are non-zero:
Mn = null(Qm');
e  = Y - Qm*(Qm\Y);
if e'*e > eps(10)
    e = e/norm(e);
else
    e = [];
end
[Qe,~,~] = svds([Mn e],k);

% A is simply a matrix that spans both the design and the residual space
A = [Qm Qe];

% Ensure the dependent data Y is in the all positive octant
if po
    y = A'*Y;
    for d = 1:3
        if y(d) < 0
            A(:,d) = -A(:,d);
        end
    end
end