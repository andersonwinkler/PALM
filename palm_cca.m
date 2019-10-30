function [cc,A,B,U,V,lW] = palm_cca(varargin)
% Do CCA via QR & SVD.
%
% Usage:
% [cc,A,B,U,V] = palm_cca(X,Y,df)
%
% Inputs:
% X and Y : Arrays with N rows. The number of columns
%           in each may differ. The ranks of X and Y
%           are not checked for speed.
% df      : Optional. Degrees of freedom (N - rank(Z)),
%           where Z is the matrix with nuisance variables
%           already regressed out from both X and Y.
%           Since both are supposed to be mean-centered,
%           df should be no larger than N-1. For speed,
%           it is not tested internally, though, so enter
%           this wisely.
%           If omited, or if equal to 0, assumes that the
%           data have NOT been residualised or mean-centered.
%           In this case, the data will be mean-centered.
%
% Outputs:
% cc       : Canonical correlations.
% A and B  : Canonical coefficients.
% U and V  : Canonical variables.
% W        : Log of Wilks' lambda.
%
% U=X*A and V=Y*B. such that each pair of columns in U and V
% are maximally correlated, under the orthonality constraint.
%
% Note: There is a faster CCA function inside palm_core.m,
%       which is used in the permutation test. This one gives
%       more outputs and is used for functionality other than 
%       what is needed for CCA permutation inference.
% 
% Based on the algorithm proposed by:
% * Bjorck A, Golub GH. Numerical methods for
%   computing angles between linear subspaces.
%   Math Comput. 1973;27(123):579-579.
%
% _____________________________________
% Anderson M. Winkler
% National Institutes of Health
% Nov/2018 (first version)
% Oct/2019 (this version)
% http://brainder.org

% Parse inputs:
X = varargin{1};
Y = varargin{2};
if nargin > 2
    df = varargin{3};
else
    df = 0;
end

% Mean center if needed:
if df == 0
    X  = X - mean(X);
    Y  = Y - mean(Y);
    df = size(X,1) - 1;
end

% Do the CCA proper:
[Qx,Rx,iX] = qr(X,0);
[Qy,Ry,iY] = qr(Y,0);
%Px = palm_idx2perm(iX');
%Py = palm_idx2perm(iY');
k  = min(rank(X),rank(Y));
[L,D,M] = svds(Qx'*Qy,k);

% Canonical correlations:
cc = min(max(diag(D(:,1:k))',0),1);

% Canonical coefficients (loadings):
if nargout > 1
    A  = Rx\L(:,1:k)*sqrt(df);
    B  = Ry\M(:,1:k)*sqrt(df);
    A(iX,:) = A;
    B(iY,:) = B;
%    A  = Rx*Px\L(:,1:k)*sqrt(df);
%    B  = Ry*Py\M(:,1:k)*sqrt(df);
end

% Canonical variables:
if nargout > 3
    U  = X*A;
    V  = Y*B;
end

% Wilks' lambda:
if nargout > 5
    lW = -fliplr(cumsum(fliplr(log(1-cc.^2))));
    %W  = 10.^-lW;
end