function [scores,coeffs,evals,extras] = palm_pca(varargin)
% Return the first p eigenvectors and eigenvalues.
%
% Usage:
% [scores,coeffs,evals,extras] = palm_pca(X,p,Z,nrm)
%
% Inputs:
% X     : 2D array
% p     : Number of eigenvectors and eigenvalues to be
%         returned. Use 0 for automatic. Default is 0.
% Z     : Nuisance variables. By default, Z = ones(N,1),
%         i.e., an intercept for mean-centering.
%         To omit removal of any nuisance, and thus not
%         even mean-center, use Z = [].
% nrm   : Boolean indicating whether variance of the
%         residualized data should be scaled to unit variance.
% Sel   : Selection matrix
% 
% Outputs:
% scores : Principal components.
% coeffs : Principal coefficients (PCA loadings).
% evals  : Eigenvalues.
% extras : A struct containing a bunch of self-explanatory
%          useful outputs, including variance explained
%          by each principal component and the semi-orthogonal
%          matrix used for residualisation.
%
% Notes:
% * If the data were residualised before calling this function
%   and Z = [], then you need to rescale the eigenvalues
%   as Evals*N/(N-Nz).
%
% * PCA is applied such that:
%      scores*coeffs' = Qz'*X*pinv(std_res)
%   Qz*scores*coeffs' = (X-Z*betas)*pinv(std_res)
%
% * Thus, the original data can be reconstructed as:
%   X = Qz*scores*coeff'*std_res + Z*betas
%
% * And the scores can be computed as:
%   scores = Qz'*X*pinv(std_res)*pinv(coeffs')
%
% This function was called "pca", but after Matlab created its
% own similarly-named "pca" function, it was renamed to "epca".
% The current version is a further improvement, hence "npca".
%
% _____________________________________
% Anderson M. Winkler
% Yale University / Institute of Living
% Nov/2010 (first version)
% Mar/2026 (this version)
% http://brainder.org

% Accept and check arguments
p   = 0;
Z   = 1;
nrm = false;
Sel = [];
if nargin < 1 || nargin > 5
    error('Incorrect number of arguments');
end
if nargin >= 1
    X   = varargin{1};
end
if nargin >= 2
    p   = varargin{2};
end
if nargin >= 3
    Z   = varargin{3};
end
if nargin >= 4
    nrm = varargin{4};
end
if nargin == 5
    Sel = varargin{5};
end
[nR,nC] = size(X);
if p > max(nR,nC)
    error('Cannot extract more eigenvalues than rows or columns.');
end

% Remove nuisance variables
if isscalar(Z)
    Z = ones(nR,1);
end
if isempty(Z)
    Q = eye(nR);
    b = [];
else
    Q = palm_semiortho(Z,Sel);
    b = Z\X;
end
X = Q'*X;
if nrm
    sX = std(X,0,1);
    X  = bsxfun(@rdivide,X,sX);
else
    sX = ones(1,size(X,2));
end

% Save some memory by working with the
% smallest possible square of X
if nR >= nC
    [~,SS,V] = svd(X'*X,0);
    rnk = simplerank(X,SS);
    if p == 0
        p = wachter(X,diag(SS),rnk);
    else
        p = min(p,rnk);
    end
    p = max(p,1);
    Vp  = V(:,1:p);
    SSp = SS(1:p,1:p);
    Up  = X*Vp;
else
    [U,SS,~] = svd(X*X',0);
    rnk = simplerank(X,SS);
    if p == 0
        p = wachter(X,diag(SS),rnk);
    else
        p = min(p,rnk);
    end
    p = max(p,1);
    Up  = U(:,1:p);
    SSp = SS(1:p,1:p);
    Vp  = X'*Up;
end

% Pick one sign
s = diag(sign(Up(1,:)));

% Eigenvectors and eigenvalues
scores = Up*s;
coeffs = Vp*s;
if isscalar(Z)
    df = nR - 1;
else
    df = size(X,1);
end
evals = diag(SSp)./df;

% Some extra outputs, not normally needed
if nargout == 4
    
    % Scaled eigenvectors
    extras.scores_scaled = X*Vp;
    extras.coeffs_scaled = Up'*X;
    
    % Unit norm eigenvectors
    if nR >= nC
        extras.scores_unit = scores/sqrt(SSp);
        extras.coeffs_unit = coeffs;
    else
        extras.scores_unit = scores;
        extras.coeffs_unit = (coeffs/sqrt(SSp))';
    end
    
    % Recovered data using the p eigenvectors
    extras.recovered_with_scores = Up*extras.coeffs_scaled;
    extras.recovered_with_coeffs = extras.scores_scaled*Vp';
    
    % Variance explained
    S = diag(SSp)./df;
    extras.variance_explained = S./sum(S);
    
    % Semiorthogonal matrix used for residualisation
    extras.semi_ortho = Q;
    
    % GLM outputs from the residualization of X by Z.
    % PCA is applied to (X-Z*betas)*pinv(std_res).
    % Thus, the original data can be reconstructed as:
    % X = scores*coeff'*std_res - Z*betas
    extras.glm_betas   = b;
    extras.glm_std_res = sX;
end

% =============================================================
function rnk = simplerank(X,SS)
dSS = diag(SS);
tol = max(size(X)) * eps(max(abs(dSS)));
rnk = sum(sqrt(dSS) > tol);

% =============================================================
function p = wachter(X,SS,rnk)
siz   = size(X);
y     = siz(1)/siz(2);
if y > 1, y = 1./y; end
Obs   = SS(1:rnk)'/rnk;
Pexp  = ((1:rnk)-.5)/rnk;
Pobs  = zeros(1,rnk);
sigsq = var(X(:),0);
for k = 1:rnk
    Pobs(k) = palm_mpcdf(Obs(k).*y,y,sigsq,true);
end
P_ratio = -log10(Pobs./Pexp);
p       = sum(P_ratio >= .5);
