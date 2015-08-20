function varargout = palm_lowrank(varargin)
% Do various tasks related to lowrank matrix completion.
% 
% Usage:
% U          = palm_lowrank(G)
% [Grec,mom] = palm_lowrank(G,U,nsel)
% Grec       = palm_lowrank(G,U,ysel,mom)
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% May/2015
% http://brainder.org

% Note that the basis U is on rowspace. This is so to minimise the number
% of matrix transposes in the code, which can be slow for large matrices.

% Common input for all cases below
G       = varargin{1};
[nP,nV] = size(G);

if nargin == 1,
    
    % Compute the basis U. This is to be done at a specific permutation.
    
    % Effective rank
    r = nP/2; % nP should always be a multiple of r.

    % Use a global mean, as each permutation is assumed to have the same
    % overall mean as the others. The average is the best estimate.
    G = G - mean(G(:),1);
    
    % Define the basis. Save some memory by working with
    % the smallest possible square of X
    if nP >= nV,
        [U,SS,~] = svd(G'*G);
        SSr = SS(1:r,1:r);
        U   = U(:,1:r)';
    else
        [~,SS,U] = svd(G*G');
        SSr = SS(1:r,1:r);
        U   = U(:,1:r)'*G;
    end
    
    % Outputs
    varargout{1} = diag(diag(SSr).^-.5)*U;
    
elseif nargin == 3,
    
    % With a basis known, reconstruct G, find the residuals,
    % and compute the moments.
    
    % Inputs
    U    = varargin{2}; % basis
    nsel = varargin{3}; % number of voxels to use

    % Reconstruct the data in the new basis. Use just a subsample, so that
    % the residuals will also include this source of variability (i.e., the
    % filling of missing data, not just low-rank basis).
    Grec = zeros(size(G));
    for p = 1:nP,
        idx = randperm(nV);
        idx = idx(1:nsel);
        Grec(p,:) = G(p,idx)*pinv(U(:,idx))*U;
    end
    
    % Compute the first three moments of the distribution of the residuals.
    [mu,s2,gamm1] = palm_moments(G(:) - Grec(:));
    
    % Add some random noise to the reconstructed data.
    Grec  = Grec + gammarnd([mu,s2,gamm1],size(Grec));
    
    % Outputs
    varargout{1} = Grec;
    varargout{2} = [mu,s2,gamm1];
    
elseif nargin == 4,
    
    % With a basis known, reconstruct a single column of G with just a few
    % known entries, and add residuals that follow a similar distribution.
    
    % Inputs
    U    = varargin{2}; % basis
    ysel = varargin{3}; % indices of selected tests 
    mom  = varargin{4}; % previously computed moments
    
    % Reconstruct a single permutation, with lots of missing values,
    % to the same basis, and add the random residuals.
    Grec  = G*pinv(U(:,ysel))*U;
    Grec  = Grec + gammarnd(mom,size(Grec));
    
    % Outputs
    varargout{1} = Grec;
end

% ==============================================================
function g = gammarnd(mom,siz)
% Generate random numbers following a Gamma distribution
% parameterised according to the first 3 moments.
kpar = 4/mom(3).^2;
tpar = sqrt(mom(2)/kpar);
sgn  = sign(mom(3));
g    = tpar*(sgn*(palm_randg(kpar,siz)-kpar)) + mom(1);
