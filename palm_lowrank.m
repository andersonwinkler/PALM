function varargout = palm_lowrank(varargin)
% Do various tasks related to lowrank matrix completion.
%
% Usage:
% U          = palm_lowrank(G)
% [Grec,mom] = palm_lowrank(G,U,nsel)
% Grec       = palm_lowrank(G,U,ysel)
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
[nJ,nV] = size(G);

if nargin == 1,
    
    % Compute the basis U. This is to be done at a specific permutation.

    % Use a global mean, as each permutation is assumed to have the same
    % overall mean as the others. The average is the best estimate.
    G = G - mean(G(:),1);
    
    % Define the basis. Save some memory by working with
    % the smallest possible square of X
    if nJ >= nV,
        [U,SS,~] = svd(G'*G);
        s   = diag(SS);
        tol = nJ * eps(max(s));
        r   = sum(s > tol);
        SSr = SS(1:r,1:r);
        U   = U(:,1:r)';
    else
        [~,SS,U] = svd(G*G');
        s   = diag(SS);
        tol = nV * eps(max(s));
        r   = sum(s > tol);
        SSr = SS(1:r,1:r);
        U   = U(:,1:r)'*G;
    end
    
    % Outputs
    varargout{1} = diag(diag(SSr).^-.5)*U;
    
elseif nargin == 3,
    
    if size(G,1) > 1,
        
        % With a basis known, reconstruct G, find the residuals,
        % and compute the moments.
        
        % Inputs
        U    = varargin{2}; % basis
        nsel = varargin{3}; % number of voxels to use
        
        % Reconstruct the data in the new basis. Use just a subsample, so that
        % the residuals will also include this source of variability (i.e., the
        % filling of missing data, not just low-rank basis).
        varargout{1} = zeros(size(G));
        for p = 1:nJ,
            idx = randperm(nV);
            idx = idx(1:nsel);
            varargout{1}(p,:) = G(p,idx)*pinv(U(:,idx))*U;
        end
        
    else
        
        % With a basis known, reconstruct a single column of G with just a few
        % known entries, and add residuals that follow a similar distribution.
        
        % Inputs
        U    = varargin{2}; % basis
        ysel = varargin{3}; % indices of selected tests
        
        % Reconstruct a single permutation with lots of missing values.
        varargout{1} = G*pinv(U(:,ysel))*U;
    end
end
