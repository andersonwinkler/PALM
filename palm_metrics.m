function varargout = palm_metrics(varargin)
% Compute some permutation metrics:
% - For permutation trees, return the ratios of entropies
% - For sets of permutations, return the average Hamming distance.
% 
% Usage:
% [lW,lW0,rl,lr,C] = palm_metrics(Ptree,X,stype)
% [Hamm,HammX]     = palm_metrics(Pset,X)
% 
% Inputs:
% - Ptree : Permutation tree.
% - X     : Design matrix (only the EVs of interest for the Freedman-Lane
%           and most methods, or the full matrix for ter Braak).
%           Note that the metrics are only meaningful if X is the same
%           used when Ptree was created originally.
% - stype : Shuffling type. It can be 'perms', 'flips' or 'both'.
% - Pset  : Set of shufflings (permutations or sign-flips).
% 
% Outputs
% - lW    : Log of the max number of permutations with the restrictions
%           imposed by the tree and the original design used to create the
%           tree.
% - lW0   : Log of the max number of permutations without the restrictions
%           imposed by the tree, but with the restrictions imposed by the
%           input design X.
% - rl    : Ratio of the logs of the number of possible permutations
%           in the restricted over the unrestricted, subtracted from 1,
%           i.e., rl = 1 - log(W)/log(W0). This is the "anisotropy".
%           Values close to 1 indicate that there is very strong structure
%           within the data. Values close to 0 indicate weak or no structure,
%           such that the data can be shuffled freely or almost freely.
% - lr    : Negative log of the ratio, i.e. lr = -log(W/W0). This is
%           self-explanatory :-)
%           For rl and lr to be meaningful, the imput X must be the same
%           used originally to create the permutation tree.
% - C     : Huberman & Hogg complexity (C) of a given tree.
%           For this to give exactly the same result as in the original
%           paper, such that it measures the relationships in the tree
%           itself, rather than the actual values found in X, the input
%           Ptree must have been constructed with X = ones(N,1) (or any
%           other constant. However, C doesn't depend on the X that is
%           input (i.e., it's not an argument needed to compute C, but
%           it's implicitly taken into account through the tree).
% - Hamm  : Average Hamming distance across the given permutation set,
%           i.e., it's the average change that a permutation cause on
%           the original indices.
% - HammX : Same as Hamm, but consider the possibility of repeated
%           elements in X. If X isn't supplied, or if X has no ties,
%           or if X is the same used originally to create the permutation
%           set, Hamm and HammX are the same.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Feb/2014
% http://brainder.org

% Take args and decide what to do
X = [];
if iscell(varargin{1}),
    dowhat = 'entropy';
    Ptree  = varargin{1};
    N      = numel(palm_permtree(Ptree,1,false,true));
    if nargin == 1,
        X     = (1:N)';
        stype = 'perms';
    elseif nargin == 2,
        X     = varargin{2};
        stype = 'perms';
    elseif nargin == 3,
        X     = varargin{2};
        stype = varargin{3};
    end
else
    dowhat = 'hamming';
    Pset   = varargin{1};
    N      = size(Pset,1);
    if nargin == 1,
        X  = (1:N)';
    elseif nargin == 2,
        X  = varargin{2};
    end
end
if isempty(X),
    X = (1:N)';
end

switch dowhat,
    case 'entropy',
        
        % Normalised entropy or anisotropy. This is computed
        % if the first input is a cell array (Ptree).
        
        % Number of permutations (log) given the data structure.
        % It is assumed that Ptree was constructed using the
        % same X as input.
        lW = palm_maxshuf(Ptree,stype,true);
        varargout{1} = lW;
        
        % Number of permutations (log) if there were no data structure:
        if nargout > 1,
            lfac = palm_factorial(N);
            [~,~,S] = unique(X,'rows');
            U   = unique(S);
            nU  = numel(U);
            if strcmpi(stype,'perms') || strcmpi(stype,'both'),
                cnt = zeros(nU,1);
                for u = 1:nU,
                    cnt(u) = sum(S == U(u));
                end
                plW0 = lfac(N+1) - sum(lfac(cnt+1));
            else
                plW0 = 0;
            end
            if strcmpi(stype,'flips') || strcmpi(stype,'both'),
                cnt = zeros(nU,1);
                for u = 1:nU,
                    cnt(u) = sum(S == U(u));
                end
                slW0 = nU*log(2);
            else
                slW0 = 0;
            end
            lW0 = plW0 + slW0;
            varargout{2} = lW0;
        end
        
        % Anisotropy of the data structure, i.e.,
        % ratio of the logs of the entropies.
        if nargout > 2,
            varargout{3} = 1 - lW/lW0;
        end
        
        % If the user wants the log of the ratios: -log(W/W0).
        if nargout > 3,
            varargout{4} = lW0 - lW;
        end
        
        % If the user wants, output also the Huberman & Hogg complexity,
        % which is computed recursively below
        if nargout > 4,
            varargout{5} = hhcomplexity(Ptree,1) - 1;
        end
        
    case 'hamming',
        
        % Average change per permutation, i.e., average
        % Hamming distance.
        varargout{1} = mean(sum(bsxfun(@ne,Pset(:,1),Pset),1),2);
        
        % Now take ties in X into account:
        XP = X(Pset);
        varargout{2} = mean(sum(bsxfun(@ne,XP(:,1),XP),1),2);
end

% ==============================================================
function D = hhcomplexity(Ptree,D)
% Computes recursively the Huberman & Hogg complexity.
% For the 1st iteration, D = 1.
for u = 1:size(Ptree,1),
    if isnan(Ptree{u,1}(1)),
        k = size(Ptree{u,3},1);
    else
        k = numel(unique(Ptree{u,1}(:,1)));
    end
    D = D * (2^k - 1);
    if size(Ptree{u,3},2) > 1,
        D = hhcomplexity(Ptree{u,3},D);
    end
end
