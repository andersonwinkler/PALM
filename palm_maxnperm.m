function maxP = palm_maxnperm(Ptree)
% Computes the maximum number of possible permutations given
% a tree that specifies the depencence between the observations.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Oct/2013
% http://brainder.org

maxP = maxpnode(Ptree,1);

% ==============================================================
function np = maxpnode(Ptree,np)
% Number of permutations per node, recursive and
% incremental.
for u = 1:size(Ptree,1),
    np = np * seq2np(Ptree{u,1});
    if ~isnan(Ptree{u,2}(2)),
        np = maxpnode(Ptree{u,3},np);
    end
end

% ==============================================================
function np = seq2np(S)
% Takes a sequence of integers and computes the 
% number of possible permutations.
U   = unique(S);
nU  = numel(U);
cnt = zeros(size(U));
for u = 1:nU,
    cnt(u) = sum(S == U(u));
end
np = factorial(numel(S))/prod(factorial(cnt));
