function P = palm_algol(varargin)
% Given a sequence of integers S (stored as a vector of length n),
% returns a set P of permutations. P can be an array n by nP of
% permutation indices, or a cell-array with nP elements, each
% being a sparse permutation matrix.
% 
% Usage:
% P = palm_algol(S,flag)
% 
% S    : Sequence of n integers (row or column vector).
% flag : Logical true or false. If true, P is a cell-array
%        with nP elements, each a sparse permutation matrix.
%        If false, P is a n by nP array with permutation
%        indices. Default is false (faster and requiring
%        less memory).
% P    : Output permutation matrix.
%
% This function is an implementation of the "Algorithm L",
% published by D. Knuth in his "The Art of Computer Programming",
% Volume 4, Fascicle 2: Generating All Tuples and Permutations.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Feb/2012 (first version with Algorithm L)
% Oct/2013 (this version)
% http://brainder.org

% Compute the largest number of unique permutations
S  = varargin{1};
n  = length(S);
U  = unique(S);
nU = zeros(size(U));
for u = 1:numel(U),
    nU(u) = sum(S == U(u));
end
maxP = factorial(n)/prod(factorial(nU));

% Algorithm L
% The 1st permutation is the sorted
[a,tmp] = sort(S);       % sort, and keep indices to permute back
[~,idxback] = sort(tmp); % to 'unsort' the permutations
p = 1;                   % current permutation number
P = zeros(n,maxP);       % init the array to store all perms
P(:,p) = a;              % 1st permutation is the identity (Step L1)

% Step L2
j = n - 1;
while j > 0 && a(j) >= a(j+1),
    j = j - 1;
end

% Interrupt if j == 0
while j > 0,
    
    % Step L3
    l = n;
    while a(j) >= a(l),
        l = l - 1;
    end
    tmp  = a(j);
    a(j) = a(l);
    a(l) = tmp;
    
    % Step L4
    k = j + 1;
    l = n;
    while k < l,
        tmp  = a(k);
        a(k) = a(l);
        a(l) = tmp;
        k = k + 1;
        l = l - 1;
    end
    p = p + 1;
    P(:,p) = a;
    
    % Step L2 again
    j = n - 1;
    while j > 0 && a(j) >= a(j+1),
        j = j - 1;
    end
end

% 'Unsort' the permutations
P = P(idxback,:);

% Convert to sparse permutation matrices
if nargin > 1 && varargin{2},
    [~,tmp] = sort(P);
    [~,idx] = sort(tmp);
    P = cell(maxP,1);
    for p = 1:maxP,
        P{p} = palm_idx2perm(idx(:,p));
    end
end
