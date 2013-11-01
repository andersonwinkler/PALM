function a = palm_nextperm(a)
% Given a sequence of integers "a", returns the next lexicographic
% permutation of this sequence. If "a" is already the last possible
% permutation, returns a vector of zeros of size(a).
% Note that to shuffle vectors, they must be supplied as
% column vectors (N by 1).
% 
% Usage:
% a1 = palm_nextperm(a)
% 
% a  : 2D array to be shuffled. Only the 1st column is
%      considered for the permutations. The rows as a whole
%      are shuffled together.
% a1 : Permuted sequence of values that corresponds to the
%      next lexicographic permutation. If a is already the last
%      possible permutation, a1 = false.
% 
% This function is an implementation of the "Algorithm L",
% published by D. Knuth in his "The Art of Computer Programming",
% Volume 4, Fascicle 2: Generating All Tuples and Permutations.
% See also palm_algol to produce all possible permutations for
% a given sequence in a single function.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Feb/2012 (first version)
% Oct/2013 (this version)
% http://brainder.org

% Algorithm L
% Step L2
n = size(a,1);
j = n - 1;
while j > 0 && a(j,1) >= a(j+1,1),
    j = j - 1;
end

% If this isn't yet the last permutation, brings up the next one.
if j > 0,
    
    % Step L3
    l = n;
    while a(j,1) >= a(l,1),
        l = l - 1;
    end
    tmp  = a(j,:);
    a(j,:) = a(l,:);
    a(l,:) = tmp;
    
    % Step L4
    k = j + 1;
    l = n;
    while k < l,
        tmp  = a(k,:);
        a(k,:) = a(l,:);
        a(l,:) = tmp;
        k = k + 1;
        l = l - 1;
    end
else
    % If the input is the last permutation, then there is no next.
    % Return then a boolean "false" that can be tested outside.
    a = false;
end