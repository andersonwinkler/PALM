function [P,maxP,lex] = palm_signflip(varargin)
% Produce a set P of nP sign-flipping matrices. P is a
% cell-array which elements are sparse arrays, each being a
% sign-flipping matrix.
% The input X can be a positive integer representing
% the dimensionality of P, or a vector of integers representing
% the blocks of observations that need to be sign-flipped
% together.
% 
% Usage:
% P = signflip(X,nP,nlx)
% 
% X   : Number of observations (a positive integer) OR
%       a vector of integers representing blocks.
% nP  : Number of sign-flipping matrices to be generated
% nlx : Force NOT using ordered sign-flips
% P   : Cell-array containing the sign-flipping matrices.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Mar/2012 (first version)
% Feb/2013 (this version)
% http://brainder.org

% Defaults
nP   = 20000;
nlx  = false;
rthr = 100;

% Accept inputs
if nargin < 1 || nargin > 3,
    error('Incorrect number of arguments');
end
if nargin >= 1, X   = varargin{1}; end
if nargin >= 2, nP  = varargin{2}; end
if nargin >= 3, nlx = varargin{3}; end

% Maximum number of permutations
if numel(X) == 1,
    n = X;
else
    [U,~,idx] = unique(X);
    n = numel(U);
end
maxP = 2^n;

% Make sure there are no more requested permutations than
% actual possible permutations
if nP > maxP || nP <= 0,
    nP = maxP;
end

% Decide whether to list all maxP sign-flips or not
if maxP/nP < rthr && ~nlx,
    lex = true;
    
    % The case P=I will be considered later, so drop the maxP case
    [~,d] = sort(rand(maxP-1,1));
    d = d(1:nP-1) - 1; % should include the 0th (all -1)
    b = de2bi(d);
else
    lex = false;
    b = randi(2,nP-1,n); % one perm less, as the P=I enters below
    b = b - 1; % put in the range 0-1
end

b(~b) = -1; % the zeros become -1
P = cell(nP,1);
if numel(X) == 1,
    P{1} = speye(n);
else
    b = b(:,idx);
    P{1} = speye(numel(idx));
end
for p = 2:nP;
    P{p} = sparse(diag(b(p-1,:)));
end
