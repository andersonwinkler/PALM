function [P,maxPall,lex] = palm_plexi(varargin)
% Given a sequence 'S' of 'n' values, some of them that
% may be repeated, output a number 'nP' of permutations.
%
% If the maximum number of possible permutations is smaller than
% 100x the requested number of permutations, then run a
% lexicographic permutation algorithm, which ensures permutation
% without replacement and takes care of the repetitions.
% If the number of possible permutations is larger than 100x nP,
% the requested number of permutations, then run a simple permutation
% algorithm without replacement. Note that in this case, although
% the permutation algorithm is without replacement, if there are
% repeated values in S, there may still be repeated permutations.
%
% Usage:
% P = plexi(S,nP,blk,pB,nlx)
%
% Inputs:
% S   : Sequence of values.
% nP  : Desired number of permutations (default: 20000).
% blk : Definition of the exchangeability blocks.
% pB  : Whole-block permutation? (default: 0)
% nlx : Force NOT using lexicographic permutations,
%       doing simple CMC instead
% 
% Outputs:
% P   : A cell array containing all the nP permutation matrices.
%       The first matrix is always I.
%
% Reference for the permutation algorithm:
% * Knuth, D. E. "Generating All Tuples and Permutations".
%   The Art of Computer Programming. Volume 4, Fascicle 2.
%   Addison-Wesley, 2005
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Feb/2012 (first version)
% Feb/2013 (this version)
% http://brainder.org

% Defaults
nP   = 20000; % Number of permutations
pB   = false; % Permute blocks
nlx  = false; % Permutations not lexicographic
rthr = 100;   % Threshold to define whether using lexicographic or not

% Accept inputs
if nargin < 1 || nargin > 5,
    error('Incorrect number of arguments');
end
if nargin >= 1, S0  = varargin{1}(:); end
if nargin >= 2, nP  = varargin{2}; end
if nargin >= 3, blk = varargin{3}; else blk = ones(numel(S0),1); end
if nargin >= 4, pB  = varargin{4}; end
if nargin == 5, nlx = varargin{5}; end

% Some preliminaries
if pB,
    [tmp,idx1] = sortrows([blk S0]);
    [~,idxback] = sort(idx1);
    S = tmp(:,2);
    nBtmp = numel(unique(blk));
    Bsize = numel(S)/nBtmp;
    S = reshape(S,[Bsize nBtmp])';
    S = mat2seq(S);
    blks = ones(size(S));
    B = 1;
else
    [blks,idx1] = sort(blk);
    [~,idxback] = sort(idx1);
    S = S0(idx1);
    B = unique(blks); % unique, sorted block indices
end
nB = numel(B);     % number of blocks
N  = zeros(nB,1);  % to store the number of elements for each block
T  = cell(nB,1);   % to contain the original values in S (or Sb) replaced for integers
Pb = cell(nB,1);   % to store temporary permutations (idx form) for each block

% Number of unique permutations
maxP = zeros(1,nB);
for b = 1:nB, % for each block
    Sb   = S(blks == B(b));
    N(b) = numel(Sb);
    T{b} = nan(size(Sb)); % init with NaN to avoid actual numbers that may cause troubles
    U    = unique(Sb);
    nU   = zeros(size(U));
    for u = 1:numel(U),
        nU(u) = sum(Sb == U(u)); % count
        T{b}(Sb == U(u)) = u;
    end
    maxP(b) = factorial(N(b))/prod(factorial(nU));
end
maxPall = prod(maxP);

% Make sure there are no more requested permutations than
% actual possible permutations
if nP > maxPall || nP <= 0,
    nP = round(maxPall);
end

% Decide whether run lexicographic or not
if maxPall/nP < rthr && ~nlx,
    lex = true;
    
    % For each block
    for b = 1:nB,
        
        % Algorithm L
        % The 1st permutation is the sorted
        n = N(b);        % number of elements in this block
        a = sort(T{b});  % sort
        p = 1;           % current permutation number
        Pb{b} = zeros(n,maxP(b));  % init the array to store all perms
        Pb{b}(:,p) = a;   % 1st permutation is the identity
        
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
            Pb{b}(:,p) = a;
            
            % Step L2 again
            j = n - 1;
            while j > 0 && a(j) >= a(j+1),
                j = j - 1;
            end
        end
    end
    
    % Once the all permutations for each block have been listed,
    % make a table with all the unique possible combinations
    Ptab = zeros(prod(maxP),nB);
    tmp = [1 maxP 1];
    for b = 2:nB+1,
        kr = kron((1:tmp(b))',ones(prod(tmp(1:b-1)),1));
        Ptab(:,b-1) = repmat(kr,[prod(tmp(b+1:end)) 1]);
    end
    
    % Get rid of some randomly selected permutations
    % until we have nP. The 1st permutation (identity)
    % is always kept.
    if maxPall ~= size(Ptab,1),
        error('Something went wrong! Check the code!');
    end
    [~,tmp] = sort(rand(maxPall-1,1));
    idxout = [true(maxPall-nP,1); false(nP-1,1)];
    Ptab([false; idxout(tmp)],:) = [];
    
    % Correct the indices that currently range 1..N(b) so that
    % they reflect the indices of the sorted overall sample
    for b = 1:nB,
        Pb{b} = Pb{b} + sum(N(1:b-1));
    end
    
    % Now assemble the permutation indices!
    Pidx = [];
    for b = 1:nB,
        Pidx = [Pidx; Pb{b}(:,Ptab(:,b))]; %#ok
    end
    
else
    lex = false;
    
    % For each block
    Pb = cell(nB,1);
    for b = 1:nB,

        % The first permutation is always 'no permutation'
        Pb{b} = (1:sum(N(b)))';
        nPtmp = 1;
        
        % Make sure there is no repetition of the permutation
        % mapping, i.e., the permutation itself is without
        % replacement, albeit if there are ties, there may be
        % effective repetitions (rare if maxP >> nP).
        % This procedure is good, but with within-block permutation
        % it may cause an infinite loop (e.g., exhaustive within-block
        % may never reach the desired nP. So, put an if inside.
        while nPtmp < nP,
            [~,Pfill] = sort(rand(N(b),nP-nPtmp));
            Pb{b} = [Pb{b} Pfill]; % not a problem to grow in the loop
            if nB == 1,
                Pb{b} = unique(Pb{b}','rows')'; % not a problem to change here either
            end
            nPtmp = size(Pb{b},2);
        end
        
        % Correct the indices that currently range 1..N(b) so that
        % they reflect the indices of the sorted overall sample
        Pb{b} = Pb{b} + sum(N(1:b-1));
    end
    
    % Now assemble the permutation indices
    Pidx = [];
    for b = 1:nB,
        Pidx = [Pidx; Pb{b}]; %#ok
    end
end

% Convert to permutation matrices. This part could have been
% integrated into the while-loops above, and it would be faster
% and use less memory. But it would make the code less clear.
[~,tmp]  = sort(Pidx);
[~,idx2] = sort(tmp);
P = cell(nP,1);
for p = 1:nP,
    P{p} = palm_idx2perm(idx2(:,p));
    if pB, P{p} = kron(P{p},eye(Bsize)); end
    P{p} = P{p}(idxback,idxback);
end
