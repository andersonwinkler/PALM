function [Bset,nB] = palm_shufsimple(varargin)
% A single function to generate a set of permutations and/or
% sign-flips. This function is a faster replacement to
% palm_shuftree.m when all observations are freely exchangeable,
% i.e., when there are no block restrictions and no tree needs
% to be constructed.
% 
% Usage
% [Bset,nB] = palm_shuftree(M,nP0,CMC,EE,ISE,idxout)
% 
% Inputs:
% - M      : Design matrix.
% - nP0    : Requested number of permutations.
% - CMC    : Use Conditional Monte Carlo.
% - EE     : Allow permutations?
% - ISE    : Allow sign-flips?
%            If you supply the EE argument, you must
%            also supply ISE argument. If one is omited,
%            the other needs to be omited too.
%            Default is true for EE, and false for ISE.
% - idxout : (Optional) If true, the output isn't a cell
%            array with permutation matrices, but an array
%            with permutation indices.
% 
% Outputs:
% - Bset   : Set of permutations and/or sign flips.
% - nB     : Number of permutations and/or sign-flips.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Jan/2014
% http://brainder.org

% Accept arguments
if nargin < 2 || nargin > 6 || nargin == 4,
    error('Incorrect number of arguments');
end
M   = varargin{1};
nP0 = varargin{2};
if nargin > 2,
    CMC = varargin{3};
else
    CMC = false;
end
if nargin > 4,
    EE  = varargin{4};
    ISE = varargin{5};
else
    EE  = true;
    ISE = false;
end
if nargin > 5,
    idxout = varargin{6};
else
    idxout = false;
end
if ~EE && ~ISE,
    error('EE and/or ISE must be enabled, otherwise there is nothing to shuffle.')
end

% Sequence of unique values to shuffle
N = size(M,1);
[~,~,seq] = unique(M,'rows');
seq = horzcat(seq,(1:N)');
seqS = sortrows(seq);

% Number of unique permutations & sign flips
maxP = 1;
maxS = 1;
if EE,
    U = unique(seqS(:,1));
    nrep = zeros(size(U));
    for u = 1:numel(U),
        nrep(u) = sum(seqS(:,1) == U(u));
    end
    maxP = factorial(N)/prod(factorial(nrep));
    if maxP > 1.7976931348623158e308,
        fprintf('Number of possible permutations is >= 1.7976931348623158e308.\n');
    else
        fprintf('Number of possible permutations is %g.\n',maxP);
    end
end
if ISE,
    maxS = 2^N;
    if maxS > 1.7976931348623158e308,
        fprintf('Number of possible sign-flips is >= 1.7976931348623158e308.\n');
    else
        fprintf('Number of possible sign-flips is %g.\n',maxS);
    end
end
maxB = maxP * maxS;

% String for the screen output below
if EE && ~ISE,
    whatshuf = 'permutations only';
elseif ISE && ~EE,
    whatshuf = 'sign-flips only';
elseif EE && ISE,
    whatshuf = 'permutations and sign-flips';
end

% This ensures that there is at least 1 permutation (no permutation)
% and 1 sign-flipping (no sign-flipping). These are modified below as
% needed.
Pset = seqS(:,2);
Sset = ones(N,1);

% Generate the Pset and Sset
if nP0 == 0 || nP0 >= maxB,
    % Run exhaustively if the user requests too many permutations.
    % Note that here CMC is irrelevant.
    fprintf('Running %g shufflings (%s).\n',maxB,whatshuf);
    if EE,
        Pset = horzcat(Pset,zeros(N,maxP-1));
        for p = 2:maxP,
            seqS = palm_nextperm(seqS);
            Pset(:,p) = seqS(:,2);
        end
    end
    if ISE,
        Sset = palm_d2b(0:maxS-1,N)';
        Sset(~~Sset) = -1;
        Sset( ~Sset) =  1;
    end
elseif nP0 < maxB,
    % Or use a subset of possible permutations. The nested conditions
    % are to avoid repetitions, and to compensate fewer flips with more
    % perms or vice versa as needed in the tight situations
    fprintf('Running %g shufflings (%s).\n',nP0,whatshuf);
    if EE,
        if nP0 >= maxP,
            Pset = horzcat(Pset,zeros(N,maxP-1));
            for p = 2:maxP,
                seqS = palm_nextperm(seqS);
                Pset(:,p) = seqS(:,2);
            end
        else
            Pset = horzcat(Pset,zeros(N,nP0-1));
            if CMC,
                for p = 1:nP0,
                    Pset(:,p) = randperm(N)';
                end
            else
                Pseq = zeros(size(Pset));
                Pseq(:,1) = seqS(:,2);
                for p = 2:nP0,
                    whiletest = true;
                    while whiletest,
                        Pset(:,p) = randperm(N)';
                        Pseq(:,p) = seqS(Pset(:,p));
                        whiletest = any(all(bsxfun(@eq,Pseq(:,p),Pseq(:,1:p-1))));
                    end
                end
            end
        end
    end
    if ISE,
        if nP0 >= maxS,
            Sset = palm_d2b(0:maxS-1,N)';
            Sset(~~Sset) = -1;
            Sset( ~Sset) =  1;
        else
            if CMC,
                Sset = double(rand(N,nP0) > .5);
                Sset(:,1) = 0;
                Sset(~~Sset) = -1;
                Sset( ~Sset) =  1;
            else
                Sset = zeros(N,nP0);
                for p = 2:nP0,
                    whiletest = true;
                    while whiletest,
                        Sset(:,p) = rand(N,1) > .5;
                        whiletest = any(all(bsxfun(@eq,Sset(:,p),Sset(:,1:p-1))));
                    end
                end
                Sset(~~Sset) = -1;
                Sset( ~Sset) =  1;
            end
        end
    end
end

% Generate the set of shufflings, mixing permutations and
% sign-flippings as needed.
nP = size(Pset,2);
nS = size(Sset,2);
if nS == 1,
    % If only 1 sign-flip is possible, ignore it.
    Bset = Pset;
elseif nP == 1,
    % If only 1 permutation is possible, ignore it.
    Bset = bsxfun(@times,Pset,Sset);
elseif nP0 == 0 || nP0 >= maxB,
    % If the user requested too many shufflings, do all
    % those that are possible.
    Bset = zeros(N,maxB);
    b = 1;
    for p = 1:size(Pset,2),
        for s = 1:size(Sset,2),
            Bset(:,b) = Pset(:,p) .* Sset(:,s);
            b = b + 1;
        end
    end
else
    % The typical case, with an enormous number of possible
    % shufflings, and the user choses a moderate number
    Bset = zeros(N,nP0);
    % 1st shuffling is no shuffling, regardless
    Bset(:,1) = (1:N)';
    if CMC,
        % If CMC, no need to take care of repetitions.
        for b = 2:nP0,
            Bset(:,b) = Pset(:,randi(nP)) .* Sset(:,randi(nS));
        end
    else
        % Otherwise, avoid them
        [~,bidx] = sort(rand(nP*nS,1));
        bidx = bidx(1:nP0);
        [pidx,sidx] = ind2sub([nP nS],bidx);
        for b = 2:nP0,
            Bset(:,b) = Pset(:,pidx(b)) .* Sset(:,sidx(b));
        end
    end
end
nB = size(Bset,2);

% Sort back to the original order
Bset = sortrows(Bset);

% If the desired outputs are permutation indices instead
% of permutation matrices
if ~ idxout,
    Bset = palm_swapfmt(Bset);
end
