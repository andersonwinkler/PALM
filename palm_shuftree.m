function [Bset,nB] = palm_shuftree(varargin)
% This is a wrapper for the palm_permtree.m and palm_fliptree.m
% that generates a sigle set of permutations. It can also generate
% only permutations with sign-flipping depending on the input
% arguments.
% 
% Usage (style 1)
% [Bset,nB] = palm_shuftree(Ptree,N,nP0,CMC,EE,ISE,idxout)
% 
% Inputs:
% - Ptree  : Permutation tree.
% - N      : Number of observations.
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
% 
% Usage (style 2, to be used by the PALM main function):
% [Bset,nB] = palm_shuftree(opts,plm)
% 
% Inputs:
% - opts   : Struct with PALM options
% - plm    : Struct with PALM data
% 
% Outputs:
% - Bset   : Set of permutations and/or sign flips.
% - nB     : Number of permutations and/or sign-flips.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Nov/2013
% http://brainder.org

% Take arguments
if nargin == 2,
    opts    = varargin{1};
    plm     = varargin{2};
    EE      = opts.EE;
    ISE     = opts.ISE;
    N       = plm.N;
    nP0     = opts.nP0;
    CMC     = opts.CMC;
    Ptree   = palm_tree(plm.tmp.EB,plm.tmp.seq);
    idxout  = false;
elseif nargin == 4 || nargin == 6 || nargin == 7,
    Ptree   = varargin{1};
    N       = varargin{2};
    nP0     = varargin{3};
    CMC     = varargin{4};
    if nargin == 6 || nargin == 7,
        EE  = varargin{5};
        ISE = varargin{6};
    else
        EE  = true;
        ISE = false;
    end
    if nargin == 7,
        idxout = varargin{7};
    else
        idxout = false;
    end
else
    error('Incorrect number of input arguments');
end
if ~EE && ~ISE,
    error('EE and/or ISE must be enabled, otherwise there is nothing to shuffle.')
end

% Maximum number of shufflings (perms, sign-flips or both)
maxP = 1;
maxS = 1;
if EE,
    maxP = palm_maxshuf(Ptree,'perms');
    if maxP > 1.7976931348623158e308,
        fprintf('Number of possible permutations is >= 1.7976931348623158e308.\n');
    else
        fprintf('Number of possible permutations is %g.\n',maxP);
    end
end
if ISE,
    maxS = palm_maxshuf(Ptree,'flips');
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
Pset{1} = speye(N);
Sset = Pset;

% Generate the Pset and Sset
if nP0 == 0 || nP0 >= maxB,
    % Run exhaustively if the user requests too many permutations.
    % Note that here CMC is irrelevant.
    fprintf('Running %g shufflings (%s).\n',maxB,whatshuf);
    if EE,
        Pset = palm_permtree(Ptree,maxP,[],false,maxP);
    end
    if ISE,
        Sset = palm_fliptree(Ptree,maxS,[],false,maxS);
    end
elseif nP0 < maxB,
    % Or use a subset of possible permutations. The nested conditions
    % are to avoid repetitions, and to compensate fewer flips with more
    % perms or vice versa as needed in the tight situations
    fprintf('Running %g shufflings (%s).\n',nP0,whatshuf);
    if EE,
        if nP0 >= maxP,
            Pset = palm_permtree(Ptree,maxP,CMC,false,maxP);
        else
            Pset = palm_permtree(Ptree,nP0,CMC,false,maxP);
        end
    end
    if ISE,
        if nP0 >= maxS,
            Sset = palm_fliptree(Ptree,maxS,CMC,false,maxS);
        else
            Sset = palm_fliptree(Ptree,nP0,CMC,false,maxS);
        end
    end
end

% Generate the set of shufflings, mixing permutations and
% sign-flippings as needed.
nP = numel(Pset);
nS = numel(Sset);
if nS == 1,
    % If only 1 sign-flip is possible, ignore it.
    Bset = Pset;
elseif nP == 1,
    % If only 1 permutation is possible, ignore it.
    Bset = Sset;
elseif nP0 == 0 || nP0 >= maxB,
    % If the user requested too many shufflings, do all
    % those that are possible.
    Bset = cell(maxB,1);
    b = 1;
    for p = 1:numel(Pset),
        for s = 1:numel(Sset),
            Bset{b} = Pset{p} * Sset{s};
            b = b + 1;
        end
    end
else
    % The typical case, with an enormous number of possible
    % shufflings, and the user choses a moderate number
    Bset = cell(nP0,1);
    % 1st shuffling is no shuffling, regardless
    Bset{1} = Pset{1} * Sset{1};
    if CMC,
        % If CMC, no need to take care of repetitions.
        for b = 2:nP0,
            Bset{b} = Pset{randi(nP)} * Sset{randi(nS)};
        end
    else
        % Otherwise, avoid them
        [~,idx] = sort(rand(nP*nS,1));
        idx = idx(1:nP0);
        [pidx,sidx] = ind2sub([nP nS],idx);
        for b = 2:nP0,
            Bset{b} = Pset{pidx(b)} * Sset{sidx(b)};
        end 
    end
end
nB = numel(Bset);

% In the draft mode, the permutations can't be in lexicographic
% order, but entirely shuffled.
if nargin == 2 && opts.draft,
    Bset2 = cell(size(Bset));
    [~,idx] = sort(rand(nB,1));
    for p = 1:nB,
        Bset2{p} = Bset(idx(p));
    end
    Bset = Bset2;
end

% If the desired outputs are permutation indices instead
% of permutation matrices
if idxout,
    Bset = palm_swapfmt(Bset);
end
