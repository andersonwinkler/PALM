function [Bset,nB] = palm_makeshuf(opts,plm)
% This is a wrapper for the palm_permtree.m and palm_fliptree.m
% that generates a sigle set of permutations based on the options
% and parameters given in 'opts' and 'plm' structs.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Nov/2013
% http://brainder.org

% Make the permutation tree, regardless of what comes next.
% The design matrix is represented by the sequence in the 2nd arg
Ptree = palm_tree(plm.EB,plm.tmp.seq);

% Maximum number of shufflings (perms, sign-flips or both)
maxP = 1;
maxS = 1;
if opts.EE,
    maxP = palm_maxshuf(Ptree,'perms');
end
if opts.ISE,
    maxS = palm_maxshuf(Ptree,'flips');
end
maxB = maxP * maxS;

% Vars for below
Pset{1} = speye(plm.N);
Sset = Pset;

if opts.nP0 >= maxB,
    
    % Run exhaustively if the user requests too many permutations.
    % Give some warning though. Note that here CMC is irrelevant.
    warning([...
            'The selected number of shufflings (%g) is larger than or\n' ...
            '         equal to the maximum possible number (%g).\n' ...
            '         Computing all %g shufflings exhaustively.'],...
            opts.nP0,maxB,maxB);
    if opts.EE,
        Pset = palm_permtree(Ptree,maxP,[],maxP);
    end
    if opts.ISE,
        Sset = palm_fliptree(Ptree,maxS,[],maxS);
    end
    
elseif opts.nP0 < maxB,
    
    % Or use a subset of possible permutations. The nested conditions
    % are to avoid repetitions, and to compensate fewer flips with more
    % perms or vice versa as needed in the tight situations
    if opts.EE,
        if opts.nP0 >= maxP,
            Pset = palm_permtree(Ptree,maxP,opts.CMC,maxP);
        else
            Pset = palm_permtree(Ptree,opts.nP0,opts.CMC,maxP);
        end
    end
    if opts.ISE,
        if opts.nP0 >= maxS,
            Sset = palm_fliptree(Ptree,maxS,opts.CMC,maxS);
        else
            Sset = palm_fliptree(Ptree,opts.nP0,opts.CMC,maxS);
        end
    end
end

% Generate the set of shufflings, mixing permutations and
% sign-flippings as needed.
nP = numel(Pset);
nS = numel(Sset);
if nS == 1,
    Bset = Pset;
elseif nP == 1,
    Bset = Sset;
elseif opts.nP0 >= maxB,
    Bset = cell(maxB,1);
    b = 1;
    for p = 1:numel(Pset),
        for s = 1:numel(Sset),
            Bset{b} = Pset{p} * Sset{s};
            b = b + 1;
        end
    end
else
    Bset = cell(opts.nP0,1);
    if opts.CMC,
        for b = 1:opts.nP0,
            Bset{b} = Pset{randi(nP)} * Sset{randi(nS)};
        end
    else
        [~,idx] = sort(rand(nP*nS,1));
        idx = idx(1:opts.nP0);
        [pidx,sidx] = ind2sub([nP nS],idx);
        for b = 1:opts.nP0,
            Bset{b} = Pset{pidx(b)} * Sset{sidx(b)};
        end 
    end
end
nB = numel(Bset);
