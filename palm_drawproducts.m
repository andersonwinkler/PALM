function B = palm_drawproducts(P,S,nB)
% Merge a set P of permutations matrices with a set S of
% sign-flipping matrices. A new set B with nB elements is
% produced with random products P*S.
% 
% Usage:
% B = mergeps(P,S,nB)
% 
% P  : Cell array with permutation matrices.
% S  : Cell array with sign-flipping matrices.
% nB : Number of random products P*S to be generated.
% B  : Cell-array with nB output matrices.
% 
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Jul/2013
% http://brainder.org

nP = numel(P);
nS = numel(S);
nBall = nP*nS;

if nB > nBall,
    nB = nBall;
end

B = cell(nB,1);
b = 1;
if nB == nBall,
    for p = 1:nP,
        for s = 1:nS,
            B{b} = P{p}*S{s};
            b = b + 1;
        end
    end
else
    idx = [];
    while numel(idx) < nB,
        idx = vertcat(idx,randi(nBall,nB-numel(idx),1)); %#ok
        idx = unique(idx);
    end
    [i,j] = ind2sub([nP nS],idx);
    for p = 1:i,
        for s = 1:j,
            B{b} = P{p}*S{s};
            b = b + 1;
        end
    end
end