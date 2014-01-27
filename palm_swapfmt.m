function Pnew = palm_swapfmt(Pset)
% Convert a set of permutation matrices to an array
% of permutation indices and vice versa.
%
%         Cell array <===> Shuffling indices
%  Shuffling indices <===> Cell array
% 
% Pnew = palm_swapfmt(Pset)
% 
% Pset : Set of permutations, sign-flips or both.
%        This can be supplied either as a cell array
%        of (sparse) permutation matrices, or an
%        array of permutation indices.
% Pnew : The converted set of permutations.
% 
% _____________________________________
% Anderson Winkler
% FMRIB / University of Oxford
% Dec/2013
% http://brainder.org

if iscell(Pset),
    Pnew = zeros(size(Pset{1},1),numel(Pset));
    for p = 1:numel(Pset),
        Pnew(:,p) = palm_perm2idx(Pset{p});
    end
    if size(unique(abs(Pnew)','rows'),1) == 1;
        Pnew = sign(Pnew);
    end
else
    Pnew = cell(size(Pset,2),1);
    for p = 1:size(Pset,2),
        sgn = sign(Pset(:,p));
        idx = abs(Pset(:,p));
        if all(true(size(idx)) == idx),
            Pnew{p} = sparse(diag(sgn));
        else
            Pnew{p} = sparse(diag(sgn))*palm_idx2perm(idx);
        end
    end
end
