function S = palm_mat2seq(M)
% Convert a matrix M into a sequence S of indices
% for each unique row. S is useful, e.g., to compute
% all possible lexicographic permutations of rows of M.
% 
% Usage:
% S = mat2seq(M)
% 
% M : A 2D array.
% S : Sequence of indices
% 
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Feb/2013
% http://brainder.org

N = size(M,1);
S = zeros(N,1);
U = unique(M,'rows');
for u = 1:size(U,1),
    uidx = all(bsxfun(@eq,U(u,:),M),2);
    S(uidx) = u;
end