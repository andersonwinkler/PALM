function P = palm_idx2perm(I)
% Convert a set of indices into a
% sparse permutation matrix.
% _____________________________________
% Anderson Winkler
% FMRIB / University of Oxford
% Feb/2012
% http://brainder.org

P = speye(numel(I));
P = P(I,:);
