function I = palm_perm2idx(P)
% Convert a permutation matrix into
% a set of indices.
% _____________________________________
% Anderson Winkler
% FMRIB / University of Oxford
% Feb/2012
% http://brainder.org

I = 1:size(P,1);
I = P*I';
