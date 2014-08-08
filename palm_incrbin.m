function B = palm_incrbin(B)
% Increment a binary number by 1. 
% 
% Usage:
% Bi = palm_incrbin(B)
%
% B  : A logical vector, with the least
%      significant digits first.
% Bi : Incremented B.
% 
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Feb/2014
% http://brainder.org

% Although this function is tiny, it's used more than
% once, so it's better to have as a separate file.

k = find(~B,1);
B(1:k) = ~B(1:k);