function P = palm_fristonp(T,plm)
% Computes the parametric p-value for the combined statistic.
% It is calculated assuming that the statistic was computed 
% usign the function for the respective method.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Aug/2013
% http://brainder.org

% Save the 1-P
P = 1-(1-T).^plm.nY;