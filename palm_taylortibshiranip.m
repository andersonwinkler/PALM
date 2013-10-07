function P = palm_taylortibshiranip(T,plm)
% Computes the parametric p-value for the combined statistic.
% It is calculated assuming that the statistic was computed 
% usign the function for the respective method.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Aug/2013
% http://brainder.org

P = normcdf(T,0,1/sqrt(plm.nY));
