function [Mr,Yr] = palm_huhjhun(P,Y,plm)
% This tiny function implements one of several
% different regression/permutation methods. See
% the function 'permglm.m' for a complete overview.
%
% Inputs:
% P   : A permutation or sign-flipping matrix.
% Y   : Data not shuffled.
% plm : A struct with the GLM data and regressors.
%       See the function 'palm.m'.
%
% Outputs:
% Mr  : Design matrix, ready for regression.
% Yr  : Data, ready for regression.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Aug/2013
% http://brainder.org

Mr = plm.tmp.Q'*plm.tmp.Rz*plm.tmp.X;
Yr = P'*plm.tmp.Q'*plm.tmp.Rz*Y;