function [Mr,Y] = palm_draperstoneman(P,Y,plm)
% This tiny function implements one of several
% different regression/permutation methods. See
% the function 'permglm.m' for a complete overview.
%
% Inputs:
% P   : A permutation or sign-flipping matrix
% Y   : Data not shuffled.
% plm : A struct with the GLM data and regressors.
%       See the function 'palm.m'
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

Mr = [P*plm.tmp.X plm.tmp.Z];
% note that Y remains unchanged in this method