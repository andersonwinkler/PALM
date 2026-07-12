function x = mpinv(P,y,sigsq,tail)
% Inverse cumulative density function of the Marcenko-Pastur
% distribution of the eigenvalues of random matrices.
% 
% Consider a matrix X (never seen by this function) with random complex
% entries with zero mean and variance sigsq, and with a ratio of
% number of rows by number columns converging to y. Let P be the
% cumulative density of some eigenvalue x of the covariance matrix
% of X. This function returns x.
% 
% Usage:
% x = mpinv(P,y,sigsq)
% 
% Inputs:
% P     : Cumulative density (i.e., 1-p, where p is the p-value).
% y     : Ratio #rows/#cols of the random matrix from which the
%         eigenvalues are characteristic.
% sigsq : Variance of the entries of X.
% tail  : Boolean indicating whether x are upper-tail probabilities
%         (p-values) or simply the CDF. Default is false.
% 
% Outputs:
% x     : Eigenvalue of the matrix X for the corresponding P.
% 
% References:
% * Bai ZD. Methodologies in spectral analysis of large dimensional
%   random matrices, a review. Statistica Sinica 1999;9(3):611?77. 
% * Marcenko VA, Pastur LA. Distribution of eigenvalues for some 
%   sets of random matrices. Math USSR Sb 1967;1(4):457?83. 
% 
% _____________________________________
% Anderson M. Winkler
% National Institutes of Health
% Mar/2021
% http://brainder.org

opts = optimset('TolX',1e-10);
a = sigsq.*(1 - y^.5).^2;
b = sigsq.*(1 + y^.5).^2;
x = fminbnd(@(x)abs(mpcdf(x,y,sigsq,tail)-P),a,b,opts);