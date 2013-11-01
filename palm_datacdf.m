function pvals = palm_datacdf(G,Gvals)
% Compute the p-values for a set of statistics G, taking
% as reference a set of observed values for G, from which
% the empirical cumulative distribution function (cdf) is
% generated, or using a custom cdf.
%
% Usage:
% pvals = datacdf(G,Gvals)
%
% Inputs:
% G     : Vector of Nx1 statistics to be converted to p-values
% Gvals : A Mx1 vector of observed values for the same statistic
%         from which the empirical cdf is build and p-values
%         obtained.
%
% Output:
% pvals : P-values.
%
% This function is a simplification of the much more generic
% 'cdftool.m' so that only the 'data' option is retained.
% To increase speed, there is no argument checking.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / Univ. of Oxford
% Jul/2012 (original version)
% Aug/2013 (this version)
% http://brainder.org

% Sort the data and compute the empirical distribution
[~,cdfG,distp] = palm_csort(Gvals(:),'descend',true);
cdfG  = unique(cdfG);
distp = flipud(unique(distp))./numel(Gvals);

% Convert the data to p-values
pvals = zeros(size(G));
for g = numel(cdfG):-1:1,
    pvals(G <= cdfG(g)) = distp(g);
end
pvals(G > cdfG(end)) = 0;
