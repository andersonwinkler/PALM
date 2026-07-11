function map = coolhot1(siz)
% Create a colormap that is a combination of
% 'cool' and 'hot'.
% 
% Usage:
% coolhot(SIZ)
% 
% - SIZ = Size of the colorbar. To ensure symmetry,
%        an even number is recommended.
%
% _____________________________________
% Anderson M. Winkler
% Yale University / Institute of Living
% Jan/2011
% http://brainder.org

if nargin < 1, siz = size(get(gcf,'colormap'),1); end

c = cool(floor(siz/2));
c(:,1) = zeros(size(c,1),1);
h = hot(ceil(siz/2));
map = [c ; h];