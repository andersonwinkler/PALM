function map = coolhot2(siz)
% Create a colormap that goes in this order:
% Cyan - Blue - Black - Red - Yellow
% 
% Usage:
% coolhot2(SIZ)
% 
% - SIZ = Size of the colorbar. To ensure symmetry,
%        an even number is recommended.
%
% _____________________________________
% Anderson M. Winkler
% Yale University / Institute of Living
% May/2011
% http://brainder.org

if nargin < 1, siz = size(get(gcf,'colormap'),1); end

p = floor(siz/4);

map = [
    zeros(1,2*p) linspace(0,1,p) ones(1,p);        % Red
    linspace(1,0,p) zeros(1,2*p) linspace(0,1,p);  % Green
    ones(1,p) linspace(1,0,p) zeros(1,2*p)]';      % Blue

