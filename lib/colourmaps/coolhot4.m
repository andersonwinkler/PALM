function map = coolhot4(siz)
% Create a colormap that goes in this order:
% Blue - White - Red
% 
% Usage:
% coolhot4(siz)
% 
% siz : Size of the colorbar. To ensure symmetry,
%       an even number is recommended.
%
% _____________________________________
% Anderson M. Winkler
% Yale University / Institute of Living
% May/2011
% http://brainder.org

if nargin < 1, siz = size(get(gcf,'colormap'),1); end

p = floor(siz/2);

map = [
    linspace(0,1,p) ones(1,p);        % Red
    linspace(0,1,p) linspace(1,0,p);  % Green
    ones(1,p) linspace(1,0,p)]';      % Blue
