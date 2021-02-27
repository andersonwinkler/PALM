function map = coolhot5(siz)
% Create a colormap that goes in this order:
% Cyan - Blue | Red - Yellow
% 
% Usage:
% coolhot5(siz)
% 
% siz : Size of the colorbar. To ensure symmetry,
%       an even number is recommended.
%
% _____________________________________
% Anderson M. Winkler
% Yale University / Institute of Living
% Aug/2011
% http://brainder.org

if nargin < 1, siz = size(get(gcf,'colormap'),1); end

p = floor(siz/2);

map = [
    zeros(1,p)      ones(1,p);        % Red
    linspace(1,0,p) linspace(0,1,p);  % Green
    ones(1,p)       zeros(1,p)]';     % Blue
