function Y2d = palm_convNto2(YNd,N)
% Convert an N-D dataset into a 2D array, choosing which dimension
% becomes the 1st dimension (rows) of the 2D array. All the remaining
% dimensions are unwrapped into the 2nd dimension (columns). For example,
% with a 4D input (x,y,z,t) and N = 4, the output is 2D of size (t, x*y*z).
% 
% Usage:
% Y2d = palm_convNto2(YNd,N);
% 
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Sep/2012 (first version)
% May/2026 (this version)
% http://brainder.org

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% PALM -- Permutation Analysis of Linear Models
% Copyright (C) 2015 Anderson M. Winkler
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% First version:
% tmp = permute(Y4d,[4 1 2 3]);
% siz = size(tmp);
% Y2d = reshape(tmp,[size(tmp,1) prod(siz(2:end))]);

% Second version:
% Y2d = reshape(YNd,...
%     numel(YNd)/size(YNd,ndims(YNd)),...
%     size(YNd,ndims(YNd)))';

% Third version, more general:
%if nargin < 2 || isempty(N)
%    N = ndims(YNd);
%end
nd   = max(ndims(YNd),N);
d    = 1:nd;
d(N) = []; % avoids setdiff (may fail in Octave)
perm = [N d];
Y2d  = reshape(permute(YNd,perm),size(YNd,N),[]);