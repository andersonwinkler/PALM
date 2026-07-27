function YNd = palm_conv2toN(Y2d,siz,N)
% Convert a 2D array back into an N-D dataset, reversing palm_convNto2.
% The 1st dimension (rows) of the 2D array is placed back at dimension N
% of the N-D output, and the 2nd dimension (columns) is rewrapped into
% all the other dimensions.
% 
% Usage:
% YNd = palm_conv2toN(Y2d,siz);
%
% Y2d : 2D data
% siz : Sizes up to dimension N-1
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
% tmp = reshape(Y2d,[size(Y2d,1) siz]);
% Y4d = permute(tmp,[2 3 4 1]);

% Second version:
% YNd = reshape(Y2d',[siz(:)' size(Y2d,1)]);

% Third version:
%if nargin < 3 || isempty(N)
%    N = numel(siz);
%end
nd            = max(numel(siz),N);
siz(end+1:nd) = 1;
d             = 1:nd;
d(N)          = [];
perm          = [N d];
YNd           = ipermute(reshape(Y2d,siz(perm)),perm);