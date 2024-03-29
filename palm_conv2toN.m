function YNd = palm_conv2toN(Y2d,siz)
% Convert a 2D array (t,x*y*z) into a N-D dataset (x,y,z,t,...)
% 
% Usage:
% YNd = conv2toN(Y2d,siz);
%
% Y2d : 2D data
% siz : Sizes up to dimension N-1
% 
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Sep/2012
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

% tmp = reshape(Y2d,[size(Y2d,1) siz]);
% Y4d = permute(tmp,[2 3 4 1]);

YNd = reshape(Y2d',[siz(:)' size(Y2d,1)]);
