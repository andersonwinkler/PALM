function [G,Gdist] = palm_prepaccel(G,Gdist)
% Prepare the test statistic and the null distribution for the
% tail and gamma approximations.
%
% Usage:
% [G,Gdist] = palm_prepaccel(G,Gdist)
% 
% - G     : Test statistic.
% - Gdist : Reference distribution. The first value is the non-permuted.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Dec/2015
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

% Make sure all values are positive.
d     = min(G,min(Gdist));
G     = G - mi + 1;
Gdist = Gdist - mi + 1;

% Remove the non-permuted case.
Gdist = Gdist(2:end);

% Do a Box-Cox transformation.
[Gdist,L] = palm_boxcox(Gdist);
G         = palm_boxcox(G,L);
