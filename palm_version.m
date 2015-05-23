function varargout = palm_version
% Read and output or print the version string.
% 
% _____________________________________
% Anderson Winkler
% FMRIB / University of Oxford
% Mar/2015
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

% Read the file with the version
fid = fopen(fullfile(fileparts(mfilename('fullpath')),'palm_version.txt'),'r');
vstr = textscan(fid,'%s');
fclose(fid);

% Assemble back as a string
vstr = sprintf('%s ',vstr{1}{:});
vstr = sprintf('%s\n',vstr(1:end-1));

% Print in the screen if no other output
if nargout == 0,
    fprintf(vstr);
else
    varargout{1} = vstr;
end