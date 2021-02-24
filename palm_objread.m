function [vtx,fac,extra] = palm_objread(filename)
% Read a surface file, in ASCII format.
% 
% [vtx,fac,extra] = palm_objread(filename);
% 
% - vtx contains the coordinates (x,y,z), one vertex per row
% - fac contains the indices for the three vertices of each face
% - extra contains additional elements from the obj file (if they exist)
% 
% The indices for the vertices start at 1, not 0.
% 
% _____________________________________
% Anderson M. Winkler
% National Institutes of Health
% Feb/2020
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

extra = [];

% Read input file
fid = fopen(filename,'r');
obj = textscan(fid,'%s');
obj = obj{1};
fclose(fid);

% Geometric vertices
idx = find(strcmp(obj,'v'));
vtx = [str2double(obj(idx+1)) str2double(obj(idx+2)) str2double(obj(idx+3))];

% Faces
idx = find(strcmp(obj,'f'));
fac = [str2double(obj(idx+1)) str2double(obj(idx+2)) str2double(obj(idx+3))];
fac = uint64(fac);

% Vertex normals
idx = find(strcmp(obj,'vn'));
if ~ isempty(idx)
   extra.vn  = [str2double(obj(idx+1)) str2double(obj(idx+2)) str2double(obj(idx+3))];
end

% Vertex normals
idx = find(strcmp(obj,'vt'));
if ~ isempty(idx)
   extra.vt  = [str2double(obj(idx+1)) str2double(obj(idx+2)) str2double(obj(idx+3))];
end

% Parameter space vertices
idx = find(strcmp(obj,'vp'));
if ~ isempty(idx)
   extra.vp  = [str2double(obj(idx+1)) str2double(obj(idx+2))];
end

% Materials library
idx = find(strcmp(obj,'mtllib'));
if ~ isempty(idx)
   extra.mtllib  = obj(idx+1);
end

% Specific materials
idx = find(strcmp(obj,'usemtl'));
if ~ isempty(idx)
   extra.usemtl  = obj(idx+1);
end
