function data = palm_hdf5read(filename,datapath)
% Read an N-D array from an HDF5 file, specified by its path within the file
% 
% Usage:
% 
% data = palm_hdf5read(filename,datapath);
% 
% filename : HDF5 filename
% datapath : Dataset path within HDF5
% data     : Data, as read
%
% See also: palm_hdf5write.m, palm_hdf5spec.m
%
% _____________________________________
% Anderson M. Winkler
% UTRGV
% May/2026
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

% The file must exist
if exist('isfile','builtin') ~= 0
    tf = isfile(filename);
else
    tf = exist(filename,'file') == 2;
end
if ~ tf
    error('File does not exist: %s',filename);
end

% Read the dataset (engine dependent)
if palm_isoctave && ~ exist('h5read') %#ok<EXIST>
    % Octave without hdf5oct: basic, memory-heavy fallback that
    % only understands Octave's own "-hdf5" files.
    warning('This Octave installation has incomplete support for HDF5 files.')
    data = octave_h5read(filename,datapath);
else
    % Octave with the hdf5oct package (shipped with PALM or 
    % installed by the user), or MATLAB
    data = h5read(filename,datapath);
end

% Ensure we have column-major (HDF5 stores as row-major)
data = permute(data,ndims(data):-1:1);

% ------------------------------------------------------------------
function raw = octave_h5read(filename,dataset)
% Fallback used when hdf5oct is unavailable. Loads the whole file with
% load('-hdf5',...) and walks the dataset path into the resulting struct
% tree. Counterpart of octave_h5write in palm_hdf5write.
S = load('-hdf5',filename);
parts = strsplit(dataset, '/');
parts = parts(~cellfun('isempty',parts));
node = S;
for k = 1:numel(parts)
    if ~isstruct(node) || ~isfield(node,parts{k})
        error('Component "%s" not found in path', parts{k});
    end
    node = node.(parts{k});
end
if isstruct(node)
    error('Path refers to a group, not a dataset');
end
raw = node;