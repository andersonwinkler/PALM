function X = palm_hdf5read(spec)
% Read an HDF5 dataset selected by a compact specification string of the
% form:
% 
%     /path/to/file.h5:/path/to/dataset[:N]
% 
% Usage:
% 
% X = palm_hdf5read(spec);
% 
% spec       : Specification string (see palm_hdf5spec)
% X.filename : HDF5 filename
% X.dataset  : Dataset path within HDF5
% X.data     : Data, as read
% X.permdim  : Permutation dimension (NaN if not supplied)
% X.readwith : Engine used to read the data
%
% See also: palm_hdf5spec.m
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

% Parse the specification string into its three components
[filename,dataset,permdim] = palm_hdf5spec(spec);

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
    % Octave without hdf5oct package
    % Basic functionality only, will use more memory and
    % will fail for some HDF5 files
    data     = octave_load(filename,dataset);
    readwith = 'octave-load-hdf5';
else
    % Octave with hdf5oct package or MATLAB
    % Complete functionality, will use less memory
    % and should work with any HDF5 file
    data     = h5read(filename,dataset);
    readwith = 'h5read';
end

% Dimensionality check (skipped when no permutation dimension was supplied)
if ~isnan(permdim) && ndims(data) < permdim
    error('Dataset "%s" has %d dimension(s), but you want to permute dimension %d.', ...
        dataset,ndims(data),permdim);
end

% Outputs
X.filename = filename;
X.dataset  = dataset;
X.permdim  = permdim;
X.data     = data';
X.readwith = readwith;

% ------------------------------------------------------------------
function raw = octave_load(filename, dataset)
% If hdf5oct package is not available, will open using "load",
% but this has incomplete support for HDF5. Fallback only.
% It will load the whole file, then walk the path into the struct tree.
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