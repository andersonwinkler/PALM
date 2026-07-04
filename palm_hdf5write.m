function palm_hdf5write(X)
% Write an HDF5 dataset from a structure X.
%
% Usage:
%
% palm_hdf5write(X);
%
% X.data     : Data to be written
% X.filename : HDF5 file to be created
% X.dataset  : Dataset path within the HDF5 file (must start with "/")
% X.permdim  : Ignored (kept for symmetry with palm_hdf5read.m)
%
% If the target file already exists it is deleted first, so that a
% fresh file is created.
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

if nargin < 1
    error('An input structure is required.');
end
if ~isstruct(X) || ~all(isfield(X,{'filename','dataset','data'}))
    error('Input must be a struct with fields "data", "filename" and "dataset".');
end

filename = strtrim(char(X.filename));
dataset  = strtrim(char(X.dataset));
if isempty(filename)
    error('The file name is empty.');
end
if isempty(dataset) || dataset(1) ~= '/'
    error('The dataset path is empty or does not start with "/": %s',dataset);
end

% If the file already exists, remove it so a clean file is created.
if exist('isfile','builtin') ~= 0
    tf = isfile(filename);
else
    tf = exist(filename,'file') == 2;
end
if tf
    delete(filename);
end

% Write the dataset (engine dependent)
data = X.data';
if palm_isoctave && ~ exist('h5write') %#ok<EXIST>
    % Octave without hdf5oct package
    % Basic functionality only, uses Octave's native HDF5 format
    % and should be read back with palm_hdf5read's own fallback.
    octave_h5write(filename,dataset,data);
else
    % Octave with hdf5oct package or MATLAB
    % Complete functionality, will work with any HDF5 tool.
    h5create(filename,dataset,size(data),'Datatype',class(data));
    h5write(filename,dataset,data);
end

% ------------------------------------------------------------------
function octave_h5write(filename, dataset, data)
% If the hdf5oct package is not available, save using Octave's native
% "-hdf5" format. This is the counterpart of octave_h5read in
% palm_hdf5read: it builds a nested struct whose fields are the
% components of the dataset path, then saves those fields so that
% load('-hdf5',...) can walk the same path back.
parts = strsplit(dataset,'/');
parts = parts(~cellfun('isempty',parts));
if isempty(parts)
    error('The dataset path is empty.');
end
node = data;
for k = numel(parts):-1:1
    tmp = struct();
    tmp.(parts{k}) = node;
    node = tmp;
end
% "node" is now a scalar struct; save its top-level fields into the file.
save('-hdf5',filename,'-struct','node');