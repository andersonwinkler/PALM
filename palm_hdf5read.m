function X = palm_hdf5read(spec)
% Read an HDF5 dataset selected by a compact specification
% string of the form:
% 
%     /path/to/file.h5:/path/to/dataset[:N]
% 
% where:
%  * /path/to/file.h5 is the path to the HDF5 file that will be read
%  * /path/to/dataset is the path within the HDF5 file that contains
%    the data
%  * N is an optional integer that is returned unchanged. In the context
%    of permutations, it is the dimension along which the data will be
%    permuted. If omitted, permdim is returned as NaN.
% 
% The spec must contain one or two ":" separators (two if N is given,
% one if it is omitted). A leading Windows drive letter (C:\... or
% C:/...) is recognised and its colon is not treated as a separator.
% 
% Usage:
% 
% X = palm_hdf5read(spec);
% 
% spec       : Specification, as described above
% X.data     : Data, as read
% X.filename : HDF5 filename
% X.dataset  : Dataset path within HDF5
% X.permdim  : Permutation dimension (NaN if not supplied)
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
    error('A specification string is required.');
end
spec = strtrim(char(spec));

% Locate the ":" separators of filename.h5:/path/to/dataset[:N]
% A leading Windows drive letter (e.g. C:\ or C:/) contains a colon that
% is not a separator, so we exclude it from the count and the split. The
% relevant separators are then always the last one or two colons in the
% string: two when the permutation dimension is supplied, one when it is not.
hasDrive = ~isempty(regexp(spec,'^[A-Za-z]:[\\/]','once'));
colons   = strfind(spec,':');
nSep     = numel(colons) - double(hasDrive);
if nSep < 1 || nSep > 2
    error(['The specification must contain one or two ":" separators, as in ' ...
        '"filename.h5:/path/to/dataset" or "filename.h5:/path/to/dataset:N" ' ...
        '(a leading Windows drive letter such as C:\\ is allowed and does not ' ...
        'count). Found %d: %s'],nSep,spec);
end

if nSep == 2
    % Permutation dimension supplied
    cPath = colons(end-1);   % between the file name and the in-file path
    cDim  = colons(end);     % between the in-file path and the dimension count
    filename = spec(1:cPath-1);
    dataset  = spec(cPath+1:cDim-1);
    dimStr   = strtrim(spec(cDim+1:end));
else
    % Permutation dimension omitted (single separator)
    cPath = colons(end);     % between the file name and the in-file path
    filename = spec(1:cPath-1);
    dataset  = spec(cPath+1:end);
    dimStr   = '';
end
if isempty(filename)
    error('The file name is empty: %s',spec);
end
if isempty(dataset) || dataset(1) ~= '/'
    error('The dataset path is empty or does not start with "/": %s',spec);
end

% The permutation dimension is optional. When supplied (two separators) it
% must be a positive integer; when omitted (one separator) it defaults to NaN.
if nSep == 1
    permdim = NaN;
else
    if isempty(regexp(dimStr,'^\d+$','once'))
        error('The permutation dimension "%s" must be a positive integer.',dimStr);
    end
    permdim = str2double(dimStr);
    if ~(isfinite(permdim) && permdim >= 1)
        error('The permutation dimension must be a positive integer.');
    end
end

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
    data = octave_h5read(filename,dataset);
else
    % Octave with hdf5oct package or MATLAB
    % Complete functionality, will use less memory
    % and should work with any HDF5 file
    data = h5read(filename,dataset);
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

% ------------------------------------------------------------------
function raw = octave_h5read(filename, dataset)
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