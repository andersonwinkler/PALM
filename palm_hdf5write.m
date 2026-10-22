function palm_hdf5write(filename,datapath,data,level)
% Write an N-D array to a dataset within an HDF5 file.
%
% Usage:
%
% palm_hdf5write(filename,datapath,data,level);
%
% filename : HDF5 file to write to (created if it does not exist)
% datapath : Path to the dataset within the HDF5 file, starting with "/"
% data     : Data to be written (numeric; logical is stored as uint8)
% level    : Optional gzip level, an integer from 0 to 9. 0, or omitting
%            the argument, stores the dataset uncompressed. A positive
%            level chunks the array and uses HDF5 deflate.
%
% If the file does not exist it is created. If it exists, the dataset at
% "datapath" is overwritten if already present, or added otherwise; any
% other datasets already in the file are left unchanged.
%
% See also: palm_hdf5read.m
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

if nargin < 3 || nargin > 4
    error('Incorrect number of arguments.');
end
if nargin < 4 || isempty(level)
    level = 0;
end
if ~(isnumeric(level) && isscalar(level) && isfinite(level) && ...
        level == fix(level) && level >= 0 && level <= 9)
    error('The compression level must be an integer between 0 and 9.');
end
level = double(level);
if isempty(datapath) || datapath(1) ~= '/'
    error('The dataset path must start with "/": %s',datapath);
end

% Swap to row-major
data = permute(data,ndims(data):-1:1);

% HDF5 numeric datasets carry no MATLAB class metadata, so logical data is
% stored (and later read back) as uint8. Non-numeric data is unsupported.
if islogical(data)
    data = uint8(data);
elseif ~ isnumeric(data)
    error('Only numeric or logical data can be written; received type (class) "%s".',class(data));
end

% Write the dataset. HDF5 support is patchy in Octave, let's adapt...
if palm_isoctave

    if exist('h5write') %#ok<EXIST>
        ds       = find_dataset(filename,datapath);
        existing = ~ isempty(ds);
        if existing
            same_shape   = isequal(ds.Dataspace.Size,size(data));
            same_type    = strcmp(class(data),h5_class_to_matlab(ds.Datatype));
            same_deflate = dataset_deflate(ds) == level;
            if same_shape && same_type && same_deflate
                % A matching dataset is already in place.
            elseif exist('h5delete') %#ok<EXIST>
                % HDF5 cannot resize, retype, or add a filter in place, so
                % drop the old link and recreate. Unlinking does not return
                % the old bytes to the file.
                h5delete(filename,datapath);
                existing = false;
            else
                error([
                    'Dataset %s already exists in file %s, but this Octave installation does not have the h5delete command.\n',...
                    'You have to either delete either the dataset or the entire file first.\n'], ...
                    datapath,filename);
            end
        end
        if ~ existing
            create_dataset(filename,datapath,data,level);
        end
        h5write(filename,datapath,data);
    else % if no hdf5oct functions
        % Octave without hdf5oct package: fallback using Octave own
        % "-hdf5" format (not interchangeable with standard HDF5 readers).
        warning('The package "hdf5oct" is not installed. Thus, there is only partial support for HDF5 files.')
        octave_h5write(filename,datapath,data);
    end

else

    % In MATLAB, we look up whether the dataset exists (and whether it 
    % matches in size and shape) to decide whether overwriting or
    % creating a dataset.
    ds       = find_dataset(filename,datapath);
    existing = ~ isempty(ds);
    if existing
        same_shape   = isequal(ds.Dataspace.Size,size(data));
        same_type    = strcmp(class(data),h5_class_to_matlab(ds.Datatype));
        same_deflate = dataset_deflate(ds) == level;
        if ~ (same_shape && same_type && same_deflate)
            % HDF5 cannot resize, retype, or add a filter in place, so
            % drop the old link and recreate. Unlinking does not return
            % the old bytes to the file.
            remove_link(filename,datapath);
            existing = false;
        end
    end
    if ~ existing
        create_dataset(filename,datapath,data,level);
    end
    % A matching dataset (or the freshly created one) is now in place, so
    % h5write simply overwrites the values.
    h5write(filename,datapath,data);
end

% ------------------------------------------------------------------
function create_dataset(filename,datapath,data,level)
% Create a dataset. A positive level gzip-compresses it. The whole array
% is one chunk unless that would be large, in which case the longest
% dimension is split into pieces of a few megabytes.
if level > 0
    h5create(filename,datapath,size(data), ...
        'Datatype',class(data), ...
        'ChunkSize',chunk_size(data), ...
        'Deflate',level);
else
    h5create(filename,datapath,size(data),'Datatype',class(data));
end

% ------------------------------------------------------------------
function chunks = chunk_size(data)
% Chunk shape for a compressed dataset. Small arrays are one chunk.
chunks = size(data);
switch class(data)
    case {'double','uint64','int64'}, nbytes = 8;
    case {'single','uint32','int32'}, nbytes = 4;
    case {'uint16','int16'},          nbytes = 2;
    otherwise,                        nbytes = 1;
end
nbytes = nbytes * numel(data);
limit  = 4 * 1024 * 1024;
if nbytes <= limit
    return
end
[~,dim] = max(chunks);
nchunk = ceil(nbytes / limit);
chunks(dim) = max(1,ceil(chunks(dim) / nchunk));

% ------------------------------------------------------------------
function level = dataset_deflate(ds)
% gzip level stored on an h5info dataset, or 0 when the dataset is not
% compressed. Octave's h5info reports Deflate; MATLAB reports Filters.
level = 0;
if isfield(ds,'Deflate') && ~isempty(ds.Deflate) && isnumeric(ds.Deflate)
    level = double(ds.Deflate);
elseif isfield(ds,'Filters') && ~isempty(ds.Filters)
    for k = 1:numel(ds.Filters)
        if isfield(ds.Filters(k),'Name') && strcmpi(ds.Filters(k).Name,'deflate') ...
                && isfield(ds.Filters(k),'Data') && ~isempty(ds.Filters(k).Data)
            level = double(ds.Filters(k).Data(1));
        end
    end
end

% ------------------------------------------------------------------
function ds = find_dataset(filename,datapath)
% Return the h5info dataset struct at datapath, or [] when it is absent.
% Uses only the file-level h5info, which (unlike h5info(filename,datapath))
% does not raise when the dataset is missing.
ds = [];
if exist(filename,'file') ~= 2
    return
end
parts = strsplit(datapath,'/');
parts = parts(~cellfun('isempty',parts));
ds    = search_node(h5info(filename),parts);

% ------------------------------------------------------------------
function ds = search_node(node,parts)
% Walk the h5info hierarchy following the path components in "parts".
ds = [];
if isempty(parts)
    return
end
if numel(parts) == 1 %#ok<ISCL>
    for i = 1:numel(node.Datasets)
        if strcmp(local_name(node.Datasets(i).Name),parts{1})
            ds = node.Datasets(i);
            return
        end
    end
else
    for i = 1:numel(node.Groups)
        if strcmp(local_name(node.Groups(i).Name),parts{1})
            ds = search_node(node.Groups(i),parts(2:end));
            return
        end
    end
end

% ------------------------------------------------------------------
function name = local_name(fullname)
% Last component of a (possibly "/"-separated) HDF5 object name. Group names
% from h5info are full paths, whereas dataset names are already local; this
% handles both.
idx = find(fullname == '/',1,'last');
if isempty(idx)
    name = fullname;
else
    name = fullname(idx+1:end);
end

% ------------------------------------------------------------------
function remove_link(filename,datapath)
% Delete an existing dataset (link) using the low-level interface, so a
% dataset of a different size or type can be recreated.
% onCleanup guarantees the file is closed even if the delete fails.
fid = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT');
closer = onCleanup(@() H5F.close(fid));
H5L.delete(fid,datapath,'H5P_DEFAULT');

% ------------------------------------------------------------------
function cls = h5_class_to_matlab(dtype)
% Map an h5info Datatype struct to the matching MATLAB class name, or ''
% when it cannot be determined confidently. Returning '' simply forces a
% safe delete-and-recreate, so an unrecognised type is never written wrongly.
cls = '';
% Prefer the concrete named type when available: it encodes size and sign
% (e.g. "H5T_STD_U32LE"), covering both endiannesses.
if isfield(dtype,'Type') && ischar(dtype.Type)
    switch dtype.Type
        case {'H5T_IEEE_F64LE','H5T_IEEE_F64BE'}, cls = 'double';
        case {'H5T_IEEE_F32LE','H5T_IEEE_F32BE'}, cls = 'single';
        case {'H5T_STD_I8LE',  'H5T_STD_I8BE'},   cls = 'int8';
        case {'H5T_STD_U8LE',  'H5T_STD_U8BE'},   cls = 'uint8';
        case {'H5T_STD_I16LE', 'H5T_STD_I16BE'},  cls = 'int16';
        case {'H5T_STD_U16LE', 'H5T_STD_U16BE'},  cls = 'uint16';
        case {'H5T_STD_I32LE', 'H5T_STD_I32BE'},  cls = 'int32';
        case {'H5T_STD_U32LE', 'H5T_STD_U32BE'},  cls = 'uint32';
        case {'H5T_STD_I64LE', 'H5T_STD_I64BE'},  cls = 'int64';
        case {'H5T_STD_U64LE', 'H5T_STD_U64BE'},  cls = 'uint64';
    end
end
% Fallback on Class + Size for floats when the named type is not exposed.
% Integer sign cannot be inferred from these two alone, so integers are left
% as '' (forcing a safe recreate) rather than guessed.
if isempty(cls) && isfield(dtype,'Class') && isfield(dtype,'Size') && strcmp(dtype.Class,'H5T_FLOAT')
    switch dtype.Size
        case 8, cls = 'double';
        case 4, cls = 'single';
    end
end

% ------------------------------------------------------------------
function octave_h5write(filename,datapath,data)
% Fallback used when hdf5oct is unavailable. Saves with Octave's native
% "-hdf5" format, building a nested struct whose fields are the components
% of the dataset path so that load('-hdf5',...) in octave_h5read can walk
% the same path back. To preserve datasets that share part of the path and
% to overwrite an existing one, the current file contents are loaded and
% merged before the whole file is rewritten.
parts = strsplit(datapath,'/');
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
if exist(filename,'file') == 2
    S = load('-hdf5',filename);
    S = merge_structs(S,node,parts);
    save('-hdf5',filename,'-struct','S');
else
    save('-hdf5',filename,'-struct','node');
end

% ------------------------------------------------------------------
function S = merge_structs(S,node,parts)
% Set S.(parts{1}).(parts{2})... to the leaf carried by "node", creating
% intermediate structs as needed and preserving existing sibling fields.
key = parts{1};
if numel(parts) == 1 %#ok<ISCL>
    S.(key) = node.(key);
else
    if isfield(S,key) && isstruct(S.(key))
        sub = S.(key);
    else
        sub = struct();
    end
    S.(key) = merge_structs(sub,node.(key),parts(2:end));
end