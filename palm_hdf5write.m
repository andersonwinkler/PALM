function palm_hdf5write(filename,datapath,data)
% Write an N-D array to a dataset within an HDF5 file.
%
% Usage:
%
% palm_hdf5write(filename,datapath,data);
%
% filename : HDF5 file to write to (created if it does not exist)
% datapath : Path to the dataset within the HDF5 file, starting with "/"
% data     : Data to be written (numeric; logical is stored as uint8)
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

if nargin ~= 3
    error('Incorrect number of arguments.');
end
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
            % As long as h5delete is not implemented, we cannot do the same
            % as in MATLAB
            error('Dataset %s already exists in file %s. Delete either the dataset or the file first.',datapath,filename)
        else
            % If the dataset doesn't already exist, we can create it and
            % write the file
            h5create(filename,datapath,size(data),'Datatype',class(data));
            h5write(filename,datapath,data);
        end
    else
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
        same_shape = isequal(ds.Dataspace.Size,size(data));
        same_type  = strcmp(class(data),h5_class_to_matlab(ds.Datatype));
        if ~ (same_shape && same_type)
            % A dataset is there but its shape or type differs. HDF5 cannot
            % resize/retype in place through the high-level interface, so
            % drop the old link and recreate. Other datasets stay intact.
            % (HDF5 does not reclaim the freed space, so replacing a dataset
            % of a different size or type grows the file a little.)
            remove_link(filename,datapath);
            existing = false;
        end
    end
    if ~ existing
        % Create new dataset
        h5create(filename,datapath,size(data),'Datatype',class(data));
    end
    % A matching dataset (or the freshly created one) is now in place, so
    % h5write simply overwrites the values.
    h5write(filename,datapath,data);
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