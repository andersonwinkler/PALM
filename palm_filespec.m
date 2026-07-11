function [filename,datapath,permdim] = palm_filespec(spec)
% Parse a compact file/dataset specification string of the form:
% 
%     /path/to/file.h5:/path/to/dataset[:N]
% 
% and return its three components separately.
% 
% where:
%  * /path/to/file.h5 is the path to the HDF5 file
%  * /path/to/dataset is the path within the HDF5 file that contains
%    the data
%  * N is an optional positive integer. In the context of permutations,
%    it is the dimension along which the data will be permuted. If
%    omitted, permdim is returned as NaN.
% 
% The spec must contain one or two ":" separators (two if N is given,
% one if it is omitted). A leading Windows drive letter (C:\... or
% C:/...) is recognised and its colon is not treated as a separator.
% 
% If the spec contains no ":" separator introducing an in-file dataset
% path (aside from a possible Windows drive letter), it is taken to refer
% to a non-HDF5 file. In that case the filename is returned unchanged and
% both dataset and permdim are returned as NaN. A caller can therefore
% test the second output (e.g. with ischar) to decide whether the spec
% refers to an HDF5 file.
% 
% Usage:
% 
% [filename, dataset, permdim] = palm_hdf5spec(spec);
% 
% spec     : Specification, as described above
% filename : Filename (returned unchanged for a non-HDF5 spec)
% datapath : Path to dataset within the file, or NaN if the spec does
%            not refer to an HDF5 or MATLAB .mat file
% permdim  : Permutation dimension (NaN if not supplied, or if the spec
%            does not refer to an HDF5 or .mat file)
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
if isempty(spec)
    error('The specification string is empty.');
end

% Locate the ":" separators of filename.h5:/path/to/dataset[:N]
% A leading Windows drive letter (e.g. C:\ or C:/) contains a colon that
% is not a separator, so we exclude it from the count and the split. The
% relevant separators are then always the last one or two colons in the
% string: two when the permutation dimension is supplied, one when it is not.
hasDrive = ~isempty(regexp(spec,'^[A-Za-z]:[\\/]','once'));
colons   = strfind(spec,':');
nSep     = numel(colons) - double(hasDrive);
if nSep < 1
    % No in-file dataset path was given, so the spec does not refer to an
    % HDF5 file. Return the filename unchanged and flag the non-HDF5 case
    % by returning both dataset and permdim as NaN.
    filename = spec;
    datapath = NaN;
    permdim  = NaN;
    return
elseif nSep > 2
    error(['The specification has too many ":" separators. Use ' ...
        '"filename.h5:/path/to/dataset" or "filename.h5:/path/to/dataset:N" ' ...
        '(a leading Windows drive letter such as C:\\ is allowed and does not ' ...
        'count). Found %d separator(s): %s'],nSep,spec);
end

if nSep == 2
    % Permutation dimension supplied
    cPath    = colons(end-1);   % between the file name and the in-file path
    cDim     = colons(end);     % between the in-file path and the dimension count
    filename = spec(1:cPath-1);
    datapath = spec(cPath+1:cDim-1);
    dimStr   = strtrim(spec(cDim+1:end));
else
    % Permutation dimension omitted (single separator)
    cPath    = colons(end);     % between the file name and the in-file path
    filename = spec(1:cPath-1);
    datapath = spec(cPath+1:end);
    dimStr   = '';
end
if isempty(filename)
    error('The file name is empty: %s',spec);
end
if isempty(datapath)
    error('The dataset path or variable name is empty: %s',spec);
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