##
##    Copyright (C) 2012 Tom Mullins
##    Copyright (C) 2015 Tom Mullins, Thorsten Liebig, Anton Starikov, Stefan Großhauser
##    Copyright (C) 2008-2013 Andrew Collette
##    Copyright (C) 2024 George Apostolopoulos
##
##    This file is part of hdf5oct.
##
##    hdf5oct is free software: you can redistribute it and/or modify
##    it under the terms of the GNU Lesser General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    hdf5oct is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU Lesser General Public License for more details.
##
##    You should have received a copy of the GNU Lesser General Public License
##    along with hdf5oct.  If not, see <http://www.gnu.org/licenses/>.
##

# -*- texinfo -*-
# @deftypefn  {Function File} { } h5disp (@var{filename})
# @deftypefnx {Function File} { } h5disp (@var{filename},@var{location})
# @deftypefnx {Function File} { } h5disp (@var{filename},@var{location},@var{mode})
#
# Display the contents of a HDF5 file.
#
# @code{h5disp(@var{filename})} displays the metadata that describes the structure of
# the specified HDF5 file.
#
# @code{h5disp(@var{filename},@var{location})} displays the metadata for 
# the specified location within the file.
#
# @code{h5disp(@var{filename},@var{location},@var{mode})}, where @var{mode} is 'min', 
# displays only the group and
# dataset names. The default for mode is 'simple', which displays the dataset
# metadata and attribute values.
#
# @seealso{h5info}
# @end deftypefn


function h5disp(filename,location,mode)

if nargin<1 || nargin>3,
    print_usage();
endif

if nargin == 1
    location = "/";
endif

if nargin < 3,
  mode = 'simple';
endif

mode = tolower(mode);
if !(strcmp(mode,'simple') || strcmp(mode,'min')),
  error("Invalid 'mode' option. Valid values: 'simple' or 'min'");
endif

info = h5info(filename,location);

disp(["HDF5 " filename]);
if strcmp(mode,'simple'),
  __disph5__(info,0);
else
  __disph5_min__(info,0);
endif


endfunction

function disp_indented(name,value,depth)
  S=char(strsplit(disp(value),'\n'));
  S=S(1:end-1,:);
  if size(S,1) == 1,
    disp([blanks(2*depth) "'" name "': " S]);
  else
    disp([blanks(2*depth) "'" name "':"]);
    disp([char(ones(size(S,1),2*(depth+1))*' ') S]);
  endif
endfunction

function s = disp_size(sz)
  if isempty(sz),
    s = "";
  elseif length(sz)==1,
    s = num2str(sz(1));
  else
    s = [num2str(sz(1)) num2str(sz(2:end),'x%d')];
  endif
endfunction

function __dispAttr__(info, depth)
n=numfields(info);
if n==0, return; endif
indent = blanks(2*depth);
disp([indent "Attributes (" num2str(n) "):"]);
names = fieldnames(info);
for i=1:n
    disp_indented(names{i},getfield(info,names{i}),depth+1)
endfor
endfunction

function __disph5__(info, depth)

indent = blanks(2*depth);

if isfield(info,"Groups"), # info is a group
    disp([indent "Group '" info.Name "'"]);
    __dispAttr__(info.Attributes, depth+1);
    datasets = info.Datasets;
    for i=1:size(datasets,1) __disph5__(datasets(i),depth+1); endfor
    groups = info.Groups;
    for i=1:size(groups,1) __disph5__(groups(i),depth+1); endfor
elseif isfield(info,"Dataspace"), # info is a dataset
    disp([indent "Dataset '" info.Name "'"]);
    __dispAttr__(info.Attributes, depth+1);
    disp([indent "  Extent: " num2str(info.Dataspace.ExtentType)]);
    if !isempty(info.Dataspace.Size), disp([indent "  Size: " disp_size(info.Dataspace.Size)]); endif
    if !isempty(info.Dataspace.MaxSize), disp([indent "  MaxSize: " disp_size(double(info.Dataspace.MaxSize))]); endif
    __disph5__(info.Datatype,depth+1);
    if !isempty(info.ChunkSize)
        disp([indent "  ChunkSize: [" num2str(info.ChunkSize) "]"]);
    endif
elseif isfield(info,"Class"), # info is a datatype
    disp([indent "Datatype"]);
    disp([indent "  Class: '" info.Class "'"]);
    disp([indent "  OctaveClass: '" info.OctaveClass "'"]);
    if isfield(info,"Size"), disp([indent "  Size: " num2str(info.Size) ]); endif
    if isfield(info,"Sign"), disp([indent "  Sign: " info.Sign]); endif
    if isfield(info,"charSet"), disp([indent "  charSet: " info.charSet]); endif
    if isfield(info,"Pading"), disp([indent "  Pading: " info.Pading]); endif
endif

endfunction

function __disph5_min__(info, depth)

indent = blanks(2*depth);

if isfield(info,"Groups"), # info is a group
    disp([indent "Group '" info.Name "'"]);
    datasets = info.Datasets;
    for i=1:size(datasets,1) __disph5_min__(datasets(i),depth+1); endfor
    groups = info.Groups;
    for i=1:size(groups,1) __disph5_min__(groups(i),depth+1); endfor
elseif isfield(info,"Dataspace"), # info is a dataset
    disp([indent "Dataset '" info.Name "'"]);
endif

endfunction








