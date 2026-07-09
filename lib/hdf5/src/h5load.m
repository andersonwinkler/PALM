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
# @deftypefn  {Function File} { @var{data} = } h5load (@var{filename})
# @deftypefnx {Function File} { @var{data} = } h5load (@var{filename},@var{location})
#
# Load an entire HDF5 file or a portion of it as a struct.
#
# @code{h5load(@var{filename})} loads all datasets in @var{filename} and returns
# the struct @var{data}. Datasets and Groups become fields of the structure reproducting
# the internal hierarchy of the HDF5 file. Attributes, datatypes and other components are
# not loaded.
#  
# @code{h5disp(@var{filename},@var{location})} loads all datasets in @var{filename} below
# the node specified by @var{location}. If @var{location} is a group, then all datasets
# and groups below it are returned in @var{data}. If @var{location} is a dataset, then only this
# dataset will be returned.
#
# This function is not provided by the MATLAB high-level HDF5 interface.
#
# @seealso{h5read}
# @end deftypefn

function data = h5load(filename,location)

if nargin<1 || nargin>2,
    print_usage();
endif

if nargin == 1
    location = "/";
endif

info = h5info(filename,location);


data = __h5load__(filename,info);

endfunction

function s = __dispAttr__(info, s)
n=numfields(info);
if n==0, return; endif
names = fieldnames(info);
for i=1:n
    disp_indented(names{i},getfield(info,names{i}),depth+1)
endfor
endfunction

function dataout = __h5load__(filename, info, datain)

if isfield(info,"Groups"), # info is a group
    g = struct(); # create a struct for the group
    datasets = info.Datasets;
    for i=1:size(datasets,1) # load all datasets as fields of the group
      g = __h5load__(filename,datasets(i),g);
    endfor
    groups = info.Groups; # load all groups as fields of the group
    for i=1:size(groups,1)
      g = __h5load__(filename,groups(i),g);
    endfor
    if nargin==3,  # means that this group is a child of another group
      dataout = datain;
      name = strsplit(info.Name,"/"){end};
      dataout.(name) = g; # store this group as a field of the parent struct
    else
      dataout = g; # else return the group directly
    endif
elseif isfield(info,"Dataspace"), # info is a dataset
    dset = h5read(filename,info.Name); # read the dataset
    if nargin==3, # means that this dataset is a child of a group
      dataout = datain;
      name = strsplit(info.Name,"/"){end}; # get just the name
      dataout.(name) = dset; # return as a field of the parent struct
    else
      dataout = dset; # else return the dataset directly
    endif
endif

endfunction