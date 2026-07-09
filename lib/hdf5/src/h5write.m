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
# @deftypefn {Function File} { } h5write (@var{filename}, @var{dsetname}, @var{data})
# @deftypefnx {Function File} { } h5write (@var{filename}, @var{dsetname}, @var{data}, @var{start}, @var{count})
# @deftypefnx {Function File} { } h5write (@var{filename}, @var{dsetname}, @var{data}, @var{start}, @var{count}, @var{stride})
#
# Write data to a HDF5 dataset.
#
# @code{h5write(filename,dsetname,data)} writes data to an entire dataset, dsetname,
# in the specified HDF5 file. If the dataset is fixed in size,
# the amount of data to be written must match the size of the dataset.
#
# @code{h5write(filename,ds,data,start,count)} writes a subset of data to a dataset,
# beginning at starting location start, and continuing for count elements.
# In a multidimensional dataset, count specifies a distance in each direction.
# @code{h5write} extends an extendable dataset along any unlimited dimensions,
# if necessary.
#
# @code{h5write(filename,ds,data,start,count,stride)} specifies the spacing between
# elements, @var{stride}, along each dimension of the dataset.
#
# @table @asis
# @item @var{filename}
# Filename of an existing HDF5 file, specified as a string
# @item @var{dsetname}
# The full path to an existing dataset, specified as a string.
# @item @var{data}
# Data to be written to the HDF5 file.
# If a numeric datatype was specified in the corresponding call to h5create,
# then data is a numeric matrix containing floating-point or integer data.
# Data must be non sparse, and must be the same size as the HDF5 dataset
# if you do not specify start or count.
# If a dimension in the dataset is unlimited, then the data to be written
# can be any size along that dimension.
#
# If "string" was specified as the datatype in the corresponding call to h5create,
# data is either a single string or a cell string array. 
# The array dimensions must match those specified
# in the call to h5create.
#
# @item @var{start}
# Starting location, specified as a numeric vector of positive integers.
# For an n-dimensional dataset, start is a vector of length n containing 1-based
# indices. The elements of start correspond, in order, to the dataset dimensions.
# If any dimension of ds is unlimited, you must specify start.
#
# If you do not specify start, then the h5write function starts writing
# to the dataset from the first index along each dimension.
#
# @item @var{count}
# Number of elements to write, specified as a numeric vector of positive integers. 
# For an n-dimensional dataset, count is a vector of length n, 
# specifying the number of elements to write to the dataset 
# along each dimension. If any dimension of ds is unlimited, 
# then count must be specified.
# @item @var{stride}
# Optional spacing between elements along each dimension. For an n-dimensional
# dataset, stride is a vector of length n. A value of 1 writes without 
# skipping elements in the corresponding dimension, a value of 
# 2 writes every other element, and so on.
# @end table
#
# @seealso{h5create}
# @end deftypefn
#

function h5write(filename,location,data,start_pos,count,stride)

# check number and types of arguments
if !((nargin == 3) || (nargin==5) || (nargin==6))
    print_usage();
endif
if (!ischar(filename))
  error("h5write: 1st argument must be a string holding the hdf5 file name");
endif
if (!isfile(filename))
  error("h5write: filename does not exist");
endif
if (!ischar(location))
  error("h5write: 2nd argument must be a string holding the dataset location");
endif

datasize = size(data);
datasize = datasize(:);

if nargin==3,
    start_pos = [];
    count = [];
    stride = [];
else
    start_pos = check_idx_vec(start_pos,datasize,'start');
    count = check_idx_vec(count,datasize,'count');
    if nargin==6,
        stride = check_idx_vec(stride,datasize,'stride');
    else
        stride = [];
    endif
endif

if ischar(data), data = cellstr(data); endif

__h5write__(filename,location,data,start_pos,count,stride);

endfunction

function j = check_idx_vec(i, sz, lbl)

if !(isvector(i) && isindex(i)),
    error(["h5write: " lbl " must be a vector of valid index values"]);
endif
j = i(:);
if length(j)==1,
  j = [j; 1];
endif
if length(j) != length(sz)
    error(["h5write: size of " lbl " incompatible with data"]);
endif

endfunction




