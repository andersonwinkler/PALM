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
# @deftypefn {Function File} {@var{data}=} h5read (@var{filename}, @var{dsetname})
# @deftypefnx {Function File} {@var{data}=} h5read (@var{filename}, @var{dsetname}, @var{start}, @var{count})
# @deftypefnx {Function File} {@var{data}=} h5read (@var{filename}, @var{dsetname}, @var{start}, @var{count}, @var{stride})
#
# Read data from a HDF5 dataset.
#
# @code{data = h5read(@var{filename}, @var{dsetname})} 
# reads all the data from the dataset @var{dsetname}
# contained in the HDF5 file @var{filename}.
#
# @code{data = h5read(@var{filename}, @var{dsetname}, @var{start}, @var{count})} 
# reads a subset of data from the dataset
# beginning at the location specified in @var{start}.
# @var{count} specifies the number of elements to read along each dimension.
#
# @code{data = h5read(@var{filename}, @var{dsetname}, @var{start}, @var{count}, @var{stride})} 
# returns a subset of data with
# the interval between the indices of each dimension of the dataset specified by stride.
#
# Input arguments:
#
# @table @asis
# @item @var{filename}
# Filename of an existing HDF5 file, specified as a string
# @item @var{dsetname}
# Dataset name, specified as a string scalar or character vector containing
# the name of an existing dataset in the HDF5 file.
# @item @var{start}
# Starting location, specified as a numeric vector of positive integers.
# For an n-dimensional dataset, start is a vector of length n containing 1-based
# indices. The elements of start correspond, in order, to the dataset dimensions.
# If any dimension of ds is unlimited, you must specify start.
#
# If you do not specify start, then the h5write function starts writing
# to the dataset from the first index along each dimension.
# @item @var{count}
# Number of elements to write, specified as a numeric vector of positive integers. 
# For an n-dimensional dataset, count is a vector of length n, 
# specifying the number of elements to write to the dataset 
# along each dimension. If any dimension of ds is unlimited, 
# then count must be specified.
# @item @var{stride}
# Optional spacing between elements along each dimension. For an n-dimensional
# dataset, stride is a vector of length n. A value of 1 reads without 
# skipping elements in the corresponding dimension, a value of 
# 2 reads every other element, and so on.
# @end table
#
# @seealso{h5create, h5write}
# @end deftypefn
#

function data = h5read(filename,location,start_pos,count,stride)

# check number and types of arguments
if !((nargin == 2) || (nargin==4) || (nargin==5))
    print_usage();
endif
if (!ischar(filename))
  error("h5read: 1st argument must be a string holding the hdf5 file name");
endif
if (!isfile(filename))
  error("h5read: filename does not exist");
endif
if (!ischar(location))
  error("h5read: 2nd argument must be a string holding the dataset location");
endif

if nargin==2,
    start_pos = [];
    count = [];
    stride = [];
else
    start_pos = check_idx_vec(start_pos,'start');
    count = check_idx_vec(count,'count');
    if nargin==5,
        stride = check_idx_vec(stride,'stride');
    else
        stride = [];
    endif
    if length(count)!=length(start_pos),
        error(["h5write: count should have the same size as start"]);
    endif
    if nargin==6 && length(stride)!=length(start_pos),
        error(["h5write: stride should have the same size as start"]);
    endif
endif

data = __h5read__(filename,location,start_pos,count,stride);

endfunction

function j = check_idx_vec(i, lbl)

if !(isvector(i) && isindex(i)),
    error(["h5write: " lbl " must be a vector of valid index values"]);
endif

j = i(:);

endfunction




