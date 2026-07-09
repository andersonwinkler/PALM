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
# @deftypefn {Function File} { } h5create (@var{filename}, @var{dsetname}, @var{size})
# @deftypefnx {Function File} { } h5create (@var{filename}, @var{dsetname}, @var{size}, @var{key}, @var{val},...)
#
# Create a HDF5 dataset.
#
# @code{h5create (@var{filename}, @var{dsetname}, @var{size})} 
# creates a dataset with name @var{dsetname} and size @var{size}
# in the HDF5 file specified by @var{filename}. 
#
# @code{h5create (@var{filename}, @var{dsetname}, @var{size}, @var{key}, @var{val},...)} 
# specifies one or more key-value arguments.
#
# Input arguments:
#
# @table @asis
# @item @var{filename}
# The path of the HDF5 file as a string. If the file does not exist it
# will be created.
# @item @var{dsetname}
# A string specifying the complete path to the dataset starting
# from the root group "/". Intermediate groups
# are created as necessary.
# @item @var{size}
# A vector specifying the dimensions of the dataset.
# @var{size} may contain one or several Inf values.
# This will lead to unlimited maximum extent of the dataset in the
# respective dimensions and 0 initial extent.
# Datasets with at least one unlimited dimension must be chunked.
# Chunking is generally recommended for large datasets.
#
# @code{size==1} or @code{size==[1 1]} results in a scalar (0-dimensional)
# dataset.
#
# Note that due to the difference in array storage layout between 
# OCTAVE (column-major)
# and HDF5 (row-major), arrays are saved transposed in the HDF5 file.
# Thus, if a @code{hdf5oct}-generated file is opened by another application
# the arrays will appear transposed. 
# @end table
#
# The list of @var{key}, @var{val} arguments allows to specify
# certain properties of the dataset. Allowed settings are:
#
# @table @asis
# @item @option{Datatype}
# one of the strings:
# @samp{double} (default) | @samp{single} | 
# @samp{double complex} | @samp{single complex} |
# @samp{uint64} | @samp{uint32} |
# @samp{uint16} | @samp{uint8} | @samp{int64} | @samp{int32} | @samp{int16} |
# @samp{int8} | @samp{logical} | @samp{string}
#
# The @samp{complex} and @samp{logical} datatypes are not supported in the MATLAB high-level HDF5 interface.
#
# @code{hdf5oct} supports these types in compatibility with @code{h5py}: 
#
# - the @samp{complex} types are mapped to compound HDF5 datatypes with 2 members, 
# @samp{r} and @samp{i}, for real and imaginary part, respectively.
#
# - the @samp{logical} type is mapped to a HDF5 Enum: @code{(FALSE = 0, TRUE = 1)}
#
# @item @option{ChunkSize}
# The value may be either a vector specifying the chunk size,
# or an empty vector [], which means no chunking (this is the default).
# @end table
#
# @seealso{h5write}
# @end deftypefn


function h5create(filename,location,sz,varargin)

# check number and types of arguments
if ((nargin < 3))
    error("h5create: you must supply at least 3 arguments");
endif
if (!ischar(filename))
  error("h5create: 1st argument must be a string holding the hdf5 file name");
endif
create_file = 1;
if isfile(filename)
  # check valid hdf5
  create_file = 0;
endif
if (!ischar(location))
  error("h5create: 2nd argument must be a string holding the dataset location");
endif
if !isreal(sz) || !isvector(sz),
  error("h5create: 3rd argument must be a dimension vector");
endif
sz = sz(:);
# sz elements must be index or Inf
for i=1:length(sz),
  if !(isindex(sz(i)) || sz(i)==inf)
    error("h5create: 3rd argument elements must be valid index values or inf");
  endif
endfor

## check options
[reg, datatype, chunksize, fillvalue] = parseparams (varargin, ...
  'Datatype', 'double',...
  'ChunkSize',[],...
  'FillValue',0);

# check datatype
if !(strcmp(datatype,'double') || ...
  strcmp(datatype,'single') || ...
  strcmp(datatype,'double complex') || ...
  strcmp(datatype,'single complex') || ...
  strcmp(datatype,'uint64') || ...
  strcmp(datatype,'int64') || ...
  strcmp(datatype,'uint32') || ...
  strcmp(datatype,'int32') || ...
  strcmp(datatype,'uint16') || ...
  strcmp(datatype,'int16') || ...
  strcmp(datatype,'uint8') || ...
  strcmp(datatype,'int8') || ...
  strcmp(datatype,'logical') || ...
  strcmp(datatype,'string'))
  error("h5create: invalid 'Datatype'");
endif
# check
if !isindex(sz) && isempty(chunksize)
  error("h5create: chunksize must be defined for inf size datasets");
endif
if !isempty(chunksize)
  if size(chunksize(:))!=size(sz),
    error("h5create: chunksize is different size from 3rd argument");
  endif
  if !isindex(chunksize),
    error("h5create: invalid chunksize");
  endif
  chunksize = chunksize(:);
endif

if size(sz,1)==1, # convert [n] to [1xn]
  sz = [sz; 1];
  if !isempty(chunksize)
    chunksize = [chunksize; 1];
  endif
endif

__h5create__(filename,create_file,location,sz,datatype,chunksize,fillvalue);

# tests for all functions in package

%!shared fname
%! fname =  tempname ();

%!function y = test1(fname,loc,x,type,size)
%!  # create, write, read
%!  h5create(fname,loc,size,'datatype',type);
%!  h5write(fname,loc,x);
%!  y=h5read(fname,loc);
%!endfunction

%!test
%! x = uint32(1:10);
%! y = test1(fname,'/T1/D1',x,'uint32',size(x));
%! assert (x, y);

%!test
%! x = int64(1:10);
%! x = x';
%! y = test1(fname,'/T1/D2',x,'int64',size(x));
%! assert (x, y);

%!test
%! x = reshape(1:150,10,15);
%! y = test1(fname,'/T1/D3',x,'double',size(x));
%! assert (x, y);

%!test
%! x = {"ένα"; "δύο"; "τρία"; "τέσσερα"; "πέντε"};
%! y = test1(fname,'/T1/D4',x,'string',size(x));
%! assert (x, y);

%!test
%! x = "Single string";
%! y = test1(fname,'/T1/D5',x,'string',1);
%! assert (x, y);

%!function y = test2(fname,loc,x,size,chunk,start,count,stride)
%!  # create, write, read
%!  h5create(fname,loc,size,'chunksize',chunk);
%!  h5write(fname,loc,x,start,count,stride);
%!  y=h5read(fname,loc,start,count,stride);
%!endfunction

%!test
%! loc = '/T2/G1/G2/D1';
%! x = reshape(1:80,10,8);
%! size = [Inf 15];
%! chunk = [5 15];
%! start = [11 1];
%! count = [10 8];
%! stride = [1 2];
%! y = test2(fname,loc,x,size,chunk,start,count,stride);
%! assert(x,y);

%!test
%! loc = '/T2/D4';
%! x = reshape(1:60,3,4,5);
%! sz = [3 4 Inf];
%! chunk = [3 4 1];
%! h5create(fname,loc,sz,'chunksize',chunk);
%! h5write(fname,loc,x,[1 1 1],size(x));
%! start = [1 1 1];
%! count = [3 4 3];
%! stride = [1 1 2];
%! y=h5read(fname,loc,start,count,stride);
%! assert(x(:,:,1:2:5),y)

%!function y = test3(fname,loc,attr,x)
%!  # write & read attr
%!  h5writeatt(fname,loc,attr,x);
%!  y=h5readatt(fname,loc,attr);
%!endfunction

%!test
%! loc = '/T2/G1/G2/D1';
%! attr = 'A1';
%! x = 1:10;
%! y = test3(fname,loc,attr,x);
%! assert(x,y);

%!test
%! loc = '/T2/G1/G2';
%! attr = 'A2';
%! x = "Χαρακτηριστικό";
%! y = test3(fname,loc,attr,x);
%! assert(x,y);

%!test
%! loc = '/T2/G1';
%! attr = 'A3';
%! x = {"1ο Χαρακτηριστικό"; "2ο Χαρακτηριστικό"};
%! y = test3(fname,loc,attr,x);
%! assert(x,y);
