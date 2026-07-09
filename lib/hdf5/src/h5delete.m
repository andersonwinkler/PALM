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
# @deftypefn {Function File} { } h5delete (@var{filename}, @var{location})
#
# Delete a dataset, group, or attribute from a HDF5 file.
#
# @code{h5delete (@var{filename}, @var{location})} removes the object
# specified by @var{location} from the HDF5 file @var{filename}.
#
# If @var{location} is the full path to an existing dataset or group,
# that object is deleted. Deleting a group also deletes all of its contents.
#
# If @var{location} does not refer to an existing dataset or group, it is
# interpreted as @code{parent/attribute}, where @var{attribute} is the name
# of an attribute attached to the HDF5 object at @var{parent}.
# For example, @code{h5delete ('data.h5', '/my_dataset/units')} deletes the
# attribute @code{units} from the dataset @code{/my_dataset}.
#
# Input arguments:
#
# @table @asis
# @item @var{filename}
# Filename of an existing HDF5 file, specified as a string.
# @item @var{location}
# Full path to the dataset, group, or attribute to be deleted,
# specified as a string with a leading @code{'/'}.
# @end table
#
# @seealso{h5create, h5writeatt}
# @end deftypefn
#

function h5delete(filename,location)

# check number and types of arguments
if nargin != 2,
    print_usage();
endif
if (!ischar(filename))
  error("h5delete: 1st argument must be a string holding the hdf5 file name");
endif
if (!isfile(filename))
  error("h5delete: filename does not exist");
endif
if (!ischar(location))
  error("h5delete: 2nd argument must be a string holding the location to delete");
endif

__h5delete__(filename,location);

endfunction

# tests for h5delete

%!shared fname
%! fname = tempname ();

%!test
%! h5create(fname, "/D1", [2 2], "datatype", "double");
%! h5write(fname, "/D1", ones(2, 2));
%! h5delete(fname, "/D1");
%! fail("h5read(fname, '/D1')");

%!test
%! h5create(fname, "/G1/D1", [1 2], "datatype", "int32");
%! h5write(fname, "/G1/D1", int32([1 2]));
%! h5delete(fname, "/G1");
%! fail("h5info(fname, '/G1')");

%!test
%! h5create(fname, "/D2", [1 1], "datatype", "double");
%! h5writeatt(fname, "/D2", "units", "meters");
%! h5delete(fname, "/D2/units");
%! fail("h5readatt(fname, '/D2', 'units')");

%!test
%! h5create(fname, "/E1", [1 1], "datatype", "double");
%! fail("h5delete(fname, '/')");

%!test
%! h5create(fname, "/E2", [1 1], "datatype", "double");
%! fail("h5delete(fname, '/nonexistent_object')");

%!error h5delete (fname)
%!error h5delete ("nonexistingfile.h5", "/foo")
