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
# @deftypefn {Function File} {@var{att_val}=} h5readatt (@var{filename}, @var{location}, @var{attr})
#
# Read a HDF5 attribute.
#
# Retrieves the value of the specified attribute named @var{attr}
# from the specified location @var{location} in the HDF5 file @var{filename}.
#
# Input arguments:
#
# @table @asis
# @item @var{filename}
# Filename of an existing HDF5 file, specified as a string
# @item @var{location}
# The path to an existing node in the HDF5 file, which can be either a Group,
# a Dataset or a Named DataType. The root group, "/", is also valid.
# @item @var{attr}
# Name of attribute, specified as a string scalar.
# @item @var{att_val}
# Value of the attribute read from the file.
# @end table
#
# @seealso{h5writeatt}
# @end deftypefn
#

function attrval = h5readatt(filename,location,attr)

# check number and types of arguments
if nargin != 3,
    print_usage();
endif
if (!ischar(filename))
  error("h5readatt: 1st argument must be a string holding the hdf5 file name");
endif
if (!isfile(filename))
  error("h5readatt: filename does not exist");
endif
if (!ischar(location))
  error("h5readatt: 2nd argument must be a string holding the HDF5 object location");
endif
if (!ischar(attr))
  error("h5readatt: 3rd argument must be a string holding the attribute name");
endif

attrval = __h5readatt__(filename,location,attr);
