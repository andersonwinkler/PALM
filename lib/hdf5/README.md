hdf5oct - a HDF5 wrapper for GNU Octave
=======================================

This is a [GNU Octave](https://octave.org) package for data serialization to/from HDF5 files. 
It provides an interface mostly compatible to MATLAB's **"High-Level Functions for HDF5 files"**.

The following functions are implemented:
```
- h5create
- h5write
- h5writeatt
- h5read
- h5readatt
- h5info
- h5disp
- h5load
- h5delete
```

The functions `h5load` (load entire file or group) and `h5delete` (delete dataset, group, or attribute) are not supported in MATLAB.

`hdf5oct` can be used to export/import multidimensional array data of class

    'double','single','double complex','single complex',
    'uint64','int64', ... 'uint8', 'int8', 
    'logical', 
    'string'

The `complex` and `logical` datatypes are not supported in the MATLAB high-level interface.

`hdf5oct` supports these types in compatibility with [`h5py`](https://www.h5py.org): 

- the `complex` types are mapped to compound HDF5 datatypes with 2 members, 
`r` and `i`, for real and imaginary part, respectively.

- the `logical` type is mapped to a HDF5 Enum: `{FALSE = 0, TRUE = 1}`

# Getting started

The following short OCTAVE code snippets show how to use the package:

```matlab
>> pkg load hdf5oct % load the package
>> data = uint32(reshape(1:10,2,5)) % create data in workspace
data =

   1   3   5   7   9
   2   4   6   8  10

% create a HDF5 file with a dataset named 'D1' of matching size and datatype
>> h5create('test.h5','/D1',size(data),'datatype','uint32')
% store the data to the dataset
>> h5write('test.h5','/D1',data)
% read back the values
>> h5read('test.h5','/D1')
ans =

   1   3   5   7   9
   2   4   6   8  10

% read the 3rd & 6th columns as a [2x2] matrix 
>> h5read('test.h5','/D1',[1 3],[2 2],[1 2])
ans =

   5   9
   6  10
```
Strings are always UTF8 encoded. A single string is written to a scalar dataset (size=1). Multiple strings must be passed as a cell array:

```matlab
>> oneliner = "This is a single string";
>> h5create('test.h5','/D2',1,'datatype','string') % scalar string dataset
>> h5write('test.h5','/D2',oneliner)
>> h5read('test.h5','/D2')
ans = This is a single string
>> str = {"one", "δύο", "три", "neljä"};
>> h5create('test.h5','/D3',size(str),'datatype','string')
>> h5write('test.h5','/D3',str)
>> h5read('test.h5','/D3')
ans =
{
  [1,1] = one
  [1,2] = δύο
  [1,3] = три
  [1,4] = neljä
}
```
The structure of the HDF5 file can be viewed with `h5disp`. `h5info` can also be used for more detail.

```matlab
>> h5disp('test.h5')
HDF5 test.h5
Group '/'
  Dataset '/D1'
    Extent: Simple
    Size: 2x5
    MaxSize: 2x5
    Datatype
      Class: 'int'
      OctaveClass: 'uint32'
      Size: 4
      Sign: unsigned
  Dataset '/D2'
    Extent: Scalar
    Datatype
      Class: 'string'
      OctaveClass: 'string'
      Size: H5T_VARIABLE
      charSet: utf8
      Pading: nullterm
  Dataset '/D3'
    Extent: Simple
    Size: 1x4
    MaxSize: 1x4
    Datatype
      Class: 'string'
      OctaveClass: 'string'
      Size: H5T_VARIABLE
      charSet: utf8
      Pading: nullterm
```

Datasets, groups, and attributes can be removed with `h5delete`:

```matlab
>> h5delete('test.h5', '/D1')              % delete a dataset
>> h5delete('test.h5', '/group_1')         % delete a group and its contents
>> h5delete('test.h5', '/D1/units')        % delete an attribute
```

# Array storage layout convention

In HDF5, arrays are stored in C-style, [row-major order](https://en.wikipedia.org/wiki/Row-_and_column-major_order). On the other hand, OCTAVE and MATLAB employ fortran-style, column-major storage order.

To avoid array transposition operations during data serialization, MATLAB employs the following convention, which is also followed by `hdf5oct`:

|      Storage Type      |         MATLAB/OCTAVE array size         |          HDF5 DataSpace dimensions           |
| :--------------------: | :--------------------------------------: | :------------------------------------------: |
|         Matrix         |              $[N \times M]$              |                $[M \times N]$                |
| Multidimensional Array | $[N_1 \times N_2 \times ... \times N_m]$ | $[N_m \times N_{m-1} \times ... \times N_1]$ |
|       Row Vector       |              $[1 \times N]$              |           $[N \times 1]$ or $[N]$            |
|     Column Vector      |              $[N \times 1]$              |                $[1 \times N]$                |

In this manner, the data is copied "as-is" from memory to disk, minimizing overhead and memory allocations.

A MATLAB user or an OCTAVE user employing `hdf5oct` to export/import data to/from HDF5 files will not notice any difference. However, when a MATLAB- or `hdf5oct`-generated file is opened by another application, or vice-versa, the arrays will appear transposed. 

# Installation

To install the latest package release run the following in OCTAVE

```matlab
    pkg install -forge hdf5oct
```

or, for the latest development snapshot

```matlab
    pkg install "https://github.com/gapost/hdf5oct/archive/master.zip"
```

After successful installation, test the package with

```matlab
    pkg test hdf5oct
```

This performs a number of basic tests on all functions in the package.

# TODO 

- support compression flags for h5create

- h5read: implement MATLAB compatible mapping to OCTAVE of the remaining HDF5 datatypes: `Bitfield, Opaque, Reference, Enum, Compound, Array`

- write more comprehensive tests.

## License

© 2012, Tom Mullins \
© 2015, Tom Mullins, Anton Starikov, Thorsten Liebig, Stefan Großhauser \
© 2008-2013, Andrew Collette \
© 2024-2025, George Apostolopoulos

[HighFive](https://github.com/BlueBrain/HighFive) is used internally as a C++ interface to the HDF5 library.

Released under the [LGPLv3+](COPYING) and [Boost Software License v1.0](src/third_party/HighFive-2.9.0/LICENSE).
