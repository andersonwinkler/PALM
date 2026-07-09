/*
 *
 *    Copyright (C) 2012 Tom Mullins
 *    Copyright (C) 2015 Tom Mullins, Thorsten Liebig, Anton Starikov, Stefan Großhauser
 *    Copyright (C) 2008-2013 Andrew Collette
 *    Copyright (C) 2024 George Apostolopoulos
 *
 *    This file is part of hdf5oct.
 *
 *    hdf5oct is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU Lesser General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    hdf5oct is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public License
 *    along with hdf5oct.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 */

#include "hdf5oct.h"

using namespace std;
namespace H5 = HighFive;
namespace h5o = hdf5oct;

// PKG_ADD: autoload("__h5read__","hdf5oct.oct")
// PKG_ADD: autoload("__h5readatt__","hdf5oct.oct")
// PKG_ADD: autoload("__h5write__","hdf5oct.oct")
// PKG_ADD: autoload("__h5writeatt__","hdf5oct.oct")
// PKG_ADD: autoload("__h5create__","hdf5oct.oct")
// PKG_ADD: autoload("__h5delete__","hdf5oct.oct")
// PKG_ADD: autoload("h5info","hdf5oct.oct")

// PKG_DEL: autoload("__h5read__","hdf5oct.oct","remove")
// PKG_DEL: autoload("__h5readatt__","hdf5oct.oct","remove")
// PKG_DEL: autoload("__h5write__","hdf5oct.oct","remove")
// PKG_DEL: autoload("__h5writeatt__","hdf5oct.oct","remove")
// PKG_DEL: autoload("__h5create__","hdf5oct.oct","remove")
// PKG_DEL: autoload("__h5delete__","hdf5oct.oct","remove")
// PKG_DEL: autoload("h5info","hdf5oct.oct","remove")

// __h5create__(fname,create_file,loc,sz,datatype,chunksize,fillvalue)
DEFUN_DLD(__h5create__, args, , "__h5create__: backend for h5create\n\
Users should not use this directly. Use h5create.m instead")
{
    if (args.length() != 7)
        error("__h5create__: wrong # of args");
    string filename = args(0).string_value();
    bool create_file = args(1).bool_value();
    string location = args(2).string_value();
    uint64NDArray size = args(3).uint64_array_value();
    string datatype = args(4).string_value();
    uint64NDArray chunksize = args(5).uint64_array_value();
    int fillvalue = args(6).int_value();
    bool has_unlimited = !chunksize.isempty();
    try
    {
        // open the hdf5 file, create it if it does not exist
        H5::File file(filename, create_file ? H5::File::Create : H5::File::ReadWrite);

        // check location
        if (!h5o::validLocation(location))
            error("h5create: %s", h5o::lastError.c_str());
        if (h5o::locationExists(file, location))
            error("h5create: location '%s' already exists", location.c_str());
        if (!h5o::canCreate(file, location))
            error("h5create: location '%s' cannot be created. "
                  "Check that intermediate nodes are of type Group",
                  location.c_str());

        // create dataspace
        H5::DataSpace fspace = H5::DataSpace::Scalar();
        size_t ndim = size.numel();
        bool is_scalar = (ndim == 1 && size_t(size(0)) == 1) ||
                         (ndim == 2 && size_t(size(0)) == 1 & size_t(size(1)) == 1);
        if (!is_scalar)
        {
            vector<size_t> dims(ndim), maxdims(ndim);
            for (int i = 0; i < ndim; i++)
            {
                bool b = size_t(size(i)) == H5::DataSpace::UNLIMITED;
                dims[ndim - 1 - i] = b ? 0 : size_t(size(i));
                maxdims[ndim - 1 - i] = size_t(size(i));
            }
            fspace = H5::DataSpace(dims, maxdims);
        }

        // Modify dataset creation properties, i.e. enable chunking.
        H5::DataSetCreateProps dscp;
        if (has_unlimited)
        {
            vector<hsize_t> chunk_dims(ndim);
            for (int i = 0; i < ndim; i++)
                chunk_dims[ndim - 1 - i] = chunksize(i);
            dscp.add(H5::Chunking(chunk_dims));
        }
        // Set fill value for the dataset
        // cparms.setFillValue( PredType::NATIVE_INT, &fill_val);

        // create the dataset
        file.createDataSet(location, fspace, h5o::h5type_from_spec(datatype), dscp);
    }
    catch (const H5::Exception &e)
    {
        error("%s", e.what());
    }
    return octave_value_list();
}

// data = h5read(filename,dsname,start,count,stride)
DEFUN_DLD(__h5read__, args, , "__h5read__: backend for h5read\n\
Users should not use this directly. Use h5read.m instead")
{
    if (args.length() != 5)
        error("__h5read__: wrong # of args");
    string filename = args(0).string_value();
    string location = args(1).string_value();
    uint64NDArray start = args(2).uint64_array_value();
    uint64NDArray count = args(3).uint64_array_value();
    uint64NDArray stride = args(4).uint64_array_value();

    try
    {
        H5::File file(filename, H5::File::ReadOnly);

        // check that location is valid, exists and that it is a dataset
        if (!h5o::validLocation(location))
            error("h5read: %s", h5o::lastError.c_str());
        if (!h5o::locationExists(file, location))
            error("h5read: location %s does not exist", location.c_str());
        else if (file.getObjectType(location) != H5::ObjectType::Dataset)
            error("h5read: location '%s' is not a Dataset", location.c_str());

        // Create file dx struct
        h5o::data_exchange dxfile;
        H5::DataSet dset = file.getDataSet(location);
        if (!dxfile.assign(&dset))
            error("h5read: dataset %s: %s", location.c_str(),
                  h5o::lastError.c_str());

        // if requested, select hyperslab
        if (!start.isempty() && !dxfile.selectHyperslab(start, count, stride, false))
            error("h5read: hyperslab selection: %s", h5o::lastError.c_str());

        return dxfile.read();
    }
    catch (const H5::Exception &e)
    {
        error("%s", e.what());
    }
}

// h5write(filename,ds,data,start,count,stride)
DEFUN_DLD(__h5write__, args, , "__h5write__: backend for h5write\n\
Users should not use this directly. Use h5write.m instead")
{
    if (args.length() != 6)
        error("__h5write__: wrong # of args");
    string filename = args(0).string_value();
    string location = args(1).string_value();
    octave_value data = args(2);
    uint64NDArray start = args(3).uint64_array_value();
    uint64NDArray count = args(4).uint64_array_value();
    uint64NDArray stride = args(5).uint64_array_value();

    try
    {
        H5::File file(filename, H5::File::ReadWrite);

        // check that location is valid, exists and that it is a dataset
        if (!h5o::validLocation(location))
            error("h5write: %s", h5o::lastError.c_str());
        if (!h5o::locationExists(file, location))
            error("h5write: location %s does not exist", location.c_str());
        else if (file.getObjectType(location) != H5::ObjectType::Dataset)
            error("h5write: location '%s' is not a Dataset", location.c_str());

        // Create file dx struct
        h5o::data_exchange dxfile;
        H5::DataSet dset = file.getDataSet(location);
        if (!dxfile.assign(&dset))
            error("h5write: dataset %s: %s", location.c_str(),
                  h5o::lastError.c_str());

        // Create octave data dx struct
        h5o::data_exchange dxmem;
        if (!dxmem.assign(data))
            error("h5write: octave data: %s", h5o::lastError.c_str());

        if (dxmem.dtype_spec == "string" && dxfile.dtype_info.size != H5T_VARIABLE)
            error("h5write: export of string data is supported only to variable size HDF5 string datasets");

        if (!start.isempty() && !dxfile.selectHyperslab(start, count, stride, true))
            error("h5write: hyperslab selection: %s", h5o::lastError.c_str());

        if (!dxmem.isCompatible(dxfile))
            error("h5write: incompatible dataset and octave data: %s", h5o::lastError.c_str());

        dxmem.write(dxfile);
    }
    catch (const H5::Exception &e)
    {
        error("%s", e.what());
    }
    return octave_value_list();
}

template <class H5Obj>
octave_value read_attr(const H5Obj &obj, const string &attrname)
{
    if (!obj.hasAttribute(attrname))
        error("h5readatt: attribute %s does not exist", attrname.c_str());
    h5o::data_exchange dx;
    H5::Attribute attr = obj.getAttribute(attrname);
    if (!dx.assign(&attr))
        error("h5readatt: attribute %s: %s", attrname.c_str(),
              h5o::lastError.c_str());
    return dx.read_attribute();
}

// h5readatt(filename,location,attr)
DEFUN_DLD(__h5readatt__, args, , "__h5readatt__: backend for h5readatt\n\
Users should not use this directly. Use h5readatt.m instead")
{
    if (args.length() != 3)
        error("__h5readatt__: wrong # of args");
    string filename = args(0).string_value();
    string location = args(1).string_value();
    string attrname = args(2).string_value();

    try
    {
        H5::File file(filename, H5::File::ReadOnly);

        // check that location is valid & exists
        if (!h5o::validLocation(location))
            error("h5readatt: %s", h5o::lastError.c_str());
        if (!h5o::locationExists(file, location))
            error("h5readatt: location '%s' does not exist", location.c_str());

        switch (file.getObjectType(location))
        {
        case H5::ObjectType::Group:
        {
            H5::Group g = file.getGroup(location);
            return read_attr(g, attrname);
        }
        case H5::ObjectType::Dataset:
        {
            H5::DataSet dset = file.getDataSet(location);
            return read_attr(dset, attrname);
        }
        // case H5::ObjectType::UserDataType:
        default:
            error("h5readatt: location has unsupported object type");
        }
    }
    catch (const H5::Exception &e)
    {
        error("%s", e.what());
    }

    return octave_value_list();
}

// h5writeatt(filename,location,attr,val)
DEFUN_DLD(__h5writeatt__, args, , "__h5writeatt__: backend for h5writeatt\n\
Users should not use this directly. Use h5writeatt.m instead")
{
    if (args.length() != 4)
        error("__h5writeatt__: wrong # of args");
    string filename = args(0).string_value();
    string location = args(1).string_value();
    string attrname = args(2).string_value();
    octave_value data = args(3);

    try
    {
        H5::File file(filename, H5::File::ReadWrite);

        // check that location is valid & exists
        if (!h5o::validLocation(location))
            error("h5writeatt: %s", h5o::lastError.c_str());
        if (!h5o::locationExists(file, location))
            error("h5writeatt: location '%s' does not exist", location.c_str());

        h5o::data_exchange dxmem;
        if (!dxmem.assign(data))
        {
            error("h5writeatt: data value val: %s", h5o::lastError.c_str());
            return octave_value();
        }

        switch (file.getObjectType(location))
        {
        case H5::ObjectType::Group:
        {
            H5::Group g = file.getGroup(location);
            if (!dxmem.write_as_attribute(g, attrname))
                error("h5writeatt: could not write attr: %s", h5o::lastError.c_str());
        }
        break;
        case H5::ObjectType::Dataset:
        {
            H5::DataSet dset = file.getDataSet(location);
            if (!dxmem.write_as_attribute(dset, attrname))
                error("h5writeatt: could not write attr: %s", h5o::lastError.c_str());
        }
        break;
        // case H5::ObjectType::UserDataType:
        default:
            error("h5writeatt: location has unsupported object type");
        }
    }
    catch (const H5::Exception &e)
    {
        error("%s", e.what());
    }

    return octave_value_list();
}

// h5delete(filename,location)
DEFUN_DLD(__h5delete__, args, , "__h5delete__: backend for h5delete\n\
Users should not use this directly. Use h5delete.m instead")
{
    if (args.length() != 2)
        error("__h5delete__: wrong # of args");
    string filename = args(0).string_value();
    string location = args(1).string_value();

    try
    {
        H5::File file(filename, H5::File::ReadWrite);

        if (!h5o::validLocation(location))
            error("h5delete: %s", h5o::lastError.c_str());
        if (location == "/")
            error("h5delete: cannot delete root group '/'");

        if (h5o::locationExists(file, location))
        {
            switch (file.getObjectType(location))
            {
            case H5::ObjectType::Group:
            case H5::ObjectType::Dataset:
            {
                string::size_type pos = location.rfind('/');
                string parent = (pos == 0) ? "/" : location.substr(0, pos);
                string name = location.substr(pos + 1);
                H5::Group g = file.getGroup(parent);
                g.unlink(name);
            }
            break;
            default:
                error("h5delete: location '%s' has unsupported object type",
                      location.c_str());
            }
        }
        else
        {
            string::size_type pos = location.rfind('/');
            string parent = (pos == 0) ? "/" : location.substr(0, pos);
            string attrname = location.substr(pos + 1);
            if (!h5o::locationExists(file, parent))
                error("h5delete: location '%s' does not exist", location.c_str());

            switch (file.getObjectType(parent))
            {
            case H5::ObjectType::Group:
            {
                H5::Group g = file.getGroup(parent);
                if (!g.hasAttribute(attrname))
                    error("h5delete: location '%s' does not exist", location.c_str());
                g.deleteAttribute(attrname);
            }
            break;
            case H5::ObjectType::Dataset:
            {
                H5::DataSet dset = file.getDataSet(parent);
                if (!dset.hasAttribute(attrname))
                    error("h5delete: location '%s' does not exist", location.c_str());
                dset.deleteAttribute(attrname);
            }
            break;
            default:
                error("h5delete: location '%s' does not exist", location.c_str());
            }
        }
    }
    catch (const H5::Exception &e)
    {
        error("%s", e.what());
    }

    return octave_value_list();
}

DEFUN_DLD(h5info, args, argout, "-*- texinfo -*- \n\
@deftypefn {Loadable Function} {@var{info}=} h5info (@var{filename}) \n\
@deftypefnx {Loadable Function} {@var{info}=} h5info (@var{filename}, @var{location}) \n\n\
Return information about the content of a HDF5 file.\n\n\
@code{info = h5info(@var{filename})} returns information about an entire HDF5 file.\n\n\
@code{info = h5info(@var{filename},@var{location})} returns \
information about the specified location in the HDF5 file.\n\n\
The return value, @var{info}, is a structure with detailed information \
on the groups, datasets, and datatypes contained in the file.\n\n\
@seealso{h5disp}\n@end deftypefn")
{
    int nargin = args.length();
    if (nargin < 1 || nargin > 2)
    {
        print_usage();
        return octave_value();
    }
    if (!(args(0).is_string()))
    {
        print_usage();
        return octave_value();
    }
    string filename = args(0).string_value();
    string location("/");
    if (nargin == 2)
    {
        if (!(args(1).is_string()))
        {
            print_usage();
            return octave_value();
        }
        location = args(1).string_value();
    }

    octave_scalar_map info;
    try
    {
        H5::File file(filename, H5::File::ReadOnly);

        // check that location is valid & exists
        if (!h5o::validLocation(location))
            error("h5info: %s", h5o::lastError.c_str());
        if (!h5o::locationExists(file, location))
            error("h5info: location '%s' does not exist", location.c_str());

        switch (file.getObjectType(location))
        {
        case H5::ObjectType::Group:
        {
            h5o::group_info_t I;
            I.assign(file.getGroup(location), location);
            info = I.oct_map();
        }
        break;
        case H5::ObjectType::Dataset:
        {
            h5o::dset_info_t I;
            I.assign(file.getDataSet(location), location);
            info = I.oct_map();
        }
        break;
        default:
            error("h5info: location '%s' has unsupported object type",
                  location.c_str());
        }
    }
    catch (const H5::Exception &e)
    {
        error("%s", e.what());
    }
    return octave_value_list(info);
}

HighFive::DataType hdf5oct::h5type_from_spec(const std::string &dtype_spec)
{
    if (dtype_spec == "double")
        return h5traits<double>::predType();
    else if (dtype_spec == "single")
        return h5traits<float>::predType();
    else if (dtype_spec == "double complex")
        return h5traits<std::complex<double>>::predType();
    else if (dtype_spec == "single complex")
        return h5traits<std::complex<float>>::predType();
    else if (dtype_spec == "uint64")
        return h5traits<uint64_t>::predType();
    else if (dtype_spec == "int64")
        return h5traits<int64_t>::predType();
    else if (dtype_spec == "uint32")
        return h5traits<uint32_t>::predType();
    else if (dtype_spec == "int32")
        return h5traits<int32_t>::predType();
    else if (dtype_spec == "uint16")
        return h5traits<uint16_t>::predType();
    else if (dtype_spec == "int16")
        return h5traits<int16_t>::predType();
    else if (dtype_spec == "uint8")
        return h5traits<uint8_t>::predType();
    else if (dtype_spec == "int8")
        return h5traits<int8_t>::predType();
    else if (dtype_spec == "logical")
        return h5traits<bool>::predType();
    else if (dtype_spec == "string")
        return h5traits<std::string>::predType();
    return HighFive::DataType();
}

string h5_concat_path(const string &path, const string &obj_name)
{
    string p(path);
    if (path.back() != '/')
        p += '/';
    p += obj_name;
    return p;
}

void hdf5oct::dspace_info_t::assign(const H5::DataSpace &ds)
{
    int ndim;
    switch (H5Sget_simple_extent_type(ds.getId()))
    {
    case H5S_SIMPLE:
        extent_type = "Simple";
        ndim = ds.getNumberDimensions();
        size = uint64NDArray(dim_vector(1, ndim));
        maxsize = uint64NDArray(dim_vector(1, ndim));
        {
            auto dims = ds.getDimensions();
            auto maxdims = ds.getMaxDimensions();
            for (octave_idx_type i = 0; i < ndim; i++)
            {
                size(i) = dims[ndim - 1 - i];
                maxsize(i) = maxdims[ndim - i - 1];
            }
        }
        break;
    case H5S_SCALAR:
        extent_type = "Scalar";
        break;
    case H5S_NULL:
        extent_type = "Empty";
        break;
    default:
        break;
    }
}
octave_scalar_map hdf5oct::dspace_info_t::oct_map() const
{
    map<string, octave_value> M;
    M["ExtentType"] = extent_type;
    M["Size"] = size;
    // Return maxsize as a double array to get the "inf" for unlimited dimensions
    RowVector dmaxsize(maxsize.numel());
    for (int i = 0; i < maxsize.numel(); i++)
    {
        dmaxsize(i) = hsize_t(maxsize(i)) == H5::DataSpace::UNLIMITED ? INFINITY : double(maxsize(i));
    }
    M["MaxSize"] = dmaxsize;
    return octave_scalar_map(M);
}
bool hdf5oct::dspace_info_t::needs_chunk() const
{
    if (size.isempty())
        return false;
    int n = size.numel();
    // HDF5 specs: chunking is needed if size!=maxsize OR maxsize==unlimited
    for (int i = 0; i < n; i++)
    {
        if (size(i) != maxsize(i) || (size_t(maxsize(i)) == H5::DataSpace::UNLIMITED))
            return true;
    }
    return false;
}

void hdf5oct::dtype_info_t::assign(const H5::DataType &dt)
{
    string s = dt.string();
    switch (dt.getClass())
    {
    case H5::DataTypeClass::Integer:
        h5class = "int";
        {
            bool isUnsigned = (H5Tget_sign(dt.getId()) == H5T_SGN_NONE);
            signed_status = isUnsigned ? "unsigned" : "signed";
            size = dt.getSize();
            if (size <= 8)
            {
                octave_class = isUnsigned ? "uint" : "int";
                octave_class += to_string(size * 8);
            }
        }
        break;
    case H5::DataTypeClass::Float:
        h5class = "float";
        {
            size = dt.getSize();
            if (size == 4)
                octave_class = "single";
            if (size == 8)
                octave_class = "double";
        }
        break;
    case H5::DataTypeClass::Time:
        h5class = "time";
        break;
    case H5::DataTypeClass::String:
        h5class = "string";
        {
            H5::StringType st = dt.asStringType();
            size = st.isVariableStr() ? H5T_VARIABLE : st.getSize();
            switch (st.getCharacterSet())
            {
            case H5::CharacterSet::Ascii:
                cset = "ascii";
                break;
            case H5::CharacterSet::Utf8:
                cset = "utf8";
                break;
            default:
                break;
            }
            switch (st.getPadding())
            {
            case H5::StringPadding::NullTerminated:
                pading = "nullterm";
                break;
            case H5::StringPadding::NullPadded:
                pading = "nullpad";
                break;
            case H5::StringPadding::SpacePadded:
                pading = "spacepad";
                break;
            default:
                break;
            }
            octave_class = "string";
        }
        break;
    case H5::DataTypeClass::BitField:
        h5class = "bitfield";
        break;
    case H5::DataTypeClass::Opaque:
        h5class = "opaque";
        break;
    case H5::DataTypeClass::Compound:
        h5class = "compound";
        if (dt == h5o::h5traits<std::complex<double>>::predType())
        {
            octave_class = "double complex";
            size = dt.getSize();
        }
        else if (dt == h5o::h5traits<std::complex<float>>::predType())
        {
            octave_class = "single complex";
            size = dt.getSize();
        }
        break;
    case H5::DataTypeClass::Reference:
        h5class = "reference";
        break;
    case H5::DataTypeClass::Enum:
        h5class = "enum";
        if (dt == h5o::h5traits<bool>::predType())
        {
            octave_class = "logical";
            size = dt.getSize();
        }
        break;
    case H5::DataTypeClass::VarLen:
        h5class = "vlen";
        break;
    case H5::DataTypeClass::Array:
        h5class = "array";
        break;
    case H5::DataTypeClass::Invalid:
    default:
        break;
    }
}
octave_scalar_map hdf5oct::dtype_info_t::oct_map() const
{
    map<string, octave_value> M;
    M["Class"] = h5class;
    M["OctaveClass"] = octave_class;
    if (size)
        M["Size"] = size == H5T_VARIABLE ? octave_value("H5T_VARIABLE") : octave_value(size);
    if (!signed_status.empty())
        M["Sign"] = signed_status;
    if (!cset.empty())
        M["charSet"] = cset;
    if (!pading.empty())
        M["Pading"] = pading;
    return octave_scalar_map(M);
}
void hdf5oct::dset_info_t::assign(const H5::DataSet &ds, const string &path)
{
    name = path;
    H5::DataType dt = ds.getDataType();
    dtype_info.assign(dt);
    dspace_info.assign(ds.getSpace());
    // check creation properties
    H5::DataSetCreateProps dscpl = ds.getCreatePropertyList();
    if (dspace_info.isSimple() && H5Pget_layout(dscpl.getId()) == H5D_CHUNKED)
    {
        int ndim = dspace_info.size.numel();
        vector<hsize_t> hdims(ndim);
        if (H5Pget_chunk(dscpl.getId(), ndim, hdims.data()) == ndim)
        {
            chunksize = uint64NDArray(dim_vector(1, ndim));
            for (int i = 0; i < ndim; i++)
                chunksize(i) = hdims[ndim - 1 - i];
        }
    }
    attributes = readAttributes(ds);
}
octave_scalar_map hdf5oct::dset_info_t::oct_map() const
{
    // all map fields mandatory !!
    map<string, octave_value> M;
    M["Name"] = name;
    M["Datatype"] = dtype_info.oct_map();
    M["Dataspace"] = dspace_info.oct_map();
    M["ChunkSize"] = chunksize;
    M["Attributes"] = attributes;
    return octave_scalar_map(M);
}
void hdf5oct::group_info_t::assign(const H5::Group &g, const string &path)
{
    name = path;
    size_t n = g.getNumberObjects();
    vector<string> group_names, ds_names, dt_names, link_names;
    for (hsize_t i = 0; i < n; i++)
    {
        string obj_name = g.getObjectName(i);
        H5::ObjectType type = g.getObjectType(obj_name);

        switch (type)
        {
        case H5::ObjectType::Group:
        {
            group_info_t info;
            info.assign(g.getGroup(obj_name),
                        h5_concat_path(path, obj_name));
            groups.push_back(info);
        }
        break;
        case H5::ObjectType::Dataset:
        {
            dset_info_t info;
            info.assign(g.getDataSet(obj_name),
                        h5_concat_path(path, obj_name));
            datasets.push_back(info);
        }
        break;
        case H5::ObjectType::UserDataType:
        {
            named_dtype_info_t info;
            info.assign(g.getDataType(obj_name),
                        h5_concat_path(path, obj_name));
            datatypes.push_back(info);
        }
        break;
        case H5::ObjectType::Other:
        default:
            break;
        }
    }
    attributes = readAttributes(g);
}
octave_scalar_map hdf5oct::group_info_t::oct_map() const
{
    // all map fields mandatory !!
    map<string, octave_value> M;
    M["Name"] = name;

    octave_map omap;

    vector<string> keys{"Name", "Groups", "Datasets",
                        "Datatypes", "Attributes"};
    omap = octave_map(dim_vector(groups.size(), 1), keys);
    for (int i = 0; i < groups.size(); i++)
        omap.fast_elem_insert(i, groups[i].oct_map());
    M["Groups"] = omap;

    keys = {"Name", "Datatype", "Dataspace", "ChunkSize", "Attributes"};
    omap = octave_map(dim_vector(datasets.size(), 1), keys);
    for (int i = 0; i < datasets.size(); i++)
        omap.fast_elem_insert(i, datasets[i].oct_map());
    M["Datasets"] = omap;

    keys = {"Name"};
    omap = octave_map(dim_vector(datatypes.size(), 1), keys);
    for (int i = 0; i < datatypes.size(); i++)
        omap.fast_elem_insert(i, datatypes[i].oct_map());
    M["Datatypes"] = omap;

    M["Attributes"] = attributes;

    return octave_scalar_map(M);
}

// compatible means:
// - same datatype spec
// - same extent type (Scalar/Simple)
// - same dimensionality
// - same # of elements in each dim
bool hdf5oct::data_exchange::isCompatible(const data_exchange &dx)
{
    if (dtype_spec != dx.dtype_spec)
    {
        lastError = "different datatype specs";
        return false;
    }
    if (dtype.getSize() != dx.dtype.getSize())
    {
        lastError = "different datatype size";
        return false;
    }

    if (dspace_info.extent_type != dx.dspace_info.extent_type)
    {
        lastError = "different dataspaces extent type";
        return false;
    }

    if (dspace_info.isSimple())
    {
        int ndim = dspace.getNumberDimensions();
        if (ndim != dx.dspace.getNumberDimensions())
        {
            lastError = "different dataspace dimensionality";
            return false;
        }
        if (ndim == 1 && (dv(1) != dx.dv(1)))
        {
            lastError = "different number of elements";
            return false;
        }
        if (ndim > 1)
        {
            for (int i = 0; i < ndim; i++)
            {
                if (dv(i) != dx.dv(i))
                {
                    lastError = "different number of elements";
                    return false;
                }
            }
        }
    }
    return true;
}

void hdf5oct::data_exchange::reset()
{
    lastError = "";
    ov = octave_value();
    dset = nullptr;
    attr = nullptr;
    dtype = H5::DataType();
    dtype_info = dtype_info_t();
    dspace = H5::DataSpace::Null();
    dspace_info = dspace_info_t();
    dtype_spec = "";
}

bool hdf5oct::data_exchange::assign(H5::DataSet *ds)
{
    reset();
    dset = ds;
    return assign(ds->getDataType(), ds->getSpace());
}

bool hdf5oct::data_exchange::assign(const H5::Attribute *aattr)
{
    reset();
    attr = aattr;
    return assign(attr->getDataType(), attr->getSpace());
}

bool hdf5oct::data_exchange::assign(const H5::DataType &t, const H5::DataSpace &s)
{
    dtype = t;
    dspace = s;

    // check datatype
    dtype_info.assign(dtype);
    if (dtype_info.h5class.empty())
    {
        lastError = "Invalid HDF5 datatype";
        return false;
    }
    if (dtype_info.octave_class.empty())
    {
        lastError = "Unsupported datatype '";
        lastError += dtype_info.h5class;
        lastError += "'";
        return false;
    }
    dtype_spec = dtype_info.octave_class;

    // check dataspace
    dspace_info.assign(dspace);
    if (dspace_info.extent_type.empty() ||
        dspace_info.isNull())
    {
        lastError = "Empty HDF5 dataspace";
        return false;
    }
    if (dspace_info.isScalar())
        dv = dim_vector(1, 1);
    if (dspace_info.isSimple())
    {
        size_t n = dspace_info.size.numel();
        if (n == 1)
            dv = dim_vector(1, dspace_info.size(0)); // 1D dataset -> [1xN] octave array
        else
        {
            dv.resize(n);
            for (size_t i = 0; i < n; i++)
                dv(i) = dspace_info.size(i);
        }
    }

    return true;
}

bool hdf5oct::data_exchange::selectHyperslab(uint64NDArray start,
                                             uint64NDArray count, uint64NDArray stride, bool tryExtend)
{
    if (dspace_info.isNull())
    {
        lastError = "Hyperslab selection not possible. Dataset is NULL";
        return false;
    }
    if (dspace_info.isScalar())
    {
        lastError = "Hyperslab selection not possible for scalar datasets";
        return false;
    }

    int ndim = dspace_info.size.numel();
    vector<size_t> hstart(ndim), hcount(ndim),
        hstride(ndim),
        fdims = dspace.getDimensions();

    bool need_extend = false;
    for (int i = 0; i < ndim; i++)
    {
        int j = ndim - 1 - i;
        hstart[j] = size_t(start(i)) - 1;
        hcount[j] = count(i);
        hstride[j] = (stride.isempty()) ? 1 : size_t(stride(i));
        size_t end = hstart[j] + 1 + (hcount[j] - 1) * hstride[j];
        if (end > fdims[j])
        {
            if (tryExtend)
            {
                if (end < size_t(dspace_info.maxsize(i)))
                {
                    fdims[j] = end;
                    need_extend = true;
                }
                else
                {
                    lastError = "Hyperslab selection superseeds size of non-extensible dataset";
                    return false;
                }
            }
            else
            {
                lastError = "hyperslab selection beyond the dataset's size";
                return false;
            }
        }
    }
    if (need_extend)
        dset->resize(fdims);
    dspace = dset->getSpace();
    H5::HyperSlab slab(H5::RegularHyperSlab(hstart, hcount, hstride));
    dspace = slab.apply(dspace);
    // adjust dim_vector
    if (ndim == 1)
        dv(1) = count(0); // special case of 1D h5 -> [1xN] octave
    else
        for (int i = 0; i < ndim; i++)
            dv(i) = count(i);
    return true;
}

bool hdf5oct::data_exchange::assign(octave_value v)
{

    reset();

    // find class
    if (v.islogical())
    {
        dtype_spec = "logical";
    }
    else if (v.isinteger())
    {
        if (v.is_uint64_type())
            dtype_spec = "uint64";
        else if (v.is_int64_type())
            dtype_spec = "int64";
        else if (v.is_uint32_type())
            dtype_spec = "uint32";
        else if (v.is_int32_type())
            dtype_spec = "int32";
        else if (v.is_uint16_type())
            dtype_spec = "uint16";
        else if (v.is_int16_type())
            dtype_spec = "int16";
        else if (v.is_uint8_type())
            dtype_spec = "uint8";
        else if (v.is_int8_type())
            dtype_spec = "int8";
    }
    else if (v.isreal())
    {
        if (v.is_double_type())
            dtype_spec = "double";
        else if (v.is_single_type())
            dtype_spec = "single";
    }
    else if (v.iscomplex())
    {
        if (v.is_double_type())
            dtype_spec = "double complex";
        else if (v.is_single_type())
            dtype_spec = "single complex";
    }
    else if (v.iscellstr())
        dtype_spec = "string";

    if (dtype_spec.empty())
    {
        lastError = "Unsupported Octave data of class '";
        lastError += v.class_name();
        lastError += "'";
        return false;
    }
    // get h5 datatype
    dtype = h5type_from_spec(dtype_spec);

    // prepare h5 memory dataspace & octave dimensions
    dspace = H5::DataSpace::Scalar();
    dv = v.dims();
    int ndim = dv.ndims();
    // find # of data elements
    hsize_t nelem = 1;
    for (int i = 0; i < ndim; i++)
        nelem *= dv(i);
    if (nelem > 1) // if not scalar (octave 1x1)
    {
        vector<size_t> dims(ndim);
        for (int i = 0; i < ndim; i++)
            dims[ndim - 1 - i] = dv(i);
        dspace = H5::DataSpace(dims);
    }
    dspace_info.assign(dspace);

    ov = v;

    return true;
}

void hdf5oct::data_exchange::write(const data_exchange &dxfile)
{
    if (dtype_spec == "double")
        write_impl<double>(dxfile);
    else if (dtype_spec == "single")
        write_impl<float>(dxfile);
    else if (dtype_spec == "double complex")
        write_impl<std::complex<double>>(dxfile);
    else if (dtype_spec == "single complex")
        write_impl<std::complex<float>>(dxfile);
    else if (dtype_spec == "uint64")
        write_impl<uint64_t>(dxfile);
    else if (dtype_spec == "int64")
        write_impl<int64_t>(dxfile);
    else if (dtype_spec == "uint32")
        write_impl<uint32_t>(dxfile);
    else if (dtype_spec == "int32")
        write_impl<int32_t>(dxfile);
    else if (dtype_spec == "uint16")
        write_impl<uint16_t>(dxfile);
    else if (dtype_spec == "int16")
        write_impl<int16_t>(dxfile);
    else if (dtype_spec == "uint8")
        write_impl<uint8_t>(dxfile);
    else if (dtype_spec == "int8")
        write_impl<int8_t>(dxfile);
    else if (dtype_spec == "logical")
        write_impl<bool>(dxfile);
    else if (dtype_spec == "string")
        write_string(dxfile);
}

void hdf5oct::data_exchange::write_string(const data_exchange &dxfile)
{
    Array<string> A = ov.cellstr_value();
    octave_idx_type n = A.numel();
    vector<const char *> p(n);
    for (octave_idx_type i = 0; i < n; i++)
        p[i] = A(i).data();
    // dxfile.dset.write(p.data(), dtype, dspace, dxfile.dspace);
    // dxfile.dset.write_raw(p.data(), dtype);
    h5write(*dxfile.dset, p.data(), dtype, dspace, dxfile.dspace);
}

octave_value hdf5oct::data_exchange::read_string()
{
    octave_idx_type n = 1;
    for (octave_idx_type i = 0; i < dv.ndims(); i++)
        n *= dv(i);

    if (dtype_info.size == H5T_VARIABLE)
    {
        vector<char *> p(n);
        Array<string> A(dv);
        H5::DataSpace memspace = from_dim_vector(dv);
        h5read(*dset, p.data(), dtype, memspace, dspace);
        for (octave_idx_type i = 0; i < n; i++)
            A(i) = p[i] ? string(p[i]) : string();
        herr_t ret = H5Dvlen_reclaim(dtype.getId(), memspace.getId(), H5P_DEFAULT, p.data());
        if (ret < 0)
        {
            lastError = "Error in call to H5Dvlen_reclaim";
            return octave_value();
        }
        return (n > 1) ? octave_value(A) : octave_value(A(0));
    }
    else
    {
        octave_idx_type sz = dtype_info.size;
        vector<char> buff(n * sz, '\0');
        H5::DataSpace memspace = from_dim_vector(dv);
        h5read(*dset, buff.data(), dtype, memspace, dspace);
        Array<string> A(dv);
        const char *p = buff.data();
        for (octave_idx_type i = 0; i < n; i++)
        {
            A(i) = string(p, sz);
            p += sz;
        }
        return n > 1 ? octave_value(A) : octave_value(A(0));
    }
}

octave_value hdf5oct::data_exchange::read_attribute()
{
    octave_value ret;
    if (dtype_spec == "double")
        ret = read_attr_impl<double>();
    else if (dtype_spec == "single")
        ret = read_attr_impl<float>();
    else if (dtype_spec == "double complex")
        ret = read_attr_impl<std::complex<double>>();
    else if (dtype_spec == "single complex")
        ret = read_attr_impl<std::complex<float>>();
    else if (dtype_spec == "uint64")
        ret = read_attr_impl<uint64_t>();
    else if (dtype_spec == "int64")
        ret = read_attr_impl<int64_t>();
    else if (dtype_spec == "uint32")
        ret = read_attr_impl<uint32_t>();
    else if (dtype_spec == "int32")
        ret = read_attr_impl<int32_t>();
    else if (dtype_spec == "uint16")
        ret = read_attr_impl<uint16_t>();
    else if (dtype_spec == "int16")
        ret = read_attr_impl<int16_t>();
    else if (dtype_spec == "uint8")
        ret = read_attr_impl<uint8_t>();
    else if (dtype_spec == "int8")
        ret = read_attr_impl<int8_t>();
    else if (dtype_spec == "logical")
        ret = read_attr_impl<bool>();
    else if (dtype_spec == "string")
        ret = read_string_attr();
    return ret;
}

octave_value hdf5oct::data_exchange::read_string_attr()
{
    if (attr == nullptr)
    {
        lastError = "Read requested from null attribute";
        return octave_value();
    }
    octave_idx_type n = 1;
    for (octave_idx_type i = 0; i < dv.ndims(); i++)
        n *= dv(i);

    Array<string> A(dv);
    if (dtype_info.size == H5T_VARIABLE)
    {
        vector<char *> p(n);
        H5::DataSpace memspace = from_dim_vector(dv);
        attr->read(p.data(), dtype);
        for (octave_idx_type i = 0; i < n; i++)
            A(i) = p[i] ? string(p[i]) : string();
        herr_t ret = H5Dvlen_reclaim(dtype.getId(), memspace.getId(), H5P_DEFAULT, p.data());
        if (ret < 0)
        {
            lastError = "Error in call to H5Dvlen_reclaim";
            throw H5::AttributeException("Failed to reclaim HDF5 internal memory");
            return octave_value();
        }
    }
    else
    {
        octave_idx_type sz = dtype_info.size;
        vector<char> buff(n * sz, '\0');
        attr->read(buff.data(), dtype);
        const char *p = buff.data();
        for (octave_idx_type i = 0; i < n; i++)
        {
            A(i) = string(p, sz);
            p += sz;
        }
    }
    return n > 1 ? octave_value(A) : octave_value(A(0));
}

octave_value hdf5oct::data_exchange::read()
{
    octave_value ret;
    if (dtype_spec == "double")
        ret = read_impl<double>();
    else if (dtype_spec == "single")
        ret = read_impl<float>();
    else if (dtype_spec == "double complex")
        ret = read_impl<std::complex<double>>();
    else if (dtype_spec == "single complex")
        ret = read_impl<std::complex<float>>();
    else if (dtype_spec == "uint64")
        ret = read_impl<uint64_t>();
    else if (dtype_spec == "int64")
        ret = read_impl<int64_t>();
    else if (dtype_spec == "uint32")
        ret = read_impl<uint32_t>();
    else if (dtype_spec == "int32")
        ret = read_impl<int32_t>();
    else if (dtype_spec == "uint16")
        ret = read_impl<uint16_t>();
    else if (dtype_spec == "int16")
        ret = read_impl<int16_t>();
    else if (dtype_spec == "uint8")
        ret = read_impl<uint8_t>();
    else if (dtype_spec == "int8")
        ret = read_impl<int8_t>();
    else if (dtype_spec == "logical")
        ret = read_impl<bool>();
    else if (dtype_spec == "string")
        ret = read_string();
    return ret;
}

H5::DataSpace hdf5oct::data_exchange::from_dim_vector(const dim_vector &dv)
{
    int ndim = dv.ndims();
    if (ndim == 2 && dv(0) == 1 && dv(1) == 1)
        return H5::DataSpace::Scalar();
    if (ndim == 2 && dv(0) == 1)
    { // [1xN] octave <-> 1D dataspace
        vector<size_t> dims(1);
        dims[0] = dv(1);
        return H5::DataSpace(dims);
    }
    vector<size_t> dims(ndim);
    for (int i = 0; i < ndim; i++)
        dims[ndim - 1 - i] = dv(i);
    return H5::DataSpace(dims);
}

bool hdf5oct::data_exchange::write_as_attribute(H5::Attribute &attr)
{
    if (dtype_spec == "double")
        write_attr_impl<double>(attr);
    else if (dtype_spec == "single")
        write_attr_impl<float>(attr);
    else if (dtype_spec == "double complex")
        write_attr_impl<std::complex<double>>(attr);
    else if (dtype_spec == "single complex")
        write_attr_impl<std::complex<float>>(attr);
    else if (dtype_spec == "uint64")
        write_attr_impl<uint64_t>(attr);
    else if (dtype_spec == "int64")
        write_attr_impl<int64_t>(attr);
    else if (dtype_spec == "uint32")
        write_attr_impl<uint32_t>(attr);
    else if (dtype_spec == "int32")
        write_attr_impl<int32_t>(attr);
    else if (dtype_spec == "uint16")
        write_attr_impl<uint16_t>(attr);
    else if (dtype_spec == "int16")
        write_attr_impl<int16_t>(attr);
    else if (dtype_spec == "uint8")
        write_attr_impl<uint8_t>(attr);
    else if (dtype_spec == "int8")
        write_attr_impl<int8_t>(attr);
    else if (dtype_spec == "logical")
        write_attr_impl<bool>(attr);
    else if (dtype_spec == "string")
        write_string_attr(attr);
    return true;
}

void hdf5oct::data_exchange::write_string_attr(H5::Attribute &attr)
{
    Array<string> A = ov.cellstr_value();
    octave_idx_type n = A.numel();
    vector<const char *> p(n);
    for (octave_idx_type i = 0; i < n; i++)
        p[i] = A(i).data();
    attr.write_raw(p.data(), dtype);
}

bool hdf5oct::locationExists(const H5::File &f, const std::string &loc)
{
    // check for intermediate groups
    string::size_type pos = loc.find_first_of('/'), last_pos = 0;
    while (pos != string::npos)
    {
        if (pos != 0)
        { // skip checking if root exists
            string iloc = loc.substr(0, pos);
            if (!f.exist(iloc))
                return false;
        }
        last_pos = pos + 1;
        pos = loc.find_first_of('/', last_pos);
    }
    if (last_pos < loc.size())
        return f.exist(loc);
    return true;
}

bool hdf5oct::canCreate(const H5::File &f, const std::string &loc)
{
    // check for intermediate groups
    string::size_type pos = loc.find_first_of('/'), last_pos = 0;
    while (pos != string::npos)
    {
        if (pos != 0)
        { // skip checking is root exists
            string iloc = loc.substr(0, pos);
            if (!f.exist(iloc))
                return true; // ok, will be created
            if (f.getObjectType(iloc) != H5::ObjectType::Group)
                return false; // error, intermediate node is not a group
        }
        last_pos = pos + 1;
        pos = loc.find_first_of('/', last_pos);
    }
    return true;
}

bool hdf5oct::validLocation(const std::string &loc)
{
    if (loc.find_first_of('/') != 0)
    {
        lastError = "dataset location should be specified "
                    "as a full pathname with leading '/'";
        return false;
    }
    return true;
}
