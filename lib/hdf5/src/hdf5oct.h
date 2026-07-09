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

#ifndef _hdf5oct_h_
#define _hdf5oct_h_

#include <octave/oct.h>

// #if defined (HAVE_HDF5) && defined (HAVE_HDF5_18)
#include <highfive/highfive.hpp>

namespace hdf5oct
{

    static std::string lastError;

    // Helper type traits template
    template <typename T>
    struct h5traits;

    // specialized for all supported types
    template <>
    struct h5traits<bool>
    {
        typedef boolNDArray OctaveArray;
        static OctaveArray toOctaveArray(octave_value v) { return v.bool_array_value(); }
        static const HighFive::DataType predType() { return HighFive::create_datatype<bool>(); }
    };
    template <>
    struct h5traits<std::complex<double>>
    {
        typedef ComplexNDArray OctaveArray;
        static OctaveArray toOctaveArray(octave_value v) { return v.complex_array_value(); }
        static const HighFive::DataType predType() { return HighFive::AtomicType<std::complex<double>>(); }
    };
    template <>
    struct h5traits<std::complex<float>>
    {
        typedef FloatComplexNDArray OctaveArray;
        static OctaveArray toOctaveArray(octave_value v) { return v.float_complex_array_value(); }
        static const HighFive::DataType predType() { return HighFive::AtomicType<std::complex<float>>(); }
    };
    template <>
    struct h5traits<double>
    {
        typedef NDArray OctaveArray;
        static OctaveArray toOctaveArray(octave_value v) { return v.array_value(); }
        static const HighFive::DataType predType() { return HighFive::AtomicType<double>(); }
    };
    template <>
    struct h5traits<float>
    {
        typedef FloatNDArray OctaveArray;
        static OctaveArray toOctaveArray(octave_value v) { return v.float_array_value(); }
        static const HighFive::DataType predType() { return HighFive::AtomicType<float>(); }
    };
    template <>
    struct h5traits<uint64_t>
    {
        typedef uint64NDArray OctaveArray;
        static OctaveArray toOctaveArray(octave_value v) { return v.uint64_array_value(); }
        static const HighFive::DataType predType() { return HighFive::AtomicType<uint64_t>(); }
    };
    template <>
    struct h5traits<uint32_t>
    {
        typedef uint32NDArray OctaveArray;
        static OctaveArray toOctaveArray(octave_value v) { return v.uint32_array_value(); }
        static const HighFive::DataType predType() { return HighFive::AtomicType<uint32_t>(); }
    };
    template <>
    struct h5traits<uint16_t>
    {
        typedef uint16NDArray OctaveArray;
        static OctaveArray toOctaveArray(octave_value v) { return v.uint16_array_value(); }
        static const HighFive::DataType predType() { return HighFive::AtomicType<uint16_t>(); }
    };
    template <>
    struct h5traits<uint8_t>
    {
        typedef uint8NDArray OctaveArray;
        static OctaveArray toOctaveArray(octave_value v) { return v.uint8_array_value(); }
        static const HighFive::DataType predType() { return HighFive::AtomicType<uint8_t>(); }
    };
    template <>
    struct h5traits<int64_t>
    {
        typedef int64NDArray OctaveArray;
        static OctaveArray toOctaveArray(octave_value v) { return v.int64_array_value(); }
        static const HighFive::DataType predType() { return HighFive::AtomicType<int64_t>(); }
    };
    template <>
    struct h5traits<int32_t>
    {
        typedef int32NDArray OctaveArray;
        static OctaveArray toOctaveArray(octave_value v) { return v.int32_array_value(); }
        static const HighFive::DataType predType() { return HighFive::AtomicType<int32_t>(); }
    };
    template <>
    struct h5traits<int16_t>
    {
        typedef int16NDArray OctaveArray;
        static OctaveArray toOctaveArray(octave_value v) { return v.int16_array_value(); }
        static const HighFive::DataType predType() { return HighFive::AtomicType<int16_t>(); }
    };
    template <>
    struct h5traits<int8_t>
    {
        typedef int8NDArray OctaveArray;
        static OctaveArray toOctaveArray(octave_value v) { return v.int8_array_value(); }
        static const HighFive::DataType predType() { return HighFive::AtomicType<int8_t>(); }
    };
    template <>
    struct h5traits<std::string>
    {
        typedef Array<std::string> OctaveArray;
        static OctaveArray toOctaveArray(octave_value v) { return v.cellstr_value(); }
        static const HighFive::DataType predType() { return HighFive::VariableLengthStringType(HighFive::CharacterSet::Utf8); }
    };

    // Translate an octave-like type spec, e.g. 'uint32', to H5 datatype
    HighFive::DataType h5type_from_spec(const std::string &dtype_spec);

    /**
     * @brief Check if a location exists in a HDF5 file.
     *
     * If the location path contains intermediate groups, these are also checked.
     *
     * @param f The HDF5 file object
     * @param loc The location path
     * @return true If the specified location exists
     * @return false If the location or any intermediate groups do not exist
     */
    bool locationExists(const HighFive::File &f, const std::string &loc);

    /**
     * @brief Check if a location can be created in a HDF5 file.
     *
     * If the location path contains intermediate objects, the function
     * checks that either they are of type Group or that they do not exist
     * and thus will be also created.
     *
     * @param f The HDF5 file object
     * @param loc The location path
     * @return true If the specified location can be created
     * @return false Otherwise
     */
    bool canCreate(const HighFive::File &f, const std::string &loc);

    bool validLocation(const std::string &loc);

    // structures with info on H5 objects (DataSpace,DataType,DataSet,Group)
    // oct_map() function returns this info as a (key,value) map
    // for reporting back to Octave in h5info

    // DataSpace
    struct dspace_info_t
    {
        std::string extent_type;
        uint64NDArray size;
        uint64NDArray maxsize;
        void assign(const HighFive::DataSpace &dspace);
        octave_scalar_map oct_map() const;
        bool needs_chunk() const;
        bool isNull() const { return extent_type == "Empty"; }
        bool isScalar() const { return extent_type == "Scalar"; }
        bool isSimple() const { return extent_type == "Simple"; }
    };
    // DataType
    struct dtype_info_t
    {
        std::string h5class;
        std::string octave_class;
        size_t size;
        std::string signed_status;
        std::string cset;
        std::string pading;
        void assign(const HighFive::DataType &dtype);
        octave_scalar_map oct_map() const;
    };
    // DataSet
    struct dset_info_t
    {
        std::string name;
        dtype_info_t dtype_info;
        dspace_info_t dspace_info;
        uint64NDArray chunksize;
        octave_value fillValue;
        std::map<std::string, octave_value> attributes;
        void assign(const HighFive::DataSet &ds, const std::string &path);
        octave_scalar_map oct_map() const;
    };
    // Named DataType (not fully implemented)
    struct named_dtype_info_t : public dtype_info_t
    {
        std::string name;
        void assign(const HighFive::DataType &dtype, const std::string &path)
        {
            name = path;
            dtype_info_t::assign(dtype);
        }
        octave_scalar_map oct_map() const
        {
            // octave_scalar_map m = dtype_info_t::oct_map(); // this needs fixing
            octave_scalar_map m;
            m.assign("Name", name);
            return m;
        }
    };
    // Group
    struct group_info_t
    {
        std::string name;
        std::vector<group_info_t> groups;
        std::vector<dset_info_t> datasets;
        std::vector<named_dtype_info_t> datatypes;
        std::map<std::string, octave_value> attributes;
        void assign(const HighFive::Group &g, const std::string &path);
        octave_scalar_map oct_map() const;
    };

    /**
     * @brief The data_exchage structure facilitates IO operations between H5 & Octave
     *
     * To write an octave value to a dataset:
     *
     * @code {.cpp}
     * octave_value ov;
     * DataSet dset(...);
     * data_exchange dxmem, dxfile;
     * if (!dxmem.assign(ov)) error ...
     * if (!dxfile.assign(&dset)) error ...
     * if (!dxmem.isCompatible(dxfile)) error ...
     * dxmem.write(dxfile);
     * @endcode
     *
     * To read a value:
     *
     * @code {.cpp}
     * DataSet dset(...);
     * data_exchange dxfile;
     * if (!dxfile.assign(&dset)) error ...
     * if (!dxfile.selectHyperslab(...)) error ...
     * octave_value ov = dxfile.read();
     * @endcode
     *
     *
     */
    struct data_exchange
    {
        octave_value ov;
        HighFive::DataSet *dset{nullptr};
        const HighFive::Attribute *attr{nullptr};
        HighFive::DataType dtype;
        dtype_info_t dtype_info;
        HighFive::DataSpace dspace{HighFive::DataSpace::Null()};
        dspace_info_t dspace_info;
        std::string dtype_spec;
        dim_vector dv;

        bool assign(octave_value v);
        bool assign(HighFive::DataSet *ds);
        bool assign(const HighFive::Attribute *attr);
        bool isCompatible(const data_exchange &dx);
        bool isValid() const
        {
            return !(dtype_spec.empty() || dtype_spec == "unsupported");
        }
        bool isScalar() const
        {
            return dv.ndims() == 2 && dv(0) == 1 && dv(1) == 1;
        }
        bool selectHyperslab(uint64NDArray start, uint64NDArray count, uint64NDArray stride, bool tryResize);

        octave_value read();
        octave_value read_attribute();
        void write(const data_exchange &dxfile);

        template <class Derivate>
        bool write_as_attribute(Derivate &obj, const std::string &name)
        {
            if (obj.hasAttribute(name))
                obj.deleteAttribute(name);
            HighFive::Attribute attr = obj.createAttribute(name, dspace, dtype);
            return write_as_attribute(attr);
        }

    private:
        void reset();

        static HighFive::DataSpace from_dim_vector(const dim_vector &dv);
        static void h5read(const HighFive::DataSet &dset,
                           void *data,
                           const HighFive::DataType &mem_type, const HighFive::DataSpace &mem_space,
                           const HighFive::DataSpace &file_space,
                           const HighFive::DataTransferProps &xfer_props = HighFive::DataTransferProps())
        {
            HighFive::detail::h5d_read(dset.getId(),
                                       mem_type.getId(),
                                       mem_space.getId(),
                                       file_space.getId(),
                                       xfer_props.getId(),
                                       data);
        }
        static void h5write(const HighFive::DataSet &dset,
                            const void *data,
                            const HighFive::DataType &mem_type, const HighFive::DataSpace &mem_space,
                            const HighFive::DataSpace &file_space,
                            const HighFive::DataTransferProps &xfer_props = HighFive::DataTransferProps())
        {
            HighFive::detail::h5d_write(dset.getId(),
                                        mem_type.getId(),
                                        mem_space.getId(),
                                        file_space.getId(),
                                        xfer_props.getId(),
                                        data);
        }

        bool assign(const HighFive::DataType &t, const HighFive::DataSpace &s);

        template <class T>
        octave_value read_impl()
        {
            typename h5traits<T>::OctaveArray A(dv);
            h5read(*dset, A.fortran_vec(), h5traits<T>::predType(), from_dim_vector(dv), dspace);
            return octave_value(A);
        }
        template <class T>
        octave_value read_attr_impl()
        {
            typename h5traits<T>::OctaveArray A(dv);
            attr->read(A.fortran_vec(), dtype);
            return octave_value(A);
        }
        octave_value read_string();
        octave_value read_string_attr();

        template <typename T>
        void write_impl(const data_exchange &dxfile)
        {
            auto A = h5traits<T>::toOctaveArray(ov);
            h5write(*dxfile.dset, A.fortran_vec(), h5traits<T>::predType(), dspace, dxfile.dspace);
        }
        template <typename T>
        void write_attr_impl(HighFive::Attribute &att)
        {
            auto A = h5traits<T>::toOctaveArray(ov);
            att.write_raw(A.fortran_vec(), h5traits<T>::predType());
        }
        void write_string(const data_exchange &dxfile);
        void write_string_attr(HighFive::Attribute &attr);
        bool write_as_attribute(HighFive::Attribute &attr);
    };

    template <class H5Obj>
    std::map<std::string, octave_value> readAttributes(const H5Obj &obj)
    {
        std::map<std::string, octave_value> M;
        for (auto &attr_name : obj.listAttributeNames())
        {
            hdf5oct::data_exchange dx;
            octave_value ov;
            HighFive::Attribute attr = obj.getAttribute(attr_name);
            if (dx.assign(&attr))
                ov = dx.read_attribute();
            M[attr_name] = ov;
        }
        return M;
    }

} // namespace hdf5oct

#endif
