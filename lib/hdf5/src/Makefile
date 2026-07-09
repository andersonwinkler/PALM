
OCTAVE ?= octave
MKOCTFILE ?= mkoctfile

OCTAVE_CMD = ${OCTAVE} -H -q -W -f --eval

H5FLAGS=$(shell ${OCTAVE_CMD} 's="";if __octave_config_info__().build_features.HDF5==1,s="-DHAVE_HDF5";end;printf(s)') \
		$(shell ${OCTAVE_CMD} 's="";if __octave_config_info__().build_features.HDF5_18==1,s="-DHAVE_HDF5_18";end;printf(s)') \
		$(shell ${OCTAVE_CMD} 's="";if __octave_config_info__().build_features.HDF5_H==1,s="-DHAVE_HDF5_H";end;printf(s)') \
		$(shell ${OCTAVE_CMD} 's="";if __octave_config_info__().build_features.HDF5_INT2FLOAT_CONVERSIONS==1,s="-DHAVE_HDF5_INT2FLOAT_CONVERSIONS";end;printf(s)') \
		$(shell ${OCTAVE_CMD} 's="";if __octave_config_info__().build_features.HDF5_UTF8==1,s="-DHAVE_HDF5_UTF8";end;printf(s)')

INCLUDES = $(shell ${OCTAVE_CMD} "printf ('%s', __octave_config_info__().build_environment.HDF5_CPPFLAGS)") \
           -I./third_party/HighFive-2.9.0/include

LIBS = $(shell ${OCTAVE_CMD} "printf ('%s %s', __octave_config_info__().build_environment.HDF5_LDFLAGS, __octave_config_info__().build_environment.HDF5_LIBS)")

TEST1 = $(findstring HAVE_HDF5,$(H5FLAGS))
TEST2 = $(findstring HAVE_HDF5_18,$(H5FLAGS))

.PHONY: clean 

all: hdf5oct.cc hdf5oct.h test_hdf5
	$(MKOCTFILE) -v $(H5FLAGS) $(INCLUDES) $(LIBS) hdf5oct.cc 

test_hdf5:
	@if [ "$(TEST1)" != "HAVE_HDF5" ]; then \
		echo "hdf5oct cannot be installed because octave has not been built with HDF5 support"; \
		false; \
	fi
	@if [ "$(TEST2)" != "HAVE_HDF5_18" ]; then \
		echo "hdf5oct cannot be installed because octave has been built with HDF5 < 1.8.0"; \
		false; \
	fi

clean:
	rm -f hdf5oct.oct 

