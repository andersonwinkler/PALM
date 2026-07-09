#!/usr/bin/env make -f
# hdf5oct build for PALM
#
# Builds hdf5oct.oct — Octave-only HDF5 backend for h5read, h5write, etc.
# MATLAB uses built-in h5* functions and needs no compilation from this tree.
#
# Targets:
#   make octave   - build hdf5oct.oct
#   make matlab   - print note (no build needed)
#   make help     - show usage
#
# Prerequisites:
#   - octave and mkoctfile in PATH
#   - Octave built with HDF5 >= 1.8, OR HDF5_PREFIX set to a system install

include Makefile.var

.PHONY: all octave matlab help clean distclean distclean-all check-tools test_hdf5

help:
	@echo "Build hdf5oct for PALM (Octave HDF5 I/O)"
	@echo ""
	@echo "Targets:"
	@echo "  make octave        Build hdf5oct.oct for Octave"
	@echo "  make matlab        No build (MATLAB has built-in h5read/h5write)"
	@echo "  make               Same as make octave"
	@echo "  make check-tools   Verify octave and mkoctfile are available"
	@echo "  make clean         Remove intermediate object files"
	@echo "  make distclean     Remove hdf5oct.oct"
	@echo ""
	@echo "Variables:"
	@echo "  OCTAVE=path        Octave executable (default: octave)"
	@echo "  MKOCTFILE=path     mkoctfile executable (default: mkoctfile)"
	@echo "  HDF5_PREFIX=path   Optional system HDF5 if Octave lacks HDF5"
	@echo "  V=1                Verbose mkoctfile output"
	@echo ""
	@echo "Prerequisite: Octave must be built with HDF5 support, unless HDF5_PREFIX"
	@echo "is set. Check with:"
	@echo "  octave --eval \"disp(__octave_config_info__().build_features.HDF5)\""
	@echo ""
	@echo "Alternative install: pkg install -forge hdf5oct"

check-tools:
	@if ! command -v $(OCTAVE) >/dev/null 2>&1; then \
		echo "Error: $(OCTAVE) not found in PATH."; exit 1; \
	fi
	@if ! command -v $(MKOCTFILE) >/dev/null 2>&1; then \
		echo "Error: $(MKOCTFILE) not found in PATH."; exit 1; \
	fi

test_hdf5:
ifndef HDF5_PREFIX
	@if [ "$(TEST_HDF5)" != "HAVE_HDF5" ]; then \
		echo "Error: Octave was not built with HDF5 support."; \
		echo "  Reinstall Octave with HDF5, set HDF5_PREFIX, or use:"; \
		echo "  pkg install -forge hdf5oct"; \
		exit 1; \
	fi
	@if [ "$(TEST_HDF5_18)" != "HAVE_HDF5_18" ]; then \
		echo "Error: Octave was built with HDF5 < 1.8.0."; \
		exit 1; \
	fi
endif

all: octave

octave: check-tools test_hdf5 $(OCTFILE)
	$(call verb,Build complete.)

matlab:
	@echo "MATLAB provides built-in h5read, h5write, h5create, h5info, etc."
	@echo "No compilation is needed from lib/hdf5."

$(OCTFILE): $(SRC)
	$(MKOCTFILE) $(MKOCTFLAGS) $(H5FLAGS) $(INCLUDES) $(LIBS) hdf5oct.cc -o $(OCTFILE)

clean:
	$(DEL) *.o

distclean: clean
	$(DEL) $(OCTFILE)

distclean-all: distclean
