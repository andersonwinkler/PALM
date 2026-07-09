#!/usr/bin/env make -f
# PALM lib/ — build all bundled compiled extensions for this platform
#
# Run from lib/:
#   make            Auto-detect MATLAB or Octave and build everything
#   make matlab       Build all MEX files for MATLAB (.mexa64, .mexmaca64, ...)
#   make octave       Build all MEX/.oct files for Octave (.mex, hdf5oct.oct)
#   make help         Show usage
#   make clean        Remove object files in all components
#   make distclean    Remove binaries for the current RUNTIME in all components
#   make distclean-all Remove every known MEX/.oct variant in all components
#
# Override runtime detection:
#   make RUNTIME=octave
#   make RUNTIME=matlab
#
# Components built:
#   @file_array/private   file2mat, mat2file, init        (NIFTI memory-map; required)
#   @gifti/private        zstream, xml_parser, base64    (GIFTI I/O; required on Octave)
#   @xmltree/private      xml_findstr                    (optional, faster XML)
#   @xmltreemod/private   xml_findstr                    (optional, CIFTI XML)
#   extras                spm_existfile                  (optional, faster file checks)
#   @nifti/private        nifti_stats                    (optional, stat intent codes)
#   hdf5                  hdf5oct.oct                    (Octave HDF5; MATLAB uses built-ins)

include Makefile.var

matlab: RUNTIME = matlab
octave: RUNTIME = octave

.PHONY: all matlab octave help clean distclean distclean-all \
        build-mex check-tools $(MEX_DIRS)

all: check-tools
ifeq (matlab,$(RUNTIME))
	@$(MAKE) matlab
else
	@$(MAKE) octave
endif

matlab: check-tools build-mex
	@$(MAKE) -C hdf5 matlab

octave: check-tools build-mex
	@$(MAKE) -C hdf5 octave

check-tools:
	@if [ "$(RUNTIME)" = "matlab" ] && ! command -v mex >/dev/null 2>&1; then \
		echo "Error: mex not found in PATH."; exit 1; \
	fi
	@if [ "$(RUNTIME)" = "octave" ] && ! command -v mkoctfile >/dev/null 2>&1; then \
		echo "Error: mkoctfile not found in PATH."; exit 1; \
	fi

build-mex: $(MEX_DIRS)

$(MEX_DIRS):
	$(call banner,Building $@ for $(RUNTIME))
	@$(MAKE) -C $@ RUNTIME=$(RUNTIME) all

help:
	@echo "PALM lib/ — build all bundled compiled extensions"
	@echo ""
	@echo "Targets:"
	@echo "  make               Auto-detect runtime and build all components"
	@echo "  make matlab        Build for MATLAB"
	@echo "  make octave        Build for Octave (includes hdf5oct.oct)"
	@echo "  make clean         Remove intermediate object files"
	@echo "  make distclean     Remove binaries for the current RUNTIME"
	@echo "  make distclean-all Remove all MEX/.oct variants"
	@echo ""
	@echo "Variables:"
	@echo "  RUNTIME=matlab|octave|auto   (default: auto; detected: $(RUNTIME))"
	@echo ""
	@echo "Components:"
	@for d in $(MEX_DIRS); do echo "  $$d"; done
	@echo "  hdf5/ (Octave only)"
	@echo ""
	@echo "After building, run PALM as usual; palm_checkprogs adds lib/ to the path."

clean:
	@for d in $(MEX_DIRS); do \
		$(MAKE) -C $$d clean; \
	done
	@$(MAKE) -C hdf5 clean

distclean:
	@for d in $(MEX_DIRS); do \
		$(MAKE) -C $$d RUNTIME=$(RUNTIME) distclean; \
	done
	@$(MAKE) -C hdf5 distclean

distclean-all:
	@for d in $(MEX_DIRS); do \
		$(MAKE) -C $$d distclean-all; \
	done
	@$(MAKE) -C hdf5 distclean-all
