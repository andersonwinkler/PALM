#!/usr/bin/env make -f
# @nifti/private MEX build for PALM
#
# Builds nifti_stats from nifti_stats.c + nifti_stats_mex.c (NIFTI intent
# stat-code conversions). Optional for PALM; basic NIFTI I/O does not need it.
#
# Targets:
#   make matlab   - build for MATLAB  (.mexa64, .mexmaca64, .mexw64, ...)
#   make octave   - build for Octave  (.mex)
#   make          - auto-detect mex or mkoctfile in PATH
#   make help     - show usage
#
# Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

include Makefile.var

SPMMEX = nifti_stats.$(SUF)
SRCS   = nifti_stats.c nifti_stats_mex.c

# MATLAB: pass -std=c99 through CFLAGS (Octave gets it from MEXOPTS in Makefile.var)
ifneq (octave,$(RUNTIME))
  ifeq (windows,$(PLATFORM))
    MEXOPTS += CFLAGS=\"$$CFLAGS -std=c99\"
  else
    MEXOPTS += CFLAGS='$$CFLAGS -std=c99'
  endif
  MEX := $(MEXBIN) $(MEXOPTS)
endif

# Math library (needed on Linux for nifti_stats.c)
ifneq (,$(filter Linux MSYS,$(PLATFORM)))
  MEXLIBS = -lm
else ifeq (octave,$(RUNTIME))
  MEXLIBS = -lm
endif

.PHONY: all matlab octave help clean distclean distclean-all check-tools

help:
	@echo "Build @nifti/private MEX files for PALM"
	@echo ""
	@echo "Targets:"
	@echo "  make matlab        Build for MATLAB (platform-specific .mex* suffix)"
	@echo "  make octave        Build for Octave  (generic .mex suffix)"
	@echo "  make               Auto-detect mex or mkoctfile"
	@echo "  make clean         Remove intermediate object files"
	@echo "  make distclean     Remove binaries for the current RUNTIME/PLATFORM"
	@echo "  make distclean-all Remove all known MEX variants in this directory"
	@echo ""
	@echo "Variables:"
	@echo "  RUNTIME=matlab|octave|auto   (default: auto)"
	@echo "  PLATFORM=Linux|Darwin|windows|MSYS|octave"
	@echo "  SUF=mexa64|mexmaca64|mexw64|mex|..."
	@echo ""
	@echo "Detected: RUNTIME=$(RUNTIME) PLATFORM=$(PLATFORM) SUF=$(SUF)"
	@echo "Compiler: $(MEX)"

check-tools:
	@if [ "$(RUNTIME)" = "matlab" ] && ! command -v mex >/dev/null 2>&1; then \
		echo "Error: mex not found in PATH."; exit 1; \
	fi
	@if [ "$(RUNTIME)" = "octave" ] && ! command -v mkoctfile >/dev/null 2>&1; then \
		echo "Error: mkoctfile not found in PATH."; exit 1; \
	fi

all: check-tools $(SPMMEX)
	$(call verb,Build complete.)

matlab:
	@$(MAKE) RUNTIME=matlab all

octave:
	@$(MAKE) RUNTIME=octave all

$(SPMMEX): $(SRCS) nifti1.h
	$(MEX) $(SRCS) $(MEXLIBS) $(MEXEND)

clean:
	$(DEL) *.o *.obj

distclean: clean
	$(DEL) $(SPMMEX) nifti_stats_mex.$(SUF)

distclean-all: clean
	$(DEL) nifti_stats.mex nifti_stats.mexglx nifti_stats.mexa64 \
	       nifti_stats.mexmaci64 nifti_stats.mexmaca64 \
	       nifti_stats.mexw32 nifti_stats.mexw64
	$(DEL) nifti_stats_mex.mex nifti_stats_mex.mexglx nifti_stats_mex.mexa64 \
	       nifti_stats_mex.mexmaci64 nifti_stats_mex.mexmaca64 \
	       nifti_stats_mex.mexw32 nifti_stats_mex.mexw64

install:
	@:

tarball: all
	$(TAR) cf nifti_stats_mex.tar $(SPMMEX)
