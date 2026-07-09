#!/usr/bin/env make -f
# @file_array MEX build for PALM
#
# Builds file2mat, mat2file, and init for memory-mapped NIFTI I/O.
#
# Targets:
#   make matlab   - build for MATLAB  (.mexa64, .mexmaca64, .mexw64, ...)
#   make octave   - build for Octave  (.mex)
#   make          - auto-detect mex or mkoctfile in PATH
#   make help     - show usage
#   make clean    - remove object files
#   make distclean - remove binaries for the selected RUNTIME/PLATFORM
#
# Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

include Makefile.var

SPMMEX = file2mat.$(SUF) mat2file.$(SUF) init.$(SUF)

.PHONY: all matlab octave help clean distclean distclean-all check-tools

help:
	@echo "Build @file_array MEX files for PALM"
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

%.$(SUF): %.c
	$(MEX) $< $(MEXEND)

clean:
	$(DEL) *.o *.obj

distclean: clean
	$(DEL) $(SPMMEX)

distclean-all: clean
	$(DEL) file2mat.mex file2mat.mexglx file2mat.mexa64 file2mat.mexmaci64 \
	       file2mat.mexmaca64 file2mat.mexw32 file2mat.mexw64
	$(DEL) mat2file.mex mat2file.mexglx mat2file.mexa64 mat2file.mexmaci64 \
	       mat2file.mexmaca64 mat2file.mexw32 mat2file.mexw64
	$(DEL) init.mex init.mexglx init.mexa64 init.mexmaci64 \
	       init.mexmaca64 init.mexw32 init.mexw64

install:
	@:

tarball: all
	$(TAR) cf file_array_mex.tar $(SPMMEX)
