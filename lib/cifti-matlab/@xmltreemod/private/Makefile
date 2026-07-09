#!/usr/bin/env make -f
# @xmltreemod MEX build for PALM / cifti-matlab
#
# Builds xml_findstr for faster XML string search (optional; pure-MATLAB
# fallback exists in xml_findstr.m).
#
# Targets:
#   make matlab   - build for MATLAB  (.mexa64, .mexmaca64, .mexw64, ...)
#   make octave   - build for Octave  (.mex)
#   make          - auto-detect mex or mkoctfile in PATH
#   make help     - show usage
#
# Copyright (C) 2008-2013 Wellcome Trust Centre for Neuroimaging

include Makefile.var

SPMMEX = xml_findstr.$(SUF)

.PHONY: all matlab octave help clean distclean distclean-all check-tools

help:
	@echo "Build @xmltreemod MEX files for PALM / cifti-matlab"
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
	$(DEL) xml_findstr.mex xml_findstr.mexglx xml_findstr.mexa64 \
	       xml_findstr.mexmaci64 xml_findstr.mexmaca64 \
	       xml_findstr.mexw32 xml_findstr.mexw64

install:
	@:

tarball: all
	$(TAR) cf xmltreemod_mex.tar $(SPMMEX)
