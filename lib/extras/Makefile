#!/usr/bin/env make -f
# lib/extras MEX build for PALM
#
# Builds spm_existfile for fast file-existence checks (optional; pure-MATLAB
# fallback exists in spm_existfile.m).
#
# Other files in this directory (spm_fileparts, spm_type, spm_flip_analyze_images)
# are pure MATLAB and require no compilation.
#
# Targets:
#   make matlab   - build for MATLAB  (.mexa64, .mexmaca64, .mexw64, ...)
#   make octave   - build for Octave  (.mex)
#   make          - auto-detect mex or mkoctfile in PATH
#   make help     - show usage
#
# Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

include Makefile.var

SPMMEX = spm_existfile.$(SUF)

.PHONY: all matlab octave help clean distclean distclean-all check-tools

help:
	@echo "Build lib/extras MEX files for PALM"
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
	$(DEL) spm_existfile.mex spm_existfile.mexglx spm_existfile.mexa64 \
	       spm_existfile.mexmaci64 spm_existfile.mexmaca64 \
	       spm_existfile.mexw32 spm_existfile.mexw64

install:
	@:

tarball: all
	$(TAR) cf extras_mex.tar $(SPMMEX)
