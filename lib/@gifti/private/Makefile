#!/usr/bin/env make -f
# @gifti MEX build for PALM
#
# Builds zstream, xml_parser, and base64 for GIFTI surface I/O.
#
# Targets:
#   make matlab   - build for MATLAB  (.mexa64, .mexmaca64, .mexw64, ...)
#   make octave   - build for Octave  (.mex)
#   make          - auto-detect mex or mkoctfile in PATH
#   make help     - show usage
#   make clean    - remove object files
#   make distclean - remove binaries for the selected RUNTIME/PLATFORM
#
# Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

include Makefile.var

SPMMEX = zstream.$(SUF) xml_parser.$(SUF) base64.$(SUF)

# MATLAB: pass -std=c99 through CFLAGS (Octave gets it from MEXOPTS in Makefile.var)
ifneq (octave,$(RUNTIME))
  ifeq (windows,$(PLATFORM))
    MEXOPTS += CFLAGS=\"$$CFLAGS -std=c99\"
  else
    MEXOPTS += CFLAGS='$$CFLAGS -std=c99'
  endif
  MEX := $(MEXBIN) $(MEXOPTS)
endif

.PHONY: all matlab octave help clean distclean distclean-all check-tools

help:
	@echo "Build @gifti MEX files for PALM"
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

# xml_parser links yxml.c (yxml XML tokenizer)
xml_parser.$(SUF): xml_parser.c yxml.c yxml.h
	$(MEX) xml_parser.c yxml.c $(MEXEND)

clean:
	$(DEL) *.o *.obj

distclean: clean
	$(DEL) $(SPMMEX)

distclean-all: clean
	$(DEL) zstream.mex zstream.mexglx zstream.mexa64 zstream.mexmaci64 \
	       zstream.mexmaca64 zstream.mexw32 zstream.mexw64
	$(DEL) xml_parser.mex xml_parser.mexglx xml_parser.mexa64 \
	       xml_parser.mexmaci64 xml_parser.mexmaca64 xml_parser.mexw32 \
	       xml_parser.mexw64
	$(DEL) base64.mex base64.mexglx base64.mexa64 base64.mexmaci64 \
	       base64.mexmaca64 base64.mexw32 base64.mexw64

install:
	@:

tarball: all
	$(TAR) cf gifti_mex.tar $(SPMMEX)
