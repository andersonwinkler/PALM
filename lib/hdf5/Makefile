#!/usr/bin/env make -f
# lib/hdf5 build wrapper for PALM (hdf5oct)
#
# Delegates to src/ where hdf5oct.oct is built.
#
# Targets:
#   make help     - show usage
#   make octave   - build hdf5oct.oct for Octave
#   make matlab   - informational (MATLAB uses built-in h5* functions)
#   make clean    - remove object files
#   make distclean - remove hdf5oct.oct

.PHONY: help octave matlab clean distclean distclean-all check-tools all

help octave matlab clean distclean distclean-all check-tools all:
	@$(MAKE) -C src $@
