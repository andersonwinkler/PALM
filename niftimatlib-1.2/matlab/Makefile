#!make -f
#
# The tardist file for the nifti1 matlab code contains precompiled
# mex files for the PCWIN, MAC, SOL2, LNX86, GLNXA64 platforms, so
# you do not need to do this build if you have one of those platforms.
# Run "uname" at your shell prompt to see your platform code.
# 
# 
# Two files need to be compiled with the MATLAB mex compiler: 
#	file2mat.c 
#	mat2file.c
# These files are in the @file_array/private/src directory.
# A successful compile will leave two "mex" files in the @file_array/private directory:
# 	file2mat.mex<platform suffix>  
#	mat2file.mex<platform suffix>
#
# To make the mex files, in this matlab directory type:
#	make all
#
# 
# Problems ?
# If you need to make platform specific tweaks to the mex call, you
# may want to look at the SPM8 (http://www.fil.ion.ucl.ac.uk/spm/)
# distribution src/Makefile and src/Makefile.var that has more detailed flag settings.
#
# For Windows you need the CygWin environment (http://sourceware.org/cygwin/).
# Install the Default configuration, along with the gcc C compiler, and make
# from the "Devel" options.
# You may want to refer to Matthew Brett's web page on mex files for Windows:
# (http://gnumex.sourceforge.net/)
#


## Defines
MEX = mex -O
CC  = cc
DIR = @file_array/private/src



## Default targets
unknown:
	@ make `uname`

all:
	@ make `uname`

help:
	@echo "make all to compile the file2mat and mat2file mex libraries"
	@echo ""


## Platforms
SunOS:	big_endian
IRIX:	big_endian
IRIX64:	big_endian
HP-UX:	big_endian
AIX:	big_endian
OSF1:	big_endian

# for MAC power pc reset to big_endian
Darwin:	little_endian

Linux:	little_endian
Linux.A64:	little_endian

CYGWIN_NT-4.0: windows
CYGWIN_NT-5.0: windows
CYGWIN_NT-5.1: windows



## Compiles
big_endian:
	(cd $(DIR); $(MEX)  -DBIGENDIAN file2mat.c; mv *mex* ..)
	(cd $(DIR); $(MEX)  -DBIGENDIAN mat2file.c; mv *mex* ..)

little_endian:
	(cd $(DIR); $(MEX)  file2mat.c; mv *mex* ..)
	(cd $(DIR); $(MEX)  mat2file.c; mv *mex* ..)

windows:
	(cd $(DIR); $(MEX)  -DSPM_WIN32 file2mat.c; mv *.dll ..)
	(cd $(DIR); $(MEX)  -DSPM_WIN32 mat2file.c; mv *.dll ..)

