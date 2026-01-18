-include ./Makefile.conf

# the following vars for libs are only defined if mkoctfile is defined
# and this file might be called without mkoctfile defined (clean, distclean) 
ifdef MKOCTFILE
	ifndef LAPACK_LIBS
		LAPACK_LIBS := $(shell $(MKOCTFILE) -p LAPACK_LIBS)
	endif
	ifndef BLAS_LIBS
		BLAS_LIBS := $(shell $(MKOCTFILE) -p BLAS_LIBS)
	endif
	ifndef FLIBS
		FLIBS := $(shell $(MKOCTFILE) -p FLIBS)
	endif
	LDFLAGS := $(shell $(MKOCTFILE) -p LDFLAGS) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)
	# other definitions
	AR := $(shell $(MKOCTFILE) -p AR)
endif

# SLICOT files, relative to src directory
SC        := slicot
SC_SUBMOD := ../slicot-reference
# subdirectories relative to $(SC), $(SC)/$(SC_SRC) is SLICOT build directory
SC_SRC    := src
SC_BUILD  := $(SC)/$(SC_SRC)
# own sources related to SLICOT relative to build directory
SC_LAPACK := lapack_aux
SRC       := ../..

PKG_CXXFLAGS := -Wall $(PKG_CXXFLAGS_APPEND)

all: __control_slicot_functions__.oct \
     __control_helper_functions__.oct

# For being able to directly call Makefile from repository instead of
# distributed file structure, copy required files into the slicot
# source directory. The file TB01ZD.f is removed afterwards because
# there is an own version of this file using a different algorithm.

# MA02ID.f : One exemplary file from SLICOT, if not present, other files
#            are most probably missing as well. Copy all files from the
#            src directories of the SLICOT submodule.
#            If we are in the dist tarball, the *.f files exist already.
$(SC_BUILD)/MA02ID.f:
	mkdir -p $(SC_BUILD)/$(SC_LAPACK)
	cp $(SC_SUBMOD)/$(SC_SRC)/*.f $(SC_BUILD)
	cp $(SC_SUBMOD)/$(SC_SRC)/$(SC_LAPACK)/*.f $(SC_BUILD)/$(SC_LAPACK)
	$(RM) $(SC_BUILD)/TB01ZD.f
	./bootstrap && $(RM) -r "autom4te.cache"

# TODO: Private oct-files for control package.

# Compilation of SLICOT library
#
#   Compile in $(SC_BUILD):  *.f                SLICOT files, copied above
#                            $(SRC)/*.f         own derived files
#                            $(SC_LAPACK)/*.f   deprecated LAPACK files
#
#   MA02ID.f is compiled individually, suppressing indexing warnings
#
slicotlibrary.a: slicot $(SC)/$(SC_SRC)/MA02ID.f
	cd $(SC_BUILD) && \
	  $(MKOCTFILE) -w -c MA02ID.f; mv MA02ID.f x && \
		$(MKOCTFILE) -c *.f $(SRC)/*.f $(SC_LAPACK)/*.f && \
		mv x MA02ID.f
	$(AR) -rc slicotlibrary.a $(SC_BUILD)/*.o

# slicot functions
__control_slicot_functions__.oct: __control_slicot_functions__.cc common.cc slicotlibrary.a
	LDFLAGS='$(LDFLAGS)' \
    $(MKOCTFILE) $(PKG_CXXFLAGS) __control_slicot_functions__.cc common.cc slicotlibrary.a

# helper functions
__control_helper_functions__.oct: __control_helper_functions__.cc
	$(MKOCTFILE) $(PKG_CXXFLAGS) __control_helper_functions__.cc

clean:
	$(RM) -r *.o core octave-core *.oct *~ $(SC_BUILD)/*.o
	$(RM) -rf $(SC_BUILD)

realclean: clean
	$(RM) -r *.a

## This should also remove any configure cache which clean should not
## remove according to GNU guidelines.
## https://www.gnu.org/prep/standards/html_node/Standard-Targets.html
distclean: clean realclean
	$(RM) -f Makefile.conf config.h config.log config.status \
					 oct-alt-includes.h configure aclocal.m4
