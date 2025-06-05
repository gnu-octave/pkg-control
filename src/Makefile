-include ./Makefile.conf

MKOCTFILE ?= mkoctfile

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

# slicot relative to src directory, SC_SRC is the build directory
SC      := slicot
SC_SRC  := $(SC)/src
SC_SMOD := ../slicot-reference

# own and auxilary sources relative to build directory
SRC := ../..
SRC_AUX := $(SRC)/src_aux

# other definitions
AR := $(shell $(MKOCTFILE) -p AR)

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
$(SC_SRC)/MA02ID.f:
	mkdir -p $(SC_SRC)
	cp $(SC_SMOD)/src/*.f $(SC_SRC)
	cp $(SC_SMOD)/src_aux/*.f $(SC_SRC)
	$(RM) $(SC_SRC)/TB01ZD.f
	./bootstrap && $(RM) -r "autom4te.cache"

# TODO: Private oct-files for control package.

# Compile SLICOT library
#
# Compile: $(SC_SRC)/*.f                       (SLICOT files, copied above)
#          $(SRC)/TG04BX.f $(SRC)/TB01ZD.f     own derived files
#          $(SRC_AUX)/*.f                      external sources like
#                                              deprecated LAPACK files
#
# MA02ID.f is compiled by an extra command
# suppressing warnings on indexing errors
#
slicotlibrary.a: slicot $(SC_SRC)/MA02ID.f
	cd $(SC_SRC) && \
	  $(MKOCTFILE) -w -c MA02ID.f; mv MA02ID.f x && \
		$(MKOCTFILE) -c *.f $(SRC)/TB01ZD.f $(SRC)/TG04BX.f $(SRC_AUX)/*.f && \
		mv x MA02ID.f
	$(AR) -rc slicotlibrary.a $(SC_SRC)/*.o

# slicot functions
__control_slicot_functions__.oct: __control_slicot_functions__.cc common.cc slicotlibrary.a
	LDFLAGS="$(LDFLAGS)" \
    $(MKOCTFILE) $(PKG_CXXFLAGS) __control_slicot_functions__.cc common.cc slicotlibrary.a

# helper functions
__control_helper_functions__.oct: __control_helper_functions__.cc
	$(MKOCTFILE) $(PKG_CXXFLAGS) __control_helper_functions__.cc

clean:
	$(RM) -r *.o core octave-core *.oct *~ $(SC_SRC)/*.o
	$(RM) -rf $(SC_SRC)

realclean: clean
	$(RM) -r *.a

## This should also remove any configure cache which clean should not
## remove according to GNU guidelines.
## https://www.gnu.org/prep/standards/html_node/Standard-Targets.html
distclean: clean realclean
	rm -f Makefile.conf config.h config.log config.status oct-alt-includes.h
