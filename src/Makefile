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
LFLAGS := $(shell $(MKOCTFILE) -p LFLAGS) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)

AR := $(shell $(MKOCTFILE) -p AR)

PKG_CXXFLAGS := -Wall -Wno-deprecated-declarations $(PKG_CXXFLAGS_APPEND)

all: __control_slicot_functions__.oct \
     __control_helper_functions__.oct

# TODO: Private oct-files for control package.

# unpack and compile SLICOT library
# Note that TG04BX is a custom routine.
# It has the extension .fortran such that
# it is not deleted by rm *.f when using
# the developer makefile makefile_control.m
# The other rotuines *.fortan are modified versions
# of the original slicot routines.
slicotlibrary.a: slicot.tar.gz
	tar -xzf slicot.tar.gz
	mkdir sltmp
	mv slicot/src/*.f ./sltmp
	mv slicot/src_aux/*.f ./sltmp
	cp SG03AD.fortran ./sltmp/SG03AD.f
	cp SG03BD.fortran ./sltmp/SG03BD.f
	cp SB04OD.fortran ./sltmp/SB04OD.f
	cp TG04BX.fortran ./sltmp/TG04BX.f
	cp TB01ZD.fortran ./sltmp/TB01ZD.f
	cd sltmp; $(MKOCTFILE) -c *.f
	$(AR) -rc slicotlibrary.a ./sltmp/*.o
	rm -rf sltmp slicot

# slicot functions
__control_slicot_functions__.oct: __control_slicot_functions__.cc common.cc slicotlibrary.a
	LFLAGS="$(LFLAGS)" \
    $(MKOCTFILE) $(PKG_CXXFLAGS) __control_slicot_functions__.cc common.cc slicotlibrary.a

# helper functions
__control_helper_functions__.oct: __control_helper_functions__.cc
	$(MKOCTFILE) $(PKG_CXXFLAGS) __control_helper_functions__.cc

clean:
	$(RM) -r *.o core octave-core *.oct *~ *.f slicot sltmp

realclean: clean
	$(RM) -r *.a

## This should also remove any configure cache which clean should not
## remove according to GNU guidelines.
## https://www.gnu.org/prep/standards/html_node/Standard-Targets.html
distclean: clean realclean
