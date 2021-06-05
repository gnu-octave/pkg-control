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
#
# MA02ID.f use in its original version is compiled by an
# extra command suppressing warnings on indexing errors
#
slicotlibrary.a: slicot.tar.gz
	tar -xzf slicot.tar.gz
	mkdir sltmp
	mv slicot/src/*.f ./sltmp
	mv slicot/src_aux/*.f ./sltmp
	if [ "$(HAVE_DGGES)" = "1" ]; then \
		echo "copy routines using DGGES"; \
		cp SB04OD.fortran ./sltmp/SB04OD.f; \
		cp SG03AD.fortran ./sltmp/SG03AD.f; \
		cp SG03BD.fortran ./sltmp/SG03BD.f; \
	fi;
	cp AB08NX.fortran ./sltmp/AB08NX.f
	cp AG08BY.fortran ./sltmp/AG08BY.f
	cp SB01BY.fortran ./sltmp/SB01BY.f
	cp SB01FY.fortran ./sltmp/SB01FY.f
	cp SB06ND.fortran ./sltmp/SB06ND.f
	cp TB01MD.fortran ./sltmp/TB01MD.f
	cp TB01ND.fortran ./sltmp/TB01ND.f
	cp TB01ZD.fortran ./sltmp/TB01ZD.f
	cp TG04BX.fortran ./sltmp/TG04BX.f
	cp ODLTZM.fortran ./sltmp/ODLTZM.f
	cp makefile.slicot ./sltmp/makefile
	cd sltmp; $(MKOCTFILE) -w -c MA02ID.f; rm MA02ID.f; $(MKOCTFILE) -c *.f;
	$(AR) -rc slicotlibrary.a ./sltmp/*.o
	rm -rf sltmp slicot

# slicot functions
__control_slicot_functions__.oct: __control_slicot_functions__.cc common.cc slicotlibrary.a
	LDFLAGS="$(LDFLAGS)" \
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
