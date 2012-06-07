MKOCTFILE ?= mkoctfile

LAPACK_LIBS := $(shell $(MKOCTFILE) -p LAPACK_LIBS)
BLAS_LIBS := $(shell $(MKOCTFILE) -p BLAS_LIBS)
FLIBS := $(shell $(MKOCTFILE) -p FLIBS)
LFLAGS := $(shell $(MKOCTFILE) -p LFLAGS) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)

all: control_slicot_functions.oct \
     is_real_scalar.oct \
     is_real_vector.oct \
     is_real_matrix.oct \
     is_real_square_matrix.oct

# TODO: Private oct-files for control package.

# unpack and compile SLICOT library
# Note that TG04BX is a custom routine.
# It has the extension .fortran such that
# it is not deleted by rm *.f when using
# the developer makefile makefile_control.m
slicotlibrary.a: slicot.tar.gz
	tar -xzf slicot.tar.gz
	mv slicot/src/*.f .
	mv slicot/src_aux/*.f .
	cp TG04BX.fortran TG04BX.f
	$(MKOCTFILE) -c *.f
	ar -rc slicotlibrary.a *.o
	rm -rf *.o *.f slicot

# slicot functions
control_slicot_functions.oct: control_slicot_functions.cc slicotlibrary.a
	LFLAGS="$(LFLAGS)" \
    $(MKOCTFILE) control_slicot_functions.cc common.cc slicotlibrary.a

# helpers
is_real_scalar.oct: is_real_scalar.cc
	$(MKOCTFILE) is_real_scalar.cc

is_real_vector.oct: is_real_vector.cc
	$(MKOCTFILE) is_real_vector.cc

is_real_matrix.oct: is_real_matrix.cc
	$(MKOCTFILE) is_real_matrix.cc

is_real_square_matrix.oct: is_real_square_matrix.cc
	$(MKOCTFILE) is_real_square_matrix.cc

clean:
	rm -rf *.o core octave-core *.oct *~ *.f slicot

realclean: clean
	rm -rf *.a