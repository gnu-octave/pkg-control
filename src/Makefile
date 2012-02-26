OCTAVE ?= octave

LAPACK_LIBS := $(shell mkoctfile -p LAPACK_LIBS)
BLAS_LIBS := $(shell mkoctfile -p BLAS_LIBS)
FLIBS := $(shell mkoctfile -p FLIBS)

all: control_slicot_functions.oct \
     is_real_scalar.oct is_real_vector.oct is_real_matrix.oct \
     is_real_square_matrix.oct \
     copy_NEWS

copy_NEWS:
	$(if $(shell $(OCTAVE) -q --eval "if (compare_versions (version (), \"3.6.0\", \"<\")) printf (\"true\"); endif"), mv ../NEWS ../doc/)

# TODO: Leading and trailing underscores for sl* functions
#       (__sl*__.oct) would be nice, but this can be an issue
#       for fortran compilers.  Private oct-files for the
#       control package would be even better.

# TODO: Create a single oct-file control_slicot.oct

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
	mkoctfile -c *.f
	ar -r slicotlibrary.a *.o
	rm -rf *.o *.f slicot

# slicot functions
control_slicot_functions.oct: control_slicot_functions.cc slicotlibrary.a
	mkoctfile control_slicot_functions.cc slicotlibrary.a \
              ${LAPACK_LIBS} ${BLAS_LIBS} ${FLIBS}

# helpers
is_real_scalar.oct: is_real_scalar.cc
	mkoctfile is_real_scalar.cc

is_real_vector.oct: is_real_vector.cc
	mkoctfile is_real_vector.cc

is_real_matrix.oct: is_real_matrix.cc
	mkoctfile is_real_matrix.cc

is_real_square_matrix.oct: is_real_square_matrix.cc
	mkoctfile is_real_square_matrix.cc

clean:
	rm -rf *.o core octave-core *.oct *~ *.f slicot

realclean: clean
	rm -rf *.a