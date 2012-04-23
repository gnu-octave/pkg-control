# LAPACK_LIBS := $(shell mkoctfile -p LAPACK_LIBS)
# BLAS_LIBS := $(shell mkoctfile -p BLAS_LIBS)
FLIBS := $(shell mkoctfile -p FLIBS)

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
	mkoctfile -c *.f
	ar -rc slicotlibrary.a *.o
	rm -rf *.o *.f slicot

lapacklibrary.a: lapack-3.4.1.tgz
	tar -xzf lapack-3.4.1.tgz
	mv lapack-3.4.1/BLAS/SRC/*.f .
	mv lapack-3.4.1/INSTALL/*.f .
	mv lapack-3.4.1/SRC/*.f .
	mkoctfile -c *.f
	ar -rc lapacklibrary.a *.o
	rm -rf *.o *.f lapack-3.4.1

# slicot functions
control_slicot_functions.oct: control_slicot_functions.cc slicotlibrary.a lapacklibrary.a
	mkoctfile control_slicot_functions.cc common.cc slicotlibrary.a lapacklibrary.a \
              ${FLIBS}

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
	rm -rf *.o core octave-core *.oct *~ *.f slicot lapack-3.4.1

realclean: clean
	rm -rf *.a