OCTAVE ?= octave

LAPACK_LIBS := $(shell mkoctfile -p LAPACK_LIBS)
BLAS_LIBS := $(shell mkoctfile -p BLAS_LIBS)
FLIBS := $(shell mkoctfile -p FLIBS)

all: slab08nd.oct slab13dd.oct slsb10hd.oct slsb10ed.oct slab13bd.oct \
     slsb01bd.oct slsb10fd.oct slsb10dd.oct slsb03md.oct slsb04md.oct \
     slsb04qd.oct slsg03ad.oct slsb02od.oct slab13ad.oct slab01od.oct \
     sltb01pd.oct slsb03od.oct slsg03bd.oct slag08bd.oct sltg01jd.oct \
     sltg01hd.oct sltg01id.oct slsg02ad.oct sltg04bx.oct sltb01id.oct \
     sltg01ad.oct slsb10id.oct slsb10kd.oct slsb10zd.oct sltb04bd.oct \
     slab04md.oct slsb10jd.oct sltd04ad.oct sltb01ud.oct \
     is_real_scalar.oct is_real_vector.oct is_real_matrix.oct \
     is_real_square_matrix.oct \
     copy_NEWS

copy_NEWS:
	$(if $(shell $(OCTAVE) -q --eval "if (compare_versions (version (), \"3.6.0\", \"<\")) printf (\"true\"); endif"), mv ../NEWS ../doc/)

# TODO: Leading and trailing underscores for sl* functions
#       (__sl*__.oct) would be nice, but this can be an issue
#       for fortran compilers.  Private oct-files for the
#       control package would be even better.

slicotlibrary.a: slicot.tar.gz
	tar -xzf slicot.tar.gz
	mv slicot/src/*.f .
	mv slicot/src_aux/*.f .
	cp TG04BX.fortran TG04BX.f
	mkoctfile *.f \
	${LAPACK_LIBS} ${BLAS_LIBS} ${FLIBS}
	ar -r slicotlibrary.a *.o
	rm -rf *.o *.f AB01MD.oct slicot

# transmission zeros of state-space models
slab08nd.oct: slab08nd.cc slicotlibrary.a
	mkoctfile slab08nd.cc slicotlibrary.a \
              ${LAPACK_LIBS} ${BLAS_LIBS}	

# L-infinity norm
slab13dd.oct: slab13dd.cc slicotlibrary.a
	mkoctfile slab13dd.cc slicotlibrary.a \
              ${LAPACK_LIBS} ${BLAS_LIBS} ${FLIBS}

# H-2 controller synthesis - continuous-time
slsb10hd.oct: slsb10hd.cc slicotlibrary.a
	mkoctfile slsb10hd.cc slicotlibrary.a \
              ${LAPACK_LIBS} ${BLAS_LIBS}

# H-2 controller synthesis - discrete-time
slsb10ed.oct: slsb10ed.cc slicotlibrary.a
	mkoctfile slsb10ed.cc slicotlibrary.a \
              ${LAPACK_LIBS} ${BLAS_LIBS}

# H-2 norm
slab13bd.oct: slab13bd.cc slicotlibrary.a
	mkoctfile slab13bd.cc slicotlibrary.a \
              ${LAPACK_LIBS} ${BLAS_LIBS}

# Pole assignment
slsb01bd.oct: slsb01bd.cc slicotlibrary.a
	mkoctfile slsb01bd.cc slicotlibrary.a \
              ${LAPACK_LIBS} ${BLAS_LIBS}

# H-inf controller synthesis - continuous-time
slsb10fd.oct: slsb10fd.cc slicotlibrary.a
	mkoctfile slsb10fd.cc slicotlibrary.a \
              ${LAPACK_LIBS} ${BLAS_LIBS}

# H-inf controller synthesis - discrete-time
slsb10dd.oct: slsb10dd.cc slicotlibrary.a
	mkoctfile slsb10dd.cc slicotlibrary.a \
              ${LAPACK_LIBS} ${BLAS_LIBS}

# Lyapunov equations
slsb03md.oct: slsb03md.cc slicotlibrary.a
	mkoctfile slsb03md.cc slicotlibrary.a \
              ${LAPACK_LIBS} ${BLAS_LIBS}

# Sylvester equations - continuous-time
slsb04md.oct: slsb04md.cc slicotlibrary.a
	mkoctfile slsb04md.cc slicotlibrary.a \
              ${LAPACK_LIBS} ${BLAS_LIBS}

# Sylvester equations - discrete-time
slsb04qd.oct: slsb04qd.cc slicotlibrary.a
	mkoctfile slsb04qd.cc slicotlibrary.a \
              ${LAPACK_LIBS} ${BLAS_LIBS}

# generalized Lyapunov equations
slsg03ad.oct: slsg03ad.cc slicotlibrary.a
	mkoctfile slsg03ad.cc slicotlibrary.a \
              ${LAPACK_LIBS} ${BLAS_LIBS}

# algebraic Riccati equations
slsb02od.oct: slsb02od.cc slicotlibrary.a
	mkoctfile slsb02od.cc slicotlibrary.a \
              ${LAPACK_LIBS} ${BLAS_LIBS}

# Hankel singular values
slab13ad.oct: slab13ad.cc slicotlibrary.a
	mkoctfile slab13ad.cc slicotlibrary.a \
              ${LAPACK_LIBS} ${BLAS_LIBS}

# staircase form using orthogonal transformations
slab01od.oct: slab01od.cc slicotlibrary.a
	mkoctfile slab01od.cc slicotlibrary.a \
              ${LAPACK_LIBS} ${BLAS_LIBS}

# minimal realization of state-space models
sltb01pd.oct: sltb01pd.cc slicotlibrary.a
	mkoctfile sltb01pd.cc slicotlibrary.a \
              $(LAPACK_LIBS) $(BLAS_LIBS)

# Cholesky factor of Lyapunov equations
slsb03od.oct: slsb03od.cc slicotlibrary.a
	mkoctfile slsb03od.cc slicotlibrary.a \
              ${LAPACK_LIBS} ${BLAS_LIBS}

# Cholesky factor of generalized Lyapunov equations
slsg03bd.oct: slsg03bd.cc slicotlibrary.a
	mkoctfile slsg03bd.cc slicotlibrary.a \
              ${LAPACK_LIBS} ${BLAS_LIBS}

# transmission zeros of descriptor state-space models
slag08bd.oct: slag08bd.cc slicotlibrary.a
	mkoctfile slag08bd.cc slicotlibrary.a \
              ${LAPACK_LIBS} ${BLAS_LIBS}

# minimal realization of descriptor state-space models
sltg01jd.oct: sltg01jd.cc slicotlibrary.a
	mkoctfile sltg01jd.cc slicotlibrary.a \
              ${LAPACK_LIBS} ${BLAS_LIBS}

# controllability staircase form of descriptor state-space models
sltg01hd.oct: sltg01hd.cc slicotlibrary.a
	mkoctfile sltg01hd.cc slicotlibrary.a \
              ${LAPACK_LIBS} ${BLAS_LIBS}

# observability staircase form of descriptor state-space models
sltg01id.oct: sltg01id.cc slicotlibrary.a
	mkoctfile sltg01id.cc slicotlibrary.a \
              ${LAPACK_LIBS} ${BLAS_LIBS}

# solution of algebraic Riccati equations for descriptor systems
slsg02ad.oct: slsg02ad.cc slicotlibrary.a
	mkoctfile slsg02ad.cc slicotlibrary.a \
              ${LAPACK_LIBS} ${BLAS_LIBS}

# gain of descriptor state-space models
sltg04bx.oct: sltg04bx.cc slicotlibrary.a
	mkoctfile sltg04bx.cc slicotlibrary.a \
              ${LAPACK_LIBS} ${BLAS_LIBS}

# scaling of state-space models
sltb01id.oct: sltb01id.cc slicotlibrary.a
	mkoctfile sltb01id.cc slicotlibrary.a \
              ${LAPACK_LIBS} ${BLAS_LIBS}

# scaling of descriptor state-space models
sltg01ad.oct: sltg01ad.cc slicotlibrary.a
	mkoctfile sltg01ad.cc slicotlibrary.a \
              ${LAPACK_LIBS} ${BLAS_LIBS}

# H-infinity loop shaping - continuous-time
slsb10id.oct: slsb10id.cc slicotlibrary.a
	mkoctfile slsb10id.cc slicotlibrary.a \
              ${LAPACK_LIBS} ${BLAS_LIBS}

# H-infinity loop shaping - discrete-time - strictly proper case
slsb10kd.oct: slsb10kd.cc slicotlibrary.a
	mkoctfile slsb10kd.cc slicotlibrary.a \
              ${LAPACK_LIBS} ${BLAS_LIBS}

# H-infinity loop shaping - discrete-time - proper case
slsb10zd.oct: slsb10zd.cc slicotlibrary.a
	mkoctfile slsb10zd.cc slicotlibrary.a \
              ${LAPACK_LIBS} ${BLAS_LIBS}

# State-space to transfer function conversion
sltb04bd.oct: sltb04bd.cc slicotlibrary.a
	mkoctfile sltb04bd.cc slicotlibrary.a \
              ${LAPACK_LIBS} ${BLAS_LIBS}

# bilinear transformation
slab04md.oct: slab04md.cc slicotlibrary.a
	mkoctfile slab04md.cc slicotlibrary.a \
              ${LAPACK_LIBS} ${BLAS_LIBS}

# descriptor to regular state-space conversion
slsb10jd.oct: slsb10jd.cc slicotlibrary.a
	mkoctfile slsb10jd.cc slicotlibrary.a \
              ${LAPACK_LIBS} ${BLAS_LIBS}

# transfer function to state-space conversion
sltd04ad.oct: sltd04ad.cc slicotlibrary.a
	mkoctfile sltd04ad.cc slicotlibrary.a \
              ${LAPACK_LIBS} ${BLAS_LIBS}

# controllable block Hessenberg realization
sltb01ud.oct: sltb01ud.cc slicotlibrary.a
	mkoctfile sltb01ud.cc slicotlibrary.a \
              ${LAPACK_LIBS} ${BLAS_LIBS}

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