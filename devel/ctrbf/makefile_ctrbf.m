mkoctfile sltb01ud.cc \
          TB01UD.f MB01PD.f MB03OY.f MB01QD.f \
          "$(mkoctfile -p LAPACK_LIBS)" \
          "$(mkoctfile -p BLAS_LIBS)"
