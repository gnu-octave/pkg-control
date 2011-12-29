## ==============================================================================
## Developer Makefile for OCT-files
## ==============================================================================
## USAGE: * fetch control from Octave-Forge by svn
##        * add control/inst, control/src and control/devel to your Octave path
##        * run makefile_*
## ==============================================================================

homedir = pwd ();
develdir = fileparts (which ("makefile_staircase"));
srcdir = [develdir, "/../src"];
cd (srcdir);

## staircase form using orthogonal transformations
mkoctfile slab01od.cc \
          AB01OD.f AB01ND.f MB03OY.f MB01PD.f MB01QD.f \
          "$(mkoctfile -p LAPACK_LIBS)" \
          "$(mkoctfile -p BLAS_LIBS)"

## controllability staircase form of descriptor state-space models
mkoctfile sltg01hd.cc \
          TG01HD.f TG01HX.f \
          "$(mkoctfile -p LAPACK_LIBS)" \
          "$(mkoctfile -p BLAS_LIBS)"

## observability staircase form of descriptor state-space models
mkoctfile sltg01id.cc \
          TG01ID.f TB01XD.f MA02CD.f AB07MD.f TG01HX.f \
          MA02BD.f \
          "$(mkoctfile -p LAPACK_LIBS)" \
          "$(mkoctfile -p BLAS_LIBS)"

system ("rm *.o");
cd (homedir);