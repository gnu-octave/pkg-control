## ==============================================================================
## Developer Makefile for OCT-files
## ==============================================================================
## USAGE: * fetch control-devel from Octave-Forge by svn
##        * add control-devel/inst, control-devel/src and control-devel/devel
##          to your Octave path (by an .octaverc file)
##        * run makefile_modred
## ==============================================================================

homedir = pwd ();
develdir = fileparts (which ("makefile_modred"));
srcdir = [develdir, "/../src"];
cd (srcdir);

mkoctfile slab09hd.cc \
          slicotlibrary.a \
          "$(mkoctfile -p LAPACK_LIBS)" \
          "$(mkoctfile -p BLAS_LIBS)"

mkoctfile slab09id.cc \
          slicotlibrary.a \
          "$(mkoctfile -p LAPACK_LIBS)" \
          "$(mkoctfile -p BLAS_LIBS)"

mkoctfile slab09jd.cc \
          slicotlibrary.a \
          "$(mkoctfile -p LAPACK_LIBS)" \
          "$(mkoctfile -p BLAS_LIBS)"

system ("rm *.o");
cd (homedir);

