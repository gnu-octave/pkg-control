## ==============================================================================
## Developer Makefile for OCT-files
## ==============================================================================
## USAGE: * fetch control from Octave-Forge by svn
##        * add control/inst, control/src and control/devel to your Octave path
##        * run makefile_*
## ==============================================================================

homedir = pwd ();
develdir = fileparts (which ("makefile_minreal"));
srcdir = [develdir, "/../src"];
cd (srcdir);

mkoctfile sltb01pd.cc \
          TB01PD.f TB01XD.f TB01ID.f AB07MD.f TB01UD.f \
          MB03OY.f MB01PD.f MB01QD.f \
          "$(mkoctfile -p LAPACK_LIBS)" \
          "$(mkoctfile -p BLAS_LIBS)"

mkoctfile sltg01jd.cc \
          TG01JD.f TG01AD.f TB01XD.f MA02CD.f TG01HX.f \
          "$(mkoctfile -p LAPACK_LIBS)" \
          "$(mkoctfile -p BLAS_LIBS)"

system ("rm *.o");
cd (homedir);
