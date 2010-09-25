## ==============================================================================
## Developer Makefile for OCT-files
## ==============================================================================
## USAGE: * fetch control from Octave-Forge by svn
##        * add control/inst and control/src to your Octave path
##        * run makefile_*
## ==============================================================================

homedir = pwd ();
srcdir = fileparts (which ("makefile_minreal"));
cd (srcdir);

mkoctfile "-Wl,-framework" "-Wl,vecLib" \
          sltb01pd.cc \
          TB01PD.f TB01XD.f TB01ID.f AB07MD.f TB01UD.f \
          MB03OY.f MB01PD.f MB01QD.f

mkoctfile "-Wl,-framework" "-Wl,vecLib" \
          sltg01jd.cc \
          TG01JD.f TG01AD.f TB01XD.f MA02CD.f TG01HX.f

cd (homedir);