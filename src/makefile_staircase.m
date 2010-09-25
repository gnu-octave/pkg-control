## ==============================================================================
## Developer Makefile for OCT-files
## ==============================================================================
## USAGE: * fetch control from Octave-Forge by svn
##        * add control/inst and control/src to your Octave path
##        * run makefile_*
## ==============================================================================

homedir = pwd ();
srcdir = fileparts (which ("makefile_staircase"));
cd (srcdir);

## staircase form using orthogonal transformations
mkoctfile "-Wl,-framework" "-Wl,vecLib" \
          slab01od.cc \
          AB01OD.f AB01ND.f MB03OY.f MB01PD.f MB01QD.f

## controllability staircase form of descriptor state-space models
mkoctfile "-Wl,-framework" "-Wl,vecLib" \
          sltg01hd.cc \
          TG01HD.f TG01HX.f

## observability staircase form of descriptor state-space models
mkoctfile "-Wl,-framework" "-Wl,vecLib" \
          sltg01id.cc \
          TG01ID.f TB01XD.f MA02CD.f AB07MD.f TG01HX.f \
          MA02BD.f

cd (homedir);