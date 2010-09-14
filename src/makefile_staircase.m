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

mkoctfile slab01od.cc \
          AB01OD.f AB01ND.f MB03OY.f MB01PD.f MB01QD.f

cd (homedir);