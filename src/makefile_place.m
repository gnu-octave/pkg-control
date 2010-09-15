## ==============================================================================
## Developer Makefile for OCT-files
## ==============================================================================
## USAGE: * fetch control from Octave-Forge by svn
##        * add control/inst and control/src to your Octave path
##        * run makefile_*
## ==============================================================================

homedir = pwd ();
srcdir = fileparts (which ("makefile_place"));
cd (srcdir);

mkoctfile slsb01bd.cc \
          SB01BD.f MB03QD.f MB03QY.f SB01BX.f SB01BY.f \
          select.f

cd (homedir);