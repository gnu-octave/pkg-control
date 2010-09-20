## ==============================================================================
## Developer Makefile for OCT-files
## ==============================================================================
## USAGE: * fetch control from Octave-Forge by svn
##        * add control/inst and control/src to your Octave path
##        * run makefile_*
## ==============================================================================

homedir = pwd ();
srcdir = fileparts (which ("makefile_lqr"));
cd (srcdir);

mkoctfile "-Wl,-framework" "-Wl,vecLib" \
          slsb02od.cc \
          SB02OD.f SB02OU.f SB02OV.f SB02OW.f SB02OY.f \
          SB02MR.f SB02MV.f

cd (homedir);