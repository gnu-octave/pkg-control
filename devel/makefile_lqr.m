## ==============================================================================
## Developer Makefile for OCT-files
## ==============================================================================
## USAGE: * fetch control from Octave-Forge by svn
##        * add control/inst, control/src and control/devel to your Octave path
##        * run makefile_*
## ==============================================================================

homedir = pwd ();
develdir = fileparts (which ("makefile_lqr"));
srcdir = [develdir, "/../src"];
cd (srcdir);

mkoctfile slsb02od.cc \
          SB02OD.f SB02OU.f SB02OV.f SB02OW.f SB02OY.f \
          SB02MR.f SB02MV.f \
          "$(mkoctfile -p BLAS_LIBS)"

mkoctfile slsg02ad.cc \
          SG02AD.f SB02OU.f SB02OV.f SB02OW.f SB02OY.f \
          MB01SD.f MB02VD.f MB02PD.f MA02GD.f \
          "$(mkoctfile -p BLAS_LIBS)"

system ("rm *.o");
cd (homedir);
