## ==============================================================================
## Developer Makefile for OCT-files
## ==============================================================================
## USAGE: * fetch control from Octave-Forge by svn
##        * add control/inst, control/src and control/devel to your Octave path
##        * run makefile_*
## ==============================================================================

homedir = pwd ();
develdir = fileparts (which ("makefile_chol"));
srcdir = [develdir, "/../src"];
cd (srcdir);

mkoctfile slsb03od.cc \
          SB03OD.f select.f SB03OU.f SB03OT.f MB04ND.f \
          MB04OD.f SB03OR.f SB03OY.f SB04PX.f MB04NY.f \
          MB04OY.f SB03OV.f
   

mkoctfile slsg03bd.cc \
          SG03BD.f SG03BV.f SG03BU.f SG03BW.f SG03BX.f \
          SG03BY.f MB02UU.f MB02UV.f

cd (homedir);