## ==============================================================================
## Developer Makefile for OCT-files
## ==============================================================================
## USAGE: * fetch control from Octave-Forge by svn
##        * add control/inst, control/src and control/devel to your Octave path
##        * run makefile_*
## ==============================================================================

homedir = pwd ();
develdir = fileparts (which ("makefile_hankel"));
srcdir = [develdir, "/../src"];
cd (srcdir);

mkoctfile slab13ad.cc \
          AB13AD.f TB01ID.f TB01KD.f AB13AX.f MA02DD.f \
          MB03UD.f TB01LD.f SB03OU.f MB03QX.f select.f \
          SB03OT.f MB03QD.f MB04ND.f MB04OD.f MB03QY.f \
          SB03OR.f SB03OY.f SB04PX.f MB04NY.f MB04OY.f \
          SB03OV.f

cd (homedir);
