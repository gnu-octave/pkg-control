## ==============================================================================
## Developer Makefile for OCT-files
## ==============================================================================
## USAGE: * fetch control from Octave-Forge by svn
##        * add control/inst, control/src and control/devel to your Octave path
##        * run makefile_*
## ==============================================================================

homedir = pwd ();
develdir = fileparts (which ("makefile_conversions"));
srcdir = [develdir, "/../src"];
cd (srcdir);

## state-space to transfer function
mkoctfile sltb04bd.cc \
          TB04BD.f MC01PY.f TB01ID.f TB01ZD.f MC01PD.f \
          TB04BX.f MA02AD.f MB02RD.f MB01PD.f MB02SD.f \
          MB01QD.f

## descriptor to regular state-space
mkoctfile slsb10jd.cc \
          SB10JD.f

## transfer function to state-space
mkoctfile sltd04ad.cc \
          TD04AD.f TD03AY.f TB01PD.f TB01XD.f AB07MD.f \
          TB01UD.f TB01ID.f MB01PD.f MB03OY.f MB01QD.f

system ("rm *.o");
cd (homedir);
