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

## descriptor to regular state-space

## transfer function to state-space

system ("rm -rf *.o sltb04bd.oct slsb10jd.oct sltd04ad.oct");
system ("make sltb04bd.oct slsb10jd.oct sltd04ad.oct");

system ("rm -rf *.o");
cd (homedir);
