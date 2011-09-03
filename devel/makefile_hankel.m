## ==============================================================================
## Developer Makefile for OCT-files
## ==============================================================================
## USAGE: * fetch control from Octave-Forge by svn
##        * add control/inst, control/src and control/devel to your Octave path
##        * run makefile_*
## NOTES: * The option "-Wl,-framework" "-Wl,vecLib" is needed for MacPorts'
##          octave-devel @3.3.52_1+gcc44 on MacOS X 10.6.4. However, this option
##          breaks other platforms. See MacPorts Ticket #26640.
## ==============================================================================

homedir = pwd ();
develdir = fileparts (which ("makefile_hankel"));
srcdir = [develdir, "/../src"];
cd (srcdir);

mkoctfile "-Wl,-framework" "-Wl,vecLib" \
          slab13ad.cc \
          AB13AD.f TB01ID.f TB01KD.f AB13AX.f MA02DD.f \
          MB03UD.f TB01LD.f SB03OU.f MB03QX.f select.f \
          SB03OT.f MB03QD.f MB04ND.f MB04OD.f MB03QY.f \
          SB03OR.f SB03OY.f SB04PX.f MB04NY.f MB04OY.f \
          SB03OV.f

system ("rm *.o");
cd (homedir);
