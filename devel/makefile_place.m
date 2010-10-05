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
develdir = fileparts (which ("makefile_place"));
srcdir = [develdir, "/../src"];
cd (srcdir);

mkoctfile "-Wl,-framework" "-Wl,vecLib" \
          slsb01bd.cc \
          SB01BD.f MB03QD.f MB03QY.f SB01BX.f SB01BY.f \
          select.f

cd (homedir);