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
develdir = fileparts (which ("makefile_staircase"));
srcdir = [develdir, "/../src"];
cd (srcdir);

## staircase form using orthogonal transformations
mkoctfile "-Wl,-framework" "-Wl,vecLib" \
          slab01od.cc \
          AB01OD.f AB01ND.f MB03OY.f MB01PD.f MB01QD.f

## controllability staircase form of descriptor state-space models
mkoctfile "-Wl,-framework" "-Wl,vecLib" \
          sltg01hd.cc \
          TG01HD.f TG01HX.f

## observability staircase form of descriptor state-space models
mkoctfile "-Wl,-framework" "-Wl,vecLib" \
          sltg01id.cc \
          TG01ID.f TB01XD.f MA02CD.f AB07MD.f TG01HX.f \
          MA02BD.f

system ("rm *.o");
cd (homedir);