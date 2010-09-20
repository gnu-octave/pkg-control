## ==============================================================================
## Developer Makefile for OCT-files
## ==============================================================================
## USAGE: * fetch control from Octave-Forge by svn
##        * add control/inst and control/src to your Octave path
##        * run makefile_*
## ==============================================================================

homedir = pwd ();
srcdir = fileparts (which ("makefile_zero"));
cd (srcdir);

## transmission zeros of state-space models
mkoctfile "-Wl,-framework" "-Wl,vecLib" \
          slab08nd.cc \
          AB08ND.f AB08NX.f TB01ID.f MB03OY.f MB03PY.f

## transmission zeros of descriptor state-space models
mkoctfile "-Wl,-framework" "-Wl,vecLib" \
          slag08bd.cc \
          AG08BD.f AG08BY.f TG01AD.f TB01XD.f MA02CD.f \
          TG01FD.f MA02BD.f MB03OY.f

cd (homedir);