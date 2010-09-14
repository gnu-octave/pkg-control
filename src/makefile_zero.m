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

mkoctfile slab08nd.cc \
          AB08ND.f AB08NX.f TB01ID.f MB03OY.f MB03PY.f

cd (homedir);