## ==============================================================================
## Developer Makefile for OCT-files
## ==============================================================================
## USAGE: * fetch control from Octave-Forge by svn
##        * add control/inst, control/src and control/devel to your Octave path
##        * run makefile_*
## ==============================================================================

homedir = pwd ();
develdir = fileparts (which ("makefile_zero"));
srcdir = [develdir, "/../src"];
cd (srcdir);

## transmission zeros of state-space models
mkoctfile slab08nd.cc \
          AB08ND.f AB08NX.f TB01ID.f MB03OY.f MB03PY.f

## transmission zeros of descriptor state-space models
mkoctfile slag08bd.cc \
          AG08BD.f AG08BY.f TG01AD.f TB01XD.f MA02CD.f \
          TG01FD.f MA02BD.f MB03OY.f

## gain of descriptor state-space models
mkoctfile sltg04bx.cc \
          TG04BX.f MB02RD.f MB02SD.f

cd (homedir);
