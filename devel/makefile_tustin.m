## ==============================================================================
## Developer Makefile for OCT-files
## ==============================================================================
## USAGE: * fetch control from Octave-Forge by svn
##        * add control/inst, control/src and control/devel to your Octave path
##        * run makefile_*
## ==============================================================================

homedir = pwd ();
develdir = fileparts (which ("makefile_tustin"));
srcdir = [develdir, "/../src"];
cd (srcdir);

## bilinear transformation
mkoctfile slab04md.cc \
          AB04MD.f

system ("rm *.o");
cd (homedir);
