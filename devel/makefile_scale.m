## ==============================================================================
## Developer Makefile for OCT-files
## ==============================================================================
## USAGE: * fetch control from Octave-Forge by svn
##        * add control/inst, control/src and control/devel to your Octave path
##        * run makefile_*
## ==============================================================================

homedir = pwd ();
develdir = fileparts (which ("makefile_scale"));
srcdir = [develdir, "/../src"];
cd (srcdir);

## scaling of state-space models
mkoctfile sltb01id.cc \
          TB01ID.f

## scaling of descriptor state-space models
mkoctfile sltg01ad.cc \
          TG01AD.f

system ("rm *.o");
cd (homedir);
