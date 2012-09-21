## ==============================================================================
## Developer Makefile for OCT-files
## ==============================================================================
## USAGE: * fetch control from Octave-Forge by svn
##        * add control/inst, control/src and control/devel to your Octave path
##        * run makefile_control
## ==============================================================================

homedir = pwd ();
develdir = fileparts (which ("makefile_control"));
srcdir = [develdir, "/../src"];
cd (srcdir);

## system ("make realclean");  # recompile slicotlibrary.a
system ("make clean");
system ("make -j1 all");
system ("rm *.o");
system ("rm *.d");

cd (homedir);
