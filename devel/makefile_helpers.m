## ==============================================================================
## Developer Makefile for OCT-files
## ==============================================================================
## USAGE: * fetch control from Octave-Forge by svn
##        * add control/inst, control/src and control/devel to your Octave path
##        * run makefile_*
## ==============================================================================

homedir = pwd ();
develdir = fileparts (which ("makefile_helpers"));
srcdir = [develdir, "/../src"];
cd (srcdir);

mkoctfile is_real_scalar.cc

mkoctfile is_real_vector.cc

mkoctfile is_real_matrix.cc

mkoctfile is_real_square_matrix.cc

system ("rm *.o");
cd (homedir);