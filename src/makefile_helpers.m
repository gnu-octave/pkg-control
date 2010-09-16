## ==============================================================================
## Developer Makefile for OCT-files
## ==============================================================================
## USAGE: * fetch control from Octave-Forge by svn
##        * add control/inst and control/src to your Octave path
##        * run makefile_*
## ==============================================================================

homedir = pwd ();
srcdir = fileparts (which ("makefile_helpers"));
cd (srcdir);

mkoctfile is_real_scalar.cc

mkoctfile is_real_vector.cc

mkoctfile is_real_matrix.cc

mkoctfile is_real_square_matrix.cc

cd (homedir);