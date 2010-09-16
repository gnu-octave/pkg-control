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

mkoctfile is_ss_mat.cc

mkoctfile is_real_scalar.cc

cd (homedir);