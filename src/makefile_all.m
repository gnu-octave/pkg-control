## ==============================================================================
## Developer Makefile for OCT-files
## ==============================================================================
## USAGE: * fetch control from Octave-Forge by svn
##        * add control/inst and control/src to your Octave path
##        * run makefile_*
## ==============================================================================

homedir = pwd ();
srcdir = fileparts (which ("makefile_all"));
cd (srcdir);

makefile_chol
makefile_h2syn
makefile_hankel
makefile_hinfsyn
makefile_lqr
makefile_lyap
makefile_minreal
makefile_norm
makefile_place
makefile_staircase
makefile_zero

cd (homedir);