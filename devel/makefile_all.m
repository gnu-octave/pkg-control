## ==============================================================================
## Developer Makefile for OCT-files
## ==============================================================================
## USAGE: * fetch control from Octave-Forge by svn
##        * add control/inst, control/src and control/devel to your Octave path
##        * run makefile_*
## ==============================================================================

homedir = pwd ();
develdir = fileparts (which ("makefile_all"));
srcdir = [develdir, "/../src"];
cd (srcdir);

makefile_chol
makefile_h2syn
makefile_hankel
makefile_helpers
makefile_hinfsyn
makefile_lqr
makefile_lyap
makefile_minreal
makefile_ncfsyn
makefile_norm
makefile_place
makefile_scale
makefile_staircase
makefile_zero

cd (homedir);