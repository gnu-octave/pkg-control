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

%{
makefile_chol
makefile_conversions
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
makefile_tustin
makefile_zero
%}

system ("make realclean");
system ("make -j4 all");
system ("rm *.o");

cd (homedir);
