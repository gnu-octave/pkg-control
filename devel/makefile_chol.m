## ==============================================================================
## Developer Makefile for OCT-files
## ==============================================================================
## USAGE: * fetch control from Octave-Forge by svn
##        * add control/inst, control/src and control/devel to your Octave path
##        * run makefile_*
## ==============================================================================

homedir = pwd ();
develdir = fileparts (which ("makefile_chol"));
srcdir = [develdir, "/../src"];
cd (srcdir);

system ("rm -rf *.o slsb03od.oct slsg03bd.oct");
system ("make slsb03od.oct slsg03bd.oct");

system ("rm -rf *.o");
cd (homedir);

%{
system ("rm -rf *.o ");
system ("make ");

system ("rm -rf *.o");
%}