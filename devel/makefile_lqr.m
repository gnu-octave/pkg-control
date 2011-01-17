## ==============================================================================
## Developer Makefile for OCT-files
## ==============================================================================
## USAGE: * fetch control from Octave-Forge by svn
##        * add control/inst, control/src and control/devel to your Octave path
##        * run makefile_*
## NOTES: * The option "-Wl,-framework" "-Wl,vecLib" is needed for MacPorts'
##          octave-devel @3.3.52_1+gcc44 on MacOS X 10.6.4. However, this option
##          breaks other platforms. See MacPorts Ticket #26640.
## ==============================================================================

homedir = pwd ();
develdir = fileparts (which ("makefile_lqr"));
srcdir = [develdir, "/../src"];
cd (srcdir);

mkoctfile "-Wl,-framework" "-Wl,vecLib" \
          slsb02od.cc \
          SB02OD.f SB02OU.f SB02OV.f SB02OW.f SB02OY.f \
          SB02MR.f SB02MV.f

mkoctfile "-Wl,-framework" "-Wl,vecLib" \
          slsg02ad.cc \
          SG02AD.f SB02OU.f SB02OV.f SB02OW.f SB02OY.f \
          MB01SD.f MB02VD.f MB02PD.f MA02GD.f

cd (homedir);