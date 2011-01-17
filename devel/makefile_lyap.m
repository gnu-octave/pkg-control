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
develdir = fileparts (which ("makefile_lyap"));
srcdir = [develdir, "/../src"];
cd (srcdir);

## Lypunov
mkoctfile "-Wl,-framework" "-Wl,vecLib" \
          slsb03md.cc \
          SB03MD.f select.f SB03MX.f SB03MY.f MB01RD.f \
          SB03MV.f SB03MW.f SB04PX.f


## Sylvester          
mkoctfile "-Wl,-framework" "-Wl,vecLib" \
          slsb04md.cc \
          SB04MD.f SB04MU.f SB04MY.f SB04MR.f SB04MW.f

mkoctfile "-Wl,-framework" "-Wl,vecLib" \
          slsb04qd.cc \
          SB04QD.f SB04QU.f SB04QY.f SB04MW.f SB04QR.f


## Generalized Lyapunov
mkoctfile "-Wl,-framework" "-Wl,vecLib" \
          slsg03ad.cc \
          SG03AD.f MB01RW.f MB01RD.f SG03AX.f SG03AY.f \
          MB02UU.f MB02UV.f

cd (homedir);