## ==============================================================================
## Developer Makefile for OCT-files
## ==============================================================================
## USAGE: * fetch control from Octave-Forge by svn
##        * add control/inst, control/src and control/devel to your Octave path
##        * run makefile_*
## ==============================================================================

homedir = pwd ();
develdir = fileparts (which ("makefile_lyap"));
srcdir = [develdir, "/../src"];
cd (srcdir);

## Lypunov
mkoctfile slsb03md.cc \
          SB03MD.f select.f SB03MX.f SB03MY.f MB01RD.f \
          SB03MV.f SB03MW.f SB04PX.f \
          "$(mkoctfile -p BLAS_LIBS)"


## Sylvester          
mkoctfile slsb04md.cc \
          SB04MD.f SB04MU.f SB04MY.f SB04MR.f SB04MW.f \
          "$(mkoctfile -p BLAS_LIBS)"

mkoctfile slsb04qd.cc \
          SB04QD.f SB04QU.f SB04QY.f SB04MW.f SB04QR.f \
          "$(mkoctfile -p BLAS_LIBS)"


## Generalized Lyapunov
mkoctfile slsg03ad.cc \
          SG03AD.f MB01RW.f MB01RD.f SG03AX.f SG03AY.f \
          MB02UU.f MB02UV.f \
          "$(mkoctfile -p BLAS_LIBS)"

system ("rm *.o");
cd (homedir);
