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
develdir = fileparts (which ("makefile_ncfsyn"));
srcdir = [develdir, "/../src"];
cd (srcdir);

## H-infinity loop shaping - continuous-time
mkoctfile "-Wl,-framework" "-Wl,vecLib" \
          slsb10id.cc \
          SB10ID.f SB02RD.f select.f SB10JD.f MB02VD.f \
          MA02GD.f SB02MS.f MA02ED.f SB02RU.f SB02SD.f \
          MB01RU.f SB02QD.f SB02MV.f SB02MW.f SB02MR.f \
          MA02AD.f MB02PD.f MB01SD.f MB01UD.f SB03SY.f \
          MB01RX.f SB03MX.f SB03SX.f MB01RY.f SB03QY.f \
          SB03QX.f SB03MY.f SB04PX.f SB03MV.f SB03MW.f

## H-infinity loop shaping - discrete-time - strictly proper case
mkoctfile "-Wl,-framework" "-Wl,vecLib" \
          slsb10kd.cc \
          SB10KD.f SB02OD.f select.f SB02OY.f SB02OW.f \
          SB02OV.f SB02MV.f SB02OU.f SB02MR.f

## H-infinity loop shaping - discrete-time - proper case
mkoctfile "-Wl,-framework" "-Wl,vecLib" \
          slsb10zd.cc \
          SB10ZD.f MA02AD.f SB02OD.f select.f MB01RX.f \
          MB02VD.f SB02OY.f SB02OW.f SB02OV.f SB02OU.f \
          SB02MR.f MA02GD.f SB02MV.f

cd (homedir);
