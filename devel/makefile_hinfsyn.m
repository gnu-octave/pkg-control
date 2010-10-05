## ==============================================================================
## Developer Makefile for OCT-files
## ==============================================================================
## USAGE: * fetch control from Octave-Forge by svn
##        * add control/inst, control/src and control/devel to your Octave path
##        * run makefile_*
## ==============================================================================

homedir = pwd ();
develdir = fileparts (which ("makefile_hinfsyn"));
srcdir = [develdir, "/../src"];
cd (srcdir);

mkoctfile "-Wl,-framework" "-Wl,vecLib" \
          slsb10fd.cc \
          SB10FD.f SB10PD.f SB10QD.f SB10RD.f SB02RD.f \
          MB01RU.f MB01RX.f MA02AD.f SB02SD.f MA02ED.f \
          SB02RU.f SB02MR.f MB01SD.f SB02MS.f SB02MV.f \
          SB02MW.f SB02QD.f MB02PD.f SB03QX.f SB03QY.f \
          MB01RY.f SB03SX.f SB03SY.f select.f SB03MX.f \
          SB03MY.f MB01UD.f SB03MV.f SB03MW.f SB04PX.f

mkoctfile "-Wl,-framework" "-Wl,vecLib" \
          slsb10dd.cc \
          SB10DD.f MB01RU.f MB01RX.f SB02SD.f SB02OD.f \
          MA02AD.f SB02OU.f SB02OV.f SB02OW.f MB01RY.f \
          SB02OY.f SB03SX.f SB03SY.f MA02ED.f select.f \
          SB03MX.f SB02MR.f SB02MV.f MB01UD.f SB03MV.f \
          SB04PX.f

cd (homedir);