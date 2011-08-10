## ==============================================================================
## Developer Makefile for OCT-files
## ==============================================================================
## USAGE: * fetch control from Octave-Forge by svn
##        * add control/inst, control/src and control/devel to your Octave path
##        * run makefile_*
## ==============================================================================

homedir = pwd ();
develdir = fileparts (which ("makefile_norm"));
srcdir = [develdir, "/../src"];
cd (srcdir);

## H-2 norm
mkoctfile slab13bd.cc \
          AB13BD.f SB08DD.f SB03OU.f SB01FY.f TB01LD.f \
          SB03OT.f MB04ND.f MB04OD.f MB03QX.f select.f \
          SB03OR.f MB04OX.f MB03QD.f SB03OY.f MA02AD.f \
          MB03QY.f SB04PX.f MB04NY.f MB04OY.f SB03OV.f

## L-inf norm
mkoctfile slab13dd.cc \
          AB13DD.f MA02AD.f MB01SD.f MB03XD.f TB01ID.f \
          TG01AD.f TG01BD.f AB13DX.f MA01AD.f MA02ID.f \
          MB03XP.f MB04DD.f MB04QB.f MB04TB.f MB03XU.f \
          MB04TS.f UE01MD.f MB02RD.f MB02SD.f MB04QC.f \
          MB04QF.f MB03YA.f MB03YD.f MB02RZ.f MB04QU.f \
          MB02SZ.f MB03YT.f

cd (homedir);
