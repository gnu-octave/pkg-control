## ==============================================================================
## Developer Makefile for OCT-files
## ==============================================================================
## USAGE: * fetch control-devel from Octave-Forge by svn
##        * add control-devel/inst, control-devel/src and control-devel/devel
##          to your Octave path (by an .octaverc file)
##        * run makefile_conred
## ==============================================================================

homedir = pwd ();
develdir = fileparts (which ("makefile_conred"));
srcdir = [develdir, "/../src"];
cd (srcdir);

mkoctfile slsb16ad.cc \
          SB16AD.f TB01ID.f SB16AY.f TB01KD.f AB09IX.f \
          MB04OD.f MB01WD.f SB03OD.f MB03UD.f AB05PD.f \
          AB09DD.f AB07ND.f TB01LD.f AB05QD.f SB03OU.f \
          MA02AD.f MB03QX.f select.f MB01YD.f MB01ZD.f \
          SB03OT.f MB04OY.f MB03QD.f MB04ND.f MB03QY.f \
          SB03OR.f SB03OY.f SB04PX.f MB04NY.f SB03OV.f \
          "$(mkoctfile -p LAPACK_LIBS)" \
          "$(mkoctfile -p BLAS_LIBS)"

mkoctfile slsb16bd.cc \
          SB16BD.f AB09AD.f AB09BD.f SB08GD.f SB08HD.f \
          TB01ID.f AB09AX.f MA02GD.f AB09BX.f TB01WD.f \
          MA02DD.f MB03UD.f select.f AB09DD.f SB03OU.f \
          MA02AD.f SB03OT.f MB04ND.f MB04OD.f SB03OR.f \
          SB03OY.f SB04PX.f MB04NY.f MB04OY.f SB03OV.f \
          "$(mkoctfile -p LAPACK_LIBS)" \
          "$(mkoctfile -p BLAS_LIBS)"

mkoctfile slsb16cd.cc \
          SB16CD.f SB16CY.f AB09IX.f SB03OD.f MB02UD.f \
          AB09DD.f MA02AD.f MB03UD.f select.f SB03OU.f \
          MB01SD.f SB03OT.f MB04ND.f MB04OD.f SB03OR.f \
          SB03OY.f SB04PX.f MB04NY.f MB04OY.f SB03OV.f \
          "$(mkoctfile -p LAPACK_LIBS)" \
          "$(mkoctfile -p BLAS_LIBS)"

system ("rm *.o");
cd (homedir);

