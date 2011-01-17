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
develdir = fileparts (which ("makefile_minreal"));
srcdir = [develdir, "/../src"];
cd (srcdir);

mkoctfile "-Wl,-framework" "-Wl,vecLib" \
          sltb01pd.cc \
          TB01PD.f TB01XD.f TB01ID.f AB07MD.f TB01UD.f \
          MB03OY.f MB01PD.f MB01QD.f

mkoctfile "-Wl,-framework" "-Wl,vecLib" \
          sltg01jd.cc \
          TG01JD.f TG01AD.f TB01XD.f MA02CD.f TG01HX.f

cd (homedir);