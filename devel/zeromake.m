homedir = pwd ();
develdir = fileparts (which ("makefile_zero"));
srcdir = [develdir, "/../src"];
cd (srcdir);

## gain of descriptor state-space models
mkoctfile "-Wl,-framework" "-Wl,vecLib" \
          sltg04bx.cc \
          TG04BX.f MB02RD.f MB02SD.f