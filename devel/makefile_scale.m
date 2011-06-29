homedir = pwd ();
develdir = fileparts (which ("makefile_scale"));
srcdir = [develdir, "/../src"];
cd (srcdir);

## scaling of state-space models
mkoctfile "-Wl,-framework" "-Wl,vecLib" \
          sltb01id.cc \
          TB01ID.f

## scaling of descriptor state-space models
mkoctfile "-Wl,-framework" "-Wl,vecLib" \
          sltg01ad.cc \
          TG01AD.f

cd (homedir);
