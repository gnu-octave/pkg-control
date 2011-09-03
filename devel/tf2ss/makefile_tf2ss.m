mkoctfile "-Wl,-framework" "-Wl,vecLib" \
          sltd04ad.cc \
          TD04AD.f TD03AY.f TB01PD.f TB01XD.f AB07MD.f \
          TB01UD.f TB01ID.f MB01PD.f MB03OY.f MB01QD.f

system ("rm *.o");
