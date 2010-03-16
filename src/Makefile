all: slab08nd.oct slab13dd.oct slsb10hd.oct slsb10ed.oct slab13bd.oct \
     slsb01bd.oct slsb10fd.oct slsb10dd.oct slsb03md.oct slsb04md.oct \
     slsb04qd.oct slsg03ad.oct slsb02od.oct slab13ad.oct

# transmission zeros of state-space models
slab08nd.oct: slab08nd.cc
	mkoctfile slab08nd.cc \
              AB08ND.f AB08NX.f TB01ID.f MB03OY.f MB03PY.f

# L-infinity norm
slab13dd.oct: slab13dd.cc
	mkoctfile slab13dd.cc \
              AB13DD.f MA02AD.f MB01SD.f MB03XD.f TB01ID.f \
              TG01AD.f TG01BD.f AB13DX.f MA01AD.f MA02ID.f \
              MB03XP.f MB04DD.f MB04QB.f MB04TB.f MB03XU.f \
              MB04TS.f UE01MD.f MB02RD.f MB02SD.f MB04QC.f \
              MB04QF.f MB03YA.f MB03YD.f MB02RZ.f MB04QU.f \
              MB02SZ.f MB03YT.f

# H-2 controller synthesis - continuous-time
slsb10hd.oct: slsb10hd.cc
	mkoctfile slsb10hd.cc \
              SB10HD.f SB10UD.f SB10VD.f SB10WD.f SB02RD.f \
              MB01RU.f SB02SD.f MA02ED.f SB02RU.f SB02MR.f \
              MB01SD.f SB02MS.f SB02MV.f SB02MW.f MA02AD.f \
              SB02QD.f MB02PD.f SB03QX.f SB03QY.f MB01RX.f \
              MB01RY.f SB03SX.f SB03SY.f select.f SB03MX.f \
              SB03MY.f MB01UD.f SB03MV.f SB03MW.f SB04PX.f

# H-2 controller synthesis - discrete-time
slsb10ed.oct: slsb10ed.cc
	mkoctfile slsb10ed.cc \
              SB10ED.f SB10SD.f SB10TD.f SB10PD.f MB01RX.f \
              SB02SD.f SB02OD.f MB01RU.f SB02OU.f SB02OV.f \
              SB02OW.f MB01RY.f SB02OY.f SB03SX.f SB03SY.f \
              MA02ED.f select.f SB03MX.f SB02MR.f SB02MV.f \
              MB01UD.f SB03MV.f SB04PX.f

# H-2 norm
slab13bd.oct: slab13bd.cc
	mkoctfile slab13bd.cc \
              AB13BD.f SB08DD.f SB03OU.f SB01FY.f TB01LD.f \
              SB03OT.f MB04ND.f MB04OD.f MB03QX.f select.f \
              SB03OR.f MB04OX.f MB03QD.f SB03OY.f MA02AD.f \
              MB03QY.f SB04PX.f MB04NY.f MB04OY.f SB03OV.f

# Pole assignment
slsb01bd.oct: slsb01bd.cc
	mkoctfile slsb01bd.cc \
              SB01BD.f MB03QD.f MB03QY.f SB01BX.f SB01BY.f \
              select.f

# H-inf controller synthesis - continuous-time
slsb10fd.oct: slsb10fd.cc
	mkoctfile slsb10fd.cc \
              SB10FD.f SB10PD.f SB10QD.f SB10RD.f SB02RD.f \
              MB01RU.f MB01RX.f MA02AD.f SB02SD.f MA02ED.f \
              SB02RU.f SB02MR.f MB01SD.f SB02MS.f SB02MV.f \
              SB02MW.f SB02QD.f MB02PD.f SB03QX.f SB03QY.f \
              MB01RY.f SB03SX.f SB03SY.f select.f SB03MX.f \
              SB03MY.f MB01UD.f SB03MV.f SB03MW.f SB04PX.f

# H-inf controller synthesis - discrete-time
slsb10dd.oct: slsb10dd.cc
	mkoctfile slsb10dd.cc \
              SB10DD.f MB01RU.f MB01RX.f SB02SD.f SB02OD.f \
              MA02AD.f SB02OU.f SB02OV.f SB02OW.f MB01RY.f \
              SB02OY.f SB03SX.f SB03SY.f MA02ED.f select.f \
              SB03MX.f SB02MR.f SB02MV.f MB01UD.f SB03MV.f \
              SB04PX.f

# Lyapunov equations
slsb03md.oct: slsb03md.cc
	mkoctfile slsb03md.cc \
              SB03MD.f select.f SB03MX.f SB03MY.f MB01RD.f \
              SB03MV.f SB03MW.f SB04PX.f

# Sylvester equations - continuous-time
slsb04md.oct: slsb04md.cc
	mkoctfile slsb04md.cc \
              SB04MD.f SB04MU.f SB04MY.f SB04MR.f SB04MW.f

# Sylvester equations - discrete-time
slsb04qd.oct: slsb04qd.cc
	mkoctfile slsb04qd.cc \
              SB04QD.f SB04QU.f SB04QY.f SB04MW.f SB04QR.f

# generalized Lyapunov equations
slsg03ad.oct: slsg03ad.cc
	mkoctfile slsg03ad.cc \
              SG03AD.f MB01RW.f MB01RD.f SG03AX.f SG03AY.f \
              MB02UU.f MB02UV.f

# algebraic Riccati equations
slsb02od.oct: slsb02od.cc
	mkoctfile slsb02od.cc \
              SB02OD.f SB02OU.f SB02OV.f SB02OW.f SB02OY.f \
              SB02MR.f SB02MV.f

# Hankel singular values
slab13ad.oct: slab13ad.cc
	mkoctfile slab13ad.cc \
              AB13AD.f TB01ID.f TB01KD.f AB13AX.f MA02DD.f \
              MB03UD.f TB01LD.f SB03OU.f MB03QX.f select.f \
              SB03OT.f MB03QD.f MB04ND.f MB04OD.f MB03QY.f \
              SB03OR.f SB03OY.f SB04PX.f MB04NY.f MB04OY.f \
              SB03OV.f

clean:
	rm *.o core octave-core *.oct *~