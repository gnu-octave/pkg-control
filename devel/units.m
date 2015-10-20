n = 3;  p = 4;  m = 2;
sys = drss (n, p, m)

stn = {'power', 'voltage', 'current'};
stu = {'VA', 'V', 'A'};

outn = {'energy', 'mass', 'force', 'frequency'};
outu = {'kJ', 'kg', 'N', 'Hz'};

inn = {'angle', 'length'};
inu = {'radians', 'm'};

set (sys, 'ts', 1, 'timeunit', 'seconds', ...
          'statename', stn, 'stateunit', stu, ...
          'outputname', outn, 'outputunit', outu, ...
          'inputname', inn, 'inputunit', inu);

sys_ss = sys
sys_tf = tf (sys)
sys_zpk = zpk (sys)
sys_frd = frd (sys, logspace (-3, -1, 10))
