sys = tf (WestlandLynx);

w = logspace (-3,4,1e4);

tic
H = freqresp (sys(1,1), w);
toc

%sigma (sys, logspace (-2,2, 20))
tic
sigma (sys, w)
toc