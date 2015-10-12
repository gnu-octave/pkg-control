A = [ -0.01357      -32.2          -46.3            0
       0.00012        0              1.214          0
      -0.0001212      0             -1.214          1
       0.00057        0             -9.1           -0.6696  ];

B = [ -0.433;
       0.1394;
      -0.1394;
      -0.1577   ];

C = [  0              0              0              1
       1              0              0              0       ];

D = [  0
       0        ];

x0 = [ 2
       2
       2
       2        ];
 
sys = ss (A, B, C, D);
sys = set (sys, 'inputname', 'u', 'outputname', 'y')


lambda_sys_desired = [ -0.5   -1+1j   -1-1j   -2  ]

F = place (sys, lambda_sys_desired)

ctrl = ss (-F);
ctrl = set (ctrl, 'inputname', 'xhat', 'outputname', 'u')


lambda_obs_desired = [ -10    -11     -12     -13 ]

K = place (sys.', lambda_obs_desired).'

obs = ss (A-K*C, [B, K]);
obs = set (obs, 'inputname', {'u', 'y1', 'y2'}, 'outputname', 'xhat')


G = ss (C);
G = set (G, 'inputname', 'xhat', 'outputname', 'yhat')


OE = sumblk ('e = y - yhat', 2)


N = connect (sys, ctrl, obs, G, OE, ':', {'e1', 'e2'})


y0 = C * x0;
x0hat = pinv (C.' * C) * C.' * y0;


figure (1)
initial (N, [x0; x0hat], 1)


