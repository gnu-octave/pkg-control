% Plant
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
 
sys = ss (A, B, [C; eye(4)], [D; zeros(4,1)]);
sys = set (sys, 'inputname', 'u', 'outputname', {'y1', 'y2', 'x1', 'x2', 'x3', 'x4'})


% Controller
lambda_sys_desired = [ -0.5   -1+1j   -1-1j   -2  ]

F = place (sys, lambda_sys_desired)

ctrl = ss (-F);
ctrl = set (ctrl, 'inputname', 'xhat', 'outputname', 'u')


% Observer
lambda_obs_desired = [ -10    -11     -12     -13 ]

K = place (sys({'y1', 'y2'}, :).', lambda_obs_desired).'

obs = ss (A-K*C, [B, K]);
obs = set (obs, 'inputname', {'u', 'y1', 'y2'}, 'outputname', 'xhat')


% Output Error
OE = sumblk ('e = x - xhat', 4)


% Entire System
N = connect (sys, ctrl, obs, OE, ':', {'e1', 'e2', 'e3', 'e4'})


% Initial Conditions
y0 = C * x0;
x0hat = pinv (C.' * C) * C.' * y0;


% Simulation
figure (1)
initial (N, [x0; x0hat], 1)

[Y, T] = initial (N, [x0; x0hat], 1);

figure (2)
plot (T, Y)
grid on

