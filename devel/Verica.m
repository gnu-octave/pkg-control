%% -*- texinfo -*-
%% Figure 7: Full-order observer estimation errors e(t) = x(t) - xhat(t) for
%% the aircraft example.
%% FIXME: The figure doesn't match the published one yet.

% ===============================================================================
% Aircraft Simulation Example             Lukas Reichlin             October 2015
% ===============================================================================
% Reference:
% Verica Radisavljevic-Gajic:
% Full- and Reduced-Order Linear Observer Implementations in Matlab/Simulink.
% IEEE Control Systems Magazine, Vol. 35, No. 5, pp. 91-101, October 2015
% Digital Object Identifier 10.1109/MCS.2015.2449691
% ===============================================================================


% Tabula Rasa
clear all, close all, clc


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
 
sys = ss (A, B, C, D);
sys = set (sys, 'inputname', 'u', 'outputname', 'y', 'statename', 'x')


% Controller
lambda_sys_desired = [ -0.5   -1+1j   -1-1j   -2  ]

F = place (sys, lambda_sys_desired)

ctrl = ss (-F);
ctrl = set (ctrl, 'inputname', 'xhat', 'outputname', 'u')


% Observer
lambda_obs_desired = [ -10    -11     -12     -13 ]

K = place (sys.', lambda_obs_desired).'

obs = ss (A-K*C, [B, K]);
obs = set (obs, 'inputname', {'u', 'y1', 'y2'}, 'outputname', 'xhat', 'statename', 'xhat')


% Entire System
N = connect (sys, ctrl, obs, 'u', {'y1', 'y2'})


% Initial Conditions
y0 = C * x0;
x0hat = pinv (C.' * C) * C.' * y0;


% Simulation
figure (1)
[Y, T, X] = initial (N, [x0; x0hat], 1);
ERR = X(:, 1:4) - X(:, 5:8);
plot (T, ERR)
grid on
title ('Figure 7: e(t) = x(t) - xhat(t)')
xlabel ('Time [s]')
ylabel ('Observation Errors Full-Order Observer')
legend ('e1', 'e2', 'e3', 'e4', 'location', 'southeast')


% ===============================================================================

% Alternative Code
Naug = augstate (N);
S = sumblk ('e = x - xhat', 4);
N2 = connect (Naug, S, 'u', 'e')
[Y2, T2] = initial (N2, [x0; x0hat], 1);

figure (2)
plot (T2, Y2)
grid on
title ('Figure 7: e(t) = x(t) - xhat(t)')
xlabel ('Time [s]')
ylabel ('Observation Errors Full-Order Observer')
legend ('e1', 'e2', 'e3', 'e4', 'location', 'southeast')
