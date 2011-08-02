% ===============================================================================
% Robust Control of a Mass-Damper-Spring System
% ===============================================================================
% Reference: Gu, D.W., Petkov, P.Hr. and Konstantinov, M.M.
%            Robust Control Design with Matlab, Springer 2005
% ===============================================================================

% Tabula Rasa
clear all, close all, clc

% ===============================================================================
% System Model
% ===============================================================================
%                +---------------+  
%                | d_m   0    0  |
%          +-----|  0   d_c   0  |<----+
%      u_m |     |  0    0   d_k |     | y_m
%      u_c |     +---------------+     | y_c
%      u_k |                           | y_k
%          |     +---------------+     |
%          +---->|               |-----+
%                |     G_nom     |
%        u ----->|               |-----> y
%                +---------------+

% Nominal Values
m_nom = 3;   % Mass
c_nom = 1;   % Damping Coefficient
k_nom = 2;   % Spring Stiffness

% Perturbations
p_m = 0.4;   % 40% uncertainty in the mass
p_c = 0.2;   % 20% uncertainty in the damping coefficient
p_k = 0.3;   % 30% uncertainty in the spring stiffness

% State-Space Representation
A =   [            0,            1
        -k_nom/m_nom, -c_nom/m_nom ];

B1 =  [            0,            0,            0
                -p_m,   -p_c/m_nom,   -p_k/m_nom ];

B2 =  [            0
             1/m_nom ];

C1 =  [ -k_nom/m_nom, -c_nom/m_nom
                   0,        c_nom
               k_nom,            0 ];

C2 =  [            1,            0 ];

D11 = [         -p_m,   -p_c/m_nom,   -p_k/m_nom
                   0,            0,            0
                   0,            0,            0 ];

D12 = [      1/m_nom
                   0
                   0 ];

D21 = [            0,            0,            0 ];

D22 = [            0 ];

inname = {'u_m', 'u_c', 'u_k', 'u'};    % Input Names
outname = {'y_m', 'y_c', 'y_k', 'y'};   % Output Names

G_nom = ss (A, [B1, B2], [C1; C2], [D11, D12; D21, D22], \
            'inputname', inname, 'outputname', outname);


% ===============================================================================
% Frequency Analysis of Uncertain System
% ===============================================================================

% Uncertainties: -1 <= delta_m, delta_c, delta_k <= 1
[delta_m, delta_c, delta_k] = ndgrid ([-1, 0, 1], [-1, 0, 1], [-1, 0, 1]);
w = logspace (-1, 1, 100);  % Frequency Vector
figure (1)

for k = 1 : numel (delta_m)
  Delta = diag ([delta_m(k), delta_c(k), delta_k(k)]);
  G_per = lft (Delta, G_nom);
  bode (G_per, w)
  subplot (2, 1, 1)
  hold on
  subplot (2, 1, 2)
  hold on
endfor


% ===============================================================================
% Mixed Sensitivity H-infinity Controller Design
% ===============================================================================
%                                    +-------+
%             +--------------------->|  W_p  |----------> e_p
%             |                      +-------+
%             |                      +-------+
%             |                +---->|  W_u  |----------> e_u
%             |                |     +-------+
%             |                |    +---------+
%             |                |  ->|         |->
%  r   +    e |   +-------+  u |    |  G_nom  |
% ----->(+)---+-->|   K   |----+--->|         |----+----> y
%        ^ -      +-------+         +---------+    |
%        |                                         |
%        +-----------------------------------------+

% Weighting Functions
s = tf ('s');
W_p = 0.95 * (s^2 + 1.8*s + 10) / (s^2 + 8.0*s + 0.01);
W_u = 10^-2;

% Mixed Sensitivity H-infinity Controller Design
G = G_nom(4, 4);    % Extract output y and input u
%K = mixsyn (G, W_p, W_u);
%[K, ~, gamma] = mixsyn (G, W_p, W_u)
[K, N, gamma] = mixsyn (G, W_p, W_u)
% g_max = 10;
% K = mixsyn (G, W_p, W_u, [], g_max)
% [K, ~, gamma] = mixsyn (G, W_p, W_u, [], g_max)

% Open Loop
L = G * K;

% Closed Loop
T = feedback (L);

figure (2)
sigma (T)

figure (3)
sigma (N)

norm (N, inf)

figure (4)
step (T)


% ===============================================================================
% H-infinity Loop-Shaping Design
% ===============================================================================

W1 = 8 * (2*s + 1) / (0.9*s);     % Precompensator
W2 = 1;                           % Postcompensator
factor = 1.1                      % Suboptimal Controller

% Compute the suboptimal positive feedback controller
K = ncfsyn (G, W1, W2, factor)

K = -K;  % negative feedback controller

% Open Loop
L = G * K;

% Closed Loop
T = feedback (L);

figure (5)
step (T)

% ===============================================================================
