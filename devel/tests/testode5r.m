%{
pkg load odepkg

function vdy = frob (t, y, varargin)
  vdy(1,1) = -0.04*y(1)+1e4*y(2)*y(3);
  vdy(2,1) =  0.04*y(1)-1e4*y(2)*y(3)-3e7*y(2)^2;
  vdy(3,1) =  y(1)+y(2)+y(3)-1;
endfunction

function vmas = fmas (vt, vy, varargin)
  vmas = [1,0,0; 0,1,0; 0,0,0];
endfunction

A = odeset ("Mass", @fmas);
B = ode5r (@frob, [0 1e8], [1 0 0], A);
%}

% Tabula Rasa
clear all, close all, clc

pkg load odepkg


function vdy = frob (t, y, A, B, C, D)

 vdy = [A, zeros(rows(A),rows(C)); C, -eye(rows(C))] * y + [B; D];

endfunction


M = blkdiag (eye (6), 0);



% ===============================================================================
% Robust Control of a Mass-Damper-Spring System     Lukas Reichlin    August 2011
% ===============================================================================
% Reference: Gu, D.W., Petkov, P.Hr. and Konstantinov, M.M.
%            Robust Control Design with Matlab, Springer 2005
% ===============================================================================


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
m_nom = 3;   % mass
c_nom = 1;   % damping coefficient
k_nom = 2;   % spring stiffness

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

inname = {'u_m', 'u_c', 'u_k', 'u'};   % input names
outname = {'y_m', 'y_c', 'y_k', 'y'};  % output names

G_nom = ss (A, [B1, B2], [C1; C2], [D11, D12; D21, D22], ...
            'inputname', inname, 'outputname', outname);

G = G_nom(4, 4);                       % extract output y and input u


% ===============================================================================
% Mixed Sensitivity H-infinity Controller Design (S over KS Method)
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
s = tf ('s');                          % transfer function variable
W_p = 0.95 * (s^2 + 1.8*s + 10) / (s^2 + 8.0*s + 0.01);  % performance weighting
W_u = 10^-2;                           % control weighting

% Synthesis
K_mix = mixsyn (G, W_p, W_u);          % mixed-sensitivity H-infinity synthesis

% Interconnections
L_mix = G * K_mix;                     % open loop
T_mix = feedback (L_mix);              % closed loop


% ===============================================================================
% H-infinity Loop-Shaping Design (Normalized Coprime Factor Perturbations)
% ===============================================================================

% Settings
W1 = 8 * (2*s + 1) / (0.9*s);          % precompensator
W2 = 1;                                % postcompensator
factor = 1.1;                          % suboptimal controller

% Synthesis
K_ncf = ncfsyn (G, W1, W2, factor);    % positive feedback controller

% Interconnections
K_ncf = -K_ncf;                        % negative feedback controller
L_ncf = G * K_ncf;                     % open loop
T_ncf = feedback (L_ncf);              % closed loop


% ===============================================================================
% Plot Results
% ===============================================================================

% Step Response
figure (1)
step (T_mix, T_ncf, 10)                % step response for 10 seconds

% ===============================================================================



[A, B, C, D, E] = ssdata (T_mix);

opt = odeset ("Mass", M);
sol = ode5r (@frob, [0, 10], zeros(1, 7), opt, A, B, C, D);

figure (2)
t = sol.x;
y = sol.y(:,7);
plot (t, y)
grid on
