% Heat Exchanger 9

% Tabula Rasa
clear all;
close all;
clc;

% Physical Parameters
m_dot_1 = 0.5;  % [kg/s]
m_dot_2 = 0.2;  % [kg/s]
V_1 = 0.01;     % [m^3]
V_2 = 0.02;     % [m^3]

rho = 1000;     % [kg/m^3]
c = 4200;       % [J/(kg K)]
Ar = 2;         % [m^2]
k = 2000;       % [W/(m^2 K)]

% Control-Oriented Parameters
tau_1 = rho * V_1 * c / (m_dot_1 * c  +  k * Ar);  % [s]
tau_2 = rho * V_2 * c / (m_dot_2 * c  +  k * Ar);  % [s]

sigma_1 = k * Ar / (m_dot_1 * c  +  k * Ar);       % [-]
sigma_2 = k * Ar / (m_dot_2 * c  +  k * Ar);       % [-]

beta_1 = m_dot_1 * c / (m_dot_1 * c  +  k * Ar);   % [-]
beta_2 = m_dot_2 * c / (m_dot_2 * c  +  k * Ar);   % [-]

% System Matrices
A = [      -1/tau_1    sigma_1/tau_1 ;
      sigma_2/tau_2         -1/tau_2 ];
      
B = [ beta_1/tau_1               0 ;
                 0    beta_2/tau_2 ];
             
C = [ 1    0 ;
      0    1 ];
      
D = [ 0    0 ;
      0    0 ];
      
% State-Space Model
P = ss (A, B, C, D)

% Stability
st = isstable (P)

% Controllability
co = isctrb (P)

% Observability
ob = isobsv (P)

% Invariant Zeros
z = zero (P)

% System Poles
p = pole (P)

% Singular Value Plot
figure (1)
sigma (P)

% Relative-Gain Array
RGA = tf (P) .* tf (inv (P)).'

RGA(0)  % FIXED by replacing minreal by sminreal in sys2tf

%{
Result from the Dark Side:

>> RGA = tf (P) .* tf (inv (P)).'

RGA =
 
  From input 1 to output...
       s^2 + 0.2029 s + 0.008368
   1:  -------------------------
       s^2 + 0.2029 s + 0.003833
 
               -0.004535
   2:  -------------------------
       s^2 + 0.2029 s + 0.003833
 
  From input 2 to output...
               -0.004535
   1:  -------------------------
       s^2 + 0.2029 s + 0.003833
 
       s^2 + 0.2029 s + 0.008368
   2:  -------------------------
       s^2 + 0.2029 s + 0.003833
 
Continuous-time transfer function.

>> freqresp (RGA, 0)

ans =

    2.1831   -1.1831
   -1.1831    2.1831

>> 
%}
