%% -*- texinfo -*-
%% Numerical optimization of a PID controller using an objective function.
%% The objective function is located in the file @command{optiPIDfun}.
%% Type @code{which optiPID} to locate, @code{edit optiPID} to open
%% and simply @code{optiPID} to run the example file.
%% In this example called @code{optiPID}, loosely based on [1], it is assumed
%% that the plant
%% @iftex
%% @tex
%% $$ P(s) = {1 \\over (s^{2} + s + 1)\\ (s + 1)^{4}} $$
%% @end tex
%% @end iftex
%% @ifnottex
%% @example
%%                   1
%% P(s) = -----------------------
%%        (s^2 + s + 1) (s + 1)^4 
%% @end example
%% @end ifnottex
%% is controlled by a PID controller with second-order roll-off
%% @iftex
%% @tex
%% $$ C(s) = k_p \\ (1 + {1 \\over T_i \\ s} + T_d \\ s) \\ {1 \\over (\\tau \\ s + 1)^{2}} $$
%% @end tex
%% @end iftex
%% @ifnottex
%% @example
%%                  1                1
%% C(s) = Kp (1 + ---- + Td s) -------------
%%                Ti s         (tau s + 1)^2
%% @end example
%% @end ifnottex
%% in the usual negative feedback structure
%% @iftex
%% @tex
%% $$ T(s) = {L(s) \\over 1 + L(s)} = {P(s) \\ C(s) \\over 1 + P(s) \\ C(s)} $$
%% @end tex
%% @end iftex
%% @ifnottex
%% @example
%%          L(s)       P(s) C(s)
%% T(s) = -------- = -------------
%%        1 + L(s)   1 + P(s) C(s)
%% @end example
%% @end ifnottex
%% The plant P(s) is of higher order but benign.  The initial values for the
%% controller parameters
%% @iftex
%% @tex
%% $k_p$, $T_i$ and $T_d$
%% @end tex
%% @end iftex
%% @ifnottex
%% Kp, Ti and Td
%% @end ifnottex
%% are obtained by applying the
%% Astroem and Haegglund rules [2].  These values are to be improved using a
%% numerical optimization as shown below.
%% As with all numerical methods, this approach can never guarantee that a
%% proposed solution is a global minimum.  Therefore, good initial guesses for
%% the parameters to be optimized are very important.
%% The Octave function @code{fminsearch} minimizes the objective function @var{J},
%% which is chosen to be
%% @iftex
%% @tex
%% $$ J(k_p, T_i, T_d) = \\mu_1 \\cdot \\int_0^{\\infty} \\! t \\ |e(t)| \\ dt \\ + \\ \\mu_2 \\cdot (|| y(t) ||_{\\infty} - 1) \\ + \\ \\mu_3 \\cdot ||S(jw)||_{\\infty} $$
%% @end tex
%% @end iftex
%% @ifnottex
%% @example
%%                     inf 
%% J(Kp, Ti, Td) = mu1 INT t |e(t)| dt  +  mu2 (||y(t)||    - 1)  +  mu3 ||S(jw)||
%%                      0                               inf                       inf
%% @end example
%% @end ifnottex
%% This particular objective function penalizes the integral of time-weighted absolute error
%% @iftex
%% @tex
%% $$ ITAE = \\int_0^{\\infty} \\! t \\ |e(t)| \\ dt $$
%% @end tex
%% @end iftex
%% @ifnottex
%% @example
%%        inf 
%% ITAE = INT t |e(t)| dt
%%         0             
%% @end example
%% @end ifnottex
%% and the maximum overshoot
%% @iftex
%% @tex
%% $$ y_{max} - 1 = || y(t) ||_{\\infty} - 1 $$
%% @end tex
%% @end iftex
%% @ifnottex
%% @example
%% y    - 1 = ||y(t)||    - 1
%%  max               inf
%% @end example
%% @end ifnottex
%% to a unity reference step
%% @iftex
%% @tex
%% $r(t) = \\varepsilon (t)$
%% @end tex
%% @end iftex
%% in the time domain. In the frequency domain, the sensitivity
%% @iftex
%% @tex
%% $$ M_s = ||S(jw)||_{\\infty} $$
%% @end tex
%% @end iftex
%% @ifnottex
%% @example
%% Ms = ||S(jw)||
%%               inf
%% @end example
%% @end ifnottex
%% is minimized for good robustness, where S(jw) denotes the @emph{sensitivity} transfer function
%% @iftex
%% @tex
%% $$ S(s) = {1 \\over 1 + L(s)} = {1 \\over 1 + P(s) \\ C(s)} $$
%% @end tex
%% @end iftex
%% @ifnottex
%% @example
%%            1            1
%% S(s) = -------- = -------------
%%        1 + L(s)   1 + P(s) C(s)
%% @end example
%% @end ifnottex
%% The constants
%% @iftex
%% @tex
%% $\\mu_1$, $\\mu_2$ and $\\mu_3$
%% @end tex
%% @end iftex
%% @ifnottex
%% mu1, mu2 and mu3
%% @end ifnottex
%% are @emph{relative weighting factors} or @guillemetleft{}tuning knobs@guillemetright{}
%% which reflect the importance of the different design goals. Varying these factors
%% corresponds to changing the emphasis from, say, high performance to good robustness.
%% The main advantage of this approach is the possibility to explore the tradeoffs of
%% the design problem in a systematic way.
%% In a first approach, all three design objectives are weigthed equally.
%% In subsequent iterations, the parameters
%% @iftex
%% @tex
%% $\\mu_1 = 1$, $\\mu_2 = 10$ and $\\mu_3 = 20$
%% @end tex
%% @end iftex
%% @ifnottex
%% mu1 = 1, mu2 = 10 and mu3 = 20
%% @end ifnottex
%% are found to yield satisfactory closed-loop performance.  This controller results
%% in a system with virtually no overshoot and a phase margin of 64 degrees.
%%
%% @*@strong{References}@*
%% [1] Guzzella, L.
%% @cite{Analysis and Design of SISO Control Systems},
%% VDF Hochschulverlag, ETH Zurich, 2007@*
%% [2] Astroem, K. and Haegglund, T.
%% @cite{PID Controllers: Theory, Design and Tuning},
%% Second Edition,
%% Instrument Society of America, 1995

% ===============================================================================
% optiPID                          Lukas Reichlin                       July 2009
% ===============================================================================
% Numerical Optimization of an A/H PID Controller
% Required OCTAVE Packages: control
% Required MATLAB Toolboxes: Control, Optimization
% ===============================================================================

% Tabula Rasa
clear all, close all, clc;

% Global Variables
global P t dt mu_1 mu_2 mu_3

% Plant
numP = [1];
denP = conv ([1, 1, 1], [1, 4, 6, 4, 1]);
P = tf (numP, denP);

% Relative Weighting Factors: PLAY AROUND WITH THESE!
mu_1 = 1;               % Minimize ITAE Criterion
mu_2 = 10;              % Minimize Max Overshoot
mu_3 = 20;              % Minimize Sensitivity Ms

% Simulation Settings: PLANT-DEPENDENT!
t_sim = 30;             % Simulation Time [s]
dt = 0.05;              % Sampling Time [s]
t = 0 : dt : t_sim;     % Time Vector [s]

% A/H PID Controller: Ms = 2.0
[gamma, phi, w_gamma, w_phi] = margin (P);

ku = gamma;
Tu = 2*pi / w_gamma;
kappa = inv (dcgain (P));

disp ('optiPID: Astrom/Hagglund PID controller parameters:');
kp_AH = ku * 0.72 * exp ( -1.60 * kappa  +  1.20 * kappa^2 )
Ti_AH = Tu * 0.59 * exp ( -1.30 * kappa  +  0.38 * kappa^2 )
Td_AH = Tu * 0.15 * exp ( -1.40 * kappa  +  0.56 * kappa^2 )

C_AH = optiPIDctrl (kp_AH, Ti_AH, Td_AH);

% Initial Values
C_par_0 = [kp_AH; Ti_AH; Td_AH];

% Optimization
warning ('optiPID: optimization starts, please be patient ...');
C_par_opt = fminsearch (@optiPIDfun, C_par_0);

% Optimized Controller
disp ('optiPID: optimized PID controller parameters:');
kp_opt = C_par_opt(1)
Ti_opt = C_par_opt(2)
Td_opt = C_par_opt(3)

C_opt = optiPIDctrl (kp_opt, Ti_opt, Td_opt);

% Open Loop
L_AH = P * C_AH;
L_opt = P * C_opt;

% Closed Loop
T_AH = feedback (L_AH, 1);
T_opt = feedback (L_opt, 1);

% A Posteriori Stability Check
disp ('optiPID: closed-loop stability check:');
st_AH = isstable (T_AH)
st_opt = isstable (T_opt)

% Stability Margins
disp ('optiPID: gain margin gamma [-] and phase margin phi [deg]:');
[gamma_AH, phi_AH] = margin (L_AH)
[gamma_opt, phi_opt] = margin (L_opt)

% Plot Step Response
figure (1)
step (T_AH, 'b', T_opt, 'r', t)
legend ('Astroem/Haegglund PID', 'Optimized PID', 'Location', 'SouthEast')

% ===============================================================================
