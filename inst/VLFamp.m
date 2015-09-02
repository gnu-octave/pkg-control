## Copyright (C) 2015   Thomas D. Dean
##
## This file is part of LTI Syncope.
##
## LTI Syncope is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## LTI Syncope is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with LTI Syncope.  If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn{Function File} {} VLFamp
## @deftypefnx{Function File} {@var{result} =} VLFamp (@var{verbose})
## Calculations on a two stage preamp for a multi-turn,
## air-core solenoid loop antenna for the reception of
## signals below 30kHz.
##
## The Octave Control Package functions are used extensively to
## approximate the behavior of operational amplifiers and passive
## electrical circuit elements.
##
## This example presents several 'screen' pages of documentation of the
## calculations and some reasoning about why.  Plots of the results are
## presented in most cases.
##
## The process is to display a 'screen' page of text followed by the
## calculation and a 'Press return to continue' message.  To proceed in
## the example, press return.  ^C to exit.
##
## At one point in the calculations, the process may seem to hang, but,
## this is because of extensive calculations.
##
## The returned transfer function is more than 100 characters long so
## will wrap in screens that are narrow and appear jumbled.
## @end deftypefn

## Author: Thomas D. Dean <tomdean@wavecable.com>
## Created: June 2015
## Version: 0.1

function retval = VLFamp (verbose = false)

  if (nargin > 1)
    print_usage ();
  endif

  clc;
  disp ("---- VLF Pre-Amplifier Design ----");
  disp ("");
  disp ("This example covers the design of a pre-amplifier for use in");
  disp ("receiving radio frequencies below 30kHz.");
  disp ("");
  disp ("See http://www.vlf.it for details of Narural Radio Sources");
  disp ("");
  disp ("The Octave Control Package functions are used extensively to");
  disp ("approximate the behavior of operational amplifiers and passive");
  disp ("electrical circuit elements.");
  disp ("");
  disp ("This example presents several 'screen' pages of documentation of the");
  disp ("calculations and some reasoning about why.  Plots of the results are");
  disp ("presented in most cases.");
  disp ("");
  disp ("Often, when multiple plots are displayed, they may be overlaid");
  disp ("on the screen.  You may use the mouse and move them for better viewing.");
  disp ("");
  disp ("The process is to display a 'screen' page of text followed by the");
  disp ("calculation and a 'Press return to continue' message. To proceed in");
  disp ("the example, press return.  ^C to exit.");
  disp ("");
  disp ("At one point in the calculations, the process may seem to hang, but,");
  disp ("this is because of extensive calculations.");
  disp ("");
  disp ("The returned transfer function is more than 100 characters long so");
  disp ("will wrap in screens that are narrow and appear jumbled.");
  disp ("");
  ##
  input ("Press Return to continue:");
  blanks ();
  ##
  disp ("");
  disp ("The amplifier consists of two AD797 op amps and a low pass filter.");
  disp ("With biasing and blocking capacitors omitted, three blocks remain.");
  disp ("");
  disp ("");
  disp ("          Gain = 10");
  disp ("        +-------------+");
  disp ("        |             |      -- Low Pass Filter --");
  disp ("     ---+ p           |");
  disp (" Loop   |   Stage 1   +--+----R3--+--R4--+--R5--+---> To Stage 2");
  disp ("        |   Amplifier |  |        |      |      |");
  disp ("     -+-+ n           |  |        C1     C2     C3");
  disp ("      | |             |  |        |      |      |");
  disp ("      | +-------------+  |       Gnd    Gnd    Gnd");
  disp ("      |                  |");
  disp ("      +----+---R2--------+");
  disp ("           |");
  disp ("           R1");
  disp ("           |");
  disp ("          Gnd");
  disp ("");
  disp ("");
  disp ("");
  disp ("                   Gain = 10");
  disp ("                 +-------------+");
  disp ("                 |             |");
  disp ("            Gnd--+ p           |");
  disp ("                 |   Stage 2   +--+----R8--+----> Output");
  disp ("                 |   Amplifier |  |        |");
  disp ("      From >---+-+ n           |  |        R9");
  disp ("     Filter    | |             |  |        |");
  disp ("               | +-------------+  |       Gnd");
  disp ("               |                  |");
  disp ("               +----+---R6--------+");
  disp ("                    |");
  disp ("                    R7");
  disp ("                    |");
  disp ("                   Gnd");
  disp ("");
  disp ("");
  disp ("R1 and R2 profide feedback to control the gain of Stage 1.");
  disp ("R3 through R5 with C1 through C3 form a low pass filter to limit the");
  disp ("   bandwidth.");
  disp ("R6 and R7 profide feedback to control the gain of Stage 2.");
  disp ("R8 and R9 provide impedance matching to the cable and/or receiver,");
  disp ("   possibly a PC sound card.");
  disp ("");
  ##
  input ("Press Return to continue:");
  blanks ();
  ##
  disp ("");
  disp ("");
  disp ("The graphs in the ad797 datasheet reveal the following parameters:");
  disp ("");

  show ("a0 = 1e7;   ## Open Loop Gain");
  show ("p1 = 55;    ## Pole (Hz)");
  show ("p2 = 1e6;   ## Pole (Hz)");
  show ("z1 = 4.3e6; ## Zero (Hz)");

  disp ("");
  disp ("The open loop transfer function of an op amp with m zeros and n");
  disp ("poles is expressed in the form:");
  disp ("  tf = open_loop_gain * zero_expressions / pole_expressions");
  disp ("where ");
  disp ("  zero_expressions = (1+s/z1) * (1+s/z2) * ... * (1+s/zm) ");
  disp ("  pole_expressions = (1+s/p1) * (1+s/p2) * ... * (1+s/pn)");
  disp ("  z1 ... zm are the m zeros");
  disp ("  p1 ... pn are the n poles");
  disp ("");

  ##
  input ("Press Return to continue:");
  blanks ();
  ##

  disp ("");
  disp ("The amplifier stages have 1 zero and 2 poles:");
  disp ("");
  show ("s = tf ('s')")
  disp ("");
  show ("TFopen = a0 * (1+s/z1) / (1+s/p1) / (1+s/p2)")
  disp ("");
  show ("TFopen_norm = minreal (TFopen)")

  disp ("");
  disp ("Note: The difference between the op amp expression and the usual");
  disp ("Zero-Pole-Gain expression is in the modification of the gain");
  disp ("parameter.  The gain argument to zpk() is modified by the zeros");
  disp ("and poles, so the derived transfer function matches actual");
  disp ("measurements.");
  disp ("");

  show ("Azpk = zpk ([-z1], [-p1, -p2], 1e7*p1*p2/z1)")

  ##
  input ("Press Return to continue:");
  blanks ();
  ##
  disp ("");
  disp ("The bode plot of these two open loop transfer functions produce");
  disp ("identical results.  And, the plots show the same shape as the");
  disp ("graphs in the datasheet.");
  disp ("");

  show ("figure 1");
  show ("bode (TFopen)");
  show ("subplot (2,1,1)");
  show ("title ('Equation Bode Diagram')");
  show ("figure 2");
  show ("bode (Azpk)");
  show ("subplot (2,1,1)");
  show ("title ('ZPK Bode Diagram')");

  disp ("");
  disp ("Two Bode Diagrams should be visible, possibly overlaid.");
  disp ("");
  ##
  input ("Press Return to close the plots and continue:");
  blanks ();
  ##

  close all;
  disp ("");
  disp ("The normalized step response of the ad797 is:");

  disp ("");
  show ("TFnorm = TFopen/dcgain(TFopen)")
  disp ("");
  show ("step (TFnorm, 'b')");
  show ("title ('AD797 Normalized Open-Loop Step Response')");
  show ("ylabel ('Normalized Amplitude')");
  disp ("");

  ##
  input ("Press Return to close the plot and continue:");
  blanks ();
  ##

  close all;
  disp ("");
  disp ("--- Design Stage 1 of the VLFamp ---");
  disp ("");
  disp ("Resistors R1 and R2 form a feedback system to control the gain of ");
  disp ("Stage 1.  This feedback system returns a portion of the output to the");
  disp ("negative input.  This is normally expressed as:");
  disp ("  Vfb = Vout * R1 / (R1 + R2)");
  disp ("So, the transfer function of the feedback network is:");
  disp ("  tf = Vfb / Vout = R1 / (R1 + R2)");
  disp ("The effects of the AD797 gain on the input and the feedback may be ");
  disp ("represented as TFstage1 = Vout/Vp = gain / (1 + dcgain * TFfeedback).");
  disp ("If dcgain is sufficiently large, this reduces to");
  disp ("   TFstage1 = 1 / TFfeedback.");
  disp ("The dcgain of the AD797 is >> 1, so, the feedback completely controls");
  disp ("the output and variations in the dcgain will not effect the Stage gain.");
  disp ("");
  disp ("The feedback is added to the AD797 using the feedback function");
  disp ("");
  show ("Gfb = 10");
  show ("b  = 1 / Gfb");
  show ("R1 = 10e3");
  show ("R2 = R1 * (1/b - 1)")
  disp ("");
  show ("TFstage1 = feedback (TFopen, b)");
  disp ("");

  show ("bodemag (TFopen, 'r', TFstage1, 'b')");
  show ("legend ('Open Loop Gain (TFopen)', 'Closed Loop Gain (TFstage1)')");
  disp ("");

  disp ("The use of negative feedback to reduce the low-frequency (LF) gain");
  disp ("has led to a corresponding increase in the system bandwidth (defined");
  disp ("as the frequency where the gain drops 3dB below its maximum value).");
  disp ("");
  disp ("With this feedback, we have a gain of 10, or 20db up to 10MHz,");
  disp ("far more than the frequency range of interest.");
  disp ("");
  ##
  input ("Press Return to close the plot and continue:");
  blanks ();
  ##
  close all;
  disp ("");
  disp ("Since the gain is now dominated by the feedback network, a useful");
  disp ("relationship to consider is the sensitivity of this gain to variation");
  disp ("in the op amp's open-loop gain.");
  disp ("");
  disp ("Before deriving the system sensitivity, however, it is useful to");
  disp ("define the loop gain, L(s)=a(s)b(s), which is the total gain a signal");
  disp ("experiences traveling around the loop:");
  disp ("");
  disp ("Sensitivity = partial(TFstage1/TFopen)*TFopen/TFstage1");
  disp ("or S(s) = 1 / (1 + TFopen(s) * TFstage1(s))");
  disp ("or S(s) = 1 / (1 + L(s)), which has the same form as feedback");
  disp ("So, use the feedback function to develop the sensitivity.");
  disp ("");
  
  show ("L = TFopen * b")
  disp ("");
  show ("Sens = feedback (1, L)")
  disp ("");
  show ("figure 1");
  show ("bodemag (TFstage1, 'b', Sens, 'g')");
  disp ("");

  disp ("The very small low-frequency sensitivity (more than -100 dB) indicates");
  disp ("a design whose closed-loop gain suffers minimally from open-loop gain");
  disp ("variation. Such variation in a(s) is common due to manufacturing");
  disp ("variability, temperature change, etc.");
  disp ("");
  ##
  input ("Press Return to close the plot and continue:");
  blanks ();
  ##
  disp ("");
  disp ("You can check the step response of A(s) using the STEP command:");
  disp ("");

  show ("figure 2");
  show ("step (TFstage1)");

  disp ("");
  disp ("The stability margin can be analyzed by plotting the loop gain, L(s)");
  disp ("with the margin function.");
  disp (" ");
  disp ("This plot may display  warning messages, you can safely ignore them.");
  disp (" "); fflush(stdout);
  show ("margin (L)");
  disp (" "); fflush(stdout); fflush(stderr);
  disp (" ");
  disp ("Two plots are displayed, possibly overlaid.");
  disp (" ");
  ##
  input ("Press Return to close the plots and continue:");
  blanks ();
  ##

  disp ("");
  disp ("The plot indicates a phase margin of less than 3 degrees.  Stage 1");
  disp ("needs to be compensated to increase this to an acceptible level,");
  disp ("more than 45 degrees, if possible.");
  disp ("");
  disp ("Feedback Lead Compensation");
  disp ("");
  disp ("A commonly used method of compensation in this type of circuit is");
  disp ("feedback lead compensation. This technique modifies b(s) by adding");
  disp ("a capacitor, C, in parallel with the feedback resistor, R2.");
  disp ("The capacitor value is chosen so as to introduce a phase lead to b(s)");
  disp ("near the crossover frequency, thus increasing the amplifier's phase");
  disp ("margin.");
  disp ("The new feedback transfer function is shown below.");
  disp ("You can approximate a value for C by placing the zero of b(s) at the");
  disp ("0dB crossover frequency of L(s):");
  disp ("");

  show ("[Gm, Pm, Wcg, Wcp] = margin (L)");
  show ("C = 1/(R2*Wcp)")
  disp ("");
  if (C < 1e-12)
    disp ("The calculated value of C is very small.");
    disp ("Now, look at a range of values.");
  endif;

  disp (" ");
  disp ("The next plots take some time...");
  disp (" ");

  ##
  input ("Press Return to close the plot and continue:");
  blanks ();
  ##
  close all;
  disp ("The next plots take some time...");
  disp ("");

  show ("K = R1/(R1+R2);");
  show ("C = [10:10:200]*1e-12;");
  show ("b_array = arrayfun (@(C) tf ([K*R2*C, K], [K*R2*C, 1]), C,'uniformoutput',false);");
  show ("A_array = cellfun (@feedback, {TFopen}, b_array, 'uniformoutput', false);");
  show ("L_array = cellfun (@mtimes, {TFopen}, b_array, 'uniformoutput', false);");
  show ("S_array = cellfun (@feedback, {1}, L_array, 'uniformoutput', false);");
  disp (" "); fflush(stdout);
  show ("[Gm, Pm, Wcg, Wcp] = cellfun (@margin, L_array);");
  disp (" ");

  close all
  show ("figure 1");
  show ("step (TFstage1, 'r', A_array{:})");
  show ("figure 2");
  show ("bode (TFstage1, A_array{:})");
  show ("figure 3");
  show ("plot (C, Pm)");
  show ("grid");
  show ("xlabel ('Compensation Capacitor, C (pF)')");
  show ("ylabel ('Phase Margin (deg)')");
  show ("figure 4");
  show ("step (A_array{C==50e-12}, 'r', A_array{C==100e-12}, 'b', A_array{C==200e-12}, 'g')");
  show ("legend ('Compensated (50 pF)', 'Compensated (100 pF)', 'Compensated (200 pF)')");

  disp (" ");
  disp ("Four plots are displayed, possibly overlaid.");
  disp (" ");
  ##
  input ("Press Return to close the plots and continue:");
  blanks ();
  ##
  close all;
  disp ("");
  disp ("");
  disp ("          Gain = 10");
  disp ("        +-------------+");
  disp ("        |             |");
  disp ("     ---+ p           |");
  disp (" Loop   |   Stage 1   +--+---->");
  disp ("        |   Amplifier |  |");
  disp ("     -+-+ n           |  |");
  disp ("      | |             |  |");
  disp ("      | +-------------+  |");
  disp ("      |                  |");
  disp ("      +----+------R2-----+");
  disp ("           |             |");
  disp ("           +-----Ccomp---+");
  disp ("           |");
  disp ("           |");
  disp ("           R1");
  disp ("           |");
  disp ("          Gnd");
  disp ("");
  disp ("The selected compensation capacitor is 100pf.");

  show ("TFcomp = A_array{C==100e-12}");
  show ("bode (TFopen, 'b', TFstage1, 'g', TFcomp, 'r')");
  show ("legend ('TFopen', 'TFstage1', 'TFcomp')");

  disp ("");
  ##
  input ("Press Return to close the plot and continue:");
  blanks ();
  ##
  close all;
  disp ("");
  disp ("--- Low Pass Filter Design ---");
  disp ("");
  disp ("The low pass filter is composed of three equal sections.");
  disp ("Develop one section and put three in series.");
  disp ("");
  
  show ("C = 20e-9");
  show ("R = 1000");
  show ("TFsection = tf ([1], [C*R, 1])");
  disp ("");
  show ("TFfilter = TFsection * TFsection * TFsection;");
  if (verbose)
    TFfilter
  endif;

  disp ("");
  disp ("---- Final Design ----");
  disp ("");
  disp ("The final configuration is:  AD797 --> LP Filter --> AD797");
  disp ("");

  show ("TFpreamp = TFcomp * TFfilter * TFcomp;");

  show ("figure 1");
  show ("bode (TFpreamp, {1, 1e5})");
  show ("figure 2");
  show ("margin (TFpreamp)");

  disp ("");
  disp ("Two plots are displayed, possibly overlaid.");
  disp ("");
  ##
  input ("Press Return to close the plots and continue:");
  blanks ();
  ##
  disp ("");
  disp ("As can be seen from the plots, the gain margin is almost 30db.");
  disp ("The phase margin is 230 degrees.");
  disp ("");
  ## disp ("Use 'close all' to close the plots.");
  ##
  close all
  blanks ();
  disp ("The resultant transfer function is over 100 characters long");
  disp ("and will appear jumbled on narrower screens.");
  disp ("");
  show ("TFpreamp")
  ##
  
  if (nargout > 0)
    retval = TFpreamp;
  endif

endfunction


## support function to display a command and then
## execute it in the caller's environment.
function show (str)

  disp ([">> ", str]);
  evalin ("caller", str);

endfunction


## support function to insert blank lines in the display
function blanks (n = 5)

  for idx = 1:n
    disp ("");
  endfor

endfunction
