%******************************************************************************
% Project     : rtGl
% Autor       : Lukas von Arx
% Filename    : 4.6.5 Regler-Auslegung
% Startdate   : 30.04.2012
%******************************************************************************
% 
%
%******************************************************************************
%
% used functions:     -
%
% input files:        - 
% 
% output files:       - 
%
%******************************************************************************
%
% global variables:  -
%
%******************************************************************************

clc; clf;
clear all;
close all;

tau=10e-3;
td=tau/2;
kp=7;
kc=2/pi();

%1. System (Tiefpass):
num1   = [kp];
denum1 = [tau 1];
Gs = tf(num1, denum1);


%2. System (Delay):
num2   = [-(td*kc) 1];
denum2 = [td*kc 1];
Gr = tf(num2, denum2);

%Serieschaltung der Systeme (open Loop):
G0=Gr*Gs;

%Führungsfrequenzgang (closed Loop):
GF=Gr*Gs/(1+Gr*Gs);

% =======================================================================
% Plotten des Bodediagramms mit zusätzlichen Eigenschaften
% Bereiche
% Frequenzachsen-----------------------------------------
N = 10000;
%fmin = 1; fmax = 1e6;
wmin=0.1; 
wmax=1e3;
base = log10(wmin); limit = log10(wmax);
w = logspace(base,limit,N);

    %omega_min    = 1; 
    %omega_max    = 1e6;
    %omega_Range  = {omega_min, omega_max};     % !brackets!
    Magni_Range  = [1e-2,2e1]; 
    Phase_Range  = [-270, 45];
    Ampl_Range   = {Magni_Range; Phase_Range}; % !brackets!

%Bodediagramm:
%subplot(2,2,1);
%bode(Gr,'b', Gs,'g', G0,'r', w);
%PL=bodeplot(Gr,'y', Gs,'g', G0,'r');
%Legend('Gr', 'Gs','G0')
%grid on;

%Nyquist-Diagramm:
subplot(1,2,1); 
impulse(Gr);
grid on;
subplot(1,2,2); 
nyquist(Gr);
grid on;

% Schrittantwort
%subplot(2,2,1);                   % Plot in die untere Bildhälfte
%T = 1e-3;                       % Zeitausschnitt vorgeben
%step(GF);                 % Schrittantwort mit Farben
%step(GF,'r' ,T);                 % Schrittantwort mit Farben
%grid;                             % Gitternetz

% Pole Nullstellen Darstellung 
%subplot(2,2,2);                   % Plot in die obere Bildhälfte
%iopzplot(GF);                     % Pole / Nullstellen plotten 
%axis([-1e4; 1e4; -1e4; 1e4]);    % Bei Bedarf Achsen skalieren

% Setzen der Einheiten und Skalierung für die Amplitude; grid on 
%setoptions(PL,'MagUnits','abs','MagScale','log');   
%setoptions(PL,'YLimMode','manual','YLim', Ampl_Range, 'Grid','on')

%figure(2);
%PL=bodeplot(Gr,'y', Gs,'g',G0,'r', omega_Range);
%Legend('Gr', 'Gs', 'G0')
%grid on;

%bode(H, w);
%nyquist(H);
%step(H);
%impulse(H);
