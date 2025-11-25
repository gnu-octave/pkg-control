% bug #46330: segfault in arrayfun
% <https://savannah.gnu.org/bugs/index.php?46330>

% Test Data
a = [-0.46E-01,            0.10681415316, 0.0,   -0.17121680433;
     -0.1675901504661613, -0.515,         1.0,    0.6420630320636088E-02;
      0.1543104215347786, -0.547945,     -0.906, -0.1521689385990753E-02;
      0.0,                 0.0,           1.0,    0.0];

b = [0.1602300107479095,      0.2111848453E-02;
     0.8196877780963616E-02, -0.3025E-01;
     0.9173594317692437E-01, -0.75283075;
     0.0,                     0.0];

c = [1.0, 0.0, 0.0, 0.0;
     0.0, 0.0, 0.0, 1.0];

d = zeros (2, 2);

e = eye (4);

w = Inf;

% Test Function
s = i * w;
H = arrayfun (@(x) c/(x*e - a)*b + d, s, "uniformoutput", false);

%{
octave:1> bug_arrayfun 
warning: matrix singular to machine precision, rcond = nan
 ** On entry to DLASCL parameter number  4 had an illegal value
panic: Segmentation fault: 11 -- stopping myself...
attempting to save variables to 'octave-workspace'...
save to 'octave-workspace' complete
Segmentation fault: 11
%}
