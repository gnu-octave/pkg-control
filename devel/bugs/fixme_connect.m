s = tf ('s');

F = 1/(20*s+1);
F.InputName = 'dy'; F.OutputName = 'dp';

% Process
%P = exp(-93.9*s) * 5.6/(40.2*s+1);
P = 1/(-93.9*s) * 5.6/(40.2*s+1);
P.InputName = 'u'; P.OutputName = 'y0';

% Prediction model
Gp = 5.6/(40.2*s+1);
Gp.InputName = 'u'; Gp.OutputName = 'yp';

%Dp = exp(-93.9*s);
Dp = 1/(-93.9*s);
Dp.InputName = 'yp'; Dp.OutputName = 'y1';

% Overall plant
S1 = sumblk ('ym = yp + dp');
S2 = sumblk ('dy = y0 - y1');


Plant = connect (P, Gp, Dp, F, S1, S2, 'u', 'ym');
% error: lti: inname 'u' is ambiguous
% FIXME: the dark side uses some kind of "vertcat" for common input names


% PI controller
Kp = 0.574;
Ti = 40.1;
C = Kp*(1 + 1/(Ti*s))
C.InputName = 'e'; C.OutputName = 'u';


% Assemble closed-loop model from [y_sp,d] to y
Sum1 = sumblk ('e = ysp - yp - dp');
Sum2 = sumblk ('y = y0 + d');
Sum3 = sumblk ('dy = y - y1');


T = connect (P, Gp, Dp, C, F, Sum1, Sum2, Sum3, {'ysp', 'd'}, 'y')
