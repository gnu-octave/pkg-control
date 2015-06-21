## http://savannah.gnu.org/bugs/?45336

a1 = tf ({[2 24], [12]; [1 48 0], [2 24]}, {[-1 24], [-1 24]; [-4 96], [-1 24]})

b1 = minreal (a1)

%{
## http://savannah.gnu.org/bugs/?45334

a1 = tf ({[1], [2]; [3], [4]}, {[1 3], [1 4]; [2 3], [1 1]}).'

tf (ss (a1))
minreal (tf (ss (a1)))

b1 = a1
%b1 = minreal (a1)
%}

%{
b2 = b1^2

b3 = b2 * b1

b4 = b2^2
b4 = b3 * b1

% b5 = b4 * b1
b5 = b3 * b2

% b6 = b4 * b2
b6 = b3^2

b7 = b6 * b1
b7 = b5 * b2
b7 = b4 * b3
%}

b2 = b1^2
b2m = minreal (b2)

b3 = b1^3
b3m = minreal (b3)

b4 = b1^4
b4m = minreal (b4)

b5 = b1^5
b5m = minreal (b5)

b6 = b1^6
b6m = minreal (b6)

b7 = b1^7
b7m = minreal (b7)
