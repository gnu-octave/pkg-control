index = [ 3     3 ];
   
dcoeff = [ 1.0   6.0  11.0   6.0
           1.0   6.0  11.0   6.0 ];
   
ucoeff = zeros (2, 2, 4);

u11 = [ 1.0   6.0  12.0   7.0 ];
u12 = [ 0.0   1.0   4.0   3.0 ];
u21 = [ 0.0   0.0   1.0   1.0 ];
u22 = [ 1.0   8.0  20.0  15.0 ];

ucoeff(1,1,:) = u11;
ucoeff(1,2,:) = u12;
ucoeff(2,1,:) = u21;
ucoeff(2,2,:) = u22;
   
[a, b, c, d] = __sl_td04ad__ (ucoeff, dcoeff, index, 0)

