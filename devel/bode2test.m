A = ss (-2, 3, 4, 5);
B = ss (-1, 1, 1, 0);
C = ss (-4, 7, 10, 12);
w = [1 2 3 4 5 6 7];

%bode2 (A, B, "b", w, C)
%bode2 (A, B, "b", C)
%bode2 (A, B, "b", C, WestlandLynx)
%bode2 (WestlandLynx)
%bode2 (A, B, "b", c2d (C, 0.375), WestlandLynx)
bode2 (A, c2d (B, 0.223), "b", c2d (C, 0.375), WestlandLynx)

