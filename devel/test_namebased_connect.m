P = Boeing707;

I = ss (-eye (2));
I.inname = P.outname;
I.outname = P.inname;

T = connect (P, I, P.inname, P.outname)

T2 = feedback (P)  # should be the same model as T

