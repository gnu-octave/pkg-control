 A=tf(20,[1 2 30]);
 bode(A)
 dcgain1=dcgain(A)
 [ mag,p,w]=bode(A);
 q=find(mag<dcgain1/2);
 index1=q(1);
 gain2=mag(index1)
 bw=w(index1)
 st1=mat2str(bw,6)
 disp(["the bandwidth of this lowpass filter is " st1 " radians per second"])