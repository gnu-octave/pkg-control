function [h,g,p,ks,ch,cg] = step_function (G)

  [z,p,k] = zpkdata (G);
  z = z{1};
  p = p{1};
  m = length (z);
  n = length (p);

  [b, a] = tfdata (G, 'vector');
  ks = b(n+1)/a(n+1);

  %  H(t=0), for  s->inf:
  %
  %  h(0)    = s*H(s) = G(s)
  %  h'(0)   = s*(G(s) - h(0));
  %  h''(0)  = s*(s*G(s) - s*h(0) - h'(0))
  %  h'''(0) = s*(s*s*G(s) - s*s*h(0) - s*h'(0) - h''(0))
  %  ...

  h0 = zeros (n,1);
  bn = b;

  for i = 1:n

    h0(i) = bn(1)/a(1);

    bn = bn - h0(i)*a;
    bn = [bn(2:end) 0];

  end

  %
  % h(t) = KS + c1*e^{p1*t} + ... cn*e^{pn*t})
  %
  %  0:  h0 - KS = c1 + c2 + ... + cn
  %  1:       h1 = c1*p1 + c2*p2 + ... + cn*pn
  %                ...
  %  n-1:     hn = c1*p1^(n-1) + c2*p2^(n-1) + ... + cn*pn^(n-1)
  %

  y = h0;
  y(1) = y(1) - ks;

  X = ones (n,n);

  for i = 2:n
    X(i,:) = X(i-1,:) .* p';
  end

  ch = inv(X)*y;
  cg = ch .* p;

  h = @(t) real (ks + ch'*exp(p*t));
  g = @(t) real (cg'*exp(p*t));

end
