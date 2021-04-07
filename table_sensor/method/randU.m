function U = randU(n) 
X = randn(n);
  [Q,R] = qr(X);
  r = sign(diag(R));
  U = bsxfun(@times, Q, r');