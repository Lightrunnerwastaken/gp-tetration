/* exp-031 prototype: Bluestein via pow2-FFT convolution (memory-lean).
   Verifies fft/fftinv semantics and equality with poly-mult bluedft. */
bluedft_pm(t) = {
  local(nn, ch, a, b, pp);
  nn = length(t);
  ch = exp(-Pi*I/nn);
  a = vector(nn, j, t[j] * ch^((j-1)^2));
  b = vector(2*nn-1, m, ch^(-(m-nn)^2));
  pp = Polrev(a) * Polrev(b);
  vector(nn, k, ch^((k-1)^2) * polcoeff(pp, (k-1)+(nn-1)));
}
bluedft_fft(t) = {
  local(nn, ch, a, b, M, w, fa, fb, pp);
  nn = length(t);
  ch = exp(-Pi*I/nn);
  a = vector(nn, j, t[j] * ch^((j-1)^2));
  b = vector(2*nn-1, m, ch^(-(m-nn)^2));
  M = 1; while (M < 3*nn-2, M = M*2);
  w = powers(exp(2*Pi*I/M), M-1);
  fa = fft(w, concat(a, vector(M-nn, i, 0)));
  fb = fft(w, concat(b, vector(M-(2*nn-1), i, 0)));
  pp = fftinv(w, vector(M, i, fa[i]*fb[i]));
  /* fftinv scaling check happens in the self-test below */
  vector(nn, k, ch^((k-1)^2) * pp[(k-1)+(nn-1)+1] / M);
}
{
  for (c = 1, 3,
    N = [100, 1000, 2304][c];
    t = vector(N, j, cos(j*j*0.37 + c) + I*sin(j*1.71 - c));
    gettime();
    g1 = bluedft_pm(t);  t1 = gettime();
    g2 = bluedft_fft(t); t2 = gettime();
    err = vecmax(vector(N, k, abs(g1[k]-g2[k]))) / vecmax(vector(N, k, abs(g1[k])));
    print("N=", N, " relerr=", precision(err,5), " polymult_ms=", t1, " fft_ms=", t2);
  );
}
