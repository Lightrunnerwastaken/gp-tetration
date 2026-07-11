/* exp-020 prototype: exact-N DFT via Bluestein chirp + polynomial
   multiplication. Goal: free the efam staylor grid from power-of-2
   rounding (sampling cost ~ N sfunc calls; the current fft() forces
   N up to 2x). Extraction budget is loose (fft: ~0.2s, rotation: ~s),
   so poly-mult Bluestein (PARI's subquadratic RgX mul) is fine.

   Contract (matches exp-011b usage): given samples t[1..N], compute
     G[k+1] = sum_{j=0..N-1} t[j+1] * w^(j*k),  k=0..N-1
   for w = om^(-1) with om = exp(2*Pi*I/N)  (i.e. plain inverse-sign DFT).
   Bluestein: j*k = (j^2 + k^2 - (k-j)^2)/2, chirp ch = exp(-Pi*I/N):
     G_k = ch^(k^2) * sum_j (t_j ch^(j^2)) * ch^(-(k-j)^2)
   which is a convolution, done as poly mult. */

bluedft(t) = {
  local(N, ch, a, b, pp);
  N = length(t);
  ch = exp(-Pi*I/N);
  /* chirp powers: ch^(j^2) for j=0..N-1 (exponents mod 2N safe: use exact) */
  a = vector(N, j, t[j] * ch^((j-1)^2));
  b = vector(2*N-1, m, ch^(-(m-N)^2));  /* index m -> shift k-j = m-N */
  pp = Polrev(a) * Polrev(b);
  vector(N, k, ch^((k-1)^2) * polcoeff(pp, (k-1)+(N-1)));
}

/* reference: naive rotation DFT, same contract */
naivedft(t) = {
  local(N, w);
  N = length(t);
  w = exp(-2*Pi*I/N);
  vector(N, k, sum(j=1, N, t[j] * w^((j-1)*(k-1))));
}

/* self-test at several N (incl. non-powers of 2) and precisions */
{
  seed = 1;
  for (c = 1, 4,
    N = [12, 100, 1000, 2100][c];
    t = vector(N, j, cos(j*j*0.37 + c) + I*sin(j*1.71 - c));  /* deterministic pseudo-random */
    gettime();
    gb = bluedft(t);
    tb = gettime();
    gn = naivedft(t);
    tn = gettime();
    err = vecmax(vector(N, k, abs(gb[k]-gn[k])));
    scale = vecmax(vector(N, k, abs(gn[k])));
    print("N=", N, " relerr=", precision(err/scale, 5), " blue_ms=", tb, " naive_ms=", tn);
  );
}
