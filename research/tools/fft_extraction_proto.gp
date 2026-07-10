\\ exp-011 Prototyp v2: staylor-Extraktion (Rotations-Loop) vs. FFT — Identitaetstest
default(realprecision, 60);

{
  \\ ---------- KOMPLEXER BRANCH (complextaylor=1): n Samples, n Terme ----------
  n = 32; r = 0.7; rinv = 1/r;
  tcrc = vector(n, t, exp(Pi*I*(-1-1/n+2*t/n)));
  fsamp = vector(n, t, exp(0.3*tcrc[t]) + 1/(2 - 0.5*tcrc[t]) + 0.2*I*tcrc[t]^2);

  test = fsamp;
  old = vector(n);
  for (s=1, n,
    tot = 0;
    for (t=1, n, test[t] = test[t]*conj(tcrc[t]); tot = tot + test[t]);
    old[s] = tot/n * rinv^s;
  );

  om = exp(2*Pi*I/n);
  c0 = exp(-Pi*I*(1+1/n));
  winv = powers(om^(-1), n-1);
  G = fft(winv, fsamp);
  new = vector(n, s, (rinv^s/n) * conj(c0)^s * om^(-s) * G[(s%n)+1]);
  err1 = vecmax(vector(n, s, abs(old[s]-new[s])));
  print("komplex: max |alt-neu| = ", err1);

  \\ ---------- REELLER BRANCH (complextaylor=0): m Samples, 2m Terme ----------
  m = 16; terms = 2*m;
  tcrc2 = vector(m, t, exp(Pi*I*(-1/(2*m)+t/m)));
  fsamp2 = vector(m, t, exp(0.3*tcrc2[t]) + 1/(2 - 0.5*tcrc2[t]));

  test2 = fsamp2;
  old2 = vector(terms);
  for (s=1, terms,
    tot = 0;
    for (t=1, m, test2[t] = test2[t]*conj(tcrc2[t]); tot = tot + real(test2[t]));
    old2[s] = tot/m * rinv^s;
  );

  mu = exp(Pi*I/m);
  c1 = exp(-Pi*I/(2*m));
  w2inv = powers(mu^(-1), 2*m-1);
  padded = concat(fsamp2, vector(m, i, 0));
  H = fft(w2inv, padded);
  new2 = vector(terms, s, (rinv^s/m) * real(conj(c1)^s * mu^(-s) * H[(s%(2*m))+1]));
  err2 = vecmax(vector(terms, s, abs(old2[s]-new2[s])));
  print("reell:   max |alt-neu| = ", err2);
}
quit
