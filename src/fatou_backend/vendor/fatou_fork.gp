\ps 21
notheta=0;

matrixradius=0;
dbgmatrix=0;
limitp=0; /* default precision can be set to 16 to speed up calculations */
throwp=4.5; /* precision thrown away to improve convergence */
throws=9.5; /* precision at which nskip gets decremented every time */
z = 1.0;
precis=precision(z);
if (precis<=38, default(format,"g0.32"));
ir  = 7/10;  /* automatically selected values: 7/10; 15/16 */
ctr = 8/10;  /* automatically selected values: 8/10;  8/9  */
ctrmul = 1;  /* exp-023 scan knob: scales the theta sampling radius */
/* thlogk= 1/10; thlogk isn't used anymore :( see c:\pari\fatou_backup....gp for more details :( */
superfr = 0.1;
isuperfr = 0.01;
fmode=0;
invabeli=1/4;
complextaylor=1;
superfk=8;
efam=0; /* exp-011b: set by loop() — base-e family gets the fft extraction */
subeta=0; /* set by loop(): real base <= eta = e^(1/e), i.e. kc <= 0 — the
   walk lives at the real attracting fixed point there; the exp-032/033/034
   cache/quantization extensions broke that regime, so sub-eta bases take
   the exact pre-exp-032 code path (they are fast anyway). */
/* exp-012: the fs/finv walk in sfunc depends only on the base map, NOT on
   the evolving ct/theta state — cache walk endpoints per sample index
   across loop iterations while the sampling grid is unchanged. */
swon=0; swidx=0; swkey=0; swz=0; swn=0; swvalid=0; swb=0;
/* exp-041b: per-stretch cache of the theta branch's decision-only y2 */
swy2=0; swy2v=0;
/* exp-021: incremental ct-Horner state */
icct=0; icvals=0; icA=0; ickey=0; icdct=0; icdig=60; icfull=1;
/* exp-048: per-point degree truncation. ct's coefficients decay like
   circr^-k, so a sample at radius rho needs only ~dig/log10(circr/rho) terms.
   Samples are bucketed by radius (cached per grid stretch, like swz) and each
   level keeps a truncated copy of the polynomial, rebuilt once per pass. */
icnlev=8; swlev=0; icdl=0; icrmax=1; icdlon=0;
/* exp-024: isuperf/isuperf2 at the raw grid point depend only on the base
   map — cache per sample index within a grid stretch (flag 1/2 = branch) */
swisf=0; swisfv=0;
/* exp-026: theta-grid quantization + per-index superf cache + incremental
   abelest for thfunc (same program as exp-020/021/024, theta side) */
thon=0; thidx=0; thskey=0; thsw=0; thswv=0; thicv=0; thicA=0; thfull=1; thdct=0; thdig=60;
/* exp-027: incremental extraction — the DFT is linear, so transform only
   the t_est diff (scale ~10^-re) at reduced precision and add to the
   previous pass's coefficients */
exkey=0; extest=0; excoef=0;
/* exp-030: chirp caches — the extraction chirps depend only on the grid,
   not on the iteration; computing them as fresh powers every pass costs
   ~11 mults/term. Two-slot cache for bluedft (staylor+theta Ns), map
   cache for the staylor coefficient mapping. */
bd1n=0; bd1p=0; bd1t=0; bd1a=0; bd1b=0; bd1w=0; bd1m=0; bd2n=0; bd2p=0; bd2t=0; bd2a=0; bd2b=0; bd2w=0; bd2m=0; bdclock=1;
/* exp-061: mixdft twiddle cache (two slots, LRU, shares bdclock). The
   twiddles depend only on (nn, realprecision) -- never on the data -- but
   were rebuilt on every call. */
mx1n=0; mx1p=0; mx1t=0; mx1w=0; mx1o=0; mx2n=0; mx2p=0; mx2t=0; mx2w=0; mx2o=0;
exmap=0; exmapkey=0; exmapprec=0;   /* exp-062: precision the map carries */
quietmode=0;
/* I added || (real(Period)>47) to handle speed for sexpinit(1.4494); takes the place of theta0lim=0.224 */
theta0lim=0.0002; /* theta0lim=0.0224; */
x2mode=0; /* possibly use in the future */
etaB= exp(1/exp(1));
prtpoly(wtaylor,t,name) = {
  local(s,z);
  if (t==0,t=60);
  if (name<>0,
    print ("{"name"=");
    z=polcoeff(wtaylor,0);
    if (real(z)<0, print("      "z), print("        " z));
    for (s=1,t-1,
      z=polcoeff(wtaylor,s);
      if (s>9, print1("+x^" s), print1("+x^ " s));
      print1("* ");
      if (imag(z)<>0, print1("("));
      if (real(z)<0, print1(z), print1(" " z));
      if (imag(z)<>0, print(")"), print(););
    );
    print("}");
  ,
    for (s=0,t-1,
      z=polcoeff(wtaylor,s);
      if (s>9, print1("a" s "= "), print1("a" s "=  "));
      if (real(z)<0, print(z), print(" " z));
    );
  );
}
strip0fx(fx) = {
  local(fz,nt);
  fx=Pol(fx);
  fz=0;
  nt=length(fx);
  for (n=1,nt, fz=fz+polcoeff(fx,n)*x^n);
  return(Ser(fz));
}

finv(z,targ) = {
  local(zn,zn2);
  if (x2mode==0,
    if (abs(z-k+1)>0.000001, zn=log(z-k+1), zn=1E10);
  ,
    zn=sqrt(z-(k-1/4))-1/2;
    if (targ<>0,
      zn2=-(zn+0.5)-0.5;
      if (abs(zn-targ)>abs(zn2-targ),zn=zn2);
    );
  );
  return(zn);
}
fs(z) = {
  local(zn);
  if (x2mode==0,
    z=exp(z)-1+k;
  ,
    z=z^2+z+k;
  );
  return(z);
}
/* exp-071: safefs() is a SEARCH aid, not an evaluator. Above 5E8 it replaces
   the imaginary part with random(), and above |Re z|>10000 it returns the
   hardcoded sentinel 1E400*exp(I*imag(z)). Inside abel() that is correct and
   cheap -- the walk only needs to know "far away, keep going". But invabel()
   used it to BUILD the value it returns, so the sentinel became the answer:
   sexp_2(5) came back as 1E400/ln(2) = 1.4426950408889634E400 instead of
   2^65536 ~ 1E19728, and sexp_2(5.5) came back NEGATIVE, which is impossible
   for a real base-2 tetration -- non-deterministic, because safefs() puts a
   random() imaginary part in first; negative in the run that was measured.
   No error, no warning. See research/METHODS.md, section 2.
   The reconstruction loop in invabel() now uses the exact fs(); PARI raises
   e_OVERFLOW when the tower genuinely leaves the representable range, which
   is a loud failure instead of a plausible-looking wrong number.
   DO NOT route abel()'s walk through fs() -- the saturation is load-bearing
   there and removing it would cost the exp-0xx keeps their fast path. */
safefs(z) = {
  if (x2mode==0,
    if (abs(z)> 5E8, z=real(z)+I*random);
    if (abs(real(z))>10000,
      if (real(z)>0, z= 1E400*exp(I*imag(z)), z=-1+k; ),
      z=fs(z);
    );
  ,
    if (abs(z)<1E400, z=fs(z), z=1E800);
  );
  return(z);
}
kfromp(p) = {
  local(y,z);
  if (x2mode==0,
    z = 2*Pi*I/p;
    y = -exp(z) + z + 1;
  ,
    /* for the x^2+x+k case                                         */
    /* L=sqrt(-k)                                                   */
    /* lambda=1+2*L                                                 */
    /* Period = 2*Pi*I/log(lambda)                                  */
    /* z = log(z-L)/L + log(z-L2)/L2 + subst(e,x,(z-circc)*circ1r); */
    /* z = log(z-L)/log(lambda1) + log(z-L2)/log(lambda2) + subst(e,x,(z-circc)*circ1r);
    /* so 1/L; 1/L2 is replaced by rlnlm, rlnlm2                    */
    z = 2*Pi*I/p; /* log(lambda) */
    z = exp(z);   /* lambda      */
    z = (z-1)/2;  /* L           */
    y = -z^2;     /* k           */
  );
  return(y);
}

/* formalfixed.gp  */
/* prtpoly(formalfixed(16),17,"xfixed"); */
fixedk(kc) = {
  local(x0,n,fixed,y0);
  if ((imag(kc)==0) && (real(kc)>0),
    complextaylor=0;
    circc=real(circc);
  ,
    complextaylor=1;
  );
  if ((x2mode==1),
    k=kc;
    L=sqrt(-k);
    if (imag(L)<0,L=-L);
    L2=-L;
    /* for the x^2+x+k case                                         */
    /* L=sqrt(-k)                                                   */
    /* lambda=1+2*L                                                 */
    /* Period = 2*Pi*I/log(lambda)                                  */
    /* z = log(z-L)/L + log(z-L2)/L2 + subst(e,x,(z-circc)*circ1r); */
    /* z = log(z-L)/log(lambda1) + log(z-L2)/log(lambda2) + subst(e,x,(z-circc)*circ1r);
    /* so 1/L; 1/L2 is replaced by rlnlm, rlnlm2                    */
    rlnlm = 1/log(1+2*L);    /* recipricol of log(lambda) */
    rlnlm2 = 1/log(1+2*L2);  /* recipricol of log(lambda) */
    Period = 2*Pi*I*rlnlm;
    Period2 = 2*Pi*I*rlnlm2;
    if (real(Period2)<0,Period2=-Period2);
    repelling=1;
    if (imag(Period)<=0, repelling=0);
    circc=(L+L2)/2;           /* center; average of the fixed points   */
    circr=(L-circc)/I;        /* radius multiplier; L/L2 map to +I/-I  */
    argc = arg(circr);        /* for complex bases                     */
    circr=abs(circr);
    lambda1=1+2*L;
    lambda2=1+2*L2;
    return(L);
  );
xfixed=
        0
+x^ 1*  1
+x^ 2* -1/6
+x^ 3*  1/36
+x^ 4* -1/270
+x^ 5*  1/4320
+x^ 6*  1/17010
+x^ 7* -139/5443200
+x^ 8*  1/204120
+x^ 9* -571/2351462400
+x^10* -281/1515591000
+x^11*  163879/2172751257600
+x^12* -5221/354648294000
+x^13*  5246819/10168475885568000
+x^14*  5459/7447614174000
+x^15* -534703531/1830325659402240000
+x^16*  91207079/1595278956070800000;
  k=kc;
  lnb = exp(k-1);
  if (imag(kc)>=0,
    L  = subst(xfixed,x,+I*sqrt(2*kc));
    L2 = subst(xfixed,x,-I*sqrt(2*kc));
  ,
    L  = subst(xfixed,x,-I*sqrt(2*kc));
    L2 = subst(xfixed,x,+I*sqrt(2*kc));
  );
  y0 = 100;
  while ((abs(k)>10^-5) && (y0>10^(10-precis)),
    /* simple Newton's method ... */
    /* nearest singularity at kc=2*Pi*I */
    y0 = (exp(L)-1+kc-L)/(exp(L)-1);
    y1 = (exp(L2)-1+kc-L2)/(exp(L2)-1);
    L = L - y0;
    L2 = L2 - y1;
    y0 = abs(y0); y1=abs(y1);
    if (y1>y0,y0=y1);
  );
  rlnlm = 1/L;    /* recipricol of log(lambda) */
  rlnlm2 = 1/L2;  /* recipricol of log(lambda) */
  Period = 2*Pi*I*rlnlm;
  Period2 = 2*Pi*I*rlnlm2;
  if (real(Period2)<0,Period2=-Period2);
  repelling=1;
  if (imag(Period)<=0, repelling=0);
  circc=(L+L2)/2;           /* center; average of the fixed points   */
  circr=(L-circc)/I;        /* radius multiplier; L/L2 map to +I/-I  */
  argc = arg(circr);        /* for complex bases                     */
  circr=abs(circr);
  lambda1=L-k+1;
  lambda2=L2-k+1;
  l  = sexp_invabel(L);
  l2 = sexp_invabel(L2);
  return(L);
}

/* theta transform mapped to unit circle */
thfunc(z,n,yin,yon) = {
  local(y,p,h,A);
  if (n<>2,
    /* exp-059: the caller knows y exactly; this log only undoes its exp */
    y = if (yon, yin, log(z)/(2*Pi*I));
    if (thon && thidx>0,
      /* exp-026: superf(zth+y) depends only on the base map at the fixed
         grid point — cache it; abelest's ct-Horner updates incrementally
         against the previous pass (diff-poly thdct at reduced precision). */
      if (thswv[thidx],
        p = thsw[thidx];
      ,
        p = superf(zth+y);
        thsw[thidx] = p; thswv[thidx] = 1;
      );
      if (thfull,
        A = rlnlm*(log(I*(p-L))-Pi*I/2) + rlnlm2*(log(-I*(p-L2))+Pi*I/2) + sfunczero;
        h = subst(ct, x, (p-circc));
        thicA[thidx] = A; thicv[thidx] = h;
      ,
        h = thicv[thidx] + subst(thdct, x, precision(p-circc, thdig));
        thicv[thidx] = h;
        A = thicA[thidx];
      );
      y = (A + h) - y - zth;
      return(y);
    );
    y = abelest(superf(zth+y),ct)-y-zth;
    return(y);
  ,
    y = if (yon, -yin, log(z)/(-2*Pi*I));
    y = abelest(superf2(ztl+y),ct)-y-ztl;
    return(y);
  );
}
/* taylor series function ussing thfunc which is defined by sfuncmode */
thtaylor(n,samples) = {
  local(s,t,x1,y,z,tot,t_est,tcrc,halfsamples,wtaylor,terms,thsc,om,c0,G,cf,pw);
  if (samples==0, samples=120);  /* no matter how many sample points, the default gie series size is 200 halfsamples */
  samples = floor(samples);
  /* exp-026: quantize the theta grid (efam, n==1) to 64-steps so it stays
     unchanged for ~20-30 iterations — prerequisite for the caches below. */
  /* exp-034: theta quantization+caches for all REAL bases (complex bases
     alternate n=1/n=2 and would thrash the single-slot cache). */
  /* exp-051: exp-026 quantized only above 64, so below it the theta grid
     moved every iteration and thfull stayed 1 -- caches reallocated and the
     incremental path never engaged. Measured at dps 300: that was 60 of 232
     iterations, the contiguous block thsamples 27..62. Quantize there too,
     with a smaller step (16) so the overshoot stays small. */
  if ((subeta==0) && (efam || (complextaylor==0)) && (n==1),
    if (samples > 64,
      samples = 64*ceil(samples/64)
    ,
      if (samples > 16, samples = 16*ceil(samples/16))));
  terms=samples-1;
  t_est    = vector (samples,i,0);
  tcrc     = vector (samples,i,0);
  wtaylor=0;
  if ((subeta==0) && (efam || (complextaylor==0)) && (n==1),
    thfull = 1;
    if ((thskey == samples) && (type(icct) == "t_POL"),
      thdct = ct - icct;
      thsc = vecmax(apply(abs, Vec(thdct)));
      if (thsc > 0,
        thdig = max(30, default(realprecision) + ceil(log(thsc)/log(10)) + 40);
        thdct = precision(thdct, thdig);
        thfull = 0;
      );
    );
    if (thfull,
      thsw = vector(samples); thswv = vector(samples);
      thicv = vector(samples); thicA = vector(samples);
    );
    thskey = samples;
    thon = 1;
  );
  /* exp-059: geometric grid, and hand thfunc the exact y = x1/2 instead of
     making it recover x1 from exp(Pi*I*x1) with a full-precision log. */
  tcrc = geoseq(exp(Pi*I*(-1-1/samples)), exp(2*Pi*I/samples), samples);
  for(s=1, samples,
    x1 = -1 + -1/(samples) + (2*s/samples);
    thidx = s;
    t_est[s]= thfunc(tcrc[s], n, x1/2, 1);
  );
  thon = 0; thidx = 0;

  /* exp-042b: replace the quadratic rotation DFT with the transform
     staylor has used since exp-011b/020, and build the polynomial once.
       coeff_s = (1/n) sum_t t_est[t] conj(tcrc[t])^s
               = (1/n) (conj(c0)/om)^s G[(s%n)+1]
     om = exp(2*Pi*I/n), c0 = exp(-Pi*I*(1+1/n)),
     G = fft(powers(om^-1, n-1), t_est); bluedft for non-power-of-two n. */
  om = exp(2*Pi*I/samples);
  c0 = exp(-Pi*I*(1+1/samples));
  /* exp-061: doubling build at both sites. `pw` is the sharper of the two --
     it is not an FFT input, it multiplies the output coefficients one by
     one, so its ~terms ulp landed directly in the theta Taylor series. */
  if (2^valuation(samples,2)==samples,
    G = fft(geopow(om^(-1), samples-1), t_est)
  ,
    G = mixdft(t_est));
  pw = geopow(conj(c0)/om, terms);
  cf = vector(terms+1, j, pw[j] * G[((j-1)%samples)+1] / samples);
  wtaylor=Polrev(cf);
  wtaylor=precision(wtaylor,precis);
  return(wtaylor);
}

/* s(f(z-L)) = lambda*s(z-L) */
/* s(lambda x)=fx(s(x)) */
/* assumes fixed point of f(x)=0, fx(x)=f(x+L)-L */
formalischroder(fx,n) = {
  local(lambda,i,j,z,absz,plo,phi,f1t,f2t,ns,f1s);
  fx=strip0fx(fx);
  lambda = polcoeff(fx,1);
  f1t=Ser(x);
  z=1.0; ns=1; absz=1;
  plo=10^(-precision(z)/2.12);
  phi=10^(precision(z)/2.12);
  if (n==0, n=floor(precision(z)*1.49999));
  i=2;
  while (((i<=n) && ((ns==0) || ((absz>plo) && (absz<phi)))),
    f1s=f1t;
    f1t=f1t+acoeff*x^i;
    f2t=0;
    f2t=subst(f1t,x,lambda*x)-subst(fx,x,f1t);
    z = polcoeff(f2t, i);
    z = subst(z,acoeff,x);
    ns=-polcoeff(z,0)/polcoeff(z,1);
    f1t=f1s+ns*x^i;
    i++;
    if (i>3,absz=abs(ns));
  );
  return(Pol(f1t));
}

/* s(fx(z)) = lambda*s(z) */
/* assumes fixed point of f(x)=0, fx(x)=f(x+L)-L */
formalschroder(fx,n) = {
  local(lambda,i,z,absz,plo,phi,f1t,f2t,ns,f1s);
  fx=strip0fx(fx);
  lambda = polcoeff(fx,1);
  f1t=Ser(x);
  z=1.0; ns=1; absz=1;
  plo=10^(-precision(z)/2.1);
  phi=10^(precision(z)/2.1);
  if (n==0, n=floor(precision(z)*1.49999));
  i=2;
  while (((i<=n) && ((ns==0) || ((absz>plo) && (absz<phi)))),
    f1s=f1t;
    f1t=f1t+acoeff*x^i;
    f2t=subst(f1t,x,fx);
    z = polcoeff(f1t*(lambda)-f2t, i);
    z = subst(z,acoeff,x);
    ns=-polcoeff(z,0)/polcoeff(z,1);
    f1t=f1s+ns*x^i;
    i++;
    if (i>3,absz=abs(ns));
  );
  return(Pol(f1t));
}
isuperf(z) = {
  local(n,y,y1,u,u1);
  /* exp-057: iterate the OFFSET u = y - L instead of y. The old form walked y
     to within isuperfr = 10^(-precis/21) of L and only then formed y-L, which
     cancels log10(1/isuperfr) digits (measured: 21.4 at dps 400) -- that
     cancellation is the engine's whole calibration loss. In u-coordinates it
     never happens, using exp(L) = lambda1 = L-k+1:
         fs(L+u)   - L = lambda1*expm1(u)
         finv(L+u) - L = log1p(u/lambda1)
     x2mode has a different finv (sqrt branch) and keeps the original path. */
  if (x2mode==0,
    u = z - L;
    n = 0;
    if (repelling,
      while (abs(u)>isuperfr,
        /* review fix: finv() guards its log with abs(z-k+1)>1e-6 and returns
           1E10 otherwise. In u-coordinates z-k+1 is lambda1+u, so the same
           guard belongs here -- log1p(u/lambda1) returns -oo exactly where the
           original returned the sentinel. Reachable in principle: for base e
           lambda1 = L, so u = -lambda1 is the sample point z = 0, which lies
           inside the sampling circle. */
        u1 = if (abs(lambda1+u) > 0.000001, log1p(u/lambda1), 1E10-L);
        if (abs(u1+2*Pi*I)<abs(u1), u = u1+2*Pi*I, u = u1);
        n++;
      );
    ,
      while (abs(u)>isuperfr,
        u = lambda1*expm1(u);
        n--;
      );
    );
    y = subst(fsl,x,u);
  ,
  y=z;
  n=0;
  if (repelling,
    while (abs(y-L)>isuperfr,
      y1=finv(y,L);
      if (abs(y1+2*Pi*I-L)<abs(y1-L), y=y1+2*Pi*I, y=y1);
      n++;
    );
  ,
    while (abs(y-L)>isuperfr,
      y=fs(y);
      n--;
    );
  );
  y = subst(fsl,x,y-L);
  );
  y = log(y)*rlnlm + n;
  while (abs(imag(y))>abs(imag(y+Period)),y=y+Period);
  while (abs(imag(y))>abs(imag(y-Period)),y=y-Period);
  return(y);
}
isuperf2(z) = {
  local(n,y,u);
  /* exp-057: same offset iteration at the lower fixed point, lambda2 = L2-k+1 */
  if (x2mode==0,
    u = z - L2;
    n = 0;
    while (abs(u)>isuperfr2,
      /* review fix: same sentinel as in isuperf, at the lower fixed point */
      u = if (abs(lambda2+u) > 0.000001, log1p(u/lambda2), 1E10-L2);
      n++;
    );
    y = subst(fsl2,x,u);
  ,
  y=z;
  n=0;
  while (abs(y-L2)>isuperfr2,
    y=finv(y,L2);
    n++;
  );
  y = subst(fsl2,x,y-L2);
  );
  y = log(y)*rlnlm2 + n;
  while (abs(imag(y))>abs(imag(y+Period2)),y=y+Period2);
  while (abs(imag(y))>abs(imag(y-Period2)),y=y-Period2);
  return(y);
}
superf(z) = {
  local(n,i,y);
  if (repelling,
    y=lambda1^(z-superfk);
    n=superfk;
    while (abs(y)>superfr, y=y/(lambda1);n++);
    y = subst(fisl,x,y);
    for (i=1,n,y=fs(y));
    return(y);
  ,
    y=lambda1^(z+superfk);
    n=superfk;
    while (abs(y)>superfr, y=y*lambda1;n++);
    y = subst(fisl,x,y);
    for (i=1,n,y=finv(y);if (imag(y)<-1.5708,y=y+2*Pi*I));
    return(y);
  );
}
superf2(z) = {
  local(n,i,y);
  y=lambda2^(z-superfk);
  n=superfk;
  while (abs(y)>superfr2, y=y/lambda2;n++);
  y = subst(fisl2,x,y);
  for (i=1,n,y=fs(y));
  return(y);
}
z0rfunc(z) = {
  local(z1);
  z=exp(z*2*Pi*I)*ircircr;
  z=abs(z)/abs(fs(z+circc)-circc)-1;
  if (z>1, z=1);
  if (z<=-1, z=-1);
  return(z);
}
invz0rfunc(i) = {
  local(z1,z2,n,tlo,thi,t1,t2,tsave,p);
  if (i==0,i=43);
  tlo=0.0;
  while ((tlo<1) && (z0rfunc(tlo)>0), tlo=tlo+0.1);
  thi=tlo+0.1;
  while ((thi<(tlo+1)) && (z0rfunc(thi)<0), thi=thi+0.1);
  p = z0rfunc(tlo)*z0rfunc(thi);
  if (p>0,
    print(k);
    if (ir==7/10,
      print("theta init failed ir=7/10; no z0r found. try ctr=19/20; try ir=0.95");
      breakpoint();
    );
    if (ir==15/16,
      print("theta init failed; no z0r found. try ctr=19/20;");
      breakpoint();
    );
    print("theta init failed; no z0r found. try ctr=19/20;");
  );
  tsave=[tlo,thi];
  for (n=1,i,
    t1 = (tlo+thi)*0.5;
    z1 = z0rfunc(t1);
    if (z1<0,tlo=t1,thi=t1);
  );
  t1 = (tlo+thi)*0.5*2*Pi;
  z1 = (circc+exp(t1*I)*ircircr);
  tlo=tsave[1]+1;
  thi=tsave[2];
  for (n=1,i,
    t2 = (tlo+thi)*0.5;
    z2 = z0rfunc(t2);
    if (z2<0,tlo=t2,thi=t2);
  );
  t2 = (tlo+thi)*0.5*2*Pi;
  z2 = (circc+exp(t2*I)*ircircr);
  return([t1,t2,z1,z2]);
}

initsch(kd,myctr,myir) = {
  local(z1,z2,y1,y2,z,myps,arp,psold);
  z = 1.0;
  precis=precision(z);
  /* exp-060: the Schroeder walk runs until |u| <= isuperfr, and
     log10(1/isuperfr) = precis/myps exactly -- so the file's `\ps 21` pins the
     walk at ~precis/21 decades and, with it, a residual precision loss of the
     same order (exp-057 removed the larger, cancellation half of it). Scaling
     the series with the working precision buys those digits back, and the
     shorter walk pays for the longer series build. Measured at dps 300:
     ps 21 -> 38.45 s / 295.2 true digits, ps 38 -> 36.34 s / 301.9.
     Capped at 64 because formalschroder is ~O(d^4). The series-truncation
     depth is the GLOBAL seriesprecision (Ser(x) inside formalschroder), so it
     has to be raised as well and restored on the way out -- raising only myps
     fails with "degree > 21". */
  psold = default(seriesprecision);
  myps = max(psold, min(64, round(precis/8)));
  default(seriesprecision, myps);
  if (kd==0, kd=1);
  L=fixedk(kd);
  sfunczero=0; sfunczero=-abelest(circc,0);
  arp = abs(arg(Period));
  if ((ctr==1) || (arp<theta0lim) || (real(Period)>47), thetamode=0,thetamode=1);
  if (notheta, thetamode=0);
  if ((abs(k)>1.8) || ((repelling==0) && (argc<1.45)) || (arp<0.12), ir=15/16, ir=7/10);
  if ((real(k)>2.82)||(imag(k)>2),ctr=8/9,ctr=4/5);
  if (myctr<>0,ctr=myctr;ir=24/25);
  if (myir<>0,ir=myir);
  ircircr=ir*ctr*circr;
  z = invz0rfunc();
  z0h = z[3];
  z0l = z[4];
/*z0l = conj(z0h);  */ /* z0l=conj(z0h) for real valued bases */
  if (thetamode,
    fisl = L+formalischroder(fs(L+x)-L,myps);
    fsl = formalschroder(fs(L+x)-L,myps);
    superfr = log(10^-precis/abs(polcoeff(fisl,length(fisl)-1)))/(length(fisl)-1);
    superfr = exp(superfr);
    isuperfr = log(10^-precis/abs(polcoeff(fsl,length(fsl)-1)))/(length(fsl)-1);
    isuperfr = exp(isuperfr);
    zth = isuperf(z0h)+0.5;
    ztl = conj(zth);
    if (complextaylor,
      fisl2 = L2+formalischroder(fs(L2+x)-L2,myps);
      fsl2 = formalschroder(fs(L2+x)-L2,myps);
      superfr2 = log(10^-precis/abs(polcoeff(fisl2,length(fisl2)-1)))/(length(fisl2)-1);
      superfr2 = exp(superfr2);
      isuperfr2 = log(10^-precis/abs(polcoeff(fsl2,length(fsl2)-1)))/(length(fsl2)-1);
      isuperfr2 = exp(isuperfr2);
      ztl = isuperf2(z0l)+0.5;
    );
  );
/*ztl = isuperf2(z0l)+0.5; */
  default(seriesprecision, psold);   /* exp-060: restore */
  return(z);
}

abelest(z,e) = {
  z = rlnlm*(log(I*(z-L))-Pi*I/2) + rlnlm2*(log(-I*(z-L2))+Pi*I/2) + sfunczero + subst(e,x,(z-circc));
  if ((imag(z)==0) && (complextaylor==0), z=real(z));
  return(z);
}

renormr(ct) = {
  local(y0,y1,z0,z1,z,m2,s2,r2,n);
  z0=z0h;
  z1=fs(z0h);
  y0=abelest(z0,ct);
  y1=abelest(z1,ct);
  /* a*z1 + k*z1^2 - a*z0 + k*z0^2 = 1 - (y1-y0);  */
  /* a*(z1-z0) + k*(z1^2-z0^2) = 1 - (y1-y0)       */
  m2=matrix(2,2);
  s2=matrix(2,1);
  m2[1,1] = (z1-circc)-(z0-circc);
  m2[1,2] = (z1-circc)^2-(z0-circc)^2;
  s2[1,1] = 1 - (y1-y0);

  z0=z0l;
  z1=fs(z0l);
  y0=abelest(z0,ct);
  y1=abelest(z1,ct);

  m2[2,1] = (z1-circc)-(z0-circc);
  m2[2,2] = (z1-circc)^2-(z0-circc)^2;
  s2[2,1] = 1 - (y1-y0);
  r2=matsolve(m2,s2);
  z = r2[1,1]*x + r2[2,1]*x^2;

  z = z - abelest(circc,ct+z);

  if (complextaylor==0, z=real(z));
  return(z);
}

/* matrix solution to simulaneous equations without theta with lctr terms */
matrix_r(B,lctr) = {
  local(kc,ar,y0,yi,ys,z0,zi,zs,z,m2,s2,r2,i,j,nt,h,rr,re);
  /* a*zi + j*zi^2 - a*z0 + j*z0^2 = 1 - (yi-y0);  */
  /* a*(zi-z0) + j*(zi^2-z0^2) = 1 - (yi-y0)       */
  if (x2mode==0, kc=log(log(B))+1, kc=B);
  thetamode=0; ctr=1;
  initsch(kc);
  if (lctr==0, lctr=128);
  if ((lctr%4), lctr=lctr+4-(lctr%4));
  m2=matrix(lctr,lctr);
  s2=matrix(lctr,1);
  h = lctr/2;
  if (complextaylor<>0,nt=lctr,nt=h);
  for (i=1,nt,
    z0 = circc+circr*exp(I*argc+(i-0.5)*2*Pi*I/lctr);
    ar = arg(z0-circc)-argc;
    y0 = abelest(z0,0);
    zi = finv(z0);
    zs = fs(z0);
    if (((abs(ar)<Pi/2) && (abs(zi-circc)<circr)) || ((abs(zi-circc)<abs(zs-circc)) && (abs(zs-circc)>circr)),
      yi = abelest(zi,0);
      for (j=1,lctr, m2[i,j] = (zi-circc)^j-(z0-circc)^j);
      s2[i,1] = -1 - (yi-y0);
    ,
      ys = abelest(zs,0);
      for (j=1,lctr, m2[i,j] = (zs-circc)^j-(z0-circc)^j);
      s2[i,1] = 1 - (ys-y0);
    );
  );
  z = -abelest(circc,0);
  if (complextaylor==0,
    for (i=1,h,
      s2[i+h,1]=imag(s2[i,1]);
      s2[i,1]  =real(s2[i,1]);
      for (j=1,lctr,
        m2[i+h,j]=imag(m2[i,j]);
        m2[i,j]  =real(m2[i,j]);
      );
    );
    z = real(z);
  );
  r2=matsolve(m2,s2);
  for (i=1,lctr, z=z+r2[i,1]*x^i);
  ct = z;
  ctd = deriv(ct);
  rslog=renormslog(ct);
  z = abs(polcoeff(ct,lctr)*circr^lctr);
  rr=renormr(ct);
  re=abs(polcoeff(rr,1));
  if (re<z,re=z);
  re = -log(abs(re))/log(10);
  re = precision(re,9);
/*if (quietmode==0, printf("%5.1f decimal digits, %4d ctsamples\n", re, lctr); ); */
  if (quietmode==0, print(re " decimal digits, "lctr" ctsamples"); );
  return([re,lctr]);
}
log10(z)=log(z)*rln10;
rlog10(z)=log10(abs(z));
matrix_ir(B,lctr,ltht,myctr,myir) = {
  local(kc,yn,zn,ys,z,n,m,i,tot,tcrc,h,nt,savect,zi,ytht,rr,re,te,vthms,vthmt,vthtmt2,m2,s2,r2);
  /* matrix solution, sampled at ctr*circr, with inner points within radius of ir */
  /* assumes initsch is done */
  /* assumes superf/isuperf all that */
  /* assumes k is initialized */
  /* if ct, tht initialized, use those for length? */
  if (x2mode==0, kc=log(log(B))+1, kc=B);
  rln10=1/log(10);
  initsch(kc,myctr,myir);
  if ((matrixradius==0)&&(myir==0)&&((ir==7/10)||(ir==15/16)||(ir==24/25)),
    z = ir*ctr;
    ir = 24/25;
    ctr = (25/24)*z; /* shrink ctr by taking advantage of ir=24/25 */
  );
  if (lctr==0,
    if (limitp==0,z=precis-throwp,z=limitp);
    z=z-7.5;
    lctr= floor((-z*log(10))/log(ctr));
  );
  if ((lctr%4), lctr=lctr+4-(lctr%4));

  if (ltht==0, ltht=18);
  tcrc = vector(ltht);   /* vector of unit circle tht sample points */
  vthms = vector(ltht);  /* vector of th taylor series matrix samples */
  vthmt = vector(ltht);  /* vector of tht taylor series terms */
  m2=matrix(lctr,lctr);
  s2=matrix(lctr,1);
  if (complextaylor,
    vthmt2 = vector(ltht); /* vector of tht2 taylor series terms for superf2 */
  );
  rr=renormr(0);
  ct=rr;
  tht=thtaylor(1,6);
  if (complextaylor,
    tht2=thtaylor(2,6);
  );

  for (n=1,ltht,
    yn = -0.5+(n-0.5)/ltht;
    tcrc[n] = exp(2*Pi*I*yn);
    yn = zth+yn;
    zn=superf(yn);
    /* now we take the Abel function of zn -- corresponds to this line below in comments */
    /* y = abelest(superf(zth+y),ct)-y-zth; */
    vthms[n] = vector(lctr+1); /* lctr+1 holds the constant term */
    for (i=1,lctr,
      vthms[n][i] = (zn-circc)^i;
    );
    vthms[n][lctr+1] = abelest(zn,0)-yn;   /* this the right place for -yn */
  );
  /* fourier/taylor series for vthmt; fourier/taylor series for n=1 is the constant term; biased by n-1 */
  for (n=1,ltht,
    tot=0;
    for (i=1,ltht,
      tot=tot+vthms[i];
      vthms[i]=vthms[i]*conj(tcrc[i]);
    );
    tot=tot/ltht;
    vthmt[n]=tot;
  );

  if (complextaylor,
    for (n=1,ltht,
      yn = -0.5+(n-0.5)/ltht;
      yn = -yn; tcrc[n]=exp(-2*Pi*I*yn); /* Reverse Order for superf2 */
      yn = ztl+yn;
      zn=superf2(yn);
      /* now we take the Abel function of zn -- corresponds to this line below in comments */
      /* y = abelest(superf(zth+y),ct)-y-zth; */
      vthms[n] = vector(lctr+1); /* lctr+1 holds the constant term */
      for (i=1,lctr,
        vthms[n][i] = (zn-circc)^i;
      );
      vthms[n][lctr+1] = abelest(zn,0)-yn;   /* this is the right place for -yn */
    );
    /* fourier/taylor series for vthmt; fourier/taylor series for n=1 is the constant term; biased by n-1 */
    for (n=1,ltht,
      tot=0;
      for (i=1,ltht,
        tot=tot+vthms[i];
        vthms[i]=vthms[i]*conj(tcrc[i]);
      );
      tot=tot/ltht;
      vthmt2[n]=tot;
    );

  );
  if (quietmode==0, print("ctr="ctr", ir="ir", finished fourier vthmt; setting up main matrix next ..."));
  h = lctr/2;
  if (complextaylor<>0,nt=lctr,nt=h);
  for (n=1,nt,
    yn = circc+circr*ctr*exp(I*argc+(n-0.5)*2*Pi*I/lctr);
    zn = sfunc(yn); /* [est,n,est_point] */
    ys = zn[3];     /* the point from which the estimation for the sfunc was made */

    if (imag(zn[2])==0,
      for (m=1,lctr, m2[n,m] = (zn[3]-circc)^m-(yn-circc)^m);
      s2[n,1] = abelest(yn,0) - (abelest(ys,0)+zn[2]);
    ,
      if (imag(zn[2])>0,
        /* ys = isuperf(yn), use superf fourier series; ys is the point for the sueprfunction approximation */
        ytht = exp((ys-zth)*2*Pi*I);
        tot=0;
        for (i=1,ltht,
          tot=tot+vthmt[i]*ytht^(i-1); /* taylor series for n=1 is the constant term */
        );
      ,
        /* ys = isuperf2(yn), use superf2 fourier series; ys is the point for the sueprfunction approximation */
        ytht = exp((ys-ztl)*-2*Pi*I);
        tot=0;
        for (i=1,ltht,
          tot=tot+vthmt2[i]*ytht^(i-1); /* taylor series for n=1 is the constant term */
        );
      );
      zi = abelest(yn,0);
      for (m=1,lctr, m2[n,m] = tot[m]-(yn-circc)^m);
      s2[n,1] = -(ys+tot[lctr+1]-zi);   /* superf(z)+theta(superf(z)) */
    );
  );
  z = 0;
  if (complextaylor==0,
    for (n=1,h,
      s2[n+h,1]=imag(s2[n,1]);
      s2[n,1]  =real(s2[n,1]);
      for (m=1,lctr,
        m2[n+h,m]=imag(m2[n,m]);
        m2[n,m]  =real(m2[n,m]);
      );
    );
    z = real(z);
  );
  if (quietmode==0, print("setup main matrix; solving matrix now ..."));
  r2=matsolve(m2,s2);
  for (i=1,lctr, z=z+r2[i,1]*x^i);
  ct = z;
  ctd = deriv(ct);
  /* tht=thtaylor(1,ltht); */
  tht=0; thti=0;
  for (i=1,ltht,
    tot=0;
    for (n=1,lctr,tot=tot+polcoeff(ct,n)*vthmt[i][n]);
    tot=tot+vthmt[i][lctr+1];
    tht=tht+tot*x^(i-1);
  );
/*tht=thtaylor(1,ltht);*/
  if (complextaylor,
    /* tht2=thtaylor(2,ltht); */
    tht2=0; tht2i=0;
    for (i=1,ltht,
      tot=0;
      for (n=1,lctr,tot=tot+polcoeff(ct,n)*vthmt2[i][n]);
      tot=tot+vthmt2[i][lctr+1];
      tht2=tht2+tot*x^(i-1);
    );
  );
  thetamode=2; /* for other parts of fatou.gp */
  rslog=renormslog(ct);
  z = abs(polcoeff(ct,lctr)*(circr*ctr)^lctr);
  rr=renormr(ct);
  re=abs(polcoeff(rr,1));
  if (re<z,re=z);
  re = -log(abs(re))/log(10);
  re = precision(re,9);
  te = abs(polcoeff(tht,ltht-1));
  te = if (te<>0, te=-log(te)/log(10), te=precis);
  te = precision(te,9);
  if (quietmode==0, print(re " decimal digits, "te " theta digits"); );
  if (dbgmatrix==0,
    return([re,lctr,ltht]);
  ,
    if (dbgmatrix==1, return(s2); );
    if (dbgmatrix==2, return(m2); );
    if (dbgmatrix==3, return(vthms); );
    if (dbgmatrix==4, return(vthmt); );
  );
}

thest(z) = {
  local(y,y1,y2,zc,k0);
  if (abs(z-circc)<circr,
    y2 = abelest(z,ct);
    if (imag(y2) > 0,
      y1 = isuperf(z);
      y1 = y1+polcoeff(tht,0);
      while (abs(y1-y2)>abs(y1-y2+Period), y1=y1+Period);
      while (abs(y1-y2)>abs(y1-y2-Period), y1=y1-Period);
      y1=y1-polcoeff(tht,0);
      y1 = y1 + subst(tht,x,exp((y1-zth)*2*Pi*I));
    ,
      if (complextaylor==0, return(conj(thest(conj(z)))));
      y1 = isuperf2(z);
      y1 = y1+polcoeff(tht2,0);
      while (abs(y1-y2)>abs(y1-y2+Period2), y1=y1+Period2);
      while (abs(y1-y2)>abs(y1-y2-Period2), y1=y1-Period2);
      y1=y1-polcoeff(tht2,0);
      y1 = y1 + subst(tht2,x,exp((y1-ztl)*-2*Pi*I));
    );
  ,
    if (abs(z-L2)>abs(z-L),
      y1 = isuperf(z);
      if (imag(y1-zth)<0,
        if (imag(Period)<0, y1=y1-Period, y1=y1+Period);
      );
      y1 = y1 + subst(tht,x,exp((y1-zth)*2*Pi*I));
    ,
      if (complextaylor==0, return(conj(thest(conj(z)))));
      y1 = isuperf2(z);
      if (imag(y1-ztl)>0,
        if (imag(Period2)>0, y1=y1-Period2, y1=y1+Period2);
      );
      y1 = y1 + subst(tht2,x,exp((y1-ztl)*-2*Pi*I));
    );
  );
  return(y1);
}

renormslog(ct) = {
  local(n,y0,y1,z,t);

  y0=-1+k;
  n=-1;
  y1=fs(y0);
  if ((repelling==0) && (thetamode>1), t=superfr, t=0.1);

  while ((abs(y1-circc)<abs(y0-circc)) && (abs(y1-L)>t),
    y0=y1;
    y1=fs(y0);
    n++;
  );

  if ((abs(y0-circc)>circr*ctr) && (repelling==0) && (thetamode>1),
    z = -thest(y0) + n;
  ,
    z = -abelest(y0,ct) + n;
  );

  rslog=z;
  return(z);
}

/* exp-021: abelest(zz,ct) for sampling — the ct-Horner dominates warm
   sample cost (1.87ms of 1.88ms at dps 220 / 2561 terms). Between passes
   ct changes only by ~10^-re, so evaluate the DIFF polynomial at reduced
   precision (additive update, no cancellation) against a cached per-sample
   value. Full recompute whenever the grid changes (ickey mismatch). */
icabel(zz) = {
  local(h, A, lv);
  if (swon && swidx>0,
    /* exp-048: the sample's radius level is a property of the cached walk
       endpoint, so derive it once per grid stretch. */
    lv = 0;
    if (icdlon,
      lv = swlev[swidx];
      if (lv == 0,
        lv = 1 + floor(icnlev * abs(zz-circc) / icrmax);
        if (lv < 1, lv = 1);
        if (lv > icnlev, lv = icnlev);
        swlev[swidx] = lv;
      );
    );
    if (icfull,
      A = rlnlm*(log(I*(zz-L))-Pi*I/2) + rlnlm2*(log(-I*(zz-L2))+Pi*I/2) + sfunczero;
      h = subst(if (lv, icdl[lv], ct), x, (zz-circc));
      icA[swidx] = A;
      icvals[swidx] = h;
    ,
      h = icvals[swidx] + subst(if (lv, icdl[lv], icdct), x, precision(zz-circc, icdig));
      icvals[swidx] = h;
      A = icA[swidx];
    );
    h = A + h;
    if ((imag(h)==0) && (complextaylor==0), h=real(h));
    return(h);
  );
  return(abelest(zz, ct));
}

sfunc(z) = {
  local(y,y1,y2,zc,k0,n);
  if ((complextaylor==0) && (imag(z)<0), return(conj(sfunc(conj(z)))));
  zc=z;
  /* exp-015: y2 (a full O(terms) series evaluation) is only used by the
     theta branch below — compute it lazily there instead of eagerly. */
  /* exp-012: reuse the cached walk endpoint when sampling an unchanged
     grid (the walk uses only fs/finv — independent of ct/theta). */
  if (swon && swidx>0 && swvalid[swidx],
    z = swz[swidx];
    n = swn[swidx];
  ,
    y = fs(z);
    n=0;
    while (abs(y-circc)<abs(z-circc), z=y;y=fs(z);n--;);
    if (n==0,
      y=finv(z);
      while (abs(y-circc)<abs(z-circc), z=y; y=finv(y); n++);
    );
    if (swon && swidx>0,
      swz[swidx]=z; swn[swidx]=n;
      /* swvalid wird erst am Ende von sfunc gesetzt (exp-016: erst wenn
         auch swb gespeichert ist) */
    );
  );
  if ((abs(z-circc)<ircircr)||(thetamode==0),
    y1 = icabel(z) + n;
  ,
    z=zc;
    /* exp-015 lazy; exp-041b: y2 only feeds discrete decisions whose
       margins are O(1), so one evaluation per grid stretch suffices. */
    if (swon && swidx>0 && swy2v[swidx],
      y2 = swy2[swidx];
    ,
      y2 = icabel(zc);
      if (swon && swidx>0, swy2[swidx]=y2; swy2v[swidx]=1);
    );
    /* use theta mapping if abs(y-circc)>(ir*circr) */
    if (imag(y2)>0,
      if (swon && swidx>0 && swisfv[swidx]==1,
        y1 = swisf[swidx];
      ,
        y1 = isuperf(z);
        if (swon && swidx>0, swisf[swidx]=y1; swisfv[swidx]=1);
      );
      y1 = y1 + polcoeff(tht,0);
      while (abs(y1-y2)>abs(y1-y2+Period), y1=y1+Period);
      while (abs(y1-y2)>abs(y1-y2-Period), y1=y1-Period);
      y1=y1-polcoeff(tht,0);
      z = y1;
      n = I;
      y1 = y1 + subst(tht,x,exp((y1-zth)*2*Pi*I));
    ,
      if (swon && swidx>0 && swisfv[swidx]==2,
        y1 = swisf[swidx];
      ,
        y1 = isuperf2(z);
        if (swon && swidx>0, swisf[swidx]=y1; swisfv[swidx]=2);
      );
      y1 = y1 + polcoeff(tht2,0);
      while (abs(y1-y2)>abs(y1-y2+Period2), y1=y1+Period2);
      while (abs(y1-y2)>abs(y1-y2-Period2), y1=y1-Period2);
      y1=y1-polcoeff(tht2,0);
      z = y1;
      n = -I;
      y1 = y1 + subst(tht2,x,exp((y1-ztl)*-2*Pi*I));
    );
  );
  /* exp-016: abelest(zc,0) uses only the base estimate (ct-independent)
     at the raw grid point — cache it alongside the walk endpoints. */
  if (swon && swidx>0 && swvalid[swidx],
    y1 = y1 - swb[swidx];
  ,
    y = abelest(zc,0);
    if (swon && swidx>0, swb[swidx]=y; swvalid[swidx]=1);
    y1 = y1 - y;
  );
  return([y1,n,z]);
}

abel(z0) = {
  local(n,z1,z0r,z1r,argz,t,z);
  subetaguard("abel");
  z0r=abs(z0-circc);
  n=0;
  /* we can use absl(fs(z0))=abel(z0)+1, or abel(finv(z0))=abel(z0)-1 */
  /* the abel function is centered at circc, with a singularity at a radius of circr */
  /* the next 20 lines of code optimize mapping z0<=fs(z0); or z0<=finv(z0) so as to optimize convergence */
  /* argc is usually<1.45 except for complex bases as we approach the two real valued fixed point case */
  if (abs(z0-circc)>1E-20, argz=arg(z0-circc)-argc, argz=-argc);
  if ((repelling==0) && (thetamode>1), t=superfr, t=0.1);

  if ((argz>(Pi/2)) || (argz<(-Pi/2)),
    z1=safefs(z0);
    z1r=abs(z1-circc);
    /* the abs(z1-L)>t is for k<0 bases with argc>1.45 */
    while ((z1r<z0r) && (((n==0) && (z0r>circr)) || (argc<1.45) || (abs(z1-L)>t)),
      z0=z1;
      z0r=z1r;
      z1=fs(z0);
      z1r=abs(z1-circc);
      n++;
    );
    if ((n==0),
      z1=finv(z0);
      z1r=abs(z1-circc);
      while ((z1r<z0r) && ((argc<1.45) || ((abs(z1-L)>t) && (abs(z1-L2)>t))),
        z0=z1;
        z0r=z1r;
        z1=finv(z0);
        z1r=abs(z1-circc);
        n--;
      );
    );
  ,
    z1=finv(z0);
    z1r=abs(z1-circc);
/*  fix for half(z) near L1 */
/*  while ((z1r<z0r) && ((argc<1.45) || (abs(z1-L2)>t)) && (((n==0) && (z0r>circr)) || (abs(z1-L)>t)), */
    while ((z1r<z0r) && ((argc<1.45) || (abs(z1-L2)>t)) && ((((n==0) || (argc<1.45)) && (z0r>circr)) || (abs(z1-L)>t)),
      z0=z1;
      z0r=z1r;
      z1=finv(z0);
      z1r=abs(z1-circc);
      n--;
    );
    if ((n==0),
      z1=fs(z0);
      z1r=abs(z1-circc);
      while ((z1r<z0r) && ((argc<1.45) || (abs(z1-L)>t)),
        z0=z1;
        z0r=z1r;
        z1=fs(z0);
        z1r=abs(z1-circc);
        n++;
      );
    );
  );

  if ((z0r>circr*ctr) && (thetamode>1),
    z = thest(z0) - n;
  ,
    z = abelest(z0,ct) - n;
  );
  return(z);
}

/* ===== exp-071b: regulaere Iteration fuer reelle Basen 1 < b < eta =====
   Kneser braucht ein komplex-konjugiertes Fixpunktpaar. Fuer reelle
   1 < b < eta = e^(1/e) gibt es das nicht: dort liegt ein REELLER
   ANZIEHENDER Fixpunkt L mit b^L = L und Multiplikator lam = L*ln b in
   (0,1). Die eindeutige reelle Tetration kommt dort aus regulaerer
   Iteration (Koenigs/Schroeder), nicht aus der Kneser-Kontur:
       sigma(z)    = lim_n (f^n(z) - L) / lam^n
       sigma^-1(u) = lim_n f^-n(L + lam^n u),   f^-1(z) = ln z / ln b
       S(x)        = sigma^-1( sigma(1) * lam^x )   =>  S(0)=1, S(x+1)=b^S(x)

   Vorher liefen diese Basen trotzdem durch die Kneser-Maschinerie und gaben
   STILL ~15 richtige Stellen zurueck, wo 60 angefordert waren (b=1.2), ab
   dps ~100 nur noch Muell -- im Original genauso wie im Fork. Der Verraeter
   war ein Imaginaerteil 7.26e-14 auf einer beweisbar reellen Groesse.
   Belegt in research/METHODS.md Abschnitt 2; unabhaengige Referenz
   research/tools/regular_subeta.py.

   DISJUNKT ZUM KNESER-PFAD: alles hier laeuft nur bei subeta==1. Jede
   Gate-Basis hat subeta==0, keine exp-0xx-Optimierung wird beruehrt.

   PRAEZISIONSFALLE: regn und Arbeitspraezision sind NICHT unabhaengig
   waehlbar. sigma^-1 braucht lam^n*u oberhalb der Aufloesung von L, sonst
   geht die Information beim Addieren verloren; und f^-1 expandiert mit
   1/lam, verstaerkt also jeden Restfehler um lam^-n. Daher
   regwp = regn*log10(1/lam) + precis + 60. Mit zu kleiner Arbeitspraezision
   kommt S(0) = 14.77 statt 1 heraus (in der Python-Fassung selbst erlebt). */

regL=0; reglam=0; regs1=0; regn=0; reglb=0; regwp=0; regon=0;

regfix(b) = {
  local(z, lb, n, nit);
  lb = log(b);
  z = 3/2;
  nit = ceil(log(precis)/log(2)) + 12;
  for (n = 1, nit, z = z - (z - exp(lb*z))/(1 - lb*exp(lb*z)));
  return(z);
}

regsigma(z) = {
  local(v, n);
  /* Auf regwp aufpolstern, BEVOR iteriert wird. Sonst konvergiert v gegen L
     nur bis zur Praezision der EINGABE, und danach ist v-L bitgleich 0 --
     log(0) statt eines Ergebnisses. Das Aufpolstern erfindet keine Stellen:
     das Resultat traegt so viele wie die Eingabe, aber die Rechnung faellt
     nicht mehr auf null zusammen. */
  v = precision(z, regwp);
  for (n = 1, regn, v = exp(reglb*v));
  return((v-regL)/reglam^regn);
}

regsigmainv(u) = {
  local(v, n);
  v = regL + reglam^regn*u;
  for (n = 1, regn, v = log(v)/reglb);
  return(v);
}

reginit(b) = {
  local(sav, shrink, e0, e1, tol, zz, need, bp, lamr, Lr);
  sav = default(realprecision);
  zz = 1.0;
  precis = precision(zz);
  /* Durchgang 1 bei aktueller Praezision: nur um lam und damit den
     Praezisionsbedarf zu kennen. */
  reglb = log(b);
  Lr = regfix(b);
  lamr = Lr*reglb;
  if ((lamr <= 0) || (lamr >= 1), error("reginit: lam=", lamr, " nicht in (0,1) -- b liegt nicht in (1,eta)"));
  shrink = -log(lamr)/log(10);
  regn = floor((precis+25)/shrink)+1;
  need = floor(regn*shrink) + precis + 60;
  /* Die Basis selbst muss `need` Stellen tragen, und need ist rund das
     DOPPELTE des Ziels (regn*shrink ~ precis). Ein Literal wie 1.2 wird bei
     der gerade gesetzten Praezision geparst und traegt daher nie genug --
     mit einem 60-stelligen b liefert sigma^-1 nach der Verstaerkung um
     lam^-n = 10^85 puren Muell. Exakte Basen (6/5, t_INT, t_FRAC) tragen
     jede gewuenschte Stellenzahl. */
  bp = if (type(b) == "t_REAL", precision(b), need);
  if (bp < need,
    error("reginit: die Basis traegt nur ", bp, " Stellen, gebraucht werden ", need,
          ". Basis EXAKT uebergeben -- 6/5 statt 1.2, allgemein bestappr(b) --",
          " ein Dezimalliteral kann das prinzipiell nicht leisten.");
  );
  regwp = need;
  default(realprecision, regwp);
  iferr(
    reglb = log(b);
    regL = regfix(b);
    reglam = regL*reglb;
    regs1 = regsigma(1);
    /* exp-071d: der frueher hier stehende Test |sigma^-1(sigma(1)) - 1| war
       TAUTOLOGISCH -- er reduziert sich algebraisch auf f^-n(f^n(1)) = 1 und
       prueft nur, dass log und exp numerisch invers sind, nicht die
       Koenigs-Linearisierung. Jetzt zwei Schritte der Funktionalgleichung:
       S(1) = b und S(2) = b^b. Beide durchlaufen lam^1 bzw. lam^2 und sind
       damit echte Proben. */
    e0 = abs(regsigmainv(regs1*reglam^2) - b^b);
    e1 = abs(regsigmainv(regs1*reglam) - b)
  , E, default(realprecision, sav); error(E));
  default(realprecision, sav);
  tol = 10.0^(-precis);
  if ((e0 > tol) || (e1 > tol), error("reginit: Selbsttests halten nicht -- S(2)-b^b=", e0, "  S(1)-b=", e1));
  regon = 1;
  if (quietmode==0, print("subeta: regulaere Iteration, L=", precision(regL,20), " lam=", precision(reglam,20), " n=", regn, " Arbeits-dps=", regwp));
  return(regL);
}

regsexp(z) = {
  local(sav, y);
  if (regon==0, error("regsexp: reginit() wurde nicht ausgefuehrt"));
  sav = default(realprecision);
  default(realprecision, regwp);
  /* exp-071c: ohne iferr bliebe die Sitzung nach einem Fehler auf regwp
     stehen und ein spaeteres sexpinit() erbte die falsche Praezision. */
  y = iferr(regsigmainv(regs1*reglam^z), E, default(realprecision, sav); error(E));
  default(realprecision, sav);
  return(precision(y, precis));
}

regslog(y) = {
  local(sav, z);
  if (regon==0, error("regslog: reginit() wurde nicht ausgefuehrt"));
  sav = default(realprecision);
  default(realprecision, regwp);
  z = iferr(log(regsigma(y)/regs1)/log(reglam), E, default(realprecision, sav); error(E));
  default(realprecision, sav);
  return(precision(z, precis));
}

/* exp-071c: sexpinit() steigt fuer sub-eta VOR loop() aus, damit kein falscher
   Kneser-Zustand entsteht. Folge: k, lnb, rslog, ct, circc, circr, L, L2 und
   argc gehoeren dann noch zur ZUVOR initialisierten Basis. sexp() und slog()
   routen an die regulaere Iteration; alles andere wuerde still die Werte der
   vorigen Basis liefern -- nach sexpinit(exp(1)); sexpinit(6/5) etwa
   bitgleich die von e, ohne Fehler. Genau die Sorte stiller Falschantwort,
   die diese Version beseitigt. Daher hier ein lauter Abbruch.
   Blosses Zuruecksetzen von k reicht nicht: abel() haengt zusaetzlich an
   circc, circr, argc, ct, L und L2. */
subetaguard(nm) = {
  if (subeta,
    error(nm, "(): fuer reelle 1 < b < eta ist nur sexp() und slog() definiert. ",
          "Die Kneser-Groessen gehoeren noch zur zuvor initialisierten Basis und ",
          "wuerden hier still deren Werte liefern."));
}

slog(z) = {
  if (subeta, return(regslog(z)));
  abel(z*lnb+k-1)+rslog;
}

/* better invabel estimate; returns [est,p,curyz]; p=3 indicates betterest failed to converge */
betterest(z,est,n) = {
  local(y,slop,lest,ly,s,sn,p,lastyz,curyz,pgoal,absest);
  p=0;
  s=0;
  y=0;
  while ((s<n) && (p<3),
    /* can remove real(Period) if n=0 */
    if (abs(est-circc)>circr, p++; est=exp(arg(est-circc)*I)*circr+circc);
    s++;
    ly = y;
    y = abelest(est,ct);
    slop = subst(ctd,x,est-circc)+rlnlm/(est-L)+rlnlm2/(est-L2);
    lest = est;
    est=est+(z-y)/slop;
    curyz=abs(y-z);
/*  if ((real(Period)>2.2) && (abs(est-circc)>circr), p++; est=exp(arg(est-circc)*I)*circr+circc); */
  );
  if (p<3,
    if (n>0,
      ly=y;
      lastyz=curyz;
      absest=abs(est-lest);
    ,
      lastyz=100;
    );
    y=abel(est);
    curyz=abs(y-z);
    if (n==0,
      lest=est+curyz*0.0001;
      ly=abel(lest);
    );
    /* exp-083: Abbruchschranke 10^(-precis+2) statt 10^(-precis/1.2). Meist
       schiesst der letzte Sekantenschritt weit ueber die Schranke hinaus; landet
       er knapp darunter, blieb invabel auf ~precis/1.2 Stellen stehen. Das waren
       die sexp-Einbrueche an bestimmten Nachkommaanteilen (dps 60: schlechtester
       Punkt 2.5e-64 statt ~4e-74, jetzt 1.3e-73; ~25 % mehr Zeit je sexp-Wert). */
    pgoal=10^(-precis+2);
    absest=100;
    sn=s;
/*  while ((curyz>pgoal) && (absest>pgoal) && ((curyz<lastyz) || (s<3), */
    while ((curyz>pgoal) && (absest>pgoal) && ((curyz<lastyz) || (s<(sn+2))),
      est=precision(est,precis);
      y=precision(y,precis);
      if (abs(est-lest)>pgoal,
        slop=(y-ly)/(est-lest);
      );
      lest=est;
      ly=y;
      est=est+(z-y)/slop;
      lastyz=curyz;
      y=abel(est);
      curyz=abs(y-z);
      s++;
    );
    if (curyz>0.1, p=3);
/*  if (curyz>0.1, print (z" "curyz " need better initial est, invabel(z)")); */
  );
  return([est,p,curyz]);
}

/* exp-092: ct am Rand seines Abtastkreises auswerten -- am Fixpunkt L bzw. L2,
   |L - circc| = circr. Die Groesse ist nur die Startschaetzung fuer betterest.
   Bei dps 40 (Grad 256) fallen die Glieder dort bis ~1e-5; bei dps 60 (Grad
   512) wachsen sie auf 1e23, exp(y/rlnlm) lief ueber, und sexp(z) brach fuer
   jedes z mit Im z > invabeli ab ("overflow in expo()", b = e, z = 0.5+0.3i;
   gefunden von der Wissensgraph-Sitzung). Ist der Schwanz klein, wie bisher die
   volle Summe; sonst die Summe bis zum kleinsten Glied der Huellkurve
   (optimale Abbruchstelle einer Reihe am Rand ihres Kreises). */
ctedge(w) = {
  my(n = poldegree(ct), t, tail, env, jmin);
  t = vector(n + 1, j, polcoeff(ct, j - 1)*w^(j - 1));
  tail = vecmax(vector(min(10, n + 1), j, abs(t[n + 2 - j])));
  if (tail < 1e-3, return(vecsum(t)));
  env = vector(n + 1, j, vecmax(vector(min(5, n + 2 - j), i, abs(t[j + i - 1]))));
  jmin = 1; for (j = 2, n + 1, if (env[j] < env[jmin], jmin = j));
  return(vecsum(t[1..jmin]));
}
invabel(z,est) = {
  local (t,tc,y,curyz,nest,estp,estn);
  subetaguard("invabel");
  if (est==0,
    t = real(z);  /* biased around zero */
    if ((x2mode && (real(k)>0.5)) , tc=-0.5, tc=0);  /* new equation added this week; not online */
    if (argc<1.45,
      if (imag(z)<0,
        tc = tc + imag(Period2)*imag(z)/real(Period2);
      ,
        tc = tc +imag(Period)*imag(z)/real(Period);
      );
    ,
      /* need better guestimate for imag(z)>imag(Period) */
      /* function becomes ~= Period*I periodic; likewise for -Period2 */
      tc = tc-10*imag(z)/abs(0.5*imag(Period)); /* guestimate */
    );
    t = floor(t+tc+0.5);
    z = (z-t);
/*  print (tc" "t" "z);  */

    if ((imag(z)>invabeli),
      /* abelest(z) = log(z+k/2-L)/L + log(z+k/2-L2)/L2 + subst(e,x,z); */
      /* abel(z) = abelest(z-k/2,ct) - 1/2                             */
      /* abel(z) = log(z-L)*rlnlm + log(z-L2)*rlnlm2 + subst(e,x,z-k/2) - 1/2   */

      y = z - log(L-L2)*rlnlm2 - ctedge(L-circc) - sfunczero;
      est = exp(y/rlnlm)+L;

      if (argc<1.45, nest=betterest(z,est,3), nest=betterest(z,est,0)); /* use old algorithm for argc>1.45, for theta mapping */
      if (nest[2]>2,  /* failed; try invabel(z+/-1) */
        estp = exp((y+1)/rlnlm)+L;
        estn = exp((y-1)/rlnlm)+L;
        nest = betterest(z+1,estp,3);
        if (nest[2]<3, z++;t--,  /* estp worked, else try estn */
          nest = betterest(z-1,estn,3);
          if (nest[2]<3, z--;t++,  /* estn worked */
            nest=betterest(z,est,0);  /* everything failed; try old algorithm w/out error limit; this helps for period<2 */
          );
        );
      );


    , if ((imag(z)<-invabeli),

      y = z - log(L2-L)*rlnlm - ctedge(L2-circc) - sfunczero;
      est = exp(y/rlnlm2)+L2;

      if (argc<1.45, nest=betterest(z,est,3), nest=betterest(z,est,0)); /* use old algorithm for argc>1.45, for theta mapping */
      if (nest[2]>2,    /* failed; try invabel(z+/-1) */
        estp = exp((y+1)/rlnlm2)+L;
        estn = exp((y-1)/rlnlm2)+L;
        nest = betterest(z+1,estp,3);
        if (nest[2]<3, z++;t--,  /* estp worked, else try estn */
          nest = betterest(z-1,estn,3);
          if (nest[2]<3, z--;t++,  /* estn worked */
            nest=betterest(z,est,0);  /* everything failed; try old algorithm w/out error limit; this helps for period<2 */
          );
        );
      );

    ,
      est=circc;
      nest=betterest(z,est,3);
    ););
  );
  est = nest[1];
  curyz=nest[3];

  if (curyz>0.1, print (z" "curyz " need better initial est, invabel(z)"));
  if (t>0,
    /* exp-071: exact map, NOT safefs -- see the note at safefs(). */
    for (n=1,t,est=fs(est));
  ,
    /* hack for x2mode symmetry */
    if ((x2mode==1) && (abs(est-circc)>circr),
      if (abs(abelest(est,ct)-(z-t))>abs(abelest(-1-est,ct)-(z-t)), est=-1-est);
    );
    for (n=1,-t, est=finv(est); );
  );
  if ((imag(z)==0) && (complextaylor==0), est=real(est));
  return(est);
}

abeltaylor( w,r,samples,invabelmode) = {
  local(rinv,s,t,x1,z,tot,t_est,tcrc,halfsamples,wtaylor,terms);
/* outputs gie taylor series, the complex taylor series for invabel, */
  if (samples==0, samples=240);  /* no matter how many sample points, the default gie series size is 200 halfsamples */
  halfsamples=samples/2;
  terms = floor(samples*200/240);
  t_est    = vector (samples,i,0);
  tcrc     = vector (samples,i,0);
  if (r==0,r=1);
  rinv = 1/r;
  wtaylor=0;
  for(s=1, samples, x1=-1/(samples)+(s/halfsamples); tcrc[s]=exp(Pi*I*x1); );

  /* uses the invabel to allow for centering anywhere in the complex plane */
  for (t=1,samples,
    z = w+r*tcrc[t];
    if (invabelmode==0, t_est[t]=abel(z), if (invabelmode==1, t_est[t]=invabel(z), t_est[t]=sexp(z)) );
  );




  for (s=0,terms-1,
    tot=0;
    for (t=1,samples,
      tot=tot+t_est[t];
      t_est[t]=t_est[t]*conj(tcrc[t]);
    );
    tot=tot/samples;
    if (s>=1, tot=tot*(rinv)^s);
    wtaylor=wtaylor+tot*x^s;
  );
  wtaylor=precision(wtaylor,precis);
  return(wtaylor);
}

invabeltaylor(w, r, samples) = {
  return(abeltaylor(w,r,samples,1));
}

sexp(z) = {
  local (y);
  if (subeta, return(regsexp(z)));  /* exp-071b */
  /* center inflection point at the origin */
  y = invabel(z-rslog);
  y = (y-k+1)/lnb;
  if ((imag(z)==0) && (complextaylor==0), y=real(y));
  return(y);
}

sexptaylor(w, r, samples) = {
  local(y);
  subetaguard("sexptaylor");
  y = abeltaylor(w-rslog,r,samples,1);
  y = (y-k+1)/lnb;
  return(y);
}

slogtaylor(w, r, samples) = {
  local(y,z0);
  subetaguard("slogtaylor");
  y = abeltaylor(w*lnb+k-1,r*lnb,samples) + rslog;
  y = subst(y,x,x*lnb);
  return(y);
}

sexp_invabel(z) = {
  local(rlnB);
  /* lnB=exp(-1+k);         */
  /* y = y - (-1+k);        */
  /* y = y/lnB;             */
  /* lnb = exp(k-1);        */
  z = (z - (-1+k))/lnb;
  return(z);
}

invabel_sexp(z) = {
  local(rlnB);
  /* lnB=exp(-1+k);         */
  /* y = y - (-1+k);        */
  /* y = y/lnB;             */
  /* lnb = exp(k-1);        */
  z = (z*lnb)+(-1+k);
  return(z);
}

/* exp-020: exact-N DFT via Bluestein chirp + poly mult, for grids that
   are not powers of 2. Verified vs naive rotation DFT to precision floor
   (research/tools/bluestein_proto.gp; N=12/100/1000/2100, ~53x faster
   than rotation at N=2100). */
/* exp-056: exact-N DFT for N = r*2^k with small odd r, by radix-r decimation
   into r power-of-two FFTs plus twiddles. Same convention as bluedft:
   X[k] = sum_j t[j] * omega^((j-1)(k-1)),  omega = exp(-2*Pi*I/N).
   Falls back to bluedft for anything that does not factor that way. */
mixdft(t) = {
  local(nn, m, r, w, sub, s, a, j, omp, res, idx, acc, sl);
  nn = #t;
  m = 2^valuation(nn, 2);
  r = nn/m;
  if ((r > 9) || (m < 8), return(bluedft(t)));
  /* exp-061: both twiddle vectors depend only on (nn, precision), so cache
     them (two slots, LRU) instead of rebuilding ~nn full-precision products
     on every call, and build them by doubling -- powers() carries ~n ulp
     (measured 4.6e-402 at nn=4608, precis 404, i.e. ~2.7 digits; doubling
     gives 1.9e-403). The exp-056 verification could not see that error:
     bluedft built its twiddles the same way, so the two agreed on a
     common-mode error. */
  sl = 0;
  if ((mx1n == nn) && (mx1p >= default(realprecision)), sl = 1);
  if ((sl == 0) && (mx2n == nn) && (mx2p >= default(realprecision)), sl = 2);
  if (sl == 0,
    w = geopow(exp(-2*Pi*I/m), m-1);
    omp = if (r == 1, 0, geopow(exp(-2*Pi*I/nn), nn-1));
    if (mx1t <= mx2t,
      mx1n=nn; mx1p=default(realprecision); mx1w=w; mx1o=omp; mx1t=bdclock++; sl=1
    ,
      mx2n=nn; mx2p=default(realprecision); mx2w=w; mx2o=omp; mx2t=bdclock++; sl=2);
  );
  if (sl == 1,
    w = mx1w; omp = mx1o; mx1t=bdclock++
  ,
    w = mx2w; omp = mx2o; mx2t=bdclock++);
  sub = vector(r, s, fft(w, vector(m, a, t[(a-1)*r + s])));
  if (r == 1, return(sub[1]));
  res = vector(nn);
  for (j=1, nn,
    idx = ((j-1) % m) + 1;
    acc = sub[1][idx];
    for (s=2, r,
      acc = acc + omp[(((j-1)*(s-1)) % nn) + 1] * sub[s][idx];
    );
    res[j] = acc;
  );
  return(res);
}

bluedft(t) = {
  local(nn, ch, M, w, fa, pp, sl, chv, fbv);
  nn = length(t);
  sl = 0;
  if ((bd1n == nn) && (bd1p >= default(realprecision)), sl = 1);
  if ((sl == 0) && (bd2n == nn) && (bd2p >= default(realprecision)), sl = 2);
  if (sl == 0,
    /* exp-031: convolution via pow2 FFT instead of poly mult — Karatsuba
       temporaries blew the PARI stack at dps 1020 / grid 2304 (anchor
       crash "not enough memory"); FFT is also 2-3x faster and the b-side
       transform is grid-fixed, so it is cached alongside the chirps. */
    ch = exp(-Pi*I/nn);
    chv = vector(nn, j, ch^((j-1)^2));
    M = 1; while (M < 3*nn-2, M = M*2);
    /* exp-061: doubling build -- M is the next power of two above 3*nn-2,
       so powers() carried up to ~4*nn ulp here, the largest of the three
       twiddle sites. Cached since exp-031, so the extra build cost is paid
       once per grid. */
    w = geopow(exp(2*Pi*I/M), M-1);
    fbv = fft(w, concat(vector(2*nn-1, m, ch^(-(m-nn)^2)), vector(M-(2*nn-1), i, 0)));
    if (bd1t <= bd2t,
      bd1n=nn; bd1p=default(realprecision); bd1a=chv; bd1b=fbv; bd1w=w; bd1m=M; bd1t=bdclock++; sl=1
    ,
      bd2n=nn; bd2p=default(realprecision); bd2a=chv; bd2b=fbv; bd2w=w; bd2m=M; bd2t=bdclock++; sl=2);
  );
  if (sl == 1,
    chv = bd1a; fbv = bd1b; w = bd1w; M = bd1m; bd1t=bdclock++
  ,
    chv = bd2a; fbv = bd2b; w = bd2w; M = bd2m; bd2t=bdclock++);
  fa = fft(w, concat(vector(nn, j, t[j] * chv[j]), vector(M-nn, i, 0)));
  pp = fftinv(w, vector(M, i, fa[i]*fbv[i]));
  /* exp-071d: Laufvariable NICHT k nennen -- k ist die Basiskonstante der
     Engine (fs(z) = e^z - 1 + k); siehe die NOTE in icbuild(). */
  vector(nn, i4, chv[i4] * pp[(i4-1)+(nn-1)+1] / M);
}

/* exp-048: build one truncated copy of pp per radius level.
   Truncation index from the actual coefficient magnitudes via exponent()
   (binary exponents, no logs): keep term k while
       exponent(c_k) + k*log2(rho)  >=  exponent(max|c|) - dig*log2(10) - 16.
   The level radius is the UPPER edge of the bucket, so every sample assigned
   to that level is covered. */
icbuild(pp, dig) = {
  /* NOTE: the loop variable must NOT be called k -- k is the base parameter
     of this engine (fs(z) = e^z - 1 + k), and using it here silently defeats
     the truncation. */
  local(v, n, lg2, tol2, j, lr2, kk, ii, mx, bv, icsq);
  icdlon = 0;
  if (type(pp) != "t_POL", return(0));
  v = Vecrev(pp);
  n = #v;
  if (n < 64, return(0));
  lg2 = vector(n, kk, if (v[kk]==0, -1000000, exponent(v[kk])));
  mx = vecmax(lg2);
  if (mx <= -1000000, return(0));
  icsq = vector(n, kk, kk-1);   /* exp-074b: einmal, fuer die Vektorform von bv */
  icdl = vector(icnlev);
  for (j=1, icnlev,
    lr2 = log(icrmax*j/icnlev)/log(2);
    /* exp-074: die Schwelle haengt am groessten BEITRAG |c_k|*rho^k DIESER
       Stufe, nicht am groessten Koeffizienten mx. Beides ist dasselbe, solange
       die Koeffizienten fallen -- dann liegt max(lg2) bei k=1 und (k-1)*lr2=0.
       Genau das setzt der Kommentar oben voraus ("coefficients decay like
       circr^-k"), und das gilt fuer circr > 1. Bei Basis e ist circr = 1.3372,
       am Shell-Thron-Rand aber ~0.133: dort WACHSEN die Koeffizienten, mx sitzt
       am hoechsten Index (gemessen 160 von 161), der groesste Beitrag bei
       Index 2, und die alte Schwelle lag ~85 Dezimalstellen zu hoch. kk fiel
       von 161 auf 45, sfunc lieferte Muell, die Iteration stagnierte sofort:
       b=1.4494 gab bei dps 150 nur 1.837 statt 147.410 Stellen, ohne Fehler.
       Belegt in research/METHODS.md Abschnitt 2; Gate danach BITGLEICH. */
    /* exp-074b: bv als VEKTOR-Operation, nicht als GP-Schleife.
       Die erste Fassung suchte das Maximum mit `for (i3=1,n, ...)` und kostete
       gemessen ~10 % auf Basis e (dps 150: 4133 -> 4531 ms; dps 300:
       36258 -> 40360 ms, Klammern 0.3-2.1 %). Mein Ueberschlag "8n Vergleiche
       gegen N^2 Auswertungen, also 0.1 %" war doppelt falsch: die N^2-Arbeit
       steckt in kompilierten PARI-Routinen, diese Schleife nicht -- zaehlbare
       Operationen sind keine vergleichbaren Operationen. Und ein Waechter
       "wenn mx vorne liegt, ist bv = mx" hilft nicht: bei Basis e liegt das
       Maximum bei Index 32 bzw. 64 von 65, die Koeffizienten wachsen dort
       also EBENFALLS. Die alte Schwelle war damit auch fuer Basis e zu scharf
       (bv = -24.2 gegen mx = -7), nur folgenlos; am Rand wurde daraus der
       Totalausfall.
       lg2 + lr2*icsq ist Vektorarithmetik und vecmax ein Builtin -- pro Stufe
       drei kompilierte O(n)-Durchgaenge statt n interpretierter Schritte. */
    bv = vecmax(lg2 + lr2*icsq);
    tol2 = bv - dig*3.3219280948873623 - 16;
    kk = n;
    while ((kk > 1) && (lg2[kk] + (kk-1)*lr2 < tol2), kk--);
    /* build the truncation explicitly from the coefficient vector: the
       polynomial Euclidean quotient does not truncate a t_REAL polynomial
       the way one expects. */
    icdl[j] = if (kk >= n, pp, Polrev(vector(kk, ii, v[ii])));
  );
  icdlon = 1;
  return(1);
}

/* exp-061: [1, mu, mu^2, ..., mu^n] -- exactly what powers(mu,n) returns,
   but built by doubling, so the error is ~log2(n) ulp instead of ~n, and
   with no division (the first form of geoseq divided by c once per element,
   and at these precisions a division costs 2-3x a multiplication). */
geopow(mu, n) = {
  local(v, half, j, blk);
  if (n < 0, return(vector(0)));   /* guard: vector(n+1) would raise */
  v = vector(n+1);
  v[1] = 1;
  if (n == 0, return(v));
  v[2] = mu;
  blk = 1;
  while (blk < n,
    half = min(blk, n-blk);
    for (j=1, half, v[blk+j+1] = v[blk+1] * v[j+1]);
    blk = blk + half;
  );
  return(v);
}

/* exp-059: tcrc[s] = c*mu^s for s = 1..n, now via geopow (no division). */
geoseq(c, mu, n) = {
  local(v);
  if (n <= 0, return(vector(0)));
  v = geopow(mu, n);
  return(vector(n, j, c*v[j+1]));
}

staylor( w,r,samples) = {
  /* exp-063: exdig/exin/exp0/exsc/icsc were assigned but not declared, so
     they leaked into the global namespace and rode along in every dumped
     engine state. They are temporaries of the extraction block. */
  local(rinv,s,t,x1,y,y0,y1,y2,st,z,tot,t_est,tcrc,halfsamples,wtaylor,terms,om,c0,mu,c1,G,coeffs,
        exdig,exin,exp0,exsc,icsc);
  if (samples==0, samples=240);  /* no matter how many sample points, the default gie series size is 200 halfsamples */
  /* exp-011b: the fft extraction needs a power-of-2 grid; the grid change
     costs the factor-2 bases their delicate sample/terms co-evolution
     (gate: 30-43 digits), so it applies ONLY to the base-e family where
     the expensive long inits live. */
  /* exp-020: 256-step grids instead of power-of-2 (waste <=12% instead
     of up to 2x in sfunc sampling, the dominant cost); the exp-012 walk
     cache stays valid within each step. Extraction routes to bluedft()
     when the grid is not a power of 2. */
  if (efam,
    if (samples <= 256,
      z=1; while (z<samples, z=z*2); samples=z
    ,
      /* exp-037 (KEEP-Kandidat): depth-adaptive quantum — above grid 3584
         cold passes (walks + full Horner/extract, ~12x a warm iteration)
         dominate; 512-steps halve the stretch changes. Verified in the
         720/733 e|700 reference pair: 1.10x normalized, full digits
         (698 proven). Below 3584 bit-identical to the 256-quantum. */
      if (samples > 3584,
        samples = 512*ceil(samples/512)
      ,
        samples = 256*ceil(samples/256)))
  ,
    /* exp-032: gentle 32-step quantization for non-e bases — stabilizes
       the grid for a few iterations so the walk/isuperf caches and the
       incremental Horner apply. Small steps to respect the delicate
       sample/terms co-evolution that pow2 rounding broke (exp-011b). */
    if ((samples > 64) && (subeta==0), samples = 32*ceil(samples/32)));
  terms=samples;
  if (complextaylor==0, samples=samples/2);
  t_est    = vector (samples,i,0);
  tcrc     = vector (samples,i,0);
  if (r==0,r=1);
  rinv = 1/r;

  /* exp-012: enable the sfunc walk cache while sampling an unchanged grid
     (exp-032: all bases except sub-eta; grid identity = [samples, w, r]). */
  if (subeta==0,
    if (swkey != [samples, w, r],
      swkey = [samples, w, r];
      swz = vector(samples);
      swn = vector(samples);
      swb = vector(samples);
      swvalid = vector(samples);
      swisf = vector(samples);
      swisfv = vector(samples);
      swlev = vector(samples);   /* exp-048: per-sample radius level */
      swy2 = vector(samples);
      swy2v = vector(samples);
    );
    swon = 1;
    /* exp-021: prepare incremental pass — diff of ct vs the previously
       sampled ct, evaluated at reduced precision (digits needed = working
       precision minus the diff's scale, +40 margin). */
    icfull = 1;
    if ((ickey == [samples, w, r]) && (type(icct) == "t_POL"),
      icdct = ct - icct;
      icsc = vecmax(apply(abs, Vec(icdct)));
      if (icsc > 0,
        icdig = max(30, default(realprecision) + ceil(log(icsc)/log(10)) + 40);
        icdct = precision(icdct, icdig);
        icfull = 0;
      );
    );
    if (icfull,
      icvals = vector(samples);
      icA = vector(samples);
    );
    ickey = [samples, w, r];
    icct = ct;
    /* exp-048: (re)build the per-level truncations for this pass */
    icrmax = r;
    icbuild(if (icfull, ct, icdct), if (icfull, precis, icdig));
  );
  if (complextaylor,
    /* exp-059: x1 is affine in s, so tcrc is geometric */
    tcrc = geoseq(exp(Pi*I*(-1-1/samples)), exp(2*Pi*I/samples), samples);
    y0 = exp(I*argc);
    for(s=1, samples,
      swidx = s;
      t_est[s]=sfunc(w+r*tcrc[s]*y0)[1];
    );
  ,

    /* exp-059: x1 is affine in s, so tcrc is geometric */
    tcrc = geoseq(exp(-Pi*I/(2*samples)), exp(Pi*I/samples), samples);
    for(s=1, samples,
      swidx = s;
      t_est[s]=sfunc(w+r*tcrc[s])[1];
    );
  );
  swon = 0;
  swidx = 0;
  wtaylor=0;
  y0=0;
  y1=0;
  y2=0;
  st=0;
  tot=0;
  /* exp-011b: for the base-e family the rotation extraction (a plain DFT,
     O(terms*samples)) is replaced by PARI's fft on the power-of-2 grid:
       complex: coeff_s = (rinv^s/n) * conj(c0)^s * om^(-s) * G[(s%n)+1]
       real:    coeff_s = (rinv^s/m) * Re(conj(c1)^s * mu^(-s) * G[(s%(2m))+1])
     Verified identical to the rotation loop to machine precision
     (research/tools/fft_extraction_proto.gp). exp-033: since exp-032 the
     grid is quantized (32-steps) for ALL bases and bluedft handles exact N,
     so every base takes this path now — the exp-011b fragility came from
     pow2 grid doubling, not from the transform. The rotation loop below
     stays live for sub-eta bases (attracting-fixed-point regime). */
  if (subeta==0,
    /* exp-027: on an unchanged grid, extract only the diff of the samples
       at reduced precision (the transform is linear); full extract at
       stretch starts. exdig from the diff's scale like exp-021. */
    exin = t_est;
    exinc = 0;
    if ((exkey == [samples, w, r]) && (type(extest) == "t_VEC"),
      exin = vector(length(t_est), i, t_est[i] - extest[i]);
      exsc = vecmax(apply(abs, exin));
      if (exsc > 0,
        exdig = max(30, default(realprecision) + ceil(log(exsc)/log(10)) + 40);
        exp0 = default(realprecision);
        default(realprecision, exdig);
        exin = vector(length(exin), i, precision(exin[i], exdig));
        exinc = 1;
      , exinc = -1);
    );
    if (exinc >= 0,
    if (complextaylor,
      om = exp(2*Pi*I/samples);
      c0 = exp(-Pi*I*(1+1/samples));
      if (2^valuation(samples,2)==samples,
        G = fft(geopow(om^(-1), samples-1), exin)   /* exp-061 */
      ,
        G = mixdft(exin));
      /* exp-062: the key must carry the working precision. The extraction
         runs at the reduced exdig, which CLIMBS while the precision ladder
         does, so a map built early in a grid stretch was reused later at a
         precision it does not carry (measured at dps 300: 68 of 223 hits,
         worst shortfall 58 digits). Reuse only downwards, as bluedft and
         mixdft already do. */
      if ((exmapkey != [samples, w, r, 1]) || (type(exmap) != "t_VEC")
          || (exmapprec < default(realprecision)),
        exmap = vector(terms, s, (rinv^s/samples) * conj(c0)^s * om^(-s));
        exmapkey = [samples, w, r, 1];
        exmapprec = default(realprecision);
      );
      coeffs = vector(terms, s, exmap[s] * G[(s%samples)+1]);
    ,
      mu = exp(Pi*I/samples);
      c1 = exp(-Pi*I/(2*samples));
      if (2^valuation(samples,2)==samples,
        G = fft(geopow(mu^(-1), 2*samples-1), concat(exin, vector(samples, i, 0)))   /* exp-061 */
      ,
        G = mixdft(concat(exin, vector(samples, i, 0))));
      /* exp-062: same for the real branch (this is the one base e takes). */
      if ((exmapkey != [samples, w, r, 2]) || (type(exmap) != "t_VEC")
          || (exmapprec < default(realprecision)),
        exmap = vector(terms, s, conj(c1)^s * mu^(-s));
        exmapkey = [samples, w, r, 2];
        exmapprec = default(realprecision);
      );
      coeffs = vector(terms, s, (rinv^s/samples) * real(exmap[s] * G[(s%(2*samples))+1]));
    );
    );
    if (exinc == 1,
      default(realprecision, exp0);
      coeffs = vector(terms, s, excoef[s] + coeffs[s]);
    );
    if (exinc == -1, coeffs = excoef);
    exkey = [samples, w, r];
    extest = t_est;
    excoef = coeffs;
    wtaylor = Polrev(concat([0], coeffs));
    for (s=1,terms,
      y2=y1;
      y1=y0;
      y0=abs(coeffs[s])*circr^s;
      if ((s>40) && (st==0) && ((y0*0.99)>y1) && ((y0*0.99)>y2),
        st=s;
      );
    );
  ,
    for (s=1,terms,
      tot=0;
      for (t=1,samples,
        t_est[t]=t_est[t]*conj(tcrc[t]);
        if (complextaylor, tot=tot+t_est[t], tot=tot+real(t_est[t]) );
      );
      tot=tot/samples;
      tot=tot*(rinv)^s;
      wtaylor=wtaylor+tot*x^s;
      y2=y1;
      y1=y0;
      y0=abs(tot)*circr^s;
      if ((s>40) && (st==0) && ((y0*0.99)>y1) && ((y0*0.99)>y2),
        st=s;
      );
    );
  );
  if (st==0, st=terms);
  wtaylor=precision(wtaylor,precis);
  if (complextaylor, wtaylor=subst(wtaylor,x,x*exp(-argc*I)));
  return([wtaylor,st]);
}

loop(kc,nlim,nskip,looplim) = {
  local(n,m,me,relast,mt,re,rr,r,z,stopterms,skipdec,ctsamples,thsamples,log10,rlog10);
  if (imag(kc)<0,
    print("bases with imag(k)<0 not supported; use loop(conj(k)); conj(abel(conj(z)));");
    if (x2mode==1, return(0));
  );
  log10=log(10);
  rlog10=1/log10;
  /* exp-014: slow convergence (~2 digits/iter) affects a NEIGHBORHOOD of
     base e (e.g. b=3, kc~1.094, capped at 62.4 true digits), not just
     kc==1. Widened window still excludes the fragile fast bases
     (b=2: kc~0.63, b=10: kc~1.83). */
  efam = (abs(kc-1) < 0.12);
  subeta = ((imag(kc)==0) && (real(kc) <= 0));
  /* exp-071b: sub-eta verlaesst die Kneser-Maschinerie ganz. Vorher lief die
     volle Kontur-Iteration und lieferte still ~15 Stellen; jetzt uebernimmt
     die regulaere Iteration am reellen anziehenden Fixpunkt. Frueher Ausstieg,
     damit kein falscher ct/theta-Zustand zurueckbleibt, den abel() spaeter
     stillschweigend benutzt. kc = 1+log(log b)  =>  b = exp(exp(kc-1)). */
  if (subeta,
    error("loop(): sub-eta (kc=", kc, ") ueber sexpinit(b) aufrufen, nicht ueber loop(kc)."
          " Der Umweg b = exp(exp(kc-1)) verliert genau die Stellen, die die"
          " regulaere Iteration braucht -- die Basis muss exakt bleiben.");
  );
  initsch(kc);
  /* exp-023: for the base-e family a smaller sampling radius trades a
     slightly slower rate (~2.05 -> ~1.7 digits/iter) for a much smaller
     grid (terms/digit ~9 -> ~6): product optimum near ctr*0.9 at dps
     150-300, ~16% faster with unchanged true digits. Fast bases keep
     their original radius (delicate sample/terms co-evolution). */
  /* exp-045: re-scan of the exp-023 radius after ~10 keeps changed the
     cost balance (exp-024/026/027/030-034/037/041b/042b). The N^2 sampling
     term now dominates enough that a smaller radius wins despite the slower
     rate (grid shrinks quadratically, iteration count only grows linearly).
     Scanned with the ctrmul knob on top of the exp-023 value: an extra 0.9
     gives 1.053x at dps 200, 1.119x at dps 300, 1.194x at dps 400 with
     unchanged true digits -- the gain grows with depth. 0.8 was faster at
     dps 300 but LOSES at dps 400 (grid quantization lands it on the same
     2304 grid as 0.9 while needing 382 instead of 299 iterations), so 0.9
     is the monotone choice. 9/10 -> 9/10 * 9/10 = 81/100. */
  if (efam, ctr = ctr*81/100);
  if (ctrmul != 1, ctr = ctr*ctrmul);  /* scan knob, default 1 */
  if (nlim==0,  nlim=70);
  if (nskip==0, nskip=6);
  if (looplim==0, if (limitp==0, looplim=precis-throwp, looplim=limitp));
  /* exp-004: looplim is a digits GOAL; callers passing dps-20 silently lose
     ~17 digits on slow-converging bases (e.g. b=e). Unless limitp requests a
     speed cap, always aim for full working precision. */
  if ((limitp==0) && (looplim < precis-throwp-2), looplim = precis-throwp-2);
  /* exp-007b: see loop body — the iteration-cap raise is rate-adaptive
     (only slow-converging bases get more iterations; fast bases degrade
     when iterated past their plateau, so they keep the caller's cap). */
  skipdec=precis-throws; /* start decrementing nskip here */
  if (thetamode, r=circr*ctr, r=circr);
  ct=0;
  re=-3;
  me=-3;
  relast=-4;
  stopterms=0;
  thlog = (7/2000);
  tht_re_mult=log10/log(thlog);
  n=0;
  rr=renormr(ct);

  while ((re<looplim) && (n<nlim) && (nskip>=0) && ((re>relast) || (nskip>0)),
    n++;
    /* exp-018: extend the iteration cap for ALL bases until the looplim
       goal is reached. exp-007d had restricted this to the e-family
       because extended runs seemed to "degrade" fast bases — that was a
       measurement illusion: they exposed the ~41-digit true error of the
       old nlim-30 references (error-vector test, journal 2026-07-10;
       sexp error == contour error, factor 1, for all bases). */
    if ((n>=nlim-1) && (limitp==0) && (re<looplim),
      nlim = nlim+20);
    /* exp-008: early iterations carry limited signal, so run them at
       reduced working precision. The state's TRUE precision is up to
       2x the contour residual re (factor-2 bases like b=2), so the
       guard must be 2*re + margin — the re+40 guard of exp-006 cut
       real digits of the fast bases. Restored to full before renormslog. */
    default(realprecision, max(48, min(precis, 2*floor(re) + 60)));
    ct=precision(ct,precis);
    ct=ct+rr;
    if (thetamode,
      thetamode++;  /* used by renormr */
      /* exp-063: tht_re_mult starts NEGATIVE (log(10)/log(7/2000)) and only
         turns positive once a nonzero theta coefficient is seen. If that
         coefficient were ever exactly zero the seed would survive, and past
         re ~ 14.7 this expression goes negative straight into vector(). 6 is
         a length the code already uses on purpose (thtaylor(1,6)). */
      thsamples=max(6, floor(re*tht_re_mult+6));
      tht=thtaylor(1,thsamples);
      if (complextaylor, tht2=thtaylor(2,thsamples));  /* might need thsamples2 ... */
      if (n<=2, m=1, m=floor(thsamples/2));
      thlog=abs(polcoeff(tht,m));
      if (thlog<>0,
        thlog=log(thlog)+0.69;           /* 2x = +0.69 buffer for tht_re_mult */
        tht_re_mult=-m*1.05*log10/thlog; /* this is the multiplier for re; 1.05x is a buffer */
      );
    );

    ctsamples = floor((stopterms+20)*1.03);
    ctsamples = 4*floor(ctsamples/4);
    if ((re<relast) || (re>=skipdec),
      nskip--;
      if (re<skipdec,
        if (ctsamples<length(ct), ctsamples=length(ct)-1+4);
      );
    );
    if (ctsamples<length(ct), ctsamples=length(ct)-1);
    z=staylor(circc,r,ctsamples); ct=z[1]; stopterms=z[2];
    rr=renormr(ct);
    relast=re;
    re=abs(polcoeff(rr,1));
    m=length(ct)-1;
    me=abs(polcoeff(ct,m)*r^m);
    if (me>re, re=me);
    re = -log(re)*rlog10;
    re=precision(re,9);
    /* exp-010: near the plateau the loop burns full-cost iterations on
       gains of <0.001 digits before nskip runs out (~3 x 10.4s at e|200).
       Marginal gains spend a stall credit immediately. */
    if ((n>3) && ((re-relast) < 0.1), nskip--);
    if (quietmode==0,
/*    printf("%4d =loopcnt %5.1f decimal digits, %4d ctsamples,%4d %4d thsamples\n", n, re, ctsamples, stopterms, thsamples);*/
      print(n "=loopcnt "re" decimal digits, "ctsamples" ctsamples, "thsamples" thsamples");
    );
  );
  default(realprecision, precis);  /* exp-008: restore full precision */
  ctd = deriv(ct);
  if ((quietmode<0) && ((n%4)<>0) ,print());
  rslog=renormslog(ct);
  if (complextaylor==0, rslog=real(rslog));
/*if (x2mode==1, abelmp5=abel(fs(-0.5))-1);*/
  re = precision(re,9);
  return([re,n,(length(ct)-1)]);
}

loop1(n,ctsamples,thsamples) = {
  local(i,m,me,mr,mt,z,re,rr,r,ct1,rlog10,stopterms);
  rlog10=1/log(10);
  if (n==0, n=1);
  if (ctsamples==0, ctsamples = length(ct)-1);
  ctsamples = 4*floor(ctsamples/4);
  if (thsamples==0, thsamples=length(tht));
  if (quietmode==0, print(ctsamples " ctsamples, "thsamples " thsamples"));
  for (i=1,n,
    ct=precision(ct,precis);
    ct1=polcoeff(ct,1);
    if (thetamode, r=circr*ctr, r=circr);
    if (thetamode,
      thetamode++;  /* used by renormr */
      m = length(tht)-1;
      tht=thtaylor(1,thsamples);
      if (complextaylor, tht2=thtaylor(2,thsamples));  /* might need thsamples2 ... */
    );
    z=staylor(circc,r,ctsamples); ct=z[1]; stopterms=z[2];
    rr=renormr(ct);
    re=abs(polcoeff(rr,1));
    m=length(ct)-1;
    me=abs(polcoeff(ct,m)*r^m);
    mr=abs(polcoeff(ct,1)-ct1);
    if (me>re, re=me);
    re = -log(re)*rlog10;
    re = precision(re,9);
    if (quietmode==0,
      if (mr<>0,
        mr=-log(mr)*rlog10;
        mr = precision(mr,9);
        print(i " =loopcnt " re " decimal digits, "mr " accord");
      ,
        print(i " =loopcnt " re " decimal digits, equations in perfect accord");
      );
    );
  );
  rslog=renormslog(ct);
  re = precision(re,9);
  return([re,n,(length(ct)-1)]);
}

sexpinit(b,nlim,nskip,looplim) = {
  local(z);
  /* exp-071b: reelle Basen 1 < b <= eta verlassen die Kneser-Maschinerie hier,
     BEVOR b ueber kc = 1+log(log b) laeuft. Der Umweg wuerde die Basis auf die
     gerade gesetzte Praezision runden, und die regulaere Iteration braucht sie
     auf etwa der doppelten Zielstellenzahl. */
  if ((imag(b)==0) && (real(b) > 1) && (real(b) <= etaB),
    /* exp-071c: subeta ERST nach erfolgreichem reginit setzen -- sonst bleibt
       nach einem gescheiterten sexpinit eine intakt geladene Kneser-Basis
       unbenutzbar, weil die Waechter greifen. */
    z = reginit(b);
    subeta = 1;
    return(z);
  );
  z=loop(log(log(b))+1,nlim,nskip,looplim);
  if (quietmode==0, print("sexp(z); slog(z); sexptaylor(0); /* sexptaylor series at 0; */"));
  return(z);
}


f(z) = {
  local(y);
  if (fmode<2,
    if (fmode==0, y=abel(z), y=invabel(z));
  ,
    if (x2mode==0,
      if (fmode==2, y=slog(z), y=sexp(z));
    ,
      if (fmode==2, if (real(z)>-0.5, y=abel(z), y=abel(-1-z)), y=invabel(z)-L;);
    );
  );
  return(y);
}

/* Complex function magnitude/phase plotter, originally from Mike */
/* To use:
*     1. Define function to graph as func(z).
*     2. Load this program.
*     3. Execute MakeGraph(width, height, x0, y0, x1, y1, filename) with the parameters given as follows:
*        width, height = width/height of image in pixels
*        x0, y0, x1, y1 = rectangle of complex plane to graph: x0 + y0i in upper-left corner to x1 + y1i in lower-right corner
*        filename = name of file to save as.
* Output is in .PPM format.
*/
/* Color conversion (HSB to RGB). */

HSB2RGB(z) = {

  local(HH,F,P,Q,T);
  local(mag,phase,H,S,B);

  mag = abs(z);
  if (z<>0, phase=arg(z), phase=0);
  H = phase/(2*Pi);
  S = 1/(1 + 0.3*log(mag + 1));
  B = 1 - 1/(1.1 + 5*log(mag + 1));

  HH = floor(6*H)%6;
  F = (6*H) - floor(6*H);
  P = B*(1 - S);
  Q = B*(1 - (F*S));
  T = B*(1 - (1-F)*S);
  if(HH == 0, return([B, T, P]));
  if(HH == 1, return([Q, B, P]));
  if(HH == 2, return([P, B, T]));
  if(HH == 3, return([P, Q, B]));
  if(HH == 4, return([T, P, B]));
  if(HH == 5, return([B, P, Q]));
}



/* This version of MakeGraph, uses brackets and commas as white space; /* instead of just spaces as deliminators.  Use "MakeGraphLegal=1; for slower legally correct version */
writeppm(filename,z) = {
/* this version is legal, writing 16 pixels at a time, but is somewhat slower than writing the entire line with brackets and commas */
  local(i,n);
  n=length(z);
  i=1;
  while (n>48, write1(filename,z[i+00]" "z[i+01]" "z[i+02]" "z[i+03]" "z[i+04]" "z[i+05]" "z[i+06]" "z[i+07]" "z[i+08]" "z[i+09]" "z[i+10]" "z[i+11]" "z[i+12]" "z[i+13]" "z[i+14]" "z[i+15]" "z[i+16]" "z[i+17]" "z[i+18]" "z[i+19]" "z[i+20]" "z[i+21]" "z[i+22]" "z[i+23]" "z[i+24]" "z[i+25]" "z[i+26]" "z[i+27]" "z[i+28]" "z[i+29]" "z[i+30]" "z[i+31]" "z[i+32]" "z[i+33]" "z[i+34]" "z[i+35]" "z[i+36]" "z[i+37]" "z[i+38]" "z[i+39]" "z[i+40]" "z[i+41]" "z[i+42]" "z[i+43]" "z[i+44]" "z[i+45]" "z[i+46]" "z[i+47]" ");n=n-48;i=i+48);
  if    (n>24, write1(filename,z[i+00]" "z[i+01]" "z[i+02]" "z[i+03]" "z[i+04]" "z[i+05]" "z[i+06]" "z[i+07]" "z[i+08]" "z[i+09]" "z[i+10]" "z[i+11]" "z[i+12]" "z[i+13]" "z[i+14]" "z[i+15]" "z[i+16]" "z[i+17]" "z[i+18]" "z[i+19]" "z[i+20]" "z[i+21]" "z[i+22]" "z[i+23]" ");n=n-24;i=i+24);
  if    (n>12, write1(filename,z[i+00]" "z[i+01]" "z[i+02]" "z[i+03]" "z[i+04]" "z[i+05]" "z[i+06]" "z[i+07]" "z[i+08]" "z[i+09]" "z[i+10]" "z[i+11]" ");n=n-12;i=i+12);
  if   (n==3,  write (filename,z[i+00]" "z[i+01]" "z[i+02])
  , if (n==6,  write (filename,z[i+00]" "z[i+01]" "z[i+02]" "z[i+03]" "z[i+04]" "z[i+05])
  , if (n==9,  write (filename,z[i+00]" "z[i+01]" "z[i+02]" "z[i+03]" "z[i+04]" "z[i+05]" "z[i+06]" "z[i+07]" "z[i+08])
  , if (n==12, write (filename,z[i+00]" "z[i+01]" "z[i+02]" "z[i+03]" "z[i+04]" "z[i+05]" "z[i+06]" "z[i+07]" "z[i+08]" "z[i+09]" "z[i+10]" "z[i+11])
  ))));
}
/* This version of MakeGraph is legal, and reasonably fast, using writeppm above */
/* Use "MakeGraphLegal=0 for even faster version that prints a line at a time */
MakeGraph(width, height, x0, y0, x1, y1, filename, n) = {
  local(xt,y,z,n2,Red,Green,Blue,xstep,ystep,xx,yy,lineg,grey);
  lineg = vector(3*width);
  xstep = (x1 - x0)/width;
  ystep = (y1 - y0)/height;
  n1 = abs(n);
  n2 = n1/2.0;
  write(filename, "P3");
  write(filename, "# fatou.gp k= "k " coordinates "x0+y0*I" "x1+y1*I);
  write(filename, width, " ", height);
  write(filename, "255");
  for(y=0, height-1,
    for(xt=0, width-1,
      xx = x0+(xstep*xt);
      yy = y0+(ystep*y);
      z = xx+yy*I;
      grey=0;
      /* grey grid lines every "n" spaces, highlighting the (0,0) origin */
      if (n1 <> 0,
        if (abs(((xx+n2)%n1)-n2)<1E-6, grey=1);
        if (abs(((yy+n2)%n1)-n2)<1E-6,
          if (xx>0, grey=1);
          if (yy>0.5, grey=1);
          if (yy<-0.5, grey=1);
          if (n<0, grey=1);
        );
      );
      if (grey,
        RGB=[128,128,128];
      ,
        RGB = floor(HSB2RGB(f(z))*255.9999);
      );
      lineg[xt*3+1]=RGB[1];
      lineg[xt*3+2]=RGB[2];
      lineg[xt*3+3]=RGB[3];
    );
    /* writeppm writes 16 pixels at a time; but is legal; MakeGraphLegal=0 writes a line at once */
    if (MakeGraphLegal, writeppm(filename,lineg), write(filename, lineg));
    print(y);
  );
}

halfsexp(z)={
  local(y,ar,m);
  subetaguard("halfsexp");
  m=0;
  z = invabel_sexp(z);
  if (abs(z-circc)>circr,
    ar = arg(z-circc)-argc;
    if ((ar<=(Pi/2)) && (ar>=-Pi/2),
      z=finv(z);
      m=1;
    ,
      z=fs(z);
      m=-1;
    );
  );
  y = invabel(abel(z)+0.5);
  if (m==1,  y=fs(y));
  if (m==-1, y=finv(y));
  y = sexp_invabel(y);
  return(y);
}

testz(z) = {
  local(y,y1,y2,ar);
  ar=arg(z-circc)-argc;
  y1=finv(z);
  y2=fs(z);
  if (((abs(ar)<Pi/2) && (abs(y1-circc)<circr)) || ((abs(y1-circc)<abs(y2-circc)) && (abs(y2-circc)>circr)),
    y = abelest(y1,ct) - abelest(z,ct) + 1;
  ,
    y = abelest(y2,ct) - abelest(z,ct) - 1;
  );
  return(y);
}

/* pentation routines begins here */
/* help_penthex() shows the help menu */
pcent=0;     /* set to ipent(1) to normalize so pent(0)=1 */
pentr=0.05;  /* constant; the pentation schroder/ischroder are very accurate withing 0.05 of the fixed point */
             /* use a smaller value of pentr for \p 57, higher than 32 digit precision */
pentfixr=0.7;

/* exp-076 (2026-10-09): Pentation mit der Praezision des Motors.
   Das Original rechnete die Poincare-Funktion mit FESTEN 21 Gliedern
   (formalischroder(tslog,21)) und lieferte damit ~30 Stellen, egal wie viele
   verlangt waren -- bei \p 60: pent(0)-1 = 6.5e-31. Jetzt waechst die
   Gliederzahl mit precis, und der Weg laeuft ueber slog() statt slogtaylor(),
   damit er auch fuer reelle 1 < b < eta (regulaere Tetration, dort ist
   slogtaylor gesperrt) und fuer komplexe Basen gilt.
   Konstruktion unveraendert: Fixpunkt L von sexp auf (-2,-1), dort ist slog
   anziehend mit mu = 1/lambdat; Psi mit Psi(mu*w) = slog(Psi(w)) ist die
   Poincare-Funktion, pent(z) = Psi(lambdat^(z+pcent)) mit pent(0) = 1. */

/* Cauchy-Taylorreihe von f um w auf |z-w| = r mit `samples` Punkten, f eine
   Closure (slog, z->pent(z), ...). Generisch: slog() leitet fuer sub-eta an
   regslog weiter. Grad samples-1; der Koeffizient k traegt den absoluten
   Fehler ~ eps/r^k, was nicht stoert: ausgewertet wird nur weit innerhalb
   von r (|Psi - L| <= pentr << r). */
pentcauchy(fn, w, r, samples) = {
  my(tc, v, p, c);
  tc = vector(samples, t, exp(2*Pi*I*(t-1)/samples));
  v = vector(samples, t, fn(w + r*tc[t]));
  p = 0;
  for (s = 0, samples-1,
    c = vecsum(v)/samples/r^s;
    p = p + c*x^s;
    v = vector(samples, t, v[t]*conj(tc[t]));
  );
  /* exp-096: die hohen Glieder kommen nach der Kuerzung mit weniger Woertern
     zurueck; im Horner-Schema reichte ein 19-stelliges Glied die kurze
     Genauigkeit bis zum Ergebnis durch (b = 1.5, Rang 7, dps 30: Fehler
     4e-20 statt 1e-32). Auf precis setzen. */
  return(precision(p, precis));
}

/* Taylorreihe von f um einen Fixpunkt nahe t0 (|t0 - L| << sr). Die
   Stuetzstellenzahl folgt aus (sr/R)^samples < eps, R = Konvergenzradius,
   erst auf 2 geschaetzt, dann am Wurzeltest der Reihe selbst nachgemessen.
   Fixpunkt per Newton auf der Reihe, Reihe darauf zentriert.
   Rueckgabe [L, sc, R] mit sc(t) = f(L+t) - L, sc(0) = 0. */
pentfixser(fn, t0, sr, isreal, eps, name) = {
  my(R, samples, S, s, t, lx, L, sc);
  R = 2;
  for (pass = 1, 3,
    samples = ceil((precis+25)*log(10)/log(R/sr));
    S = pentcauchy(fn, t0, sr, samples);
    if (isreal, S=real(S));
    /* exp-096: die winzigen hohen Glieder kommen aus der Kuerzung mit weniger
       Woertern zurueck (absolut genau, relativ kurz); pentpsi multipliziert sie
       weiter, und der Verlust wuchs bis "Psi-Glied 32 traegt nur 19 Stellen"
       (b = 1.5, Rang 7, dps 45). Wie ct in loop(): auf precis setzen. */
    S = precision(S, precis);
    s = max(8, samples\4);
    R = min(abs(polcoeff(S, s))^(-1/s), abs(polcoeff(S, s\2))^(-2/s));
    if (R <= 1.25*sr, error(name, ": Reihe hat um den Fixpunkt nur Radius ", precision(R,5), " -- sr kleiner waehlen"));
    if ((sr/R)^samples < eps, break);
  );
  if ((sr/R)^samples >= eps, error(name, ": Abtastung erreicht die Praezision nicht (R = ", precision(R,5), ")"));
  t = 0;
  for (i = 1, 200,
    lx = (subst(S, x, t) - t0 - t)/(subst(deriv(S), x, t) - 1);
    t = t - lx;
    if (abs(lx) < 10.0^(-precis-5), break);
  );
  L = t0 + t;
  sc = subst(S, x, x + t) - L;
  sc = sc - polcoeff(sc, 0);
  return([L, sc, R]);
}

/* Poincare-Funktion des zentrierten slog sc (sc(0) = 0, sc'(0) = mu):
   Psi(mu*w) = sc(Psi(w)), Psi = w + O(w^2). Newton-Verdopplung: ist Psi
   mod w^m bekannt, so loest Delta = O(w^m) das Dreieckssystem
     (mu^j - mu) d_j = e_j + sum_{m<=i<j} g_{j-i} d_i,
   E = sc(Psi) - Psi(mu*w), G = sc'(Psi). Kosten O(K*M(N)*log N) statt der
   O(N^4) von formalischroder, das jedes Glied mit voller Komposition neu
   ansetzt. Stopp, wenn der neueste Block, gewichtet mit rr^j, unter eps liegt
   (fruehestens ab Grad 8), sonst bei nmax. Rueckgabe [Psi, ok]. */
pentpsi(sc, rr, eps, nmax) = {
  my(mu, scd, P, m, m2, Pw, E, G, d, acc, tail, ok);
  mu = polcoeff(sc, 1);
  scd = deriv(sc);
  P = x;
  m = 2;
  ok = 0;
  while (m <= nmax,
    m2 = min(2*m, nmax+1);
    Pw = P + O(x^m2);
    E = subst(truncate(sc + O(x^m2)), x, Pw) - subst(P, x, mu*x);
    G = subst(truncate(scd + O(x^m2)), x, Pw);
    d = vector(m2, j, 0);
    tail = 0;
    for (j = m, m2-1,
      acc = polcoeff(E, j);
      for (i = m, j-1, acc = acc + polcoeff(G, j-i)*d[i]);
      d[j] = acc/(mu^j - mu);
      tail = max(tail, abs(d[j])*rr^j);
    );
    P = P + sum(j = m, m2-1, d[j]*x^j);
    m = m2;
    if ((m >= 8) && (tail < eps), ok = 1; break);
  );
  return([P, ok]);
}

/* Reihe fuer die Auswertung bei |w| <= rr herrichten. Die hohen Glieder
   verlieren im Newton-Loeser relative Stellen (Ausloeschung in E, bei Basis
   100 sind es 19 von 57 Stellen bei Grad 135) -- legitim, denn ihr Beitrag
   |c|*rr^n liegt weit unter eps. PARI fuehrt die Praezision im Horner-Schema
   aber RELATIV mit: ein einziges solches Glied drueckte pent(0)-1 bei Basis
   100 auf 3.6e-21. Daher: jedes Glied, dessen absoluter Fehlerbeitrag
   |c|*10^-prec(c)*rr^n unter 10^(1-precis) liegt, auf precis auffuellen
   (wie regsigma);
   liegt er darueber, laut abbrechen. Glieder hinter dem letzten, das noch
   ueber eps*1e-10 beitraegt, entfallen. */
pentclean(P, rr, eps, name) = {
  my(c, n0, Q, e);
  n0 = 0;
  for (n = 0, poldegree(P),
    c = polcoeff(P, n);
    if ((c != 0) && (abs(c)*rr^n > eps*1e-10), n0 = n));
  Q = 0;
  for (n = 0, n0,
    c = polcoeff(P, n);
    if (c != 0,
      e = if ((type(c)=="t_INT") || (type(c)=="t_FRAC"), 0, abs(c)*rr^n*10.0^(-precision(c)));
      if (e > 10.0^(1-precis), error("pentinit: ", name, "-Glied ", n, " traegt nur ", precision(c), " Stellen, Fehlerbeitrag ", precision(e,3)));
      c = precision(c, precis));
    Q = Q + c*x^n;
  );
  return(Q);
}

/* pentinit(B, skip, rr, sr): rr = Auswerteradius der Reihen (Default 0.05
   wie im Original, gedeckelt durch den Konvergenzradius ~ lambdat*|L+2| von
   Psi), sr = Radius der slog-Abtastung (Default 1/2). */
pentinit(B,skip,rr,sr) = {
  my(s,lx,t0,sc,res,eps,nmax,z1,z2,errp,errz,zz);
  errp=0;
  hexreg=0;           /* exp-077: eine neue Basis entwertet jede Hexation */
  if (skip==0,
    zz=sexpinit(B);   /* initialize sexp using fatou.gp program */
    if (type(zz)=="t_VEC", errp=zz[1]);   /* sub-eta liefert nur L, keine Schaetzung */
  );
  if (sr==0, sr=1/2);
  eps = 10.0^(-precis-8);
  /* Fixpunkt von sexp auf (-2,-1): fuer slog anziehend, Rate 1/lambdat.
     Grob iterieren, fein per Newton auf der Reihe. */
  lx=100; t0=-1.97; s=0;
  while ((abs(lx-t0) > 1e-8) && (s < 2000), lx = t0; t0 = slog(t0); s++);
  if (s >= 2000, error("pentinit: slog-Iteration zum Fixpunkt konvergiert nicht"));
  res = pentfixser(slog, t0, sr, imag(B)==0, eps, "pentinit");
  tfixed = res[1];
  sc = res[2];
  tslog = sc + tfixed;                 /* slog um tfixed, wie im Original */
  lambdat = 1/polcoeff(sc, 1);         /* slope at fixed point of sexp */
  if (rr==0, rr = min(0.05, abs(lambdat*(tfixed+2))/8));
  pentr = rr;
  nmax = poldegree(sc) - 1;
  res = pentpsi(sc, rr, eps, nmax);
  if (res[2]==0, error("pentinit: Poincare-Reihe bis Grad ", nmax, " nicht konvergiert -- kleineres rr oder groesseres sr"));
  ipentz = Pol(serreverse(res[1] + O(x^(poldegree(res[1])+1))));
  if (abs(polcoeff(ipentz, poldegree(ipentz)))*rr^poldegree(ipentz) > eps,
    error("pentinit: Umkehrreihe konvergiert bei pentr = ", rr, " nicht ausreichend"));
  pentz  = pentclean(res[1], rr, eps, "Psi") + tfixed;
  ipentz = pentclean(ipentz, rr, eps, "sigma");

  pperiod=2*Pi*I/log(lambdat);
  pcent=0;              /* set to zero before we normalize so pent(0)=1 */
  pcent=ipent(1); /* normalize pentation so pent(0)=1, pcent=pent^{-1}(1) */
  /* Fehlerschaetzung: dieselben Werte ueber zwei Reduktionstiefen (anderer
     Reihenpunkt, andere sexp-Kette), dazu pent(0) = 1 (Psi gegen sigma) und
     ein Rundweg. Das Original verglich mit `z1<1` und lief deshalb nur fuer
     reelle Basen. Der Motor selbst bleibt dabei unsichtbar: seine Stellen
     begrenzen die Pentation zusaetzlich (errp). */
  errz = abs(pent(0)-1);
  errz = max(errz, abs(ipent(pent(1/3))-1/3));
  for (j = 1, 2,
    zz = if (j==1, 1/2, -1/2);
    z1 = pent(zz); z2 = pent(zz, 1);
    errz = max(errz, abs(z1-z2)/abs(z1));
  );
  errz = if (errz==0, precis, min(precis, -log(errz)/log(10)));
  errz = precision(errz,9);
  if ((errp>0) && (errp<errz), errz=errp); /* if pentation precision<tetration precision ... */
  psing=ipent(-2)+1;  /* this is the pentation singularity, where pent(z-1)=-2, pent(z)=singularity */
  if (quietmode==0,
    print (" ");
    print ("pentation base              " sexp(1));
    print ("pentation(-0.5)             " pent(-0.5));
    print ("sexp fixed point    tfixed  " tfixed);
    print ("sexp slope at fixed lambdat " lambdat);
    print ("pentation period    pperiod " pperiod);
    print ("pentation singularity psing ");
    print(" "psing);
    print ("pentation precision:        " errz " decimal digits");
  );
  return(errz);
}
pent(z,extra)={
  /* generate pentation approximation, pentinit first */
  /* generated from from lambdat, the slope at the fixed point, and the pentz ischroder function */
  /* exp-076: extra = zusaetzliche Reduktionsschritte (fuer die Fehlerschaetzung) */
  local(n,i,y);
  y=lambdat^(z+pcent);
  n=0;
  while (abs(y)>pentr, y=y/lambdat;n++); /* iterate n times until y<pentr, pentr=0.05 */
  for (i=1,extra, y=y/lambdat;n++);
  y = subst(pentz,x,y);
  for (i=1,n,y=sexp(y));                 /* iterate sexp n times */
  return(y);
}
ipent(z)={
  /* generate the inverse pentation approximation, pentinit first */
  /* generated from from lambdat, the slope at the fixed point, and the ipentz schroder function */
  local(n,i,y);
  n=0;
  y=z;
  while (abs(y-tfixed)>pentr, y=slog(y);n++; /* iterate slog n times until y is within pentr of tfixed, pentr=0.05 */
    /* exp-081: fuer manche komplexen z erreicht die slog-Kette den Fixpunkt
       nie -- ohne Obergrenze hing knhexinit bei b = e in der Bogenplanung */
    if (n > 60, error("ipent: slog-Kette erreicht den Fixpunkt nicht")));
  y=subst(ipentz,x,(y-tfixed));
  y=log(y)/log(lambdat)-pcent+n;  /* add n to the schroder function approximation for ipent(z) */
  return(y);
}
pentaylor( w,r,samples,inv_pent) = {
  local(rinv,s,t,x1,y,z,tot,t_est,tcrc,halfsamples,wtaylor,terms);
/* outputs gie taylor series, the complex taylor series for sexp, */
/* default halfsamples=240 */
  if (samples==0, samples=240);  /* no matter how many sample points, the default gie series size is 200 halfsamples */
  halfsamples=samples/2;
  terms = floor(samples*200/240);
  t_est    = vector (samples,i,0);
  tcrc     = vector (samples,i,0);
  if (r==0,r=1);
  rinv = 1/r;
  wtaylor=0;
  for(s=1, samples, x1=-1/(samples)+(s/halfsamples); tcrc[s]=exp(Pi*I*x1); );
  for (t=1,samples,
    if (inv_pent==0, t_est[t]=pent(w+r*tcrc[t]), t_est[t]=ipent(w+r*tcrc[t]));
  );
  for (s=0,terms-1,
    tot=0;
    for (t=1,samples,
      tot=tot+t_est[t];
      t_est[t]=t_est[t]*conj(tcrc[t]);
    );
    tot=tot/samples;
    if (s>=1, tot=tot*(rinv)^s);
    wtaylor=wtaylor+tot*x^s;
  );
  wtaylor=precision(wtaylor,precis);
  return(wtaylor);
}
sexpupfixed(limitlo,limithi) = {
  local (s,z,duplim,dlolim,limittst,evalsexp,evallo,evalhi,pencnt,firstlo,firsthi,plim,tsexp,dsexp);
  quietmode=1;
  if ((limitlo==0) || (limithi==0), limitlo=1.5; limithi=2);
  evallo=-1.0;
  evalhi=1.0;
  evalsexp=1.0;
  firstlo=1;
  firsthi=1;
  pencnt=0;
  precis=precision(evalsexp);
  plim=10^(-precis+8);
  while (abs(evalsexp)>plim,
    ratio=-evallo/(evalhi-evallo);
    if (ratio<0.1,ratio=ratio*2);
    if (ratio>0.9,ratio=ratio*2-1;);
    limittst=limitlo+ratio*(limithi-limitlo);
    limittst=precision(limittst,precis); /* make sure base has full precision */
    sexpinit(limittst);
    tsexp = real(sexptaylor(3,1));
    /* derivative of sexp(x)-x */
    dsexp = precision(deriv(tsexp-x),precis);
    dlolim=2.0;
    duplim=4.0;
    s=1;
    while ((duplim-dlolim)>plim && (s<200),
      upfixed=(duplim+dlolim)/2;
      z = upfixed-3;
      if (subst(dsexp,x,z)>0, duplim=upfixed, dlolim=upfixed);
      s++;
    );
    evalsexp=sexp(upfixed)-upfixed;
    print ("base " sexp(1)" "evalsexp);
    if (evalsexp>0, limithi=limittst;evalhi=evalsexp;firsthi=0,
                    limitlo=limittst;evallo=evalsexp;firstlo=0);
    /* first loop initializations */
    if (firsthi, evalhi=-evalsexp); /* binary searching until evaluated hi */
    if (firstlo, evallo=-evalsexp); /* binary searching until evaluated lo */
    pencnt++;
  );
  pentinit(limittst);
  tsexp = real(sexptaylor(upfixed,1));
  print (" ");
  print ("sexp base, sexp(upfixed)=upfixed ", real(sexp(1)));
  print ("sexpinit(B) iterations required  ", pencnt);
  print ("upfixed, parabolic fixed point   ", upfixed);
  print ("sexp'(upfixed) base B            ", real(polcoeff(tsexp,1)));
  print ("sexp(upfixed)-upfixed error      ", real(polcoeff(tsexp,0))-upfixed);
  quietmode=0;
  upfixed;
}



/* abel for hexation approximation to pentation simulaneous equations without theta */
hex_abel(z,e) = {
  z = (1/log(lambda_pf))*(log(I*(z-pentf))-Pi*I/2) + conj(1/log(lambda_pf))*(log(-I*(z-conj(pentf)))+Pi*I/2) + subst(e,x,(z-real(pentf)));
  if (imag(z)==0, z=real(z));
  return(z);
}
/* matrix solution to hexation simulaneous equations without theta */
matrix_h(lctr) = {
  local(ar,y0,yi,ys,z0,zi,zs,z,m2,s2,r2,i,j,nt,h,lf,r,tot);
  i=0; ;pentf=0.5; lf=2;
  while ((lf<>pentf) && (i<1000), lf=pentf;pentf=ipent(pentf);i++;);
  r = abs(pentf-psing)*0.7;
  pt = pentaylor(pentf,r,0,1);
  lambda_pf = 1/polcoeff(pt,1);
  print ("pentf ");
  print (pentf);
  print ("lambda_pf ");
  print (lambda_pf);
  hex_period=2*Pi*I/log(lambda_pf);
  print ("hex pseudo period");
  print (hex_period);
  /* a*zi + j*zi^2 - a*z0 + j*z0^2 = 1 - (yi-y0);  */
  /* a*(zi-z0) + j*(zi^2-z0^2) = 1 - (yi-y0)       */
  if ((lctr%4), lctr=lctr+4-(lctr%4));
  m2=matrix(lctr,lctr);
  s2=matrix(lctr,1);
  h = lctr/2;
  nt=h;
  for (i=1,nt,
    z0 = real(pentf)+imag(pentf)*exp((i-0.5)*2*Pi*I/lctr);
    ar = arg(z0-real(pentf));
    y0 = hex_abel(z0,0);
    zi = ipent(z0);
    zs = pent(z0);
    if (abs(ar)<Pi/2,
      yi = hex_abel(zi,0);
      for (j=1,lctr, m2[i,j] = (zi-real(pentf))^j-(z0-real(pentf))^j);
      s2[i,1] = -1 - (yi-y0);
    ,
      ys = hex_abel(zs,0);
      for (j=1,lctr, m2[i,j] = (zs-real(pentf))^j-(z0-real(pentf))^j);
      s2[i,1] = 1 - (ys-y0);
    );
  );
  for (i=1,h,
    s2[i+h,1]=imag(s2[i,1]);
    s2[i,1]  =real(s2[i,1]);
    for (j=1,lctr,
      m2[i+h,j]=imag(m2[i,j]);
      m2[i,j]  =real(m2[i,j]);
    );
  );
  z = real(z);
  r2=matsolve(m2,s2);
  for (i=1,lctr, z=z+r2[i,1]*x^i);
  pt = z;
  ptd = deriv(pt);
  hex_bias=hex_abel(-2,pt)+3;   /* hex_bias so ihex(1)=0 */
  /* calculate decimal precision */
  tot = 0;
  z0 = real(pentf)+imag(pentf);
  zi = ipent(z0);
  tot = tot + abs(hex_abel(z0,pt)-hex_abel(zi,pt)-1);
  z0 = real(pentf)-imag(pentf);
  zs = pent(z0);
  tot = tot + abs(hex_abel(z0,pt)-hex_abel(zs,pt)+1);
  tot = -log(tot/2)/log(10);
  tot = precision(tot,9);
  print ("hex(z); ihex(z); ihex_deriv(z);"); /* returns ihex normalized so ihex(1)=0 */
  print ("hex precision "tot " decimal digits");
  return(tot);
}
ihex(z) = {
  local(n,y,zn,ar);
  if (hexreg==1, return(regihex(z)));   /* exp-077 */
  if (hexreg==2, return(knihex(z)));    /* exp-079 */
  return(ihex_deriv(z)[1]);
}
ihex_deriv(z) = {
  local(n,y,yd,zn,ar);
  /* returns hex_log(z) and its derivative */
  n=0;
  /* decide whether to iterate pent(z) or ipent(z) */
  /* before calculating hex_log, based on arg(z) */
  ar = z-real(pentf);
  if ((ar<>0) && (abs(arg(ar))<=Pi/2),
    while (abs(z-real(pentf))>imag(pentf),n++;z=ipent(z));
  ,
    while (abs(z-real(pentf))>imag(pentf),n--;z=pent(z));
  );
  y = hex_abel(z,pt)+n-hex_bias; /* hex_bias sets log(1)=0 */
  yd = (1/log(lambda_pf))/(z-pentf) +
       conj(1/log(lambda_pf))/(z-conj(pentf)) +
       subst(ptd,x,(z-real(pentf)));
  if (imag(z)==0, yd=real(yd));
  return([y,yd]);
}
/* a pretty good version of hexation function; inverse of ihex(z) */
hex(z) = {
  local(n,y,i,tc,est,iadd);
  if (hexreg==1, return(reghex(z)));    /* exp-077 */
  if (hexreg==2, return(knhex(z)));     /* exp-079 */
  /* center on -3.5 .. -2.5 */
  if (imag(z)<0,
    tc = tc - imag(hex_period)*imag(z)/real(hex_period);
  ,
    tc = tc + imag(hex_period)*imag(z)/real(hex_period);
  );
  n = -floor(real(z)+3.5+tc);
  z = z+n;
  if ((imag(z)>0.2),
    /* abel(z) = log(z-lamda_pf)*(1/log(lambda_pf)) + log(z-conj(pentf))*conj(1/log(lambda_pf)) + subst(pt,x,z-real(lamda_pf)) + hex_bias  */
    /* inverse of abel as z approaches lambda approaches this ...  which is a good enough initial approximation */
    y = z - log(2*imag(pentf))*conj(1/log(lambda_pf)) - subst(pt,x,I*imag(pentf)) - hex_bias;
    est = exp(y/(1/log(lambda_pf)))+pentf;
    est = est;
  , if ((imag(z)<-0.2),
    y = z - log(-2*imag(pentf))*(1/log(lambda_pf)) - subst(pt,x,-I*imag(pentf)) - hex_bias;
    est = exp(y/conj(1/log(lambda_pf)))+conj(pentf);
    est = est;
  ,
    est= -2.5;  /* this is the initial estimate for |imag(z)|<0.2 */
  ));
  i = 0;
  y = [0,1];
  while ((i<13) && (abs(z-y[1])>1E-15),
    i++;
    y = ihex_deriv(est);
    est = est + (z-y[1])/y[2]
  );
  if (i>12, print("hex didn't converge");print(i" "z); print(n" "est););
  if (n>0, for (i=1,n,est=ipent(est)));
  if (n<0, for (i=1,-n,est=pent(est)));
  return(est);
}
hexinit(i,n,skip)={
  my(e);
  if (n==0, n=256);
  e = pentinit(i,skip);
  return(hexsetup(e));
}
/* Hexation auf der schon initialisierten Pentation; e = deren Stellenzahl */
hexsetup(e)={
  /* exp-077: hat die Pentation einen anziehenden reellen Fixpunkt (reelle
     1 < b < eta5 = 1.7373559...), ist die Hexation dort regulaer und traegt
     die volle Praezision. Sonst Sheldons Naeherung am komplexen Fixpunkt. */
  if (reghexinit(), return(min(e, hexerr)));
  /* exp-079: sonst Kneser am komplexen Fixpunkt der Pentation; Sheldons
     matrix_h (~8 Stellen) nur noch, wenn das scheitert */
  if (iferr(knhexinit(), E, if (quietmode==0, print("knhexinit: ", E)); 0), return(min(e, hexerr)));
  if (quietmode==0, print());
  return(min(e, matrix_h(256)));
}

/* exp-080 (2026-10-09): Einstiege fuer die Python-Schnittstelle. Sie
   initialisieren Pentation bzw. Hexation nur, wenn Basis oder Praezision
   gewechselt haben -- so kann jeder Ausdruck pentx(z) fuer sich stehen, auch
   wenn ein Pool die Ausdruecke auf mehrere Worker verteilt. */
pentstate=[0,0]; hexstate=[0,0]; pentdigits=0; hexdigits=0;
pentensure()={
  if ((pentstate[1] != sexp(1)) || (pentstate[2] != precis),
    pentdigits = pentinit(sexp(1), 1);
    pentstate = [sexp(1), precis];
    hexstate = [0, 0]);
  return(pentdigits);
}
hexensure()={
  pentensure();
  if ((hexstate[1] != sexp(1)) || (hexstate[2] != precis),
    hyatclear(6);   /* exp-097: Atlas ab Rang 6 gehoert zur alten Hexation */
    hexdigits = hexsetup(pentdigits);
    hexstate = [sexp(1), precis]);
  return(hexdigits);
}
/* exp-083: jenseits des Ueberlaufs liefert die Kette manchmal eine Null ohne
   gueltige Stellen (b = 100: pent(1.96875) = 0.E495855885102019914, d. h. nur
   "Betrag unter 10^(5e20)"); Python las das als 0. Dann Fehler statt Wert. */
levelchk(y, name)={
  if (y == 0 && exponent(y) > 0, error(name, ": Ueberlauf, keine gueltigen Stellen"));
  return(y);
}
pentx(z)={ pentensure(); return(levelchk(pent(z), "pentx")); }
ipentx(z)={ pentensure(); return(levelchk(ipent(z), "ipentx")); }
hexx(z)={ hexensure(); return(levelchk(hex(z), "hexx")); }
ihexx(z)={ hexensure(); return(levelchk(ihex(z), "ihexx")); }

/* exp-088 (2026-10-13): Hyperoperationen beliebigen ganzzahligen Rangs r.
   H_0(b,x) = x+1, H_1 = b+x, H_2 = b*x, H_3 = b^x, H_4 = sexp, H_5 = pent,
   H_6 = hex; ab r = 7 dieselbe Rekursion H_r(b,0) = 1,
   H_r(b,x+1) = H_{r-1}(b, H_r(b,x)), und zwar
   - ungerade r: regulaere Iteration am linken reellen Fixpunkt von H_{r-1},
     der fuer H_{r-1} abstossend ist -- wie die Pentation am Fixpunkt von
     sexp (b = 2, r = 7: p = -3.6994, Multiplikator 9.49),
   - gerade r: regulaere Iteration am reellen anziehenden Fixpunkt von
     H_{r-1}, solange es ihn gibt (b < eta_{r-1}: eta_6 = 1.7991,
     eta_7 = 1.8398, eta_8 = 1.8681); darueber braucht es Kneser wie bei der
     Hexation (noch nicht fuer r >= 8).
   Bauweise nach dem Wissensgraphen (forschung/eta_hyperoperationen/hyper.gp,
   12.10.2026), hier mit Radiensuche, Reihenpflege (pentclean) und
   Fehlerschaetzung wie bei pentinit/reghexinit. Zustand je Rang in Vektoren,
   damit sich die Raenge beim Verschachteln nicht ueberschreiben. */
/* exp-097: stabile Auswertung eines Polynoms an komplexer Stelle. PARIs subst
   nimmt fuer reelle Polynome an komplexer Stelle einen Weg ueber das reelle
   Quadrat (Goertzel-artig), der nahe |u| = 1 Stellen verliert: Grad 127,
   dps 30: |u| = 0.9 Fehler 6e-37, |u| = 0.99 1e-27 statt 1e-38. Das
   Skalarprodukt mit den Potenzen von u (geopow, Verdopplung wie bei den
   FFT-Drehfaktoren) ist exakt. Nur im Atlas-Code (er wertet genau auf |u| = 1
   aus); die Kneser-Auswertung behaelt subst -- mit cpe dort lag H_7 bei b = 2,
   dps 16 um 0.8 Stellen schlechter (H_6 selbst etwas besser). */
/* dieselben Verdopplungsprodukte wie geopow (v[blk+1]*v[j+1], gleich bis aufs
   letzte Bit der Rundung), aber je Verdopplung EIN Vektorbefehl statt einer
   Schleife im Interpreter: n = 255 doppelt so schnell. cpe wird im Atlas
   hunderttausendfach gerufen. geopow selbst bleibt, der Tetrationskern haengt
   daran. */
geopowv(mu, n) = {
  my(v = [1, mu], m, h);
  if (n == 0, return([1]));
  while (#v < n + 1, m = #v; h = min(m - 1, n - m + 1); v = concat(v, v[m]*v[2..h+1]));
  return(v);
}
/* nur nahe dem Kreisrand (|u| >= 0.85): weiter innen ist subst genau und die
   volle komplexe Arithmetik hier etwas schlechter (b = 2, dps 16: H_7(1/2) mit cpe
   ueberall 1e-12.7, mit subst 1e-13.5; H_6 dagegen mit cpe genauer) */
cpe(P, u) = {
  my(v);
  if ((type(u) != "t_COMPLEX") || (type(P) != "t_POL") || (abs(u) < 0.85), return(subst(P, x, u)));
  v = Vecrev(P);
  return(v*geopowv(u, #v - 1)~);
}
hymax = 32;
hyp = vector(hymax); hylam = vector(hymax); hypsi = vector(hymax);
hysig = vector(hymax); hyr = vector(hymax); hyc = vector(hymax);
hystate = vector(hymax, j, [0, 0]); hydig = vector(hymax);
/* exp-089: hykind[r] = 1 regulaer, 2 Kneser (gerade r ueber eta_{r-1});
   knS[r] = Zustand der Kneser-Loesung von Rang r (knS[6] = die Hexation),
   knSd[r] ihre Stellen; knrank = der Rang, den knhexinit gerade baut, und
   knf/knfi = dessen untere Stufe H_{r-1} und ihre Umkehrung */
hykind = vector(hymax); knS = vector(hymax); knSd = vector(hymax);
knrank = 6; knslow = 1;
/* untere Stufe ohne Waechter -- hyensure hat sie schon aufgebaut.
   exp-090: hycall[r] = [vorwaerts, rueckwaerts] zaehlt die Aufrufe -- ein von
   der Rechnerlast unabhaengiges Mass fuer den Aufwand */
hycall = vector(hymax, j, [0, 0]);
hyncalls() = sum(j = 1, hymax, hycall[j][1] + hycall[j][2]);
knf(z) = hylo(knrank - 1, z);
knfi(y) = hyloinv(knrank - 1, y);
/* exp-097: Atlas je Rang. Ein Wert von H_r ist eine Kette von Werten von
   H_{r-1} entlang der Bahn vom Fixpunkt p_r bis zum Ziel, jeder davon wieder
   eine Kette -- die Kosten wachsen je Rang um ein Vielfaches (b = 1.5, dps 100:
   Rang 7/8/9 bauen in 8/13/149 s; H_9(1/2) bei dps 30: 8x H_8, 35x H_7^-1,
   169x H_6^-1). Ein Fenster um x = 1/2 (exp-096) half nicht: die Arbeit liegt
   auf der Bahn, z.B. bei p_9 = -5.56.
   Deshalb ein Raster von Kreisen auf der reellen Achse (Mitten 1/2 + i/2,
   Radius rho = 0.6, benutzt bis 0.75 rho), jeder mit einer Taylorreihe von H_r
   aus N = 2^j Kettenwerten (FFT, Verdoppeln behaelt die alten Werte, bis die
   obersten 20 % der Glieder unter der Genauigkeit der Kette liegen). Gebaut
   wird ein Kreis erst, wenn in seiner Zelle schon N Kettenwerte gerechnet
   wurden -- dann ist sein Aufbau hoechstens das Doppelte dessen, was dort
   ohnehin anfiel, und er entsteht genau da, wo die Arbeit liegt (Bahnen,
   Abtastkreise um Fixpunkte, Kneser-Bilder nahe der Achse), fuer regulaere wie
   Kneser-Raenge. hylo/hyloinv fragen zuerst den Atlas (Horner, Umkehrung per
   Newton), sonst die Kette. Der oberste Rang bekommt keinen Atlas: er wird nur
   von aussen ausgewertet, nicht ueber hylo.
   P = [xc, rho, good, T(u), T'(u), ylo, yhi], u = (x - xc)/rho */
/* erst ab 30 Arbeitsstellen (zwei Maschinenwoerter): bei b = 2, dps 16 kostete die
   Toleranz der Kreise 0.7-0.9 Stellen (H_9-H_11), und noetig ist der Atlas dort nicht */
hyaton = 1; hyatmin = 6; hyatx0 = 1/2; hyath = 1/2; hyatfreeze = 0;
hyat = vector(hymax, j, Map()); hyathit = vector(hymax, j, Map()); hyatlist = vector(hymax, j, List());
/* Schutz um den Fixpunkt p_{r+1}, von dem die Kette des naechsten Rangs ausgeht:
   jeder Schritt dort verstaerkt absolute Fehler um 1/lambda bzw. lambda, und die
   Kreise sind absolut nur so genau wie die Kette insgesamt (b = 6/5, Rang 8:
   lambda = 0.0065, Kreise von H_7 auf 1e-35 statt 1e-38 -> 24.4 statt 31.4
   Stellen). Dort rechnet die Kette. */
hyatprot = vector(hymax, j, List());
hyatnear(r, z) = { for (j = 1, #hyatprot[r], if (abs(z - hyatprot[r][j][1]) < hyatprot[r][j][2], return(1))); return(0); }
hyatbusy = vector(hymax); hyatN = vector(hymax); hyatuse = vector(hymax, j, [0, 0, 0]);
hyatclear(r0) = {
  for (j = max(1, r0), hymax,
    hyat[j] = Map(); hyathit[j] = Map(); hyatlist[j] = List(); hyatfixl[j] = List(); hyatprot[j] = List(); hyatN[j] = 0; hyatbusy[j] = 0; hyatdone[j] = 0);
}
hyatcell(z) = round((real(z) - hyatx0)/hyath);
hyatdig(r) = if (r == 5, pentdigits, if (r == 6, hexdigits, hydig[r]));
hybase(r, z) = if (r == 4, sexp(z), if (r == 5, pent(z), if (r == 6, hex(z), hyval(r, z))));
hybaseinv(r, y) = if (r == 4, slog(y), if (r == 5, ipent(y), if (r == 6, ihex(y), hyinv(r, y))));
hyatget(r, z) = {
  my(i, P);
  for (j = 1, #hyatfixl[r], P = hyatfixl[r][j]; if (abs(z - P[1]) <= P[3], return(P)));
  if (abs(imag(z)) > 0.6, return(0));
  i = hyatcell(z);
  foreach ([i, i + 1, i - 1], j,
    if (mapisdefined(hyat[r], j, &P) && (type(P) == "t_VEC") && (abs(z - P[1]) <= P[3]), return(P)));
  return(0);
}
/* Umkehrung: Newton im Kreis, dessen reeller Wertebereich Re y enthaelt */
hyatnewton(P, y) = {
  my(t = P[1] - P[3] + 2*P[3]*(real(y) - P[6])/(P[7] - P[6]), f, dt, tol = 10.0^(-precis + 2));
  for (it = 1, 100,
    f = cpe(P[4], (t - P[1])/P[2]) - y;
    dt = f/(cpe(P[5], (t - P[1])/P[2])/P[2]);
    t -= dt;
    if (abs(t - P[1]) > 1.05*P[3], return(oo));
    if (abs(dt) < tol, return(t)));
  return(oo);
}
hyatinvget(r, y) = {
  my(t, ry = real(y));
  if (abs(imag(y)) > 0.6, return(oo));
  for (j = 1, #hyatlist[r],
    my(P = hyatlist[r][j]);
    if ((ry >= P[6]) && (ry <= P[7]), t = hyatnewton(P, y); if (t != oo, return(t))));
  return(oo);
}
/* ein Kettenwert in Zelle i: mitzaehlen, bei N/2 den Kreis bauen */
hyatthr(r) = if (hyatN[r], hyatN[r], 2^ceil(log(1.2*(precis + 25))/log(2)))/hyatq;
hyatq = 1;
hyatnote(r, z) = {
  my(i, c = 0);
  if (hyatbusy[r] || (abs(imag(z)) > 0.6), return(0));
  i = hyatcell(z);
  if (mapisdefined(hyat[r], i), return(0));
  if (mapisdefined(hyathit[r], i), c = mapget(hyathit[r], i));
  c++; mapput(hyathit[r], i, c);
  if ((c >= hyatthr(r)) && !hyatfreeze, hyatbuild(r, i));
  return(0);
}
hyatpts(r, xc, rho, N, v) = {
  my(vn);
  if (!v, return(vector(N, t, hybase(r, xc + rho*exp(2*Pi*I*(t - 1)/N)))));
  vn = vector(N);
  for (t = 1, N/2, vn[2*t - 1] = v[t]; vn[2*t] = hybase(r, xc + rho*exp(2*Pi*I*(2*t - 1)/N)));
  return(vn);
}
hyatcoef(v) = {
  my(N = #v, G = fft(geopow(exp(-2*Pi*I/N), N - 1), v));
  return(precision(Polrev(vector(N, k, real(G[k])/N)), precis));
}
hyatpatch(r, xc, rho) = {
  my(good = 0.9*rho, N = 64, v, T, tgt = hyatdig(r), eps, nmax = 8*(precis + 25), yl, yr, d);
  eps = 10.0^(-tgt - 3);
  v = hyatpts(r, xc, rho, N, 0);
  /* frueh aufgeben: aus dem Abfall der Glieder 8..15 gegen 40..51 die noetige
     Punktzahl schaetzen (b = 1.5, Rang 8 bei -5.5 neben der Singularitaet: ein
     vergeblicher Aufbau bis 504 Punkte kostete 98 s, der ganze Rang 9 sonst 10 s) */
  T = hyatcoef(v);
  my(lo = vecmax(vector(8, j, abs(polcoeff(T, 7 + j)))), hi = vecmax(vector(12, j, abs(polcoeff(T, 39 + j)))), q, nest);
  if (hi >= 100*eps,
    q = (hi/lo)^(1/36);
    if (q >= 0.97, error("Atlas: Glieder fallen nicht"));
    nest = 51 + 1.25*log(100*eps/hi)/log(q);
    if (nest > nmax, error("Atlas: braeuchte ~", ceil(nest), " Punkte")));
  while (1,
    T = hyatcoef(v);
    if (vecmax(vector(N - floor(0.8*N), j, abs(polcoeff(T, floor(0.8*N) + j - 1)))) < 100*eps, break);
    if (2*N > nmax, error("Atlas: mehr als ", nmax, " Punkte"));
    N = 2*N; v = hyatpts(r, xc, rho, N, v));
  yl = cpe(T, -good/rho); yr = cpe(T, good/rho);
  if ((abs(imag(yl)) > 1e-20) || (abs(imag(yr)) > 1e-20) || (real(yr) <= real(yl)), error("Atlas: nicht reell steigend"));
  d = max(abs(cpe(T, 0.3) - hybase(r, xc + 0.3*rho)), abs(cpe(T, -0.4 + 0.3*I) - hybase(r, xc + rho*(-0.4 + 0.3*I))));
  if (d > 10.0^(2 - tgt), error("Atlas: Probe ", precision(d, 3)));
  hyatN[r] = N;
  return([xc, rho, good, T, deriv(T), real(yl), real(yr)]);
}
hyatbuild(r, i) = {
  my(xc = hyatx0 + i*hyath, P = 0, tb = getabstime(), msg = "");
  hyatbusy[r] = 1;
  foreach ([0.6, 0.4], rho,
    P = iferr(hyatpatch(r, xc, rho), E, hyatbusy[r] = 0; if (errname(E) == "e_ALARM", error(E)); msg = Str(E); 0);
    hyatbusy[r] = 1;
    if (P, break));
  hyatbusy[r] = 0;
  mapput(hyat[r], i, if (P, P, -1));
  if (P, listput(hyatlist[r], P));
  if (kndebug >= 2, print("hy: Atlas Rang ", r, " Zelle ", xc, ": ", if (P, Str("Grad ", poldegree(P[4]), ", rho ", P[2]), Str("keiner (", msg, ")")), "  ", getabstime() - tb, " ms"));
  return(P);
}
/* Fixpunkt-Kreis: die Reihe von H_r um den Fixpunkt p von H_r, die hyinit(r+1)
   ohnehin abtastet (pentfixser), als Kreis mit eigener Mitte -- er deckt
   gerade die Stelle neben der Singularitaet, an der Gitterkreise scheitern
   (b = 1.5: p_9 = -5.56). S(t) = H_r(p + t) - p, Radius rho. */
hyatfixl = vector(hymax, j, List());
hyatfix(r, p, S, rho) = {
  my(T, P, good = 0.75*rho, yl, yr, d, tgt = hyatdig(r));
  if (!hyaton || (r < hyatmin) || (precis < 30) || (imag(p) != 0) || (rho <= 0), return(0));
  T = precision(p + subst(S, x, rho*x), precis);
  yl = real(subst(T, x, -0.75)); yr = real(subst(T, x, 0.75));
  if (yr <= yl, return(0));
  d = abs(subst(T, x, 0.5) - hybase(r, p + 0.5*rho));
  if (d > 10.0^(2 - tgt), return(0));
  P = [p, rho, good, T, deriv(T), yl, yr];
  listput(hyatfixl[r], P); listput(hyatlist[r], P);
  if (kndebug >= 2, print("hy: Atlas Rang ", r, " Fixpunkt-Kreis um ", precision(p, 6), ", rho ", precision(rho, 4), ", Grad ", poldegree(T)));
  return(P);
}
/* alle Gitterzellen von lo bis hi belegen (ohne Schwelle) */
hyatcover(r, lo, hi) = {
  my(tb = getabstime(), n0 = #hyatlist[r]);
  if (!hyaton || (r < hyatmin) || hyatfreeze, return(0));
  for (i = ceil((lo - hyatx0)/hyath), floor((hi - hyatx0)/hyath),
    if (!mapisdefined(hyat[r], i), hyatbuild(r, i)));
  if (kndebug, print("hy: Atlas Rang ", r, " auf [", precision(lo, 4), ", ", precision(hi, 4), "]: ", #hyatlist[r] - n0, " neue Kreise, ", getabstime() - tb, " ms"));
  return(0);
}
/* Bahnbereich von Rang j fuer den Rang darueber: zwischen dessen Fixpunkt und [0, 1] */
hyatrange(j) = {
  my(q = if (hykind[j+1] == 2, real(knS[j+1][3]), real(hyp[j+1])));
  return([min(q, 0) - 0.75, max(q, 1) + 0.75]);
}
hyatcov = 99;
hyatfromr = 11; hyatdone = vector(hymax);
hyattol = 1;  /* vollstaendige Belegung: aus, s. exp-097 */
/* Atlas per Fortpflanzung (exp-097, "Verkettung"): die Kreise eines regulaeren
   Rangs r werden nicht ueber die Kette abgetastet, sondern ueber die
   Funktionalgleichung aus dem Nachbarkreis im Abstand 1:
     H_r(z) = H_{r-1}(H_r(z - 1))      (ungerade r, von links)
     H_r(z) = H_{r-1}^-1(H_r(z + 1))   (gerade r, von rechts)
   -- auf den Abtastpunkten Horner im Nachbarkreis und Horner/Newton in einem
   Kreis von H_{r-1}. Begonnen wird dort, wo H_r = p + psi(c lambda^z) ohne
   Kette gilt (|c lambda^z| <= rr/2 auf dem ganzen Kreis), je ein Kreis fuer
   beide Raster-Paritaeten. Gewandert wird bis 2.5 (ungerade) bzw. bis ein
   Kreis scheitert (gerade, links liegt die Singularitaet). */
hyatptsf(f, xc, rho, N, v) = {
  my(vn);
  if (!v, return(vector(N, t, f(xc + rho*exp(2*Pi*I*(t - 1)/N)))));
  vn = vector(N);
  for (t = 1, N/2, vn[2*t - 1] = v[t]; vn[2*t] = f(xc + rho*exp(2*Pi*I*(2*t - 1)/N)));
  return(vn);
}
/* Kreis aus einer Wertefunktion f auf dem Kreis um xc; Probe gegen die Kette */
hyatfrom(r, f, xc, rho, probe) = {
  my(good = 0.9*rho, N = 64, v, T, tgt = hyatdig(r), eps, nmax = 8*(precis + 25), yl, yr, d, lo, hi, q, nest);
  eps = 10.0^(-tgt - 3);
  v = hyatptsf(f, xc, rho, N, 0);
  T = hyatcoef(v);
  lo = vecmax(vector(8, j, abs(polcoeff(T, 7 + j)))); hi = vecmax(vector(12, j, abs(polcoeff(T, 39 + j))));
  if (hi >= 100*eps,
    q = (hi/lo)^(1/36);
    if (q >= 0.97, error("Atlas: Glieder fallen nicht"));
    nest = 51 + 1.25*log(100*eps/hi)/log(q);
    if (nest > nmax, error("Atlas: braeuchte ~", ceil(nest), " Punkte")));
  /* bis die obersten 20 % der Glieder unter 100 eps liegen -- oder am Rauschboden
     der Werte stehen bleiben: nach links verstaerkt H_{r-1}^-1 am anziehenden
     Fixpunkt jeden Fehler um 1/lambda (b = 1.5, Rang 8: 17.5, 1.2 Stellen je
     Schritt, wie in der Kette auch). Dann gilt der Kreis mit seinem Boden
     (P[8]), hoechstens 5 Stellen unter der Kette. */
  my(tl = hyattail(T, N), tlast = oo, off);
  while (1,
    if (tl >= 100*eps,
      /* weiter verdoppeln, solange der Schwanz deutlich faellt und Punkte frei
         sind; sonst annehmen, wenn er hoechstens tgt/6 Stellen (1 bis 5) unter
         der Kette liegt (Fehler steht in P[8]) -- nahe den Singularitaeten fallen
         die Glieder langsam (b = 1.5, Rang 10 bei 0: 504 Punkte, Schwanz 9e-31).
         Mit fest 5 Stellen trug H_11 bei b = 2, dps 16 nur 6.95 Stellen, und
         Rang 12 (braucht 8) ging nicht mehr. */
      if ((4*tl <= tlast) && (2*N <= nmax),
        tlast = tl;
        N = 2*N; v = hyatptsf(f, xc, rho, N, v); T = hyatcoef(v); tl = hyattail(T, N);
        next);
      if (tl >= 10.0^(hyattol*min(5, max(1, tgt/6)) - tgt), error("Atlas: Schwanz ", precision(tl, 3), " bei ", N, " Punkten")));
    /* Zwischen den Abtastpunkten pruefen: der Schwanz kann bei Unterabtastung
       klein aussehen (b = 1.5, Rang 8, Startkreis um 5/2: Schwanz 1e-39 bei 128
       Punkten, zwischen zwei Punkten 1e-27 daneben) */
    off = vecmax(vector(4, j, my(th = 2*Pi*((j - 1)*N/4 + 0.5 + 0.13*j)/N);
                         abs(cpe(T, exp(I*th)) - f(xc + rho*exp(I*th)))));
    if (off <= 100*max(tl, eps), break);
    if (2*N > nmax, error("Atlas: zwischen den Punkten ", precision(off, 3), " daneben"));
    N = 2*N; v = hyatptsf(f, xc, rho, N, v); T = hyatcoef(v); tl = hyattail(T, N); tlast = oo);
  yl = cpe(T, -good/rho); yr = cpe(T, good/rho);
  if ((abs(imag(yl)) > 1e-20) || (abs(imag(yr)) > 1e-20) || (real(yr) <= real(yl)), error("Atlas: nicht reell steigend"));
  if (probe,
    d = abs(cpe(T, 0.3) - hyval(r, xc + 0.3*rho));
    if (d > max(10.0^(2 - tgt), 100*tl), error("Atlas: Probe ", precision(d, 3))));
  hyatN[r] = max(hyatN[r], N);
  return([xc, rho, good, T, deriv(T), real(yl), real(yr), max(tl, eps)]);
}
hyattail(T, N) = vecmax(vector(N - floor(0.8*N), j, abs(polcoeff(T, floor(0.8*N) + j - 1))));
hyatput(r, i, P) = { mapput(hyat[r], i, if (P, P, -1)); if (P, listput(hyatlist[r], P)); }
hyatprop(r) = {
  my(la = hylam[r], up, c = hyc[r], zs, z0, i0, dir, rho = 0.6, n0 = #hyatlist[r], tb = getabstime(), P, Pn, i, j, xcj, f, nfail = 0);
  if (!hyaton || (hykind[r] != 1) || (precis < 30), return(0));
  up = (abs(la) > 1); dir = if (up, 1, -1);
  zs = log(hyr[r]/(2*abs(c)))/log(abs(la));
  z0 = if (up, zs - rho, zs + rho);
  i0 = if (up, floor((z0 - hyatx0)/hyath), ceil((z0 - hyatx0)/hyath));
  hyatbusy[r] = 1;
  for (par = 0, 1,
    i = i0 - dir*par;
    my(xc = hyatx0 + i*hyath);
    f = (z -> hyp[r] + cpe(hypsi[r], c*la^z));
    P = iferr(hyatfrom(r, f, xc, rho, 0), E, if (kndebug >= 2, print("hy:   Rang ", r, " Startkreis ", xc, ": ", E)); 0, errname(E) != "e_ALARM");
    if (!P, nfail++; next);
    if (!mapisdefined(hyat[r], i), hyatput(r, i, P));
    for (step = 1, 40,
      j = i + 2*dir; xcj = hyatx0 + j*hyath;
      if (up && (xcj > 2.5), break);
      if (!up && (xcj < -12), break);
      f = if (up, (z -> hylo(r - 1, cpe(P[4], (z - 1 - P[1])/P[2]))),
                  (z -> hyloinv(r - 1, cpe(P[4], (z + 1 - P[1])/P[2]))));
      /* der Radius darf beim Wandern schrumpfen (der neue Kreis liegt dann im
         alten), nicht wachsen: bei zu langsamem Abfall 2/3 davon, bis 0.26 */
      my(rr2 = min(rho, P[2]));
      Pn = 0;
      while (!Pn && (rr2 >= 0.26),
        Pn = iferr(hyatfrom(r, f, xcj, rr2, (step % 4 == 0)), E, if (kndebug >= 2, print("hy:   Rang ", r, " Kreis ", xcj, " rho ", precision(rr2, 3), ": ", E)); 0, errname(E) != "e_ALARM");
        if (!Pn, rr2 = 2*rr2/3));
      if (!Pn, nfail++; break);
      if (!mapisdefined(hyat[r], j), hyatput(r, j, Pn));
      P = Pn; i = j));
  hyatbusy[r] = 0;
  if (kndebug, print("hy: Atlas Rang ", r, " fortgepflanzt: ", #hyatlist[r] - n0, " Kreise (", nfail, " Abbrueche), Mitten ",
                     precision(vecmin(vector(#hyatlist[r], k, real(hyatlist[r][k][1]))), 3), "..",
                     precision(vecmax(vector(#hyatlist[r], k, real(hyatlist[r][k][1]))), 3), ", ", getabstime() - tb, " ms"));
  return(#hyatlist[r] - n0);
}
hylo(r, z) = {
  my(P, v);
  hycall[r][1]++;
  if (!hyaton || (r < hyatmin) || (precis < 30) || hyatnear(r, z), return(hybase(r, z)));
  P = hyatget(r, z);
  if (P, hyatuse[r][1]++; return(cpe(P[4], (z - P[1])/P[2])));
  v = hybase(r, z); hyatnote(r, z);
  return(v);
}
hyloinv(r, y) = {
  my(t, v);
  hycall[r][2]++;
  if (!hyaton || (r < hyatmin) || (precis < 30) || hyatnear(r, y), return(hybaseinv(r, y)));
  t = hyatinvget(r, y);
  if (t != oo, hyatuse[r][2]++; return(t));
  v = hybaseinv(r, y); hyatnote(r, v);
  return(v);
}
hyatinfo(r) = [#hyatlist[r], hyatuse[r], if (#hyatlist[r], vecsort(vector(#hyatlist[r], j, hyatlist[r][j][1])), [])];
hyval(r, z, extra) = {
  my(la = hylam[r], w, k = 0, y, up = (abs(hylam[r]) > 1));
  if (hykind[r] == 2, return(knhexS(knS[r], z)));
  w = hyc[r]*la^z;
  while (abs(w) > hyr[r],
    w = if (up, w/la, w*la); k++;
    if (k > 400, error("hyper: Rang ", r, ", Kette laenger als 400")));
  for (j = 1, extra, w = if (up, w/la, w*la); k++);
  y = hyp[r] + subst(hypsi[r], x, w);
  if (up, for (j = 1, k, y = hylo(r-1, y)), for (j = 1, k, y = hyloinv(r-1, y)));
  return(y);
}
hyinv(r, y) = {
  my(la = hylam[r], k = 0, up = (abs(hylam[r]) > 1));
  if (hykind[r] == 2, return(knihexS(knS[r], y)));
  while (abs(y - hyp[r]) > hyr[r],
    y = if (up, hyloinv(r-1, y), hylo(r-1, y)); k++;
    if (k > 400, error("ihyper: Rang ", r, ", Kette erreicht den Fixpunkt nicht")));
  return(log(subst(hysig[r], x, y - hyp[r])/hyc[r])/log(la) + if (up, k, -k));
}
hyinit(r) = {
  my(odd = r % 2, y = 1, y2, bad, ok = 0, res = 0, sc, la, R, sr, rr, ps, sig, eps, target,
     k, fn, isreal, err, zz, z1, z2, lowd, b = sexp(1));
  lowd = if (r == 7, hexdigits, hydig[r-1]);
  /* komplexe Basen: die Hexation ist dort nur Sheldons Naeherung (b = 2+i/2:
     0.8 Stellen) -- darauf laesst sich kein hoeherer Rang bauen */
  if (lowd < 8, error("hyper: Rang ", r, ": H_", r-1, " traegt nur ", precision(lowd, 3),
                      " Stellen, darauf baut kein hoeherer Rang"));
  target = max(8, min(precis - 10, floor(lowd) - 5));
  eps = 10.0^(-target - 3);
  hykind[r] = 0;
  /* Fixpunkt: ungerade r von 1 aus mit der Umkehrung (sie zieht zum linken
     abstossenden Fixpunkt), gerade r vorwaerts (zum anziehenden). Waechst
     die Vorwaertsbahn davon, gibt es keinen reellen: b > eta_{r-1}. */
  for (i = 1, 2000,
    bad = 0;
    y2 = iferr(if (odd, hyloinv(r-1, y), hylo(r-1, y)), E, bad = 1; 0);
    if (bad || (abs(y2) > 1e6), break);
    if ((imag(b) == 0) && (abs(imag(y2)) > 1e-10*(1 + abs(y2))), break);
    if (abs(y2 - y) < 1e-12*(1 + abs(y)), y = y2; ok = 1; break);
    y = y2);
  if (!ok,
    if (odd, error("hyper: Rang ", r, ": kein linker Fixpunkt von H_", r-1, " gefunden"));
    /* exp-089: gerades r ueber eta_{r-1} -- Kneser am komplexen Fixpunkt von
       H_{r-1}, dieselbe Konstruktion wie die Hexation auf der Pentation */
    if (imag(b) != 0, error("hyper: Rang ", r, ": kein anziehender Fixpunkt, und Kneser nur fuer reelle Basen"));
    knhexinit(r);
    hykind[r] = 2;
    return(precision(min(knSd[r], lowd), 9)));
  isreal = (imag(b) == 0) && (abs(imag(y)) < 1e-10);
  if (isreal, y = real(y));
  if (hyaton, listput(hyatprot[r-1], [y, 0.25]));
  /* exp-097: ab Rang hyatcov den Atlas der unteren Stufen auf ihren Bahnen
     vollstaendig belegen, von unten nach oben -- dann kostet jeder Kettenwert
     von H_{r-1} nur Horner-Schritte */
  if (hyaton && (r >= hyatcov),
    for (j = hyatmin, r - 2, my(rg = hyatrange(j)); hyatcover(j, rg[1], rg[2]));
    hyatcover(r - 1, min(real(y), 0) - 0.75, max(real(y), 1) + 0.75));
  fn = if (odd, (t -> hyloinv(r-1, t)), (t -> hylo(r-1, t)));
  sr = 1/8;
  for (i = 1, 10,
    res = iferr(pentfixser(fn, y, sr, isreal, eps, Str("hyper Rang ", r)), E, 0);
    if (type(res) == "t_VEC", break);
    sr = sr/2);
  if (type(res) != "t_VEC", error("hyper: Rang ", r, ": H_", r-1, " am Fixpunkt nicht entwickelbar"));
  hyp[r] = res[1]; sc = res[2]; R = res[3];
  /* exp-097: die Reihe am Fixpunkt als Kreis von H_{r-1} (ungerade r: sc ist
     die Umkehrung, vorwaerts per serreverse, Radius ~ R/lambda) */
  if (hyaton && (r - 1 >= hyatmin),
    iferr(hyatfix(r - 1, hyp[r], if (odd, Pol(serreverse(sc + O(x^(poldegree(sc) + 1)))), sc),
                  if (odd, 0.4*R*abs(polcoeff(sc, 1)), 0.5*R)), E, 0, errname(E) != "e_ALARM"));
  la = polcoeff(sc, 1);
  if (odd, la = 1/la);
  if (odd && (abs(la) <= 1), error("hyper: Rang ", r, ": Fixpunkt nicht abstossend, |lambda| = ", precision(abs(la), 6)));
  if (!odd && (abs(la) >= 1), error("hyper: Rang ", r, ": Fixpunkt nicht anziehend, |lambda| = ", precision(abs(la), 6)));
  rr = R/4;
  for (i = 1, 8,
    ps = pentpsi(sc, rr, eps, poldegree(sc) - 1);
    if (ps[2], break);
    rr = rr/2);
  if (!ps[2], error("hyper: Rang ", r, ": Poincare-Reihe konvergiert nicht"));
  sig = Pol(serreverse(ps[1] + O(x^(poldegree(ps[1])+1))));
  if (abs(polcoeff(sig, poldegree(sig)))*rr^poldegree(sig) > eps,
    error("hyper: Rang ", r, ": Umkehrreihe konvergiert bei r = ", precision(rr, 5), " nicht"));
  hypsi[r] = pentclean(ps[1], rr, eps, "Psi");
  hysig[r] = pentclean(sig, rr, eps, "sigma");
  hylam[r] = la; hyr[r] = rr; hykind[r] = 1;
  /* Normierung H_r(0) = 1: von 1 aus in den Radius, dann zurueckskalieren */
  y = 1; k = 0;
  while (abs(y - hyp[r]) > rr,
    y = if (odd, hyloinv(r-1, y), hylo(r-1, y)); k++;
    if (k > 400, error("hyper: Rang ", r, ": Normierungskette laenger als 400")));
  hyc[r] = subst(hysig[r], x, y - hyp[r])*la^(if (odd, k, -k));
  /* Fehlerschaetzung wie bei pentinit: H(0) = 1, H(1) = b, Rundweg und
     dieselben Werte ueber zwei Reduktionstiefen */
  err = max(abs(hyval(r, 0) - 1), abs(hyval(r, 1) - b)/abs(b));
  err = max(err, abs(hyinv(r, hyval(r, 1/3)) - 1/3));
  for (j = 1, 2,
    zz = if (j == 1, 1/2, -1/2);
    z1 = hyval(r, zz); z2 = hyval(r, zz, 1);
    err = max(err, abs(z1 - z2)/abs(z1)));
  err = if (err == 0, precis, min(precis, -log(err)/log(10)));
  if (quietmode == 0,
    print("hyper Rang ", r, ": ", if (odd, "abstossender", "anziehender"), " Fixpunkt p = ", precision(hyp[r], 15),
          ", lambda = ", precision(la, 12), ", ", precision(min(err, lowd), 4), " Stellen"));
  return(precision(min(err, lowd), 9));
}
hyensure(r) = {
  if (r <= 4, return(precis));
  if (r == 5, return(pentensure()));
  if (r == 6, return(hexensure()));
  if (r > hymax, error("hyper: Rang ", r, " > hymax = ", hymax));
  if (hystate[r] != [sexp(1), precis],
    hyensure(r - 1);
    hyatclear(r); hyatfreeze = 0;
    /* exp-097: ab Rang hyatfromr = 11 bekommen alle regulaeren Raenge darunter
       ihren Atlas (von unten nach oben, einmal) -- darunter kostet er mehr, als er
       spart (b = 1.5, dps 30: Rang 9/10 ohne 9/58 s, mit 50/138 s; Rang 11 ohne
       ~480 s, mit 68 s; b = 6/5, Rang 8: Atlas von H_7 75 s, der Rang 8 s) */
    if (hyaton && (r >= hyatfromr),
      for (j = 7, r - 1, if ((hykind[j] == 1) && !hyatdone[j], hyatprop(j); hyatdone[j] = 1)));
    hydig[r] = hyinit(r);
    hystate[r] = [sexp(1), precis]);
  return(hydig[r]);
}
hyrank(r) = {
  if ((type(r) != "t_INT") || (r < 0), error("hyper: der Rang muss eine ganze Zahl >= 0 sein, nicht ", r));
  return(r);
}
hyper(r, z) = {
  hyrank(r);
  if (r == 0, return(z + 1));
  if (r == 1, return(sexp(1) + z));
  if (r == 2, return(sexp(1)*z));
  if (r == 3, return(sexp(1)^z));
  if (r == 4, return(levelchk(sexp(z), "hyper")));
  if (r == 5, return(pentx(z)));
  if (r == 6, return(hexx(z)));
  hyensure(r);
  return(levelchk(hyval(r, z), "hyper"));
}
ihyper(r, y) = {
  hyrank(r);
  if (r == 0, return(y - 1));
  if (r == 1, return(y - sexp(1)));
  if (r == 2, return(y/sexp(1)));
  if (r == 3, return(log(y)/log(sexp(1))));
  if (r == 4, return(levelchk(slog(y), "ihyper")));
  if (r == 5, return(ipentx(y)));
  if (r == 6, return(ihexx(y)));
  hyensure(r);
  return(levelchk(hyinv(r, y), "ihyper"));
}
hyperdigits(r) = { hyrank(r); return(hyensure(r)); }
/* exp-091: Tabellen fuer den Rechner -- H_r auf dem Gitter x0 + (i-1)*step,
   Punkt fuer Punkt in aufsteigender Reihenfolge abgefragt (i = 1, 2, ...).
   Ist der Punkt eine Periode frueher reell bekannt, reicht ein Schritt der
   unteren Stufe, H_r(x + 1) = H_{r-1}(H_r(x)), statt einer eigenen Kette
   (b = 1.5, r = 9: 0.57 s je Wert direkt). Nicht reelle Werte (gerade Raenge
   links) werden nicht weitergereicht; dort rechnet der naechste Punkt direkt. */
hygtab = [0, 0, 0, []];
hypergridval(r, x0, step, i) = {
  my(per = round(1/step), p = 0, y);
  hyrank(r);
  if ((i == 1) || (hygtab[1] != r) || (hygtab[2] != x0) || (hygtab[3] != step) || (i > #hygtab[4]),
    if (i == 1, hygtab = [r, x0, step, vector(8192)], error("hypergridval: Punkte der Reihe nach ab i = 1 abfragen")));
  if (r >= 5, hyensure(r));
  if (i > per, p = hygtab[4][i - per]);
  y = if ((r >= 5) && (type(p) == "t_REAL"), levelchk(hylo(r - 1, p), "hyper"), hyper(r, x0 + (i - 1)*step));
  if (abs(imag(y)) <= 1e-9*(1 + abs(y)), hygtab[4][i] = real(y)*1.);
  return(y);
}

/* exp-077 (2026-10-09): regulaere Hexation fuer reelle 1 < b < eta5.
   Dort schneidet die Pentation y = x zweimal (b = 1.70: bei ~2.44 und ~5.69);
   der untere Schnitt a ist anziehend, P'(a) = mu in (0,1). Wie die Tetration
   unter eta (regsexp) ist die Hexation dort die regulaere Iteration:
   Psi_H(mu*w) = P(Psi_H(w)), hex(z) = Psi_H(kappa*mu^z), kappa = sigma_H(1),
   also hex(0) = 1, hex(z+1) = pent(hex(z)), hex(z) -> a fuer z -> oo.
   Sheldons hex/ihex (ohne theta, Newton auf 1E-15) liefern ~8 Stellen und
   bleiben fuer b > eta5 zustaendig. */
hexreg=0;
knforce=0;
/* exp-094: knlsqpass = 1 (ein MGS-Durchgang, rechte Seite nacheinander, s. knlsq).
   knsolvep wirkt nicht: PARI rechnet mit der Genauigkeit der Operanden, nicht mit
   realprecision, und die Matrix hat precis Stellen. */
knptsf=0.8;  /* Bogenpunkte je Unbekannte in knsolve */
knlsqpass=1; knsolvep=2; knlsqcond=0;
/* exp-095: knmode = 1 loest jeden Bogen per Wiederholung (knisolve, wie loop()
   in fatou.gp) statt per Ausgleichsrechnung (knsolvelsq); knmaxki ersetzt dann
   knmaxk. Regler: Kreispunkte je Glied, Grenze |Bild|/R fuer die
   Funktionalgleichung, Anderson-Tiefe, hoechstens so viele Schritte. */
knmode=1; knfast=0; knmaxki=3000; knis=0.7; kniqmax=0.92; kniumax=0.9; knimem=20; knimaxit=300;
kpsi=0; kur=0; kntarget=0; knistat=0;
knreachtol=1.02;
kndebug=0;
knRmeas=0;
knmaxk=600;
reghexinit(sr) = {
  my(y, y2, bad, ok, res, eps, rr, sc, ps, R, n, z1, z2, zz);
  hexreg = 0;
  if (sr==0, sr=1/2);
  eps = 10.0^(-precis-8);
  if (imag(sexp(1)) != 0, return(0));          /* komplexe Basis: kein reeller Fixpunkt */
  /* anziehender Fixpunkt: von 1 aus iterieren. Waechst P^n(1) davon oder
     laeuft ueber, gibt es keinen (b > eta5). */
  y = 1; ok = 0;
  for (i = 1, 5000,
    bad = 0;
    y2 = iferr(pent(y), E, bad = 1; 0);
    if (bad || (abs(y2) > 1e6), break);
    if (abs(y2 - y) < 1e-8, y = y2; ok = 1; break);
    y = y2;
  );
  if (!ok, return(0));
  res = pentfixser(z->pent(z), y, sr, 1, eps, "reghexinit");
  hexa = res[1]; sc = res[2]; R = res[3];
  hexmu = polcoeff(sc, 1);
  if (!((hexmu > 0) && (hexmu < 1)), error("reghexinit: P'(a) = ", precision(hexmu,10), " liegt nicht in (0,1)"));
  /* Auswerteradius so gross wie moeglich: jeder fehlende Faktor kostet einen
     ipent-Schritt, und der verstaerkt um 1/mu (b = 6/5: ~40). */
  rr = R/4;
  for (k = 1, 6,
    ps = pentpsi(sc, rr, eps, poldegree(sc) - 1);
    if (ps[2], break);
    rr = rr/2;
  );
  if (!ps[2], error("reghexinit: Poincare-Reihe der Pentation konvergiert nicht"));
  hexr = rr;
  hexsig = Pol(serreverse(ps[1] + O(x^(poldegree(ps[1])+1))));
  if (abs(polcoeff(hexsig, poldegree(hexsig)))*rr^poldegree(hexsig) > eps,
    error("reghexinit: Umkehrreihe konvergiert bei r = ", precision(rr,5), " nicht ausreichend"));
  hexpsi = pentclean(ps[1], rr, eps, "Psi_H");
  hexsig = pentclean(hexsig, rr, eps, "sigma_H");
  /* kappa = sigma_H(1): von 1 aus vorwaerts in den Radius, dann zurueckskalieren */
  y = 1; n = 0;
  while (abs(y - hexa) > rr, y = pent(y); n++);
  hexkappa = subst(hexsig, x, y - hexa)/hexmu^n;
  hexreg = 1;
  /* Fehlerschaetzung wie bei pentinit */
  hexerr = abs(hex(0) - 1);
  hexerr = max(hexerr, abs(ihex(hex(1/3)) - 1/3));
  for (j = 1, 2,
    zz = if (j==1, 1/2, -1/2);
    z1 = reghex(zz); z2 = reghex(zz, 1);
    hexerr = max(hexerr, abs(z1-z2)/abs(z1));
  );
  hexerr = if (hexerr==0, precis, min(precis, -log(hexerr)/log(10)));
  hexerr = precision(hexerr, 9);
  if (quietmode==0,
    print("hexation regulaer am anziehenden Fixpunkt der Pentation");
    print("  a = ", hexa, "   P'(a) = ", hexmu);
    print("  hex precision ", hexerr, " decimal digits");
  );
  return(1);
}
reghex(z, extra) = {
  my(w, m, y);
  w = hexkappa*hexmu^z;
  m = 0;
  while (abs(w) > hexr, w = w*hexmu; m++);
  for (i = 1, extra, w = w*hexmu; m++);
  y = hexa + subst(hexpsi, x, w);
  for (i = 1, m, y = ipent(y));
  return(y);
}
regihex(y) = {
  my(n);
  n = 0;
  while (abs(y - hexa) > hexr, y = pent(y); n++; if (n > 500, error("regihex: pent-Kette erreicht a nicht")));
  return(log(subst(hexsig, x, y - hexa)/hexkappa)/log(hexmu) - n);
}

/* exp-079 (2026-10-09): Kneser-Hexation fuer reelle b > eta5 (2, e, 10, ...).
   Dort hat die Pentation keinen reellen Fixpunkt, aber ein konjugiertes Paar
   p, conj(p), abstossend (b = 2: p = -1.99295 + 1.44513 i, P'(p) = 1.369 +
   2.488 i). Wie Kneser fuer exp: die Abel-Funktion alpha der Hexation
   (alpha(pent(z)) = alpha(z) + 1, reell auf der reellen Achse, alpha(1) = 0)
   ist oben alpha_p + theta(alpha_p), alpha_p = log(sigma_p)/log(lambda) die
   regulaere Abel-Funktion an p, theta 1-periodisch und zu p hin abklingend.
   Sheldons matrix_h laesst theta weg (~8 Stellen). Hier EIN lineares
   Ausgleichsproblem: reelle Taylorreihe A um c = Re p, auf einem oberen Bogen
     A(z) = w + sum_n c_n e^(2 pi i s n w), w = alpha_p(z)   (nahe p),
     A(z) - A(ipent z) = 1  bzw.  A(pent z) - A(z) = 1       (nahe der Achse),
   und eine Normierung; die Konjugationssymmetrie steckt in den reellen
   Koeffizienten. Eindeutig bis auf eine reelle Konstante (Maximumprinzip auf
   dem Halbzylinder), die die Normierung festlegt. Basis 2, Prototyp exp-078:
   ~26 Stellen bei \p 30, Sheldon 7e-9 daneben. */

/* min |M x - b| per modifiziertem Gram-Schmidt -- ohne die quadrierte
   Kondition der Normalgleichungen, an denen Basis e scheiterte. PARIs matqr
   nimmt nur quadratische Matrizen.
   exp-094: EIN Durchgang (knlsqpass = 1) und Q^T b NACHEINANDER wie eine
   weitere Spalte (c_j = q_j.w, w -= c_j q_j). Ein MGS-Durchgang laesst Q um
   ~ eps*Kondition schief, aber so angewandt ist es Bjoercks rueckwaertsstabile
   MGS-Loesung; erst als Produkt Q^T b braucht es den zweiten Durchgang.
   Gemessen b = 2, dps 100, 1275x644: zwei Durchgaenge + Produkt 88 s, Rest
   |b - Mx| 1.7e-102; ein Durchgang nacheinander 47 s, Rest 1.4e-109, x gleich
   bis 1e-95. Die Zeit steckt in den k-1 Vektorbefehlen je Spalte im
   Interpreter, kaum in der Stellenzahl (38/57/77/115 Stellen: 28/32/36/47 s);
   gekuerzte Zerlegung + Nachbesserung stockt (38 St.: x bis 1e-49, 57 St.:
   1e-83; Rest ~ eps_f^2) und lohnt nicht. knlsqpass = 2 bleibt als Schalter. */
knlsq(M, b, passes = knlsqpass) = {
  my(m = matsize(M)[1], n = matsize(M)[2], Q = vector(n), R = matrix(n, n), q, r, xs, cb, w);
  for (k = 1, n,
    q = M[,k];
    for (pass = 1, passes, for (j = 1, k-1, r = Q[j]~ * q; R[j,k] += r; q -= r*Q[j]));
    R[k,k] = sqrt(norml2(q));
    if (R[k,k] == 0, error("knlsq: Spalte ", k, " linear abhaengig"));
    Q[k] = q/R[k,k];
  );
  cb = vector(n); w = b;
  for (j = 1, n, cb[j] = Q[j]~ * w; w -= cb[j]*Q[j]);
  xs = vector(n)~;
  forstep (k = n, 1, -1, xs[k] = (cb[k] - if (k < n, R[k, k+1..n]*xs[k+1..n], 0))/R[k,k]);
  /* Groessenordnung der Kondition aus der R-Diagonale (Spalten sind auf
     Betrag 1 normiert); nur eine untere Schranke, fuer die Ausgabe */
  knlsqcond = vecmax(vector(n, j, abs(R[j,j])))/vecmin(vector(n, j, abs(R[j,j])));
  return(xs);
}

/* regulaere Abel-Funktion an p, Hauptzweig; der Aufrufer macht die Zweige
   entlang einer Punktfolge stetig (Mehrdeutigkeit kper = 2 pi i/log lambda). */
knalphap(z) = {
  my(m = 0);
  while (abs(z - kp) > krs, z = knfi(z); m++; if (m > 60, error("knalphap: ipent laeuft nicht gegen p")));
  return(log(subst(ksig, x, z - kp))/klog + m);
}
knwalk(pts) = {
  my(W = vector(#pts), v, best, d, tw);
  for (j = 1, #pts,
    tw = getabstime();
    if (kndebug >= 3, print("kn:     alpha_p ", j, "/", #pts, " bei ", precision(pts[j], 8)));
    W[j] = knalphap(pts[j]);
    if ((kndebug >= 2) && (getabstime() - tw > 1000),
      print("kn:     alpha_p langsam: ", getabstime() - tw, " ms bei ", precision(pts[j], 8)));
    if (j > 1,
      best = W[j]; d = abs(W[j] - W[j-1]);
      for (k = -3, 3, v = W[j] + k*kper; if (abs(v - W[j-1]) < d, d = abs(v - W[j-1]); best = v));
      W[j] = best));
  return(W);
}
knarc(rad, M) = vector(M, j, knc + rad*exp(I*(0.05 + (Pi - 0.1)*(j-1)/(M-1))));
knen(w, n) = exp(2*Pi*I*kns*n*w);
/* Huellkurve der Koeffizienten |a_k| .. |a_k+9| (gegen zufaellig kleine
   Glieder). Keine anonyme Funktion in knhexinit: dort schluckte "(k) -> ..."
   den ganzen Rest des Funktionsrumpfs, knhexinit gab die Closure zurueck. */
knenv(A, k) = vecmax(vector(10, j, abs(polcoeff(A, k + j - 1))));
/* exp-082: Abfallgesetz log|a_k| ~ a + b*k + c*log(k) aus der Huellkurve auf
   k in [K/4, 3K/4]. Die reine Zweipunkt-Steigung ueberschaetzte den Abfall,
   wenn ein Potenzfaktor k^c mitspielt (Randsingularitaeten), und K kam zu
   klein heraus (b = e, dps 45: Rest 6e-43 statt ~1e-54). Rueckgabe [a, b, c]. */
knfit(A, K) = {
  my(ks = vector(13, j, K\4 + ((j - 1)*(K\2))\12), X, y);
  X = matrix(#ks, 3, i, j, [1, ks[i], log(ks[i])][j]);
  y = vector(#ks, i, log(knenv(A, ks[i])))~;
  return(matsolve(X~*X, X~*y)~);
}
/* kleinstes K mit a + b*K + c*log(K) <= ziel (b < 0) */
knkfor(f, ziel, K0) = {
  my(K = K0);
  for (it = 1, 50, my(Kn = ceil((ziel - f[1] - f[3]*log(K))/f[2])); if (Kn == K, break); K = max(Kn, 8));
  return(K);
}
knser(z) = subst(knA, x, (z - knc)/knrho);
knserd(z) = subst(knAd, x, (z - knc)/knrho)/knrho;

/* Vorlauf fuer einen Bogenradius: Hoehenprofil s*Im alpha_p, Bilder der
   achsennahen Punkte, geschaetzte Glieder. Rueckgabe [kosten, K, N, rad, hmin,
   htop, reichweite, mittelpunkt]. K setzt voraus, dass alpha um den
   Mittelpunkt knc bis |knc - p| analytisch ist, N eine Abklinghoehe
   0.6*(htop - hmin) fuer theta -- beides Faustregeln (theta hat an den Ecken
   der reellen Kurve Randsingularitaeten, die Koeffizienten fallen langsamer;
   Befund der Wissensgraph-Sitzung 09.10.2026). Deshalb misst knhexinit die
   erreichten Stellen per Fehler-Vektor statt sie zu glauben. */
knplan(rad, target) = {
  /* exp-089: ab Rang 8 16 statt 30 Probepunkte -- jeder kostet dort eine
     Kette in H_{r-1} (b = 2, r = 8, dps 16: Planung 88 von 90 s) */
  my(pts = knarc(rad, if (knfast, 16, 30)), W, h, hmin, htop, hc2, g, reach, K, N, hc, xl, tp = getabstime());
  W = knwalk(pts);
  h = vector(#W, j, kns*imag(W[j]));
  hmin = vecmin(h); htop = vecmax(h);
  hc = 0.6*(htop - hmin); hc2 = hmin + 0.85*(htop - hmin);
  if (hc < 0.15, return([oo, 0, 0, rad, hmin, htop, rad, knc]));
  reach = rad;
  for (j = 1, #pts, if (h[j] < hc2,
    /* e_ALARM durchlassen: sonst schluckt dieses iferr das alarm(30) der
       Bogenplanung, und ein langsamer Kandidat laeuft ohne Grenze weiter
       (b = e, c = -0.875, rad 0.75: > 5 min) */
    g = iferr(if (real(pts[j]) > knc, knfi(pts[j]), knf(pts[j])), E, knc + knr, errname(E) != "e_ALARM");
    reach = max(reach, abs(g - knc))));
  if (reach >= 0.95*abs(knc - kp), return([oo, 0, 0, rad, hmin, htop, reach, knc]));
  /* langsame Bogen verwerfen: normal brauchen die 30 Probepunkte 15-50 ms.
     Mehr heisst lange ipent-Ketten (Punkte weit weg von p), und das Loesen
     mit ~10x so vielen Punkten lief dann bei b = e ueber 10 min */
  if (getabstime() - tp > 400*knslow, return([oo, 0, 0, rad, hmin, htop, reach, knc]));
  /* Bilder ausserhalb des Bogens: A wird dort nur extrapoliert, die hohen
     Koeffizienten sind aber nur bis eps*rad^-k bestimmt -- Verstaerkung
     (reach/rad)^K (Befund der Wissensgraph-Sitzung: alte e-Grenze bei 26 Stellen) */
  if (reach > knreachtol*rad, return([oo, 0, 0, rad, hmin, htop, reach, knc]));
  /* auf der reellen Achse muss eine volle Periode [x, H_{r-1}(x)] im Bogen liegen */
  xl = knc - 0.97*rad;
  if (real(knf(xl)) > knc + 0.97*rad, return([oo, 0, 0, rad, hmin, htop, reach, knc]));
  K = ceil((target + 5)/log(abs(knc - kp)/reach)*log(10));
  /* exp-082: wirksame Abklinghoehe von theta = (hcut - hmin) + 0.3. Die Ecken der
     reellen Kurve liegen 0.16-0.23 unter hmin (Wissensgraph-Sitzung), und die
     c_n fallen zusaetzlich. Zu viele Harmonische schaden: b = e, dps 45, gleicher
     Bogen: N = 80 Rest 9e-24, N = 49 2e-42, N = 35 7e-47, N = 25 8e-47. */
  N = ceil((target + 10)/(2*Pi*(hc + 0.3))*log(10)) + 2;
  return([K + 2*N, K, N, rad, hmin, htop, reach, knc]);
}

/* Ein Ausgleichsproblem auf dem Bogen vom Radius rad. JEDER Bogenpunkt
   bekommt eine Gleichung (theta nahe p, Funktionalgleichung nahe der Achse),
   und die Reihe ist bis zum weitesten Bildpunkt skaliert. Laesst man Punkte
   aus, ist die Loesung nicht mehr eindeutig: bei b = 2 erfuellte dann eine
   um 2e-21 falsche Reihe alle Gleichungen auf 1e-32.
   Rueckgabe [A(u), A'(u), theta-Koeff., rho, max. Restfehler, #theta, #FE]. */
knsolve(rad, K, N, hcut, hcut2) = {
  if (!knmode, return(knsolvelsq(rad, K, N, hcut, hcut2)));
  /* exp-095: scheitert die Wiederholung auf diesem Bogen, die Ausgleichsrechnung */
  return(iferr(knisolve(rad, K), E,
               if (kndebug, print("kn:   Wiederholung scheitert (", E, "), Ausgleichsrechnung"));
               knsolvelsq(rad, min(K, knmaxk), N, hcut, hcut2), errname(E) != "e_ALARM"));
}
knsolvelsq(rad, K, N, hcut, hcut2) = {
  my(pts, W, eqs, U, rows, scl, Mt, bv, cn, sol, A, cc, rho, tp, u1, u2, uk1, uk2, val, e, nx, nk, g, res, nth = 0, nfe = 0, sav, ts = getabstime());
  pts = knarc(rad, ceil(knptsf*(K + 1 + 2*(N + 1))));
  W = knwalk(pts);
  if (kndebug, print("kn:   solve K=", K, " N=", N, " Punkte ", #pts, ": alpha_p ", getabstime() - ts, " ms"));
  eqs = List();
  rho = rad;
  for (j = 1, #pts,
    if (kns*imag(W[j]) >= hcut, listput(eqs, [1, pts[j], W[j]]); nth++);
    if (kns*imag(W[j]) < hcut2,
      tp = getabstime();
      if (real(pts[j]) > knc, g = knfi(pts[j]); listput(eqs, [2, pts[j], g, 1]),
                              g = knf(pts[j]);  listput(eqs, [2, pts[j], g, -1]));
      if ((kndebug >= 2) && (getabstime() - tp > 1000),
        print("kn:     Bild langsam: ", getabstime() - tp, " ms bei ", precision(pts[j], 8)));
      rho = max(rho, abs(g - knc)); nfe++));
  if (kndebug, print("kn:   Bilder ", getabstime() - ts, " ms (", nth, " theta, ", nfe, " FE)"));
  U = K + 1 + 2*(N + 1);
  rows = 2*#eqs + 1;
  scl = vector(N + 1, n, vecmax(concat([1e-300], vector(#eqs, j, if (eqs[j][1] == 1, abs(knen(eqs[j][3], n-1)), 0)))));
  Mt = matrix(rows, U); bv = vector(rows)~;
  for (j = 1, #eqs,
    tp = eqs[j][1]; u1 = (eqs[j][2] - knc)/rho; uk1 = 1;
    if (tp == 1,
      for (k = 1, K + 1, Mt[2*j-1, k] = real(uk1); Mt[2*j, k] = imag(uk1); uk1 *= u1);
      for (n = 1, N + 1, e = knen(eqs[j][3], n-1)/scl[n];
        Mt[2*j-1, K+1+2*n-1] = -real(e); Mt[2*j, K+1+2*n-1] = -imag(e);
        Mt[2*j-1, K+1+2*n]   =  imag(e); Mt[2*j, K+1+2*n]   = -real(e));
      bv[2*j-1] = real(eqs[j][3]); bv[2*j] = imag(eqs[j][3]);
    ,
      u2 = (eqs[j][3] - knc)/rho; uk2 = 1;
      for (k = 1, K + 1, val = uk1 - uk2;
        Mt[2*j-1, k] = real(val); Mt[2*j, k] = imag(val); uk1 *= u1; uk2 *= u2);
      bv[2*j-1] = eqs[j][4]; bv[2*j] = 0;
    );
  );
  /* Normierung alpha(1) = 0: von 1 mit ipent in den Bogen, alpha(ipent^k 1) = -k.
     exp-090: bis 60 statt 6 Schritte -- der Bogen von Rang r liegt bei
     Re ~ -(r - 5) (b = 2, r = 12: um -6.99, die Bahn 1, 0, -1, ... braucht ~8) */
  nx = 1; nk = 0;
  while ((abs(nx - knc) > 0.9*rad) && (nk < 60), nx = knfi(nx); nk++);
  if ((abs(nx - knc) > 0.9*rad) || (abs(imag(nx)) > 1e-20), error("knhexinit: kein reeller Normierungspunkt im Bogen"));
  u1 = (real(nx) - knc)/rho; for (k = 1, K + 1, Mt[rows, k] = u1^(k-1)); bv[rows] = -nk;
  cn = vector(U, k, vecmax(abs(Mt[,k])) + 1e-300);
  Mt = matrix(rows, U, i, k, Mt[i,k]/cn[k]);
  Mt[rows,] = Mt[rows,]*1e6; bv[rows] *= 1e6;
  if (kndebug, print("kn:   Matrix ", rows, "x", U, " ", getabstime() - ts, " ms"));
  sav = default(realprecision);
  default(realprecision, ceil(knsolvep*precis));
  sol = knlsq(Mt, bv);
  default(realprecision, sav);
  if (kndebug, print("kn:   geloest ", getabstime() - ts, " ms, log10 R-Diagonale ", precision(log(knlsqcond)/log(10), 4)));
  sol = vector(U, k, precision(sol[k]/cn[k], precis));
  A = sum(k = 1, K + 1, sol[k]*x^(k-1));
  cc = vector(N + 1, n, (sol[K+1+2*n-1] + I*sol[K+1+2*n])/scl[n]);
  res = 0;
  for (j = 1, #eqs,
    tp = eqs[j][1];
    val = if (tp == 1, subst(A, x, (eqs[j][2] - knc)/rho) - eqs[j][3] - sum(n = 1, N + 1, cc[n]*knen(eqs[j][3], n-1)),
          subst(A, x, (eqs[j][2] - knc)/rho) - subst(A, x, (eqs[j][3] - knc)/rho) - eqs[j][4]);
    res = max(res, abs(val)));
  if (kndebug, print("kn:   Rest ", precision(res, 3), "  rho/rad ", precision(rho/rad, 5), "  ", getabstime() - ts, " ms"));
  return([A, deriv(A), cc, rho, res, nth, nfe, knc]);
}

/* exp-095: Kneser als Wiederholung, nach loop()/staylor()/thtaylor() in fatou.gp.
   A(u) = sum a_k u^k, u = (z - c)/R, reell. Ein Schritt:
     theta  aus A an Mt Punkten EINER Periode w = wre + j/Mt + i s hth im Kreis
            (z = alpha_p^-1(w)): t_j = A(z_j) - w_j, DFT, nur n >= 0 behalten;
     Werte  auf dem oberen Halbkreis |u| = 1: nahe der Achse A(g) +- 1 mit dem Bild
            g = ipent/pent z, sofern |g - c| <= kniqmax R; sonst W + theta(W),
            W = alpha_p(z);
     Glieder per DFT der Kreiswerte (unten konjugiert), dann A(1) = 0.
   Alles, was nicht von A abhaengt (Kreispunkte, alpha_p, Bilder, theta-Punkte,
   DFT-Matrizen), wird einmal aufgebaut; ein Schritt sind nur Polynomwerte und
   zwei Matrix-Vektor-Produkte. Die Abbildung ist affin, Anderson-Beschleunigung
   (Tiefe knimem, entspricht GMRES) statt einfacher Wiederholung: b = 2, dps 40,
   20 statt 101 Schritte. Die Reihe wird nach jedem Schritt auf precis Stellen
   gesetzt (wie ct in loop()): sonst speichert PARI die winzigen hohen Glieder mit
   weniger Woertern, und der Verlust steckt die naechsten Schritte an (gemessen
   57 -> 38 -> 19 Stellen, Boden 3e-22 statt 1e-51).
   Gemessen b = 2, dps 100, R = 0.65 Im p, K 500: 30 Schritte, 17.5 s (davon 12.8 s
   Aufbau) gegen 206 s fuer die Ausgleichsrechnung, beide ~110.5 Stellen. */
knainv(w) = {
  my(m = max(0, ceil((real(klog*w) - log(kur))/real(klog))), z);
  z = kp + subst(kpsi, x, exp(klog*(w - m)));
  for (i = 1, m, z = knf(z));
  return(z);
}
knimap(P, a) = {
  my(A = Polrev(a), S = P[3], v = vector(S), tj, thv, An);
  tj = vector(#P[10], j, subst(A, x, P[10][j]) - P[11][j]);
  thv = P[14]*(P[13]*tj~);
  for (s = 1, S, if (P[5][s] == 2, v[s] = subst(A, x, P[6][s]) + P[7][s]));
  for (i = 1, #P[8], v[P[8][i]] = P[9][i] + thv[i]);
  An = real(P[12]*v~)/S;
  An[1] -= subst(Polrev(An), x, P[15]) + P[16];
  return(precision(An, precis));
}
/* Urbilder der theta-Linie w = wre + j/Mt + i s hth; 0, wenn eines scheitert
   oder ausserhalb von |u| <= kniumax liegt */
knitline(c, R, wre, hth, Mt) = {
  my(z = vector(Mt));
  for (j = 1, Mt,
    z[j] = iferr(knainv(wre + (j - 1)/Mt + I*kns*hth), E, oo, errname(E) != "e_ALARM");
    if ((z[j] == oo) || (abs(z[j] - c)/R > kniumax), return(0)));
  return(z);
}
knisolve(rad, K) = {
  my(c = knc, R = rad, S, phi, zs, Ws, hs, hmin, htop, typ, ug, sg, g, q, nfe = 0, Wth = List(), ith = List(),
     hth, wre = 0, found = 0, tlo, thi, nt, us, mx, best, Mt, Nt, wj, zj, uj, nx, nk, E, eth, eW, P, a, gv, f, fn,
     Fs = List(), Gs = List(), dF, dG, gam, fbest = oo, nst = 0, tol, it, A, tj, cc, res, val, ts = getabstime(), tg,
     tgt = if (kntarget, kntarget, precis));
  S = max(32, ceil(knis*(K + 1)));
  phi = vector(S, s, Pi*(s - 1/2)/S);
  zs = vector(S, s, c + R*exp(I*phi[s]));
  Ws = knwalk(zs);
  hs = vector(S, s, kns*imag(Ws[s])); hmin = vecmin(hs);
  typ = vector(S); ug = vector(S); sg = vector(S);
  for (s = 1, S,
    g = iferr(if (real(zs[s]) > c, knfi(zs[s]), knf(zs[s])), E, oo, errname(E) != "e_ALARM");
    q = if (g == oo, oo, abs(g - c)/R);
    if (q <= kniqmax, typ[s] = 2; ug[s] = (g - c)/R; sg[s] = if (real(zs[s]) > c, 1, -1); nfe++,
                      typ[s] = 1; listput(Wth, Ws[s]); listput(ith, s)));
  if ((#Wth == 0) || (nfe == 0), error("knisolve: Kreis ohne theta- oder FE-Teil (", #Wth, "/", nfe, ")"));
  htop = vecmin(vector(#ith, i, hs[ith[i]]));
  tg = getabstime() - ts;
  /* theta-Linie: eine Periode auf der Hoehe hth knapp unter dem theta-Teil des
     Kreises, deren Urbilder im Kreis liegen (|u| <= kniumax); das Fenster mit
     dem kleinsten groessten |u| */
  tlo = vecmin(vector(#Wth, i, real(Wth[i]))) - 1; thi = vecmax(vector(#Wth, i, real(Wth[i]))) + 1;
  nt = ceil(16*(thi - tlo)) + 1;
  hth = htop - 0.1;
  /* je Hoehe zuerst die Periode um die Mitte des theta-Teils, sonst das beste
     Fenster aus einem 1/16-Raster; angenommen wird erst, wenn ALLE Mt Urbilder
     gelingen und im Kreis liegen (b = 2, Rang 8: ein einzelnes H_7 an der
     Linie scheiterte zwischen den Probepunkten und brach den Aufbau ab) */
  while ((!found) && (hth > hmin + 0.3),
    /* theta faellt wie e^(-2 pi n (h - hs)), hs ~ hmin - 0.25 (Ecken der reellen Kurve, exp-082) */
    Mt = 8*ceil(2*ceil((tgt + 10)*log(10)/(2*Pi*(hth - hmin + 0.25)))/8);
    wre = (tlo + thi)/2 - 1/2;
    zj = knitline(c, R, wre, hth, Mt);
    if (zj, found = 1,
      us = vector(nt, i, iferr(abs(knainv(tlo + (i - 1)/16 + I*kns*hth) - c)/R, E, oo, errname(E) != "e_ALARM"));
      best = oo;
      for (i = 1, nt - 16, mx = vecmax(us[i..i+16]); if (mx < best, best = mx; wre = tlo + (i - 1)/16));
      if (best <= kniumax, zj = knitline(c, R, wre, hth, Mt); if (zj, found = 1)));
    if (!found, hth -= 0.1));
  if (!found, error("knisolve: keine theta-Periode im Kreis"));
  best = vecmax(vector(Mt, j, abs(zj[j] - c)/R));
  Nt = Mt/2 - 1;
  wj = vector(Mt, j, wre + (j - 1)/Mt + I*kns*hth);
  uj = vector(Mt, j, (zj[j] - c)/R);
  nx = 1; nk = 0;
  while ((abs(nx - c) > 0.9*R) && (nk < 60), nx = knfi(nx); nk++);
  if ((abs(nx - c) > 0.9*R) || (abs(imag(nx)) > 1e-20), error("knisolve: kein reeller Normierungspunkt im Kreis"));
  /* Potenzen per Verdopplung (geopow, Fehler ~log K) statt K*S exp() -- b = 2,
     dps 100, K 484: 8 s des Aufbaus */
  E = vector(S, s, geopow(exp(-I*phi[s]), K));
  E = matrix(K + 1, S, k, s, E[s][k]);
  eth = matrix(Nt + 1, Mt, n, j, exp(-2*Pi*I*kns*(n - 1)*((j - 1)/Mt + wj[1]))/Mt);
  eW = matrix(#Wth, Nt + 1, i, n, exp(2*Pi*I*kns*(n - 1)*Wth[i]));
  P = [c, R, S, K, typ, ug, sg, Vec(ith), Vec(Wth), uj, wj, E, eth, eW, (nx - c)/R, nk];
  if (kndebug, print("kn:   iter K=", K, " S=", S, ": ", nfe, " FE, ", #Wth, " theta, theta-Linie h ", precision(hth, 4),
                     " (Kreis ", precision(hmin, 4), "..", precision(htop, 4), "), |u| <= ", precision(best, 3), ", Mt ", Mt,
                     "  Aufbau ", getabstime() - ts, " ms (Kreis ", tg, ")"));
  /* Anderson-Wiederholung ab A = 0 */
  a = vectorv(K + 1); tol = 10.^(1 - precis);
  for (i = 1, knimaxit,
    it = i;
    gv = knimap(P, a); f = gv - a; fn = vecmax(apply(abs, f));
    if (fn < tol, a = gv; break);
    if (2*fn < fbest, fbest = fn; nst = 0, nst++);
    if (nst >= 6, a = gv; break);
    listput(Gs, gv); listput(Fs, f);
    if (#Fs > knimem + 1, listpop(Gs, 1); listpop(Fs, 1));
    if (#Fs > 1,
      dF = matrix(#f, #Fs - 1, i1, j1, Fs[j1+1][i1] - Fs[j1][i1]);
      dG = matrix(#f, #Fs - 1, i1, j1, Gs[j1+1][i1] - Gs[j1][i1]);
      gam = knlsq(dF, f, 2);
      a = precision(gv - dG*gam, precis),
      a = gv));
  knistat = [it, fn];
  A = Polrev(a);
  tj = vector(Mt, j, subst(A, x, uj[j]) - wj[j]);
  cc = (eth*tj~)~;
  /* Rest: die Gleichungen auf dem Kreis mit dem fertigen A */
  res = 0;
  for (s = 1, S,
    val = if (typ[s] == 2, subst(A, x, exp(I*phi[s])) - subst(A, x, ug[s]) - sg[s],
              subst(A, x, exp(I*phi[s])) - Ws[s] - sum(n = 1, Nt + 1, cc[n]*knen(Ws[s], n - 1)));
    res = max(res, abs(val)));
  if (kndebug, print("kn:   iter ", it, " Schritte, |G(a) - a| ", precision(fn, 3), ", Rest ", precision(res, 3), "  ", getabstime() - ts, " ms"));
  return([A, deriv(A), cc, R, res, #Wth, nfe, c]);
}

knhexinit(r) = {
  my(eps = 10.0^(-precis-8), target = precis - 3, z, z2, res = 0, sr, R, sc, ps, W, pl, best,
     pl2, s1, s2, tst, d, hf, t0, kd, tq, tslow, nafter);
  /* exp-089: Rang r >= 8 baut auf H_{r-1} (knf/knfi); ohne Angabe Rang 6 auf
     der Pentation. Ziel nach den Stellen der unteren Stufe. */
  if (!r, r = 6);
  knrank = r; knslow = 1;
  /* exp-097: Atlas ab Rang r gehoert zur alten Loesung; waehrend der Planung keine
     Kreise bauen (sie liefe unter alarm) */
  hyatclear(r); hyatfreeze = 1;
  my(kmax = if (knmode, knmaxki, knmaxk));
  /* exp-095: die schnelle Planung (16 Probepunkte, Mittelpunkte von rechts, frueher
     Schluss) auch fuer Rang 6 in der Wiederholung ab precis > 100 -- b = 2, dps 200:
     Planung 264 von 1033 s, dps 100: 35 von 94 s; gewaehlt wurde der rechteste
     Mittelpunkt */
  knfast = (r > 6) || (knmode && (precis > 100));
  if (r > 6,
    target = max(8, min(precis - 3, floor(hydig[r-1]) - 3));
    eps = 10.0^(-target - 8));
  if (r == 6, hexreg = 0); knRmeas = 0;
  kntarget = target;
  if (imag(sexp(1)) != 0, error("knhexinit: nur reelle Basen"));
  t0 = getabstime();
  /* komplexer Fixpunkt der Pentation: fuer ipent anziehend */
  z = 0.5;
  for (i = 1, 2000, z2 = knfi(z); if (abs(z2 - z) < 1e-15, z = z2; break); z = z2);
  if (imag(z) < 0, z = conj(z));
  if (imag(z) < 1e-3, error("knhexinit: kein komplexer Fixpunkt der Pentation"));
  if (kndebug, print("kn: p0 ", precision(z, 8), "  ", getabstime() - t0, " ms"));
  sr = 1/4;
  /* exp-084: bis 1/2048 statt 1/128 -- die Reihe um p hat bei b = 50 den
     Radius 0.013, bei b = 100 nur 0.0075 (geht ab sr = 1/256) */
  for (k = 1, 10,
    res = iferr(pentfixser(t->knf(t), z, sr, 0, eps, "knhexinit"), E, 0);
    if (type(res) == "t_VEC", break);
    sr = sr/2);
  if (type(res) != "t_VEC", error("knhexinit: Pentation um p nicht entwickelbar"));
  kp = res[1]; sc = res[2]; R = res[3];
  if (kndebug, print("kn: Reihe um p, R = ", precision(R, 5), "  ", getabstime() - t0, " ms"));
  klam = polcoeff(sc, 1);
  if (abs(klam) <= 1, error("knhexinit: p ist nicht abstossend, |P'(p)| = ", precision(abs(klam),6)));
  ps = pentpsi(sc, R/4, eps, poldegree(sc) - 1);
  if (!ps[2], error("knhexinit: Poincare-Reihe an p konvergiert nicht"));
  krs = R/5;
  ksig = Pol(serreverse(ps[1] + O(x^(poldegree(ps[1])+1))));
  if (abs(polcoeff(ksig, poldegree(ksig)))*krs^poldegree(ksig) > eps,
    error("knhexinit: Schroeder-Reihe an p konvergiert bei r = ", precision(krs,5), " nicht"));
  ksig = pentclean(ksig, krs, eps, "sigma_p");
  /* exp-095: Poincare-Funktion psi = sigma^-1 fuer alpha_p^-1 (knainv), gueltig bis |u| = R/4 */
  kpsi = Pol(ps[1]); kur = R/4;
  klog = log(klam); kper = 2*Pi*I/klog;
  if (kndebug, print("kn: sigma_p fertig  ", getabstime() - t0, " ms"));
  knc = real(kp); knr = imag(kp);
  /* Abklingrichtung von theta: s*Im alpha_p waechst zu p hin */
  W = knwalk([knc + 0.6*knr*exp(0.05*I), knc + 0.6*knr*I]);
  kns = if (imag(W[2]) > imag(W[1]), 1, -1);
  /* exp-089: die Zeitgrenzen der Planung (alarm, 400 ms je Kandidat) sind fuer
     Rang 6 gemessen, wo ein pent ~0.5 ms kostet; hoehere Raenge skalieren sie
     mit den gemessenen Kosten der unteren Stufe */
  /* exp-095: auch fuer Rang 6 messen -- bei dps 200 kostet ein pent ein
     Vielfaches, alle Kandidaten liefen in alarm/400 ms ("kein brauchbares
     Bogenpaar"); und schon bei dps 100 hing die Wahl an der Last: unter Last
     fielen die rechten Mittelpunkte (lange ipent-Ketten) durch die 400 ms, und
     b = 2 bekam c = Re p + 0.5 Im p mit K 526 statt Re p + 1.5 Im p. Das
     Verhaeltnis der Messung gleicht die Last aus; bei dps 45 bleibt knslow ~1. */
  if (1,
    tq = getabstime();
    for (j = 1, 6,
      iferr(knf(knc + 0.8*knr*exp(I*Pi*j/7)), E, 0);
      iferr(knfi(knc + 0.8*knr*exp(I*Pi*j/7)), E, 0));
    knslow = max(1, (getabstime() - tq)/12/0.5);
    if (kndebug, print("kn: Rang ", r, ": ", precision((getabstime() - tq)/12., 3), " ms je Wert der unteren Stufe -> Zeitgrenzen x", precision(knslow, 3))));
  /* Bogenradius: kleinster geschaetzter Aufwand */
  /* exp-081: auch der Mittelpunkt ist frei. Re p ist nur dann gut, wenn p
     weit ueber der Achse liegt; bei b = 10 (Im p = 0.87) muss der Kreis um
     Re p bis 0.85 Im p reichen, um eine Periode der reellen Achse samt Bildern
     zu fassen (K = 740). Weiter rechts ist p weiter weg, die Periode gleich
     lang: Verhaeltnis Reichweite/|c - p| 0.88 -> ~0.5 (gemessen). */
  /* Geplant werden PAARE benachbarter Radien um denselben Mittelpunkt: der
     Fehler-Vektor braucht einen zweiten Bogen, und bei b = e lag der beste
     Einzelbogen isoliert zwischen unbrauchbaren Nachbarn. Kosten des Paars =
     die des teureren Bogens; K und N gelten fuer beide. */
  my(rads = [0.45, 0.55, 0.65, 0.75, 0.85, 1.0, 1.15, 1.3], P, c2);
  best = [oo]; pl2 = [oo];
  /* erst streng (Bilder <= 1.02 rad), sonst einmal mit 1.15 rad: b = 50 hat
     streng keinen Bogen; der Fehler-Vektor zeigt dann, was es kostet */
  /* exp-087: der gelockerte Durchgang (nur grosse Basen) sucht auch bei
     1.45 und 1.6 Im p und paart gleiche Radien benachbarter Mittelpunkte;
     der zweite Bogen wird um seinen eigenen Mittelpunkt (pl2[8]) gerechnet.
     b = 100, dps 20, K 300: Versatz 1.5 Radius 1.45 Rest 3e-31, der alte
     Hauptbogen (Versatz 1.0, Radius 1.3) 1e-13 -- er lag ausserhalb der Liste
     und hatte keinen Nachbarradius. */
  my(Pprev, prs = List());
  for (pass = 1, 2,
  if ((pass == 2) && (best[1] != oo), break);
  knreachtol = if (pass == 1, 1.02, 1.15);
  if (pass == 2, rads = concat(rads, [1.45, 1.6]));
  Pprev = 0; nafter = -1;
  for (ic = 1, 7,
    /* exp-090: ab Rang 8 die Mittelpunkte von rechts nach links und Schluss
       einen Mittelpunkt nach dem ersten mit brauchbarem Paar -- gewaehlt wurde
       bisher immer Versatz 1.5 oder 1.25 (Raenge 6, 8, 10, 12 bei b = 2), und
       die Planung trug 74-87 % aller Aufrufe der unteren Stufen */
    if (knfast && (nafter >= 1), break);
    knc = real(kp) + [0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5][if (knfast, 8 - ic, ic)]*knr;
    P = vector(#rads); tslow = 0;
    for (j = 1, #rads,
      /* exp-089: ab Rang 8 nach dem ersten zu langsamen Radius den Rest dieses
         Mittelpunkts auslassen -- groessere Radien haben laengere Ketten (b = 2,
         r = 8, dps 16: 11 Kandidaten liefen in alarm, 31 von 43 s) */
      if (tslow && knfast, P[j] = [oo]; next);
      /* alarm: ein haengender Kandidat kostet hoechstens 1 s, ein normaler
         braucht 15-50 ms (alarm gibt dabei ein t_ERROR zurueck statt zu werfen) */
      if (kndebug, print1("kn: plan c = ", precision(knc, 5), " rad = ", rads[j], " ... "));
      tq = getabstime();
      pl = iferr(alarm(ceil(knslow), knplan(rads[j]*knr, target)), E, [oo]);
      if (kndebug, print(if (type(pl) == "t_ERROR", "alarm", pl[1]), "  ", getabstime() - t0, " ms"));
      if ((type(pl) == "t_ERROR") || (getabstime() - tq > 400*knslow), tslow = 1);
      if (type(pl) == "t_ERROR", pl = [oo]);
      P[j] = pl);
    for (j = 1, #rads - 1,
      if ((P[j][1] != oo) && (P[j+1][1] != oo),
        c2 = max(P[j][1], P[j+1][1]);
        listput(prs, [c2, P[j], P[j+1]]);
        if ((best[1] == oo) || (c2 < max(best[1], pl2[1])),
          best = P[j]; pl2 = P[j+1])));
    if ((pass == 2) && (Pprev != 0),
      for (j = 1, #rads,
        if ((P[j][1] != oo) && (Pprev[j][1] != oo),
          c2 = max(P[j][1], Pprev[j][1]);
          listput(prs, [c2, P[j], Pprev[j]]);
          if ((best[1] == oo) || (c2 < max(best[1], pl2[1])),
            best = P[j]; pl2 = Pprev[j]))));
    if (nafter >= 0, nafter++, if (best[1] != oo, nafter = 0));
    Pprev = P));
  knreachtol = 1.02; hyatfreeze = 0;
  /* knforce = [Bogenanteil, theta-Schwellenanteil, K, N, Mittelpunkt-Versatz] (0 = automatisch), nur fuer Experimente */
  if (type(knforce) == "t_VEC" && knforce[1],
    knc = real(kp) + if (#knforce >= 5, knforce[5], 0)*knr;
    best = knplan(knforce[1]*knr, target);
    pl2 = knplan((knforce[1] + 0.05)*knr, target);
    prs = List([[0, best, pl2]]));
  if ((best[1] == oo) || (pl2[1] == oo), error("knhexinit: kein brauchbares Bogenpaar"));
  /* exp-091: alle brauchbaren Paare der Kosten nach; scheitert die Vorab-Loesung
     auf einem (b = 1.85 und 1.9, Rang 8: ein Bildpunkt verlangt die Hexation
     ausserhalb ihres Kreises, "Newton bleibt nicht im Bogen"), das naechste --
     hoechstens vier */
  /* zuerst genau das bisher gewaehlte Paar (Rang 6 unveraendert), dann die
     uebrigen nach Kosten */
  prs = concat([[0, best, pl2]], select(e -> (e[2] != best) || (e[3] != pl2), vecsort(Vec(prs), 1)));
  my(K0, s0 = 0, Kn, ft, Kneed, hc1, hc2, tsw, s3, pl3, c0, cm, rr);
  for (ip = 1, min(#prs, 4),
    best = prs[ip][2]; pl2 = prs[ip][3];
    knc = best[8];
    best[2] = max(best[2], pl2[2]); best[3] = max(best[3], pl2[3]);
    if (type(knforce) == "t_VEC" && knforce[3], best[2] = knforce[3]);
    if (type(knforce) == "t_VEC" && knforce[4], best[3] = knforce[4]);
    /* b = 10 (alter Aufbau um Re p): K = 740, 4.3 min bei \p 20 fuer 13.5
       Stellen -- ueber knmaxk nimmt hexinit Sheldons schnelle Naeherung */
    /* exp-084: der Plan darf bis 2 knmaxk gehen -- er schaetzt K aus |knc - p|
       und liegt oft doppelt zu hoch (b = 28: Plan 575, K0 = 287 genuegte); die
       Vorab-Loesung mit K0 = Plan/2 <= knmaxk misst den wirklichen Bedarf
       (b = 100, dps 20: Plan 1094) */
    if (best[2] > 2*kmax,
      if (ip == 1, error("knhexinit: zu teuer (K = ", best[2], " > 2 knmaxk = ", 2*kmax, ")"));
      next);
    hf =if (type(knforce) == "t_VEC" && knforce[2], knforce[2], 0.6);
    if (kndebug, print("kn: Plan K = ", best[2], " N = ", best[3], " rad = ", precision(best[4]/knr, 3), " c = ", precision(best[8], 5), "  ", getabstime() - t0, " ms, Aufrufe ", hyncalls()));
    /* exp-083: K aus einer kleinen Vorab-Loesung (halbes K, ~1/8 der Kosten)
       statt aus einer vollen ersten Loesung, die bei b = 2 und e ohnehin noch
       einmal mit groesserem K gerechnet werden musste. Der Planer nimmt |knc - p|
       als Konvergenzradius; alpha hat naehere Singularitaeten (b = e um
       knc = -0.18 bei 0.33 +- 1.37 i, Abstand 1.46 statt ~2.5). Gemessen wird
       das Abfallgesetz log|a_k| ~ a + b k + c log k (knfit), K mit 10 % Reserve. */
    K0 = max(60, best[2]\2); Kn = best[2];
    hc1 = best[5] + hf*(best[6] - best[5]); hc2 = best[5] + max(0.85, hf + 0.1)*(best[6] - best[5]);
    s0 = iferr(knsolve(best[4], K0, best[3], hc1, hc2),
               E, if (kndebug, print("kn: Paar ", ip, " scheitert: ", E)); 0, errname(E) != "e_ALARM");
    if (type(s0) == "t_VEC", break));
  if (type(s0) != "t_VEC", error("knhexinit: kein Bogenpaar laesst sich loesen (", min(#prs, 4), " versucht)"));
  if (kndebug, print("kn: Vorab-Loesung K0 = ", K0, "  ", getabstime() - t0, " ms, Aufrufe ", hyncalls()));
  ft = knfit(s0[1], K0);
  if (ft[2] < 0,
    Kneed = ceil(1.1*knkfor(ft, -(target + 5)*log(10), K0));
    knRmeas = s0[4]*exp(-ft[2]);
    if (kndebug, print("kn: Fit b = ", precision(ft[2], 5), " c = ", precision(ft[3], 4), " -> K ", Kneed, " (Plan ", best[2], ")"));
    Kn = max(Kneed, K0),
    /* kein Abfall auf [K0/4, 3K0/4]: die Glieder liegen dort schon am Boden,
       mehr Glieder bringen nichts -- K0 reicht (b = 28, dps 16: Plan K = 575,
       Vorab K0 = 287 Rest 9.7e-14, volle Loesungen 1.2e-13 und 1.7e-13 in je
       40 s). Ob das stimmt, misst danach der Fehler-Vektor. */
    Kn = K0);
  /* der Plan lag unter knmaxk, die Vorab-Loesung lief: dann mit knmaxk
     weiterrechnen statt auf Sheldons Naeherung zu fallen (b = 50, dps 20:
     Fit K = 603; die Naeherung gibt dort 0.2 Stellen) */
  Kn = min(Kn, max(K0, kmax));
  best[2] = Kn;
  s1 = if (Kn == K0, s0, knsolve(best[4], Kn, best[3], hc1, hc2));
  /* grosses K kann die Ausgleichsrechnung verderben: Bogenpunkte bei
     |u| = rad/rho < 1, Bilder bis |u| = 1, Kondition ~ (rho/rad)^K (b = 50,
     dps 20, rho/rad = 1.19: K0 = 221 Rest 6e-32, K = 600 Rest 1e-13) */
  if (s1[5] > 100*s0[5],
    if (kndebug, print("kn: K = ", Kn, " Rest ", precision(s1[5], 3), " schlechter als K0 = ", K0, " -> K0"));
    s1 = s0; Kn = K0; best[2] = K0);
  if (kndebug, print("kn: erste Loesung  ", getabstime() - t0, " ms, Aufrufe ", hyncalls()));
  /* exp-095: in der Wiederholung zuerst den Geschwisterkreis (gleicher Radius,
     Mittelpunkt + 0.125 Im p) als zweite Loesung -- der groessere Nachbarradius
     braucht mehr Glieder (b = 2, dps 100 und 200: Rest 1.5e-98 bzw. 9e-179,
     danach rechnete doch der Geschwisterkreis) */
  if (knmode,
    knc = best[8] + 0.125*knr;
    pl3 = iferr(alarm(2*ceil(knslow), knplan(best[4], target)), E, [oo]);
    if ((type(pl3) != "t_ERROR") && (pl3[1] != oo), pl2 = pl3);
    knc = best[8]);
  knc = pl2[8];
  /* exp-089: scheitert der zweite Bogen (ab Rang 8: Punkte nahe dem komplexen
     Fixpunkt der Stufe darunter lassen sich dort nicht umkehren), gilt er als
     unbrauchbar -- Rest oo, die Ausweichwege unten versuchen andere Boegen */
  s2 = iferr(knsolve(pl2[4], Kn, best[3], pl2[5] + hf*(pl2[6] - pl2[5]), pl2[5] + max(0.85, hf + 0.1)*(pl2[6] - pl2[5])),
             E, if (kndebug, print("kn: zweiter Bogen scheitert: ", E)); [0, 0, [0], oo, oo, 0, 0, pl2[8]],
             errname(E) != "e_ALARM");
  knc = best[8];
  /* exp-084: die Loesung mit dem kleineren Rest wird Hauptloesung, der
     Fehler-Vektor nahm sonst den schlechteren Rest (b = 10, dps 45:
     Hauptbogen 1.1e-46, Nachbarbogen 1.2e-53). Gleicher Mittelpunkt, K, N. */
  if (s2[5] < s1[5],
    tsw = s1; s1 = s2; s2 = tsw;
    pl2[2] = best[2]; pl2[3] = best[3]; tsw = best; best = pl2; pl2 = tsw;
    if (kndebug, print("kn: Nachbarbogen rad = ", precision(best[4]/knr, 3), " wird Hauptloesung")));
  /* exp-084: liegen Bilder ausserhalb eines Bogens, waechst sein Rest wie
     Boden * (rho/rad)^K (b = 10, dps 45, Bogen 1.0: rho/rad = 1.06, K 352
     Rest 1e-46, K 460 7e-44; Bogen 1.15: rho/rad = 1, 1e-53). Der Planer
     sieht das mit 30 Punkten nicht (Reichweite 1.0014). Ist der zweite Bogen
     viel schlechter, ihn mit halbem K nachrechnen (1/8 der Kosten) und den
     besseren nehmen. Ein dritter Bogen weiter aussen half nicht: er braucht
     mehr K (b = 10, rad 1.32, K 352: Rest 4e-32). */
  /* exp-085: zuerst einen gleich grossen Bogen um einen leicht
     verschobenen Mittelpunkt rechnen -- eine unabhaengige zweite Loesung,
     deren Bilder innen bleiben (b = 10, dps 45, Hauptbogen 1.15 um
     Re p + 1.25 Im p: Versatz +0.125 Im p Rest 1e-53, +0.25 1e-53,
     -0.125 2e-44, -0.25 1e-34). Hoechstens drei Versuche. */
  c0 = s1[8];
  for (io = 1, 3,
    if (s2[5] <= 100*s1[5], break);
    /* exp-095: +0.125 hat die Wiederholung schon als zweiten Kreis gerechnet */
    knc = c0 + if (knmode, [0.25, -0.125, -0.25], [0.125, 0.25, -0.125])[io]*knr;
    pl3 = iferr(alarm(2*ceil(knslow), knplan(best[4], target)), E, [oo]);
    if (type(pl3) == "t_ERROR" || pl3[1] == oo, next);
    s3 = iferr(knsolve(pl3[4], best[2], best[3], pl3[5] + hf*(pl3[6] - pl3[5]), pl3[5] + max(0.85, hf + 0.1)*(pl3[6] - pl3[5])),
               E, [0, 0, [0], oo, oo, 0, 0, knc], errname(E) != "e_ALARM");
    if (kndebug, print("kn: Bogen um c ", precision(knc, 6), ": Rest ", precision(s3[5], 3), " statt ", precision(s2[5], 3)));
    if (s3[5] < s2[5], s2 = s3; pl3[2] = best[2]; pl3[3] = best[3]; pl2 = pl3));
  knc = c0;
  /* exp-087: Geschwisterbogen vor eigenem/halbem K -- b = 100, dps 20: K 600
     auf dem zweiten Bogen kostete 64 s fuer 7e-29, der Geschwisterbogen danach
     gab 2e-32. */
  /* exp-086: das gilt nur, wenn die Bilder deutlich ausserhalb liegen. Liegen
     sie fast innen (rho/rad <= 1.05), fehlen dem zweiten Bogen Glieder: er
     bekommt sein eigenes K aus seinem eigenen Abfall (b = 50, dps 20: Bogen
     1.3, rho/rad 1.015, K 221 Rest 6e-23, K 330 1e-32; Bogen 1.15, rho/rad
     1.18, K 221 6e-32, K 330 5e-31). Danach noch einmal tauschen. */
  if (s2[5] > 100*s1[5],
    tsw = 0;
    if (s2[4] <= 1.05*pl2[4],
      ft = knfit(s2[1], poldegree(s2[1]));
      if (ft[2] < 0, tsw = min(ceil(1.1*knkfor(ft, -(target + 5)*log(10), Kn)), max(Kn, kmax)));
      if (tsw <= Kn, tsw = 0),
      tsw = max(60, Kn\2));
    if (tsw,
      knc = pl2[8];
      s3 = iferr(knsolve(pl2[4], tsw, best[3], pl2[5] + hf*(pl2[6] - pl2[5]), pl2[5] + max(0.85, hf + 0.1)*(pl2[6] - pl2[5])),
                 E, [0, 0, [0], oo, oo, 0, 0, pl2[8]], errname(E) != "e_ALARM");
      knc = best[8];
      if (kndebug, print("kn: zweiter Bogen mit K = ", tsw, ": Rest ", precision(s3[5], 3), " statt ", precision(s2[5], 3)));
      if (s3[5] < s2[5], s2 = s3)));
  if (s2[5] < s1[5],
    tsw = s1; s1 = s2; s2 = tsw;
    tsw = best; best = pl2; pl2 = tsw;
    if (kndebug, print("kn: zweiter Bogen rad = ", precision(best[4]/knr, 3), " wird Hauptloesung")));
  if (s2[5] == oo, error("knhexinit: kein brauchbarer zweiter Bogen fuer den Fehler-Vektor"));
  best[2] = poldegree(s1[1]); best[3] = #s1[3] - 1;
  knc = s1[8];
  knA = s1[1]; knAd = s1[2]; kncc = s1[3]; knrho = s1[4]; knres = s1[5];
  kngood = 0.97*best[4];
  /* Fehler-Vektor: beide Loesungen an Punkten im gemeinsamen Gebiet; jede
     wird um ihren eigenen Mittelpunkt (s[8]) ausgewertet */
  d = 0;
  cm = (s1[8] + s2[8])/2; rr = 0.9*min(best[4], pl2[4]) - abs(s1[8] - s2[8]);
  tst = concat(vector(5, j, cm + (j - 3)/2*rr), [cm + 0.5*rr*exp(Pi*I/3), cm + 0.5*rr*exp(2*Pi*I/3)]);
  for (j = 1, #tst, d = max(d, abs(subst(s1[1], x, (tst[j] - s1[8])/s1[4]) - subst(s2[1], x, (tst[j] - s2[8])/s2[4]))));
  if (kndebug, print("kn: Fehler-Vektor |s1 - s2| = ", precision(d, 3), ", Rest s1 = ", precision(s1[5], 3), ", s2 = ", precision(s2[5], 3), ", Aufrufe ", hyncalls()));
  d = max(d, knres);
  kd = if (d == 0, precis, min(precis, precision(-log(d)/log(10), 9)));
  if (kd < 6, error("knhexinit: nur ", precision(kd,3), " Stellen -- Fehler-Vektor zu gross"));
  knt0 = knser(knc);
  /* reelles Fenster fuer knhex: A steigt darauf und ueberstreicht mehr als
     eine Periode (b = e: -3.624 .. -2.591) */
  knxl = knc - kngood; knxr = knc + kngood;
  kntl = real(knser(knxl)); kntr = real(knser(knxr));
  if (kntr - kntl <= 1, error("knhexinit: Fenster deckt keine volle Periode (", precision(kntr - kntl, 5), ")"));
  knpar = [best[2], best[3], best[4]/knr, best[5], best[6], s1[6] + s1[7], getabstime() - t0, s1[6], s1[7], pl2[4]/knr, best[8], knRmeas];
  knS[r] = knpack(r); knSd[r] = kd;
  if (r == 6, hexerr = kd; hexreg = 2);
  if (quietmode == 0,
    print("Rang ", r, " (Kneser) am komplexen Fixpunkt von H_", r-1, ": p = ", precision(kp, 20));
    print("  H_", r-1, "'(p) = ", precision(klam, 15), "   K = ", best[2], "  N = ", best[3], "  Bogen ", precision(best[4]/knr,3), " Im p");
    print("  ", kd, " Stellen (Fehler-Vektor zweier Boegen)"));
  return(1);
}

/* exp-089: Auswertung aus einem Zustand je Rang, damit hoehere gerade
   Raenge (Kneser auf H_{r-1}) die Hexation nicht ueberschreiben.
   S = [A, A', Mittelpunkt, rho, Fensterradius, xl, xr, tl, tr, Rang]. */
knpack(r) = [knA, knAd, knc, knrho, kngood, knxl, knxr, kntl, kntr, r];
knserS(S, z) = subst(S[1], x, (z - S[3])/S[4]);
knserdS(S, z) = subst(S[2], x, (z - S[3])/S[4])/S[4];
/* alpha(y): mit der Umkehrung (rechts) bzw. H_{r-1} (links) in den Bogen, dann A */
knihexS(S, y) = {
  my(k = 0, n = 0, r1 = S[10] - 1);
  /* exp-089: alle Schritte zaehlen, nicht nur k -- nahe dem komplexen
     Fixpunkt springt y zwischen links und rechts hin und her, ohne das Fenster
     zu treffen, und k bleibt bei 0 (b = 2, Rang 8 bei -4.086 + 0.707i: hing
     beim Aufbau von Rang 10 ohne Ende) */
  while (abs(y - S[3]) > S[5],
    if (real(y) > S[3], y = hyloinv(r1, y); k++, y = hylo(r1, y); k--);
    n++;
    if (n > 200, error("knihex: Punkt laesst sich nicht in den Bogen abbilden (Rang ", S[10], ")")));
  return(knserS(S, y) + k);
}
/* H_r(z): k so, dass Re(z) - k im reellen Fenster [tl, tr] liegt; dort
   reelle Loesung per Newton mit Bisektions-Rueckfall (A steigt), fuer
   komplexes z gedaempfter Newton ab dieser Loesung; dann H_{r-1}^k bzw. ^-k.
   Ein ungesicherter Newton aus der Fenstermitte sprang bei b = e schon im
   ersten Schritt ueber den Rand, obwohl die Loesung innen lag. */
knhexS(S, z) = {
  my(k, t, tr, a, b, y, fy, dy, yn, tol = 10.0^(-precis+2), r1 = S[10] - 1);
  k = ceil(real(z) - S[9]);
  t = z - k; tr = real(t);
  a = S[6]; b = S[7]; y = S[3];
  for (it = 1, 400,
    fy = real(knserS(S, y)) - tr;
    if (fy > 0, b = y, a = y);
    yn = y - fy/real(knserdS(S, y));
    if ((yn <= a) || (yn >= b), yn = (a + b)/2);
    dy = yn - y; y = yn;
    if (abs(dy) < tol, break));
  if (imag(t) != 0,
    for (it = 1, 200,
      dy = (knserS(S, y) - t)/knserdS(S, y);
      yn = y - dy;
      while (abs(yn - S[3]) > S[5], dy = dy/2; yn = y - dy; if (abs(dy) < tol, error("knhex: Newton bleibt nicht im Bogen")));
      y = yn;
      if (abs(dy) < tol, break)));
  if (k > 0, for (i = 1, k, y = hylo(r1, y)), for (i = 1, -k, y = hyloinv(r1, y)));
  return(y);
}
knihex(y) = knihexS(knS[6], y);
knhex(z) = knhexS(knS[6], z);


/* ecalle routine is used for base eta cheta and sexpeta functions */
/* works for other cases, for fz or fz^n with parabolic multiplier */
ecallelog=0;
ecalle(fz,n) = {
  local(i,z,ns,rem,logfzx);
  kecalle=0;
  ecallelog=0;
  kfz=fz; /* save fz */
  logfzx=Ser(log(fz/x));
  negterms=1;
  while (abs(polcoeff(fz,negterms+1))<1E-35,negterms++);
/*print("terms with negative coeffients= "negterms);*/
  for (i=-negterms,n,
    if (i==0, ecallelog=acoeff, kecalle=kecalle+acoeff*x^i);
    rem = Ser(subst(kecalle,x,fz) - kecalle + ecallelog*logfzx - 1);
    z=polcoeff(rem,i+negterms);
    z=subst(z,acoeff,x);
    ns=-polcoeff(z,0)/polcoeff(z,1);
    kecalle=subst(kecalle,acoeff,ns);
    ecallelog=subst(ecallelog,acoeff,ns);
  );
  return([kecalle,ecallelog]);
}
ecalleu(z) = {
  local(y);
  y=vector(2);
  y[1]=subst(kecalle,x,z)+ecallelog*log(z);
  y[2]=subst(deriv(kecalle),x,z)+ecallelog/z;
  return(y);
}
ecallel(z) = {
  local(y);
  y=vector(2);
  y[1]=subst(kecalle,x,z)+ecallelog*log(-z);
  y[2]=subst(deriv(kecalle),x,z)+ecallelog/z;
  return(y);
}

invcheta(z) = {
  local(n,y);
  if (ecallelog==0, initeta());
  z = z/exp(1)-1;
  n=0;
  while (abs(z)>0.05,z=log(z+1);n++);
  y = ecalleu(z)[1]+n;
}
slogeta(z) = {
  local(n,y);
  if (ecallelog==0, initeta());
  z = z/exp(1)-1;
  n=0;
  while (abs(z)>0.05,z=exp(z)-1;n++);
  y = ecallel(z)[1]-n;
  y = y- etarnorm;
  return(y);
}
cheta(z) = {
  local (y,dy,n,lastyz,curyz,est,pgoal,k1x,m,i);
  if (ecallelog==0, initeta());
  i=floor(real(z)+42);
  if (i<0,i=0);
  z=z-i;
  k1x=polcoeff(kecalle,-1);
  est=1/(z/k1x + (ecallelog/k1x)*log(z/k1x));
  lastyz=100;
  y=ecalleu(est);
  curyz=abs(y[1]-z);
  pgoal=10^(-precis/1.3);
  n=1;
  while ((curyz>pgoal) && ((curyz<lastyz) || (n<3)),
    est=precision(est,precis);
    y=precision(y,precis);
    est=est+(z-y[1])/y[2];
    lastyz=curyz;
    y=ecalleu(est);
    curyz=abs(y[1]-z);
    n++;
  );
  for (m=1,i,est=exp(est)-1);
  est=(est+1)*exp(1);
  return(est);
}
sexpeta(z) = {
  local (y,dy,n,lastyz,curyz,est,pgoal,k1x,m,i);
  if (ecallelog==0, initeta());
  z=z+etarnorm;
  i=floor(-real(z)+40);
  if (i<0,i=0);
  z=z+i;
  k1x=polcoeff(kecalle,-1);
  est=1/(z/k1x + (ecallelog/k1x)*log(-z/k1x));
  lastyz=100;
  y=ecallel(est);
  curyz=abs(y[1]-z);
  pgoal=10^(-precis/1.1);
  n=1;
  while ((curyz>pgoal) && ((curyz<lastyz) || (n<3)),
    est=precision(est,precis);
    y=precision(y,precis);
    est=est+(z-y[1])/y[2];
    lastyz=curyz;
    y=ecallel(est);
    curyz=abs(y[1]-z);
    n++;
  );
  for (m=1,i,est=log(est+1));
  est=(est+1)*exp(1);
  return(est);
}
initeta(i) = {
  local(n,expxm,z);
  etaB=exp(1/exp(1));
  /* ecalle(exp(x)-1,18); */ /* this is the best we can do with \ps 21 default */
  if (i==0,i=40);  /* good enough for 50 digits of precision, */ /* 1/x term + constant term + 40 terms + klog term, */
  z=1.0; precis=precision(z);
  expxm=0; for (n=1,i+3,expxm=expxm+1.0*x^n/n!); expxm=expxm+O(x^(i+4));
  ecalle(expxm,i);  /* i=40, 1/x term + constant term + 40 terms + klog term, 43 terms total */
  etarnorm=0;
  etarnorm=slogeta(1); /* renorm so slogeta(1)=0 */
}

help(w) = {
  print ("help(); help2(); andrewjay(); help_penthex() help_eta(); uniquechart(); ");
  print ("fatou.gp generates the Abel function for iterating z <= exp(z)-1+k; f(z)");
  print ("\\p 38                                 default \\p 38 ~=33 decimal digits;");
  print ("loop(k,nlim,nskip,looplim);           main routine to generation abel function");
  print ("quietmode=1;                          disables output during loop iterations");
  print ("limitp=16;                            limits precision to 16 digits to speedup");
  print ("                                      this is the same as using looplim=16;");
  print ("sexpinit(b,nlim,nkip,looplim);        calls loop(log(log(b))+1);");
  print ("sexpinit(exp(1));  loop(1);           two examples for tetration for base e");
  print ("prtpoly(wtaylor,t,name);              pretty prints a polynomial");
  print ("slog(z); sexp(z);                     tetration and its inverse");
  print ("abel(z); invabel(z,est);              abel functions are used internally");
  print ("sexptaylor(center,radius,samples);    taylor series for sexp");
  print ("slogtaylor(center,radius,samples);    taylor series for slog");
  print ("invabeltaylor(c,r,s); abeltaylor(.);  taylor series for abel/invabel");
  print ("x2mode=1;                             if x2mode=1; iterate z^2+z+k instead");
  print ("x2mode=0;                             default iterates exp(x)-1+k");
  print ("fmode                                 defines f(z) for MakeGraph ");
  print ("fmode=2:slog(z) fmode=3:sexp(z)       fmode=0:abel(z); fmode=1:invabel(z)");
  print ("MakeGraph(width,height,x0,y0,x1,y1,filename, n);  /* f(z) defined by fmode */");
}
help2(w) = {
  print ("fatou.gp help2() menu;  help() for main help() menu");
  print ("ir=7/10; ctr=8/10;                 typical sexpinit setting; matrix_ir override");
  print ("notheta=0;                         loop without generating theta mapping");
  print ("throwp=4.5; throws=9.5;            throwaway precision,  throws skip precision");
  print ("k; kfromp(p); bfromp(p);           generate k, b from period of p");
  print ("invabel_sexp(z); sexp_invabel(z);  generate from invabel from sexp and reverse");
  print ("finv(z); fs(z); testz(z);          fs=exp(z)-1+k; finv=log(z+1-k)");
  print ("abelest(z,ct);                     abel function estimate without theta mapping");
  print ("thfunc(z,n); thtaylor(n,samples);  tht; tht2; upper and lower theta series");
  print ("ct=staylor(circc,circr,samples);   used in loop; taylor series for sfunc(z)");
  print ("formalischroder(fx,n);             fisl; also works for any f(x)");
  print ("formalschroder(fx,n);              fsl;  also works for any f(x)");
  print ("halfsexp(z);                       f(f(z))=sexp(z);");
  print ("serreverse(strip0fx(Ser(z)));      formal series reverse");
  print ("superf(z);  isuperf(z);            using upper formal schroder/ischroder");
  print ("superf2(z); isuperf2(z);           using lower formal schroder/ischroder");
  print ("z0rfunc(z); invz0rfunc(i);         used to generate z0h; z0l; for theta mapping");
  print ("loop1(i,ctsamples,thsamples);      iterate to get \"exact\" sol'n for ctsamples");
  print ("thest(z);          thetamode indicates theta taylor series was generated");
  print ("initsch(k);                        initializes L; L2; superf; superf2, z0h/z0l");
  print ("matrix_ir(B,lctr,ltht,myctr,myir); matrix solution with theta mapping! ");
  print ("  matrixradius=1   causes matrix_ir to mimics slower sexpinit lctr, ctr setting");
  print ("matrix_r(B,n);     matrix solution for B, n samples without theta mapping");
}
help_penthex(w) = {
  print ("fatou.gp help menu for pentation routines; help_pentation; help() for main help");
  print ("*PENTATION ROUTINES*    pentaylor(w,r)  pent polynomial centered at w, radius r");
  print ("pentinit(B);            calls initsexp(B), and then generates pentation base B");
  print ("pent(z)                 evaluates pentation at z");
  print ("ipent(z)                evaluates the inverse of pentation at z");
  print ("tfixed;                 tetration fixed point");
  print ("lambdat; pperiod;       sexp multiplier and pent period at the fixed point");
  print ("sexpupfixed()           generate sexp base with parabolic upper fixed point");
  print ("pentr=0.05;             use 0.01 value of pentr for \\p 57, high precision");
  print ("pentfixr=0.7;           default pentation fixed point sample radius");
  print ("help()                  fatou help menu");
  print ("help_pentation()        pentatonp help menu");
  print ("                                                   ");
  print ("*HEXATION ROUTINES*     will add taylor series later");
  print ("hexinit(i,n);           init hexation for base(i); calls matrix_h(n)");
  print ("hex_abel(z,pt);         internal hex_abel function evaluation     ");
  print ("pentf;                  pentation complex fixed point used by hex ");
  print ("lambda_pf;              pentation multiplier at pentf used by hex ");
  print ("hex_period;             hexation pseudo period                    ");
/*print ("pt; ptd;                pentation taylor series, and derivative   "); */
/*print ("hex_bias;               bias required so hex(1)=0                 "); */
  print ("hex(z); ihex(z)         hexation approxmation for z, inverse hex  ");
  print ("ihex_deriv(z);          inverse of hexation and its derivative    ");
}
help_eta() = {
  print ("fatou.gp help menu for eta/cheta routines; help_eta; help() for main help");
  print ("etaB=exp(1/exp(1));         The parabolic base eta uses ecalle");
  print ("ecalle(fz,n); initeta(i));  generate kecalle series for f(z); initeta: exp(x)-1");
  print ("kecalle, ecallelog;         ecallelog is a constant;");
  print ("ecalleu(z); ecallel(z)      evaluate upper and lower ecalle abel series");
  print ("invcheta(z); slogeta(z);    via upper and lower ecalle functions for exp(x)-1");
  print ("cheta(z);    sexpeta(z)     via the inverse o fthe ecalle functions");
}
{
   help();
}


/***********************************************/
/*                                             */
/*  andrew's slog and jay's accelerated slog   */
/*                                             */
/***********************************************/
/*
slog(z) ....
slog(log(x+1))+1 = slog(x+1)
slog'(y-1) = slog'(log(y))/(y)

Jay's equation uses this times n!
slog(exp(x)) = slog(x) + 1
slog(exp(x))-slog(x) = 1


Jay's implementation of Andrew's slog
http://math.eretrandre.org/tetrationforum/showthread.php?tid=33
aQ = matrix(80, 80, r, c, (c^(r-1))/(c!)-if(c-r+1,0,1))
aZ = matrix(80, 80, r, c, (c^(r-1))-if(c-r+1,0,c!))
b = vector(80, n, if(n-1,0,1))
cQ = matsolve(aQ, b~)
##
cZ = matsolve(aZ, b~);
##

accelerated slog
http://math.eretrandre.org/tetrationforum/showthread.php?tid=63

Andrew's slog
http://iteror.org/big/Source/articles/TetrationSuperlog_Robbins.pdf
*/

andrewjay() = {
  print ("computeandrewslog(72,0); /* returns andrewtaylor; m=1 realmode; works faster */");
  print ("allocatemem(); /* computeandrewslog needs a lot of memory */");
  print ("\\p 77  /* computeandrewslog(140,1); real mode needs higher of precision */");
  print ("loop(1);  /* initialize internal slog for comparison */");
  print ("quickfourier(13); /* gt=gtaylor(w,r,samples); using gfunc(z) */ ");
  print ("andrewslog(z) = subst(andrewtaylor,x,z)-1;");
  print ("andrewsexp(z,est);");
  print ("slogtheta(z); /* andrewslog/jayslog approximation from slog(z)+theta(slog(z))*/");
  print ("slogthetadelta(z); /* vs theta mapping */");
  print ("andrewslogerr(z);");
  print ("\\p  67 /* required precision for computejayslog(100) */");
  print ("computejayslog(100); /* returns jaytaylor; */ quickfourier(10); ");
  print ("jayslog(z);");
  print ("jayslogerr(z);");
  print ("logt(z);  /* log_10(abs(z)) for err plots */");
}

logt(z) = {
  z=abs(z);
  if (z>1E-65,
    z=log(z)/log(10);
  ,
    z=-67;
  );
  return(z);
}
gtaylor( w,r,samples) = {
  local(rinv,s,t,x1,y,z,tot,t_est,tcrc,halfsamples,wtaylor,terms);
/* outputs gie taylor series, the complex taylor series for sexp, */
/* default halfsamples=240 */
  if (samples==0, samples=240);  /* no matter how many sample points, the default gie series size is 200 halfsamples */
  halfsamples=samples/2;
  terms = floor(samples*200/240);
  t_est    = vector (samples,i,0);
  tcrc     = vector (samples,i,0);
  if (r==0,r=1);
  rinv = 1/r;
  wtaylor=0;
  for(s=1, samples, x1=-1/(samples)+(s/halfsamples); tcrc[s]=exp(Pi*I*x1); );

  /* uses the gfunc to allow for centering anywhere in the complex plane */
  for (t=1,samples, t_est[t] = gfunc(w+r*tcrc[t]); );

  for (s=0,terms-1,
    tot=0;
    for (t=1,samples,
      tot=tot+t_est[t];
      t_est[t]=t_est[t]*conj(tcrc[t]);
    );
    tot=tot/samples;
    if (s>=1, tot=tot*(rinv)^s);
    wtaylor=wtaylor+tot*x^s;
  );
  wtaylor=precision(wtaylor,precis);
  return(wtaylor);
}
computeandrewslog(i,m) = {
  local(aZ,b,cZ);
  modegfunc=1;  /* used for quickfourier */
  /* aQ = matrix(i, i, r, c, (c^(r-1))/(c!)-if(c-r+1,0,1)) */
  if (m==0, /* integer mode */
    aZ = matrix(i, i, r, c, (c^(r-1))-if(c-r+1,0,c!));
  ,
    aZ = matrix(i, i, r, c, (1.0/(r-1)!)*(c^(r-1)-if(c-r+1,0,c!)));
  );
  b = vector(i, n, if(n-1,0,1));
  /* cQ = matsolve(aQ, b~)  */
  /*                        */
  cZ = matsolve(aZ, b~);
  /*                        */

  andrewtaylor = 0; for (n=1,i,andrewtaylor=andrewtaylor+x^n*cZ[n]);
  return("returns andrewtaylor;");
}
andrewslog(z) = subst(andrewtaylor,x,z)-1;
andrewsexp(z,est) = {
  local (y,s,slop,lastyz,curyz,lest,ly,pgoal,pplus);
/* inverse sexp(z), using "est" as initial estimate */
/* est is an optional parameter, important near the singularity */
/* the 0.01 radius used for the slope means won't converge within 0.01 of singularity */
  lastyz=100;
  pplus=0; if (real(z)>0,z=z-1; pplus=1);
  y=andrewslog(est);
  curyz=abs(y-z);
  lest=est+curyz*0.05;
  ly=andrewslog(lest);
  pgoal=10^(-precis/1.3);
  /* generate the fixed point for pentation by iteration slog */
  s=1;
  while ((curyz>pgoal) && ((curyz<lastyz) || (s<3)),
    est=precision(est,precis);
    y=precision(y,precis);
    slop=(y-ly)/(est-lest);
    lest=est;
    ly=y;
    est=est+(z-y)/slop;
    lastyz=curyz;
    y=andrewslog(est);
    curyz=abs(y-z);
    s++;
  );
  if (curyz>0.1, print (curyz " bad result, need better initial est, slog(z,est)"));
  if (pplus==1, est=exp(est));
  return(est);
}

/* could add other 1-cyclic real valued function to generate quickfourier series transform for */
gfunc(z) = {
  local(y);
  y=log(z)/(2*Pi*I);
  if (real(y)>0, y--);
  if (modegfunc==1,
    y=andrewslog(sexp(y))-y;
  , if (modegfunc==2,
    y=jayslog(sexp(y))-y;
  , if (modegfunc==3,
    y=slog(andrewsexp(y))-y;
  )));
  return(y);
}

/* quickie fourier transform! */
quickfourier(n) = {
  local(z,i);
  gt=gtaylor(0,1);
  z=gt;gt=0;
  for (i=0,n,gt=gt+polcoeff(z,i)*x^i);
  for (i=-n,-1,gt=gt+conj(polcoeff(z,-i))*x^(i));
  return(2*abs(polcoeff(gt,1)));
}

slogtheta(z) = {
  local(y);
  y=slog(z);
  y=y+subst(gt,x,exp(y*2*Pi*I));
  return(y);
}
slogthetadelta(z) = {
  local(y);
  if (modegfunc==1, y=andrewslog(z), y=jayslog(z));
  y=andrewslog(z)-slogtheta(z);
}
andrewslogerr(z) = {
  local(y);
  y=andrewslog(z)-slog(z);
  if (imag(z)==0, y=real(y));
  return(y);
}


/* z = rlnlm*log(z-L) + rlnlm2*log(z-L2) + subst(e,x,(z-circc));  */
/* slog(exp(x))-slog(x) = 1                                       */
/* accelerated equation:                                          */
/* slog(exp(x)) - slog(x) + (rlnlm*log(exp(x)-L)+conj(rlnlm)*log(exp(x)-conj(L))) - (rlnlm*log(x-L)+conj(rlnlm)*log(x-conj(L))) = 1 */
/* slog(exp(x)) - slog(x) = 1 - (rlnlm*log(exp(x)-L)+conj(rlnlm)*log(exp(x)-conj(L))) + (rlnlm*log(x-L)+conj(rlnlm)*log(x-conj(L))) */

computejayslog(i) = {
  local(aZ,b,cZ,zacc);
  default(seriesprecision,i);
  modegfunc=2;  /* used for quickfourier */
  L=fixedk(1);
  /* aQ = matrix(i, i, r, c, (c^(r-1))/(c!)-if(c-r+1,0,1)) */
  aZ = matrix(i, i, r, c, (1.0/(r-1)!)*(c^(r-1)-if(c-r+1,0,c!)));
  b = vector(i, n, if(n-1,0,1));
  zacc = 1 - (rlnlm*log(exp(x)-L)+conj(rlnlm)*log(exp(x)-conj(L))) + (rlnlm*log(x-L)+conj(rlnlm)*log(x-conj(L)));
  zacc = real(zacc);
  for (n=1,i, b[n]=polcoeff(zacc,n-1); );
  /* cQ = matsolve(aQ, b~)  */
  /*                        */
  cZ = matsolve(aZ, b~);
  /*                        */

  jaytaylor = 0; for (n=1,i,jaytaylor=jaytaylor+x^n*cZ[n]);
  jtk = 0;
  jtk = -jayslog(0)-1;
  default(seriesprecision,21);
  return("returns jaytaylor;");
}
jayslog(z) = { subst(jaytaylor,x,z) + (rlnlm*log(z-L)+conj(rlnlm)*log(z-conj(L))) + jtk; }
jayslogerr(z) = {
  local(y);
  y=jayslog(z)-slog(z);
  if (imag(z)==0, y=real(y));
  return(y);
}

/* functions from uniquechart.gp code */
bfromp(z)=exp(exp(kfromp(z)-1));
plotn(z,bz,plotsexp) = {
  local(v,i,n,zi);
  n = length(z);
  v=vector(2*n);
  for (i=1,n,
    zi=z[i];
    if ((plotsexp==1) && (zi<>0), zi=sexp_invabel(z[i]));
    v[i*2-1]=real(zi);
    v[i*2  ]=imag(zi);
    v[i*2-1]=Strprintf ("%18.15f", v[i*2-1]);
    v[i*2  ]=Strprintf ("%18.15f", v[i*2  ]);
    if (bz && (zi==0),
      v[i*2-1]="";
      v[i*2  ]="";
    );
  );
  return(v);
};
writecsv(filename,z) = {
/* this version includes 8x speedup for performance of writing csv's vectors without brackets */
  local(i,n);
  n=length(z);
  i=1;
  while (n>8, write1(filename,z[i+0],",",z[i+1],",",z[i+2],",",z[i+3],",",z[i+4],",",z[i+5],",",z[i+6],",",z[i+7],",");n=n-8;i=i+8);
  if   (n==1, write(filename,z[i]);
  , if (n==2, write(filename,z[i+0],",",z[i+1]);
  , if (n==3, write(filename,z[i+0],",",z[i+1],",",z[i+2]);
  , if (n==4, write(filename,z[i+0],",",z[i+1],",",z[i+2],",",z[i+3]);
  , if (n==5, write(filename,z[i+0],",",z[i+1],",",z[i+2],",",z[i+3],","z[i+4]);
  , if (n==6, write(filename,z[i+0],",",z[i+1],",",z[i+2],",",z[i+3],","z[i+4],",",z[i+5]);
  , if (n==7, write(filename,z[i+0],",",z[i+1],",",z[i+2],",",z[i+3],","z[i+4],",",z[i+5],",",z[i+6]);
  , if (n==8, write(filename,z[i+0],",",z[i+1],",",z[i+2],",",z[i+3],","z[i+4],",",z[i+5],",",z[i+6],",",z[i+7]);
  ))))))));
}
sickle(filename,fbias,plottet,imagtet) = {
  local(lave,zave,t,y,z,v,i,n,z1,endat);
  lave=(sexp_invabel(L)+sexp_invabel(L2))/2;
  zave=slog(lave);
  zr  =(sexp_invabel(L)-sexp_invabel(L2))/2;
  v = vector(10);
  z = vector(5);
  z[1] = slog(lave-zr*0.80)-0.5+fbias;
  z[2] = slog(lave-zr*0.40)-0.5+fbias;
  z[3] = slog(lave+zr*0.00)-0.5+fbias;
  z[4] = slog(lave+zr*0.40)-0.5+fbias;
  z[5] = slog(lave+zr*0.80)-0.5+fbias;
  forstep(t=0,1,1/49,
    v[1] =(z[1]+t);
    v[2] =(z[2]+t);
    v[3] =(z[3]+t);
    v[4] =(z[4]+t);
    v[5] =(z[5]+t);
    yt=t*1999/1000-999.5/1000;
    y=slog(yt*zr+lave)+fbias;
    v[6] =(y-0.50);
    v[7] =(y-0.25);
    v[8] =(y+0.00);
    v[9] =(y+0.25);
    v[10]=(y+0.50);
    writecsv(filename,plotn(v));
  );
  write(filename);
  forstep(t=0,1,1/49,
    v[1] =sexp(z[1]+t);
    v[2] =sexp(z[2]+t);
    v[3] =sexp(z[3]+t);
    v[4] =sexp(z[4]+t);
    v[5] =sexp(z[5]+t);
    yt=t*1.9996-0.9998;
    y=slog(yt*zr+lave)+fbias;
    v[6] =sexp(y-0.50);
    v[7] =sexp(y-0.25);
    v[8] =sexp(y+0.00);
    v[9] =sexp(y+0.25);
    v[10]=sexp(y+0.50);
    writecsv(filename,plotn(v));
  );
  write(filename);
  z=vector(4);
  forstep(t=0,1,1/49,
    z[1]=lave+zr*exp(t*Pi*I);
    z[2]=lave+zr*exp(-t*Pi*I);
    z[3]=exp(z[1]*lnb);
    z[4]=log(z[2])/lnb;
    writecsv(filename,plotn(z));
  );
  z=vector(2);
  write(filename);
  if (plottet<>0,
    endat=plottet;
    for(t=0,500,
      z1 = (t/500)*(endat+1.98)-1.98;
      z[1] = z1;
      z[2] = sexp(z1+I*imagtet);
      writecsv(filename,plotn(z,1));
    );
  );
}
circchart(filename,longtheta,lctr,ltht) = {
  local(t,ti,y,y1,y2,z,v,w,n,startat,stopat,endat,endpt);
  startat=0;
  stopat=0;
  v=vector(10);
  w=vector(8);
  if (lctr==0, lctr=length(ct)-1);
  if (ltht==0, ltht=length(tht));
  if (thetamode<>0, for (t=1,lctr*2,
    tz=exp((t-0.5)*2*Pi*I/lctr);
    if (t<=lctr+1,
      v[1] = circc+circr*tz;
      v[8] = circc+ircircr*tz;
    ,
      v[1] = 0;
      v[8] = 0;
    );
    z    = circc+circr*ctr*tz;
    v[2] = z;
    y1 = fs(z);
    y2 = finv(z);
    if (abs(y1-circc)>abs(y2-circc),
      if (abs(y2-circc)>ircircr, v[3]=0, v[3]=y2);
    ,
      if (abs(y1-circc)>ircircr, v[3]=0, v[3]=y1);
    );
    if ((v[3]==0) && (t>lctr), stopat=1);
    if ((v[3]==0) && (t<lctr) && (startat==0), startat=1);
    if ((v[3]<>0) && (t<lctr) && (startat==1), startat=2);
    if ((v[3]==0) && (t<lctr) && (startat==2), startat=3);
    if ((v[3]<>0) && (startat==3), startat=4);
    if ((v[3]==0) && (startat==1), v[9]=z,  v[9]=0);
    if ((v[3]==0) && (startat==3), v[10]=z, v[10]=0);
    if (t<=ltht,
      v[4]=superf(zth-0.5+(t-0.5)/ltht);
      if (complextaylor,
        v[5]=superf2(ztl-0.5+(t-0.5)/ltht);
      ,
        v[5]=conj(v[4]);
      );
    ,
      v[4]=0;
      v[5]=0;
    );
    v[6]=0;
    v[7]=0;
    if (t>lctr,
      v[4]=0;
      v[5]=0;
    );
    if ((v[3]<>0) && (startat>2),
      v[6]=v[2];
      v[7]=v[3];
      v[2]=0;
      v[3]=0;
    );
    if ((t<lctr) && (startat<>2),
      v[2]=0;
      v[3]=0;
    );
    if (v[3]==0, v[2]=0);
    if (t<=lctr,writecsv(filename,plotn(v,1,1)));
    if (t>lctr,
      if (stopat==0, writecsv(filename,plotn(v,1,1)));
      if ((stopat<>0) & (t<=513), write(filename));
    );
  );
  );
  if (thetamode<>0, for (t=(2*lctr),513, write(filename)));
  if (thetamode==0, for (t=0,513, write(filename)));
  startat=0;
  endat=128;
  bias=0;
  if (longtheta==0,endpt=endat,endpt=longtheta);
  for(t=startat,endpt,
    if (t<=endat,
      tz=exp((t-bias)*Pi*I/endat+argc*I-Pi*I/2);
      z    = circc+circr*tz;
      w[1] = z;
      w[3] = finv(z);
      tz=exp((t-bias)*Pi*I/endat+argc*I+Pi*I/2);
      z    = circc+circr*tz;
      w[2] = z;
      w[4] = fs(z);
    );
    if (repelling,ti=t,ti=t+ltht+endpt+1);
    if (thetamode<>0,
      w[5] = superf(zth-0.5+(ti-endpt-0.5)/ltht);
/*    y=(ti-endpt-0.5)/ltht; print(y" "w[5]); */
      if (t>endat,
        w[1] = 0;
        w[2] = 0;
        w[3] = 0;
        w[4] = 0;
        w[6] = 0;
      ,
        if (complextaylor,
          w[6] = superf2(ztl-0.5+(t-endat-0.5)/ltht);
        ,
          w[6] = conj(w[5]);
        );
      );
    );
    writecsv(filename,plotn(w,1,1))
  );
}
uniquechart(w) = {
print ("requires fatou.gp; uniquechart() for this menu");
print ("/* more stuff in uniquechart.gp; ask Sheldon */");
print ("plotn(z);");
print ("bfromp(p);");
print ("writecsv(filename,z);");
print ("sickle(filename,fbias);");
print ("circchart(filename);");
print ("matrix_ir(32000+3000*I,,,8/9);");
print ("matrix_ir(0.2*I,400,90,14/15,45/46);");
print ("matrix_ir(0.15,400,90,14/15,45/46);");
print ("matrix_ir(0.1+I*1E-30,400,250,14/15,45/46);");
print ("sexp(z) invabel doesn't work too well for some of these bases");
print ("sickle finds bugs in invabel for base e!  invabel(-1.8206806279041703166094448554879 - 6.4801614304879953859598867303807*I);");
print ("");
}

