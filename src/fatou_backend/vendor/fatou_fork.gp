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
/* thlogk= 1/10; thlogk isn't used anymore :( see c:\pari\fatou_backup....gp for more details :( */
superfr = 0.1;
isuperfr = 0.01;
fmode=0;
invabeli=1/4;
complextaylor=1;
superfk=8;
efam=0; /* exp-011b: set by loop() — base-e family gets the fft extraction */
/* exp-012: the fs/finv walk in sfunc depends only on the base map, NOT on
   the evolving ct/theta state — cache walk endpoints per sample index
   across loop iterations while the sampling grid is unchanged. */
swon=0; swidx=0; swkey=0; swz=0; swn=0; swvalid=0; swb=0;
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
thfunc(z,n) = {
  local(y);
  if (n<>2,
    y=log(z)/(2*Pi*I);
    y = abelest(superf(zth+y),ct)-y-zth;
    return(y);
  ,
    y=log(z)/(-2*Pi*I);
    y = abelest(superf2(ztl+y),ct)-y-ztl;
    return(y);
  );
}
/* taylor series function ussing thfunc which is defined by sfuncmode */
thtaylor(n,samples) = {
  local(s,t,x1,y,z,tot,t_est,tcrc,halfsamples,wtaylor,terms);
  if (samples==0, samples=120);  /* no matter how many sample points, the default gie series size is 200 halfsamples */
  samples = floor(samples);
  terms=samples-1;
  t_est    = vector (samples,i,0);
  tcrc     = vector (samples,i,0);
  wtaylor=0;
  for(s=1, samples,
    x1 = -1 + -1/(samples) + (2*s/samples); /* -Pi to Pi */
    tcrc[s] = exp(Pi*I*x1);
    t_est[s]= thfunc(tcrc[s],n);
  );

  for (s=0,terms,
    tot=0;
    for (t=1,samples,
      tot=tot+t_est[t];
      t_est[t]=t_est[t]*conj(tcrc[t]);
    );
    tot=tot/samples;
/*  if (s>=1, tot=tot*(rinv)^s );   */
    wtaylor=wtaylor+tot*x^s;
  );
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
  local(n,y,y1);
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
  y = log(y)*rlnlm + n;
  while (abs(imag(y))>abs(imag(y+Period)),y=y+Period);
  while (abs(imag(y))>abs(imag(y-Period)),y=y-Period);
  return(y);
}
isuperf2(z) = {
  local(n,y);
  y=z;
  n=0;
  while (abs(y-L2)>isuperfr2,
    y=finv(y,L2);
    n++;
  );
  y = subst(fsl2,x,y-L2);
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
  local(z1,z2,y1,y2,z,myps,arp);
  z = 1.0;
  precis=precision(z);
  myps = default(seriesprecision);
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
    y1 = abelest(z,ct) + n;
  ,
    z=zc;
    y2 = abelest(zc,ct); /* exp-015: lazy — only the theta path needs it */
    /* use theta mapping if abs(y-circc)>(ir*circr) */
    if (imag(y2)>0,
      y1 = isuperf(z)+polcoeff(tht,0);
      while (abs(y1-y2)>abs(y1-y2+Period), y1=y1+Period);
      while (abs(y1-y2)>abs(y1-y2-Period), y1=y1-Period);
      y1=y1-polcoeff(tht,0);
      z = y1;
      n = I;
      y1 = y1 + subst(tht,x,exp((y1-zth)*2*Pi*I));
    ,
      y1 = isuperf2(z)+polcoeff(tht2,0);
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

slog(z) = {
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
    pgoal=10^(-precis/1.2);
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

invabel(z,est) = {
  local (t,tc,y,curyz,nest,estp,estn);
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

      y = z - log(L-L2)*rlnlm2 - subst(ct,x,(L-circc)) - sfunczero;
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

      y = z - log(L2-L)*rlnlm - subst(ct,x,(L2-circc)) - sfunczero;
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
    for (n=1,t,est=safefs(est));
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
  /* center inflection point at the origin */
  y = invabel(z-rslog);
  y = (y-k+1)/lnb;
  if ((imag(z)==0) && (complextaylor==0), y=real(y));
  return(y);
}

sexptaylor(w, r, samples) = {
  local(y);
  y = abeltaylor(w-rslog,r,samples,1);
  y = (y-k+1)/lnb;
  return(y);
}

slogtaylor(w, r, samples) = {
  local(y,z0);
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

staylor( w,r,samples) = {
  local(rinv,s,t,x1,y,y0,y1,y2,st,z,tot,t_est,tcrc,halfsamples,wtaylor,terms,om,c0,mu,c1,G,coeffs);
  if (samples==0, samples=240);  /* no matter how many sample points, the default gie series size is 200 halfsamples */
  /* exp-011b: the fft extraction needs a power-of-2 grid; the grid change
     costs the factor-2 bases their delicate sample/terms co-evolution
     (gate: 30-43 digits), so it applies ONLY to the base-e family where
     the expensive long inits live. */
  if (efam, z=1; while (z<samples, z=z*2); samples=z);
  terms=samples;
  if (complextaylor==0, samples=samples/2);
  t_est    = vector (samples,i,0);
  tcrc     = vector (samples,i,0);
  if (r==0,r=1);
  rinv = 1/r;

  /* exp-012: enable the sfunc walk cache while sampling an unchanged grid
     (e-family only; grid identity = [samples, w, r]). */
  if (efam,
    if (swkey != [samples, w, r],
      swkey = [samples, w, r];
      swz = vector(samples);
      swn = vector(samples);
      swb = vector(samples);
      swvalid = vector(samples);
    );
    swon = 1;
  );
  if (complextaylor,
    for(s=1, samples,
      x1=-1+-1/(samples)+(2*s/samples);
      tcrc[s]=exp(Pi*I*x1); /* -Pi to Pi */
      swidx = s;
      t_est[s]=sfunc(w+r*tcrc[s]*exp(I*argc))[1];
    );
  ,

    for(s=1, samples,
      x1=-1/(2*samples)+(s/samples);
      tcrc[s]=exp(Pi*I*x1);
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
     (research/tools/fft_extraction_proto.gp). All other bases keep the
     original rotation loop (their accuracy depends on the unrounded grid). */
  if (efam,
    if (complextaylor,
      om = exp(2*Pi*I/samples);
      c0 = exp(-Pi*I*(1+1/samples));
      G = fft(powers(om^(-1), samples-1), t_est);
      coeffs = vector(terms, s, (rinv^s/samples) * conj(c0)^s * om^(-s) * G[(s%samples)+1]);
    ,
      mu = exp(Pi*I/samples);
      c1 = exp(-Pi*I/(2*samples));
      G = fft(powers(mu^(-1), 2*samples-1), concat(t_est, vector(samples, i, 0)));
      coeffs = vector(terms, s, (rinv^s/samples) * real(conj(c1)^s * mu^(-s) * G[(s%(2*samples))+1]));
    );
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
  initsch(kc);
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
    /* exp-007d: extend the iteration cap incrementally, but ONLY for the
       base-e family (kc=log(log(b))+1=1). Measured: for b=e true accuracy
       tracks the contour residual 1:1 and extension lifts it massively
       (e|200: 64->193 true digits). For other bases truth ~ 2x contour-re
       and iterating past the caller cap DESTROYS accuracy (2|80: 80->41.6
       true) — they keep the caller's behavior exactly. */
    if ((n>=nlim-1) && (limitp==0) && (re<looplim) && efam,
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
      thsamples=floor(re*tht_re_mult+6);
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
pentinit(B,skip) = {
  local(s,lx,tsexp,z1,z2,n,errp,errz,pgoal);
  if (skip==0, errp=sexpinit(B)[1]); /* initialize sexp using fatou.gp program */
  s=0;
  lx=100;
  tfixed=-1.97;
  pgoal=10^(4-precis);
  /* generate the fixed point for pentation by iteration slog */
  while (abs(lx-tfixed)>pgoal && (s<500) && (lx<>tfixed),
    lx=tfixed;
    tfixed=slog(tfixed);
    s++;
  );
  for (i=1,10,tfixed=slog(tfixed));
  /* generate taylor series with radius, 0.7* radius of convergence */
  /* this will give us an accurate slope and taylor series at tfixed */
  /* the slope and taylor series will be used for the pentation schroder function solution */

  tslog = slogtaylor(tfixed,1);  /* slog has larger radis of convergence at fixed point than sexp */
  if (imag(B)==0, tslog=real(tslog));
  lambdat = 1/polcoeff(tslog,1);  /* slope at fixed point of slog */
  pentz  = formalischroder(tslog,21)+tfixed;
  ipentz = formalschroder(tslog,21);

  pperiod=2*Pi*I/log(lambdat);
  pcent=0;              /* set to zero before we normalize so pent(0)=1 */
  pcent=ipent(1); /* normalize pentation so pent(0)=1, pcent=pent^{-1}(1) */
/*generate an error estimate for pentation ischroder series */
  z1 = subst(pentz,x,pentr);
  z2 = subst(pentz,x,pentr/lambdat);
  n=0;while (z1<1,z1=sexp(z1);n++);
  n++;while (n>0, z2=sexp(z2);n--;);
  errz = -log(abs(z1-z2))/log(10);
  errz = precision(errz,9);
  if (errp<errz, errz=errp); /* if pentation precision<tetration precision ... */
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
pent(z)={
  /* generate pentation approximation, pentinit first */
  /* generated from from lambdat, the slope at the fixed point, and the pentz ischroder function */
  local(n,i,y);
  y=lambdat^(z+pcent);
  n=0;
  while (abs(y)>pentr, y=y/lambdat;n++); /* iterate n times until y<pentr, pentr=0.05 */
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
  while (abs(y-tfixed)>pentr, y=slog(y);n++); /* iterate slog n times until y is within pentr of tfixed, pentr=0.05 */
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
hexinit(i,n)={
  if (n==0, n=256);
  pentinit(i);
  print();
  matrix_h(256);
}


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

