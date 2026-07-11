/* micro-benchmarks of sfunc building blocks after a converged init.
   All times: gettime() ms for the stated repetition count. */
{
  zc0 = circc + circr*ctr*exp(Pi*I*0.37);
  print("MICRO terms_ct=", poldegree(ct)+1, " terms_tht=", poldegree(tht)+1);
  gettime();
  for(i=1,1000, y=fs(zc0));
  print("MICRO fs_x1000_ms=", gettime());
  for(i=1,1000, y=log(I*(zc0-L)));
  print("MICRO log_x1000_ms=", gettime());
  for(i=1,1000, y=exp(zc0));
  print("MICRO exp_x1000_ms=", gettime());
  w0 = zc0-circc;
  for(i=1,100, y=subst(ct,x,w0));
  print("MICRO ct_horner_x100_ms=", gettime());
  for(i=1,100, y=subst(tht,x,exp((zc0-zth)*2*Pi*I)));
  print("MICRO tht_subst_x100_ms=", gettime());
  for(i=1,100, y=abelest(zc0,ct));
  print("MICRO abelest_x100_ms=", gettime());
  for(i=1,100, y=sfunc(zc0));
  print("MICRO sfunc_coldwalk_x100_ms=", gettime());
  ctl = precision(ct, 60); wl = precision(w0, 60);
  print("MICRO coeff_prec_check=", precision(polcoeff(ctl, 10)));
  for(i=1,100, y=subst(ctl,x,wl));
  print("MICRO ct_horner_low_x100_ms=", gettime());
}
