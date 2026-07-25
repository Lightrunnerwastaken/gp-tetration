/* Driver for the sfunc branch probe.
   Prepend a line defining brdps (working precision), brbasestr (base as a
   string, evaluated AFTER realprecision is raised) and brquiet. */
default(realprecision, brdps);
\r fatou_branchprobe.gp
default(realprecision, brdps);
quietmode = brquiet;
brbase = eval(brbasestr);
brreset();
gettime();
sexpinit(brbase, brdps, 4, 0);
brt = gettime();
print("BRPROBE dps=", brdps, " time_ms=", brt);
brreport();
print("BRPROBE anchor=", sexp(0.5));
quit();
