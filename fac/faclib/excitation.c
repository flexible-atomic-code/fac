#include "excitation.h"

static int n_usr = 0;
static double usr_egrid[MAX_USR_EGRID];
static double log_usr_egrid[MAX_USR_EGRID];
static int usr_egrid_type = 0;

static int n_egrid = 0;
static double egrid[MAX_EGRID];
static double log_egrid[MAX_EGRID];
static int n_tegrid = 0;
static double tegrid[MAX_TEGRID];
static double log_te[MAX_TEGRID];

static EXCIT_TIMING timing = {0, 0, 0};

static struct { 
  int qr;
  int max_kl;
  int max_kl_1;
  int max_kl_f;
  double tolerence;
  int nkl0;
  int nkl;
  int nkappa;
  int nklp;
  int nkappap;
  int kappa0[2*MAX_NKAPPA];
  int kappa1[2*MAX_NKAPPA];  
  double *pk;
  double pk1[MAX_NKAPPA];
  double pk2[MAX_NKAPPA];
  double kl[MAX_NKL];
  double qk[MAX_NKL];
  double qk_y2[MAX_NKL];
} pw_scratch = {1, MAX_KL, 15, 2, EPS3, 0, 0, 0};

static MULTI *pk_array;

int SetTEGrid(int n, double emin, double emax) {
  int i, ie;
  double del, log_del;

  if (n < 1) {
    n_tegrid = 0;
    tegrid[0] = -1.0;
    return 0;
  }

  if (emin < 0.0) {
    tegrid[0] = emin;
    return 0;
  }

  if (n > MAX_TEGRID) {
    printf("Max # of grid points reached \n");
    return -1;
  }

  if (n == 1) {
    n_tegrid = 1;
    tegrid[0] = emin;
    log_te[0] = log(emin);
    return 0;
  }

  if (n == 2) {
    n_tegrid = 2;
    tegrid[0] = emin;
    tegrid[1] = emax;
    log_te[0] = log(emin);
    log_te[1] = log(emax);
    return 0;
  }

  if (emax < emin) {
    printf("emin must > 0 and emax < emin\n");
    return -1;
  }
  
  n_tegrid = n;
  
  del = emax - emin;
  del /= n-1.0;
  tegrid[0] = emin;
  log_te[0] = log(emin);
  for (i = 1; i < n; i++) {
    tegrid[i] = tegrid[i-1] + del;
    log_te[i] = log(tegrid[i]);
  }
  
  return 0;
}


int SetCEPWOptions(int qr, int max, int max_1, int max_f, double tolerence) {
  pw_scratch.qr = qr;
  if (max > MAX_KL) {
    printf("The maximum partial wave reached in Excitation: %d\n", MAX_KL);
    abort();
  }
  pw_scratch.max_kl = max;
  pw_scratch.max_kl_1 = max_1;
  pw_scratch.max_kl_f = max_f;
  if (tolerence >= EPS16) pw_scratch.tolerence = tolerence;
  else tolerence = EPS3;
  pw_scratch.nkl0 = 1;
  pw_scratch.kl[0] = 0;
  pw_scratch.nkl = 0;
  return 0;
}

int AddCEPW(int n, int step) {
  int i;
  for (i = pw_scratch.nkl0; i < n+pw_scratch.nkl0; i++) {
    if (i >= MAX_NKL) {
      printf("Maximum partial wave grid points reached in Excitation: %d\n", 
	     MAX_NKL);
      abort();
    }
    pw_scratch.kl[i] = pw_scratch.kl[i-1] + step;
    if ((int) (pw_scratch.kl[i]) > pw_scratch.max_kl) break;
  }
  pw_scratch.nkl0 = i;
}

int SetCEPWGrid(int ns, int *n, int *step) {
  int i, m, k, j;
  if (ns > 0) {
    for (i = 0; i < ns; i++) {
      AddCEPW(n[i], step[i]);
    }
    k = step[ns-1]*2;
    j = 2;
  } else {
    ns = -ns;
    if (ns == 0) ns = 5;
    AddCEPW(ns, 1);
    k = 2;
    j = 2;
  } 

  m = pw_scratch.kl[pw_scratch.nkl0-1];
  while (m+k <= pw_scratch.max_kl) {
    AddCEPW(j, k);
    m = pw_scratch.kl[pw_scratch.nkl0-1];
    k *= 2;
  }
  pw_scratch.nkl = pw_scratch.nkl0;
}
  

int SetEGridDetail(int n, double *xg) {
  int i;
 
  n_egrid = n;
  for (i = 0; i < n; i++) {
    egrid[i] = xg[i];
    log_egrid[i] = log(egrid[i]);
  }

  return 0;
}

int SetEGrid(int n, double emin, double emax, int type) {
  double del;
  int i;

  if (n < 1) {
    n_egrid = 0;
    return -1;
  }
  if (emin < 0.0) {
    egrid[0] = emin;
    return 0;
  }
  if (n > MAX_EGRID) {
    printf("Max # of grid points reached \n");
    return -1;
  }

  if (emax < emin) {
    printf("emin must > 0 and emax < emin\n");
    return -1;
  }

  egrid[0] = emin;
  log_egrid[0] = log(egrid[0]);
  n_egrid = n;
  if (type == 0) {
    del = emax - emin;
    del /= n-1.0;
    for (i = 1; i < n; i++) {
      egrid[i] = egrid[i-1] + del;
      log_egrid[i] = log(egrid[i]);
    }
  } else {
    del = log(emax) - log(emin);
    del = del/(n-1.0);
    del = exp(del);
    for (i = 1; i < n; i++) {
      egrid[i] = egrid[i-1]*del;
      log_egrid[i] = log(egrid[i]);
    }
  }

  return 0;
}

int SetUsrEGridDetail(int n, double *xg, int type) {
  int i; 

  if (n < 1) {
    printf("Grid points must be at least 1\n");
    return -1;
  }
  if (n > MAX_USR_EGRID) {
    printf("Max # of grid points reached \n");
    return -1;
  }

  n_usr = n;
  usr_egrid_type = type;
  for (i = 0; i < n; i++) {
    usr_egrid[i] = xg[i];
    log_usr_egrid[i] = log(usr_egrid[i]);
  }
  return 0;
}

int SetUsrEGrid(int n, double emin, double emax, int type) {
  double del;
  int i;

  if (n < 1) {
    printf("Grid points must be at least 1\n");
    return -1;
  }
  if (n > MAX_USR_EGRID) {
    printf("Max # of grid points reached \n");
    return -1;
  }

  if (emin < EPS10 || emax < emin) {
    printf("emin must > 0 and emax < emin\n");
    return -1;
  }

  usr_egrid[0] = emin;
  log_usr_egrid[0] = log(usr_egrid[0]);
  n_usr = n;
  usr_egrid_type = type;

  del = log(emax) - log(emin);
  del = del/(n-1.0);
  del = exp(del);
  for (i = 1; i < n; i++) {
    usr_egrid[i] = usr_egrid[i-1]*del;
    log_usr_egrid[i] = log(usr_egrid[i]);
  }
  return 0;
}

int CEContinuaKappas(int ie, int k, int *nkl, int *nkappa,
		     int *kappa0, int *kappa1) {
  double e1;
  int i, kl_max, kl, kl2, m, km0, km1, jmin, jmax, j0, j1;
  double z, rmax;
  int last_kl;
#ifdef PERFORM_STATISTICS
  clock_t start, stop;
#endif

#ifdef PERFORM_STATISTICS
  start = clock();
#endif
  
  e1 = egrid[ie];
  if (pw_scratch.max_kl_f > 0) {
    kl_max = pw_scratch.max_kl_f*(e1/tegrid[0]);
    kl_max = Max(kl_max, pw_scratch.max_kl_1);
    kl_max = Min(pw_scratch.max_kl, kl_max); 
  } else {
    kl_max = pw_scratch.max_kl;
  }
  if (k > 2) kl_max = Min(kl_max, pw_scratch.max_kl_1);

  m = 0;
  last_kl = 0;
  for (i = 0; i <= pw_scratch.nkl0; i++) {
    if (i == pw_scratch.nkl0) {
      kl = pw_scratch.kl[i-1]+1;
      last_kl = 1;
    } else {
      kl = pw_scratch.kl[i];
      if (kl > kl_max) {
	kl = pw_scratch.kl[i-1]+1;
	last_kl = 1;
      }
    }
    kl2 = 2*kl;
    j0 = kl2 - 1;
    km0 = GetKappaFromJL(j0, kl2);
    jmin = j0 - k;
    jmax = j0 + k;
    for (j1 = jmin; j1 <= jmax; j1 += 2) {
      km1 = GetKappaFromJL(j1, j1-1);
      kappa0[m] = km0;
      kappa1[m++] = km1;
      km1 = GetKappaFromJL(j1, j1+1);
      kappa0[m] = km0;
      kappa1[m++] = km1;
    }

    j0 = kl2 + 1;
    km0 = GetKappaFromJL(j0, kl2);
    jmin = j0 - k;
    jmax = j0 + k;
    for (j1 = jmin; j1 <= jmax; j1 += 2) {
      km1 = GetKappaFromJL(j1, j1-1);
      kappa0[m] = km0;
      kappa1[m++] = km1;
      km1 = GetKappaFromJL(j1, j1+1);
      kappa0[m] = km0;
      kappa1[m++] = km1;
    }
    if (last_kl) break;
  }

  *nkl = i;
  *nkappa = m;
#ifdef PERFORM_STATISTICS
  stop = clock();
  timing.set_kappa += stop-start;
#endif

}

    
int CERadialPk(int ie, int k0, int k1, int k, 
	       int nkappa, int *kappa0, int *kappa1) {
  ORBITAL *orb0, *orb1;
  int i, m, ite;
  int type, km0, km1, kl0, kl1;
  int js[4], ks[4];
  int j0, j1, c0, c1, mode[MAX_TEGRID];
  double e0, e1, sd, se, x;
  int index[5];
  double **p;
#ifdef PERFORM_STATISTICS
  clock_t start, stop;
#endif

#ifdef PERFORM_STATISTICS
  start = clock();
#endif
  type = 4;
  if (k == 2) {
    orb0 = GetOrbital(k0);
    orb1 = GetOrbital(k1);
    i = GetLFromKappa(orb0->kappa);
    m = GetLFromKappa(orb1->kappa);
    if (IsOdd((i+m)/2)) {
	type = 1;
    }
  } else if (k == 0) {
    if (k0 == k1) {
      type = 0;
    }
  }

  index[0] = ie;
  index[1] = k0;
  index[2] = k1;
  index[3] = k/2;

  p = (double **) MultiSet(pk_array, index, NULL);
  if (*p) {
    pw_scratch.pk = *p ;
#ifdef PERFORM_STATISTICS
  stop = clock();
  timing.rad_pk += stop-start;
#endif
    return type;
  } 

  *p = (double *) malloc(sizeof(double)*n_tegrid*nkappa);
  pw_scratch.pk = *p;

  e1 = egrid[ie];
  for (ite = 0; ite < n_tegrid; ite++) mode[ite] = 1;

  for (i = 0; i < nkappa; i++) {
    km0 = kappa0[i];
    km1 = kappa1[i];

    if (km0 == 0 || km1 == 0) {
      for (ite = 0; ite < n_tegrid; ite++) {
	pw_scratch.pk[i*n_tegrid + ite] = 0.0;
      }
      continue;
    }

    js[0] = 0;
    js[2] = 0;
    GetJLFromKappa(km0, &j0, &kl0);
    kl0 /= 2;
    if (kl0 < pw_scratch.qr) {
      js[1] = 0;
    } else {
      js[1] = j0;
      if (IsOdd(kl0)) {
	if (km0 < 0) km0 = -km0-1;
      } else {
	if (km0 > 0) km0 = -km0-1;
      }
    }
    GetJLFromKappa(km1, &j1, &kl1);
    kl1 /= 2;
    if (kl1 < pw_scratch.qr) {
      js[3] = 0;
    } else {
      js[3] = j1;
      if (IsOdd(kl1)) {
	if (km1 < 0) km1 = -km1-1;
      } else {
	if (km1 > 0) km1 = -km1-1;
      }
    }
    
    for (ite = 0; ite < n_tegrid; ite++) {
      e0 = e1 + tegrid[ite];
      c0 = OrbitalIndex(0, km0, e0);
      c1 = OrbitalIndex(0, km1, e1);
    
      ks[0] = k0;
      ks[1] = c0;
      ks[2] = k1;
      ks[3] = c1;
                
      if (kl1 >= pw_scratch.qr && 
	  kl0 >= pw_scratch.qr && mode[ite] == 1) {
	SlaterTotal(&sd, &se, js, ks, k, -mode[ite]); 
      } else {	
	SlaterTotal(&sd, &se, js, ks, k, mode[ite]);
      }
      if (mode[ite] == 1 && kl0 >= 5 && kl1 >= 5) { 
	x = (sd+1.0 != 1.0)? (fabs(se/sd)):1.0; 
	if (x+1.0 == 1.0) { 
	  if (fabs(sd) < EPS10) { 
	    mode[ite] = 2; 
	  } 
	} else { 
	  if (x < pw_scratch.tolerence) { 
	    mode[ite] = 2; 
	  } 
	} 
      } 
      m = i*n_tegrid + ite;
      pw_scratch.pk[m] = sd + se;
    }
  }

#ifdef PERFORM_STATISTICS
  stop = clock();
  timing.rad_pk += stop-start;
#endif
  return type;
}

double CERadialQk(int ie, double te, int k0, int k1, int k2, int k3, int k) {
  int type1, type2;
  int i, j, kl0, kl1, nk4, kl, km;
  double *pk, r, s, ratio, a, b, z;
#ifdef PERFORM_STATISTICS
  clock_t start, stop;
#endif

#ifdef PERFORM_STATISTICS
  start = clock();
#endif

  CEContinuaKappas(ie, k, &pw_scratch.nkl, &pw_scratch.nkappa, 
		   pw_scratch.kappa0, pw_scratch.kappa1);

  type1 = CERadialPk(ie, k0, k1, k, pw_scratch.nkappa,
		     pw_scratch.kappa0, pw_scratch.kappa1);

  pk = pw_scratch.pk;
  for (i = 0; i < pw_scratch.nkappa; i++) {
    pw_scratch.pk1[i] = InterpolatePk(te, type1, pk);
    pk += n_tegrid;
  }

  if (k2 == k0 && k3 == k1) {
    type2 = type1;
    for (i = 0; i < pw_scratch.nkappa; i++) {
      pw_scratch.pk2[i] = pw_scratch.pk1[i];
    } 
  } else {
    type2 = CERadialPk(ie, k2, k3, k, pw_scratch.nkappa,
		       pw_scratch.kappa0, pw_scratch.kappa1);
    pk = pw_scratch.pk;
    for (i = 0; i < pw_scratch.nkappa; i++) {
      pw_scratch.pk2[i] = InterpolatePk(te, type2, pk);
      pk += n_tegrid;
    }
  }

  for (i = 0; i < pw_scratch.nkl; i++) {
    pw_scratch.qk[i] = 0.0;
  }

  nk4 = 4*(k+1);
  a = 0.0;
  b = 0.0;
  for (j = 0; j < pw_scratch.nkappa; j++) {
    i = j / nk4;
    s =  pw_scratch.pk1[j]*pw_scratch.pk2[j];
    if (i < pw_scratch.nkl) pw_scratch.qk[i] += s;
    if (type1 == 1 && type2 == 1) {
      if (i == pw_scratch.nkl-1) {
	kl = pw_scratch.kl[i];
	km = pw_scratch.kappa1[j];
	km = GetLFromKappa(km)/2;
	if (km > kl) {
	  a += s;
	}
      } else if (i == pw_scratch.nkl) {
	km = pw_scratch.kappa1[j];
	km = GetLFromKappa(km)/2;
	if (km == kl) {
	  b += s;
	}
      }
    }
  }
 
  spline(pw_scratch.kl, pw_scratch.qk, pw_scratch.nkl, 
	 1.0E30, 1.0E30, pw_scratch.qk_y2);
  
  r = pw_scratch.qk[0];
  for (i = 1; i < pw_scratch.nkl; i++) {
    r += pw_scratch.qk[i];
    kl0 = pw_scratch.kl[i-1];
    kl1 = pw_scratch.kl[i];
         
    for (j = kl0+1; j < kl1; j++) {      
      splint(pw_scratch.kl, pw_scratch.qk, pw_scratch.qk_y2, 
	     pw_scratch.nkl, (double)j, &s);
      r += s;
    }     
  }   
  
  if (type1 == 1 && type2 == 1) {
    z = GetResidualZ(1);
    z *= z*0.5;
    z /= (1.0+kl1)*(1.0+kl1);
    s = egrid[ie] + te + z;    
    s /= te;
    s *= (b-a);
    r += s;
  }
  
  r *= 4.0 / (k+1.0);

#ifdef PERFORM_STATISTICS
  stop = clock();
  timing.rad_qk += stop-start;
#endif

  return r;
}

int CERadialQkMSub(double *rq, int ie, double te, int k0, int k1,
		   int k2, int k3, int k, int kp, 
		   int nq, int *q) {
  int type1, type2;
  int i, j, n, kl0, klp0, klp1, kl0_2, klp0_2, kl1;
  int km0, km1, j0, jp0, j1, nk4, kmp0, kmp1, km0_m, kmp0_m;
  int mi, mf, t, c0, cp0, dkl;
  double *pk, r, e0, e1, s;
  double pha0, phap0;
  double s3j1, s3j2, s3j3, s3j4;
  double cos_p[MAX_TEGRID], y2[MAX_TEGRID];
  int ite, iq;
  double **qk;

#ifdef PERFORM_STATISTICS
  clock_t start, stop;
#endif     

#ifdef PERFORM_STATISTICS
  start = clock();
#endif
  
  CEContinuaKappas(ie, k, &pw_scratch.nkl, &pw_scratch.nkappa,
		   pw_scratch.kappa1, pw_scratch.kappa0);
  type1 = CERadialPk(ie, k0, k1, k, pw_scratch.nkappa,
		     pw_scratch.kappa0, pw_scratch.kappa1);
  pk = pw_scratch.pk;
  for (i = 0; i < pw_scratch.nkappa; i++) {
    pw_scratch.pk1[i] = InterpolatePk(te, type1, pk);
    pk += n_tegrid;
  }
  
  CEContinuaKappas(ie, kp, &pw_scratch.nklp, &pw_scratch.nkappap,
		   pw_scratch.kappa1+MAX_NKAPPA, 
		   pw_scratch.kappa0+MAX_NKAPPA);
  if (kp == k && k2 == k0 && k3 == k1) {
    type2 = type1;
    for (i = 0; i < pw_scratch.nkappa; i++) {
      pw_scratch.pk2[i] = pw_scratch.pk1[i];
    }
  } else {
    type2 = CERadialPk(ie, k2, k3, kp, pw_scratch.nkappap,
		       pw_scratch.kappa0+MAX_NKAPPA, 
		       pw_scratch.kappa1+MAX_NKAPPA);
    pk = pw_scratch.pk;
    for (i = 0; i < pw_scratch.nkappap; i++) {
      pw_scratch.pk2[i] = InterpolatePk(te, type2, pk);
      pk += n_tegrid;
    }
  }

  qk = (double **) malloc(sizeof(double *)*nq);
  for (i = 0; i < nq; i++) {
    qk[i] = (double *) malloc(pw_scratch.nkl*sizeof(double));
    for (j = 0; j < pw_scratch.nkl; j++) {
      qk[i][j] = 0.0;
    }  
  }

  nk4 = 4*(k + 1);
  for (j = 0; j < pw_scratch.nkappa; j++) {
    i = j / nk4;
    if (i == pw_scratch.nkl) continue;
    km0 = pw_scratch.kappa0[j];
    km1 = pw_scratch.kappa1[j];
    if (km0 == 0 || km1 == 0) continue;
    GetJLFromKappa(km0, &j0, &kl0);
    kl0_2 = kl0/2;
    GetJLFromKappa(km1, &j1, &kl1);
    for (t = 0; t < pw_scratch.nkappap; t++) {
      kmp0 = pw_scratch.kappa0[t+MAX_NKAPPA];
      kmp1 = pw_scratch.kappa1[t+MAX_NKAPPA];
      if (kmp0 == 0 || kmp1 == 0) continue;
      if (kmp1 != km1) continue;

      GetJLFromKappa(kmp0, &jp0, &klp0);
      klp0_2 = klp0/2;
      
      s = pw_scratch.pk1[j] * pw_scratch.pk2[t];       
      s *= sqrt((j0+1.0)*(jp0+1.0)*(kl0+1.0)*(klp0+1.0));
      
      if (km0 != kmp0) {
	dkl = (kl0_2 - klp0_2);
	km0_m = km0;
	kmp0_m = kmp0;
	if (kl0_2 >= pw_scratch.qr) {
	  if (IsOdd(kl0_2)) {
	    if (km0 < 0) km0_m = -km0 - 1;
	  } else {
	    if (km0 > 0) km0_m = -km0 - 1;
	  }
	}

	if (klp0_2 >= pw_scratch.qr) {
	  if (IsOdd(klp0_2)) {
	    if (kmp0 < 0) kmp0_m = -kmp0 - 1;
	  } else {
	    if (kmp0 > 0) kmp0_m = -kmp0 - 1;
	  }
	}
	
	for (ite = 0; ite < n_tegrid; ite++) {
	  e0 = e1 + tegrid[ite];
	  c0 = OrbitalIndex(0, km0_m, e0);
	  cp0 = OrbitalIndex(0, kmp0_m, e0);
	  pha0 = GetPhaseShift(c0, 0);
	  phap0 = GetPhaseShift(cp0, 0);
	  cos_p[ite] = phap0-pha0;
	}
	if (n_tegrid == 1) r = cos_p[ite]; 
	else {
	  spline(tegrid, cos_p, n_tegrid, 1E30, 1E30, y2);
	  splint(tegrid, cos_p, y2, n_tegrid, te, &r);
	}
	s *= cos(r);
	if (IsOdd(dkl)) s = -s;        
      }
      
      for (iq = 0; iq < nq; iq++) {
	rq[iq] = 0.0;
      }
      for (mi = -1; mi <= 1; mi += 2) {
	s3j1 = W3j(j0, 1, kl0, -mi, mi, 0);
	s3j2 = W3j(jp0, 1, klp0, -mi, mi, 0);
	for (iq = 0; iq < nq; iq++) {
	  mf = mi + q[iq];
	  s3j3 = W3j(j0, k, j1, -mi, -q[iq], mf);
	  s3j4 = W3j(jp0, kp, j1, -mi, -q[iq], mf);
	  rq[iq] += s3j1*s3j2*s3j3*s3j4;
	}
      }

      for (iq = 0; iq < nq; iq++) {
	qk[iq][i] += s*rq[iq];
      }
    }
  }
   
  for (iq = 0; iq < nq; iq++) {
    spline(pw_scratch.kl, qk[iq], pw_scratch.nkl, 
	   1.0E30, 1.0E30, pw_scratch.qk_y2);
    
    r = qk[iq][0];
    for (i = 1; i < pw_scratch.nkl; i++) {
      r += qk[iq][i];
      kl0 = pw_scratch.kl[i-1];
      kl1 = pw_scratch.kl[i];       
      for (j = kl0+1; j < kl1; j++) {      
	splint(pw_scratch.kl, qk[iq], pw_scratch.qk_y2, 
	       pw_scratch.nkl, (double)j, &s);
	r += s;
      }     
    }   
    r *= 4.0;
    rq[iq] = r;
  }

  for (i = 0; i < nq; i++) {
    free(qk[i]);
  }
  free(qk);

#ifdef PERFORM_STATISTICS
  stop = clock();
  timing.rad_qk += stop-start;
#endif
  
  return 0;
}
 
     
double InterpolatePk(double te, int type, double *pk) {
  double *x, y2[MAX_TEGRID];
  double r;

  if (n_tegrid == 1) {
    r = pk[0];
    return r;
  }
  if (type == 0 || type == 1) {
    x = log_te;
    te = log(te);
  } else {
    x = tegrid;
  }

  spline(x, pk, n_tegrid, 1.0E30, 1.0E30, y2);
  splint(x, pk, y2, n_tegrid, te, &r);
  
  return r;
}

int CollisionStrength(double *s, double *e, int lower, int upper, int msub) {
  int i, j, t, h, p, m;  
  LEVEL *lev1, *lev2;
  double te, c, r, s3j;
  ANGULAR_ZMIX *ang;
  int nz, type, m1, m2, j1, j2;
  int nq, q[MAX_MSUB];
  double rq[MAX_MSUB][MAX_EGRID];
  int ie;
  double rqk_tmp[MAX_MSUB];
  double rqk[MAX_EGRID], y2[MAX_EGRID], qy2[MAX_MSUB][MAX_EGRID];

  lev1 = GetLevel(lower);
  if (lev1 == NULL) return -1;
  lev2 = GetLevel(upper);
  if (lev2 == NULL) return -1;
  te = lev2->energy - lev1->energy;
  if (te <= 0) return -1;
  *e = te;

  if (msub) {
    j1 = GetHamilton(lev1->ham_index)->pj;
    j2 = GetHamilton(lev2->ham_index)->pj;
    DecodePJ(j1, NULL, &j1);
    DecodePJ(j2, NULL, &j2);
    j = 0;
    for (i = -j1-j2; i <= 0; i += 2) {
      q[j++] = i;
    }
    nq = j;
    p = 0;
    for (t = -j1; t <= 0; t += 2) {
      for (h = -j2; h <= j2; h += 2) {
	for (ie = 0; ie < n_usr; ie++) {
	  s[p++] = 0.0;
	}
      }
    }
  } else {
    p = 0;
    for (ie = 0; ie < n_usr; ie++) {
      s[ie] = 0.0;
    }
  }


  nz = AngularZMix(&ang, lower, upper, -1, -1);
  for (i = 0; i < nz; i++) {
    for (j = i; j < nz; j++) {
      c = ang[i].coeff * ang[j].coeff;
      if (fabs(c) < EPS10) continue;
      if (!msub) {
	if (ang[i].k != ang[j].k) continue;
	if (ang[i].k > MAX_K) continue;
	for (ie = 0; ie < n_egrid; ie++) {
	  rqk[ie] = CERadialQk(ie, te, ang[i].k0, ang[i].k1, 
			       ang[j].k0, ang[j].k1, ang[i].k); 
	}
	if (n_usr > n_egrid) {
	  spline(egrid, rqk, n_egrid, 1.0E30, 1.0E30, y2);
	}
	for (ie = 0; ie < n_usr; ie++) {
	  if (n_usr > n_egrid) {
	    splint(egrid, rqk, y2, n_egrid, usr_egrid[ie], &r);
	  } else {
	    r = rqk[ie];
	  }
	  if (i != j) r *= 2.0; 
	  s[ie] += c*r;
	}
      } else {
	if (ang[i].k > MAX_K || ang[j].k > MAX_K) continue;
	for (ie = 0; ie < n_egrid; ie++) {
	  CERadialQkMSub(rqk_tmp, ie, te, ang[i].k0, ang[i].k1,
			 ang[j].k0, ang[j].k1,
			 ang[i].k, ang[j].k, nq, q);
	  for (t = 0; t < nq; t++) {
	    rq[t][ie] = rqk_tmp[t];
	    if (i != j) rq[t][ie] *= 2.0;
	  }
	}
	for (t = 0; t < nq; t++) {
	  if (n_usr > n_egrid) {
	    spline(egrid, rq[t], n_egrid, 1.0E30, 1.0E30, qy2[t]);
	  }
	}
	p = 0;
	for (t = -j1; t <= 0; t += 2) {
	  for (h = -j2; h <= j2; h += 2) {
	    m = (-abs(t-h)-q[0])/2;
	    s3j = W3j(j1, ang[i].k, j2, -t, t-h, h);
	    if (ang[j].k != ang[i].k) {
	      s3j *= W3j(j1, ang[j].k, j2, -t, t-h, h);
	    } else {
	      s3j *= s3j;
	    }
	    for (ie = 0; ie < n_usr; ie++) {
	      if (n_usr > n_egrid) {
		splint(egrid, rq[m], qy2[m], 
		       n_egrid, usr_egrid[ie], &r);
	      } else {
		r = rq[m][ie];
	      }
	      s[p++] += c*r*s3j;
	    }
	  }
	}
      }
    }
  }

  if (nz > 0) {
    free(ang);
  }

  if (!msub) {
    for (ie = 0; ie < n_usr; ie++) {
      s[ie] *= 2.0;     
    }
    return 1;
  } else {
    p = 0; 
    for (t = -j1; t <= 0; t += 2) {
      for (h = -j2; h <= j2; h += 2) {
	for (ie = 0; ie < n_usr; ie++) {
	  s[p++] *= 2;
	}
      }
    }
    return (p/n_usr);
  }
}

int CEQkTable(char *fn, int k, double te) {
  FILE *f;
  double x;
  int i, j, p, m, t, n;
  ORBITAL *orb1, *orb2;
  int n1, n2, kl1, kl2, j1, j2, k1, k2;

  f = fopen(fn, "w");
  if (!f) return -1;

  if (pw_scratch.nkl0 == 0) {
    SetCEPWOptions(1, 50, 15, 2, 0.0);
  }
  if (pw_scratch.nkl == 0) {
    SetCEPWGrid(0, NULL, NULL);
  }
  
  n = GetNumBounds();
  t = 0;
  n = 4;
  for (i = 1; i < 2; i++) {
    for (j = 2; j < 3; j++) {
      orb1 = GetOrbital(i);
      n1 = orb1->n;
      k1 = orb1->kappa;
      kl1 = GetLFromKappa(k1);
      j1 = GetJFromKappa(k1);
      orb2 = GetOrbital(j);
      n2 = orb2->n;
      k2 = orb2->kappa;
      kl2 = GetLFromKappa(k2);
      j2 = GetJFromKappa(k2); 
      fprintf(f, "## %-3d ##Orbitals: (%d %d %d), (%d %d %d)\n", 
	      t++, n1, kl1, j1, n2, kl2, j2);
      for (p = 0; p < n_egrid; p++) {
	x = CERadialQk(p, te, i, j, i, j, k);
	for (m = 0; m < pw_scratch.nkl; m++) {
	  fprintf(f, " %10.3E %10.3E %3.0f %10.3E %10.3E %10.3E\n",
		  egrid[p], te,	pw_scratch.kl[m], pw_scratch.qk[m], 
		  pw_scratch.qk_y2[m], x);
	}
	fprintf(f, "\n\n");
	for (m = 0; m < pw_scratch.nkappa; m++) {
	  fprintf(f, "%-3d %-3d %10.3E %10.3E\n",
		  pw_scratch.kappa0[m], pw_scratch.kappa1[m], 
		  pw_scratch.pk1[m], pw_scratch.pk1[m]);
	}
	fprintf(f, "\n\n");
      }

      fprintf(f, "\n\n");
    }
  }
  
  fclose(f);
}

void spline(double *x, double *y, int n, 
	    double yp1, double ypn, double *y2) {
  int i, k;
  double p, qn, sig, un, *u;
  
  u = malloc(sizeof(double)*n);
  
  if (yp1 > 0.99E30) {
    y2[0] = u[0] = 0.0;
  } else {
    y2[0] = -0.5;
    u[0] = (3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
  }

  for (i = 1; i < n-1; i++) {
    sig = (x[i] - x[i-1])/(x[i+1]-x[i-1]);
    p = sig*y2[i-1]+2.0;
    y2[i] = (sig-1.0)/p;
    u[i] = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i] = (6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }
  if (ypn > 0.99E30) {
    qn = un = 0.0;
  } else {
    qn = 0.5;
    un = (3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
  }

  y2[n-1] = (un-qn*u[n-2])/(qn*y2[n-2]+1.0);
  for (k = n-2; k >= 0; k--) {
    y2[k] = y2[k]*y2[k+1]+u[k];
  }

  free(u);

}
  
int splint(double *xa, double *ya, double *y2a, 
	   int n, double x, double *y) {
  int klo, khi, k;
  double h, b, a;
  
  *y = 0.0;
  klo = 0;
  khi = n-1;
  while (khi-klo > 1) {
    k = (khi+klo) >> 1;
    if (xa[k] > x) khi = k;
    else klo = k;
  }
  h = xa[khi] - xa[klo];
  if (h == 0.0) return -1;
  a = (xa[khi] - x) / h;
  b = (x - xa[klo]) / h;
  *y = a*ya[klo] + b*ya[khi] + ((a*a*a-a)*y2a[klo] +
				(b*b*b-b)*y2a[khi]) * (h*h)/6.0;
  return 0;
}

void splie2(double *x1a, double *x2a, double **ya, 
	    int m, int n, double **y2a) {
  int j;

  for (j = 0; j < m; j++) {
    spline(x2a, ya[j], n, 1.0E30, 1.0E30, y2a[j]);
  }
}

void splin2(double *x1a, double *x2a, double **ya, double **y2a,
	    int m, int n, double x1, double x2, double *y) {
  int j;
  double *ytmp, *yytmp;

  ytmp = malloc(sizeof(double)*m);
  yytmp = malloc(sizeof(double)*m);
  for (j = 0; j < m; j++) {
    splint(x2a, ya[j], y2a[j], n, x2, &yytmp[j]);
  }
  spline(x1a, yytmp, m, 1.0E30, 1.0E30, ytmp);
  splint(x1a, yytmp, ytmp, m, x1, y);
  free(ytmp);
  free(yytmp);

}
 
int SaveExcitation(int nlow, int *low, int nup, int *up, int msub, char *fn) {
  int i, j, k, n, m, ie;
  int j1, j2;

#ifdef PERFORM_STATISTICS
  STRUCT_TIMING structt;
  ANGULAR_TIMING angt;
  RECOUPLE_TIMING recouplet;
  RAD_TIMING radt;
#endif

  FILE *f;
  double s[MAX_MSUB*MAX_USR_EGRID], energy;
  int *alev;
  LEVEL *lev1, *lev2;
  double emin, emax, e;

  if (n_usr == 0) {
    printf("No collisional energy specified\n");
    return -1;
  }
  f = fopen(fn, "w");
  if (!f) return -1;

  n = 0;
  if (nlow == 0 || nup == 0) {
    n = GetNumLevels();
    if (n <= 0) return -1;
    alev = malloc(sizeof(int)*n);
    if (!alev) return -1;
    
    for (i = 0; i < n; i++) alev[i] = i;

    if (nlow == 0) {
      nlow = n; 
      low = alev;
    }
    if (nup == 0) {
      nup = n;
      up = alev;
    }
  }
  
  if (n_tegrid == 0) {
    n_tegrid = 3;
  } 
  if (tegrid[0] < 0.0) {
    emin = 1E10;
    emax = 1E-10;
    for (i = 0; i < nlow; i++) {
      lev1 = GetLevel(low[i]);
      for (j = 0; j < nup; j++) {
	lev2 = GetLevel(up[j]);
	e = lev2->energy - lev1->energy;
	if (e < emin && e > 0) emin = e;
	if (e > emax) emax = e;
      }
    }
    if ((emax - emin) < EPS3) {
      SetTEGrid(1, 0.5*(emin+emax), emax);
    } else {
      SetTEGrid(n_tegrid, emin, emax);
    }
  }

  if (n_egrid == 0) {
    n_egrid = 5;
  }
  if (egrid[0] < 0.0) {
    if (n_usr <= 5) {
      SetEGridDetail(n_usr, usr_egrid);
    } else {
      SetEGrid(n_egrid, usr_egrid[0], usr_egrid[n_usr-1], 0);
    }
  }
  if (pw_scratch.nkl0 == 0) {
    SetCEPWOptions(1, 40, 15, 2, 1E-3);
  }
  if (pw_scratch.nkl == 0) {
    SetCEPWGrid(0, NULL, NULL);
  }

  fprintf(f, " TEGRID:   ");
  for (i = 0; i < n_tegrid; i++) {
    fprintf(f, "%10.4E ", tegrid[i]*HARTREE_EV);
  }
  fprintf(f, "\n");

  fprintf(f, " EGRID:    ");
  for (i = 0; i < n_egrid; i++) {
    fprintf(f, "%10.4E ", egrid[i]*HARTREE_EV);
  }
  fprintf(f, "\n\n");

  fprintf(f, "low  2J\tup   2J\tDelta_E\n");
  for (i = 0; i < nlow; i++) {
    j1 = LevelTotalJ(low[i]);
    for (j = 0; j < nup; j++) {
      j2 = LevelTotalJ(up[j]);
      k = CollisionStrength(s, &e, low[i], up[j], msub); 
      if (k < 0) continue;
      fprintf(f, "%-4d %-2d\t%-4d %-2d\t%10.4E\n", 
	      low[i], j1, up[j], j2, e*HARTREE_EV);
      for (ie = 0; ie < n_usr; ie++) {
	fprintf(f, "%-10.3E ", usr_egrid[ie]*HARTREE_EV);
	for (m = 0; m < k; m++) {
	  fprintf(f, "%-10.3E ", s[m*n_usr+ie]);
	}
	fprintf(f, "\n");
      }
      fprintf(f, "\n");
    }
  }
  fprintf(f, "\n");

#ifdef PERFORM_STATISTICS
  GetStructTiming(&structt);
  fprintf(f, "AngZMix: %6.1E, AngZFB: %6.1E, AngZxZFB: %6.1E, SetH: %6.1E DiagH: %6.1E\n",
	  ((double) (structt.angz_mix))/CLOCKS_PER_SEC,
	  ((double) (structt.angz_fb))/CLOCKS_PER_SEC,
	  ((double) (structt.angzxz_fb))/CLOCKS_PER_SEC,
	  ((double) (structt.set_ham))/CLOCKS_PER_SEC,
	  ((double) (structt.diag_ham))/CLOCKS_PER_SEC);
  fprintf(f, "AngZS: %6.1E, AngZFBS: %6.1E, AngZxZFBS: %6.1E, AddZ: %6.1E, AddZxZ: %6.1E\n",
	  ((double) (structt.angz_states))/CLOCKS_PER_SEC,
	  ((double) (structt.angzfb_states))/CLOCKS_PER_SEC,
	  ((double) (structt.angzxzfb_states))/CLOCKS_PER_SEC,
	  ((double) (structt.add_angz))/CLOCKS_PER_SEC,
	  ((double) (structt.add_angzxz))/CLOCKS_PER_SEC);

  GetAngularTiming(&angt);
  fprintf(f, "W3J: %6.1E, W6J: %6.1E, W9J: %6.1E\n", 
	  ((double)angt.w3j)/CLOCKS_PER_SEC, 
	  ((double)angt.w6j)/CLOCKS_PER_SEC, 
	  ((double)angt.w9j)/CLOCKS_PER_SEC);
  GetRecoupleTiming(&recouplet);
  fprintf(f, "AngZ: %6.1E, AngZxZ: %6.1E, Interact: %6.1E\n",
	  ((double)recouplet.angz)/CLOCKS_PER_SEC,
	  ((double)recouplet.angzxz)/CLOCKS_PER_SEC,
	  ((double)recouplet.interact)/CLOCKS_PER_SEC);
  GetRadTiming(&radt);
  fprintf(f, "Dirac: %d, %6.1E, 1E: %6.1E, 2E: %6.1E\n", 
	  GetNumContinua(),
	  ((double)radt.dirac)/CLOCKS_PER_SEC, 
	  ((double)radt.radial_1e)/CLOCKS_PER_SEC,
	  ((double)radt.radial_2e)/CLOCKS_PER_SEC);
  fprintf(f, "RadPk: %6.1E, SetKappa: %6.1E, RadQk: %6.1E\n",
	  ((double)timing.rad_pk)/CLOCKS_PER_SEC, 
	  ((double)timing.set_kappa)/CLOCKS_PER_SEC, 
	  ((double)timing.rad_qk)/CLOCKS_PER_SEC);

  fprintf(f, "\n");
#endif /* PERFORM_STATISTICS */

  fclose(f);
  if (n > 0) free(alev);
  return 0;
}

void _FreeExcitationPk(void *p) {
  double *dp;
  dp = *((double **) p);
  free(dp);
}

int FreeExcitationPk(int ie) {
  ARRAY *b;

  b = pk_array->array;
  if (b == NULL) return 0;
  if (ie < 0) {
    MultiFreeData(b, pk_array->ndim, _FreeExcitationPk);
  } else {
    b = (ARRAY *) ArrayGet(b, ie);
    if (b) {
      MultiFreeData(b, pk_array->ndim - 1, _FreeExcitationPk);
    }
  }

  return 0;
}

int InitExcitation() {
  int blocks[4] = {5, 10, 10, 5};
  int ndim = 4;

  pk_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(pk_array, sizeof(double *), ndim, blocks);
  n_egrid = 0;
  n_tegrid = 0;
  egrid[0] = -1.0;
  tegrid[0] = -1.0;
}

  
  




