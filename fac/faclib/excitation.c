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

static CEPW_SCRATCH pw_scratch = {1, MAX_KL, 1E-1, 1E-1, 1E-3, 0, 0};

static MULTI *pk_array;
static MULTI *kappa0_array;
static MULTI *kappa1_array;

CEPW_SCRATCH *GetCEPWScratch() {
  return &pw_scratch;
}

int SetTEGridDetail(int n, double *x) {
  int i;
  
  n_tegrid = n;
  for (i = 0; i < n; i++) {
    tegrid[i] = x[i];
    log_te[i] = log(tegrid[i]);
  }
  return 0;
}

int SetTEGrid(int n, double emin, double emax) {
  int i;
  double del;

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


int SetCEPWOptions(int qr, int max, double eps_dipole, 
		   double eps_allowed, double eps_forbidden) {
  pw_scratch.qr = qr;
  if (max > MAX_KL) {
    printf("The maximum partial wave reached in Excitation: %d > %d\n", 
	   max, MAX_KL);
    abort();
  }
  pw_scratch.max_kl = max;
  pw_scratch.eps_dipole = eps_dipole;
  pw_scratch.eps_allowed = eps_allowed;
  pw_scratch.eps_forbidden = eps_forbidden;
  pw_scratch.nkl0 = 1;
  pw_scratch.kl[0] = 0;
  pw_scratch.log_kl[0] = -100.0;
  pw_scratch.nkl = 0;
  return 0;
}

int AddCEPW(int n, int step) {
  int i;
  for (i = pw_scratch.nkl0; i < n+pw_scratch.nkl0; i++) {
    if (i >= MAX_NKL) {
      printf("Maximum partial wave grid points reached in Excitation: "); 
      printf("%d > %d in constructing grid\n",  i, MAX_NKL);
      abort();
    }
    pw_scratch.kl[i] = pw_scratch.kl[i-1] + step;
    pw_scratch.log_kl[i] = log(pw_scratch.kl[i]);
    if ((int) (pw_scratch.kl[i]) > pw_scratch.max_kl) break;
  }
  pw_scratch.nkl0 = i;
  return 0;
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
  pw_scratch.kl[pw_scratch.nkl] = pw_scratch.max_kl+1;
  return 0;
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
  if (type < 0) {
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

int CERadialPk(int *nkappa, int *nkl, double **pk, 
	       short **kappa0, short **kappa1, int ie,
	       int k0, int k1, int k, int mode) {
  int type, ko2, i, j, m, t, q;
  int kf0, kf1, kpp0, kpp1, km0, km1;
  int kl0, kl1, kl0p, kl1p;
  int j0, j1, kl_max, j1min, j1max;
  ORBITAL *orb0, *orb1;
  int index[4];
  double te, e0, e1, z;
  double qkt, qkl, r, rp;  
  int js1, js3, js[4], ks[4];
  double sd, se;
  int last_kl0, second_last_kl0;
  double **p, s;
  short **kp0, **kp1;
  double eps_dipole, eps_allowed, eps_forbidden;

#ifdef PERFORM_STATISTICS
  clock_t start, stop;
#endif

#ifdef PERFORM_STATISTICS
  start = clock();
#endif

  ko2 = k/2;  
  index[0] = ie;
  index[1] = k0;
  index[2] = k1;
  index[3] = ko2;

  kp0 = (short **) MultiSet(kappa0_array, index, NULL);
  kp1 = (short **) MultiSet(kappa1_array, index, NULL);
  p = (double **) MultiSet(pk_array, index, NULL);
  if (*kp0) {
    type = (*p)[0];
    *nkappa = (*kp0)[0];
    *nkl = (*kp1)[0];
    if (mode == 0) {
      *kappa0 = *kp0 + 1;
      *kappa1 = *kp1 + 1;
    } else {
      *kappa0 = *kp1 + 1;
      *kappa1 = *kp0 + 1;
    }
    *pk = *p + 1;
    return type;
  } 
  *nkappa = (MAX_NKL+1)*(GetMaxRank()+1)*4;
  *kp0 = (short *) malloc(sizeof(short)*(*nkappa));
  *kp1 = (short *) malloc(sizeof(short)*(*nkappa));
  *p = (double *) malloc(sizeof(double)*((*nkappa)*n_tegrid));
  
  type = -1;
  if (k > 0) {
    orb0 = GetOrbital(k0);
    orb1 = GetOrbital(k1);
    kl0 = GetLFromKappa(orb0->kappa);
    kl1 = GetLFromKappa(orb1->kappa);
    kl0 = kl0/2;
    kl1 = kl1/2;
    if (IsEven(kl0 + kl1 + ko2)) {
      if (k == 2) type = 1;
      else type = 2;
    }
  } else {
    if (k0 == k1) type = 0;
  }
  e1 = egrid[ie];
  r = GetRMax();
  z = GetResidualZ();
  kl_max = r*sqrt(e1+2.0*z/r);
  kl_max = Min(kl_max, pw_scratch.max_kl);

  js[0] = 0;
  ks[0] = k0;
  js[2] = 0;
  ks[2] = k1;		

  last_kl0 = 0;
  second_last_kl0 = 0;
  qkt = 0.0;
  q = 1;
  m = 1;
  z = (tegrid[0] / e1);
  if (z > 1.0) z = 1.0;
  eps_dipole = pw_scratch.eps_dipole * z;
  eps_allowed = pw_scratch.eps_allowed * z;
  eps_forbidden = pw_scratch.eps_forbidden;
  if (mode == 0) {
    js1 = 1;
    js3 = 3;
  } else {
    js1 = 3;
    js3 = 1;
  }
  for (t = 0; !last_kl0; t++) {
    if (second_last_kl0) {
      last_kl0 = 1;
      kl0 += 1;
    } else {
      kl0 = pw_scratch.kl[t];
      if (kl0 > 10) { 
	rp = fabs(qkl/qkt);
	if (type == 1) {
	  if (rp < eps_dipole) kl_max = kl0;
	} else if (type > 0) {
	  if (rp < eps_allowed) kl_max = kl0;
	} else {
	  if (rp < eps_forbidden) kl_max = kl0;
	}
      }
      if (pw_scratch.kl[t+1] > kl_max) {      
	if (type == 1 && mode == 0) second_last_kl0 = 1;
	else last_kl0 = 1;
      }
    }
    qkl = 0.0;
    kl0p = 2*kl0;
    for (j0 = abs(kl0p-1); j0 <= kl0p+1; j0 += 2) {
      kpp0 = GetKappaFromJL(j0, kl0p); 
      km0 = kpp0;
      if (kl0 < pw_scratch.qr) {
	js[js1] = 0;
      } else {
	js[js1] = j0;
	if (IsOdd(kl0)) {
	  if (kpp0 < 0) km0 = -kpp0 - 1;
	} else {
	  if (kpp0 > 0) km0 = -kpp0 - 1;
	}
      }
      j1min = abs(j0 - k);
      j1max = j0 + k;
      if (mode == 1) {
	kf1 = OrbitalIndex(0, km0, e1);
	ks[3] = kf1;
      }
      for (j1 = j1min; j1 <= j1max; j1 += 2) {
	for (kl1p = j1 - 1; kl1p <= j1 + 1; kl1p += 2) {	
	  kl1 = kl1p/2;
	  if (last_kl0 && type == 1 && kl1 >= kl0) {
	    continue; 
	  }
	  kpp1 = GetKappaFromJL(j1, kl1p);
	  km1 = kpp1;
	  if (kl1 < pw_scratch.qr) {
	    js[js3] = 0;
	  } else {
	    js[js3] = j1;
	    if (IsOdd(kl1)) {
	      if (kpp1 < 0) km1 = -kpp1 - 1;
	    } else {
	      if (kpp1 > 0) km1 = -kpp1 - 1;
	    }
	  }
	  if (mode == 0) {
	    kf1 = OrbitalIndex(0, km1, e1);
	    ks[3] = kf1;
	  }
	  for (i = 0; i < n_tegrid; i++) {
	    te = tegrid[i];
	    e0 = e1 + te;
	    if (mode == 0) {
	      kf0 = OrbitalIndex(0, km0, e0);
	      ks[1] = kf0;
	    } else {
	      kf0 = OrbitalIndex(0, km1, e0);
	      ks[1] = kf0;	      
	    }
	    
	    if (kl1 >= pw_scratch.qr &&
		kl0 >= pw_scratch.qr) {
	      SlaterTotal(&sd, &se, js, ks, k, -1);
	    } else {
	      SlaterTotal(&sd, &se, js, ks, k, 1);
	    } 
	    r = sd + se;
	    if (i == 0 && (!sd && !se)) break;
	    (*p)[q++] = r;
	    r = r*r;
	    if (i == 0) {
	      qkl += r;
	    }
	  }
	  if (i > 0) {
	    (*kp0)[m] = kpp0;
	    (*kp1)[m] = kpp1;
	    *nkl = t;
	    m++;
	  }
	}
      }
    }    
    if (t == 0) qkt += qkl;
    else qkt += qkl*(kl0-pw_scratch.kl[t-1]);
  }

  if (type != 1) *nkl += 1;
  *nkappa = m-1;
  *p = realloc(*p, sizeof(double)*q);
  (*p)[0] = type;
  *pk = *p + 1;
  *kp0 = realloc(*kp0, sizeof(short)*m);
  *kp1 = realloc(*kp1, sizeof(short)*m);
  (*kp0)[0] = *nkappa;
  (*kp1)[0] = *nkl;
  if (mode == 0) {
    *kappa0 = *kp0 + 1;
    *kappa1 = *kp1 + 1;
  } else {
    *kappa0 = *kp1 + 1;
    *kappa1 = *kp0 + 1;
  }
    
#ifdef PERFORM_STATISTICS
  stop = clock();
  timing.rad_qk += stop-start;
#endif
  return type;
}

double CERadialQk(int ie, double te, int k0, 
		  int k1, int k2, int k3, int k) {
  int type, t, swap;
  int i, j, kl0, kl1, kl, nkappa, nkl, nkappap, nklp;
  short *kappa0, *kappa1, *kappa0p, *kappa1p, *tmp;
  double *pk, *pkp, y2[MAX_NKL], r, s, a, b, c, z;
  double *pk1, *pk2;  

#ifdef PERFORM_STATISTICS
  clock_t start, stop;
#endif

#ifdef PERFORM_STATISTICS
  start = clock();
#endif

  pk2 = NULL;
  type = CERadialPk(&nkappa, &nkl, &pk, &kappa0, &kappa1, ie,
		    k0, k1, k, 0);
  pk1 = malloc(sizeof(double)*nkappa);
  for (i = 0; i < nkappa; i++) {
    pk1[i] = InterpolatePk(te, type, pk);
    pk += n_tegrid;
  }
  nklp = nkl;
  if (k2 != k0 || k3 != k1) {
    type = CERadialPk(&nkappap, &nklp, &pkp, &kappa0p, &kappa1p, ie,
		      k2, k3, k, 0);
    pk2 = malloc(sizeof(double)*nkappap);
    for (i = 0; i < nkappap; i++) {
      pk2[i] = InterpolatePk(te, type, pkp);
      pkp += n_tegrid;
    }
  }
  swap = 0;
  if (nkl > nklp) {
    pk = pk1;
    pk1 = pk2;
    pk2 = pk;
    i = nkappa;
    nkappa = nkappap;
    nkappap = i;
    tmp = kappa0;
    kappa0 = kappa0p;
    kappa0p = tmp;
    tmp = kappa1;
    kappa1 = kappa1p;
    kappa1p = tmp;
    nkl = nklp;
    swap = 1;
  } else if (nkl < nklp) {
    swap = 1;
  }
  nklp = nkl-1;
  for (i = 0; i < nkl; i++) {
    pw_scratch.qk[i] = 0.0;
  }
  t = -1;
  kl0 = -1;
  a = 0.0;
  b = 0.0;
  c = 0.0;

  for (i = 0; i < nkappa; i++) {
    kl = GetLFromKappa(kappa0[i]);
    if (kl != kl0) {
      t++;
      kl0 = kl;
    }
    if (k2 == k0 && k3 == k1) {
      s = pk1[i]*pk1[i];
      if (t < nkl) {
	pw_scratch.qk[t] += s;
	if (type == 1 && t == nklp) {
	  kl1 = GetLFromKappa(kappa1[i]);
	  if (kl1 > kl0) a += s;
	}
      } else if (type == 1) {
	b += s;
      }
    } else {
      s = 0.0;
      for (j = 0; j < nkappap; j++) {
	if (kappa0[i] == kappa0p[j] && kappa1[i] == kappa1p[j]) {
	  s = pk1[i]*pk2[j];
	  break;
	}
      }
      if (t < nkl) {
	pw_scratch.qk[t] += s;
	if (type == 1 && t == nklp) {
	  if (j < nkappap) {
	    kl1 = GetLFromKappa(kappa1[i]);
	    if (kl1 > kl0) a += s;
	  }
	  if (swap) {
	    c += pk1[i]*pk1[i];
	  }
	}
      } else if (type == 1) {
	if (j < nkappap) b += s;
	else b += pk1[i]*pk1[i];
      } 
    }
  }

  spline(pw_scratch.log_kl, pw_scratch.qk, nkl, 1E30, 1E30, y2); 
  r = pw_scratch.qk[0];
  for (i = 1; i < nkl; i++) {
    r += pw_scratch.qk[i];
    kl0 = pw_scratch.kl[i-1];
    kl1 = pw_scratch.kl[i];
    for (j = kl0+1; j < kl1; j++) {
      splint(pw_scratch.log_kl, pw_scratch.qk, y2, nkl, 
	     LnInteger(j), &s);
      r += s;
    }
  }

  if (type == 1) {
    if (fabs(b) > fabs(a)) {
      if (swap && c != 0.0) b *= pw_scratch.qk[nklp]/c;
      z = GetResidualZ();
      z *= z*0.5;
      kl1 += 1;
      z /= kl1*kl1;
      s = egrid[ie] + te + z;
      s /= te;
      s *= b - a;
      r += s;
    }
  } else if (type > 0) { 
    t = nklp-1;
    a = pw_scratch.qk[nklp]/pw_scratch.qk[t];
    if (a < 1.0 && a > 0.0) {
      b = pw_scratch.kl[nklp] - pw_scratch.kl[t];
      b = pow(a, 1.0/b);
      s = pw_scratch.qk[t] * b/(1.0-b);
      r += s;
    }
  }
 
  free(pk1);
  if (pk2) {
    free(pk2);
  }

#ifdef PERFORM_STATISTICS
  stop = clock();
  timing.rad_qk += stop-start;
#endif
  return r;
}

  
int CERadialQkMSub(double *rq, int ie, double te, int k0, int k1,
		   int k2, int k3, int k, int kp, 
		   int nq, int *q) {
  int type1, type2, kl;
  int i, j, kl0, klp0, kl0_2, klp0_2, kl1;
  int nkappa, nkappap, nkl, nklp;
  short  *kappa0, *kappa1, *kappa0p, *kappa1p;
  double *pk1, *pk2;
  double *pk, *pkp, qy2[MAX_NKL];
  int km0, km1, j0, jp0, j1, nk4, kmp0, kmp1, km0_m, kmp0_m;
  int mi, mf, t, c0, cp0;
  double r, e0, e1, s;
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

  pk2 = NULL;
  type2 = 0;
  type1 = CERadialPk(&nkappa, &nkl, &pk, &kappa0, &kappa1, ie, 
		     k0, k1, k, 1);
  pk1 = malloc(sizeof(double)*nkappa);
  for (i = 0; i < nkappa; i++) {
    pk1[i] = InterpolatePk(te, type1, pk);
    pk += n_tegrid;
  }
  if (kp == k && k2 == k0 && k3 == k1) {
    pk2 = pk1;
    nkappap = nkappa;
    nklp = nkl;
    kappa0p = kappa0;
    kappa1p = kappa1;
  } else {
    type2 = CERadialPk(&nkappap, &nklp, &pkp, &kappa0p, &kappa1p, ie,
		       k2, k3, kp, 1);
    pk2 = malloc(sizeof(double)*nkappap);
    for (i = 0; i < nkappap; i++) {
      pk2[i] = InterpolatePk(te, type2, pkp);
      pkp += n_tegrid;
    }
    type2 = 1;
    if (nklp < nkl) nkl = nklp;
  }
  qk = (double **) malloc(sizeof(double *)*nq); 
  for (i = 0; i < nq; i++) { 
    qk[i] = (double *) malloc(nkl*sizeof(double)); 
    for (j = 0; j < nkl; j++) { 
      qk[i][j] = 0.0; 
    }   
  } 
    
  kl = -1;
  i = -1;
  e1 = egrid[ie];
  for (j = 0; j < nkappa; j++) {
    km0 = kappa0[j];
    km1 = kappa1[j];
    GetJLFromKappa(km0, &j0, &kl0);
    GetJLFromKappa(km1, &j1, &kl1);
    if (kl1 != kl) {
      i++;
      kl = kl1;
      if (i >= nkl) break;
    }
    kl0_2 = kl0/2;
    for (t = 0; t < nkappap; t++) {
      kmp0 = kappa0p[t];
      kmp1 = kappa1p[t];
      if (kmp1 != km1) continue;
      
      GetJLFromKappa(kmp0, &jp0, &klp0);
      klp0_2 = klp0/2;
      
      s = pk1[j]*pk2[t];
      s *= sqrt((j0+1.0)*(jp0+1.0)*(kl0+1.0)*(klp0+1.0));
      if (km0 != kmp0) { 
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
	  pha0 = GetPhaseShift(c0); 
	  phap0 = GetPhaseShift(cp0); 
	  cos_p[ite] = pha0-phap0; 
	} 
	if (n_tegrid == 1) r = cos_p[0];  
	else { 
	  spline(tegrid, cos_p, n_tegrid, 1E30, 1E30, y2); 
	  splint(tegrid, cos_p, y2, n_tegrid, te, &r); 
	} 
	s *= cos(r);
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
    spline(pw_scratch.log_kl, qk[iq], nkl, 1.0E30, 1.0E30, qy2);     
    r = qk[iq][0];
    for (i = 1; i < nkl; i++) { 
      r += qk[iq][i]; 
      kl0 = pw_scratch.kl[i-1]; 
      kl1 = pw_scratch.kl[i];        
      for (j = kl0+1; j < kl1; j++) {       
	splint(pw_scratch.log_kl, qk[iq], qy2, 
	       nkl, LnInteger(j), &s); 
	r += s; 
      }      
    }    
    rq[iq] = r; 
  } 

  for (i = 0; i < nq; i++) { 
    free(qk[i]); 
  } 
  free(qk); 
  free(pk1);
  if (type2) {
    free(pk2);
  }


#ifdef PERFORM_STATISTICS
  stop = clock();
  timing.rad_qk += stop-start;
#endif
  
  return type1;
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
  int nz, j1, j2;
  int nq, q[MAX_MSUB];
  double rq[MAX_MSUB][MAX_EGRID];
  int ie, type;
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
    j1 = lev1->pj;
    j2 = lev2->pj;
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
      if (!msub) {
	if (ang[i].k != ang[j].k) continue;
	for (ie = 0; ie < n_egrid; ie++) {
	  rqk[ie] = CERadialQk(ie, te, ang[i].k0, ang[i].k1, 
			       ang[j].k0, ang[j].k1, ang[i].k); 
	}
	
	if (n_usr > n_egrid) {
	  spline(log_egrid, rqk, n_egrid, 1.0E30, 1.0E30, y2);
	}
	for (ie = 0; ie < n_usr; ie++) {
	  if (n_usr > n_egrid) {
	    splint(log_egrid, rqk, y2, n_egrid, log_usr_egrid[ie], &r);
	  } else {
	    r = rqk[ie];
	  }
	  if (i != j) r *= 2.0; 
	  r /= (ang[i].k + 1.0);
	  s[ie] += c*r;
	}
      } else {
	for (ie = 0; ie < n_egrid; ie++) {
	  type = CERadialQkMSub(rqk_tmp, ie, te, ang[i].k0, ang[i].k1,
				ang[j].k0, ang[j].k1,
				ang[i].k, ang[j].k, nq, q);
	  for (t = 0; t < nq; t++) {
	    rq[t][ie] = rqk_tmp[t];
	    if (i != j) rq[t][ie] *= 2.0;
	  }
	}
	if (n_usr > n_egrid) {
	  for (t = 0; t < nq; t++) {
	    spline(log_egrid, rq[t], n_egrid, 1.0E30, 1.0E30, qy2[t]);
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
		splint(log_egrid, rq[m], qy2[m], 
		       n_egrid, log_usr_egrid[ie], &r);
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

  /* there is a factor of 4 coming from normalization and the 2 
     from the formula */
  if (!msub) {
    for (ie = 0; ie < n_usr; ie++) {
      s[ie] *= 8.0;     
    }
    return 1;
  } else {
    p = 0; 
    for (t = -j1; t <= 0; t += 2) {
      for (h = -j2; h <= j2; h += 2) {
	for (ie = 0; ie < n_usr; ie++) {
	  s[p++] *= 8.0;
	}
      }
    }
    return (p/n_usr);
  }
}

int CEQkTable(char *fn, int k, double te) {
  FILE *f;
  double x;
  int i, j, p, m, t;
  ORBITAL *orb1, *orb2;
  int n1, n2, kl1, kl2, j1, j2, k1, k2;

  f = fopen(fn, "w");
  if (!f) return -1;

  if (pw_scratch.nkl0 == 0) {
    SetCEPWOptions(1, 40, 1E-1, 1E-2, 1E-3);
  }
  if (pw_scratch.nkl == 0) {
    SetCEPWGrid(0, NULL, NULL);
  }
  
  t = 0;

  for (i = 3; i < 4; i++) {
    for (j = 6; j < 7; j++) {
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
	for (m = 0; m <= n_tegrid; m++) {
	  fprintf(f, "%d %10.3E %10.3E %10.3E %10.3E\n",
		  t, egrid[p], tegrid[m], pw_scratch.qk[m], x);
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
  double s[MAX_MSUB*MAX_USR_EGRID];
  int *alev;
  LEVEL *lev1, *lev2;
  double emin, emax, e, c, b;

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
    m = 0;
    for (i = 0; i < nlow; i++) {
      lev1 = GetLevel(low[i]);
      for (j = 0; j < nup; j++) {
	lev2 = GetLevel(up[j]);
	e = lev2->energy - lev1->energy;
	if (e > 0) m++;
	if (e < emin && e > 0) emin = e;
	if (e > emax) emax = e;
      }
    }
    e = 2.0*(emax-emin)/(emax+emin);
    if (m == 2) {
      SetTEGrid(2, emin, emax);
    } else if (e < EPS3) {
      SetTEGrid(1, 0.5*(emin+emax), emax);
    } else if (e < 0.2) {
      SetTEGrid(2, emin, emax);
    } else {
      SetTEGrid(n_tegrid, emin, emax);
    }
  }

  if (n_egrid == 0) {
    n_egrid = 6;
  }
  if (egrid[0] < 0.0) {
    if (n_usr <= 6) {
      SetEGridDetail(n_usr, usr_egrid);
    } else {
      SetEGrid(n_egrid, usr_egrid[0], usr_egrid[n_usr-1], 0);
    }
  }
  if (pw_scratch.nkl0 == 0) {
    if (msub) {
      SetCEPWOptions(0, 100, 1E-2, 1E-2, 1E-3);
    } else {
      SetCEPWOptions(0, 100, 1E-1, 1E-1, 1E-3);
    }
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
	b = 1.0+0.5*FINE_STRUCTURE_CONST2*(usr_egrid[ie]+e);
	for (m = 0; m < k; m++) {
	  c = s[m*n_usr+ie];
	  fprintf(f, "%-10.3E", c);
	  c *= PI * AREA_AU20/(2*(usr_egrid[ie]+e)*b*(j1+1));
	  fprintf(f, "%-10.3E  ", c);
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
  fprintf(f, "Dirac: %d, %6.1E, 1E: %6.1E, Slater: %6.1E, 2E: %6.1E\n", 
	  GetNumContinua(),
	  ((double)radt.dirac)/CLOCKS_PER_SEC, 
	  ((double)radt.radial_1e)/CLOCKS_PER_SEC,
	  ((double)radt.radial_slater)/CLOCKS_PER_SEC,
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

void _FreeExcitationKappa(void *p) {
  short *kp;
  kp = *((short **) p);
  free(kp);
}

int FreeExcitationPk(int ie) {
  ARRAY *b;

  b = pk_array->array;
  if (b == NULL) return 0;
  if (ie < 0) {
    MultiFreeData(b, pk_array->ndim, _FreeExcitationPk);
    b = kappa0_array->array;
    MultiFreeData(b, kappa0_array->ndim, _FreeExcitationKappa);
    b = kappa1_array->array;
    MultiFreeData(b, kappa1_array->ndim, _FreeExcitationKappa);
  } else {
    b = (ARRAY *) ArrayGet(b, ie);
    if (b) {
      MultiFreeData(b, pk_array->ndim - 1, _FreeExcitationPk);
    }
    b = kappa0_array->array;
    b = (ARRAY *) ArrayGet(b, ie);
    if (b) {
      MultiFreeData(b, kappa0_array->ndim - 1, _FreeExcitationKappa);
    }
    b = kappa1_array->array;
    b = (ARRAY *) ArrayGet(b, ie);
    if (b) {
      MultiFreeData(b, kappa1_array->ndim - 1, _FreeExcitationKappa);
    }
  }

  return 0;
}

int InitExcitation() {
  int blocks[] = {6, 10, 10, 6};
  int ndim = 4;

  pk_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(pk_array, sizeof(double *), ndim, blocks);
  kappa0_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(kappa0_array, sizeof(short *), ndim, blocks);
  kappa1_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(kappa1_array, sizeof(short *), ndim, blocks);
  
  n_egrid = 0;
  n_tegrid = 0;
  egrid[0] = -1.0;
  tegrid[0] = -1.0;
}

  
  




