#include "excitation.h"

static char *rcsid="$Id: excitation.c,v 1.38 2002/08/21 22:01:31 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

#define MAXMSUB  32
#define NPARAMS  4

static int qk_mode;
static double qk_fit_tolerance;

static int egrid_type = -1;
static int usr_egrid_type = -1;
static int pw_type = -1;

static int n_usr = 0;
static double usr_egrid[MAXNUSR];
static double log_usr[MAXNUSR];
static double xusr[MAXNUSR];
static double log_xusr[MAXNUSR];

static int n_egrid = 0;
static double egrid[MAXNE];
static double log_egrid[MAXNE];
static double egrid_min;
static double egrid_max;
static int egrid_limits_type = 0;

static int n_tegrid = 0;
static double tegrid[MAXNTE];
static double log_te[MAXNTE];

#define NGOSK 60
static double kgrid[NGOSK];
static double log_kgrid[NGOSK];
static double gos[NGOSK];
static double kint[NGOSK];
static double log_kint[NGOSK];
static double gosint[NGOSK];
static int n_born;
static double born_egrid[MAXNE];

#ifdef PERFORM_STATISTICS
static EXCIT_TIMING timing = {0, 0, 0};
#endif

static CEPW_SCRATCH pw_scratch = {1, MAXKL, 100, 5E-2, 0, 0, 10};

static MULTI *pk_array;
static MULTI *kappa0_array;
static MULTI *kappa1_array;
static MULTI *qk_array;

void uvip3p_(int *np, int *ndp, double *x, double *y, 
	     int *n, double *xi, double *yi);


CEPW_SCRATCH *GetCEPWScratch(void) {
  return &pw_scratch;
}

int SetCEQkMode(int m, double tol) {
  if (m == QK_DEFAULT) qk_mode = QK_EXACT;
  else qk_mode = m;
  if (tol > 0.0) qk_fit_tolerance = tol;
  return 0;
}

int SetCEEGridLimits(double min, double max, int type) {
  if (min <= 0) egrid_min = 0.05;
  else egrid_min = min;
  if (max <= 0) egrid_max = 8.0;
  else egrid_max = max;
  egrid_limits_type = type;

  return 0;
}

int SetCEEGridType(int type) {
  if (type >= 0) egrid_type = type;
  return 0;
}

int SetUsrCEEGridType(int type) {
  if (type >= 0) usr_egrid_type = type;
  return 0;
}

int SetCEPWGridType(int type) {
  if (type >= 0) pw_type = type;
  return 0;
}

int SetCETEGridDetail(int n, double *x) {
  n_tegrid = SetTEGridDetail(tegrid, log_te, n, x);
  return n_tegrid;
}

int SetCETEGrid(int n, double emin, double emax) {
  n_tegrid = SetTEGrid(tegrid, log_te, n, emin, emax);
  return n_tegrid;
}

int SetCEEGridDetail(int n, double *xg) {
  n_egrid = SetEGridDetail(egrid, log_egrid, n, xg);
  return n_egrid;
}

int SetCEEGrid(int n, double emin, double emax, double eth) {
  n_egrid = SetEGrid(egrid, log_egrid, n, emin, emax, eth);
  return n_egrid;
}

int SetUsrCEEGridDetail(int n, double *xg) {
  if (n > MAXNUSR) {
    printf("Max # of grid points reached \n");
    return -1;
  }
  n_usr = SetEGridDetail(usr_egrid, log_usr, n, xg);
  return n_usr;
}
 
int SetUsrCEEGrid(int n, double emin, double emax, double eth) {
  if (n > MAXNUSR) {
    printf("Max # of grid points reached \n");
    return -1;
  }
  n_usr = SetEGrid(usr_egrid, log_usr, n, emin, emax, eth);
  return n_usr;
}

int SetCEPWOptions(int qr, int max, int kl_cb, double tol) {
  pw_scratch.qr = qr;
  if (max > MAXKL) {
    printf("The maximum partial wave reached in Excitation: %d > %d\n", 
	   max, MAXKL);
    exit(1);
  }
  pw_scratch.max_kl = max;
  pw_scratch.kl_cb = kl_cb;
  pw_scratch.tolerance = tol;
  pw_scratch.nkl0 = 1;
  pw_scratch.kl[0] = 0;
  pw_scratch.log_kl[0] = -100.0;
  pw_scratch.nkl = 0;
  return 0;
}

int SetCEPWGrid(int ns, int *n, int *step) {
  if (pw_scratch.nkl0 <= 0) SetCEPWOptions(0, 256, 64, 5E-2);
  pw_scratch.nkl = SetPWGrid(&(pw_scratch.nkl0),
			     pw_scratch.kl,
			     pw_scratch.log_kl,
			     pw_scratch.max_kl,
			     &ns, n, step);
  pw_scratch.ns = ns;  
  return 0;
}

int CERadialPk(int *nkappa, int *nkl, double **pk, 
	       short **kappa0, short **kappa1, int ie,
	       int k0, int k1, int k) {
  int type, ko2, i, m, t, q;
  int kf0, kf1, kpp0, kpp1, km0, km1;
  int kl0, kl1, kl0p, kl1p;
  int j0, j1, kl_max, max0, j1min, j1max;
  ORBITAL *orb0, *orb1;
  int index[4];
  double te, e0, e1, z;
  double qkt, qkl, qkl0, h, a, r, rp;  
  int js1, js3, js[4], ks[4];
  double sd, se;
  int last_kl0, second_last_kl0;
  double **p;
  short **kp0, **kp1;
  double eps;

#ifdef PERFORM_STATISTICS
  clock_t start, stop;
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
    if (pw_type == 0) {
      *kappa0 = *kp0 + 1;
      *kappa1 = *kp1 + 1;
    } else {
      *kappa0 = *kp1 + 1;
      *kappa1 = *kp0 + 1;
    }
    *pk = *p + 1;
    return type;
  } 
  *nkappa = (MAXNKL+1)*(GetMaxRank()+1)*4;
  *kp0 = (short *) malloc(sizeof(short)*(*nkappa));
  *kp1 = (short *) malloc(sizeof(short)*(*nkappa));
  *p = (double *) malloc(sizeof(double)*((*nkappa)*n_tegrid));
  
  type = -1;
  orb0 = GetOrbital(k0);
  orb1 = GetOrbital(k1);
  kl0 = GetLFromKappa(orb0->kappa);
  kl1 = GetLFromKappa(orb1->kappa);
  kl0 = kl0/2;
  kl1 = kl1/2;
  if (IsEven(kl0 + kl1 + ko2)) {
    type = ko2;
  }

  if (egrid_type == 0) {
    e0 = egrid[ie];
    e1 = e0 - tegrid[0];
  } else {
    e1 = egrid[ie];
  }

  max0 = 2*Max(orb0->n, orb1->n);
  max0 *= sqrt(2.0-e1/Min(orb0->energy, orb1->energy));
  max0 = Max(16, max0);
  if (type >= 0 && type < CBMULTIPOLES) {
    kl_max = pw_scratch.kl_cb;
  } else {
    z = GetResidualZ();
    r = GetRMax();
    kl_max = r*sqrt(e1+2.0*z/r);
    kl_max = Min(kl_max, pw_scratch.max_kl);
  }

  js[0] = 0;
  ks[0] = k0;
  js[2] = 0;
  ks[2] = k1;		

  last_kl0 = 0;
  second_last_kl0 = 0;
  qkt = 0.0;
  q = 1;
  m = 1;
  *nkl = -1;
  eps = pw_scratch.tolerance;
  if (type >= CBMULTIPOLES) {
    z = GetCoulombBetheAsymptotic(tegrid[0], e1);
  }
  if (pw_type == 0) {
    js1 = 1;
    js3 = 3;
  } else {
    js1 = 3;
    js3 = 1;
  }

  for (t = 0; !last_kl0; t++) {
    kl0 = pw_scratch.kl[t];
    if (second_last_kl0) last_kl0 = 1;
    else {
      if (kl0 > max0) { 	
	if (type < 0) {
	  rp *= qkl;
	  rp = rp/qkt;
	  if (rp < eps) last_kl0 = 1;
	} else if (type >= CBMULTIPOLES) {
	  h = z*qkl;
	  h = h/(h + qkt);
	  if (h < eps) last_kl0 = 1;
	  else {
	    rp = fabs(1.0 - rp/z);
	    rp *= h;
	    if (rp < eps) last_kl0 = 1;
	  }
	} else {
	  z = (GetCoulombBethe(0, 0, ie, type, 1))[t-1];
	  h = z*qkl;
	  h = h/(h+qkt);
	  if (h < eps) last_kl0 = 1;
	  else {
	    z = (GetCoulombBethe(0, 0, ie, type, 0))[t-1];
	    z = z/(1.0-z);
	    rp = fabs(1.0 - rp/z);
	    rp *= h;
	    if (rp < eps) last_kl0 = 1;
	  }
	}
	if (pw_scratch.kl[t+1] > kl_max) {      
	  second_last_kl0 = 1;
	}
      }
    }
   
    qkl0 = qkl;
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
      if (pw_type == 1 && egrid_type == 1) {
	kf1 = OrbitalIndex(0, km0, e1);
	ks[3] = kf1;
      } else if (pw_type == 0 && egrid_type == 0) {
	kf1 = OrbitalIndex(0, km0, e0);
	ks[1] = kf1;
      }
      for (j1 = j1min; j1 <= j1max; j1 += 2) {
	for (kl1p = j1 - 1; kl1p <= j1 + 1; kl1p += 2) {	
	  kl1 = kl1p/2;
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
	  if (pw_type == 0 && egrid_type == 1) {
	    kf1 = OrbitalIndex(0, km1, e1);
	    ks[3] = kf1;
	  } else if (pw_type == 1 && egrid_type == 0) {
	    kf1 = OrbitalIndex(0, km1, e0);
	    ks[1] = kf1;
	  }
	  for (i = 0; i < n_tegrid; i++) {
	    te = tegrid[i];
	    if (egrid_type == 0) {
	      e1 = e0 - te;
	      if (pw_type == 0) {
		kf0 = OrbitalIndex(0, km1, e1);
		ks[3] = kf0;
	      } else {
		kf0 = OrbitalIndex(0, km0, e1);
		ks[3] = kf0;
	      }
	    } else {
	      e0 = e1 + te;
	      if (pw_type == 0) {
		kf0 = OrbitalIndex(0, km0, e0);
		ks[1] = kf0;
	      } else {
		kf0 = OrbitalIndex(0, km1, e0);
		ks[1] = kf0;	      
	      }
	    }
	    
	    if (kl1 >= pw_scratch.qr &&
		kl0 >= pw_scratch.qr) {
	      SlaterTotal(&sd, &se, js, ks, k, -1);
	    } else {
	      SlaterTotal(&sd, &se, js, ks, k, 1);
	    } 
	    r = sd + se;
	    if (i == 0 && (sd == 0.0 && se == 0.0)) break;
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
    if (t == 0) {
      qkt += qkl;
    } else if (!last_kl0 && !second_last_kl0) {
      if (qkl + 1.0 == 1.0) rp = 0.0;
      else rp = qkl/qkl0;
      h = kl0 - pw_scratch.kl[t-1];
      if (h > 1) {
	a = (1.0 - rp)*qkl0;
	rp = pow(rp, 1.0/h);
	rp = rp/(1.0-rp);
	a *= rp;
	qkt += a;
      } else {
	qkt += qkl;
	rp = rp/(1.0-rp);
      }
    }
  }
  *nkl += 1;
  *nkappa = m-1;
  *p = realloc(*p, sizeof(double)*q);
  (*p)[0] = type;
  *pk = *p + 1;
  *kp0 = realloc(*kp0, sizeof(short)*m);
  *kp1 = realloc(*kp1, sizeof(short)*m);
  (*kp0)[0] = *nkappa;
  (*kp1)[0] = *nkl;
  if (pw_type == 0) {
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

double *CERadialQkTable(int k0, int k1, int k2, int k3, int k) {
  int type, t, ie, ite, ipk, ipkp, nqk;
  int i, j, kl0, kl1, kl, nkappa, nkl, nkappap, nklp;
  short *kappa0, *kappa1, *kappa0p, *kappa1p, *tmp;
  double *pk, *pkp, *ptr, r, s, b;
  double *qk, *y2, rq[MAXNTE][MAXNE], e1, te;
  double *rqc, **p;
  int index[6];

#ifdef PERFORM_STATISTICS
  clock_t start, stop;
  start = clock();
#endif

  index[0] = 0;
  index[1] = k/2;
  index[2] = k0;
  index[3] = k1;
  index[4] = k2;
  index[5] = k3;
  
  
  p = (double **) MultiSet(qk_array, index, NULL);
  if (*p) {
    return *p;
  }   
  
  qk = pw_scratch.qk;
  y2 = pw_scratch.y2;

  nqk = n_tegrid*n_egrid;
  *p = (double *) malloc(sizeof(double)*(nqk+1));
  rqc = *p;

  for (ie = 0; ie < n_egrid; ie++) {
    e1 = egrid[ie];
    type = CERadialPk(&nkappa, &nkl, &pk, &kappa0, &kappa1, 
		      ie, k0, k1, k);
    nklp = nkl;
    if (k2 != k0 || k3 != k1) {
      type = CERadialPk(&nkappap, &nklp, &pkp, &kappa0p, &kappa1p,
			ie, k2, k3, k);
    }
    if (nkl > nklp) {
      ptr = pk;
      pk = pkp;
      pkp = ptr;
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
    } 
    nklp = nkl-1;
    for (ite = 0; ite < n_tegrid; ite++) {
      te = tegrid[ite];
      if (egrid_type == 0) e1 = egrid[ie] - te;      
      for (i = 0; i < nkl; i++) {
	qk[i] = 0.0;
      }
      t = -1;
      kl0 = -1;
      
      ipk = ite;
      for (i = 0; i < nkappa; i++) {
	if (pw_type == 0) kl = GetLFromKappa(kappa0[i]);
	else kl = GetLFromKappa(kappa1[i]);
	if (kl != kl0) {
	  t++;
	  kl0 = kl;
	}
	if (k2 == k0 && k3 == k1) {
	  s = pk[ipk]*pk[ipk];
	  qk[t] += s;
	} else {
	  s = 0.0;
	  ipkp = ite;
	  for (j = 0; j < nkappap; j++) {
	    if (kappa0[i] == kappa0p[j] && kappa1[i] == kappa1p[j]) {
	      s = pk[ipk]*pkp[ipkp];
	      break;
	    }
	    ipkp += n_tegrid;
	  }
	  qk[t] += s;
	}
	ipk += n_tegrid;
      }
      
      for (i = 0; i < nkl; i++) {
	if (!(qk[i] <= 0) && !(qk[i] > 0)) break;
      }
      nkl = i;
      nklp = nkl-1;
      if (nkl == 0) {
	r = 1E-31;
      } else if (nkl > 1) {
	spline(pw_scratch.log_kl, qk, nkl, 1E30, 1E30, y2); 
	r = qk[0];
	for (i = 1; i < nkl; i++) {
	  r += qk[i];
	  kl0 = pw_scratch.kl[i-1];
	  kl1 = pw_scratch.kl[i];
	  for (j = kl0+1; j < kl1; j++) {
	    splint(pw_scratch.log_kl, qk, y2, nkl, 
		 LnInteger(j), &s);
	    r += s;
	  }
	}
      } else {
	r = qk[0];
      }
      if (nkl > 0) {
	if (type >= CBMULTIPOLES) {
	  if (nkl > 10) {
	    b = qk[nklp]/qk[nklp-1];
	    if (b < 1.0 && b > 0.0) {
	      b = pow(b, 1.0/(pw_scratch.kl[nklp] - pw_scratch.kl[nklp-1]));
	      b = b/(1.0 - b);
	    } else {
	      b = GetCoulombBetheAsymptotic(te, e1);
	    }
	  } else {
	    b = GetCoulombBetheAsymptotic(te, e1);
	  }
	  s = qk[nklp]*b;
	  r = r + s;
	} else if (type >= 0) {
	  b = (GetCoulombBethe(0, ite, ie, type, 1))[nklp];
	  if (b > 0) {
	    s = qk[nklp]*b;
	    r = r + s;
 	  }
	}      
      }
      rq[ite][ie] = r;
    }
  }

  ptr = rqc;
  for (ite = 0; ite < n_tegrid; ite++) {
    for (ie = 0; ie < n_egrid; ie++) {
      ptr[ie] = rq[ite][ie];
    }
    ptr += n_egrid;
  }
  rqc[nqk] = type;

#ifdef PERFORM_STATISTICS
  stop = clock();
  timing.rad_qk += stop-start;
#endif
  return rqc;
}

  
double *CERadialQkMSubTable(int k0, int k1, int k2, int k3, 
			    int k, int kp, int nq, int *q) {
  int type1, type2, kl, nqk;
  int i, j, kl0, klp0, kl0_2, klp0_2, kl1;
  int nkappa, nkappap, nkl, nklp;
  short  *kappa0, *kappa1, *kappa0p, *kappa1p;
  double *pk, *pkp, *ptr, *qy2;
  int km0, km1, j0, jp0, j1, kmp0, kmp1, km0_m, kmp0_m;
  int mi, mf, t, c0, cp0;
  double r, e0, e1, te, s, b;
  double pha0, phap0;
  double s3j1, s3j2, s3j3, s3j4;
  int ie, ite, iq, ipk, ipkp;
  double qk[MAXMSUB][MAXNKL];
  double rq[MAXMSUB][MAXNTE][MAXNE];
  double rqt[MAXMSUB];
  double *rqc, **p;
  int index[6];

#ifdef PERFORM_STATISTICS
  clock_t start, stop;
  start = clock();
#endif

  index[0] = kp/2;
  index[1] = k/2;
  index[2] = k0;
  index[3] = k1;
  index[4] = k2;
  index[5] = k3;
  
  p = (double **) MultiSet(qk_array, index, NULL);
  if (*p) {
    return *p;
  } 
  
  qy2 = pw_scratch.y2;
  pkp = NULL;
  nqk = nq*n_tegrid*n_egrid;
  *p = (double *) malloc(sizeof(double)*(nqk+1));
  rqc = *p;

  for (ie = 0; ie < n_egrid; ie++) {
    type1 = CERadialPk(&nkappa, &nkl, &pk, &kappa0, &kappa1,
		       ie, k0, k1, k);
    if (kp == k && k2 == k0 && k3 == k1) {
      pkp = pk;
      nkappap = nkappa;
      nklp = nkl;
      kappa0p = kappa0;
      kappa1p = kappa1;
      type2 = type1;
    } else {
      type2 = CERadialPk(&nkappap, &nklp, &pkp, &kappa0p, &kappa1p,
			 ie, k2, k3, kp);
      if (nklp < nkl) nkl = nklp;
    }

    for (ite = 0; ite < n_tegrid; ite++) {
      te = tegrid[ite];
      if (egrid_type == 0) {
	e0 = egrid[ie];
	e1 = e0 - te;
      } else {
	e1 = egrid[ie];
	e0 = e1 + te;
      }
      
      for (i = 0; i < nq; i++) { 
	for (j = 0; j < nkl; j++) { 
	  qk[i][j] = 0.0; 
	}   
      } 
      
      kl = -1;
      i = -1;
      ipk = ite;
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
	ipkp = ite;
	for (t = 0; t < nkappap; t++) {
	  kmp0 = kappa0p[t];
	  kmp1 = kappa1p[t];
	  if (kmp1 != km1) {
	    ipkp += n_tegrid;
	    continue;
	  }
	  GetJLFromKappa(kmp0, &jp0, &klp0);
	  klp0_2 = klp0/2;
	  
	  s = pk[ipk]*pkp[ipkp];
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
	    c0 = OrbitalIndex(0, km0_m, e0); 
	    cp0 = OrbitalIndex(0, kmp0_m, e0);
	    pha0 = GetPhaseShift(c0); 
	    phap0 = GetPhaseShift(cp0);  
	    r = pha0 - phap0;
	    s *= cos(r);
	  }
	
	  for (iq = 0; iq < nq; iq++) { 
	    rqt[iq] = 0.0; 
	  } 
	  for (mi = -1; mi <= 1; mi += 2) { 
	    s3j1 = W3j(j0, 1, kl0, -mi, mi, 0); 
	    s3j2 = W3j(jp0, 1, klp0, -mi, mi, 0); 
	    for (iq = 0; iq < nq; iq++) { 
	      mf = mi + q[iq]; 
	      s3j3 = W3j(j0, k, j1, -mi, -q[iq], mf); 
	      s3j4 = W3j(jp0, kp, j1, -mi, -q[iq], mf); 
	      rqt[iq] += s3j1*s3j2*s3j3*s3j4; 
	    } 
	  } 
	  
	  for (iq = 0; iq < nq; iq++) { 
	    qk[iq][i] += s*rqt[iq]; 
	  } 

	  ipkp += n_tegrid;
	}
	ipk += n_tegrid;
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
	rq[iq][ite][ie] = r; 
      } 
      
      i = nkl - 1;
      if (type1 == type2) {
	if (type1 >= CBMULTIPOLES) {
	  for (iq = 0; iq < nq; iq++) {
	    if (nkl > 10) {
	      b = qk[iq][i]/qk[iq][i-1];
	      if (b < 1.0 && b > 0.0) {
		b = pow(b, 1.0/(pw_scratch.kl[i] - pw_scratch.kl[i-1]));
		b = b/(1.0 - b);
	      } else {
		b = GetCoulombBetheAsymptotic(te, e1);
	      }
	    } else {
	      b = GetCoulombBetheAsymptotic(te, e1);
	    }
	    s = qk[iq][i]*b;
	    rq[iq][ite][ie] += s;
	  }
	} else if (type1 >= 0) {
	  for (iq = 0; iq < nq; iq++) {
	    b = (GetCoulombBethe(0, ite, ie, type1, iq+1))[i];
	    if (b > 0) {
	      s = qk[iq][i]*b;
	      rq[iq][ite][ie] += s;
	    }
	  }
	}
      } else if (type1 >= 0) {
	for (iq = 0; iq < nq; iq++) {
	  if (nkl > 10) {
	    b = qk[iq][i]/qk[iq][i-1];
	    if (b < 1.0 && b > 0.0) {
	      b = pow(b, 1.0/(pw_scratch.kl[i] - pw_scratch.kl[i-1]));
	      b = b/(1.0 - b);
	    } else {
	      b = GetCoulombBetheAsymptotic(te, e1);
	    }
	  } else {
	    b = GetCoulombBetheAsymptotic(te, e1);
	  }
	  s = qk[iq][i]*b;
	  rq[iq][ite][ie] += s;
	}
      }
    }
  }

  ptr = rqc;
  for (iq = 0; iq < nq; iq++) {
    for (ite = 0; ite < n_tegrid; ite++) {
      for (ie = 0; ie < n_egrid; ie++) {
	ptr[ie] = rq[iq][ite][ie];
      }
      ptr += n_egrid;
    }
  }    
  rqc[nqk] = type1;
  if (type2 != 1) rqc[nqk] = type2;

#ifdef PERFORM_STATISTICS
  stop = clock();
  timing.rad_qk += stop-start;
#endif
  
  return rqc;
} 	
  
int CERadialQk(double *rqc, double te, int k0, int k1, int k2, int k3, int k) {
  int i, np, nd, type;
  int j, m;
  double *rqe, rq[MAXNTE];
  double *xte, x0;
  
  rqe = CERadialQkTable(k0, k1, k2, k3, k);
  if (n_tegrid == 1) {
    for (i = 0; i < n_egrid; i++) {
      rqc[i] = rqe[i];
    }
    type = rqe[i];
  } else {
    np = 3;
    nd = 1;
    type = rqe[n_tegrid*n_egrid];
    if (type == 0 || type == 1) {
      xte = log_te;
      x0 = log(te);
    } else {
      xte = tegrid;
      x0 = te;
    }
    for (i = 0; i < n_egrid; i++) {
      j = i;
      for (m = 0; m < n_tegrid; m++) {
	rq[m] = rqe[j];
	if (type == 0) {
	  rq[m] = log(fabs(rq[m]));
	}
	j += n_egrid;
      }
      uvip3p_(&np, &n_tegrid, xte, rq, &nd, &x0, &rqc[i]);
      if (type == 0) {
	rqc[i] = exp(rqc[i]);
	if (rqe[i] < 0.0) rqc[i] = -rqc[i];
      }
    }
  }

  return type;
}

int CERadialQkMSub(double *rqc, double te, int k0, int k1, int k2, int k3, 
		   int k, int kp, int nq, int *q) {
  int i, np, nd, iq, n;
  int j, m, type;
  double *rqe, rq[MAXNTE];
  double *xte, x0;
  
  rqe = CERadialQkMSubTable(k0, k1, k2, k3, k, kp, nq, q);
  if (n_tegrid == 1) {
    for (iq = 0; iq < nq; iq++) {
      for (i = 0; i < n_egrid; i++) {
	rqc[i] = rqe[i];
      }
      rqc += n_egrid;
      rqe += n_egrid;
    }
    type = rqe[0];
  } else {
    np = 3;
    nd = 1;
    n = n_tegrid*n_egrid;
    type = rqe[nq*n];
    if (type == 1) {
      xte = log_te;
      x0 = log(te);
    } else {
      xte = tegrid;
      x0 = te;
    }
    for (iq = 0; iq < nq; iq++) {
      for (i = 0; i < n_egrid; i++) {
	j = i;
	for (m = 0; m < n_tegrid; m++) {
	  rq[m] = rqe[j];
	  j += n_egrid;
	}
	uvip3p_(&np, &n_tegrid, xte, rq, &nd, &x0, &rqc[i]);
      }
      rqe += n;
      rqc += n_egrid;
    }
  }  
  return type;
}

void CERadialQkFromFit(int np, double *p, int n, double *x, double *logx,
		       double *y, double *dy, int ndy, void *extra) {
  double a, b, c, d, D;
  int i, k;

  D = *((double *) extra);
  if (D >= 0.0) {
    if (ndy <= 0) {
      for (i = 0; i < n; i++) {
	a = p[1]*p[1];
	b = 1.0 - x[i];
	c = pow(x[i], a);
	d = pow(b, p[3]*p[3]);
	y[i] = p[0]*c + p[2]*d;
      }
    } else {
      for (i = 0; i < n; i++) {
	a = p[1]*p[1];
	b = 1.0 - x[i];
	c = pow(x[i], a);
	d = pow(b, p[3]*p[3]);
	k = i;
	dy[k] = c;
	k += ndy;
	c *= p[0]*logx[i];
	dy[k] = 2.0*c*p[1];
	k += ndy;
	dy[k] = d;
	k += ndy;
	d *= p[2]*log(b);
	dy[k] = 2.0*d*p[3];
      }
    }
  } else {
    if (ndy <= 0) {
      for (i = 0; i < n; i++) {
	a = 1.0/(x[i]+p[3]);
	b = a*a;
	c = -2.0 + p[1]*a + p[2]*b;
	y[i] = p[0]*pow(x[i], c);
      }
    } else {
      for (i = 0; i < n; i++) {
	a = 1.0/(p[3] + x[i]);
	b = a*a;
	c = -2.0 + p[1]*a + p[2]*b;
	d = pow(x[i], c);
	k = i;
	dy[k] = d;
	k += ndy;
	d = p[0]*d*logx[i];
	dy[k] = d*a;
	k += ndy;
	dy[k] = d*b;
	k += ndy;
	dy[k] = -d*(p[1]*b + 2.0*p[2]*b*a);
      }
    }
  }
}

int CollisionStrength(double *qkt, double *params, double *e, double *bethe,
		      int lower, int upper, int msub) {
  int i, j, t, h, p, m, type, ty;  
  LEVEL *lev1, *lev2;
  double te, c, r, s3j;
  ANGULAR_ZMIX *ang;
  int nz, j1, j2, ie, np;
  int nq, q[MAXMSUB];
  double rq[MAXMSUB*(MAXNE+1)], qkc[MAXMSUB*(MAXNE+1)];
  double *rqk, tol;
  double c1, c2, *g1, *g2;
  double x1, x2, x3, x1s, x2s, x3s;
  int ierr, ipvt[NPARAMS];
  int lwa=5*NPARAMS+MAXNE;
  double wa[5*NPARAMS+MAXNE];
  double fvec[MAXNE], fjac[MAXNE*NPARAMS];

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
    rqk = qkc;
    for (t = -j1; t <= 0; t += 2) {
      for (h = -j2; h <= j2; h += 2) {
	for (ie = 0; ie < n_egrid; ie++) {
	  rqk[ie] = 0.0;
	}
	rqk += n_egrid;
      }
    }
  } else {
    rqk = qkc;
    for (ie = 0; ie < n_egrid; ie++) {
      rqk[ie] = 0.0;
    }
  }
 
  nz = AngularZMix(&ang, lower, upper, -1, -1);
  if (nz <= 0) return -1;
  type = -1;
  for (i = 0; i < nz; i++) {
    for (j = i; j < nz; j++) {
      c = ang[i].coeff * ang[j].coeff;
      if (i != j) c *= 2.0;
      if (!msub) {
	if (ang[i].k != ang[j].k) continue;
	c /= ang[i].k + 1.0;
	ty = CERadialQk(rq, te, ang[i].k0, ang[i].k1,
			ang[j].k0, ang[j].k1, ang[i].k);
	if (ty > type) type = ty;	  
	for (ie = 0; ie < n_egrid; ie++) {
	  qkc[ie] += c*rq[ie];
	}
      } else {
	ty = CERadialQkMSub(rq, te, ang[i].k0, ang[i].k1,
			    ang[j].k0, ang[j].k1,
			    ang[i].k, ang[j].k, nq, q);
	if (ty > type) type = ty;
	rqk = qkc;
	for (t = -j1; t <= 0; t += 2) {
	  for (h = -j2; h <= j2; h += 2) {
	    m = (-abs(t-h)-q[0])/2;
	    m *= n_egrid;
	    s3j = W3j(j1, ang[i].k, j2, -t, t-h, h);
	    if (ang[j].k != ang[i].k) {
	      s3j *= W3j(j1, ang[j].k, j2, -t, t-h, h);
	    } else {
	      s3j *= s3j;
	    }
	    for (ie = 0; ie < n_egrid; ie++) {
	      rqk[ie] += c*rq[m+ie]*s3j;
	    }
	    rqk += n_egrid;
	  }
	}
      }
    }
  }

  for (t = 0; t < NGOSK; t++) {
    gos[t] = 0.0;
  }
  r = 0;
  if (type >= 0) {
    for (i = 0; i < nz; i++) {
      p = GetOrbital(ang[i].k0)->kappa;
      p = GetJFromKappa(p);
      m = GetOrbital(ang[i].k1)->kappa;
      m = GetJFromKappa(m);
      c1 = ReducedCL(p, ang[i].k, m);
      if (c1 == 0) continue;
      g1 = GeneralizedMoments(NGOSK, kgrid, 
			      ang[i].k0, ang[i].k1, ang[i].k/2);
      for (j = i; j < nz; j++) {
	if (ang[j].k != ang[i].k) continue;
	p = GetOrbital(ang[j].k0)->kappa;
	p = GetJFromKappa(p);
	m = GetOrbital(ang[j].k1)->kappa;
	m = GetJFromKappa(m);
	c2 = ReducedCL(p, ang[j].k, m);
	if (c2 == 0) continue;
	g2 = GeneralizedMoments(NGOSK, kgrid, 
				ang[j].k0, ang[j].k1, ang[j].k/2);
	c = ang[i].coeff*ang[j].coeff;
	if (i != j) c *= 2.0;
	c *= 2.0*c1*c2;
	if (ang[i].k == 2) {
	  r += c * (RadialMoments(1, ang[i].k0, ang[i].k1)*
		    RadialMoments(1, ang[j].k0, ang[j].k1));
	}
	c *= ang[i].k + 1.0;
	for (t = 0; t < NGOSK; t++) {
	  gos[t] += (c/(kgrid[t]*kgrid[t]))*g1[t]*g2[t];
	}
      }
    }
    r /= 3.0;
    bethe[0] = r*2.0;
    for (t = 0; t < NGOSK; t++) {
      gos[t] = log(gos[t]);
    }
    for (i = 0; i < n_born; i++) {
      c1 = sqrt(2.0*te*born_egrid[i]);
      c2 = sqrt(2.0*te*(born_egrid[i]-1.0));
      kint[0] = c1 - c2;
      kint[NGOSK-1] = c1 + c2;
      log_kint[0] = log(kint[0]);
      log_kint[NGOSK-1] = log(kint[NGOSK-1]);
      c1 = (log_kint[NGOSK-1] - log_kint[0])/(NGOSK-1);
      for (t = 1; t < NGOSK-1; t++) {
	log_kint[t] = log_kint[t-1] + c1;
	kint[t] = exp(log_kint[t]);
      }
      np = 3;
      j = NGOSK;
      uvip3p_(&np, &j, log_kgrid, gos, &j, log_kint, gosint);
      for (t = 0; t < NGOSK; t++) {
	gosint[t] = exp(gosint[t]);
      }
      fvec[i] = Simpson(gosint, 0, NGOSK-1);
      fvec[i] *= 4.0*c1;
      if (bethe[0]) {
	fvec[i] -= bethe[0]*log(born_egrid[i]);
      }
    }
    x1 = 1.0/born_egrid[0];
    x2 = 1.0/born_egrid[1];
    x3 = 1.0/born_egrid[2];
    x1s = x1*x1;
    x2s = x2*x2;
    x3s = x3*x3;
    bethe[2] = (((fvec[0]*x2s - fvec[1]*x1s)*(x3s - x1s) -
		 (fvec[0]*x3s - fvec[2]*x1s)*(x2s - x1s))/
		((x1*x2s - x2*x1s)*(x3s - x1s) - 
		 (x1*x3s - x3*x1s)*(x2s - x1s)));
    bethe[1] = fvec[0]*x2s - fvec[1]*x1s - bethe[2]*(x1*x2s - x2*x1s);
    bethe[1] /= (x2s - x1s);
    bethe[3] = (fvec[0] - bethe[1] - bethe[2]*x1)/x1s;
  } else {
    bethe[0] = -1.0;
    bethe[1] = 0.0;
    bethe[2] = 0.0;
    bethe[3] = 0.0;
  }

  free(ang);
  
  /* there is a factor of 4 coming from normalization and the 2 
     from the formula */
  if (!msub) {
    if (qk_mode == QK_FIT) {
      for (ie = 0; ie < n_egrid; ie++) {
	qkc[ie] = 8.0*qkc[ie];
	qkt[ie] = qkc[ie];
	xusr[ie] = egrid[ie]/te;
	if (egrid_type == 1) xusr[ie] += 1.0;
	log_xusr[ie] = log(xusr[ie]);
      }	
      if (qkt[0] < EPS16 && qkt[n_egrid-1] < EPS16) {
	params[0] = 0.0;
	params[1] = 0.0;
	params[2] = 0.0;
	params[3] = 0.0;
      } else {
	tol = qk_fit_tolerance;
	if (*bethe < 0.0) {
	  params[0] = qkc[0];
	  params[1] = 0.0;
	  params[2] = 0.0;
	  params[3] = 0.0;
	} else {
	  for (ie = 0; ie < n_egrid; ie++) {
	    qkc[ie] -= (*bethe)*log_xusr[ie];
	    xusr[ie] = 1.0/xusr[ie];
	    log_xusr[ie] = -log_xusr[ie];
	  }
	  params[0] = qkc[0];
	  params[1] = 1.0;
	  params[2] = qkc[n_egrid-1];
	  params[3] = 0.1;
	}
	np = NPARAMS;
	ierr = NLSQFit(np, params, tol, ipvt, fvec, fjac, MAXNE, wa, lwa,
		       n_egrid, xusr, log_xusr, qkc, qkt, 
		       CERadialQkFromFit, (void *) bethe);
	if (ierr > 3) params[0] = ierr*1E30;
	if (*bethe >= 0.0) {
	  params[1] *= params[1];
	  params[3] *= params[3];
	}
      }
    } else if (qk_mode == QK_INTERPOLATE) {
      np = 3;
      if (type != 1) {
	for (ie = 0; ie < n_egrid; ie++) {
	  qkc[ie] *= 8.0;
	  qkc[ie] = log(qkc[ie]);
	}
	uvip3p_(&np, &n_egrid, log_egrid, qkc, &n_usr, log_usr, qkt);
	for (ie = 0; ie < n_usr; ie++) {
	  qkt[ie] = exp(qkt[ie]);
	}
      } else {
	for (ie = 0; ie < n_egrid; ie++) {
	  qkc[ie] *= 8.0;
	}
	uvip3p_(&np, &n_egrid, log_egrid, qkc, &n_usr, log_usr, qkt);
      }
    } else if (qk_mode == QK_EXACT) {
      for (ie = 0; ie < n_usr; ie++) {
	qkt[ie] = 8.0*qkc[ie];
      }
    }
    return 1;
  } else {
    rqk = qkc;
    p = 0;
    if (qk_mode == QK_FIT) {
      for (t = -j1; t <= 0; t += 2) {
	for (h = -j2; h <= j2; h += 2) {
	  for (ie = 0; ie < n_egrid; ie++) {
	    qkt[ie] = 8.0*rqk[ie];
	  }	
	  if (qkt[0] < EPS10 && qkt[n_egrid-1] < EPS10) {
	    continue;
	  }
	  p++;
	  rqk += n_egrid;
	  qkt += n_usr;
	}
      }
    } else if (qk_mode == QK_INTERPOLATE) {
      np = 3;
      if (type != 1) {	
	for (t = -j1; t <= 0; t += 2) {
	  for (h = -j2; h <= j2; h += 2) {
	    for (ie = 0; ie < n_egrid; ie++) {
	      rqk[ie] *= 8.0;
	      rqk[ie] = log(rqk[ie]);
	    }
	    uvip3p_(&np, &n_egrid, log_egrid, rqk, &n_usr, log_usr, qkt);
	    for (ie = 0; ie < n_usr; ie++) {
	      qkt[ie] = exp(qkt[ie]);
	    }
	    p++;
	    rqk += n_egrid;
	    qkt += n_usr;
	  }
	} 
      } else {
	for (t = -j1; t <= 0; t += 2) {
	  for (h = -j2; h <= j2; h += 2) {
	    for (ie = 0; ie < n_egrid; ie++) {
	      rqk[ie] *= 8.0;
	    }
	    uvip3p_(&np, &n_egrid, log_egrid, rqk, &n_usr, log_usr, qkt);
	    p++;
	    rqk += n_egrid;
	    qkt += n_usr;
	  }
	}
      }
    } else if (qk_mode == QK_EXACT) {
      for (t = -j1; t <= 0; t += 2) {
	for (h = -j2; h <= j2; h += 2) {	
	  for (ie = 0; ie < n_egrid; ie++) {
	    qkt[ie] = 8.0*rqk[ie];
	  }
	  p++;
	  rqk += n_egrid;
	  qkt += n_usr;
	}
      }
    }  
    return p;
  }
}

int SaveExcitation(int nlow, int *low, int nup, int *up, int msub, char *fn) {
#ifdef PERFORM_STATISTICS
  STRUCT_TIMING structt;
  ANGULAR_TIMING angt;
  RECOUPLE_TIMING recouplet;
  RAD_TIMING radt;
#endif  
  int i, j, k, n, m, ie, ip;
  FILE *f;
  double qkc[MAXMSUB*MAXNUSR];
  double params[MAXMSUB*NPARAMS];
  int *alev;
  int nsub;
  LEVEL *lev1, *lev2;
  CE_RECORD r;
  CE_HEADER ce_hdr;
  F_HEADER fhdr;
  ARRAY subte;
  int isub, n_tegrid0, n_egrid0, n_usr0;
  int te_set, e_set, usr_set;
  double emin, emax, e, c, e0, e1, te0;
  double rmin, rmax, bethe[4];

  n = 0;
  alev = NULL;
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
  if (m == 0) {
    return 0;
  }

  if (tegrid[0] < 0) {
    te_set = 0;
  } else {
    te_set = 1;
  }
  if (egrid[0] < 0) {
    e_set = 0;
  } else {
    e_set = 1;
  }
  if (usr_egrid[0] < 0) {
    usr_set = 0;
  } else {
    usr_set = 1;
  }
  n_tegrid0 = n_tegrid;
  n_egrid0 = n_egrid;
  n_usr0 = n_usr;

  ArrayInit(&subte, sizeof(double), 128);
  ArrayAppend(&subte, &emin);
  c = 1.0/TE_MIN_MAX;
  if (!e_set || !te_set) {
    e = c*emin;
    while (e < emax) {
      ArrayAppend(&subte, &e);
      e *= c;
    }
  }
  ArrayAppend(&subte, &emax);
 
  if (msub) {
    pw_type = 1;
    qk_mode = QK_EXACT;
  } else {
    if (pw_type < 0) pw_type = 0;
  }
  egrid_type = 1;
  usr_egrid_type = 1;
  if (pw_scratch.nkl == 0) {
    SetCEPWGrid(0, NULL, NULL);
  }

  e = (emin + emax)*0.5;
  if (egrid_limits_type == 0) {
    rmin = egrid_min;
    rmax = egrid_max;
  } else {
    rmin = egrid_min/e;
    rmax = egrid_max/e;
  }
  te0 = emax;

  born_egrid[0] = 15.0;
  born_egrid[1] = 30.0;
  born_egrid[2] = 60.0;
  n_born = 3;

  e0 = emin;
  fhdr.type = DB_CE;
  strcpy(fhdr.symbol, GetAtomicSymbol());
  fhdr.atom = GetAtomicNumber();
  f = OpenFile(fn, &fhdr);
  for (isub = 1; isub < subte.dim; isub++) {
    e1 = *((double *) ArrayGet(&subte, isub));
    if (isub == subte.dim-1) e1 = e1*1.001;
    if (!te_set) {
      emin = e0;
      emax = e1;
      if (emin < 1.0/HARTREE_EV) {
	emin = 1.0/HARTREE_EV;
	emax = 3.0*emin;
	if (emax < e1) emax = e1;
      }
      e = emax/emin;  
      if (e < 1.1) {
	SetCETEGrid(1, 0.5*(emax+emin), emax);
      } else if (e < 2.0) {
	SetCETEGrid(2, emin, emax);
      } else {
	if (n_tegrid0 == 0) {
	  n_tegrid = 3;
	}
	SetCETEGrid(n_tegrid, emin, emax);
      }
    }

    c = 10.0*e1;
    if (te0 > c) ce_hdr.te0 = c;
    else ce_hdr.te0 = te0;
    emin = rmin*ce_hdr.te0;
    emax = rmax*ce_hdr.te0;
    
    if (qk_mode == QK_EXACT) {
      if (n_egrid0 <= 0) {
	if (n_usr0 <= 0) n_usr = 6;
	if (!usr_set) {
	  SetUsrCEEGrid(n_usr, emin, emax, ce_hdr.te0);
	  usr_egrid_type = 1;
	}  
	SetCEEGridDetail(n_usr, usr_egrid);
      } else {
	if (!e_set) {
	  SetCEEGrid(n_egrid, emin, emax, ce_hdr.te0);
	  usr_egrid_type = 1;
	}
	SetUsrCEEGridDetail(n_egrid, egrid);
      }
    } else {
      if (n_egrid0 == 0) {
	n_egrid = 6;
      }
      if (!e_set) {
	SetCEEGrid(n_egrid, emin, emax, ce_hdr.te0);
      }
      if (qk_mode == QK_INTERPOLATE) {
	if (n_usr0 <= 0) {
	  SetUsrCEEGridDetail(n_egrid, egrid);
	  usr_egrid_type = 1;
	} else if (!usr_set) {
	  SetUsrCEEGrid(n_usr, emin, emax, ce_hdr.te0);
	  usr_egrid_type = 1;
	}
      } else if (qk_mode == QK_FIT) {
	SetUsrCEEGridDetail(n_egrid, egrid);
	usr_egrid_type = 1;
      }
    }
    if (qk_mode == QK_FIT && n_egrid <= NPARAMS) {
      printf("n_egrid must > %d to use QK_FIT mode\n", NPARAMS);
      return -1;
    }
    if (qk_mode == QK_INTERPOLATE) {
      for (i = 0; i < n_egrid; i++) {
	log_egrid[i] = egrid[i];
	if (egrid_type == 1) log_egrid[i] += ce_hdr.te0;
	log_egrid[i] = log(log_egrid[i]);
      }
      for (i = 0; i < n_egrid; i++) {
	log_usr[i] = usr_egrid[i];
	if (usr_egrid_type == 1) log_usr[i] += ce_hdr.te0;
	log_usr[i] = log(log_usr[i]);
      }
    }

    emax = sqrt(2.0*e1*born_egrid[2]);
    emin = sqrt(2.0*(e1*born_egrid[2] - e1));
    rmin = emax - emin;
    rmax = emax + emin;
    kgrid[0] = rmin;
    kgrid[NGOSK-1] = rmax;
    log_kgrid[0] = log(rmin);
    log_kgrid[NGOSK-1] = log(rmax);
    e = (log_kgrid[NGOSK-1] - log_kgrid[0])/(NGOSK-1);
    for (i = 1; i < NGOSK-1; i++) {
      log_kgrid[i] = log_kgrid[i-1] + e;
      kgrid[i] = exp(log_kgrid[i]);
    }

    e = 0.0;
    c = GetResidualZ();
    PrepCoulombBethe(1, n_tegrid, n_egrid, c, &e, tegrid, egrid,
		     pw_scratch.nkl, pw_scratch.kl, egrid_type, 
		     pw_type, msub);

    ce_hdr.nele = GetNumElectrons(low[0]);
    ce_hdr.qk_mode = qk_mode;
    if (qk_mode == QK_FIT) 
      ce_hdr.nparams = NPARAMS;
    else
      ce_hdr.nparams = 0;
    ce_hdr.pw_type = pw_type;
    ce_hdr.n_tegrid = n_tegrid;
    ce_hdr.n_egrid = n_egrid;
    ce_hdr.egrid_type = egrid_type;
    ce_hdr.n_usr = n_usr;
    ce_hdr.usr_egrid_type = usr_egrid_type;
    ce_hdr.msub = msub;
    ce_hdr.tegrid = tegrid;
    ce_hdr.egrid = egrid;
    ce_hdr.usr_egrid = usr_egrid;

    InitFile(f, &fhdr, &ce_hdr);  
    nsub = 1;
    if (qk_mode == QK_FIT) {
      m = ce_hdr.nparams * nsub;
      r.params = (float *) malloc(sizeof(float)*m);
    }
    m = ce_hdr.n_usr * nsub;
    r.strength = (float *) malloc(sizeof(float)*m);
    
    for (i = 0; i < nlow; i++) {
      lev1 = GetLevel(low[i]);
      for (j = 0; j < nup; j++) {
	lev2 = GetLevel(up[j]);
	e = lev2->energy - lev1->energy;
	if (e < e0 || e >= e1) continue;
	k = CollisionStrength(qkc, params, &e, bethe, low[i], up[j], msub); 
	if (k < 0) continue;
	for (m = 0; m < 4; m++) {
	  r.bethe[m] = bethe[m];
	}
	r.lower = low[i];
	r.upper = up[j];
	r.nsub = k;
	if (r.nsub > nsub) {
	  if (qk_mode == QK_FIT) {
	    m = ce_hdr.nparams * r.nsub;
	    r.params = (float *) realloc(r.params, sizeof(float)*m);
	  }
	  m = ce_hdr.n_usr * r.nsub;
	  r.strength = (float *) realloc(r.strength, sizeof(float)*m);
	  nsub = r.nsub;
	}

	if (qk_mode == QK_FIT) {
	  ip = 0;
	  for (m = 0; m < r.nsub; m++) {
	    for (ie = 0; ie < ce_hdr.nparams; ie++) {
	      r.params[ip] = (float) params[ip];
	      ip++;
	    }
	  }
	}
      
	ip = 0;
	for (m = 0; m < r.nsub; m++) {
	  for (ie = 0; ie < ce_hdr.n_usr; ie++) {
	    r.strength[ip] = (float) qkc[ip];
	    ip++;
	  }
	}
	WriteCERecord(f, &r);
      }
    }
    if (qk_mode == QK_FIT) free(r.params);
    free(r.strength);
    DeinitFile(f, &fhdr);
    e0 = e1;
    FreeExcitationQk();
    ReinitRadial(1);
  }

  ArrayFree(&subte, NULL);
  if (alev) free(alev);
  CloseFile(f, &fhdr);

#ifdef PERFORM_STATISTICS
  GetStructTiming(&structt);
  fprintf(perform_log, "AngZMix: %6.1E, AngZFB: %6.1E, AngZxZFB: %6.1E, SetH: %6.1E DiagH: %6.1E\n",
	  ((double) (structt.angz_mix))/CLOCKS_PER_SEC,
	  ((double) (structt.angz_fb))/CLOCKS_PER_SEC,
	  ((double) (structt.angzxz_fb))/CLOCKS_PER_SEC,
	  ((double) (structt.set_ham))/CLOCKS_PER_SEC,
	  ((double) (structt.diag_ham))/CLOCKS_PER_SEC);
  fprintf(perform_log, "AngZS: %6.1E, AngZFBS: %6.1E, AngZxZFBS: %6.1E, AddZ: %6.1E, AddZxZ: %6.1E\n",
	  ((double) (structt.angz_states))/CLOCKS_PER_SEC,
	  ((double) (structt.angzfb_states))/CLOCKS_PER_SEC,
	  ((double) (structt.angzxzfb_states))/CLOCKS_PER_SEC,
	  ((double) (structt.add_angz))/CLOCKS_PER_SEC,
	  ((double) (structt.add_angzxz))/CLOCKS_PER_SEC);

  GetAngularTiming(&angt);
  fprintf(perform_log, "W3J: %6.1E, W6J: %6.1E, W9J: %6.1E\n", 
	  ((double)angt.w3j)/CLOCKS_PER_SEC, 
	  ((double)angt.w6j)/CLOCKS_PER_SEC, 
	  ((double)angt.w9j)/CLOCKS_PER_SEC);
  GetRecoupleTiming(&recouplet);
  fprintf(perform_log, "AngZ: %6.1E, AngZxZ: %6.1E, Interact: %6.1E\n",
	  ((double)recouplet.angz)/CLOCKS_PER_SEC,
	  ((double)recouplet.angzxz)/CLOCKS_PER_SEC,
	  ((double)recouplet.interact)/CLOCKS_PER_SEC);
  GetRadTiming(&radt);
  fprintf(perform_log, "Dirac: %d, %6.1E, 1E: %6.1E, Slater: %6.1E, 2E: %6.1E\n", 
	  GetNumContinua(),
	  ((double)radt.dirac)/CLOCKS_PER_SEC, 
	  ((double)radt.radial_1e)/CLOCKS_PER_SEC,
	  ((double)radt.radial_slater)/CLOCKS_PER_SEC,
	  ((double)radt.radial_2e)/CLOCKS_PER_SEC);
  fprintf(perform_log, "RadPk: %6.1E, SetKappa: %6.1E, RadQk: %6.1E\n",
	  ((double)timing.rad_pk)/CLOCKS_PER_SEC, 
	  ((double)timing.set_kappa)/CLOCKS_PER_SEC, 
	  ((double)timing.rad_qk)/CLOCKS_PER_SEC);

  fprintf(perform_log, "\n");
#endif /* PERFORM_STATISTICS */

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

int FreeExcitationQk(void) {
  ARRAY *b;
  b = qk_array->array;
  if (b == NULL) return 0;
  MultiFreeData(b, qk_array->ndim, _FreeExcitationPk);
  FreeExcitationPk(-1);
  
  return 0;
}
  
int InitExcitation(void) {
  int blocks1[] = {4, 8, 8, 4};
  int blocks2[] = {4, 4, 4, 4, 4, 4};
  int ndim;

  ndim = 4;
  pk_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(pk_array, sizeof(double *), ndim, blocks1);
  kappa0_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(kappa0_array, sizeof(short *), ndim, blocks1);
  kappa1_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(kappa1_array, sizeof(short *), ndim, blocks1);

  ndim = 6;
  qk_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(qk_array, sizeof(double *), ndim, blocks2);
  
  n_egrid = 0;
  n_tegrid = 0;
  n_usr = 0;
  egrid[0] = -1.0;
  SetCEEGridLimits(0.05, 8.0, 0);
  usr_egrid[0] = -1.0;
  tegrid[0] = -1.0;  
  SetCEQkMode(QK_DEFAULT, 1E-3);

  return 0;
}

int ReinitExcitation(int m) {
  
  if (m < 0) return 0;
  FreeExcitationQk();  
  n_egrid = 0;
  n_tegrid = 0;
  n_usr = 0;
  egrid[0] = -1.0;
  SetCEEGridLimits(0.05, 8.0, 0);
  usr_egrid[0] = -1.0;
  tegrid[0] = -1.0;  
  SetCEQkMode(QK_DEFAULT, 1E-3);

  return 0;
}
  




