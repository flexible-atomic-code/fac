#include "ionization.h"

static char *rcsid="$Id: ionization.c,v 1.29 2002/02/28 16:55:04 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

#define NINT 15
#define MAXNQK (MAXNE*(MAXNE+1)/2)
#define NPARAMS 4
#define BUFSIZE 128
#define NCBOMAX 6

static double cbo_params[(NCBOMAX+1)*NCBOMAX/2][NPARAMS+1] = {
  /* 1s */
  {1.130,	4.41,   -2.00,	3.80},
  
  /* 2s 2p */
  {0.823,	3.69,	0.62,	1.79},
  {0.530,	5.07,	1.20,	2.50},

  /* 3s 3p 3d */
  {0.652,	3.83,	0.64,	2.10},
  {0.551,	4.38,	1.83,	1.90},
  {0.280,	5.70,	2.21,	2.65},

  /* 4s 4p 4d 4f */
  {0.549,	3.93,	0.75,	2.51},
  {0.512,	4.29,	1.86,	1.68},
  {0.373,	4.87,	2.85,	1.90},
  {0.151,	5.94,	3.12,	2.38},

  /* 5s 5p 5d 5f 5g */
  {0.495,	4.02,	1.20,	2.31},
  {0.475,	4.48,	0.96,	2.80},
  {0.390,	4.80,	2.69,	1.75},
  {0.240,	5.46,	2.87,	2.62},
  {0.084,	6.11,	3.57,	2.37},

  /* 6s 6p 6d 6f 6g 6h */
  {0.475,	4.28,	1.05,	2.49},
  {0.450,	4.59,	0.75,	3.27},
  {0.390,	5.02,	1.66,	2.67},
  {0.285,	5.28,	3.02,	2.22},
  {0.145,	5.91,	2.84,	3.12},
  {0.047,	6.21,	3.90,	2.34},
};

static int egrid_type = -1;
static int usr_egrid_type = -1;
static int pw_type = -1;
static int qk_mode;
static double qk_fit_tolerance;

static double yegrid0[NINT];

static int n_usr = 0;
static double usr_egrid[MAXNUSR];
static double log_usr[MAXNUSR];
static double xusr[MAXNUSR];
static double log_xusr[MAXNUSR];

static int n_egrid = 0;
static double egrid[MAXNE];
static double log_egrid[MAXNE];
static double egrid_min = 0.05;
static double egrid_max = 8.0;
static int egrid_limits_type = 0;
static double sigma[MAXNE];
static double xegrid[MAXNTE][MAXNE];
static double log_xegrid[MAXNTE][MAXNE];

static int n_tegrid = 0;
static double tegrid[MAXNTE];
static double log_te[MAXNTE];

static double _dwork[17*MAXNQK];
static int _iwork[39*MAXNQK];

void uvip3p_(int *np, int *ndp, double *x, double *y, 
	     int *n, double *xi, double *yi);
void sdbi3p_(int *md, int *ndp, double *xd, double *yd, 
	     double *zd, int *nip, double *xi, double *yi,
	     double *zi, int *ierr, double *wk, int *iwk);
void rgbi3p_(int *md, int *nxp, int *nyd, double *xd, double *yd, 
	     double *zd, int *nip, double *xi, double *yi,
	     double *zi, int *ierr, double *wk);

static struct {
  int max_k;
  int qr;
  int max_kl;
  int max_kl_eject;
  int kl_cb;
  double tolerance;
  int nkl0;
  int nkl;
  int ns;
  double kl[MAXNKL+1];
  double log_kl[MAXNKL];
} pw_scratch = {8, 0, MAXKL, 8, 0, 1E-2, 0, 0};

static MULTI *qk_array;

int SetIEGrid(int n, double emin, double emax) {
  n_tegrid = SetTEGrid(tegrid, log_te, n, emin, emax);
  return n_tegrid;
}

int SetIEGridDetail(int n, double *x) {
  n_tegrid = SetTEGridDetail(tegrid, log_te, n, x);
  return n_tegrid;
}

int SetCIMaxK(int k) {
  pw_scratch.max_k = GetMaxRank();
  if (k >= 0) pw_scratch.max_k = Min(k, pw_scratch.max_k);
  return 0;
}

int SetCIPWOptions(int qr, int max, int max_eject, int kl_cb, double tol) {
  pw_scratch.qr = qr;
  if (max > MAXKL) {
    printf("The maximum partial wave reached in Ionization: %d\n", MAXKL);
    exit(1);
  }
  pw_scratch.max_kl = max;
  pw_scratch.max_kl_eject = max_eject;
  pw_scratch.kl_cb = kl_cb;
  pw_scratch.tolerance = tol;
  pw_scratch.nkl0 = 1;
  pw_scratch.kl[0] = 0;
  pw_scratch.log_kl[0] = -100.0;
  pw_scratch.nkl = 0;
  return 0;
}  

int SetCIEGridLimits(double min, double max, int type) {
  if (min <= 0) egrid_min = 0.05;
  else egrid_min = min;
  if (max <= 0) egrid_max = 8.0;
  else egrid_max = max;
  egrid_limits_type = type;

  return 0;
}

int SetCIEGridDetail(int n, double *xg) {
  n_egrid = SetEGridDetail(egrid, log_egrid, n, xg);
  return n_egrid;
}

int SetUsrCIEGridType(int type) {
  if (type >= 0) usr_egrid_type = type;
  return 0;
}

int SetCIEGrid(int n, double emin, double emax, double eth) {
  n_egrid = SetEGrid(egrid, log_egrid, n, emin, emax, eth);
  return n_egrid;
}

int SetUsrCIEGridDetail(int n, double *xg) {
  if (n > MAXNUSR) {
    printf("Max # of grid points reached \n");
    return -1;
  }
  n_usr = SetEGridDetail(usr_egrid, log_usr, n, xg);
  return n_usr;
}

int SetUsrCIEGrid(int n, double emin, double emax, double eth) {
  if (n > MAXNUSR) {
    printf("Max # of grid points reached \n");
    return -1;
  }
  n_usr = SetEGrid(usr_egrid, log_usr, n, emin, emax, eth);
  return n_usr;
}

int SetCIPWGrid(int ns, int *n, int *step) {
  if (pw_scratch.nkl0 <= 0) SetCIPWOptions(0, 256, 8, 64, 5E-2);
  pw_scratch.nkl = SetPWGrid(&(pw_scratch.nkl0),
			     pw_scratch.kl,
			     pw_scratch.log_kl,
			     pw_scratch.max_kl,
			     &ns, n, step);
  pw_scratch.ns = ns;
  return 0;
}

int SetCIQkMode(int m, double tol) {
  if (m == QK_DEFAULT) qk_mode = QK_BED;
  else qk_mode = m;
  if (tol > 0.0) qk_fit_tolerance = tol;
  return 0;
}

int CIRadialQk(double *qk, int ie1, int ie2, int kb, int kbp, int k) {
  double pk[MAXNTE][MAXNKL];
  double y2[MAXNKL];
  double e1, e2, e0, te;
  ORBITAL *orb;
  int kappab, jb, klb, i, j, t, kl, klp;
  int js[4], ks[4];
  int jmin, jmax, ko2;
  int kf, kappaf, kf0, kappa0, kf1, kappa1;
  int kl0, kl0p, kl1, kl1p;
  int j0, j1, j1min, j1max;
  double z, z2, r, rp, sd, se, s;
  double qkjt, qkj, qklt, qkl, qkl0;
  double eps, a, b, h, jb1;
  int kl_max0, kl_max1, kl_max2, max0;
  int type, last_kl0, second_last_kl0;

  ko2 = k/2;

  for (i = 0; i < n_tegrid; i++) {
    for (j = 0; j < pw_scratch.nkl0; j++) {
      pk[i][j] = 0.0;
    }
  }
  for (i = 0; i < n_tegrid; i++) qk[i] = 0.0;
  
  e1 = egrid[ie1];
  r = GetRMax();
  z = GetResidualZ();
  z2 = z*z;
  t = r*sqrt(e1+2.0*z/r);
  kl_max0 = pw_scratch.max_kl;
  kl_max0 = Min(kl_max0, t);
  max0 = pw_scratch.ns;

  e2 = egrid[ie2];
  t = r*sqrt(e2+2.0*z/r);
  kl_max2 = pw_scratch.max_kl_eject;
  kl_max2 = Min(kl_max2, t);

  orb = GetOrbital(kb);
  kappab = orb->kappa;
  GetJLFromKappa(kappab, &jb, &klb);
  klb /= 2;
  jmin = abs(k - jb);
  jmax = k + jb;
  jb1 = jb + 1.0;

  qkjt = 0.0;
  for (j = jmin; j <= jmax; j += 2) {
    if (qkjt && fabs(qkj/qkjt) < eps) break;
    qkj = 0.0;
    for (klp = j - 1; klp <= j + 1; klp += 2) {
      kappaf = GetKappaFromJL(j, klp);
      kl = klp/2;
      if (kl > kl_max2) break;
      if (IsEven(kl + klb + ko2)) {
	type = ko2;
      } else {
        type = -1;
      }
      if (kl < pw_scratch.qr) {
	js[2] = 0;
      } else {
	js[2] = j;
	if (IsOdd(kl)) {
	  if (kappaf < 0) kappaf = -kappaf - 1;
	} else {
	  if (kappaf > 0) kappaf = -kappaf - 1;
	}
      }
      kf = OrbitalIndex(0, kappaf, e2);	
      ks[2] = kf;  

      if (type >= 0 && type < CBMULTIPOLES) {
	kl_max1 = pw_scratch.kl_cb;
      } else {
	kl_max1 = kl_max0;
      }
      eps = pw_scratch.tolerance;
      if (type >= CBMULTIPOLES) {
	z = GetCoulombBetheAsymptotic(tegrid[0]+e2, e1);
      }

      last_kl0 = 0;
      second_last_kl0 = 0; 
      qklt = 0.0;
      for (t = 0; ; t++) {
        if (second_last_kl0) {
          last_kl0 = 1;
        } else {
	  kl0 = pw_scratch.kl[t];
          if (kl0 > max0) {
	    if (type < 0) {
	      rp *= qkl;
	      rp = rp/qklt;
	      if (rp < eps) last_kl0 = 1;
	    } else if (type >= CBMULTIPOLES) {
	      h = z*qkl;
	      h = h/(h+qklt);
	      if (h < eps) last_kl0 = 1;
	      else {
		rp = fabs(1.0 - rp/z);
		rp *= h;
		if (rp < eps) last_kl0 = 1;
	      }
	    } else {
	      z = (GetCoulombBethe(ie2, 0, ie1, type, 1))[t-1];
	      h = z*qkl;
	      h = h/(h+qklt);
	      if (h < eps) last_kl0 = 1;
	      else {
		z = (GetCoulombBethe(ie2, 0, ie1, type, 0))[t-1];
		z = z/(1.0-z);
		rp = fabs(1.0 - rp/z);
		rp *= h;
		if (rp < eps) last_kl0 = 1;
	      }
	    }
	    if (pw_scratch.kl[t+1] > kl_max1) {
	      second_last_kl0 = 1;
	    }  
	  } 
	}	  
	kl0p = 2*kl0;
	qkl0 = qkl;
        qkl = 0.0;
	for (j0 = abs(kl0p - 1); j0 <= kl0p + 1; j0 += 2) {
	  kappa0 = GetKappaFromJL(j0, kl0p);
	  if (kl0 < pw_scratch.qr) {
	    js[1] = 0;
	  } else {
	    js[1] = j0;
	    if (IsOdd(kl0)) {
	      if (kappa0 < 0) kappa0 = -kappa0 - 1;
	    } else {
	      if (kappa0 > 0) kappa0 = -kappa0 - 1;
	    }
	  }
	  j1min = abs(j0 - k);
	  j1max = j0 + k;
	  for (j1 = j1min; j1 <= j1max; j1 += 2) {
	    for (kl1p = j1 - 1; kl1p <= j1 + 1; kl1p += 2) {
	      kappa1 = GetKappaFromJL(j1, kl1p);
	      kl1 = kl1p/2;
	      if (kl1 < pw_scratch.qr) {
		js[3] = 0;
	      } else {
		js[3] = j1;
		if (IsOdd(kl1)) {
		  if (kappa1 < 0) kappa1 = -kappa1 - 1;
		} else {
		  if (kappa1 > 0) kappa1 = -kappa1 - 1;
		}
	      }
	      kf1 = OrbitalIndex(0, kappa1, e1);
	      ks[3] = kf1;
	      for (i = 0; i < n_tegrid; i++) {
		te = tegrid[i];		
		e0 = e1 + e2 + te;
		kf0 = OrbitalIndex(0, kappa0, e0);
		ks[1] = kf0;

		js[0] = 0;
		ks[0] = kb;
		if (kl1 >= pw_scratch.qr &&
		    kl0 >= pw_scratch.qr) {
		  SlaterTotal(&sd, &se, js, ks, k, -1);
		} else {
		  SlaterTotal(&sd, &se, js, ks, k, 1);
		} 
		r = sd + se;
                if (!r) break;

		if (kbp != kb) {
		  js[0] = 0;
		  ks[0] = kbp;
		  if (kl1 >= pw_scratch.qr &&
		      kl0 >= pw_scratch.qr) {
		    SlaterTotal(&sd, &se, js, ks, k, -1);
		  } else {
		    SlaterTotal(&sd, &se, js, ks, k, 1);
		  } 
		  rp = sd + se;
                  if (!rp) break;
		} else {
		  rp = r;
		}
		r = r*rp;
		if (i == 0) qkl += r;
		pk[i][t] += r;
	      }
	    }
	  }
	}
        if (t == 0) {
	  qklt += qkl;
	} else if (!last_kl0 && !second_last_kl0) {
	  if (qkl + 1.0 == 1.0) rp = 0.0;
	  else rp = qkl / qkl0;
	  h = kl0 - pw_scratch.kl[t-1];
	  if (h > 1) {
	    a = (1.0 - rp) * qkl0;
	    rp = pow(rp, 1.0/h);
	    rp = rp/(1.0 - rp);
	    a *= rp;
	    qklt += a;
	  } else {
	    qklt += qkl;
	    rp = rp/(1.0 - rp);
	  }
	} else if (last_kl0) {
	  if (type >= CBMULTIPOLES) {
	    for (i = 0; i < n_tegrid; i++) {
	      b = GetCoulombBetheAsymptotic(tegrid[i]+e2, e1);
	      r = qkl * b;
	      qk[i] += r;  
	    }
	  } else if (type >= 0) {
	    for (i = 0; i < n_tegrid; i++) {
	      b = (GetCoulombBethe(ie2, i, ie1, type, 1))[t];
	      r = qkl*b;
	      qk[i] += r;
	    }    
          } 
	  break;
	}
      }
      qkj += qklt;
    }  
    qkjt += qkj;
  }
  
  t = pw_scratch.nkl;
  for (i = 0; i < n_tegrid; i++) {
    spline(pw_scratch.log_kl, pk[i], t, 1E30, 1E30, y2);
    r = pk[i][0];
    for (j = 1; j < t; j++) {
      r += pk[i][j];
      kl0 = pw_scratch.kl[j-1];
      kl1 = pw_scratch.kl[j];
      for (kl = kl0+1; kl < kl1; kl++) {
	splint(pw_scratch.log_kl, pk[i], y2, t, LnInteger(kl), &s);
	r += s;
      }
    } 
    qk[i] += r;
    qk[i] /= jb1;
  }

  return 0;
}

int CIRadialQkIntegrated(double *qke, double te, int kb, int kbp) {
  int np, i, j, k, nd;
  double *qk, qkc[MAXNTE];

  qk = CIRadialQkIntegratedTable(kb, kbp);
  if (qk == NULL) {
    return -1;
  }

  if (n_tegrid > 1) {
    nd = 1;
    np = 3;
    for (i = 0; i < NPARAMS; i++) {
      k = i;
      for (j = 0; j < n_tegrid; j++) {
	qkc[j] = qk[k];
	k += NPARAMS;
      }
      uvip3p_(&np, &n_tegrid, tegrid, qkc, &nd, &te, &qke[i]);
    }
  } else {
    for (i = 0; i < NPARAMS; i++) {
      qke[i] = qk[i];
    }
  }
  if (qk_mode == QK_BED) {
    k = n_tegrid*NPARAMS;
    qke[NPARAMS] = qk[k];
    qke[NPARAMS+1] = qk[k+1];
    return qke[NPARAMS];
  } else {
    return 0;
  }
}

void CIRadialQkBasis(int npar, double *yb, double x, double log_x) {
  double y1, y2;

  y1 = 1.0/x;
  y2 = 1.0 - y1;
  yb[0] = y2*y2;
  yb[1] = y2*y1 ;
  yb[2] = yb[1]*y1;
}

void CIRadialQkFromFit(int np, double *p, int n, 
		       double *x, double *logx, double *y) {
  double a, c, d, e, f, g;
  int i;

  if (qk_mode == QK_BED) {
    for (i = 0; i < n; i++) {
      a = 3.5 + p[4];
      c = 2.0/(1.0+x[i]);
      d = p[2]+sqrt(c);
      e = p[2]+1.0;
      f = pow(c,a);
      g = pow(e/d, p[1]);
      y[i] = p[0]*(1.0 - f*g)*logx[i];
      a = 1.0 - 1.0/x[i] - logx[i]/(1.0 + x[i]);
      y[i] += p[3]*a;
    }
  } else {
    for (i = 0; i < n; i++) {
      a = 1.0/x[i];
      d = 1.0 - a;
      c = d*a;
      y[i] = p[0]*logx[i];
      y[i] += p[1]*d*d;
      y[i] += p[2]*c;
      y[i] += p[3]*c*a;
    }
  }
}

void BEDFromFit(int np, double *p, int n, double *x, double *logx,
		double *y, double *dy, int ndy, void *extra) {
  int kl, i, k;
  double a, c, d, e, f, g, h;

  kl = *((int *) extra);
  if (ndy <= 0) {
    for (i = 0; i < n; i++) {
      a = 3.5 + kl;
      c = sqrt(x[i]);
      d = p[1] + c;
      e = p[1] + 1;
      f = logx[i]*a;
      h = log(e/d);
      g = h*p[0];
      y[i] = (1.0 - exp(f+g));
    }
  } else {
    for (i = 0; i < n; i++) {
      k = i;
      a = 3.5 + kl;
      c = sqrt(x[i]);
      d = p[1] + c;
      e = p[1] + 1;
      f = logx[i]*a;
      h = log(e/d);
      g = h*p[0];
      a = -exp(f+g);
      dy[k] = a*h;
      k += ndy;
      dy[k] = a*p[0]*(1.0/e - 1.0/d);
    }
  }
}

int LoadCIRadialQkIntegrated(int n, char *s) {
  int kappa, i, nd, np, kl, kl2, j, k, ite;
  FILE *f;
  ORBITAL *orb;
  double e, e0, a, **p;
  double xte[MAXNTE], xte0[MAXNTE];
  double qkt[NPARAMS][MAXNTE];
  int index[3];
  double emin, emax;
  char buf[BUFSIZE];
  char *sb;
  int nz, nza, nqkz, nqkc, jp, jz, md, ierr;
  int iz, iza, nqk;
  double *z, *za, *qkc, *dtmp;
  double z0, za0;

  if (n <= 0) {
    return 0;
  }

  emin = 1E10;
  emax = 0.0;
  if (s == NULL) {
    if (n > 6) {
      return 0;
    }
    iz = (n-1)*n/2;    
  }
  for (kl = 0; kl < n; kl++) {
    kl2 = 2*kl;
    for (j = kl2-1; j <= kl2+1; j += 2) {
      if (j < 0) continue;
      kappa = GetKappaFromJL(j, kl2);
      i = OrbitalExists(n, kappa, 0.0);
      if (i < 0) continue;
      orb = GetOrbital(i);
      e0 = -(orb->energy);
      if (s == NULL) {
	index[0] = 0;
	index[1] = i;
	index[2] = i;
	p = MultiSet(qk_array, index, NULL);
	if (*p) free(*p);
	*p = (double *) malloc(sizeof(double)*NPARAMS);
	for (np = 0; np < NPARAMS; np++) {
	  (*p)[np] = (cbo_params[iz+kl][np]/e0) * (PI/32.0);
	}
      }
      if (e0 > emax) emax = e0;
      if (e0 < emin) emin = e0;
    }
  }

  if (s == NULL) {
    n_tegrid = 1;
    emin = 0.5*(emin+emax);
    SetIEGrid(n_tegrid, emin, emax);
    qk_mode = QK_CB;
    return 0;
  }

  f = fopen(s, "r");
  
  sb = fgets(buf, BUFSIZE, f);
  if (sb == NULL) return -1;
  sscanf(sb, "%d", &qk_mode);
  sb = fgets(buf, BUFSIZE, f);
  if (sb == NULL) return -1;
  sscanf(sb, "%d", &n_tegrid);
  if (n_tegrid == 1) emin = 0.5*(emin+emax);
  SetIEGrid(n_tegrid, emin, emax);

  sb = fgets(buf, BUFSIZE, f);
  if (sb == NULL) return -1;
  sscanf(sb, "%d", &nz);
  z = (double *) malloc(sizeof(double)*nz);
  sb += 6;
  for (i = 0; i < nz; i++) {
    sscanf(sb, "%lf", &(z[i]));
    sb += 6;
  }
  sb = fgets(buf, BUFSIZE, f);
  if (sb == NULL) {
    free(z);
    return -1;
  }
  sscanf(sb, "%d", &nza);
  za = (double *) malloc(sizeof(double)*nza);
  sb += 6;
  for (i = 0; i < nza; i++) {
    sscanf(sb, "%le", &(za[i]));
    sb += 12;
  }

  nqk = n_tegrid*NPARAMS;
  if (qk_mode == QK_BED) nqk += 1;

  nqkz = nz*nza;
  nqkc = n_tegrid*nqkz;
  qkc = (double *) malloc(sizeof(double)*nqkc*NPARAMS);
  if (nz > 1 && nza > 1) {
    dtmp = (double *) malloc(sizeof(double)*3*nqkz);
  } else {
    dtmp = NULL;
  }
  z0 = GetAtomicNumber();
  za0 = GetResidualZ();
  za0 = za0/z0;
  
  while ((sb = fgets(buf, BUFSIZE, f))) {
    sscanf(buf, "%d%d%d", &np, &kl2, &j);
    if (np != n) {
      for (ite = 0; ite < nqkc; ite++) 
	fgets(buf, BUFSIZE, f);
    } else {
      kappa = GetKappaFromJL(j, kl2);
      k = OrbitalExists(n, kappa, 0.0);
      if (k < 0) {
	for (i = 0; i < nqkc; i++) 
	  fgets(buf, BUFSIZE, f);
	continue;
      }
      index[0] = 0;
      index[1] = k;
      index[2] = k;
      p = MultiSet(qk_array, index, NULL);
      if (*p) free(*p);
      *p = (double *) malloc(sizeof(double)*nqk);
      
      orb = GetOrbital(k);
      e = -(orb->energy);
      for (ite = 0; ite < n_tegrid; ite++) {
	xte[ite] = tegrid[ite]/e;
      }
      jz = 0;
      for (iz = 0; iz < nz; iz++) {
	for (iza = 0; iza < nza; iza++) {
          for (ite = 0; ite < n_tegrid; ite++) {
	    sb = fgets(buf, BUFSIZE, f);
	    sscanf(sb, "%lE%lE", &e0, &(xte0[ite]));
	    sb += 24;
	    for (i = 0; i < NPARAMS; i++) {
	      sscanf(sb, "%lE", &(qkt[i][ite]));
	      sb += 12;
	    }
          }
	  jp = jz;
	  if (n_tegrid == 1) {
	    for (i = 0; i < NPARAMS; i++) {
	      if (qk_mode == QK_BED && (i == 1 || i == 2)) {
		qkc[jp] = qkt[i][0];
	      } else {
		qkc[jp] = qkt[i][0]/e;
	      }
	      jp += nqkz;
	    }
	  } else {
	    np = 3; 
	    nd = 1;
	    for (ite = 0; ite < n_tegrid; ite++) {
	      for (i = 0; i < NPARAMS; i++) {
		uvip3p_(&np, &n_tegrid, xte0, qkt[i], &nd, &(xte[ite]), &a);
		if (qk_mode == QK_BED && (i == 1 || i == 2)) {
		  qkc[jp] = a;
		} else {
		  qkc[jp] = a/e;
		}
		jp += nqkz;
	      }
	    }
	  }
	  jz++;
	}
      }
      np = 3;
      nd = 1;
      jp = 0;
      jz = 0;
      if (nqkz == 1) {
	for (ite = 0; ite < n_tegrid; ite++) {
	  for (i = 0; i < NPARAMS; i++) {
	    (*p)[jp] = qkc[jp];
	    jp++;
	  }
	}
	if (qk_mode == QK_BED) {
	  (*p)[jp] = kl2/2;
	}
      } else if (nz == 1) {
	for (ite = 0; ite < n_tegrid; ite++) {
	  for (i = 0; i < NPARAMS; i++) {
	    uvip3p_(&np, &nza, za, &(qkc[jz]), &nd, &za0, &a);
	    (*p)[jp] = a;
	    jp++;
	    jz += nqkz;
	  }
	}
	if (qk_mode == QK_BED) {
	  (*p)[jp] = kl2/2;
	}
      } else if (nza == 1) {
	for (ite = 0; ite < n_tegrid; ite++) {
	  for (i = 0; i < NPARAMS; i++) { 
	    uvip3p_(&np, &nz, z, &(qkc[jz]), &nd, &z0, &a);
	    (*p)[jp] = a;
	    jp++;
	    jz += nqkz;
	  }
	}
	if (qk_mode == QK_BED) {
	  (*p)[jp] = kl2/2;
	}
      } else {
	md = 1;
	for (ite = 0; ite < n_tegrid; ite++) {
	  for (i = 0; i < NPARAMS; i++) { 
	    rgbi3p_(&md, &nza, &nz, za, z, &(qkc[jz]), &nd, 
		    &za0, &za0, &a, &ierr, dtmp);
	    (*p)[jp] = a;
	    jp++;
	    jz += nqkz;
	  }
	}
	if (qk_mode == QK_BED) {
	  (*p)[jp] = kl2/2;
	}
      }
    }
  }

  if (dtmp) free(dtmp);
  free(qkc);
  free(z);
  free(za);

  fclose(f);
  return 0;
}

int PrepCIRadialQkIntegrated(int nz, double *z, int na, double *a,
			     int np, int *n, int nte, 
			     double emin, double emax, char *s) {
#define MAXNSHELLS 128
  FILE *f;
  int i, ite, ie, k, kl, kl2, j;
  int jz, jp, iz, ia, in;
  int an[MAXNSHELLS], akappa[MAXNSHELLS];
  double mnq, tnq;
  double anq[MAXNSHELLS];
  double e0, *qkc, delta, rz;
  double *qk[MAXNTE][NPARAMS], *e, *xte[MAXNTE]; 
  ORBITAL *orb;
  int nshells, maxshells, kappa;

  k = 0;
  for (i = 1; i <= 10; i++) {
    for (kl = 0; kl < i; kl++) {
      kl2 = 2*kl;
      for (j = kl2-1; j <= kl2+1; j += 2) {
	if (j < 0) continue;
	an[k] = i;
	akappa[k] = GetKappaFromJL(j, kl2);
	k++;
      }
    }
  }
  maxshells = k;

  i = 0;
  for (in = 0; in < np; in++) {
    i += 2*n[in]-1;
  }
  i *= nz*na;
  e = (double *) malloc(sizeof(double)*i);
  for (ite = 0; ite < nte; ite++) {
    xte[ite] = (double *) malloc(sizeof(double)*i);
    for (j = 0; j < NPARAMS; j++) {
      qk[ite][j] = (double *) malloc(sizeof(double)*i);
    }
  }
   
  egrid_type = 1;
  pw_type = 0; 
  yegrid0[0] = 0.0;
  delta = 0.5/(NINT-1.0);
  for (i = 1; i < NINT; i++) {
    yegrid0[i] = yegrid0[i-1] + delta;
  }	    
  if (pw_scratch.nkl == 0) {
    SetCIPWGrid(0, NULL, NULL);
  }	    

  jz = 0;
  for (iz = 0; iz < nz; iz++) {    
    SetAtom(" ", z[iz], -1.0);
    for (ia = 0; ia < na; ia++) {
      tnq = (1.0-a[ia])*z[iz] + 1.0;
      for (i = 0; i < maxshells; i++) {
	j = GetJFromKappa(akappa[i]);
	mnq = j+1.0;
	if (tnq > mnq) {
	  anq[i] = mnq;
	  tnq -= mnq;
	} else if (tnq > 0) {
	  anq[i] = tnq;
	  tnq = -1.0;
	} else {
	  break;
	}
      }
      nshells = i;
         
      FreeIonizationQk();
      FreeSlaterArray();
      FreeResidualArray();  
     
      ClearOrbitalTable(1);
      SetAverageConfig(nshells, an, akappa, anq);
      SetRadialGrid(-1.0, -1.0);
      OptimizeRadial(0, NULL, NULL);
      rz = GetResidualZ();
      jp = jz;
      for (in = 0; in < np; in++) {
	for (kl = 0; kl < n[in]; kl++) {
	  kl2 = 2*kl;
	  for (j = kl2-1; j <= kl2+1; j += 2) {
	    if (j < 0) continue;
	    kappa = GetKappaFromJL(j, kl2);
	    k = OrbitalIndex(n[in], kappa, 0.0);
	    printf("%10.3E %10.3E %d %d %d\n", z[iz], a[ia], n[in], kl2, j);
	    orb = GetOrbital(k);
	    e0 = -(orb->energy);	
	    if (qk_mode != QK_CB) {
	      SetTransitionOptions(2, 1, 4, 4);
	      SetRRTEGrid(1, e0, e0);
	      SetPEGridLimits(-1, -1, 0);
	      SetRecQkMode(QK_FIT, 0.1);
	      SetPEGrid(0, -1.0, -1.0, 1.0);
	      PrepRREGrids(e0);
	    }
	    SetIEGrid(nte, emin*e0, emax*e0);
	    if (kl == 0) SetCIEGrid(6, 0.05*e0, 8.0*e0, e0);
	    for (ie = 0; ie < n_egrid; ie++) {
	      for (i = 0; i < n_tegrid; i++) {
		xegrid[i][ie] = egrid[ie]/tegrid[i];
		if (egrid_type == 1) xegrid[i][ie] += 1.0;
		log_xegrid[i][ie] = log(xegrid[i][ie]);
	      }
	      sigma[ie] = 1.0/sqrt(xegrid[0][ie]);
	    }
	    PrepCoulombBethe(n_egrid, n_tegrid, n_egrid, 
			     rz, egrid, tegrid, egrid,
			     pw_scratch.nkl, pw_scratch.kl, 
			     egrid_type, pw_type, 0);
	    qkc = CIRadialQkIntegratedTable(k, k);
	    e[jp] = e0;
	    k = 0;
	    for (ite = 0; ite < n_tegrid; ite++) {
	      xte[ite][jp] = tegrid[ite]/e0;
	      for (i = 0; i < NPARAMS; i++) {
		if (qk_mode == QK_BED && (i == 1 || i == 2)) {
		  qk[ite][i][jp] = qkc[k];		
		} else {
		  qk[ite][i][jp] = qkc[k]*e0;
		}
		k++;
	      }
	    }
	    jp += nz*na;
	  }
	}
      }
      jz++;
    }
  }
  
  f = fopen(s, "w");  
  fprintf(f, "%-5d\n", qk_mode);
  fprintf(f, "%-5d\n", n_tegrid);
  fprintf(f, "%-5d ", nz);
  for (iz = 0; iz < nz; iz++) {
    fprintf(f, "%-5d ", (int)(z[iz]));
  }
  fprintf(f, "\n");
  fprintf(f, "%-5d ", na);
  for (ia = 0; ia < na; ia++) {
    fprintf(f, "%11.4E ", a[ia]);
  }
  fprintf(f, "\n"); 

  jp = 0;
  for (in = 0; in < np; in++) {
    for (kl = 0; kl < n[in]; kl++) {
      kl2 = 2*kl;
      for (j = kl2-1; j <= kl2+1; j += 2) {
	if (j < 0) continue;
	fprintf(f, "%-2d %-2d %-2d\n", n[in], kl2, j);
	for (iz = 0; iz < nz; iz++) {
	  for (ia = 0; ia < na; ia++) {
	    for (ite = 0; ite < n_tegrid; ite++) {
	      fprintf(f, "%11.4E %11.4E ", e[jp], xte[ite][jp]);
	      for (i = 0; i < NPARAMS; i++) {
		fprintf(f, "%11.4E ", qk[ite][i][jp]);
	      }
	      fprintf(f, "\n");
	    }
	    jp++;
	  }
	}
      }
    }
  }

  free(e);
  for (ite = 0; ite < n_tegrid; ite++) {
    free(xte[ite]);
    for (i = 0; i < NPARAMS; i++) {
      free(qk[ite][i]);
    }
  }

  fclose(f);
  return 0;

#undef MAXNSHELLS
}

int SaveCIRadialQkIntegrated(int n, char *s) {
  int kappa, i, kl, kl2, j, k, ite;
  ORBITAL *orb;
  FILE *f;
  double e0, z, rz, a, *qkc;

  z = GetAtomicNumber();
  rz = GetResidualZ();
  rz = rz/z;
  
  f = fopen(s, "w");  
  
  fprintf(f, "%-5d\n", qk_mode);
  fprintf(f, "%-5d\n", n_tegrid);
  fprintf(f, "%-5d %-5d\n", 1, (int)z);
  fprintf(f, "%-5d %11.4E\n", 1, rz);
  for (kl = 0; kl < n; kl++) {
    kl2 = 2*kl;
    for (j = kl2-1; j <= kl2+1; j += 2) {
      if (j < 0) continue;
      kappa = GetKappaFromJL(j, kl2);
      i = OrbitalExists(n, kappa, 0.0);
      if (i < 0) continue;
      orb = GetOrbital(i);
      e0 = -(orb->energy);
      qkc = CIRadialQkIntegratedTable(i, i);
      if (qkc == NULL) continue;
      fprintf(f, "%-2d %-2d %-2d\n", n, kl2, j); 
      k = 0;
      for (ite = 0; ite < n_tegrid; ite++) {
	a = tegrid[ite]/e0;
	fprintf(f, "%11.4E %11.4E ", e0, a);
	for (i = 0; i < NPARAMS; i++) {
	  if (qk_mode == QK_BED && (i == 1 || i == 2)) {
	    a = qkc[k];
	  } else {
	    a = qkc[k]*e0;
	  }
	  fprintf(f, "%11.4E ", a);
	  k++;
	}
	fprintf(f, "\n");
      }
    }
  }
  
  fclose(f);
  return 0;
}     
 
int CIRadialQkBED(double *dp, double *b, int kl, double *p) {
  double integrand[NINT];
  double x[NINT], y[NINT], t[NINT], logt[NINT], d;
  int i, j, n;
  int ipvt[2];
  int ierr;
  int lwa=10+NINT;
  double wa[10+NINT];
  double fvec[NINT], fjac[2*NINT];
  
  for (i = 1; i < NINT; i++) {
    x[i] = 1.0/(2.0*yegrid0[i]);
    t[i] = log(x[i]);
  }
  RRRadialQkFromFit(3, p, NINT-1, &(x[1]), &(t[1]),
		    &(integrand[1]), NULL, 0, &kl);
  integrand[0] = 0.0;
  for (i = 1; i < NINT; i++) {
    integrand[i] *= x[i]*x[i];
  }

  d = 4.0*yegrid0[1];
  if (qk_mode == QK_DW) {
    dp[0] = Simpson(integrand, 0, NINT-1);
    dp[0] *= d;
  } else {
    n = (NINT+1)/2;
    j = n-1;
    t[j] = 1.0;
    y[j--] = 0.0;
    for (i = NINT-1; i > 0; i -= 2, j--) {
      y[j] = Simpson(integrand, i-2, i);
      y[j] *= d;
      y[j] += y[j+1];
      t[j] = 2.0*yegrid0[i-2];
    }
    dp[0] = y[0];
    dp[1] = p[1];
    dp[2] = 1.0/p[2];
    for (j = 1; j < n-1; j++) {
      y[j] = y[j]/dp[0];
      logt[j] = log(t[j]);
    }
    ierr = NLSQFit(2, &(dp[1]), qk_fit_tolerance, ipvt, fvec, fjac, NINT,
		   wa, lwa, n-2, &(t[1]), &(logt[1]), &(y[1]), &(y[1]),
		   BEDFromFit, &kl);
    for (i = 1; i < NINT; i++) {
      integrand[i] *= x[i];
    }
    *b = Simpson(integrand, 0, NINT-1);
    *b *= d;
  }
  return 0;
}

double *CIRadialQkIntegratedTable(int kb, int kbp) {
  ORBITAL *orb;
  int index[3], ie1, ie2, ite, i, j, kl;
  double **p, *qkc, *qk, qkt[MAXNE], tol;
  double xr[MAXNQK], yr[MAXNQK], zr[MAXNE];
  double x, y, integrand[NINT];
  int np, nd, ierr, nqk;
  double dp[NPARAMS], a, b, c, params[10], te;

  if (qk_mode == QK_CB && kbp != kb) return NULL;

  index[0] = 0;
  index[1] = kb;
  index[2] = kbp;
  
  p = (double **) MultiSet(qk_array, index, NULL);
  if (*p) {
    return (*p);
  } 
  nqk = n_tegrid*NPARAMS;
  if (qk_mode == QK_BED) nqk += 2;
  *p = (double *) malloc(sizeof(double)*nqk);
  qkc = *p;

  if (qk_mode == QK_CB) {
    orb = GetOrbital(kb);
    i = orb->n;
    j = orb->kappa;
    kl = GetLFromKappa(j);
    kl = kl/2;
    i = (i-1)*i/2;
    for (j = 0; j < NPARAMS; j++) {
      (*p)[j] = (cbo_params[i+kl][j]/(-orb->energy)) * (PI/32.0);
    }
    return *p;
  }

  te = 0.5*(tegrid[0]+tegrid[n_tegrid-1]);
  kl = RRRadialQk(params, te, kb, kbp, -1);
  CIRadialQkBED(dp, &b, kl, params);
  a = 1.0/(16.0*te);
  if (qk_mode == QK_BED) {
    for (ite = 0; ite < n_tegrid; ite++) {
      c = a*tegrid[ite];
      qkc[0] = dp[0]*c;
      qkc[1] = dp[1];
      qkc[2] = dp[2];
      qkc[3] = 4.0*PI/(16.0*tegrid[ite]) - b*c;
      qkc += NPARAMS;
    }
    *qkc = kl;
    *(qkc+1) = params[4];
    return (*p);
  }
      
  qk = CIRadialQkTable(kb, kbp);
  if (qk == NULL) return NULL;

  nd = 1;
  for (ite = 0; ite < n_tegrid; ite++) {
    c = a*tegrid[ite];
    for (ie1 = 0; ie1 < n_egrid; ie1++) {
      zr[ie1] = log(1.0 + egrid[ie1]/tegrid[ite]);
    }
    j = 0;
    for (ie1 = 0; ie1 < n_egrid; ie1++) {
      for (ie2 = 0; ie2 <= ie1; ie2++) {
	xr[j] = zr[ie1];
	yr[j] = zr[ie2];
	j++;
      }
    }
    nqk = j;
    np = 0;
    qkc[0] = dp[0]*c;
    tol = qk_fit_tolerance;
    for (ie1 = 0; ie1 < n_egrid; ie1++) {
      for (i = 0; i < NINT; i++) {	
	x = log(1.0 + egrid[ie1]*(1.0-yegrid0[i])/tegrid[ite]);
	y = log(1.0 + egrid[ie1]*yegrid0[i]/tegrid[ite]);
	if (np == 0) np = 1;
	else np = 3;
	sdbi3p_(&np, &nqk, xr, yr, qk, &nd, &x, &y, &integrand[i],
		&ierr, _dwork, _iwork);
	if (kb == kbp) {
	  integrand[i] = exp(integrand[i]);
	}
      }   
      x = Simpson(integrand, 0, NINT-1);
      x *= yegrid0[1]*egrid[ie1];
      qkt[ie1] = x - qkc[0]*log_xegrid[ite][ie1];
    }
    SVDFit(NPARAMS-1, qkc+1, NULL, tol, n_egrid, xegrid[ite], 
	   log_xegrid[ite], qkt, sigma, CIRadialQkBasis);
    qk += nqk;
    qkc += NPARAMS;
  }

  /*
  qkc = *p;
  for (j = 0; j < NPARAMS; j++) {
    i = j;
    for (ite = 0; ite < n_tegrid; ite++) {
      printf("%d %d %10.3E %10.3E\n", 
	     kb, kbp, tegrid[ite], qkc[i]);
      i += NPARAMS;
    }
    printf("\n\n");
  }
  printf("===============\n");
  */

  return (*p);
}

double *CIRadialQkTable(int kb, int kbp) {
  int index[3], j, k, t;
  int nqk, ie1, ie2, ite;
  double **p, *qk, qi[MAXNTE];

  index[0] = 1;
  index[1] = kb;
  index[2] = kbp;

  p = (double **) MultiSet(qk_array, index, NULL);
  if (*p) {
    return (*p);
  }

  nqk = n_egrid*(n_egrid+1)/2;
  *p = (double *) malloc(sizeof(double)*n_tegrid*nqk);
  qk = *p;

  for (j = 0; j < nqk*n_tegrid; j++) qk[j] = 0.0;
  j = 0;
  for (ie1 = 0; ie1 < n_egrid; ie1++) {
    for (ie2 = 0; ie2 <= ie1; ie2++) {
      for (k = 0; k <= pw_scratch.max_k; k += 2) {	
	CIRadialQk(qi, ie1, ie2, kb, kbp, k);
	for (ite = 0; ite < n_tegrid; ite++) {
	  t = nqk*ite + j;
	  qi[ite] /= k + 1.0;
	  qk[t] += qi[ite];
	}
	if (k > 2 && fabs(qi[0]/qk[j]) < pw_scratch.tolerance) break;
      }
      j++;
    }
  }

  if (kb == kbp) {
    t = 0;
    for (ite = 0; ite < n_tegrid; ite++) {
      for (ie1 = 0; ie1 < n_egrid; ie1++) {
	for (ie2 = 0; ie2 <= ie1; ie2++) {
	  qk[t] = log(qk[t]);
	  t++;
	}
      }
    }
  }

  return qk;
}
  
int IonizeStrength(double *qku, int *nqkc, double **qkc, double *te, 
		   int b, int f) {
  int i, ip, j, ierr;
  LEVEL *lev1, *lev2;
  ANGULAR_ZFB *ang;
  double c, qke[MAXNUSR];
  int nz, j0, j0p, kb, kbp, nq, nqk;
  double *p;

  lev1 = GetLevel(b);
  lev2 = GetLevel(f);
  *te = lev2->energy - lev1->energy;
  if (*te <= 0) return -1;
  
  nz = AngularZFreeBound(&ang, f, b);
  if (nz <= 0) return -1;

  nq = 0;
  nqk = NPARAMS;
  if (qk_mode == QK_BED) {
    nqk += 2;
    for (i = 0; i < n_usr; i++) {
      qku[i] = 0.0;
    }
  }
  for (i = 0; i < n_usr; i++) {
    xusr[i] = usr_egrid[i]/(*te);
    if (usr_egrid_type == 1) xusr[i] += 1.0;
    log_xusr[i] = log(xusr[i]);
  }

  p = *qkc;
  for (i = 0; i < nz; i++) {
    kb = ang[i].kb;
    j0 = GetOrbital(kb)->kappa;
    j0 = GetJFromKappa(j0);
    for (ip = 0; ip <= i; ip++) {
      kbp = ang[ip].kb;
      j0p = GetOrbital(kbp)->kappa;
      j0p = GetJFromKappa(j0p);
      if (j0p != j0) continue;
      c = 16.0*ang[i].coeff*ang[ip].coeff;
      if (ip != i) {
	c *= 2.0;
      }
      if (nq == *nqkc) {
	*nqkc *= 2;
	*qkc = (double *) realloc(*qkc, sizeof(double)*(*nqkc));
	p = *qkc + nq*nqk;
      }
      ierr = CIRadialQkIntegrated(p, (*te), kb, kbp);
      if (ierr < 0) continue;
      if (qk_mode == QK_BED) {
	p[0] *= c;
	p[3] *= c;
      } else {
	for (j = 0; j < NPARAMS; j++) {
	  p[j] *= c;
	}
      }
      nq++;
      p += nqk;
    }
  }
  free(ang);

  if (qk_mode == QK_BED) {
    p = *qkc;  
    for (i = 0; i < nq; i++) {
      CIRadialQkFromFit(NPARAMS, p, n_usr, xusr, log_xusr, qke);
      for (j = 0; j < n_usr; j++) {
	qku[j] += qke[j]*xusr[j]/(xusr[j]+1.0);
      }
      p += nqk;
    }
  } else {
    p = *qkc+nqk;
    for (i = 1; i < nq; i++) {
      for (j = 0; j < NPARAMS; j++) {
	(*qkc)[j] += p[j];
      }
      p += nqk;
    }
    CIRadialQkFromFit(NPARAMS, *qkc, n_usr, xusr, log_xusr, qku);
    nq = 1;
  }

  return nq;
}

int SaveIonization(int nb, int *b, int nf, int *f, char *fn) {
  int i, j, k;
  int ie, ip;
  FILE *file;
  LEVEL *lev1, *lev2;
  CI_RECORD r;
  CI_HEADER ci_hdr;
  F_HEADER fhdr;
  double delta, emin, emax, e;
  double *qk, qku[MAXNUSR];
  int nqkc, nq, nqk, nshells;  

  emin = 1E10;
  emax = 1E-10;
  k = 0;
  for (i = 0; i < nb; i++) {
    lev1 = GetLevel(b[i]);
    for (j = 0; j < nf; j++) {
      lev2 = GetLevel(f[j]);
      e = lev2->energy - lev1->energy;
      if (e > 0) k++;
      if (e < emin && e > 0) emin = e;
      if (e > emax) emax = e;
    }
  }
  if (emin < TE_MIN_MAX*emax) {
    emin = TE_MIN_MAX*emax;
  }
  if (k == 0) {
    printf("No ionizations can occur\n");
    return 0;
  }
  if (qk_mode == QK_CB) {
    SetIEGrid(1, 0.5*(emin+emax), emax);
  } else {
    if (n_tegrid == 0) {
      n_tegrid = 3;
    }
    if (tegrid[0] < 0.0) {
      e = 2.0*(emax-emin)/(emax+emin);
      if (e < 0.1) {
	SetIEGrid(1, 0.5*(emin+emax), emax);
      } else if (e < 0.5) {
	SetIEGrid(2, emin, emax);
      } else {
	if (k == 2) n_tegrid = 2;
	SetIEGrid(n_tegrid, emin, emax);
      }
    }
  }

  e = 0.5*(emin + emax);
  if (qk_mode != QK_CB) {
    SetTransitionOptions(2, 1, 4, 4);
    SetRRTEGrid(1, e, e);
    SetPEGridLimits(-1, -1, 0);
    SetRecQkMode(QK_FIT, 0.1);
    SetPEGrid(0, -1.0, -1.0, 1.0);
    PrepRREGrids(e);
  }
  if (egrid_limits_type == 0) {
    emin = egrid_min*e;
    emax = egrid_max*e;
  } else {
    emin = egrid_min;
    emax = egrid_max;
  }
  egrid_type = 1;
  pw_type = 0;
  if (usr_egrid_type < 0) usr_egrid_type = 1;

  if (n_usr > 0 && usr_egrid[0] < 0.0) {
    SetUsrCIEGrid(n_usr, emin, emax, e);
  }    

  if (n_egrid == 0) {    
    n_egrid = 6;
  } else if (n_egrid <= NPARAMS) {
    printf("n_egrid should be at least %d\n", NPARAMS+1);
    exit(1);
  }
  if (egrid[0] < 0.0) {
    SetCIEGrid(n_egrid, emin, emax, e);
  }

  if (n_usr <= 0) {
    SetUsrCIEGridDetail(n_egrid, egrid);
    usr_egrid_type = 1;
  }
		  
  for (ie = 0; ie < n_egrid; ie++) {
    for (i = 0; i < n_tegrid; i++) {
      xegrid[i][ie] = egrid[ie]/tegrid[i];
      if (egrid_type == 1) xegrid[i][ie] += 1.0;
      log_xegrid[i][ie] = log(xegrid[i][ie]);
    }
    sigma[ie] = 1.0/sqrt(xegrid[0][ie]);
  }
  yegrid0[0] = 0.0;
  delta = 0.5/(NINT-1.0);
  for (i = 1; i < NINT; i++) {
    yegrid0[i] = yegrid0[i-1] + delta;
  }

  if (pw_scratch.nkl == 0) {
    SetCIPWGrid(0, NULL, NULL);
  }
  
  e = GetResidualZ();
  PrepCoulombBethe(n_egrid, n_tegrid, n_egrid, e, egrid, tegrid, egrid,
		   pw_scratch.nkl, pw_scratch.kl, egrid_type, pw_type, 0);  

  nqkc = 10;
  nqk = NPARAMS;
  if (qk_mode == QK_BED) nqk += 2;
  qk = (double *) malloc(sizeof(double)*nqkc*nqk);

  fhdr.type = DB_CI;
  strcpy(fhdr.symbol, GetAtomicSymbol());
  fhdr.atom = GetAtomicNumber();
  ci_hdr.nele = GetNumElectrons(b[0]);
  ci_hdr.qk_mode = qk_mode;
  ci_hdr.nparams = nqk;
  ci_hdr.pw_type = pw_type;
  ci_hdr.n_tegrid = n_tegrid;
  ci_hdr.n_egrid = n_egrid;
  ci_hdr.egrid_type = egrid_type;
  ci_hdr.n_usr = n_usr;
  ci_hdr.usr_egrid_type = usr_egrid_type;
  ci_hdr.tegrid = tegrid;
  ci_hdr.egrid = egrid;
  ci_hdr.usr_egrid = usr_egrid;
  file = OpenFile(fn, &fhdr);
  InitFile(file, &fhdr, &ci_hdr);

  nshells = 1;
  r.params = (float *) malloc(sizeof(float)*nqk);
  r.strength = (float *) malloc(sizeof(float)*n_usr);

  for (i = 0; i < nb; i++) {
    for (j = 0; j < nf; j++) {
      nq = IonizeStrength(qku, &nqkc, &qk, &e, b[i], f[j]);
      if (nq < 0) continue;
      r.b = b[i];
      r.f = f[j];
      r.nshells = nq;
      if (r.nshells > nshells) {
	r.params = (float *) realloc(r.params, sizeof(float)*nqk*r.nshells);
	nshells = r.nshells;
      }
      
      ip = 0;
      for (k = 0; k < nq; k++) {
	for (ie = 0; ie < nqk; ie++) {
	  r.params[ip] = (float) qk[ip];
	  ip++;
	}
      }
      
      for (ie = 0; ie < n_usr; ie++) {
	r.strength[ie] = (float) qku[ie];
      }

      WriteCIRecord(file, &r);
    }
  }

  DeinitFile(file, &fhdr);
  CloseFile(file, &fhdr);
  free(r.params);
  free(r.strength);
  free(qk);

  return 0;
}

void _FreeIonizationQk(void *p) {
  double *dp;
  dp = *((double **) p);
  free(dp);
}

int FreeIonizationQk(void) {
  ARRAY *b;
  b = qk_array->array;
  if (b == NULL) return 0;
  MultiFreeData(b, qk_array->ndim, _FreeIonizationQk);
  return 0;
}

int InitIonization(void) {
  int blocks[3] = {2, 10, 10};
  int ndim = 3;
  
  qk_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(qk_array, sizeof(double *), ndim, blocks);
  
  SetCIMaxK(8);
  n_egrid = 0;
  n_tegrid = 0;
  n_usr = 0;
  egrid[0] = -1.0;
  SetCIEGridLimits(-1.0, -1.0, 0);
  tegrid[0] = -1.0;
  usr_egrid[0] = -1.0;
  SetCIQkMode(QK_DEFAULT, 1E-3);
  return 0;
}

int ReinitIonization(int m) {

  if (m < 0) return 0;

  FreeIonizationQk(); 

  SetCIMaxK(8);
  n_egrid = 0;
  n_tegrid = 0;
  n_usr = 0;
  egrid[0] = -1.0;
  SetCIEGridLimits(-1.0, -1.0, 0);
  tegrid[0] = -1.0;
  usr_egrid[0] = -1.0;
  SetCIQkMode(QK_DEFAULT, 1E-3);

  return 0;
}  

