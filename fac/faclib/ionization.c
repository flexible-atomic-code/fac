#include "ionization.h"
#include "cf77.h"

static char *rcsid="$Id: ionization.c,v 1.40 2003/03/11 21:22:45 mfgu Exp $";
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
} pw_scratch = {IONMAXK, IONLQR, IONLMAX, IONLEJEC, 
		IONLCB, IONTOL, 0, 0};

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
  if (pw_scratch.nkl0 <= 0) SetCIPWOptions(IONLQR, IONLMAX, 
					   IONLEJEC, IONLCB, IONTOL);
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
  int np = 3, one = 1;
  double logj;

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
	      if (b < 10.0) {
		r = qkl * b;
	      } else {
		r = qkl * 10.0;
	      }
	      qk[i] += r;  
	    }
	  } else if (type >= 0) {
	    for (i = 0; i < n_tegrid; i++) {
	      b = (GetCoulombBethe(ie2, i, ie1, type, 1))[t];
	      if (b < 50.0) {
		r = qkl*b;
	      } else {
		r = qkl*50.0;
	      }
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
    r = pk[i][0];
    for (j = 1; j < t; j++) {
      r += pk[i][j];
      kl0 = pw_scratch.kl[j-1];
      kl1 = pw_scratch.kl[j];
      for (kl = kl0+1; kl < kl1; kl++) {
	logj = LnInteger(kl);
	UVIP3P(np, t, pw_scratch.log_kl, pk[i], one, &logj, &s);
	r += s;
      }
    } 
    qk[i] += r;
    qk[i] /= jb1;
  }

  return 0;
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
  double a, c, d;
  int i;

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

int CIRadialQkBED(double *dp, double *bethe, double *b, int kl,
		  double *logxe, double *q, double *p, double te) {
  double integrand[NINT];
  double x[NINT], y[NINT], t[NINT], s[NINT], d;
  double x0[MAXNUSR];
  int i, j, n, np, n1;
  
  for (i = 1; i < NINT; i++) {
    x[i] = 1.0/(2.0*yegrid0[i]);
    t[i] = log(x[i]);
  }

  for (i = 0; i < n_egrid; i++) {
    q[i] = log(q[i]);
  }
  for (i = 1; i < NINT; i++) {
    if (t[i] < logxe[i]) break;
  }
  if (i > 1) {
    d = te/p[3];
    for (j = 1; j < i; j++) {
      y[j] = (x[j]-1.0)*d + 1.0;
      s[j] = log(y[j]);
    }
    RRRadialQkFromFit(3, p, i-1, &(y[1]), &(s[1]), 
		      &(integrand[1]), NULL, 0, &kl);
    for (j = 1; j < i; j++) {
      integrand[j] *= x[j]*d;
    }    
  }
  n = NINT-i;
  np = 3;
  n1 = n_egrid;
  UVIP3P(np, n1, logxe, q, n, &(t[i]), &(integrand[i]));
  for (; i < NINT; i++) {
    integrand[i] = exp(integrand[i]);
  }
  
  integrand[0] = 0.0;
  for (i = 1; i < NINT; i++) {
    integrand[i] *= x[i];
  }

  d = 4.0*yegrid0[1];
  if (qk_mode == QK_DW) {
    (*bethe) = Simpson(integrand, 0, NINT-1);
    (*bethe) *= d;
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
    *bethe = y[0];
    for (i = 1; i < NINT; i++) {
      integrand[i] *= x[i];
    }
    (*b) = Simpson(integrand, 0, NINT-1);
    (*b) *= d;
    for (i = 0; i < n_usr; i++) {
      x0[i] = 1.0/xusr[i];
    }
    UVIP3P(np, n, t, y, n_usr, x0, dp);
  }
  return 0;
}

double *CIRadialQkIntegratedTable(int kb, int kbp) {
  int index[3], ie1, ie2, ite, i, j;
  double **p, *qkc, *qk;
  double xr[MAXNQK], yr[MAXNQK], zr[MAXNE];
  double x, y, integrand[NINT];
  int np, nd, ierr, nqk;

  index[0] = 0;
  index[1] = kb;
  index[2] = kbp;
  
  p = (double **) MultiSet(qk_array, index, NULL);
  if (*p) {
    return (*p);
  } 
  nqk = n_tegrid*n_egrid;
  *p = (double *) malloc(sizeof(double)*nqk);
  qkc = *p;

  qk = CIRadialQkTable(kb, kbp);
  if (qk == NULL) return NULL;

  nd = 1;
  for (ite = 0; ite < n_tegrid; ite++) {
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
    for (ie1 = 0; ie1 < n_egrid; ie1++) {
      for (i = 0; i < NINT; i++) {	
	x = log(1.0 + egrid[ie1]*(1.0-yegrid0[i])/tegrid[ite]);
	y = log(1.0 + egrid[ie1]*yegrid0[i]/tegrid[ite]);
	if (np == 0) np = 1;
	else np = 3;
	SDBI3P(np, nqk, xr, yr, qk, nd, &x, &y, &integrand[i],
		&ierr, _dwork, _iwork);
	if (kb == kbp) {
	  integrand[i] = exp(integrand[i]);
	}
      }   
      x = Simpson(integrand, 0, NINT-1);
      x *= yegrid0[1]*egrid[ie1];
      qkc[ie1] = 16.0*x;
    }
    qk += nqk;
    qkc += n_egrid;
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
 
int CIRadialQkIntegrated(double *qke, double te, int kb, int kbp) {
  int np, i, j, k, nd, nq;
  double *qk, qkc[MAXNTE];

  qk = CIRadialQkIntegratedTable(kb, kbp);
  if (qk == NULL) {
    return -1;
  }

  nq = n_egrid;
  if (n_tegrid > 1) {
    nd = 1;
    np = 3;
    for (i = 0; i < nq; i++) {
      k = i;
      for (j = 0; j < n_tegrid; j++) {
	qkc[j] = qk[k];
	k += n_egrid;
      }
      UVIP3P(np, n_tegrid, tegrid, qkc, nd, &te, &qke[i]);
    }
  } else {
    for (i = 0; i < nq; i++) {
      qke[i] = qk[i];
    }
  }
  return 0;
}
 
int IonizeStrength(double *qku, double *qkc, double *te, 
		   int b, int f) {
  int i, ip, j, ierr;
  LEVEL *lev1, *lev2;
  ORBITAL *orb;
  ANGULAR_ZFB *ang;
  double bethe, b0, c, cmax, qke[MAXNUSR], sigma[MAXNUSR];
  int nz, j0, j0p, kl0, kl, kb, kbp, nq, nqk;
  double tol, x[MAXNE], logx[MAXNE], es;

  c = GetResidualZ()-1.0;
  es = GetAtomicNumber();
  es = (es - c)/(es + c);
  cmax = 0.0;
  if (qk_mode == QK_CB) {
    nqk = NPARAMS;
    lev1 = GetLevel(b);
    lev2 = GetLevel(f);
    *te = lev2->energy - lev1->energy;
    if (*te <= 0) return -1;
    nz = AngularZFreeBound(&ang, f, b);
    if (nz <= 0) return -1;
    for (i = 0; i < n_usr; i++) {
      xusr[i] = usr_egrid[i]/(*te);
      if (usr_egrid_type == 1) xusr[i] += 1.0;
      log_xusr[i] = log(xusr[i]);
      qku[i] = 0.0;
    }
    for (j = 0; j < nqk; j++) {
      qkc[j] = 0.0;
    }
    for (i = 0; i < nz; i++) {
      kb = ang[i].kb;
      orb = GetOrbital(kb);
      kb = orb->kappa;
      GetJLFromKappa(kb, &j0, &kl);
      kl /= 2;
      c = ang[i].coeff*ang[i].coeff;	
      if (c > cmax) {
	kl0 = kl;
	nq = orb->n;
	cmax = c;
      }
      ip = (orb->n - 1)*orb->n/2;
      for (j = 0; j < nqk; j++) {
	qke[j] = (cbo_params[ip+kl][j]/(*te)) * (PI/2.0);
	qkc[j] += c*qke[j];
      }
    }
    CIRadialQkFromFit(NPARAMS, qkc, n_usr, xusr, log_xusr, qku);
    free(ang);
    return kl0;
  } else {
    kl0 = BoundFreeOS(qke, qkc, te, b, f, -1);
    if (kl0 < 0) return kl0;
    nz = AngularZFreeBound(&ang, f, b);
    for (i = 0; i < n_usr; i++) {
      xusr[i] = usr_egrid[i]/(*te);
      if (usr_egrid_type == 1) xusr[i] += 1.0;
      log_xusr[i] = log(xusr[i]);
    }
    for (i = 0; i < n_egrid; i++) {
      x[i] = (*te + egrid[i])/(*te);
      logx[i] = log(x[i]);
    }
    CIRadialQkBED(qku, &bethe, &b0, kl0, logx, qke, qkc, *te);
    if (qk_mode == QK_BED) {
      for (i = 0; i < nz; i++) {
	kb = ang[i].kb;
	j0 =GetOrbital(kb)->kappa;
	j0 = GetJFromKappa(j0);
	for (ip = 0; ip <= i; ip++) {
	  kbp = ang[ip].kb;
	  j0p = GetOrbital(kbp)->kappa;
	  j0p = GetJFromKappa(j0p);
	  if (j0p != j0) continue;
	  c = ang[i].coeff*ang[ip].coeff;
	  if (ip != i) {
	    c *= 2.0;
	  }
	  cmax += c;
	}
      }
      b0 = ((4.0*PI)/(*te))*cmax - b0;
      for (j = 0; j < n_usr; j++) {
	qku[j] = qku[j]*log_xusr[j] + 
	  b0*(1.0-1.0/xusr[j]-log_xusr[j]/(1.0+xusr[j]));
	qku[j] *= xusr[j]/(es+xusr[j]);
      }
      for (i = 0; i < n_usr; i++) {
	qke[i] = qku[i] - bethe*log_xusr[i];
	sigma[i] = qke[i];
      }
      qkc[0] = bethe;
      for (i = 1; i < NPARAMS; i++) {
	qkc[i] = 0.0;
      }
      tol = qk_fit_tolerance;
      SVDFit(NPARAMS-1, qkc+1, NULL, tol, n_usr, xusr, log_xusr, 
	     qke, sigma, CIRadialQkBasis);
    } else {  
      for (i = 0; i < n_egrid; i++) {
	qku[i] = 0.0;
      }
      for (i = 0; i < nz; i++) {
	kb = ang[i].kb;
	j0 =GetOrbital(kb)->kappa;
	j0 = GetJFromKappa(j0);
	for (ip = 0; ip <= i; ip++) {
	  kbp = ang[ip].kb;
	  j0p = GetOrbital(kbp)->kappa;
	  j0p = GetJFromKappa(j0p);
	  if (j0p != j0) continue;
	  c = ang[i].coeff*ang[ip].coeff;
	  if (ip != i) {
	    c *= 2.0;
	  } 
	  ierr = CIRadialQkIntegrated(qke, (*te), kb, kbp);
	  if (ierr < 0) continue;
	  for (j = 0; j < n_egrid; j++) {
	    qku[j] += c*qke[j];
	  }
	}
      }
      for (i = 0; i < n_egrid; i++) {
	qke[i] = qku[i] - bethe*logx[i];
	sigma[i] = qke[i];
      }
      qkc[0] = bethe;
      for (i = 1; i < NPARAMS; i++) {
	qkc[i] = 0.0;
      }
      tol = qk_fit_tolerance;
      SVDFit(NPARAMS-1, qkc+1, NULL, tol, n_egrid, x, logx, 
	     qke, sigma, CIRadialQkBasis);
      CIRadialQkFromFit(NPARAMS, qkc, n_usr, xusr, log_xusr, qku);
    }
    free(ang);
    
    return kl0;
  }
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
  double qk[MAXNE], qku[MAXNUSR];
  int nq, nqk;  

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

  if (qk_mode != QK_CB) {
    SetTransitionOptions(2, 1, 4, 4);
    SetRRTEGrid(1, e, e);
    SetPEGridLimits(egrid_min, egrid_max, egrid_limits_type);
    SetPEGridDetail(n_egrid, egrid);
    PrepRREGrids(e);
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
    sigma[ie] = 1.0;
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

  nqk = NPARAMS;
  r.params = (float *) malloc(sizeof(float)*nqk);
  r.strength = (float *) malloc(sizeof(float)*n_usr);
  
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

  for (i = 0; i < nb; i++) {
    for (j = 0; j < nf; j++) {
      nq = IonizeStrength(qku, qk, &e, b[i], f[j]);
      if (nq < 0) continue;
      r.b = b[i];
      r.f = f[j];
      r.kl = nq;
     
      for (ip = 0; ip < nqk; ip++) {
	r.params[ip] = (float) qk[ip];
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

  return 0;
}

void FreeIonizationQkData(void *p) {
  double *dp;
  dp = *((double **) p);
  free(dp);
  *((double **) p) = NULL;
}

int FreeIonizationQk(void) {
  ARRAY *b;
  b = qk_array->array;
  if (b == NULL) return 0;
  MultiFreeData(b, qk_array->ndim, FreeIonizationQkData);
  return 0;
}

int InitIonization(void) {
  int blocks[3] = {MULTI_BLOCK3,MULTI_BLOCK3,MULTI_BLOCK3};
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

