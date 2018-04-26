/*
 *   FAC - Flexible Atomic Code
 *   Copyright (C) 2001-2015 Ming Feng Gu
 * 
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 * 
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *   GNU General Public License for more details.
 * 
 *   You should have received a copy of the GNU General Public License
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "ionization.h"
#include "recombination.h"
#include "cf77.h"

static char *rcsid="$Id$";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

#define MAXMSUB 32
#define NINT0 12
#define NINT 128
#define YEG0 5E-4
#define YEG1 25.0
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
static int usr_different = 0;

static int n_tegrid = 0;
static double tegrid[MAXNTE];
static double log_te[MAXNTE];

static int xborn = 0;

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

void SetCIBorn(int x) {
  if (x >= 0) {
    xborn = x;
    SetCIQkMode(QK_DW, 1E-3);
  } else {
    xborn = 0;
  }
}

void FreeIonizationQkData(void *p) {
  double *dp;
  dp = *((double **) p);
  free(dp);
  *((double **) p) = NULL;
}

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

void SetCILQR(int m) {
  pw_scratch.qr = m;
}

void SetCILMax(int m) {
  pw_scratch.max_kl = m;
}

void SetCILMaxEject(int m) {
  pw_scratch.max_kl_eject = m;
}

void SetCILCB(int m) {
  pw_scratch.kl_cb = m;
}

void SetCITol(double t) {
  pw_scratch.tolerance = t;
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
  double bms, bte;

  BornFormFactorTE(&bte);
  bms = BornMass();
  eth = (eth + bte)/bms;
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
  double bms, bte;

  if (n > MAXNUSR) {
    printf("Max # of grid points reached \n");
    return -1;
  }
  BornFormFactorTE(&bte);
  bms = BornMass();
  eth = (eth + bte)/bms;
  n_usr = SetEGrid(usr_egrid, log_usr, n, emin, emax, eth);
  return n_usr;
}

int SetCIPWGrid(int ns, int *n, int *step) {
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

int CIRadialQk(double *qk, double e1, double e2, int kb, int kbp, int k) {
  double pk[MAXNTE][MAXNKL], pkp[MAXNTE][MAXNKL];
  double e0, te;
  ORBITAL *orb;
  int kappab, jb, klb, i, j, t, kl, klp;
  int js[4], ks[4];
  int jmin, jmax, ko2;
  int kf, kappaf, kf0, kappa0, kf1, kappa1;
  int kl0, kl0p, kl1, kl1p;
  int j0, j1, j1min, j1max;
  double z, z2, r, rp, sd, se, s;
  double eps, a, b, h, jb1;
  int kl_max0, kl_max1, kl_max2, kl_min2, max0;
  int type;
  int np = 3, one = 1;
  double logj;

  ko2 = k/2;

  for (i = 0; i < n_tegrid; i++) {
    for (j = 0; j < pw_scratch.nkl0; j++) {
      pk[i][j] = 0.0;
    }
  }
  for (i = 0; i < n_tegrid; i++) {
    qk[i] = 0.0;
  }
  
  r = GetRMax();
  z = GetResidualZ();
  z2 = z*z;
  t = r*sqrt(e1+2.0*z/r);
  kl_max0 = pw_scratch.max_kl;
  kl_max0 = Min(kl_max0, t);
  max0 = 12;

  orb = GetOrbital(kb);
  kappab = orb->kappa;
  GetJLFromKappa(kappab, &jb, &klb);
  klb /= 2;
  jmin = abs(k - jb);
  jmax = k + jb;
  jb1 = jb + 1.0;

  t = r*sqrt(e2+2.0*z/r);
  kl_max2 = pw_scratch.max_kl_eject + klb;
  kl_max2 = Min(kl_max2, t);
  kl_min2 = klb - pw_scratch.max_kl_eject;
  if (kl_min2 < 0) kl_min2 = 0;

  for (j = jmin; j <= jmax; j += 2) {
    for (klp = j - 1; klp <= j + 1; klp += 2) {
      kappaf = GetKappaFromJL(j, klp);
      kl = klp/2;
      if (kl > kl_max2) break;
      if (kl < kl_min2) continue;
      if (IsEven(kl + klb + ko2)) {
	type = ko2;
      } else {
        type = -1;
      }
      /*
      if (kl < pw_scratch.qr) {
	js[2] = 0;
      } else {
	js[2] = j;
	if (kappaf > 0) kappaf = -kappaf - 1;
      }
      */
      js[2] = 0;
      kf = OrbitalIndex(0, kappaf, e2);	
      ks[2] = kf;  

      if (xborn) {
	for (i = 0; i < n_tegrid; i++) {
	  CERadialQkBorn(kb, kf, kbp, kf, k, tegrid[i]+e2, e1, &r, 0);
	  qk[i] += r;
	}
	continue;
      }

      if (type >= 0 && type < CBMULT) {
	kl_max1 = pw_scratch.kl_cb;
      } else {
	kl_max1 = kl_max0;
      }

      for (i = 0; i < n_tegrid; i++) {
	for (t = 0; t < pw_scratch.nkl0; t++) {
	  pkp[i][t] = 0.0;
	}
      }
      for (t = 0; ; t++) {
	kl0 = pw_scratch.kl[t];
	if (kl0 > kl_max1) {
	  break;
	}	  
	kl0p = 2*kl0;
	for (j0 = abs(kl0p - 1); j0 <= kl0p + 1; j0 += 2) {
	  kappa0 = GetKappaFromJL(j0, kl0p);
	  if (kl0 < pw_scratch.qr) {
	    js[1] = 0;
	  } else {
	    js[1] = j0;
	    if (kappa0 > 0) kappa0 = -kappa0 - 1;
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
		if (kappa1 > 0) kappa1 = -kappa1 - 1;
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
		pk[i][t] += r;
		pkp[i][t] += r;
	      }
	    }
	  }
	}
      }
    }  
  }
 
  if (xborn == 0) {
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
  } else {
    for (i = 0; i < n_tegrid; i++) {
      qk[i] /= jb1;
    }
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

void CIRadialQkBasis0(int npar, double *yb, double x, double log_x) {
  double y1, y2;

  y1 = 1.0/x;
  y2 = 1.0 - y1;
  yb[0] = log_x;
  yb[1] = y2*y2;
  yb[2] = y2*y1 ;
  yb[3] = yb[1]*y1;
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
  
  for (i = 0; i < NINT; i++) {
    x[i] = 1.0/(2.0*yegrid0[i]);
    t[i] = log(x[i]);
  }

  for (i = 0; i < n_egrid; i++) {
    q[i] = log(q[i]);
  }
  for (i = 0; i < NINT; i++) {
    if (t[i] < logxe[n_egrid-1]) break;
  }
  if (i > 1) {
    d = te/p[3];
    for (j = 0; j < i; j++) {
      y[j] = (x[j]-1.0)*d + 1.0;
      s[j] = log(y[j]);
    }
    RRRadialQkFromFit(3, p, i, y, s, integrand, NULL, 0, &kl);
    for (j = 0; j < i; j++) {
      integrand[j] *= x[j]*d;
    }    
  }
  if (i < NINT) {
    n = NINT-i;
    np = 3;
    n1 = n_egrid;
    UVIP3P(np, n1, logxe, q, n, &(t[i]), &(integrand[i]));
    for (; i < NINT; i++) {
      integrand[i] = exp(integrand[i]);
    }
  }

  for (i = 0; i < NINT; i++) {
    integrand[i] *= x[i]*yegrid0[i];
  }
  
  d = 4.0*log(yegrid0[1]/yegrid0[0]);
  if (qk_mode == QK_DW) {
    (*bethe) = Simpson(integrand, 0, NINT-1);
    (*bethe) *= d;
  } else {
    n = (NINT+1)/2;
    j = n-1;
    t[j] = 1.0;
    y[j--] = 0.0;
    for (i = NINT-1; i > 1; i -= 2, j--) {
      y[j] = Simpson(integrand, i-2, i);
      y[j] *= d;
      y[j] += y[j+1];
      t[j] = 2.0*yegrid0[i-2];
    }
    *bethe = y[0];
    for (i = 0; i < NINT; i++) {
      integrand[i] *= x[i];
    }
    (*b) = Simpson(integrand, 0, NINT-1);
    (*b) *= d;
    for (i = 0; i < n_egrid; i++) {
      x0[i] = 1.0/xusr[i];
    }
    UVIP3P(np, n, t, y, n_egrid, x0, dp);
  }
  return 0;
}

double *CIRadialQkIntegratedTable(int kb, int kbp) {
  int index[2], ie, ite, i, j, k, nqk, qlog;
  double **p, *qkc, e1, e2;
  double yint[NINT], integrand[NINT];
  double ymin, ymax, dy, y, bte, bms;
  double yegrid[NINT0], qi[MAXNTE], qt[MAXNTE][NINT0];

  index[0] = kb;
  index[1] = kbp;
  LOCK *lock = NULL;
  int locked = 0;
  p = (double **) MultiSet(qk_array, index, NULL, &lock,
			   InitPointerData, FreeIonizationQkData);
  if (lock && !(*p)) {
    SetLock(lock);
    locked = 1;
  }
  if (*p) {
    if (locked) ReleaseLock(lock);
    return (*p);
  } 

  nqk = n_tegrid*n_egrid;
  double *pd = (double *) malloc(sizeof(double)*nqk);
  qkc = pd;
  
  bms = BornMass();
  BornFormFactorTE(&bte);
  for (ie = 0; ie < n_egrid; ie++) {
    for (ite = 0; ite < n_tegrid; ite++) {
      for (i = 0; i < NINT0; i++) {
	qt[ite][i] = 0.0;
      }
    }
    ymin = ((bte+tegrid[0])*YEG0)/egrid[ie];
    ymax = ((bte+tegrid[n_tegrid-1])*YEG1)/egrid[ie];
    if (ymax >= 0.99*bms) ymax = 0.99*bms;
    if (fabs(bms-1) < EPS3 && ymax > 0.5) ymax = 0.5;
    ymin = log(ymin);
    ymax = log(ymax);
    dy = (ymax - ymin)/(NINT0-1.0);
    yegrid[0] = ymin;
    for (i = 1; i < NINT0; i++) {
      yegrid[i] = yegrid[i-1] + dy;
    }
    dy = (ymax - ymin)/(NINT-1.0);
    yint[0] = ymin;
    for (i = 1; i < NINT; i++) {
      yint[i] = yint[i-1] + dy;
    }
    for (i = 0; i < NINT0; i++) {
      y = exp(yegrid[i]);
      e2 = egrid[ie]*y;
      e1 = egrid[ie]*(1.0-y/bms);
      for (k = 0; k <= pw_scratch.max_k; k += 2) {
	CIRadialQk(qi, e1, e2, kb, kbp, k);
	for (ite = 0; ite < n_tegrid; ite++) {
	  qi[ite] /= (k + 1.0);
	  qt[ite][i] += qi[ite];
	}
      }
    }
    for (ite = 0; ite < n_tegrid; ite++) {
      for (i = 0; i < NINT0; i++) {
	y = exp(yegrid[i]);
	qt[ite][i] *= y;
      }
      qlog = 1;
      for (i = 0; i < NINT0; i++) {
	if (qt[ite][i] <= 0.0) {
	  qlog = 0;
	  break;
	}
      }
      if (qlog) {
	for (i = 0; i < NINT0; i++) {
	  qt[ite][i] = log(qt[ite][i]);
	}
      }
      UVIP3P(3, NINT0, yegrid, qt[ite], NINT, yint, integrand);
      if (qlog) {
	for (i = 0; i < NINT; i++) {
	  integrand[i] = exp(integrand[i]);
	}
      }
      y = Simpson(integrand, 0, NINT-1);
      y *= dy*egrid[ie];
      i = ite*n_egrid + ie;
      qkc[i] = 16.0*y;
    }
  }
  *p = pd;
  if (locked) ReleaseLock(lock);
#pragma omp flush
  ReinitRadial(1);
  return (*p);
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
 
double BEScale(int k, double e) {
  double z, a, b, c;
  ORBITAL *orb;

  z = GetAtomicNumber();
  a = MeanPotential(k, k);
  b = RadialMoments(-1, k, k);
  c = -a/b;
  orb = GetOrbital(k);
  a = orb->energy - a;
  if (c >= z) c = z;  
  b = (1.0 + a/e);
  b /= orb->n;
  b *= 1 - 0.75*c/z;
  return b;
}

int IonizeStrengthUTA(double *qku, double *qkc, double *te, 
		      int b, int f) {
  INTERACT_DATUM *idatum;
  LEVEL *lev1, *lev2;
  int jb, kb, qb, klb, nq, ierr;
  double bethe, b0, c, cmax, qke[MAXNUSR], sigma[MAXNUSR];
  int ns, nqk, ip, i, j;
  double tol, x[MAXNE], logx[MAXNE], es;
    
  lev1 = GetLevel(b);
  lev2 = GetLevel(f);
  *te = lev2->energy - lev1->energy;
  if (*te <= 0) return -1;
  
  idatum = NULL;
  ns = GetInteract(&idatum, NULL, NULL, lev2->iham, lev1->iham,
		   lev2->pb, lev1->pb, 0, 0, 1);  
  if (ns <= 0) return -1;
  if (idatum->s[1].index < 0 || idatum->s[3].index >= 0) {
    free(idatum->bra);
    free(idatum);
    return -1;
  }
  jb = idatum->s[1].j;
  qb = idatum->s[1].nq_ket;

  if (qk_mode == QK_CB) {
    klb = GetLFromKappa(idatum->s[1].kappa);
    klb /= 2;
    nq = idatum->s[1].n;
    nqk = NPARAMS;
    for (j = 0; j < nqk; j++) {
      qkc[j] = 0.0;
    }
    ip = (nq - 1)*nq/2;
    for (j = 0; j < nqk; j++) {
      qkc[j] = (cbo_params[ip+klb][j]/(*te)) * (PI/2.0);
      qkc[j] *= qb*(lev1->ilev+1.0);
    }
    for (i = 0; i < n_usr; i++) {
      xusr[i] = usr_egrid[i]/(*te);
      if (usr_egrid_type == 1) xusr[i] += 1.0;
      log_xusr[i] = log(xusr[i]);
      qku[i] = 0.0;
    }
    CIRadialQkFromFit(NPARAMS, qkc, n_usr, xusr, log_xusr, qku);
    free(idatum->bra);
    free(idatum);
    return klb;
  } else {
    klb = BoundFreeOSUTA(qke, qkc, te, b, f, -1);
    for (i = 0; i < n_egrid; i++) {
      x[i] = (*te + egrid[i])/(*te);
      logx[i] = log(x[i]);
      xusr[i] = x[i];
    }
    CIRadialQkBED(qku, &bethe, &b0, klb, logx, qke, qkc, *te);
    if (qk_mode == QK_BED) {
      kb = OrbitalIndex(idatum->s[1].n, idatum->s[1].kappa, 0.0);
      es = BEScale(kb, *te);
      b0 = ((4.0*PI)/(*te))*qb*(lev1->ilev+1.0) - b0;      
      for (j = 0; j < n_egrid; j++) {
	qku[j] = qku[j]*logx[j] + 
	  b0*(1.0-1.0/x[j]-logx[j]/(1.0+x[j]));
	qku[j] *= (x[j]/(es+x[j]));
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
      if (usr_different) {
	for (i = 0; i < n_usr; i++) {
	  xusr[i] = usr_egrid[i]/(*te);
	  if (usr_egrid_type == 1) xusr[i] += 1.0;
	  log_xusr[i] = log(xusr[i]);
	}
	CIRadialQkFromFit(NPARAMS, qkc, n_usr, xusr, log_xusr, qku);
      }
    } else { 
      for (i = 0; i < n_egrid; i++) {
	qku[i] = 0.0;
      } 
      kb = OrbitalIndex(idatum->s[1].n, idatum->s[1].kappa, 0.0);
      ierr = CIRadialQkIntegrated(qke, *te, kb, kb);
      for (j = 0; j < n_egrid; j++) {
	qku[j] = qke[j]*qb*(lev1->ilev+1.0);
      }      
      ip = BornFormFactorTE(NULL);
      if (ip >= 0) {
	bethe = 0.0;
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
      if (ip < 0) {
	SVDFit(NPARAMS-1, qkc+1, NULL, tol, n_egrid, x, logx, 
	       qke, sigma, CIRadialQkBasis);
      } else {
	SVDFit(NPARAMS, qkc, NULL, tol, n_egrid, x, logx, 
	       qke, sigma, CIRadialQkBasis0);
      }
      if (usr_different) {
	for (i = 0; i < n_usr; i++) {
	  xusr[i] = usr_egrid[i]/(*te);
	  if (usr_egrid_type == 1) xusr[i] += 1.0;
	  log_xusr[i] = log(xusr[i]);
	}
	CIRadialQkFromFit(NPARAMS, qkc, n_usr, xusr, log_xusr, qku);
      }
    }
    free(idatum->bra);
    free(idatum);
    return klb;
  }
}
 
int IonizeStrength(double *qku, double *qkc, double *te, 
		   int b, int f) {
  int i, ip, j, ierr;
  LEVEL *lev1, *lev2;
  ORBITAL *orb;
  ANGULAR_ZFB *ang;
  double bethe, b0, c, c0, cmax, qke[MAXNUSR], sigma[MAXNUSR];
  int nz, j0, j0p, kl0, kl, kb, kbp, nq, nqk, kb0;
  double tol, x[MAXNE], logx[MAXNE], es;

  cmax = 0.0;
  c0 = 0.0;
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
    for (i = 0; i < n_egrid; i++) {
      x[i] = (*te + egrid[i])/(*te);
      logx[i] = log(x[i]);
      xusr[i] = x[i];
    }
    CIRadialQkBED(qku, &bethe, &b0, kl0, logx, qke, qkc, *te);
    if (qk_mode == QK_BED) {
      kb0 = -1;
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
	  } else if (c > c0) {
	    c0 = c;
	    kb0 = kb;
	  }
	  cmax += c;
	}
      }
      es = BEScale(kb0, *te);
      b0 = ((4.0*PI)/(*te))*cmax - b0;
      for (j = 0; j < n_egrid; j++) {
	c = b0*(1.0-1.0/x[j]-logx[j]/(1.0+x[j]));
	qku[j] = qku[j]*logx[j] + c;
	qku[j] *= x[j]/(es+x[j]);
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
      if (usr_different) {
	for (i = 0; i < n_usr; i++) {
	  xusr[i] = usr_egrid[i]/(*te);
	  if (usr_egrid_type == 1) xusr[i] += 1.0;
	  log_xusr[i] = log(xusr[i]);
	}
	CIRadialQkFromFit(NPARAMS, qkc, n_usr, xusr, log_xusr, qku);
      }
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
      ip = BornFormFactorTE(NULL);
      if (ip >= 0) {
	bethe = 0.0;
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
      if (ip < 0) {
	SVDFit(NPARAMS-1, qkc+1, NULL, tol, n_egrid, x, logx, 
	       qke, sigma, CIRadialQkBasis);
      } else {
	SVDFit(NPARAMS, qkc, NULL, tol, n_egrid, x, logx, 
	       qke, sigma, CIRadialQkBasis0);
      }
      if (usr_different) {
	for (i = 0; i < n_usr; i++) {
	  xusr[i] = usr_egrid[i]/(*te);
	  if (usr_egrid_type == 1) xusr[i] += 1.0;
	  log_xusr[i] = log(xusr[i]);
	}
	CIRadialQkFromFit(NPARAMS, qkc, n_usr, xusr, log_xusr, qku);
      }
    }
    free(ang);
    
    return kl0;
  }
}

int SaveIonization(int nb, int *b, int nf, int *f, char *fn) {
  int i, j, k;
  int ie, ip;
  TFILE *file;
  LEVEL *lev1, *lev2;
  CI_RECORD r;
  CI_HEADER ci_hdr;
  F_HEADER fhdr;
  double delta, emin, emax, e, emax0;
  double qk[MAXNE], qku[MAXNUSR];
  int nq, nqk;  
  ARRAY subte;
  int isub, n_tegrid0, n_egrid0, n_usr0;
  int te_set, e_set, usr_set, iuta;
  double c, e0, e1;

  iuta = IsUTA();

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
  if (k == 0) {
    return 0;
  }
  /*
#if USE_MPI == 2
  int mr, nr;
  mr = MPIRank(&nr);
  if (nr > 1) {
    RandIntList(nb, b);
    RandIntList(nf, f);
  }
#endif
  */
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

  if (egrid_limits_type == 0) {
    emax0 = 0.5*(emin + emax)*egrid_max;
  } else {
    emax0 = egrid_max;
  }

  ArrayInit(&subte, sizeof(double), 128);
  ArrayAppend(&subte, &emin, NULL);
  c = 1.0/TE_MIN_MAX;
  if (!e_set || !te_set) {
    e = c*emin;
    while (e < emax) {
      ArrayAppend(&subte, &e, NULL);
      e *= c;
    }
  }
  ArrayAppend(&subte, &emax, NULL);

  egrid_type = 1;
  pw_type = 0;
  if (usr_egrid_type < 0) usr_egrid_type = 1;
  nqk = NPARAMS;
    
  fhdr.type = DB_CI;
  strcpy(fhdr.symbol, GetAtomicSymbol());
  fhdr.atom = GetAtomicNumber();
  ci_hdr.nele = GetNumElectrons(b[0]);
  ci_hdr.qk_mode = qk_mode;
  ci_hdr.nparams = nqk;
  ci_hdr.pw_type = pw_type;
  ci_hdr.egrid_type = egrid_type;
  ci_hdr.usr_egrid_type = usr_egrid_type;
  file = OpenFile(fn, &fhdr);

  e0 = emin*0.999;
  for (isub = 1; isub < subte.dim; isub++) {
    e1 = *((double *) ArrayGet(&subte, isub));
    if (isub == subte.dim-1) e1 = e1*1.001;
    emin = e1;
    emax = e0;
    k = 0;
    for (i = 0; i < nb; i++) {
      lev1 = GetLevel(b[i]);
      for (j = 0; j < nf; j++) {
	lev2 = GetLevel(f[j]);
	e = lev2->energy - lev1->energy;
	if (e < e0 || e >= e1) continue;
	if (e < emin) emin = e;
	if (e > emax) emax = e;
	k++;
      }
    }
    if (k == 0) {
      e0 = e1;
      continue;
    }

    if (qk_mode == QK_CB) {
      SetIEGrid(1, 0.5*(emin+emax), emax);
    } else {
      if (n_tegrid0 == 0) {
	n_tegrid = 3;
      }
      if (!te_set) {
	e = 2.0*(emax-emin)/(emax+emin);
	if (e < EPS3) {
	  SetIEGrid(1, 0.5*(emin+emax), emax);
	} else if (e < 0.5) {
	  SetIEGrid(2, emin, emax);
	} else {
	  if (k == 2) n_tegrid = 2;
	  SetIEGrid(n_tegrid, emin, emax);
	}
      }
    }

    n_egrid = n_egrid0;
    n_usr = n_usr0;
    if (!usr_set) usr_egrid[0] = -1.0;
    if (!e_set) egrid[0] = -1.0;
    
    e = 0.5*(emin + emax);
    if (egrid_limits_type == 0) {
      emin = egrid_min*e;
      emax = egrid_max*e;
    } else {
      emin = egrid_min;
      emax = egrid_max;
    }
    if (emax < emax0) {
      emax = 50.0*e;
      if (emax > emax0) emax = emax0;
    }
    
    if (n_egrid <= 0) {    
      n_egrid = 6;
    }
    if (egrid[0] < 0.0) {
      SetCIEGrid(n_egrid, emin, emax, e);
    }
    
    usr_different = 1;
    if (n_usr > 0 && usr_egrid[0] < 0.0) {
      SetUsrCIEGrid(n_usr, emin, emax, e);
      usr_egrid_type = 1;
      usr_different = 0;
    }   
    if (n_usr <= 0) {
      SetUsrCIEGridDetail(n_egrid, egrid);
      usr_egrid_type = 1;
      usr_different = 0;
    } 
    
    if (qk_mode != QK_CB) {
      SetTransitionOptions(2, 1, 4, 4);
      SetRRTEGrid(1, e, e);
      SetPEGridLimits(egrid_min, egrid_max, egrid_limits_type);
      SetPEGridDetail(n_egrid, egrid);
      PrepRREGrids(e, emax0);
    }    
		  
    for (ie = 0; ie < n_egrid; ie++) {
      for (i = 0; i < n_tegrid; i++) {
	xegrid[i][ie] = egrid[ie]/tegrid[i];
	if (egrid_type == 1) xegrid[i][ie] += 1.0;
	log_xegrid[i][ie] = log(xegrid[i][ie]);
      }
      sigma[ie] = 1.0;
    }
    yegrid0[0] = log(1E-5);
    delta = (log(0.5) - yegrid0[0])/(NINT-1.0);
    for (i = 1; i < NINT; i++) {
      yegrid0[i] = yegrid0[i-1] + delta;
    }
    for (i = 0; i < NINT; i++) {
      yegrid0[i] = exp(yegrid0[i]);
    }

    if (pw_scratch.nkl == 0) {
      SetCIPWGrid(0, NULL, NULL);
    }

    ci_hdr.n_tegrid = n_tegrid;
    ci_hdr.n_egrid = n_egrid;
    ci_hdr.n_usr = n_usr;
    ci_hdr.tegrid = tegrid;
    ci_hdr.egrid = egrid;
    ci_hdr.usr_egrid = usr_egrid;
    InitFile(file, &fhdr, &ci_hdr);
    ResetWidMPI();
#pragma op parallel default(shared) private(i, j, lev1, lev2, e, nq, qku, qk, r, ip, ie)
    {
    r.strength = (float *) malloc(sizeof(float)*n_usr);
    r.params = (float *) malloc(sizeof(float)*nqk);
    for (i = 0; i < nb; i++) {
      lev1 = GetLevel(b[i]);
      for (j = 0; j < nf; j++) {
	lev2 = GetLevel(f[j]);
	e = lev2->energy - lev1->energy;
	if (e < e0 || e >= e1) continue;
	int skip = SkipMPI();
	if (skip) continue;
	if (iuta) {
	  nq = IonizeStrengthUTA(qku, qk, &e, b[i], f[j]);
	} else {
	  nq = IonizeStrength(qku, qk, &e, b[i], f[j]);
	}
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
    free(r.params);
    free(r.strength);
    }
    DeinitFile(file, &fhdr);

    ReinitRadial(1);
    FreeRecQk();
    FreeRecPk();
    FreeIonizationQk();
    
    e0 = e1;
  }

  ReinitRecombination(1);
  ReinitIonization(1);

  ArrayFreeLock(&subte, NULL);
  CloseFile(file, &fhdr);

  return 0;
}

double CIRadialQkMSub(int J0, int M0, int J1, int M1, int k0, int k1, 
		      double e1, double e2, double e0) {
  ORBITAL *orb0, *orb1;
  int jk0, jk1, t, kl1, j1, kappa1, k, ks[4];
  int Jmax, Jmin, kp, Jpmax, Jpmin, J, Jp;
  int j2, kl2, kappa2, j0max, j0min, j0pmax, j0pmin;
  int j0, kl0, kappa0, j0p, kl0p, kappa0p, m0, q, M, kmax, m1, m2;
  int j2max, j2min, j2pmax, j2pmin;
  double y[MAXNKL], r, rp, d, ph0, ph0p, sd, se, c, d1, d2;
  double w6j1, w6j2, w3j1, w3j2, w3j3, w3j4, w3j5, w3j6, w3j7, w3j8;

  orb0 = GetOrbital(k0);
  jk0 = GetJFromKappa(orb0->kappa);
  orb1 = GetOrbital(k1);
  jk1 = GetJFromKappa(orb1->kappa);

  for (t = 0; t < pw_scratch.nkl; t++) {
    y[t] = 0.0;
    kl1 = 2*pw_scratch.kl[t];
    for (j1 = kl1-1; j1 <= kl1+1; j1 += 2) {
      if (j1 < 0) continue;
      kappa1 = GetKappaFromJL(j1, kl1);
      ks[3] = OrbitalIndex(0, kappa1, e1);
      for (k = 0; k <= pw_scratch.max_k; k += 2) {
	for (kp = 0; kp <= pw_scratch.max_k; kp += 2) {
	  kmax = Min(k, kp);
	  j2max = jk0 + k;
	  j2min = abs(jk0 - k);
	  j2pmax = jk1 + kp;
	  j2pmin = abs(jk1 - kp);
	  j2max = Min(j2max, j2pmax);
	  j2min = Max(j2min, j2pmin);
	  for (j2 = j2min; j2 <= j2max; j2 += 2) {
	    for (kl2 = j2-1; kl2 <= j2+1; kl2 += 2) {
	      if (kl2/2 > pw_scratch.max_kl_eject) continue;
	      kappa2 = GetKappaFromJL(j2, kl2);
	      ks[2] = OrbitalIndex(0, kappa2, e2);
	      j0max = j1 + k;
	      j0min = abs(j1 - k);
	      j0pmax = j1 + kp;
	      j0pmin = abs(j1 - kp);
	      for (j0 = j0min; j0 <= j0max; j0 += 2) {
		for (kl0 = j0 - 1; kl0 <= j0 + 1; kl0 += 2) {
		  kappa0 = GetKappaFromJL(j0, kl0);
		  ks[1] = OrbitalIndex(0, kappa0, e0);
		  ph0 = GetPhaseShift(ks[1]);
		  ks[0] = k0;
		  SlaterTotal(&sd, &se, NULL, ks, k, 1);
		  r = sd + se;
		  for (j0p = j0pmin; j0p <= j0pmax; j0p += 2) {
		    for (kl0p = j0p - 1; kl0p <= j0p + 1; kl0p += 2) {
		      kappa0p = GetKappaFromJL(j0p, kl0p);
		      ks[1] = OrbitalIndex(0, kappa0p, e0);
		      ph0p = GetPhaseShift(ks[1]);
		      ks[0] = k1;
		      SlaterTotal(&sd, &se, NULL, ks, kp, 1);
		      rp = sd + se;
		      Jmin = abs(J0 - k);
		      Jmax = J0 + k;
		      j2pmin = abs(J1 - j2);
		      j2pmax = J1 + j2;
		      Jmin = Max(Jmin, j2pmin);
		      Jmax = Min(Jmax, j2pmax);
		      Jpmin = abs(J0 - kp);
		      Jpmax = J0 + kp;
		      Jpmin = Max(Jpmin, j2pmin);
		      Jpmax = Min(Jpmax, j2pmax);
		      c = 0.0;
		      d1 = sqrt((kl0+1.0)*(j0+1.0)*(kl0p+1.0)*(j0p+1.0));
		      for (J = Jmin; J <= Jmax; J += 2) {
			w6j1 = W6j(J1, j2, J, k, J0, jk0);
			for (Jp = Jpmin; Jp <= Jpmax; Jp += 2) {
			  w6j2 = W6j(J1, j2, Jp, kp, J0, jk1);
			  d2 = (J+1.0)*(Jp+1.0)*w6j1*w6j2;
			  for (m0 = -1; m0 <= 1; m0 += 2) {
			    w3j1 = W3j(j0, 1, kl0, -m0, m0, 0);
			    w3j2 = W3j(j0p, 1, kl0p, -m0, m0, 0);	
			    for (q = -kmax; q <= kmax; q += 2) {
			      m1 = m0 + q;
			      if (m1 > j1) continue;
			      M = M0 - q;
			      if (M > J || M > Jp) continue;
			      m2 = M-M1;
			      if (m2 > j2) continue;
			      w3j3 = W3j(j0, k, j1, -m0, -q, m1);
			      w3j5 = W3j(J0, k, J, -M0, q, M);
			      w3j7 = W3j(J1, j2, J, M1, m2, -M);
			      w3j4 = W3j(j0p, kp, j1, -m0, -q, m1);
			      w3j6 = W3j(J0, kp, Jp, -M0, q, M);
			      w3j8 = W3j(J1, j2, Jp, M1, m2, -M);
			      d = w3j1*w3j2*w3j3*w3j4*w3j5*w3j6*w3j7*w3j8;
			      d *= d1*d2;
			      if (IsOdd(abs(jk0-jk1)/2)) d = -d;
			      c += d;
			    }
			  }
			}
		      }
		      y[t] += c*r*rp*cos(ph0-ph0p);
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }

  r = y[0];
  for (t = 1; t < pw_scratch.nkl; t++) {
    r += y[t];
    kl0 = pw_scratch.kl[t-1];
    kl1 = pw_scratch.kl[t];
    for (kl2 = kl0+1; kl2 < kl1; kl2++) {
      rp = LnInteger(kl2); 
      UVIP3P(3, pw_scratch.nkl, pw_scratch.log_kl, y, 1, &rp, &d);
      r += d;
    }
  }

  r *= 16.0;
  
  return r;
}

double CIRadialQkIntegratedMSub(int j1, int m1, int j2, int m2,
				int k0, int k1, double te, double e12) {
  double e0, e1, e2, d, r;
  double x[NINT0], y[NINT0], xi[NINT], yi[NINT];
  double ymin, ymax, bms, bte;
  int i, qlog;

  bms = BornMass();
  BornFormFactorTE(&bte);
  e0 = te + e12;
  ymin = (YEG0*(bte+te))/e12;
  ymax = (YEG1*(bte+te))/e12;
  if (ymax > 0.99*bms) ymax = 0.99*bms;
  if (fabs(bms-1) < EPS3 && ymax > 0.5) ymax = 0.5;
  ymin = log(ymin);
  ymax = log(ymax);

  x[0] = ymin;
  d = (ymax - ymin)/(NINT0-1);
  for (i = 1; i < NINT0; i++) {
    x[i] = x[i-1] + d;
  }
  
  for (i = 0; i < NINT0; i++) {
    r = exp(x[i]);
    e2 = e12*r;
    e1 = e12 - e2/bms;
    y[i] = CIRadialQkMSub(j1, m1, j2, m2, k0, k1, e1, e2, e0);
    y[i] *= r;
  }

  xi[0] = ymin;
  d = (ymax - ymin)/(NINT-1.0);
  for (i = 1; i < NINT; i++) {
    xi[i] = xi[i-1] + d;
  }
  qlog = 1;
  for (i = 0; i < NINT0; i++) {
    if (y[i] <= 0) {
      qlog = 0;
      break;
    }
  }
  if (qlog) {
    for (i = 0; i < NINT0; i++) {
      y[i] = log(y[i]);
    }
  }
  UVIP3P(3, NINT0, x, y, NINT, xi, yi);
  if (qlog) {
    for (i = 0; i < NINT; i++) {
      yi[i] = exp(yi[i]);
    }
  }
  r = Simpson(yi, 0, NINT-1);
  r *= d*e12;
  
  ReinitRadial(1);
  return r;
}

int IonizeStrengthMSub(double *qku, double *te, int b, int f) {
  LEVEL *lev1, *lev2;
  ANGULAR_ZFB *ang;
  double c, d, x[MAXNE], logx[MAXNE];
  double qkc[MAXMSUB*MAXNE], *rqk;
  int j1, j2, m1, m2, nz, i, ip, j, ie, kb, kbp;
  
  lev1 = GetLevel(b);
  lev2 = GetLevel(f);
  *te = lev2->energy - lev1->energy;
  if (*te <= 0) return -1;
  
  j1 = lev1->pj;
  j2 = lev2->pj;
  DecodePJ(j1, NULL, &j1);
  DecodePJ(j2, NULL, &j2);

  for (i = 0; i < n_usr; i++) {
    xusr[i] = usr_egrid[i]/(*te);
    if (usr_egrid_type == 1) xusr[i] += 1.0;
    log_xusr[i] = log(xusr[i]);
  }
  for (i = 0; i < n_egrid; i++) {
    x[i] = (*te + egrid[i])/(*te);
    logx[i] = log(x[i]);
  }

  rqk = qkc;
  for (m1 = -j1; m1 <= 0; m1 += 2) {
    for (m2 = -j2; m2 <= j2; m2 += 2) {
      for (ie = 0; ie < n_egrid; ie++) {
	rqk[ie] = 0.0;
      }
      rqk += n_egrid;
    }
  }


  nz = AngularZFreeBound(&ang, f, b);
  for (i = 0; i < nz; i++) {
    kb = ang[i].kb;
    for (ip = 0; ip <= i; ip++) {
      kbp = ang[ip].kb;
      c = ang[i].coeff*ang[ip].coeff;
      if (ip != i) c *= 2.0;
      rqk = qkc;
      for (m1 = -j1; m1 <= 0; m1 += 2) {
	for (m2 = -j2; m2 <= j2; m2 += 2) {
	  for (ie = 0; ie < n_egrid; ie++) {
	    d = CIRadialQkIntegratedMSub(j1, m1, j2, m2,
					 kb, kbp, *te, egrid[ie]);
	    rqk[ie] += c*d;
	  }
	  rqk += n_egrid;
	}
      }
    }
  }
  if (nz <= 0) {
    return -1;
  }

  free(ang);
  rqk = qkc;
  i = 0;
  for (m1 = -j1; m1 <= 0; m1 += 2) {
    for (m2 = -j2; m2 <= j2; m2 += 2) {
      if (n_egrid > 1) {
	UVIP3P(3, n_egrid, logx, rqk, n_usr, log_xusr, qku);
      } else {
	for (ie = 0; ie < n_usr; ie++) {
	  qku[ie] = rqk[0];
	}
      }
      rqk += n_egrid;
      qku += n_usr;
      i++;
    }
  }
  
  return i;
}

int SaveIonizationMSub(int nb, int *b, int nf, int *f, char *fn) {
  TFILE *file;
  LEVEL *lev1, *lev2;
  CIM_RECORD r;
  CIM_HEADER ci_hdr;
  F_HEADER fhdr;
  double qku[MAXNUSR*MAXMSUB];
  double delta, emin, emax, e, emax0;
  int nq, i, j, k, ie;

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
  if (k == 0) {
    return 0;
  }
  /*
#if USE_MPI == 2
  int mr, nr;
  mr = MPIRank(&nr);
  if (nr > 1) {
    RandIntList(nb, b);
    RandIntList(nf, f);
  }
#endif
  */
  if (egrid_limits_type == 0) {
    emax0 = 0.5*(emin + emax)*egrid_max;
  } else {
    emax0 = egrid_max;
  }

  e = 0.5*(emin + emax);
  if (egrid_limits_type == 0) {
    emin = egrid_min*e;
    emax = egrid_max*e;
  } else {
    emin = egrid_min;
    emax = egrid_max;
  }
  if (emax < emax0) {
    emax = 50.0*e;
    if (emax > emax0) emax = emax0;
  }

  egrid_type = 1;
  if (usr_egrid_type < 0) usr_egrid_type = 1;
  if (n_egrid <= 0) n_egrid = 6;
  if (egrid[0] < 0.0) {
    SetCIEGrid(n_egrid, emin, emax, e);
  }
  if (n_usr > 0 && usr_egrid[0] < 0.0) {
    SetUsrCIEGrid(n_usr, emin, emax, e);
    usr_egrid_type = 1;
  }    
  if (n_usr <= 0) {
    SetUsrCIEGridDetail(n_egrid, egrid);
    usr_egrid_type = 1;
  }  
    
  fhdr.type = DB_CIM;
  strcpy(fhdr.symbol, GetAtomicSymbol());
  fhdr.atom = GetAtomicNumber();
  ci_hdr.nele = GetNumElectrons(b[0]);
  ci_hdr.egrid_type = egrid_type;
  ci_hdr.usr_egrid_type = usr_egrid_type;
  file = OpenFile(fn, &fhdr);

  yegrid0[0] = log(1E-5);
  delta = (log(0.5)-yegrid0[0])/(NINT-1.0);
  for (i = 1; i < NINT; i++) {
    yegrid0[i] = yegrid0[i-1] + delta;
  }
  for (i = 0; i < NINT; i++) {
    yegrid0[i] = exp(yegrid0[i]);
  }

  if (pw_scratch.nkl == 0) {
    SetCIPWGrid(0, NULL, NULL);
  }
    
  ci_hdr.n_egrid = n_egrid;
  ci_hdr.n_usr = n_usr;    
  ci_hdr.egrid = egrid;
  ci_hdr.usr_egrid = usr_egrid;
  InitFile(file, &fhdr, &ci_hdr);
  ResetWidMPI();
#pragma omp parallel default(shared) private(i, j, lev1, lev2, e, nq, r, qku, ie)
  {
  for (i = 0; i < nb; i++) {
    lev1 = GetLevel(b[i]);
    for (j = 0; j < nf; j++) {
      int skip = SkipMPI();
      if (skip) continue;
      lev2 = GetLevel(f[j]);
      e = lev2->energy - lev1->energy;
      nq = IonizeStrengthMSub(qku, &e, b[i], f[j]);
      if (nq < 0) continue;
      r.b = b[i];
      r.f = f[j];
      r.nsub = nq;
      r.strength = (float *) malloc(sizeof(float)*nq*n_usr);
      for (ie = 0; ie < nq*n_usr; ie++) {
	r.strength[ie] = qku[ie];
      }
      WriteCIMRecord(file, &r);
      free(r.strength);
    }
  }
  }
  DeinitFile(file, &fhdr);
  CloseFile(file, &fhdr);

  ReinitIonization(1);
  return 0;
}

int FreeIonizationQk(void) {
  MultiFreeData(qk_array, FreeIonizationQkData);
  return 0;
}

int InitIonization(void) {
  int blocks[2] = {MULTI_BLOCK2,MULTI_BLOCK2};
  int ndim = 2;
  
  qk_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(qk_array, sizeof(double *), ndim, blocks, "iqk_array");
  
  SetCIMaxK(IONMAXK);
  n_egrid = 0;
  n_tegrid = 0;
  n_usr = 0;
  egrid[0] = -1.0;
  SetCIEGridLimits(-1.0, -1.0, 0);
  tegrid[0] = -1.0;
  usr_egrid[0] = -1.0;
  SetCIQkMode(QK_DEFAULT, 1E-3);
  SetCIPWOptions(IONLQR, IONLMAX, IONLEJEC, IONLCB, IONTOL);

  return 0;
}

int ReinitIonization(int m) {

  if (m < 0) return 0;

  FreeIonizationQk(); 

  n_egrid = 0;
  n_tegrid = 0;
  n_usr = 0;
  egrid[0] = -1.0;
  SetCIEGridLimits(-1.0, -1.0, 0);
  tegrid[0] = -1.0;
  usr_egrid[0] = -1.0;
  /*
  SetCIQkMode(QK_DEFAULT, 1E-3);
  SetCIPWOptions(IONLQR, IONLMAX, IONLEJEC, IONLCB, IONTOL);
  */
  return 0;
}  

