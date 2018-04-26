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

#include "excitation.h"
#include "transition.h"
#include "cf77.h"

static char *rcsid="$Id$";
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
static int n_egrid1 = 0;
static double egrid[MAXNE+2];
static double log_egrid[MAXNE+2];
static double egrid_min;
static double egrid_max;
static int egrid_limits_type = 0;

static int n_tegrid = 0;
static double tegrid[MAXNTE];
static double log_te[MAXNTE];

static int n_thetagrid = 0;
static double thetagrid[MAXNTHETA];
static int n_phigrid = 0;
static double phigrid[MAXNPHI];

#define NKINT 256
static double kgrid[NKINT];
static double log_kgrid[NKINT];
static double kint[NKINT];
static double log_kint[NKINT];
static double gos1[NKINT];
static double gos2[NKINT];
static double gost[NKINT];
static double gosint[NKINT];
static double xborn = XBORN;
static double xborn0 = XBORN0;
static double xborn1 = XBORN1;
static double eborn = EBORN;
#pragma omp threadprivate (kgrid, log_kgrid, kint, log_kint, gos1, gos2, gost, gosint)

static FILE *fpw=NULL;

#ifdef PERFORM_STATISTICS
static EXCIT_TIMING timing = {0, 0, 0};
#endif

static CEPW_SCRATCH pw_scratch = {1, MAXKL, 100, 5E-2, 0, 0, 10};

static int maxcecache = MAXCECACHE;
static CECACHE cecache = {0};

static MULTI *pk_array;
static MULTI *qk_array;
static MULTI *qkm_array;

static void InitCEPK(void *p, int n) {
  CEPK *d;
  int i;
  
  d = (CEPK *) p;
  for (i = 0; i < n; i++) {
    d[i].nkl = -1;
  }
}

void SetMaxCECache(int n) {
  if (n > 0) {
    maxcecache = n;
  } else {
    maxcecache = MAXCECACHE;
  }
  FreeCECache(0);
}

void AllocCECache(int msub) {
  cecache.msub = msub;
  cecache.low = malloc(sizeof(int)*maxcecache);
  cecache.up = malloc(sizeof(int)*maxcecache);
  cecache.nz = malloc(sizeof(int)*maxcecache);
  cecache.az = malloc(sizeof(ANGULAR_ZMIX *)*maxcecache);
  cecache.nmk = malloc(sizeof(int)*maxcecache);
  cecache.mbk = malloc(sizeof(double *)*maxcecache);
  int i;
  for (i = 0; i < maxcecache; i++) {
    cecache.low[i] = -1;
    cecache.up[i] = -1;
    cecache.nz[i] = 0;
    cecache.az[i] = NULL;
    cecache.nmk[i] = 0;
    cecache.mbk[i] = NULL;
  }
  cecache.ks = NULL;
  cecache.kd0 = NULL;
  cecache.kd1 = NULL;
  cecache.pk = NULL;
  cecache.qk = NULL;
}

void FreeCEQKK(CEQKK *qk, int msub) {
  int nk, i;
  nk = qk->kmax-qk->kmin+1;
  if (msub) nk *= qk->kmaxp-qk->kminp+1;
  for (i = 0; i < nk; i++) {
    free(qk->qk[i]);
  }
  free(qk->qk);
  free(qk);
}

void FreeCEPKK(CEPKK *pk) {
  int nk, i;
  CEPK *dp;
  nk = pk->kmax-pk->kmin+1;
  for (i = 0; i < nk; i++) {
    dp = pk->pk[i];
    if (dp->nkl > 0) {
      free(dp->kappa0);
      free(dp->kappa1);
      free(dp->pkd);
      free(dp->pke);
    }
    free(dp);
  }
  free(pk->pk);
  free(pk);
}

void FreeCECache(int m) {
  int i;
  for (i = 0; i < cecache.nc; i++) {
    if (cecache.nz[i] > 0) {
      if (cecache.az[i]) {
	free(cecache.az[i]);
	cecache.az[i] = NULL;
      }
      cecache.nz[i] = 0;
    }
    if (cecache.nmk[i] > 0) {
      if (cecache.mbk[i]) {
	free(cecache.mbk[i]);
	cecache.mbk[i] = NULL;
      }
      cecache.nmk[i] = 0;
    }
  }
  if (cecache.ks) {
    FreeIdxAry(cecache.ks, 0);
    free(cecache.ks);
    cecache.ks = NULL;
  }
  if (cecache.kd0) {
    FreeIdxAry(cecache.kd0, 0);
    free(cecache.kd0);
    cecache.kd0 = NULL;
  }
  if (cecache.kd1) {
    FreeIdxAry(cecache.kd1, 0);
    free(cecache.kd1);
    cecache.kd1 = NULL;
  }
  if (cecache.pk) {
    FreeCEPKK(cecache.pk);
    cecache.pk = NULL;
  }
  if (cecache.qk) {
    FreeCEQKK(cecache.qk, cecache.msub);
    cecache.qk = NULL;
  }
  if (m == 0) {
    free(cecache.low);
    free(cecache.up);
    free(cecache.nz);
    free(cecache.az);
    free(cecache.nmk);
    free(cecache.mbk);
  }
  cecache.nc = 0;
}

void FreeExcitationPkData(void *p) {
  CEPK *dp;
  
  dp = (CEPK *)p;
  if (dp->nkl > 0) {
    free(dp->kappa0);
    free(dp->kappa1);
    free(dp->pkd);
    free(dp->pke);
  }
  dp->nkl = -1;
}

void FreeExcitationQkData(void *p) {
  double *dp;

  dp = *((double **) p);
  free(dp);
  *((double **) p) = NULL;
}

CEPW_SCRATCH *GetCEPWScratch(void) {
  return &pw_scratch;
}

int SetCEPWFile(char *fn) {
  fpw = fopen(fn, "w");
  return 0;
}

int SetCEQkMode(int m, double tol) {
  if (m == QK_DEFAULT) qk_mode = QK_EXACT;
  else qk_mode = m;
  if (tol > 0.0) qk_fit_tolerance = tol;
  return 0;
}

int SetCEBorn(double eb, double x, double x1, double x0) {
  xborn = x;
  xborn1 = x1;
  if (x0 > 0) {
    xborn0 = x0;
  } else {
    xborn0 = XBORN0;
  }
  if (eb > 0) {
    eborn = eb;
  } else if (eb < 0) {
    eborn = eb/HARTREE_EV;
  } else {
    eborn = EBORN;
  }
  
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
  double bms, bte;

  BornFormFactorTE(&bte);
  bms = BornMass();
  eth = (eth + bte)/bms;
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

int SetAngleGridDetail(int m, int n, double *xg) {
  int i;

  if (m == 0) {
    if (n > MAXNTHETA) {
      printf("Max # of grid points reached\n");
      return -1;
    }
    
    n_thetagrid = n;
    for (i = 0; i < n; i++) {
      thetagrid[i] = xg[i];
    }
  } else {
    if (n > MAXNPHI) {
      printf("Max # of grid points reached\n");
      return -1;
    }
    
    n_phigrid = n;
    for (i = 0; i < n; i++) {
      phigrid[i] = xg[i];
    }
  }
  return n;
}

int SetAngleGrid(int m, int n, double xmin, double xmax) {
  if (m == 0) {
    n_thetagrid = SetLinearGrid(thetagrid, n, xmin, xmax);    
  } else {
    n_phigrid = SetLinearGrid(phigrid, n, xmin, xmax);
  }
  
  return n;
}

void SetCELQR(int m) {
  pw_scratch.qr = m;
}

void SetCELMax(int m) {
  pw_scratch.max_kl = m;
}

void SetCELCB(int m) {
  pw_scratch.kl_cb = m;
}

void SetCETol(double t) {
  pw_scratch.tolerance = t;
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
  pw_scratch.nkl = SetPWGrid(&(pw_scratch.nkl0),
			     pw_scratch.kl,
			     pw_scratch.log_kl,
			     pw_scratch.max_kl,
			     &ns, n, step);
  pw_scratch.ns = ns;  
  return 0;
}

int CERadialPk(CEPK **pk, int ie, int k0, int k1, int k, int trylock) {
  int type, ko2, i, m, t, q;
  int kf0, kf1, kpp0, kpp1, km0, km1;
  int kl0, kl1, kl0p, kl1p;
  int j0, j1, kl_max, j1min, j1max;
  ORBITAL *orb0, *orb1;
  int index[3];
  double te, e0, e1, sd, se;
  double a, tdi[MAXNTE], tex[MAXNTE];
  int js1, js3, js[4], ks[4];
  int nkappa, noex[MAXNTE];
  short *kappa0, *kappa1;
  double *pkd, *pke;

#ifdef PERFORM_STATISTICS
  clock_t start, stop;
  start = clock();
#endif

  ko2 = k/2;  
  index[0] = ko2*MAXNE + ie;
  index[1] = k0;
  index[2] = k1;    
  
  type = -1;
  orb0 = GetOrbital(k0);
  orb1 = GetOrbital(k1);
  GetJLFromKappa(orb0->kappa, &j0, &kl0);
  GetJLFromKappa(orb1->kappa, &j1, &kl1);
  kl0 = kl0/2;
  kl1 = kl1/2;
  if (IsEven(kl0 + kl1 + ko2) && Triangle(j0, j1, k)) {
    type = ko2;
  }
  LOCK *lock = NULL;
  int locked = 0;
  *pk = (CEPK *) MultiSet(pk_array, index, NULL, &lock,
			  InitCEPK, FreeExcitationPkData);
  if (lock && (*pk)->nkl < 0) {
    if (trylock) {
      if (0 == TryLock(lock)) {
	locked = 1;
      } else {
	return -9999;
      }
    } else {
      SetLock(lock);
      locked = 1;
    }
  }
  if ((*pk)->nkl >= 0) {
    if (locked) ReleaseLock(lock);
    return type;
  }

  nkappa = (MAXNKL)*(GetMaxRank()+1)*4;
  kappa0 = (short *) malloc(sizeof(short)*nkappa);
  kappa1 = (short *) malloc(sizeof(short)*nkappa);
  pkd = (double *) malloc(sizeof(double)*(nkappa*n_tegrid));
  pke = (double *) malloc(sizeof(double)*(nkappa*n_tegrid));

  e1 = egrid[ie];
  if (type > 0 && type <= CBMULT) {
    kl_max = pw_scratch.kl_cb;
  } else {
    kl_max = pw_scratch.max_kl;
  }
  
  js[0] = 0;
  ks[0] = k0;
  js[2] = 0;
  ks[2] = k1;		

  q = 0;
  m = 0;
  if (pw_type == 0) {
    js1 = 1;
    js3 = 3;
  } else {
    js1 = 3;
    js3 = 1;
  }

  for (i = 0; i < n_tegrid; i++) {
    noex[i] = 0;
    tdi[i] = 0.0;
    tex[i] = 0.0;
  }
  for (t = 0; t < pw_scratch.nkl; t++) {
    kl0 = pw_scratch.kl[t];
    if (pw_scratch.kl[t] > kl_max) break;
    kl0p = 2*kl0;
    for (i = 0; i < n_tegrid; i++) {
      if (noex[i] == 0) {
	if (1+tex[i] != 1) {
	  a = fabs(1.0 - tdi[i]/tex[i]);
	  if (a < EPS10) noex[i] = 1;
	}
      }
      if (noex[i] == 0) {
	tdi[i] = 0.0;
	tex[i] = 0.0;
      }
    }
    for (j0 = abs(kl0p-1); j0 <= kl0p+1; j0 += 2) {
      kpp0 = GetKappaFromJL(j0, kl0p); 
      km0 = kpp0;
      if (kl0 < pw_scratch.qr) {
	js[js1] = 0;
      } else {
	js[js1] = j0;
	if (kpp0 > 0) km0 = -kpp0 - 1;
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
	    if (kpp1 > 0) km1 = -kpp1 - 1;
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
	    e0 = e1 + te;
	    if (pw_type == 0) {
	      kf0 = OrbitalIndex(0, km0, e0);
	      ks[1] = kf0;
	    } else {
	      kf0 = OrbitalIndex(0, km1, e0);
	      ks[1] = kf0;	      
	    }
	    
	    if (noex[i] == 0) {
	      if (kl1 >= pw_scratch.qr &&
		  kl0 >= pw_scratch.qr) {
		SlaterTotal(&sd, &se, js, ks, k, -1);
	      } else {
		SlaterTotal(&sd, &se, js, ks, k, 1);
	      }
	      if (i == 0) {
		if (1.0+sd == 1.0 && 1.0+se == 1.0) {
		  break;
		}
	      }
	      tdi[i] += sd*sd;
	      tex[i] += (sd+se)*(sd+se);
	    } else {
	      se = 0.0;
	      if (kl1 >= pw_scratch.qr &&
		  kl0 >= pw_scratch.qr) {
		SlaterTotal(&sd, NULL, js, ks, k, -1);
	      } else {
		SlaterTotal(&sd, NULL, js, ks, k, 1);
	      }
	      if (i == 0) {
		if (1.0+sd == 1.0) {
		  break;
		}
	      }
	    }
	    pkd[q] = sd;
	    pke[q] = se;
	    q++;
	  }
	  if (i > 0) {
	    kappa0[m] = kpp0;
	    kappa1[m] = kpp1;
	    m++;
	  }
	}
      }
    }
  }

  (*pk)->nkappa = m;
  if (pw_type == 0) {
    (*pk)->kappa0 = ReallocNew(kappa0, sizeof(short)*m);
    (*pk)->kappa1 = ReallocNew(kappa1, sizeof(short)*m);
  } else {
    (*pk)->kappa0 = ReallocNew(kappa1, sizeof(short)*m);
    (*pk)->kappa1 = ReallocNew(kappa0, sizeof(short)*m);
  }
  (*pk)->pkd = ReallocNew(pkd, sizeof(double)*q);
  (*pk)->pke = ReallocNew(pke, sizeof(double)*q);  
  (*pk)->nkl = t;
  
  if (locked) ReleaseLock(lock);
#pragma omp flush
#ifdef PERFORM_STATISTICS
  stop = clock();
  timing.rad_qk += stop-start;
#endif
  return type;
}


static void InterpolateGOS(int n, double *x, double *g, 
			   int ni, double *xi, double *gi) {
  int t;

  UVIP3P(3, n, x, g, ni, xi, gi);
  for (t = 0; t < ni; t++) {
    if (xi[t] < x[0]) {
      gi[t] = g[0];
    } else {
      break;
    }
  }
  for (t = ni-1; t >= 0; t--) {
    if (xi[t] > x[n-1]) {
      gi[t] = 0.0;
    } else {
      break;
    }
  }
}

double BornFormFactorK(double q, FORM_FACTOR *bform) {
  double a, r;

  r = 1.0;
  if (bform->nk == 0) {
    a = 1.0/(1.0 + 0.25*q*q);
    a *= a;
    if (bform->te > 0) {
      r = 1.0 - a*a;
    } else {
      r = (1-a) * (1-a);
    }
  } else if (bform->nk > 0) {
    if (q <= bform->k[0]) r = bform->fk[0];
    else if (q >= bform->k[bform->nk-1]) r = bform->fk[bform->nk-1];
    else {   
      q = log(q);
      UVIP3P(3, bform->nk, bform->logk, bform->fk, 1, &q, &r);
    }
  }

  return r;
}

int CERadialQkBorn(int k0, int k1, int k2, int k3, int k, 
		   double te, double e1, double *qk, int m) {
  int p0, p1, p2, p3;
  int m0, m1, m2, m3;
  int j0, j1, j2, j3;
  int ko2, t, nk, ty, bnk;
  double r, c0, c1, dk;
  double x, b, d, c, a, h, g0, b0, b1, a0, a1;
  double *g1, *g2, *x1, *x2;
  double bte, bms;
  FORM_FACTOR *bform;

  ko2 = k/2;  
  ty = ko2;
  *qk = 0.0;
  p0 = GetOrbital(k0)->kappa;
  GetJLFromKappa(p0, &j0, &m0);
  p1 = GetOrbital(k1)->kappa;
  GetJLFromKappa(p1, &j1, &m1);
  if (!Triangle(j0, k, j1)) {
    return -1;
  }
  if (IsOdd((m0+m1+k)/2)) {
    ty = -1;
  }
  p2 = GetOrbital(k2)->kappa;
  GetJLFromKappa(p2, &j2, &m2);
  p3 = GetOrbital(k3)->kappa;
  GetJLFromKappa(p3, &j3, &m3);
  if (!Triangle(j2, k, j3)) {
    return -1;
  }
  if (IsOdd((m2+m3+k)/2)) {
    ty = -1;
  }
  if (ty < 0) return -1;

  r = ReducedCL(j0, k, j1) * ReducedCL(j2, k, j3);
  r *= (k+1.0)*(k+1.0);
  g1 = GeneralizedMoments(k0, k1, ko2);
  x1 = g1 + NGOSK;
  g2 = GeneralizedMoments(k2, k3, ko2);
  x2 = g2 + NGOSK;

  bnk = BornFormFactorTE(&bte);
  bms = BornMass();
  c0 = e1 + (te + bte)/bms;
  if (m <= 0) {
    a0 = FINE_STRUCTURE_CONST2*c0;
    a1 = FINE_STRUCTURE_CONST2*e1;
    b0 = 1 + a0;
    b1 = 1 + a1;
    a0 = 1.0 + 0.5*a0;
    a1 = 1.0 + 0.5*a1;
    c0 = 2.0*c0*a0;
    c1 = 2.0*e1*a1;
  } else {
    c0 = 2.0*c0;
    c1 = 2.0*e1;
  }
  c0 = bms * sqrt(c0);
  c1 = bms * sqrt(c1);
  nk = NKINT-1;
  kint[0] = c0 - c1;
  kint[nk] = c0 + c1;
  log_kint[0] = log(kint[0]);
  log_kint[nk] = log(kint[nk]);
  x = Min(x1[NGOSK-1], x2[NGOSK-1]);  
  if (x < log_kint[nk] && x > log_kint[0]) {
    log_kint[nk] = x;
    kint[nk] = exp(x);
  }
  dk = (log_kint[nk] - log_kint[0])/nk;
  for (t = 1; t < nk; t++) {
    log_kint[t] = log_kint[t-1] + dk;
    kint[t] = exp(log_kint[t]);
  }

  nk = NKINT;
  InterpolateGOS(NGOSK, x1, g1, nk, log_kint, gos1);
  InterpolateGOS(NGOSK, x2, g2, nk, log_kint, gos2);
  
  for (t = 0; t < nk; t++) {
    gosint[t] = r*gos1[t]*gos2[t];
  }
  if (m <= 0) {
    a = a0*a1;
    d = 0.25*FINE_STRUCTURE_CONST2/a;
    for (t = 0; t < nk; t++) {
      c = 0.5*(c0*c0 + c1*c1 - kint[t]*kint[t]);
      h = c/(c0*c1);
      c /= bms*bms;
      c = 1 + d*c;
      c *= c;
      h = 1.0 - h;
      x = (c0/bms)*(c1/bms);
      h *= x*x;
      h *= d*d;
      gosint[t] *= a*(c + h);
    }
  }
  if (bnk >= 0) {
    bform = BornFormFactor();
    for (t = 0; t < nk; t++) {
      a = BornFormFactorK(kint[t], bform);
      gosint[t] *= a;
    }
  }
  a = dk*Simpson(gosint, 0, nk-1);
  *qk += a;
    
  return ty;
}
  
int CERadialQkBornMSub(int k0, int k1, int k2, int k3, int k, int kp,
		       double te, double e1,
		       int nq, int *q, double *qk, int m) {
  int p0, p1, p2, p3;
  int m0, m1, m2, m3;
  int j0, j1, j2, j3;
  int ko2, ko2p, t, nk;
  int nudiff, mu1, mu2, ierr, ipqa[MAXMSUB];
  int kkp, iq, bnk;
  double xc, theta, dnu1, pqa[MAXMSUB];
  double r, c0, c1, c01, dk, a0, a1;
  double x, b, d, c, a, h, g0, b0, b1, bte, bms;  
  double *g1, *g2, *x1, *x2;
  double gosm1[MAXMSUB][NKINT];
  double gosm2[MAXMSUB][NKINT];
  FORM_FACTOR *bform;
  
  for (iq = 0; iq < nq; iq++) {
    qk[iq] = 0.0;
  }
  p0 = GetOrbital(k0)->kappa;
  GetJLFromKappa(p0, &j0, &m0);
  p1 = GetOrbital(k1)->kappa;
  GetJLFromKappa(p1, &j1, &m1);
  if (IsOdd((m0+m1+k)/2) || !Triangle(j0, k, j1)) {
    return -1;
  }
  p2 = GetOrbital(k2)->kappa;
  GetJLFromKappa(p2, &j2, &m2);
  p3 = GetOrbital(k3)->kappa;
  GetJLFromKappa(p3, &j3, &m3);
  if (IsOdd((m2+m3+kp)/2) || !Triangle(j2, kp, j3)) {
    return -1;
  }
  
  ko2 = k/2;
  ko2p = kp/2;  
  kkp = (ko2 + ko2p)%4;
  if (kkp == 1 || kkp == 3) {
    return Max(ko2, ko2p);
  }

  r = ReducedCL(j0, k, j1) * ReducedCL(j2, kp, j3);
  r *= (k+1.0)*(kp+1.0);
  g1 = GeneralizedMoments(k0, k1, ko2);
  x1 = g1 + NGOSK;
  g2 = GeneralizedMoments(k2, k3, ko2p);
  x2 = g2 + NGOSK;

  bnk = BornFormFactorTE(&bte);
  bms = BornMass();  
  c0 = e1 + (te + bte)/bms; 
  if (m <= 0) {
    a0 = FINE_STRUCTURE_CONST2*c0;
    a1 = FINE_STRUCTURE_CONST2*e1;
    b0 = 1 + a0;
    b1 = 1 + a1;
    a0 = 1.0 + 0.5*a0;
    a1 = 1.0 + 0.5*a1;
    c0 = 2.0*c0*a0;
    c1 = 2.0*e1*a1;
  } else {
    c0 = 2.0*c0;
    c1 = 2.0*e1;
  }
  c01 = (c0 - c1)*bms;
  c0 = bms * sqrt(c0);
  c1 = bms * sqrt(c1);
  nk = NKINT-1;
  kint[0] = c0 - c1;
  kint[nk] = c0 + c1;
  log_kint[0] = log(kint[0]);
  log_kint[nk] = log(kint[nk]);
  x = Min(x1[NGOSK-1], x2[NGOSK-1]);
  if (x < log_kint[nk] && x > log_kint[0]) {
    log_kint[nk] = x;
    kint[nk] = exp(x);
  }
  dk = (log_kint[nk] - log_kint[0])/nk;
  for (t = 1; t < nk; t++) {
    log_kint[t] = log_kint[t-1] + dk;
    kint[t] = exp(log_kint[t]);
  }

  nk = NKINT;
  InterpolateGOS(NGOSK, x1, g1, nk, log_kint, gos1);
  InterpolateGOS(NGOSK, x2, g2, nk, log_kint, gos2);

  for (t = 0; t < nk; t++) {
    gost[t] = r*gos1[t]*gos2[t];
  }
  if (m <= 0) {
    a = a0*a1;
    d = 0.25*FINE_STRUCTURE_CONST2/a;
    for (t = 0; t < nk; t++) {
      c = 0.5*(c0*c0 + c1*c1 - kint[t]*kint[t]);
      h = c/(c0*c1);
      c /= bms*bms;
      c = 1 + d*c;
      c *= c;
      h = 1.0 - h;
      x = (c0/bms)*(c1/bms);
      h *= x*x;
      h *= d*d;
      gost[t] *= a*(c + h);
    }
  }

  if (bnk >= 0) {
    bform = BornFormFactor();
    for (t = 0; t < nk; t++) {
      a = BornFormFactorK(kint[t], bform);
      gost[t] *= a;
    }
  }

  nudiff = 0;
  mu1 = 0;
  mu2 = q[nq-1]/2;
  for (t = 0; t < nk; t++) {
    xc = (c01+kint[t]*kint[t])/(2.0*c0*kint[t]);
    if (xc < 0.0) xc = 0.0;
    if (xc > 1.0) xc = 1.0;
    theta = acos(xc);
    if (theta < EPS10) theta = EPS10;
    dnu1 = ko2;
    DXLEGF(dnu1, nudiff, mu1, mu2, theta, 3, pqa, ipqa, &ierr);
    for (iq = 0; iq < nq; iq++) {
      gosm1[iq][t] = pqa[iq]*	
	  exp(0.5*(LnFactorial(ko2-q[iq]/2)-LnFactorial(ko2+q[iq]/2)));
    }
    if (kp != k) {
      dnu1 = ko2p;
      DXLEGF(dnu1, nudiff, mu1, mu2, theta, 3, pqa, ipqa, &ierr);
      for (iq = 0; iq < nq; iq++) {
	gosm2[iq][t] = pqa[iq]*	
	  exp(0.5*(LnFactorial(ko2p-q[iq]/2)-LnFactorial(ko2p+q[iq]/2)));
      }
    } else {
      for (iq = 0; iq < nq; iq++) {
	gosm2[iq][t] = gosm1[iq][t];
      }
    }
  }

  for (iq = 0; iq < nq; iq++) {
    for (t = 0; t < nk; t++) {
      gosint[t] = gost[t]*gosm1[iq][t]*gosm2[iq][t];
      if (IsOdd(ko2p+kkp/2)) gosint[t] = -gosint[t];
    }
    qk[iq] = dk*Simpson(gosint, 0, nk-1);
  }  
  
  return Max(ko2, ko2p);
}

double *CERadialQkTable(int k0, int k1, int k2, int k3, int k, int trylock) {
  int type, t, ie, ite, ipk, ipkp, nqk, nopb;
  int i, j, kl0, kl1, kl, nkappa, nkl, nkappap, nklp;
  CEPK *cepk, *cepkp, *tmp;
  short *kappa0, *kappa1, *kappa0p, *kappa1p;
  double *pkd, *pke, *pkdp, *pkep, *pd;
  double r, rd, s, b, a, c;
  double qk[MAXNKL], dqk[MAXNKL];
  double rq[MAXNTE][MAXNE+1], e1, te, te0;
  double drq[MAXNTE][MAXNE+1], *rqc, **p, *ptr;
  int index[5], mb, mk;
  int np = 3, one = 1, ieb[MAXNTE];
  double logj, xb;

#ifdef PERFORM_STATISTICS
  clock_t start, stop;
  start = clock();
#endif

  index[0] = k/2;
  index[1] = k0;
  index[2] = k1;
  index[3] = k2;
  index[4] = k3;
  LOCK *lock = NULL;
  int locked = 0;
  p = (double **) MultiSet(qk_array, index, NULL, &lock,
			   InitPointerData, FreeExcitationQkData);
  if (lock && !(*p)) {
    if (trylock) {
      if (0 == TryLock(lock)) {
	locked = 1;
      } else {
	return NULL;
      }
    } else {
      SetLock(lock);
      locked = 1;
    }
  }
  if (*p) {
    if (locked) ReleaseLock(lock);
    return *p;
  }   
  if (xborn == 0 || xborn < -1E30) {
    for (ie = 0; ie < n_egrid1; ie++) {
      if (ie == n_egrid) mb = 1;
      else {
	if (xborn == 0) {
	  mb = 0;
	} else {
	  mb = -1;
	}
      }
      e1 = egrid[ie];
      for (ite = 0; ite < n_tegrid; ite++) {
	te = tegrid[ite];
	type = CERadialQkBorn(k0, k1, k2, k3, k,
			      te, e1, &(rq[ite][ie]), mb);
	drq[ite][ie] = rq[ite][ie];
      }
    }
  } else {
    te0 = -GetOrbital(k0)->energy;
    te = -GetOrbital(k1)->energy;
    te0 = Max(te0, te);
    te = -GetOrbital(k2)->energy;
    te0 = Max(te0, te);
    te = -GetOrbital(k3)->energy;
    te0 = Max(te0, te);
    ie = n_egrid;
    for (ie = n_egrid; ie < n_egrid1; ie++) {
      e1 = egrid[ie];
      for (ite = 0; ite < n_tegrid; ite++) {
	te = tegrid[ite];
	type = CERadialQkBorn(k0, k1, k2, k3, k,
			      te, e1, &(rq[ite][ie]), 1);
	drq[ite][ie] = rq[ite][ie];
      }    
    }
    for (ite = 0; ite < n_tegrid; ite++) {
      ieb[ite] = 0;
    }
    for (ie = 0; ie < n_egrid; ie++) {
      e1 = egrid[ie];
      t = 0;
      for (ite = 0; ite < n_tegrid; ite++) {
	if (ieb[ite] == 2) {
	  te = tegrid[ite];
	  type = CERadialQkBorn(k0, k1, k2, k3, k,
				te, e1, &(rq[ite][ie]), 0);
	  drq[ite][ie] = rq[ite][ie];
	} else {
	  t = 1;
	}
      }
      if (t == 0) continue;
      type = CERadialPk(&cepk, ie, k0, k1, k, 0);
      if (k2 != k0 || k3 != k1) {
	type = CERadialPk(&cepkp, ie, k2, k3, k, 0);
      } else {
	cepkp = cepk;
      }
      if (cepk->nkl > cepkp->nkl) {
	tmp = cepk;
	cepk = cepkp;
	cepkp = tmp;
      }
      kappa0 = cepk->kappa0;
      kappa1 = cepk->kappa1;
      nkappa = cepk->nkappa;
      kappa0p = cepkp->kappa0;
      kappa1p = cepkp->kappa1;
      nkappap= cepkp->nkappa;
      pkd = cepk->pkd;
      pke = cepk->pke;
      pkdp = cepkp->pkd;
      pkep = cepkp->pke;
      nkl = cepk->nkl;
      nklp = nkl-1;
      for (ite = 0; ite < n_tegrid; ite++) {
	if (ieb[ite] == 2) continue;
	te = tegrid[ite];
	for (i = 0; i < nkl; i++) {
	  qk[i] = 0.0;
	  dqk[i] = 0.0;
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
	    s = (pkd[ipk]+pke[ipk])*(pkd[ipk]+pke[ipk]);
	    qk[t] += s;
	    s = pkd[ipk]*pkd[ipk];
	    dqk[t] += s;
	  } else {
	    s = 0.0;
	    ipkp = ite;
	    for (j = 0; j < nkappap; j++) {
	      if (kappa0[i] == kappa0p[j] && kappa1[i] == kappa1p[j]) {
		s = (pkd[ipk]+pke[ipk])*(pkdp[ipkp]+pkep[ipkp]);
		qk[t] += s;
		s = pkd[ipk]*pkdp[ipkp];
		dqk[t] += s;
		break;
	      }
	      ipkp += n_tegrid;
	    }
	  }
	  ipk += n_tegrid;
	}
	if (fpw) {
	  for (i = 0; i < nkl; i++) {
	    kl0 = pw_scratch.kl[i];
	    fprintf(fpw, "%12.5E %12.5E %3d %3d %12.5E %12.5E\n",
		    te*HARTREE_EV, e1*HARTREE_EV, i, kl0, qk[i], dqk[i]);
	  }
	}
	r = qk[0];
	rd = dqk[0];
	for (i = 1; i < nkl; i++) {
	  r += qk[i];
	  rd += dqk[i];
	  kl0 = pw_scratch.kl[i-1];
	  kl1 = pw_scratch.kl[i];
	  for (j = kl0+1; j < kl1; j++) {
	    logj = LnInteger(j);
	    UVIP3P(np, nkl, pw_scratch.log_kl, qk, 
		   one, &logj, &s);
	    r += s;
	    UVIP3P(np, nkl, pw_scratch.log_kl, dqk, 
		   one, &logj, &s);
	    rd += s;
	  }
	}
	nklp = nkl-1;
	if (ieb[ite] == 0) {
	  s = 0.0;	  
	  if (1+dqk[nklp] != 1) {
	    a = dqk[nklp]/rd;
	    if (type == 0 || type > CBMULT) {
	      c = dqk[nklp]/dqk[nklp-1];
	      c = pow(c, 1.0/(pw_scratch.kl[nklp]-pw_scratch.kl[nklp-1]));
	      if (c > 0 && c < 1) {
		b = c/(1.0-c);
	      } else {
		b = GetCoulombBetheAsymptotic(te, e1);	    
		if (b*a > xborn0) b = -1.0;
		else b = 0.0;
	      }
	      xb = xborn;
	    } else if (type > 0) {
	      b = (GetCoulombBethe(0, ite, ie, k/2, 0))[nklp];
	      if (b < 0 || IsNan(b)) {
		c = dqk[nklp]/dqk[nklp-1];
		c = pow(c, 1.0/(pw_scratch.kl[nklp]-pw_scratch.kl[nklp-1]));
		if (c > 0 && c < 1) {
		  b = c/(1.0-c);
		} else {
		  b = GetCoulombBetheAsymptotic(te, e1);	    
		  if (b*a > xborn0) b = -1.0;
		  else b = 0.0;
		}
	      }
	      xb = xborn1;
	    } else {
	      b = 0.0;
	      xb = xborn;
	    }
	  } else { 
	    b = 0.0;
	    xb = xborn;
	  }
	  if (b < 0) {
	    ieb[ite] = 1;
 	  } else {
  	    s = dqk[nklp]*b;	  
	    if (ite == 0 &&
	        ((xb < 0 && rd && -xb < s/rd) ||
	         (xb > 0 && xb < e1/te0))) {	      
	      ieb[ite] = 1;
	    } else {
	      rq[ite][ie] = r + s;
	      drq[ite][ie] = rd + s;
	    }
	  }
	}
	if (ieb[ite]) {
	  type = CERadialQkBorn(k0, k1, k2, k3, k,
				te, e1, &(rq[ite][ie]), 0);
	  drq[ite][ie] = rq[ite][ie];
	  rq[ite][ie] += r-rd;
	  if (rq[ite][ie]) {
	    if (fabs((r-rd)/rq[ite][ie]) < 0.05) ieb[ite] = 2;
	  }
	}
      }
    }
  }  

  nqk = n_tegrid*n_egrid1;
  t = nqk + 1;
  if (type >= 0 && k > 0) {
    mk = GetMaxKMBPT();
    if (k/2 <= mk) t = nqk*2 + 1;
  }
  pd = (double *) malloc(sizeof(double)*t);
  rqc = pd;

  ptr = rqc;
  for (ite = 0; ite < n_tegrid; ite++) {
    for (ie = 0; ie < n_egrid1; ie++) {
      ptr[ie] = rq[ite][ie];
    }
    ptr += n_egrid1;
  }
  rqc[nqk] = type;
  
  if (t > nqk+1) {
    ptr = &(rqc[nqk+1]);
    for (ite = 0; ite < n_tegrid; ite++) {
      for (ie = 0; ie < n_egrid1; ie++) {
	ptr[ie] = drq[ite][ie];
      }
      ptr += n_egrid1;
    }
  }

#ifdef PERFORM_STATISTICS
  stop = clock();
  timing.rad_qk += stop-start;
#endif
  if (p) *p = pd;
  if (locked) ReleaseLock(lock);
#pragma omp flush
  return *p;
}

double *CERadialQkMSubTable(int k0, int k1, int k2, int k3, int k, int kp,
			    int trylock) {
  int type1, type2, kl, nqk;
  int i, j, kl0, klp0, kl0_2, klp0_2, kl1;
  CEPK *cepk, *cepkp;
  int nkappa, nkappap, nkl, nklp;
  short  *kappa0, *kappa1, *kappa0p, *kappa1p;
  double *pkd, *pke, *pkdp, *pkep, *ptr;
  int km0, km1, j0, jp0, j1, kmp0, kmp1, km0_m, kmp0_m;
  int mi, mf, t, c0, cp0;
  double r, rd, e0, e1, te, s, sd, b, te0;
  double pha0, phap0, xb, c, a;
  double s3j1, s3j2, s3j3, s3j4;
  int ie, ite, q[MAXMSUB], nq, iq, ipk, ipkp, ieb[MAXNTE];
  double qk[MAXMSUB][MAXNKL], dqk[MAXMSUB][MAXNKL];
  double rq[MAXMSUB][MAXNTE][MAXNE+2], drq[MAXMSUB][MAXNTE][MAXNE+2];
  double rqt[MAXMSUB];
  double *rqc, **p;
  int index[5], mb;
  int np = 3, one = 1;
  double logj;
#ifdef PERFORM_STATISTICS
  clock_t start, stop;
  start = clock();
#endif

  index[0] = (k/2)*(1+(GetMaxRank()/2)) + kp/2;
  index[1] = k0;
  index[2] = k1;
  index[3] = k2;
  index[4] = k3;
  LOCK *lock = NULL;
  int locked = 0;
  p = (double **) MultiSet(qk_array, index, NULL, &lock,
			   InitPointerData, FreeExcitationQkData);
  if (lock && !(*p)) {
    if (trylock) {
      if (0 == TryLock(lock)) {
	locked = 1;
      } else {
	return NULL;
      }
    } else {
      SetLock(lock);
      locked = 1;
    }
  }
  if (*p) {
    if (locked) ReleaseLock(lock);
    return *p;
  }

  nq = Min(k, kp)/2 + 1;
  q[0] = 0;
  for (iq = 1; iq < nq; iq++) {
    q[iq] = q[iq-1] + 2;
  }  
  nqk = nq*n_tegrid*n_egrid1;
  double *pd = (double *) malloc(sizeof(double)*(nqk+1));
  rqc = pd;
  if (xborn == 0) {
    for (ie = 0; ie < n_egrid1; ie++) {
      e1 = egrid[ie];
      if (ie == n_egrid) mb = 1;
      else {
	if (xborn == 0) {
	  mb = 0;
	} else {
	  mb = -1;
	}
      }
      for (ite = 0; ite < n_tegrid; ite++) {
	te = tegrid[ite];
	type1 = CERadialQkBornMSub(k0, k1, k2, k3, k, kp, te, e1, 
				   nq, q, rqt, mb);	
	for (iq = 0; iq < nq; iq++) {
	  rq[iq][ite][ie] = rqt[iq];
	}
      }
    }
    type2 = type1;
  } else {
    pkdp = NULL;
    pkep = NULL;    
    for (ie = n_egrid; ie < n_egrid1; ie++) {
      e1 = egrid[ie];
      for (ite = 0; ite < n_tegrid; ite++) {
	te = tegrid[ite];
	type1 = CERadialQkBornMSub(k0, k1, k2, k3, k, kp, te, e1, 
				   nq, q, rqt, 1);	
	for (iq = 0; iq < nq; iq++) {
	  rq[iq][ite][ie] = rqt[iq];
	}
      }
    }
    te0 = -GetOrbital(k0)->energy;
    te = -GetOrbital(k1)->energy;
    te0 = Max(te0, te);
    te = -GetOrbital(k2)->energy;
    te0 = Max(te0, te);
    te = -GetOrbital(k3)->energy;
    te0 = Max(te0, te);
    for (ite = 0; ite < n_tegrid; ite++) {
      ieb[ite] = 0;
    }
    for (ie = 0; ie < n_egrid; ie++) {
      e1 = egrid[ie];
      type1 = CERadialPk(&cepk, ie, k0, k1, k, 0);
      nkl = cepk->nkl;
      nkappa = cepk->nkappa;
      kappa0 = cepk->kappa0;
      kappa1 = cepk->kappa1;
      pkd = cepk->pkd;
      pke = cepk->pke;
      if (kp == k && k2 == k0 && k3 == k1) {
	cepkp = cepk;
	type2 = type1;
      } else {
	type2 = CERadialPk(&cepkp, ie, k2, k3, kp, 0);
      }
      nklp = cepkp->nkl;
      if (nklp < nkl) nkl = nklp;
      nkappap = cepkp->nkappa;
      kappa0p = cepkp->kappa0;
      kappa1p = cepkp->kappa1;
      pkdp = cepkp->pkd;
      pkep = cepkp->pke;

      for (ite = 0; ite < n_tegrid; ite++) {
	te = tegrid[ite];
	e0 = e1 + te;
      
	for (i = 0; i < nq; i++) { 
	  for (j = 0; j < nkl; j++) { 
	    qk[i][j] = 0.0; 
	    dqk[i][j] = 0.0;
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
	    GetJLFromKappa(kmp0, &jp0, &klp0);
	    if (kmp1 != km1) {
	      ipkp += n_tegrid;
	      continue;
	    }
	    GetJLFromKappa(kmp0, &jp0, &klp0);
	    klp0_2 = klp0/2;
	  
	    s = (pkd[ipk]+pke[ipk])*(pkdp[ipkp]+pkep[ipkp]);
	    sd = pkd[ipk]*pkdp[ipkp];
	    b = sqrt((j0+1.0)*(jp0+1.0)*(kl0+1.0)*(klp0+1.0));
	    s *= b;
	    sd *= b;
	    if (km0 != kmp0) { 
	      km0_m = km0; 
	      kmp0_m = kmp0; 
	      if (kl0_2 >= pw_scratch.qr) { 
		if (km0 > 0) km0_m = -km0 - 1; 
	      } 
	    
	      if (klp0_2 >= pw_scratch.qr) { 
		if (kmp0 > 0) kmp0_m = -kmp0 - 1; 
	      } 
	      c0 = OrbitalIndex(0, km0_m, e0); 
	      cp0 = OrbitalIndex(0, kmp0_m, e0);
	      pha0 = GetPhaseShift(c0); 
	      phap0 = GetPhaseShift(cp0);	      
	      r = cos(pha0 - phap0);
	      s *= r;
	      sd *= r;
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
	      dqk[iq][i] += sd*rqt[iq];
	    } 

	    ipkp += n_tegrid;
	  }
	  ipk += n_tegrid;
	}
      
	for (iq = 0; iq < nq; iq++) { 
	  r = qk[iq][0];
	  rd = qk[iq][0];
	  for (i = 1; i < nkl; i++) { 
	    r += qk[iq][i]; 
	    rd += qk[iq][i];
	    kl0 = pw_scratch.kl[i-1]; 
	    kl1 = pw_scratch.kl[i]; 
	    for (j = kl0+1; j < kl1; j++) {       
	      logj = LnInteger(j);
	      UVIP3P(np, nkl, pw_scratch.log_kl, qk[iq],
		     one, &logj, &s);
	      r += s;
	      UVIP3P(np, nkl, pw_scratch.log_kl, dqk[iq],
		     one, &logj, &s);
	      rd += s;
	    }      
	  }    
	  rq[iq][ite][ie] = r;
	  drq[iq][ite][ie] = rd;
	} 

	if (ieb[ite] == 0) {
	  i = nkl - 1;
	  r = 0.0;
	  if (type1 > 0 && type1 <= CBMULT) xb = xborn1;
	  else xb = xborn;
	  for (iq = 0; iq < nq; iq++) {
	    if (1+dqk[iq][i] != 1) {
	      a = dqk[iq][i]/rq[iq][ite][ie];
	      if (k != kp || type1 == 0 || type1 > CBMULT) {
		c = dqk[iq][i]/dqk[iq][i-1];
		c = pow(c, 1.0/(pw_scratch.kl[i]-pw_scratch.kl[i-1]));
		if (c > 0 && c < 1) {
		  b = c/(1.0-c);
		} else {
		  b = GetCoulombBetheAsymptotic(te, e1);
		  if (b*a > xborn0) b = -1.0;
		  else b = 0.0;
		}
	      } else if (type1 >= 0) {
		b = (GetCoulombBethe(0, ite, ie, k/2, abs(q[iq])/2))[i];
		if (b < 0 || IsNan(b)) {
		  c = dqk[iq][i]/dqk[iq][i-1];
		  c = pow(c, 1.0/(pw_scratch.kl[i]-pw_scratch.kl[i-1]));
		  if (c > 0 && c < 1) {
		    b = c/(1.0-c);
		  } else {
		    b = GetCoulombBetheAsymptotic(te, e1);
		    if (b*a > xborn0) b = -1.0;
		    else b = 0.0;
		  }
		}
	      } else {	  
		b = 0.0;
	      }
	    } else {
	      b = 0.0;
	    }
	    if (b < 0) {
	      ieb[ite] = 1;
	      break;
	    } else {
	      s = dqk[iq][i]*b;
	      rqt[iq] = s;
	      if (xb < 0 && drq[iq][ite][ie]) {
	        s /= drq[iq][ite][ie];
	        if (s > r) r = s;
	      }
	    }
	  }
	  if (ieb[ite] == 0) {
	    if ((xb < 0 && -xb < r) ||
		(xb > 0 && xb < e1/te0)) {
	      ieb[ite] = 1;
	    } else {
	      for (iq = 0; iq < nq; iq++) {
		rq[iq][ite][ie] += rqt[iq];
	      }
	    }
	  }
        } 
	if (ieb[ite]) {
	  type1 = CERadialQkBornMSub(k0, k1, k2, k3, k, kp, te, e1, 
				     nq, q, rqt, 0);
	  for (iq = 0; iq < nq; iq++) {
	    rq[iq][ite][ie] += rqt[iq] - drq[iq][ite][ie];
	  }
	}
      }
    }
  }

  ptr = rqc;
  for (iq = 0; iq < nq; iq++) {
    for (ite = 0; ite < n_tegrid; ite++) {
      for (ie = 0; ie < n_egrid1; ie++) {
	ptr[ie] = rq[iq][ite][ie];
      }
      ptr += n_egrid1;
    }
  }    
  rqc[nqk] = type1;
  if (type2 != 1) rqc[nqk] = type2;

#ifdef PERFORM_STATISTICS
  stop = clock();
  timing.rad_qk += stop-start;
#endif
  *p = pd;
  if (locked) ReleaseLock(lock);
#pragma omp flush
  return rqc;
} 	
  
int CERadialQk(double *rqc, double te, int k0, int k1, int k2, int k3,
	       int k, int trylock) {
  int i, np, nd, type;
  int j, m, mk;
  double *rqe, rq[MAXNTE];
  double *xte, x0;
  
  rqe = CERadialQkTable(k0, k1, k2, k3, k, trylock);
  if (rqe == NULL) return -9999;
  if (n_tegrid == 1) {
    for (i = 0; i < n_egrid1; i++) {
      rqc[i] = rqe[i];
    }
    type = rqe[i];
  } else {
    np = 3;
    nd = 1;
    type = rqe[n_tegrid*n_egrid1];
    if (type == 0 || type == 1) {
      xte = log_te;
      x0 = log(te);
    } else {
      xte = tegrid;
      x0 = te;
    }
    for (i = 0; i < n_egrid1; i++) {
      j = i;
      for (m = 0; m < n_tegrid; m++) {
	rq[m] = rqe[j];
	j += n_egrid1;
      }
      UVIP3P(np, n_tegrid, xte, rq, nd, &x0, &rqc[i]);
    }
  }

  if (type >= 0 && k > 0) {
    mk = GetMaxKMBPT();
    if (k/2 <= mk) {
      rqe += n_tegrid*n_egrid1+1;
      rqc += n_egrid1;
      if (n_tegrid == 1) {
	for (i = 0; i < n_egrid1; i++) {
	  rqc[i] = rqe[i];
	}
      } else {
	np = 3;
	nd = 1;
	if (type == 0 || type == 1) {
	  xte = log_te;
	  x0 = log(te);
	} else {
	  xte = tegrid;
	  x0 = te;
	}
	for (i = 0; i < n_egrid1; i++) {
	  j = i;
	  for (m = 0; m < n_tegrid; m++) {
	    rq[m] = rqe[j];
	    j += n_egrid1;
	  }
	  UVIP3P(np, n_tegrid, xte, rq, nd, &x0, &rqc[i]);
	}
      }
    }
  }
 
  return type;
}

int CERadialQkMSub(double *rqc, double te, int k0, int k1, int k2, int k3, 
		   int k, int kp, int trylock) {
  int i, np, nd, iq, n;
  int j, m, type, nq;
  double *rqe, rq[MAXNTE];
  double *xte, x0;
  
  rqe = CERadialQkMSubTable(k0, k1, k2, k3, k, kp, trylock);
  if (rqe == NULL) return -9999;
  nq = Min(k, kp)/2 + 1;

  if (n_tegrid == 1) {
    for (iq = 0; iq < nq; iq++) {
      for (i = 0; i < n_egrid1; i++) {
	rqc[i] = rqe[i];
      }
      rqc += n_egrid1;
      rqe += n_egrid1;
    }
    type = rqe[0];
  } else {
    np = 3;
    nd = 1;
    n = n_tegrid*n_egrid1;
    type = rqe[nq*n];
    if (type == 0 || type == 1) {
      xte = log_te;
      x0 = log(te);
    } else {
      xte = tegrid;
      x0 = te;
    }
    for (iq = 0; iq < nq; iq++) {
      for (i = 0; i < n_egrid1; i++) {
	j = i;
	for (m = 0; m < n_tegrid; m++) {
	  rq[m] = rqe[j];
	  j += n_egrid1;
	}
	UVIP3P(np, n_tegrid, xte, rq, nd, &x0, &rqc[i]);
      }
      rqe += n;
      rqc += n_egrid1;
    }
  }  
  return type;
}

void BornFromFit(int np, double *p, int n, double *x, double *logx,
		 double *y, double *dy, int ndy, void *extra) {
  double a, b;

  int i, k;

  if (ndy <= 0) {
    for (i = 0; i < n; i++) {
      a = (1.0 + p[1])/(x[i] + p[1]);
      y[i] = p[0] + p[2]*a + p[3]*a*a;
    }
  } else {
    for (i = 0; i < n; i++) {
      a = x[i] + p[1];
      b = (x[i] - 1.0)/(a*a);
      a = (1.0 + p[1])/a;
      k = i;
      dy[k] = 1.0;
      k += ndy;
      dy[k] = p[2]*b + 2.0*p[3]*b*a;
      k += ndy;
      dy[k] = a;
      k += ndy;
      dy[k] = a*a;
    }
  }
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
 
void RelativisticCorrection(int m, double *s, double *p, double te, double b) {
  int i, j, k;
  double a, c, b1, b0;

  if (b <= 0.0) return;
  for (j = 0; j < n_usr; j++) {
    a = usr_egrid[j];
    c = FINE_STRUCTURE_CONST2*a;
    b1 = 1.0 + c;
    a = usr_egrid[j] + te;
    c = FINE_STRUCTURE_CONST2*a;
    b0 = 1.0 + c;
    a = 2.0*a*(1.0 + 0.5*c)*FINE_STRUCTURE_CONST2;
    c = a/(1.0+a);
    c = -b0*b1*b*(log(1.0-c) + c);
    if (m <= 0) {
      s[j] += c;
    } else {
      for (i = 0; i < m; i++) {
	k = i*n_usr + j;
	s[k] += c*s[k]/p[j];
      }
    }
  }
}

int CollisionStrengthUTA(double *qkt, double *params, double *e, double *bethe,
			 int lower, int upper) {  
  INTERACT_DATUM *idatum;
  LEVEL *lev1, *lev2;
  int p1, p2, j1, j2, k0, k1, type, ty;
  int ns, q1, q2, ie, kmin, kmax, k, np;
  double te, *rqk, tol;
  double rq[MAXMSUB*(MAXNE+1)], qkc[MAXMSUB*(MAXNE+1)];
  int ierr, ipvt[NPARAMS];
  int lwa=5*NPARAMS+MAXNE;
  double wa[5*NPARAMS+MAXNE];
  double fvec[MAXNE], fjac[MAXNE*NPARAMS];
  double born_egrid, born_cross, c, d, r;
  double bte, bms;
  FORM_FACTOR *bform;
  
  lev1 = GetLevel(lower);
  if (lev1 == NULL) return -1;
  lev2 = GetLevel(upper);
  if (lev2 == NULL) return -1;
  te = lev2->energy - lev1->energy;
  if (te <= 0) return -1;
  *e = te;
  p1 = lev1->pj;
  p2 = lev2->pj;

  rqk = qkc;
  for (ie = 0; ie < n_egrid1; ie++) {
    rqk[ie] = 0.0;
  }

  idatum = NULL;
  ns = GetInteract(&idatum, NULL, NULL, lev1->iham, lev2->iham,
		   lev1->pb, lev2->pb, 0, 0, 0);
  if (ns <= 0) return -1;
  if (idatum->s[0].index < 0 || idatum->s[3].index >= 0) {
    free(idatum->bra);
    free(idatum);
    return -1;
  }
  if (idatum->s[0].nq_bra > idatum->s[0].nq_ket) {
    j1 = idatum->s[0].j;
    j2 = idatum->s[1].j;
    q1 = idatum->s[0].nq_bra;
    q2 = idatum->s[1].nq_bra;
    k0 = OrbitalIndex(idatum->s[0].n, idatum->s[0].kappa, 0.0);
    k1 = OrbitalIndex(idatum->s[1].n, idatum->s[1].kappa, 0.0);
  } else {
    j1 = idatum->s[1].j;
    j2 = idatum->s[0].j;
    q1 = idatum->s[1].nq_bra;
    q2 = idatum->s[0].nq_bra;
    k1 = OrbitalIndex(idatum->s[0].n, idatum->s[0].kappa, 0.0);
    k0 = OrbitalIndex(idatum->s[1].n, idatum->s[1].kappa, 0.0);
  }
    
  type = -1;
  kmin = abs(j1-j2);
  kmax = j1 + j2;
  for (k = kmin; k <= kmax; k += 2) {
    ty = CERadialQk(rq, te, k0, k1, k0, k1, k, 0);
    if (ty > type) type = ty;
    for (ie = 0; ie < n_egrid1; ie++) {
      qkc[ie] += rq[ie]/(k+1.0);
    }
  }

  d = (lev1->ilev+1.0)*q1*(j2+1.0-q2)/((j1+1.0)*(j2+1.0));
  for (ie = 0; ie < n_egrid1; ie++) {
    qkc[ie] *= d;
  }
  bte = te;
  if (type >= 0) {
    r = 0.0;
    if (Triangle(j1, j2, 2) && IsOdd(p1+p2)) {
      r = MultipoleRadialNR(-1, k0, k1, G_BABUSHKIN);
    }
    if (fabs(r) > 0.0) {
      r = OscillatorStrength(-1, te, r, NULL);
      bethe[0] = d*2.0*r/te;
      BornFormFactorTE(&bte);
      bms = BornMass();
      bte = (te + bte)/bms;
    } else {
      bethe[0] = 0.0;
    }
    ie = n_egrid;
    born_cross = qkc[ie]*8.0;
    if (born_cross > 0) {
      c = egrid[ie]+bte;
      born_egrid = c/te;
      if (bethe[0] > 0) bethe[1] = born_cross - bethe[0]*log(born_egrid);
      else bethe[1] = born_cross;
      bethe[2] = egrid[ie]; 
    } else {
      bethe[0] = -1.0;
      bethe[1] = 0.0;
      bethe[2] = 0.0;
    }
  } else {
    bethe[0] = -1.0;
    bethe[1] = 0.0;
    bethe[2] = 0.0;
  }
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
      if (bethe[0] < 0.0) {
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
      UVIP3P(np, n_egrid, log_egrid, qkc, n_usr, log_usr, qkt);
      for (ie = 0; ie < n_usr; ie++) {
	qkt[ie] = exp(qkt[ie]);
      }
    } else {
      for (ie = 0; ie < n_egrid; ie++) {
	qkc[ie] *= 8.0;
      }
      UVIP3P(np, n_egrid, log_egrid, qkc, n_usr, log_usr, qkt);
    }
  } else if (qk_mode == QK_EXACT) {
    for (ie = 0; ie < n_usr; ie++) {
      qkt[ie] = 8.0*qkc[ie];
    }
  }
  RelativisticCorrection(0, qkt, params, bte, bethe[0]);
  free(idatum->bra);
  free(idatum);
  return 1;
}

int CollisionStrengthEB(double *qkt, double *e, double *bethe,
			int lower, int upper) {
  LEVEL *lev1, *lev2, *plev1, *plev1p, *plev2, *plev2p;
  double te, a, ap, c, cp, r, s[3];
  double rq[MAXNE+1], qkc[MAXNE+1];
  double born_egrid, born_cross, *mbk;
  int ie, i1, i2, i1p, i2p, p1, p2, p1p, p2p;
  int j1, j2, j1p, j2p, mlev1, mlev2, mlev1p, mlev2p;
  int ilev1, ilev2, ilev1p, ilev2p, i, ip, nz, nzp, k, nmk;
  ANGULAR_ZMIX *ang, *angp;
  double bte, bms;
  FORM_FACTOR *bform;
        
  lev1 = GetEBLevel(lower);
  if (lev1 == NULL) return -1;
  lev2 = GetEBLevel(upper);
  if (lev2 == NULL) return -1;
  te = lev2->energy - lev1->energy;
  if (te <= 0) return -1;
  *e = te;
  
  for (ie = 0; ie < n_egrid1; ie++) {
    qkc[ie] = 0.0;
  }

  for (i1 = 0; i1 < lev1->n_basis; i1++) {
    DecodeBasisEB(lev1->basis[i1], &ilev1, &mlev1);
    plev1 = GetLevel(ilev1);
    DecodePJ(plev1->pj, &p1, &j1);
    for (i2 = 0; i2 < lev2->n_basis; i2++) {
      DecodeBasisEB(lev2->basis[i2], &ilev2, &mlev2);
      plev2 = GetLevel(ilev2);
      DecodePJ(plev2->pj, &p2, &j2);
      c = lev1->mixing[i1]*lev2->mixing[i2];      
      if (fabs(c) < EPS10) continue;
      nz = AngularZMix(&ang, ilev1, ilev2, -1, -1, &nmk, &mbk);
      for (i = 0; i < nz; i++) {
	a = W3j(j1, ang[i].k, j2, -mlev1, mlev1-mlev2, mlev2);
	if (IsOdd((j1-mlev1)/2)) a = -a;
	a *= c*ang[i].coeff;
	if (fabs(a) < EPS10) continue;
	for (i1p = 0; i1p < lev1->n_basis; i1p++) {
	  DecodeBasisEB(lev1->basis[i1p], &ilev1p, &mlev1p);
	  plev1p = GetLevel(ilev1p);
	  DecodePJ(plev1p->pj, &p1p, &j1p);
	  for (i2p = 0; i2p < lev2->n_basis; i2p++) {
	    DecodeBasisEB(lev2->basis[i2p], &ilev2p, &mlev2p);
	    if (mlev1p-mlev2p != mlev1-mlev2) continue;
	    plev2p = GetLevel(ilev2p);
	    DecodePJ(plev2p->pj, &p2p, &j2p);
	    cp = lev1->mixing[i1p]*lev2->mixing[i2p];
	    if (fabs(cp) < EPS10) continue;
	    nzp = AngularZMix(&angp, ilev1p, ilev2p, -1, -1, &nmk, &mbk);
	    for (ip = 0; ip < nzp; ip++) {
	      if (angp[ip].k != ang[i].k) continue;
	      ap = W3j(j1p, angp[ip].k, j2p, -mlev1p, mlev1p-mlev2p, mlev2p);
	      if (IsOdd((j1p-mlev1p)/2)) ap = -ap;
	      ap *= cp*angp[ip].coeff;
	      if (fabs(ap) < EPS10) continue;
	      r = a*ap/(ang[i].k + 1.0);
	      k = CERadialQk(rq, te, ang[i].k0, ang[i].k1, 
			     angp[ip].k0, angp[ip].k1, ang[i].k, 0);
	      for (ie = 0; ie < n_egrid1; ie++) {
		qkc[ie] += r*rq[ie];
		/*
		printf("%d %d %d %d %d %d %d %d %d %10.3E %10.3E\n",
		       lower, upper, i1, i2, i, i1p, i2p, ip, ie, r, rq[ie]);
		*/
	      }
	    }
	    if (nzp) free(angp);
	  }
	}	
      }
      if (nz) free(ang);
    }
  }

  BornFormFactorTE(&bte);
  bms = BornMass();
  bte = (te + bte)/bms;
  SetTransitionMode(M_NR);
  SetTransitionGauge(G_BABUSHKIN);
  k = TRMultipoleEB(s, &te, -1, lower, upper);
  if (k != 0) bethe[0] = 0;
  else {
    r = 0.0;
    for (k = 0; k < 3; k++) {
      r += OscillatorStrength(-1, te, s[k], NULL);
    }
    bethe[0] = 2.0*r/te;
  }
  ie = n_egrid;
  born_cross = qkc[ie]*8.0;
  if (born_cross > 0) {
    c = egrid[ie] + bte;
    born_egrid = c/te;
    if (bethe[0] > 0) bethe[1] = born_cross-bethe[0]*log(born_egrid);
    else bethe[1] = born_cross;
    bethe[2] = egrid[ie];
  } else {
    bethe[0] = -1.0;
    bethe[1] = 0.0;
    bethe[2] = 0.0;
  }
  
  for (ie = 0; ie < n_egrid; ie++) {
    qkt[ie] = 8.0*qkc[ie];
  }

  RelativisticCorrection(0, qkt, NULL, bte, bethe[0]);

  return 1;
}

int CollisionStrengthEBD(double *qkt, double *e, double *bethe, double *born,
			 int lower, int upper) {
  LEVEL *lev1, *lev2, *plev1, *plev1p, *plev2, *plev2p;
  double te, a, ap, c, cp, r, s[3];
  double rq[(MAXNE+2)*MAXMSUB];
  double d, d1, d2, rs;
  int q, nq, kkp, qb, qbp, ith, iph, m, ka, nmk;
  double born_egrid, born_cross, *mbk;
  int ie, i1, i2, i1p, i2p, p1, p2, p1p, p2p;
  int j1, j2, j1p, j2p, mlev1, mlev2, mlev1p, mlev2p;
  int ilev1, ilev2, ilev1p, ilev2p, i, ip, nz, nzp, k;
  ANGULAR_ZMIX *ang, *angp;
  double bte, bms;
  FORM_FACTOR *bform;
      
  lev1 = GetEBLevel(lower);
  if (lev1 == NULL) return -1;
  lev2 = GetEBLevel(upper);
  if (lev2 == NULL) return -1;
  te = lev2->energy - lev1->energy;
  if (te <= 0) return -1;
  *e = te;
  
  m = n_egrid1*n_thetagrid*n_phigrid;
  for (ie = 0; ie < m; ie++) {
    qkt[ie] = 0.0;
  }

  for (i1 = 0; i1 < lev1->n_basis; i1++) {
    DecodeBasisEB(lev1->basis[i1], &ilev1, &mlev1);
    plev1 = GetLevel(ilev1);
    DecodePJ(plev1->pj, &p1, &j1);
    for (i2 = 0; i2 < lev2->n_basis; i2++) {
      DecodeBasisEB(lev2->basis[i2], &ilev2, &mlev2);
      plev2 = GetLevel(ilev2);
      DecodePJ(plev2->pj, &p2, &j2);
      c = lev1->mixing[i1]*lev2->mixing[i2];      
      if (fabs(c) < EPS10) continue;
      nz = AngularZMix(&ang, ilev1, ilev2, -1, -1, &nmk, &mbk);
      for (i = 0; i < nz; i++) {
	a = W3j(j1, ang[i].k, j2, -mlev1, mlev1-mlev2, mlev2);
	if (IsOdd((j1-mlev1)/2)) a = -a;
	a *= c*ang[i].coeff;
	if (fabs(a) < EPS10) continue;
	for (i1p = 0; i1p < lev1->n_basis; i1p++) {
	  DecodeBasisEB(lev1->basis[i1p], &ilev1p, &mlev1p);
	  plev1p = GetLevel(ilev1p);
	  DecodePJ(plev1p->pj, &p1p, &j1p);
	  for (i2p = 0; i2p < lev2->n_basis; i2p++) {
	    DecodeBasisEB(lev2->basis[i2p], &ilev2p, &mlev2p);
	    if (mlev1p-mlev2p != mlev1-mlev2) continue;
	    plev2p = GetLevel(ilev2p);
	    DecodePJ(plev2p->pj, &p2p, &j2p);
	    cp = lev1->mixing[i1p]*lev2->mixing[i2p];
	    if (fabs(cp) < EPS10) continue;
	    nzp = AngularZMix(&angp, ilev1p, ilev2p, -1, -1, &nmk, &mbk);
	    for (ip = 0; ip < nzp; ip++) {	      
	      ap = W3j(j1p, angp[ip].k, j2p, -mlev1p, mlev1p-mlev2p, mlev2p);
	      if (IsOdd((j1p-mlev1p)/2)) ap = -ap;
	      ap *= cp*angp[ip].coeff;
	      if (fabs(ap) < EPS10) continue;	      
	      r = a*ap;
	      nq = Min(ang[i].k, angp[ip].k)/2;
	      kkp = (ang[i].k + angp[ip].k)/2;
	      qb = mlev1 - mlev2;
	      qbp = mlev1p - mlev2p;
	      k = CERadialQkMSub(rq, te, ang[i].k0, ang[i].k1, 
				 angp[ip].k0, angp[ip].k1,
				 ang[i].k, angp[ip].k, 0);
	      for (q = -nq; q <= nq; q++) {
		m = abs(q);	
		ka = 0;
		if (q < 0 && IsOdd(kkp)) rs = -1.0;
		else rs = 1.0;
		for (ith = 0; ith < n_thetagrid; ith++) {
		  d1 = WignerDMatrix(thetagrid[ith], ang[i].k, qb, 2*q);
		  d2 = WignerDMatrix(thetagrid[ith], angp[ip].k, qbp, 2*q);
		  for (iph = 0; iph < n_phigrid; iph++) {
		    d = cos(0.5*(qb-qbp)*phigrid[iph]);		
		    d *= d1*d2;
		    for (ie = 0; ie < n_egrid1; ie++) {
		      qkt[ie + ka*n_egrid1] += r*rs*rq[ie + m*n_egrid1]*d;
		    }
		    ka++;
		  }
		}
	      }
	    }
	    if (nzp) free(angp);
	  }
	}	
      }
      if (nz) free(ang);
    }
  }

  BornFormFactorTE(&bte);
  bms = BornMass();
  bte = (te + bte)/bms;
  m = n_egrid1*n_thetagrid*n_phigrid;
  for (ie = 0; ie < m; ie++) {
    qkt[ie] *= 8.0;
  }
  m = n_thetagrid*n_phigrid;
  ie = n_egrid1-1;
  for (i = 0; i < m; i++) {
    born_cross = qkt[ie + i*n_egrid1];
    d = qkt[ie-1 + i*n_egrid1];
    c = egrid[ie] + bte;
    born_egrid = c/te;
    d1 = log(born_egrid);
    c = egrid[ie-1] + bte;
    d2 = log(c/te);
    if (born_cross > d) {
      bethe[i] = (born_cross-d)/(d1-d2);
      born[i] = born_cross - bethe[i]*d1;
    } else {
      bethe[i] = 0.0;
      born[i] = born_cross;
    }
  }
  born[m] = egrid[ie];

  RelativisticCorrection(0, qkt, NULL, bte, bethe[0]);
  return 1;
}

double AngZCorrection(int nmk, double *mbk, ANGULAR_ZMIX *ang, int t) {
  double r;

  if (nmk < t) return 0.0;
  if (t <= 0) return 0.0;
  r = mbk[t-1];
  if (1.0 + r == 1.0) return 0.0;
  r /= MultipoleRadialNR(-t, ang->k0, ang->k1, G_BABUSHKIN);
  if (ang->coeff < 0) r = -r;
  r /= mbk[nmk+t-1];

  return r;
}

int CollisionStrength(double *qkt, double *params, double *e, double *bethe,
		      int lower, int upper, int msub) {
  int i, j, t, h, p, m, type, ty, p1, p2, gauge;  
  LEVEL *lev1, *lev2;
  double te, c, r, s3j, c1, c2, *mbk, aw;
  int j1, j2, ie, nk, np, nq, kkp, nmk;
  double rq[MAXMSUB*(MAXNE+1)];
  double qkc[MAXMSUB*(MAXNE+1)];
  double *rqk, *rqkt, tol;
  int ierr, ipvt[NPARAMS];
  int lwa=5*NPARAMS+MAXNE;
  double wa[5*NPARAMS+MAXNE];
  double fvec[MAXNE], fjac[MAXNE*NPARAMS];
  double born_egrid, born_cross, bt, ubt[MAXNUSR];
  double bte, bms;
  FORM_FACTOR *bform;
  int nz;
  ANGULAR_ZMIX *ang;

  /*
  nz = cecache.nz[ic];
  ang = cecache.az[ic];
  nmk = cecache.nmk[ic];
  mbk = cecache.mbk[ic];
  if (nz <= 0) {
    if (nmk > 0) {
      free(mbk);
      cecache.nmk[ic] = 0;
      cecache.mbk[ic] = NULL;
    }
    cecache.az[ic] = NULL;
    return -1;
  }
  */
  lev1 = GetLevel(lower);
  if (lev1 == NULL) return -1;
  lev2 = GetLevel(upper);
  if (lev2 == NULL) return -1;
  te = lev2->energy - lev1->energy;
  if (te <= 0) return -1;
  *e = te;
  aw = te*FINE_STRUCTURE_CONST;

  j1 = lev1->pj;
  j2 = lev2->pj;
  DecodePJ(j1, &p1, &j1);
  DecodePJ(j2, &p2, &j2);

  if (msub) {  
    j = 0;
    rqk = qkc;
    for (t = -j1; t <= 0; t += 2) {
      for (h = -j2; h <= j2; h += 2) {
	for (ie = 0; ie < n_egrid1; ie++) {
	  rqk[ie] = 0.0;
	}
	rqk += n_egrid1;
      }
    }
  } else {
    for (ie = 0; ie < n_egrid1; ie++) {
      qkc[ie] = 0.0;
      qkc[ie+n_egrid1] = 0.0;
    }    
  }
  gauge = GetTransitionGauge();
  nz = AngularZMix(&ang, lower, upper, -1, -1, &nmk, &mbk);
  if (nz <= 0) {
    if (nmk > 0) free(mbk);
    return -1;
  }

  type = -1;
  short **ia;
  ia = malloc(sizeof(short *)*nz);
  for (i = 0; i < nz; i++) {
    ia[i] = malloc(sizeof(short)*nz);
    for (j = i; j < nz; j++) {
      ia[i][j] = 0;
    }
  }
  //int iter = 0;
  while (1) {
    //iter++;
    int nleft = 0;
    for (i = 0; i < nz; i++) {
      for (j = i; j < nz; j++) {
	if (ia[i][j]) continue;
	c = ang[i].coeff * ang[j].coeff;
	if (fabs(c) < EPS30) {
	  ia[i][j] = 1;
	  continue;
	}
	if (i != j) c *= 2.0;
	if (!msub) {
	  if (ang[i].k != ang[j].k) {
	    ia[i][j] = 1;
	    continue;
	  }
	  c /= ang[i].k + 1.0;
	  if (fpw) {
	    fprintf(fpw, "# %3d %3d %12.5E %12.5E %2d %2d %2d\n", 
		    lower, upper, te*HARTREE_EV, 8.0*c,
		    ang[i].k, n_tegrid, n_egrid);
	  }
	  ty = CERadialQk(rq, te, ang[i].k0, ang[i].k1,
			  ang[j].k0, ang[j].k1, ang[i].k, 1);
	  if (ty == -9999) {
	    nleft++;
	    continue;
	  }
	  ia[i][j] = 1;
	  t = ang[i].k/2;
	  if (ty >= 0 && t > 0 && t <= nmk) {
	    c1 = AngZCorrection(nmk, mbk, ang+i, t);
	    if (i == j) c2 = c1;
	    else {
	      c2 = AngZCorrection(nmk, mbk, ang+j, t);
	    }
	    if (gauge == G_COULOMB) {
	      c1 /= aw;
	      c2 /= aw;
	    }
	    c1 = c*(1.0+c1)*(1.0+c2);
	  } else {
	    c1 = c;
	  }
	  if (fpw) {
	    fprintf(fpw, "\n\n");
	  }
	  if (ty > type) type = ty;	  
	  if (ty >= 0 && t > 0 && t <= nmk) {
	    for (ie = 0; ie < n_egrid1; ie++) {
	      qkc[ie] += c1*(rq[ie+n_egrid1]) + c*(rq[ie]-rq[ie+n_egrid1]);
	      qkc[ie+n_egrid1] += c*rq[ie];
	    }
	  } else {
	    for (ie = 0; ie < n_egrid1; ie++) {
	      qkc[ie] += c*rq[ie];
	      qkc[ie+n_egrid1] += c*rq[ie];
	    }
	  }
	} else {
	  ty = CERadialQkMSub(rq, te, ang[i].k0, ang[i].k1,
			      ang[j].k0, ang[j].k1, ang[i].k, ang[j].k, 1);
	  if (ty == -9999) {
	    nleft++;
	    continue;
	  }
	  nq = Min(ang[i].k, ang[j].k);
	  kkp = (ang[i].k + ang[j].k)/2;
	  if (ty > type) type = ty;
	  rqk = qkc;
	  for (t = -j1; t <= 0; t += 2) {
	    for (h = -j2; h <= j2; h += 2) {
	      m = t-h;
	      if (abs(m) <= nq) {
		s3j = W3j(j1, ang[i].k, j2, -t, m, h);
		if (ang[j].k != ang[i].k) {
		  s3j *= W3j(j1, ang[j].k, j2, -t, m, h);
		} else {
		  s3j *= s3j;
		}
		if (m < 0 && IsOdd(kkp)) s3j = -s3j;
		m = (abs(m)*n_egrid1)/2;
		for (ie = 0; ie < n_egrid1; ie++) {
		  rqk[ie] += c*rq[m+ie]*s3j;
		}
	      }
	      rqk += n_egrid1;
	    }
	  }
	}
      }
    }
    if (nleft == 0) break;
  }
  for (i = 0; i < nz; i++) {
    free(ia[i]);
  }
  free(ia);
  BornFormFactorTE(&bte);
  bms = BornMass();
  bte = (te + bte)/bms;

  if (msub) {
    for (t = 0; t < MAXMSUB; t++) {
      params[t] = 0.0;
    }
  }

  if (type >= 0) {
    t = 0;
    if (!msub) {
      for (ie = 0; ie < n_egrid1; ie++) {
	if (qkc[ie] <= 0 || fabs(qkc[ie]/qkc[ie+n_egrid1]-1.0) >= 0.75) {
	  t = 1;
	  break;
	}
      }
      if (t) {
	for (ie = 0; ie < n_egrid1; ie++) {
	  qkc[ie] = qkc[ie+n_egrid1];
	}
      }
    }
    r = 0.0;
    if (Triangle(j1, j2, 2) && IsOdd(p1+p2)) {
      for (i = 0; i < nz; i++) {
	if (ang[i].k != 2) continue;
	c = MultipoleRadialNR(-1, ang[i].k0, ang[i].k1, G_BABUSHKIN);
	c1 = ang[i].coeff;
	r += c1*c;
      }
      if (nmk >= 1 && t == 0) {
	c1 = mbk[0];
	if (c1 + 1.0 != 1.0) {
	  if (gauge == G_COULOMB) c1 /= aw;
	  r += c1;
	}
      }
    }    
    if (fabs(r) > 0.0) {
      r = OscillatorStrength(-1, te, r, NULL);
      bethe[0] = 2.0*r/te;
    } else {
      bethe[0] = 0.0;
    }
    ie = n_egrid;
    if (!msub) {
      born_cross = qkc[ie]*8.0;
    } else {
      rqk = qkc;
      bt = 0.0;
      for (t = -j1; t <= 0; t += 2) {
	for (h = -j2; h <= j2; h += 2) {
	  bt += rqk[ie];
	  if (t != 0) bt += rqk[ie];
	  rqk += n_egrid1;
	}
      }
      if (bt > 0) {
	p = 0;
	rqk = qkc;
	for (t = -j1; t <= 0; t += 2) {
	  for (h = -j2; h <= j2; h += 2) {
	    params[p] = rqk[ie]/bt;
	    p++;
	    rqk += n_egrid1;
	  }
	}
      }
      born_cross = bt*8.0;      
    }
    if (born_cross > 0) {
      c = egrid[ie] + bte;
      born_egrid = c/te;
      if (bethe[0] > 0) bethe[1] = born_cross - bethe[0]*log(born_egrid);
      else bethe[1] = born_cross;
      bethe[2] = egrid[ie]; 
    } else {
      bethe[0] = -1.0;
      bethe[1] = 0.0;
      bethe[2] = 0.0;
    }
  } else {
    bethe[0] = -1.0;
    bethe[1] = 0.0;
    bethe[2] = 0.0;
  }

  free(ang);
  /*
  cecache.az[ic] = NULL;
  cecache.nz[ic] = 0;
  if (nmk > 0) {
    free(mbk);
    cecache.nmk[ic] = 0;
    cecache.mbk[ic] = NULL;
  }
  */
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
	if (bethe[0] < 0.0) {
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
	UVIP3P(np, n_egrid, log_egrid, qkc, n_usr, log_usr, qkt);
	for (ie = 0; ie < n_usr; ie++) {
	  qkt[ie] = exp(qkt[ie]);
	}
      } else {
	for (ie = 0; ie < n_egrid; ie++) {
	  qkc[ie] *= 8.0;
	}
	UVIP3P(np, n_egrid, log_egrid, qkc, n_usr, log_usr, qkt);
      }
    } else if (qk_mode == QK_EXACT) {
      for (ie = 0; ie < n_usr; ie++) {
	qkt[ie] = 8.0*qkc[ie];
      }
    }
    RelativisticCorrection(0, qkt, NULL, bte, bethe[0]);
    return 1;
  } else {
    rqk = qkc;
    rqkt = qkt;
    p = 0;
    if (qk_mode == QK_FIT) {
      for (t = -j1; t <= 0; t += 2) {
	for (h = -j2; h <= j2; h += 2) {
	  for (ie = 0; ie < n_egrid; ie++) {
	    rqkt[ie] = 8.0*rqk[ie];
	  }	
	  if (rqkt[0] < EPS10 && rqkt[n_egrid-1] < EPS10) {
	    continue;
	  }
	  p++;
	  rqk += n_egrid1;
	  rqkt += n_usr;
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
	    UVIP3P(np, n_egrid, log_egrid, rqk, n_usr, log_usr, rqkt);
	    for (ie = 0; ie < n_usr; ie++) {
	      rqkt[ie] = exp(rqkt[ie]);
	    }
	    p++;
	    rqk += n_egrid1;
	    rqkt += n_usr;
	  }
	} 
      } else {
	for (t = -j1; t <= 0; t += 2) {
	  for (h = -j2; h <= j2; h += 2) {
	    for (ie = 0; ie < n_egrid; ie++) {
	      rqk[ie] *= 8.0;
	    }
	    UVIP3P(np, n_egrid, log_egrid, rqk, n_usr, log_usr, rqkt);
	    p++;
	    rqk += n_egrid1;
	    rqkt += n_usr;
	  }
	}
      }
    } else if (qk_mode == QK_EXACT) {
      for (t = -j1; t <= 0; t += 2) {
	for (h = -j2; h <= j2; h += 2) {	
	  for (ie = 0; ie < n_egrid; ie++) {
	    rqkt[ie] = 8.0*rqk[ie];
	  }
	  p++;
	  rqk += n_egrid1;
	  rqkt += n_usr;
	}
      }
    }  
    for (ie = 0; ie < n_usr; ie++) {
      ubt[ie] = 0.0;
    }
    rqkt = qkt;
    for (t = -j1; t <= 0; t += 2) {
      for (h = -j2; h <= j2; h += 2) {
	for (ie = 0; ie < n_egrid; ie++) {
	  ubt[ie] += rqkt[ie];
	  if (t != 0) ubt[ie] += rqkt[ie];
	}
	rqkt += n_usr;
      }
    }
    RelativisticCorrection(p, qkt, ubt, bte, bethe[0]);
    return p;
  }
}

void ProcessCECache(int msub, int iuta, TFILE *f) {
  int ic, iz, jz, ie, skip, ilow, iup, nparams;

  if (qk_mode == QK_FIT) 
    nparams = NPARAMS;
  else
    nparams = 0;

  if (!iuta) {
    /*
    int nb = GetNumBounds();
    int i, *kb;
    kb = malloc(sizeof(int)*nb);
    for (i = 0; i < nb; i++) kb[i] = 0;
    */
    ResetWidMPI();    
#pragma omp parallel default(shared) private(ic, iz, ilow, iup, skip)
    {
      for (ic = 0; ic < cecache.nc; ic++) {
	ilow = cecache.low[ic];
	iup = cecache.up[ic];
	skip = SkipMPI();
	if (skip) continue;
	cecache.nz[ic] = AngularZMix(&cecache.az[ic], ilow, iup,
				     -1, -1,
				     &cecache.nmk[ic],
				     &cecache.mbk[ic]);
	/*
	for (iz = 0; iz < cecache.nz[ic]; iz++) {
	  if (cecache.az[ic][iz].k0 != cecache.az[ic][iz].k1) {
	    kb[cecache.az[ic][iz].k0] |= 1;
	    kb[cecache.az[ic][iz].k1] |= 2;
	  } else {
	    kb[cecache.az[ic][iz].k0] |= 4;
	  }
	}
	*/
      }
    }
    /*
    int nks, nkd0, nkd1;
    nks = 0;
    nkd0 = 0;
    nkd1 = 0;
    for (i = 0; i < nb; i++) {
      if (kb[i] == 4) nks++;
      else {
	if (kb[i] & 1) nkd0++;
	if (kb[i] & 2) nkd1++;
      }
    }
    int *ks = NULL, *kd0 = NULL, *kd1 = NULL;
    if (nks > 0) {
      ks = malloc(sizeof(int)*nks);
    }
    if (nkd0 > 0) {
      kd0 = malloc(sizeof(int)*nkd0);
    }
    if (nkd1 > 0) {
      kd1 = malloc(sizeof(int)*nkd1);
    }    
    int is = 0, i0 = 0, i1 = 0;
    for (i = 0; i < nb; i++) {
      if (kb[i] == 4) ks[is++] = i;
      else {
	if (kb[i] & 1) kd0[i0++] = i;
	if (kb[i] & 2) kd1[i1++] = i;
      }
    }
    IDXARY *iks = NULL, *ikd0 = NULL, *ikd1 = NULL;
    if (nks > 0) {
      iks = malloc(sizeof(IDXARY));
      InitIdxAry(iks, nks, ks);
    }
    if (nkd0 > 0) {
      ikd0 = malloc(sizeof(IDXARY));
      InitIdxAry(ikd0, nkd0, kd0);
    }
    if (nkd1 > 0) {
      ikd1 = malloc(sizeof(IDXARY));
      InitIdxAry(ikd1, nkd1, kd1);
    }

    CEPKK *pk = NULL;
    int npk = nks + nkd0*nkd1;
    pk = malloc(sizeof(CEPKK)*npk);
    for (i = 0; i < npk; i++) {
      pk[i].kmin = (short)5000;
      pk[i].kmax = 0;
    }
    ORBITAL *orb0, *orb1;
    int j0, j1, jmin, jmax, k, j, m;
    for (i = 0; i < nks; i++) {
      orb0 = GetOrbital(ks[i]);
      j0 = GetJFromKappa(orb0->kappa);
      pk[i].kmin = 0;
      pk[i].kmax = Max(j0, pk[i].kmax);
    }
    for (i = 0; i < nkd0; i++) {
      orb0 = GetOrbital(kd0[i]);
      j0 = GetJFromKappa(orb0->kappa);
      for (j = 0; j < nkd1; j++) {
	orb1 = GetOrbital(kd1[j]);
	j1 = GetJFromKappa(orb1->kappa);
	jmin = abs(j0-j1)/2;
	jmax = (j0+j1)/2;
	k = nks + i0*nkd1+i1;
	pk[k].kmin = Max(pk[k].kmin, jmin);
	pk[k].kmax = Min(pk[k].kmax, jmax);
      }
    }
    for (i = 0; i < npk; i++) {
      pk[i].pk = malloc(sizeof(CEPK *)*(pk[i].kmax-pk[i].kmin+1));
      m = 0;
      for (k = pk[i].kmin; k <= pk[i].kmax; k++, m++) {
	pk[i].pk[m] = NULL;
      }
    }
    CEQKK *qk = NULL;
    int nks2 = nks*(1+nks)/2;
    int nkd01 = nkd0*nkd1;
    int nksd = nks*nkd01;
    int nkd2 = nkd01*(1+nkd01)/2;
    int nqk = nks*(1+nks)/2 + nks*nkd01 + nkd2;
    int p, q;
    qk = malloc(sizeof(CEQKK)*nqk);
    for (i = 0; i < nks; i++) {
      for (j = 0; j <= i; j++) {
	m = i*(i+1)/2 + j;
	qk[m].kmin = pk[i].kmin;
	qk[m].kmax = pk[i].kmax;
	qk[m].kminp = pk[j].kmin;
	qk[m].kmaxp = pk[j].kmax;
      }
      for (i0 = 0; i0 < nkd01; i0++) {
	m = nks2+i*nkd01+i0;
	qk[m].kmin = pk[i].kmin;
	qk[m].kmax = pk[i].kmax;
	qk[m].kminp = pk[i0].kmin;
	qk[m].kmaxp = pk[i0].kmax;
      }
    }
    for (i0 = 0; i < nkd01; i0++) {
      for (i1 = 0; i1 <= i0; i1++) {
	m = nks2+nks*nkd01 + i0*(i0+1)/2 + i1;
	qk[m].kmin = pk[i0].kmin;
	qk[m].kmax = pk[i0].kmax;
	qk[m].kminp = pk[i1].kmin;
	qk[m].kmaxp = pk[i1].kmax;
      }
    }

    for (i = 0; i < nqk; i++) {
      if (msub) {
	p = qk[i].kmax-qk[i].kmin+1;
	q = qk[i].kmaxp-qk[i].kminp+1;
	qk[i].qk = malloc(sizeof(double *)*p*q);
	m = 0;
	for (k = qk[i].kmin; k <= qk[i].kmax; k++) {
	  for (j = qk[i].kminp; j <= qk[i].kmaxp; j++, m++) {
	    qk[i].qk[m] = NULL;
	  }
	}
      } else {
	qk[i].kmin = Max(qk[i].kmin, qk[i].kminp);
	qk[i].kmax = Min(qk[i].kmax, qk[i].kmaxp);
	qk[i].qk = malloc(sizeof(double *)*(qk[i].kmax-qk[i].kmin+1));
	m = 0;
	for (k = qk[i].kmin; k <= qk[i].kmax; k++, m++) {
	  qk[i].qk[m] = NULL;
	}
      }
    }

    ResetWidMPI();
#pragma omp parallel default(shared) private(ic, skip, iz, ie)
    {
      for (ic = 0; ic < cecache.nc; ic++) {
	skip = SkipMPI();
	if (skip) continue;
	for (iz = 0; iz < cecache.nz[ic]; iz++) {
	  for (ie = 0; ie < n_egrid; ie++) {
	    CEPK *cepk;
	    int type = CERadialPk(&cepk, ie,
				  cecache.az[ic][iz].k0,
				  cecache.az[ic][iz].k1,
				  cecache.az[ic][iz].k, 1);
	  }
	}
      }
    }

    ResetWidMPI();
#pragma omp parallel default(shared) private(ic, iz, jz, skip)
    {
      for (ic = 0; ic < cecache.nc; ic++) {
	skip = SkipMPI();
	if (skip) continue;
	for (iz = 0; iz < cecache.nz[ic]; iz++) {
	  for (jz = iz; jz < cecache.nz[ic]; jz++) {
	    if (msub) {
	      double *rqe = CERadialQkMSubTable(cecache.az[ic][iz].k0,
						cecache.az[ic][iz].k1,
						cecache.az[ic][jz].k0,
						cecache.az[ic][jz].k1,
						cecache.az[ic][iz].k,
						cecache.az[ic][jz].k, 1);
	    } else {
	      if (cecache.az[ic][iz].k != cecache.az[ic][jz].k) continue;
	      skip = SkipMPI();
	      if (skip) continue;
	      double *rqe = CERadialQkTable(cecache.az[ic][iz].k0,
					    cecache.az[ic][iz].k1,
					    cecache.az[ic][jz].k0,
					    cecache.az[ic][jz].k1,
					    cecache.az[ic][iz].k, 1);
	    }
	  }
	}
      }
    }    
    */
  }

  ResetWidMPI();
#pragma omp parallel default(shared) private(ic, ilow, iup, skip, ie)
  {
    int nsub, m, k, ip, iempty;
    double bethe[3];
    double params[MAXMSUB*NPARAMS];
    double qkc[MAXMSUB*MAXNUSR];
    double e;
    CE_RECORD r;
    nsub = 1;
    if (msub) {
      r.params = (float *) malloc(sizeof(float)*nsub);
    } else if (qk_mode == QK_FIT) {
      m = nparams * nsub;
      r.params = (float *) malloc(sizeof(float)*m);
    }
    m = n_usr * nsub;
    r.strength = (float *) malloc(sizeof(float)*m);
    for (ic = 0; ic < cecache.nc; ic++) {
      skip = SkipMPI();
      if (skip) continue;
      ilow = cecache.low[ic];
      iup = cecache.up[ic];
      if (iuta) {
	k = CollisionStrengthUTA(qkc, params, &e, bethe, ilow, iup);
      } else {
	k = CollisionStrength(qkc, params, &e, bethe, ilow, iup, msub);
      }
      if (k < 0) continue;
      r.bethe = bethe[0];
      r.born[0] = bethe[1];
      r.born[1] = bethe[2];
      r.lower = ilow;
      r.upper = iup;
      r.nsub = k;
      if (r.nsub > nsub) {
	r.params = (float *) realloc(r.params, sizeof(float)*r.nsub);
	m = n_usr * r.nsub;
	r.strength = (float *) realloc(r.strength, sizeof(float)*m);
	nsub = r.nsub;
      }      
      if (msub) {
	for (m = 0; m < r.nsub; m++) {
	  r.params[m] = (float) params[m];
	}
      } else if (qk_mode == QK_FIT) {
	ip = 0;
	for (m = 0; m < r.nsub; m++) {
	  for (ie = 0; ie < nparams; ie++) {
	    r.params[ip] = (float) params[ip];
	    ip++;
	  }
	}
      }
      
      ip = 0;
      iempty = 1;
      for (m = 0; m < r.nsub; m++) {
	for (ie = 0; ie < n_usr; ie++) {
	  r.strength[ip] = (float) qkc[ip];
	  if (r.strength[ip]) iempty = 0;
	  ip++;
	}
      }
      if (iempty == 0) {
	WriteCERecord(f, &r);
      }
    }
    if (msub || qk_mode == QK_FIT) free(r.params);
    free(r.strength);
  }
}

int SaveExcitation(int nlow, int *low, int nup, int *up, int msub, char *fn) {
#ifdef PERFORM_STATISTICS
  STRUCT_TIMING structt;
  ANGULAR_TIMING angt;
  RECOUPLE_TIMING recouplet;
  RAD_TIMING radt;
#endif
  SYMMETRY *sym;
  STATE *st;
  CONFIG *cfg;
  int i, j, k, n, m, ie, ip, iempty, nsub;
  TFILE *f;
  double qkc[MAXMSUB*MAXNUSR];
  double params[MAXMSUB*NPARAMS];
  int *alev;
  LEVEL *lev1, *lev2;
  CE_RECORD r;
  CE_HEADER ce_hdr;
  F_HEADER fhdr;
  ARRAY subte;
  int isub, n_tegrid0, n_egrid0, n_usr0;
  int te_set, e_set, usr_set, iuta;
  double emin, emax, e, c;
  double e0, e1, te0, ei;
  double rmin, rmax, bethe[3];
  int nc, ilow, iup;

  iuta = IsUTA();
  if (iuta && msub) {
    printf("cannot call CETableMSub in UTA mode\n");
    return -1;
  }
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

  nc = OverlapLowUp(nlow, low, nup, up);

  emin = 1E10;
  emax = 1E-10;
  m = 0;
  for (i = 0; i < nlow; i++) {
    lev1 = GetLevel(low[i]);
    for (j = 0; j < nup; j++) {
      lev2 = GetLevel(up[j]);
      e = lev2->energy - lev1->energy;
      if (i < nlow-nc || j < nup-nc) e = fabs(e);
      if (e > 0) m++;
      if (e < emin && e > 0) emin = e;
      if (e > emax) emax = e;
    }
  }
  if (m == 0) {
    return 0;
  }

  /*
  if (!iuta) {
    AllocCECache(msub);
    PrepAngZStates(nlow, low, nup, up);
  }
  */
  ei = 1E31;
  if (iuta) {
    for (j = 0; j < nup; j++) {
      lev2 = GetLevel(up[j]);
      cfg = GetConfigFromGroup(lev2->iham, lev2->pb);
      k = OrbitalIndex(cfg->shells[0].n, cfg->shells[0].kappa, 0.0);
      e = -(GetOrbital(k)->energy);
      if (e < ei) ei = e;
    }
  } else {
    for (j = 0; j < nup; j++) {
      lev2 = GetLevel(up[j]);
      sym = GetSymmetry(lev2->pj);
      st = (STATE *) ArrayGet(&(sym->states), lev2->pb);
      if (st->kgroup < 0) {
	k = st->kcfg;
      } else {
	cfg = GetConfig(st);
	k = OrbitalIndex(cfg->shells[0].n, cfg->shells[0].kappa, 0.0);
      }
      e = -(GetOrbital(k)->energy);
      if (e < ei) ei = e;
    }
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
 
  if (msub) {
    pw_type = 1;
    qk_mode = QK_EXACT;
  } else {
    pw_type = 0;
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

  e0 = emin*0.999;
  fhdr.type = DB_CE;
  strcpy(fhdr.symbol, GetAtomicSymbol());
  fhdr.atom = GetAtomicNumber();
  f = OpenFile(fn, &fhdr);
  for (isub = 1; isub < subte.dim; isub++) {
    e1 = *((double *) ArrayGet(&subte, isub));
    if (isub == subte.dim-1) e1 = e1*1.001;
    emin = e1;
    emax = e0;
    m = 0;
    for (i = 0; i < nlow; i++) {
      lev1 = GetLevel(low[i]);
      for (j = 0; j < nup; j++) {
	lev2 = GetLevel(up[j]);
	e = lev2->energy - lev1->energy;
	if (i < nlow-nc || j < nup-nc) e = fabs(e);
	if (e < e0 || e >= e1) continue;
	if (e < emin) emin = e;
	if (e > emax) emax = e;
	m++;
      }
    }
    if (m == 0) {
      e0 = e1;
      continue;
    }
    if (!te_set) {
      e = emax/emin;  
      if (e < 1.001) {
	SetCETEGrid(1, 0.5*(emax+emin), emax);
      } else if (e < 1.5) {
	SetCETEGrid(2, emin, emax);
      } else if (e < 5.0) {
	if (m == 2) n_tegrid = 2; 
	else if (n_tegrid0 == 0) n_tegrid = 3;
	SetCETEGrid(n_tegrid, emin, emax);
      } else {
	if (m == 2) n_tegrid = 2;
	else if (n_tegrid0 == 0) n_tegrid = 4;
	SetCETEGrid(n_tegrid, emin, emax);
      }
    }

    e = 0.5*(emin + emax);
    emin = rmin*e;
    if (te0 > ei) {
      emax = rmax*te0;
      ce_hdr.te0 = te0;
    } else {
      emax = rmax*te0*3.0;
      ce_hdr.te0 = te0;
    }
    
    if (n_egrid0 == 0) {
      n_egrid = 6;
    }
    if (!e_set) {
      SetCEEGrid(n_egrid, emin, emax, ce_hdr.te0);
    }
    if (n_usr0 <= 0) {
      SetUsrCEEGridDetail(n_egrid, egrid);
      usr_egrid_type = 1;
    } else if (!usr_set) {
      SetUsrCEEGrid(n_usr, emin, emax, ce_hdr.te0);
      usr_egrid_type = 1;
    }
    if (n_egrid > MAXNE) {
      printf("n_egrid exceeded MAXNE=%d\n", MAXNE);
      return -1;
    }
    n_egrid1 = n_egrid + 1;
    ie = n_egrid;
    if (eborn > 0.0) {
      egrid[ie] = Max(te0,ei)*eborn;
    } else {
      egrid[ie] = -eborn;
    }
    if (egrid[ie] < 2*egrid[ie-1]) egrid[ie] = 2*egrid[ie-1];
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

    e = 0.0;
    c = GetResidualZ();
    if (xborn+1.0 != 1.0) {
      PrepCoulombBethe(1, n_tegrid, n_egrid, c, &e, tegrid, egrid,
		       pw_scratch.nkl, pw_scratch.kl, msub);
    }
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
    ResetWidMPI();
#pragma omp parallel default(shared) private(i, j, lev1, lev2, e, ilow, iup, k, qkc, params, bethe, r, m, ip, nsub, ie, iempty)
    {
    nsub = 1;
    if (msub) {
      r.params = (float *) malloc(sizeof(float)*nsub);
    } else if (qk_mode == QK_FIT) {
      m = ce_hdr.nparams * nsub;
      r.params = (float *) malloc(sizeof(float)*m);
    }
    m = ce_hdr.n_usr * nsub;
    r.strength = (float *) malloc(sizeof(float)*m);    
    //ic = 0;
    for (i = 0; i < nlow; i++) {
      lev1 = GetLevel(low[i]);
      for (j = 0; j < nup; j++) {
	lev2 = GetLevel(up[j]);	
	e = lev2->energy - lev1->energy;	
	ilow = low[i];
	iup = up[j];
	if (i < nlow-nc || j < nup-nc) {
	  if (e < 0) {
	    ilow = up[j];
	    iup = low[i];
	    e = -e;
	  }
	}	    
	if (e < e0 || e >= e1) continue;
	int skip = SkipMPI();
	if (skip) continue;
	if (iuta) {
	  k = CollisionStrengthUTA(qkc, params, &e, bethe, ilow, iup);
	} else {
	  k = CollisionStrength(qkc, params, &e, bethe, ilow, iup, msub); 
	}
	if (k < 0) continue;
	r.bethe = bethe[0];
	r.born[0] = bethe[1];
	r.born[1] = bethe[2];
	r.lower = ilow;
	r.upper = iup;
	r.nsub = k;
	if (r.nsub > nsub) {
	  r.params = (float *) realloc(r.params, sizeof(float)*r.nsub);
	  m = ce_hdr.n_usr * r.nsub;
	  r.strength = (float *) realloc(r.strength, sizeof(float)*m);
	  nsub = r.nsub;
	}

	if (msub) {
	  for (m = 0; m < r.nsub; m++) {
	    r.params[m] = (float) params[m];
	  }
	} else if (qk_mode == QK_FIT) {
	  ip = 0;
	  for (m = 0; m < r.nsub; m++) {
	    for (ie = 0; ie < ce_hdr.nparams; ie++) {
	      r.params[ip] = (float) params[ip];
	      ip++;
	    }
	  }
	}
      
	ip = 0;
	iempty = 1;
	for (m = 0; m < r.nsub; m++) {
	  for (ie = 0; ie < ce_hdr.n_usr; ie++) {
	    r.strength[ip] = (float) qkc[ip];
	    if (r.strength[ip]) iempty = 0;
	    ip++;
	  }
	}
	if (iempty == 0) {
	  WriteCERecord(f, &r);
	}
      }
    }
    if (msub || qk_mode == QK_FIT) free(r.params);
    free(r.strength);
    }
  /*
    cecache.low[ic] = ilow;
    cecache.up[ic] = iup;
    ic++;
    if (ic == maxcecache) {
    cecache.nc = ic;
    ic = 0;
    ProcessCECache(msub, iuta, f);
    }
    }
    }
    if (ic > 0) {
      cecache.nc = ic;
      ProcessCECache(msub, iuta, f);
    }
  */
    DeinitFile(f, &fhdr);
    e0 = e1;
    FreeExcitationQk();
    ReinitRadial(2);
  }
  
  ReinitExcitation(1);
  //FreeCECache(0);
  ArrayFreeLock(&subte, NULL);
  if (alev) free(alev);
  CloseFile(f, &fhdr);

  if (fpw) {
    fclose(fpw);
    fpw = NULL;
  }

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

int SaveExcitationEB(int nlow0, int *low0, int nup0, int *up0, char *fn) {
  int nlow, *low, nup, *up;
  int i, j, k, n, m, ie;
  CEF_RECORD r;
  CEF_HEADER ce_hdr;
  F_HEADER fhdr;
  LEVEL *lev1, *lev2;
  SYMMETRY *sym;
  STATE *st;
  CONFIG *cfg;
  ARRAY subte;
  int isub, n_tegrid0, n_egrid0;
  int te_set, e_set, iempty;
  double emin, emax, e, c;
  double e0, e1, te0, ei;
  double rmin, rmax, bethe[3];
  int nc, ilow, iup;
  TFILE *f;
  double qkc[MAXNE+1];
 
  if (GetLowUpEB(&nlow, &low, &nup, &up, nlow0, low0, nup0, up0) == -1)
    return 0;

  nc = OverlapLowUp(nlow, low, nup, up);

  emin = 1E10;
  emax = 1E-10;
  m = 0;
  for (i = 0; i < nlow; i++) {
    lev1 = GetEBLevel(low[i]);
    for (j = 0; j < nup; j++) {
      lev2 = GetEBLevel(up[j]);
      e = lev2->energy - lev1->energy;
      if (i < nlow-nc || j < nup-nc) e = fabs(e);
      if (e > 0) m++;
      if (e < emin && e > 0) emin = e;
      if (e > emax) emax = e;
    }
  }
  if (m == 0) {
    return 0;
  }

  ei = 1E31;
  for (j = 0; j < nup0; j++) {
    lev2 = GetLevel(up0[j]);    
    sym = GetSymmetry(lev2->pj);
    st = (STATE *) ArrayGet(&(sym->states), lev2->pb);
    if (st->kgroup < 0) {
      k = st->kcfg;
    } else {
      cfg = GetConfig(st);
      k = OrbitalIndex(cfg->shells[0].n, cfg->shells[0].kappa, 0.0);
    }
    e = -(GetOrbital(k)->energy);
    if (e < ei) ei = e;
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

  n_tegrid0 = n_tegrid;
  n_egrid0 = n_egrid;

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
 
  qk_mode = QK_EXACT;
  pw_type = 0;
  egrid_type = 1;

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
  
  e0 = emin*0.999;
  fhdr.type = DB_CEF;
  strcpy(fhdr.symbol, GetAtomicSymbol());
  fhdr.atom = GetAtomicNumber();
  f = OpenFile(fn, &fhdr);
  for (isub = 1; isub < subte.dim; isub++) {
    e1 = *((double *) ArrayGet(&subte, isub));
    if (isub == subte.dim-1) e1 = e1*1.001;
    emin = e1;
    emax = e0;
    m = 0;
    for (i = 0; i < nlow; i++) {
      lev1 = GetEBLevel(low[i]);
      for (j = 0; j < nup; j++) {
	lev2 = GetEBLevel(up[j]);
	e = lev2->energy - lev1->energy;
	if (i < nlow-nc || j < nup-nc) e = fabs(e);
	if (e < e0 || e >= e1) continue;
	if (e < emin) emin = e;
	if (e > emax) emax = e;
	m++;
      }
    }
    if (m == 0) {
      e0 = e1;
      continue;
    }
    if (!te_set) {
      e = emax/emin;  
      if (e < 1.001) {
	SetCETEGrid(1, 0.5*(emax+emin), emax);
      } else if (e < 1.5) {
	SetCETEGrid(2, emin, emax);
      } else if (e < 5.0) {
	if (m == 2) n_tegrid = 2; 
	else if (n_tegrid0 == 0) n_tegrid = 3;
	SetCETEGrid(n_tegrid, emin, emax);
      } else {
	if (m == 2) n_tegrid = 2;
	else if (n_tegrid0 == 0) n_tegrid = 4;
	SetCETEGrid(n_tegrid, emin, emax);
      }
    }

    e = 0.5*(emin + emax);
    emin = rmin*e;
    if (te0 > ei) {
      emax = rmax*te0;
      ce_hdr.te0 = te0;
    } else {
      emax = rmax*te0*3.0;
      ce_hdr.te0 = te0;
    }
    
    if (n_egrid0 == 0) {
      n_egrid = 6;
    }
    if (!e_set) {
      SetCEEGrid(n_egrid, emin, emax, ce_hdr.te0);
    }
    if (n_egrid > MAXNE) {
      printf("n_egrid exceeded MAXNE=%d\n", MAXNE);
      return -1;
    }
    n_egrid1 = n_egrid + 1;
    ie = n_egrid;
    if (eborn > 0.0) {
      egrid[ie] = Max(te0,ei)*eborn;
    } else {
      egrid[ie] = -eborn;
    }
    if (egrid[ie] < 2*egrid[ie-1]) egrid[ie] = 2*egrid[ie-1];

    e = 0.0;
    c = GetResidualZ();
    if (xborn+1.0 != 1.0) {
      PrepCoulombBethe(1, n_tegrid, n_egrid, c, &e, tegrid, egrid,
		       pw_scratch.nkl, pw_scratch.kl, 0);
    }
    ce_hdr.nele = GetNumElectrons(low0[0]);
    ce_hdr.n_tegrid = n_tegrid;
    ce_hdr.n_egrid = n_egrid;
    ce_hdr.tegrid = tegrid;
    ce_hdr.egrid = egrid;
    GetFields(&ce_hdr.bfield, &ce_hdr.efield, &ce_hdr.fangle);
    InitFile(f, &fhdr, &ce_hdr);  
    m = ce_hdr.n_egrid;

    ResetWidMPI();
#pragma omp parallel default(shared) private(i, j, r, lev1, lev2, e, ilow, iup, k, iempty, ie, qkc, bethe)
    {
    r.strength = (float *) malloc(sizeof(float)*m);    
    for (i = 0; i < nlow; i++) {
      lev1 = GetEBLevel(low[i]);
      for (j = 0; j < nup; j++) {
	lev2 = GetEBLevel(up[j]);
	e = lev2->energy - lev1->energy;	
	ilow = low[i];
	iup = up[j];
	if (i < nlow-nc || j < nup-nc) {
	  if (e < 0) {
	    ilow = up[j];
	    iup = low[i];
	    e = -e;
	  }
	}	    
	if (e < e0 || e >= e1) continue;
	int skip = SkipMPI();
	if (skip) continue;
	k = CollisionStrengthEB(qkc, &e, bethe, ilow, iup); 
	if (k < 0) continue;

	r.bethe = bethe[0];
	r.born[0] = bethe[1];
	r.born[1] = bethe[2];
	r.lower = ilow;
	r.upper = iup;

	iempty = 1;
	for (ie = 0; ie < ce_hdr.n_egrid; ie++) {
	  r.strength[ie] = (float) qkc[ie];
	  if (r.strength[ie]) iempty = 0;
	}
	if (iempty == 0) {
	  WriteCEFRecord(f, &r);
	}
      }
    }
    free(r.strength);
    }
    DeinitFile(f, &fhdr);
    e0 = e1;
    FreeExcitationQk();
    ReinitRadial(2);
  }

  ReinitExcitation(1);

  ArrayFreeLock(&subte, NULL);
  free(low);
  free(up);

  CloseFile(f, &fhdr);

  if (fpw) {
    fclose(fpw);
    fpw = NULL;
  }

  return 0;
}

int SaveExcitationEBD(int nlow0, int *low0, int nup0, int *up0, char *fn) {
  int nlow, *low, nup, *up;
  int i, j, k, n, m, ie;
  CEMF_RECORD r;
  CEMF_HEADER ce_hdr;
  F_HEADER fhdr;
  LEVEL *lev1, *lev2;
  SYMMETRY *sym;
  STATE *st;
  CONFIG *cfg;
  ARRAY subte;
  int isub, n_tegrid0, n_egrid0;
  int te_set, e_set, iempty;
  double emin, emax, e, c;
  double e0, e1, te0, ei;
  double rmin, rmax;
  int nc, ilow, iup;
  TFILE *f;
  double *qkc;
  double *bethe, *born;
 
  if (GetLowUpEB(&nlow, &low, &nup, &up, nlow0, low0, nup0, up0) == -1)
    return 0;

  nc = OverlapLowUp(nlow, low, nup, up);

  emin = 1E10;
  emax = 1E-10;
  m = 0;
  for (i = 0; i < nlow; i++) {
    lev1 = GetEBLevel(low[i]);
    for (j = 0; j < nup; j++) {
      lev2 = GetEBLevel(up[j]);
      e = lev2->energy - lev1->energy;
      if (i < nlow-nc || j < nup-nc) e = fabs(e);
      if (e > 0) m++;
      if (e < emin && e > 0) emin = e;
      if (e > emax) emax = e;
    }
  }
  if (m == 0) {
    return 0;
  }

  ei = 1E31;
  for (j = 0; j < nup0; j++) {
    lev2 = GetLevel(up0[j]);    
    sym = GetSymmetry(lev2->pj);
    st = (STATE *) ArrayGet(&(sym->states), lev2->pb);
    if (st->kgroup < 0) {
      k = st->kcfg;
    } else {
      cfg = GetConfig(st);
      k = OrbitalIndex(cfg->shells[0].n, cfg->shells[0].kappa, 0.0);
    }
    e = -(GetOrbital(k)->energy);
    if (e < ei) ei = e;
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

  n_tegrid0 = n_tegrid;
  n_egrid0 = n_egrid;

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
 
  qk_mode = QK_EXACT;
  pw_type = 1;
  egrid_type = 1;

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
  
  m = n_thetagrid*n_phigrid;
  qkc = (double *) malloc(sizeof(double)*(MAXNE+2)*m);
  bethe = (double *) malloc(sizeof(double)*m);
  born = (double *) malloc(sizeof(double)*(m+1));
  e0 = emin*0.999;
  fhdr.type = DB_CEMF;
  strcpy(fhdr.symbol, GetAtomicSymbol());
  fhdr.atom = GetAtomicNumber();
  f = OpenFile(fn, &fhdr);
  for (isub = 1; isub < subte.dim; isub++) {
    e1 = *((double *) ArrayGet(&subte, isub));
    if (isub == subte.dim-1) e1 = e1*1.001;
    emin = e1;
    emax = e0;
    m = 0;
    for (i = 0; i < nlow; i++) {
      lev1 = GetEBLevel(low[i]);
      for (j = 0; j < nup; j++) {
	lev2 = GetEBLevel(up[j]);
	e = lev2->energy - lev1->energy;
	if (i < nlow-nc || j < nup-nc) e = fabs(e);
	if (e < e0 || e >= e1) continue;
	if (e < emin) emin = e;
	if (e > emax) emax = e;
	m++;
      }
    }
    if (m == 0) {
      e0 = e1;
      continue;
    }
    if (!te_set) {
      e = emax/emin;  
      if (e < 1.001) {
	SetCETEGrid(1, 0.5*(emax+emin), emax);
      } else if (e < 1.5) {
	SetCETEGrid(2, emin, emax);
      } else if (e < 5.0) {
	if (m == 2) n_tegrid = 2; 
	else if (n_tegrid0 == 0) n_tegrid = 3;
	SetCETEGrid(n_tegrid, emin, emax);
      } else {
	if (m == 2) n_tegrid = 2;
	else if (n_tegrid0 == 0) n_tegrid = 4;
	SetCETEGrid(n_tegrid, emin, emax);
      }
    }

    e = 0.5*(emin + emax);
    emin = rmin*e;
    if (te0 > ei) {
      emax = rmax*te0;
      ce_hdr.te0 = te0;
    } else {
      emax = rmax*te0*3.0;
      ce_hdr.te0 = te0;
    }
    
    if (n_egrid0 == 0) {
      n_egrid = 6;
    }
    if (!e_set) {
      SetCEEGrid(n_egrid, emin, emax, ce_hdr.te0);
    }
    if (n_egrid > MAXNE) {
      printf("n_egrid exceeded MAXNE=%d\n", MAXNE);
      return -1;
    }
    n_egrid1 = n_egrid + 2;
    ie = n_egrid+1;
    if (eborn > 0.0) {
      egrid[ie] = Max(te0,ei)*eborn;
    } else {
      egrid[ie] = -eborn;
    }
    if (egrid[ie] < 2*egrid[ie-2]) egrid[ie] = 2*egrid[ie-2];
    egrid[ie-1] = 0.7*egrid[ie];

    e = 0.0;
    c = GetResidualZ();
    if (xborn+1.0 != 1.0) {
    PrepCoulombBethe(1, n_tegrid, n_egrid, c, &e, tegrid, egrid,
		     pw_scratch.nkl, pw_scratch.kl, 1);
    }
    ce_hdr.nele = GetNumElectrons(low0[0]);
    ce_hdr.n_tegrid = n_tegrid;
    ce_hdr.n_egrid = n_egrid;
    ce_hdr.n_thetagrid = n_thetagrid;
    ce_hdr.n_phigrid = n_phigrid;
    ce_hdr.tegrid = tegrid;
    ce_hdr.egrid = egrid;
    ce_hdr.thetagrid = thetagrid;
    ce_hdr.phigrid = phigrid;
    GetFields(&ce_hdr.bfield, &ce_hdr.efield, &ce_hdr.fangle);
    InitFile(f, &fhdr, &ce_hdr);  
    m = n_thetagrid*n_phigrid;
    r.strength = (float *) malloc(sizeof(float)*m*n_egrid);
    r.born = (float *) malloc(sizeof(float)*m);
    r.bethe = (float *) malloc(sizeof(float)*(m+1));
    
    for (i = 0; i < nlow; i++) {
      lev1 = GetEBLevel(low[i]);
      for (j = 0; j < nup; j++) {
	lev2 = GetEBLevel(up[j]);
	e = lev2->energy - lev1->energy;	
	ilow = low[i];
	iup = up[j];
	if (i < nlow-nc || j < nup-nc) {
	  if (e < 0) {
	    ilow = up[j];
	    iup = low[i];
	    e = -e;
	  }
	}	    
	if (e < e0 || e >= e1) continue;
	k = CollisionStrengthEBD(qkc, &e, bethe, born, ilow, iup); 
	if (k < 0) continue;
	
	for (ie = 0; ie < m; ie++) {
	  r.bethe[ie] = bethe[ie];
	}
	for (ie = 0; ie <= m; ie++) {
	  r.born[ie] = born[ie];
	}
	r.lower = ilow;
	r.upper = iup;
	iempty = 1;
	for (k = 0; k < m; k++) {	  
	  for (ie = 0; ie < n_egrid; ie++) {
	    r.strength[ie+k*n_egrid] = (float) qkc[ie+k*n_egrid1];
	    if (r.strength[ie+k*n_egrid]) iempty = 0;
	  }
	}
	if (iempty == 0) {
	  k = WriteCEMFRecord(f, &r);
	}
      }
    }
    free(r.strength);
    free(r.bethe);
    free(r.born);
    DeinitFile(f, &fhdr);
    e0 = e1;
    FreeExcitationQk();
    ReinitRadial(2);
  }    			   

  free(bethe);
  free(born);
  free(qkc);
  ReinitExcitation(1);

  ArrayFree(&subte, NULL);
  free(low);
  free(up);

  CloseFile(f, &fhdr);

  if (fpw) {
    fclose(fpw);
    fpw = NULL;
  }

  return 0;
}

int FreeExcitationQk(void) {
  MultiFreeData(qk_array, FreeExcitationQkData);
  MultiFreeData(pk_array, FreeExcitationPkData);
  return 0;
}

  
int InitExcitation(void) {
  int blocks1[] = {MULTI_BLOCK3,MULTI_BLOCK3,MULTI_BLOCK3};
  int blocks2[] = {MULTI_BLOCK5,MULTI_BLOCK5,MULTI_BLOCK5,
		   MULTI_BLOCK5,MULTI_BLOCK5};
  int ndim;

  ndim = 3;
  pk_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(pk_array, sizeof(CEPK), ndim, blocks1, "pk_array");

  ndim = 5;
  qk_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(qk_array, sizeof(double *), ndim, blocks2, "qk_array");

  if (fpw) {
    fclose(fpw);
    fpw = NULL;
  }
  n_egrid = 0;
  n_tegrid = 0;
  n_usr = 0;
  egrid[0] = -1.0;
  SetCEEGridLimits(0.05, 8.0, 0);
  usr_egrid[0] = -1.0;
  tegrid[0] = -1.0;  
  SetCEQkMode(QK_DEFAULT, 1E-3);
  SetCEPWOptions(EXCLQR, EXCLMAX, EXCLCB, EXCTOL);

  SetAngleGrid(0, 10, 0.0, PI);
  SetAngleGrid(1, 20, 0.0, TWO_PI);

  SetMaxCECache(-1);
  return 0;
}

int ReinitExcitation(int m) {
  
  if (m < 0) return 0;
  FreeExcitationQk();  
  if (fpw) {
    fclose(fpw);
    fpw = NULL;
  }
  n_egrid = 0;
  n_tegrid = 0;
  n_usr = 0;
  egrid[0] = -1.0;
  SetCEEGridLimits(0.05, 8.0, 0);
  usr_egrid[0] = -1.0;
  tegrid[0] = -1.0;  
  SetCEQkMode(QK_DEFAULT, 1E-3);
  SetCEPWOptions(EXCLQR, EXCLMAX, EXCLCB, EXCTOL);

  return 0;
}
