#include "excitation.h"
#include "cf77.h"

static char *rcsid="$Id: excitation.c,v 1.63 2003/10/14 21:26:16 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

#define MAXMSUB  16
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

#define NKINT 128
static double kgrid[NKINT];
static double log_kgrid[NKINT];
static double kint[NKINT];
static double log_kint[NKINT];
static double gos1[NKINT];
static double gos2[NKINT];
static double gost[NKINT];
static double gosint[NKINT];
static double xborn = XBORN;

#ifdef PERFORM_STATISTICS
static EXCIT_TIMING timing = {0, 0, 0};
#endif

static CEPW_SCRATCH pw_scratch = {1, MAXKL, 100, 5E-2, 0, 0, 10};

static MULTI *pk_array;
static MULTI *qk_array;

static void InitCEPK(void *p, int n) {
  CEPK *d;
  int i;
  
  d = (CEPK *) p;
  for (i = 0; i < n; i++) {
    d[i].nkl = -1;
  }
}

CEPW_SCRATCH *GetCEPWScratch(void) {
  return &pw_scratch;
}

int SetCEQkMode(int m, double tol) {
  if (m == QK_DEFAULT) qk_mode = QK_EXACT;
  else qk_mode = m;
  if (tol > 0.0) qk_fit_tolerance = tol;
  return 0;
}

int SetCEBorn(double x) {
  xborn = x;
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

int CERadialPk(CEPK **pk, int ie, int k0, int k1, int k) {
  int type, ko2, i, m, t, q;
  int kf0, kf1, kpp0, kpp1, km0, km1;
  int kl0, kl1, kl0p, kl1p;
  int j0, j1, kl_max, j1min, j1max;
  ORBITAL *orb0, *orb1;
  int index[4];
  double te, e0, e1, sd, se;
  int js1, js3, js[4], ks[4];
  int nkappa;
  short *kappa0, *kappa1;
  double *pkd, *pke;

#ifdef PERFORM_STATISTICS
  clock_t start, stop;
  start = clock();
#endif

  ko2 = k/2;  
  index[0] = ie;
  index[1] = k0;
  index[2] = k1;
  index[3] = ko2;
  
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

  *pk = (CEPK *) MultiSet(pk_array, index, NULL, InitCEPK);
  if ((*pk)->nkl >= 0) {
    return type;
  }

  nkappa = (MAXNKL)*(GetMaxRank()+1)*4;
  kappa0 = (short *) malloc(sizeof(short)*nkappa);
  kappa1 = (short *) malloc(sizeof(short)*nkappa);
  pkd = (double *) malloc(sizeof(double)*(nkappa*n_tegrid));
  pke = (double *) malloc(sizeof(double)*(nkappa*n_tegrid));

  e1 = egrid[ie];
  if (type >= 0 && type < CBMULTIPOLES) {
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

  for (t = 0; t < pw_scratch.nkl; t++) {
    kl0 = pw_scratch.kl[t];
    if (pw_scratch.kl[t] > kl_max) break;
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
	    e0 = e1 + te;
	    if (pw_type == 0) {
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
	    if (i == 0) {
	      if (1.0+sd == 1.0 && 1.0+se == 1.0) {
		break;
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

  (*pk)->nkl = t;
  (*pk)->nkappa = m;
  if (pw_type == 0) {
    (*pk)->kappa0 = realloc(kappa0, sizeof(short)*m);
    (*pk)->kappa1 = realloc(kappa1, sizeof(short)*m);
  } else {
    (*pk)->kappa0 = realloc(kappa1, sizeof(short)*m);
    (*pk)->kappa1 = realloc(kappa0, sizeof(short)*m);
  }
  (*pk)->pkd = realloc(pkd, sizeof(double)*q);
  (*pk)->pke = realloc(pke, sizeof(double)*q);  
    
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

int CERadialQkBorn(int k0, int k1, int k2, int k3, int k, 
		   double te, double e1, double *qk) {
  int p0, p1, p2, p3;
  int m0, m1, m2, m3;
  int j0, j1, j2, j3;
  int ko2, t, nk;
  double r, c0, c1, dk;
  double *g1, *g2, *x1, *x2;

  *qk = 0.0;
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
  if (IsOdd((m2+m3+k)/2) || !Triangle(j2, k, j3)) {
    return -1;
  }

  ko2 = k/2;
  r = ReducedCL(j0, k, j1) * ReducedCL(j2, k, j3);
  r *= (k+1.0)*(k+1.0);
  g1 = GeneralizedMoments(k0, k1, ko2);
  x1 = g1 + NGOSK;
  g2 = GeneralizedMoments(k2, k3, ko2);
  x2 = g2 + NGOSK;

  c0 = sqrt(2.0*(te + e1));
  c1 = sqrt(2.0*e1);
  nk = NKINT-1;
  kint[0] = c0 - c1;
  kint[nk] = c0 + c1;
  log_kint[0] = log(kint[0]);
  log_kint[nk] = log(kint[nk]);
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

  *qk = dk*Simpson(gosint, 0, nk-1);
  
  return ko2;
}
  
int CERadialQkBornMSub(int k0, int k1, int k2, int k3, int k, int kp,
		       double te, double e1, 
		       int nq, int *q, double *qk) {
  int p0, p1, p2, p3;
  int m0, m1, m2, m3;
  int j0, j1, j2, j3;
  int ko2, ko2p, t, nk;
  int nudiff, mu1, mu2, ierr, ipqa[MAXMSUB];
  int kkp, iq;
  double xc, theta, dnu1, pqa[MAXMSUB];
  double r, c0, c1, c01, dk;
  double *g1, *g2, *x1, *x2;
  double gosm1[MAXMSUB][NKINT];
  double gosm2[MAXMSUB][NKINT];
  
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

  c0 = 2.0*(te+e1);
  c1 = 2.0*e1;
  c01 = c0 - c1;
  c0 = sqrt(c0);
  c1 = sqrt(c1);
  nk = NKINT-1;
  kint[0] = c0 - c1;
  kint[nk] = c0 + c1;
  log_kint[0] = log(kint[0]);
  log_kint[nk] = log(kint[nk]);
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

double *CERadialQkTable(int k0, int k1, int k2, int k3, int k) {
  int type, t, ie, ite, ipk, ipkp, nqk, ieb;
  int i, j, kl0, kl1, kl, nkappa, nkl, nkappap, nklp;
  CEPK *cepk, *cepkp, *tmp;
  short *kappa0, *kappa1, *kappa0p, *kappa1p;
  double *pkd, *pke, *pkdp, *pkep;
  double r, rd, s, b;
  double qk[MAXNKL], dqk[MAXNKL];
  double rq[MAXNTE][MAXNE], e1, te, te0;
  double *rqc, **p, *ptr;
  int index[6];
  int np = 3, one = 1;
  double logj;

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
  
  p = (double **) MultiSet(qk_array, index, NULL, InitPointerData);
  if (*p) {
    return *p;
  }   
  
  nqk = n_tegrid*n_egrid;
  *p = (double *) malloc(sizeof(double)*(nqk+1));
  rqc = *p;

  te0 = -GetOrbital(k0)->energy;
  te = -GetOrbital(k1)->energy;
  te0 = Max(te0, te);
  te = -GetOrbital(k2)->energy;
  te0 = Max(te0, te);
  te = -GetOrbital(k3)->energy;
  te0 = Max(te0, te);
  ieb = 0;
  for (ie = 0; ie < n_egrid; ie++) {
    e1 = egrid[ie];
    type = CERadialPk(&cepk, ie, k0, k1, k);
    if (k2 != k0 || k3 != k1) {
      type = CERadialPk(&cepkp, ie, k2, k3, k);
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
      
      if (ieb == 0) {
	nklp = nkl-1;
	if (type >= CBMULTIPOLES) {
	  b = GetCoulombBetheAsymptotic(te, e1);
	} else if (type >= 0) {
	  b = (GetCoulombBethe(0, ite, ie, type, 1))[nklp];
	  if (b < 0 || IsNan(b)) b = GetCoulombBetheAsymptotic(te, e1);
	} else {
	  b = 0.0;
	}
	s = dqk[nklp]*b;
	if (ite == 0 &&
	    ((xborn < 0 && rd && -xborn < s/rd) ||
	     (xborn > 0 && xborn < e1/te0))) {
	  ieb = 1;
	} else {
	  rq[ite][ie] = r + s;
	}
      }
      if (ieb) {
	type = CERadialQkBorn(k0, k1, k2, k3, k,
			      te, e1, &(rq[ite][ie]));
	rq[ite][ie] += r-rd;
      }
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

double *CERadialQkMSubTable(int k0, int k1, int k2, int k3, int k, int kp) {
  int type1, type2, kl, nqk;
  int i, j, kl0, klp0, kl0_2, klp0_2, kl1;
  CEPK *cepk, *cepkp;
  int nkappa, nkappap, nkl, nklp;
  short  *kappa0, *kappa1, *kappa0p, *kappa1p;
  double *pkd, *pke, *pkdp, *pkep, *ptr;
  int km0, km1, j0, jp0, j1, kmp0, kmp1, km0_m, kmp0_m;
  int mi, mf, t, c0, cp0;
  double r, rd, e0, e1, te, s, sd, b, te0;
  double pha0, phap0;
  double s3j1, s3j2, s3j3, s3j4;
  int ie, ite, q[MAXMSUB], nq, iq, ipk, ipkp, ieb;
  double qk[MAXMSUB][MAXNKL], dqk[MAXMSUB][MAXNKL];
  double rq[MAXMSUB][MAXNTE][MAXNE], drq[MAXMSUB][MAXNTE][MAXNE];
  double rqt[MAXMSUB];
  double *rqc, **p;
  int index[6];
  int np = 3, one = 1;
  double logj;

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

  p = (double **) MultiSet(qk_array, index, NULL, InitPointerData);
  if (*p) {
    return *p;
  }
 
  nq = Min(k, kp)/2 + 1;
  q[0] = 0;
  for (iq = 1; iq < nq; iq++) {
    q[iq] = q[iq-1] + 2;
  }
  pkdp = NULL;
  pkep = NULL;
  nqk = nq*n_tegrid*n_egrid;
  *p = (double *) malloc(sizeof(double)*(nqk+1));
  rqc = *p;

  te0 = -GetOrbital(k0)->energy;
  te = -GetOrbital(k1)->energy;
  te0 = Max(te0, te);
  te = -GetOrbital(k2)->energy;
  te0 = Max(te0, te);
  te = -GetOrbital(k3)->energy;
  te0 = Max(te0, te);
  ieb = 0;
  for (ie = 0; ie < n_egrid; ie++) {
    e1 = egrid[ie];

    type1 = CERadialPk(&cepk, ie, k0, k1, k);
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
      type2 = CERadialPk(&cepkp, ie, k2, k3, kp);
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
	      if (-q[iq] <= k && -q[iq] <= kp) {
		mf = mi + q[iq]; 
		s3j3 = W3j(j0, k, j1, -mi, -q[iq], mf); 
		s3j4 = W3j(jp0, kp, j1, -mi, -q[iq], mf); 
		rqt[iq] += s3j1*s3j2*s3j3*s3j4; 
	      } 
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

      if (ieb == 0) {
	i = nkl - 1;
	r = 0.0;
	for (iq = 0; iq < nq; iq++) {
	  if (type1 >= CBMULTIPOLES) {
	    b = GetCoulombBetheAsymptotic(te, e1);	  
	  } else if (type1 >= 0) {
	    if (abs(q[iq]) == 2) {
	      b = (GetCoulombBethe(0, ite, ie, type1, 1))[i];
	    } else if (q[iq] == 0) {
	      b = (GetCoulombBethe(0, ite, ie, type1, 2))[i];
	    }
	  } else {	  
	    b = 0.0;
	  }
	  s = dqk[iq][i]*b;
	  rqt[iq] = s;
	  if (xborn < 0 && drq[iq][ite][ie]) {
	    s /= drq[iq][ite][ie];
	    if (s > r) r = s;
	  }
	}
	if (ite == 0 && 
	    ((xborn < 0 && -xborn < r) ||
	     (xborn > 0 && xborn < e1/te0))) {
	  ieb = 1;
	} else {
	  for (iq = 0; iq < nq; iq++) {
	    rq[iq][ite][ie] += rqt[iq];
	  }
	}
      }
      
      if (ieb) {
	type1 = CERadialQkBornMSub(k0, k1, k2, k3, k, kp, te, e1, 
				   nq, q, rqt);
	for (iq = 0; iq < nq; iq++) {
	  rq[iq][ite][ie] += rqt[iq] - drq[iq][ite][ie];
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
	j += n_egrid;
      }
      UVIP3P(np, n_tegrid, xte, rq, nd, &x0, &rqc[i]);
    }
  }

  return type;
}

int CERadialQkMSub(double *rqc, double te, int k0, int k1, int k2, int k3, 
		   int k, int kp) {
  int i, np, nd, iq, n;
  int j, m, type, nq;
  double *rqe, rq[MAXNTE];
  double *xte, x0;
  
  rqe = CERadialQkMSubTable(k0, k1, k2, k3, k, kp);
  nq = Min(k, kp)/2 + 1;

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
    if (type == 0 || type == 1) {
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
	UVIP3P(np, n_tegrid, xte, rq, nd, &x0, &rqc[i]);
      }
      rqe += n;
      rqc += n_egrid;
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

int CollisionStrength(double *qkt, double *params, double *e, double *bethe,
		      int lower, int upper, int msub) {
  int i, j, t, h, p, m, type, ty;  
  LEVEL *lev1, *lev2;
  double te, c, r, s3j;
  ANGULAR_ZMIX *ang;
  int nz, j1, j2, ie, nk, np, nq, kkp, q[MAXMSUB];
  double rq[MAXMSUB*(MAXNE+1)], qkc[MAXMSUB*(MAXNE+1)];
  double *rqk, tol;
  double c1, c2, **g, *ck, dk, kmin, kmax;
  double *g1, *g2, *x1, *x2;
  int ierr, ipvt[NPARAMS];
  int lwa=5*NPARAMS+MAXNE;
  double wa[5*NPARAMS+MAXNE];
  double fvec[MAXNE], fjac[MAXNE*NPARAMS];
  double born_egrid, born_cross;

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
			    ang[j].k0, ang[j].k1, ang[i].k, ang[j].k);
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
	      m = abs(m)*n_egrid/2;
	      for (ie = 0; ie < n_egrid; ie++) {
		rqk[ie] += c*rq[m+ie]*s3j;
	      }
	    }
	    rqk += n_egrid;
	  }
	}
      }
    }
  }

  for (t = 0; t < NKINT; t++) {
    gost[t] = 0.0;
  }
  r = 0;
  if (msub) {
    for (t = 0; t < MAXMSUB; t++) {
      params[t] = 0.0;
    }
  }
  if (type >= 0) {
    nk = NGOSK-1;
    kmin = 1E30;
    kmax = -1E30;
    g = (double **) malloc(sizeof(double *)*nz);
    ck = (double *) malloc(sizeof(double)*nz);
    for (i = 0; i < nz; i++) {
      p = GetOrbital(ang[i].k0)->kappa;
      GetJLFromKappa(p, &t, &h);
      m = GetOrbital(ang[i].k1)->kappa;
      GetJLFromKappa(m, &np, &ie);
      if (IsOdd((h+ang[i].k+ie)/2)) {
	ck[i] = 0.0;
	g[i] = NULL;
      } else {
	ck[i] = ReducedCL(t, ang[i].k, np);
	g[i] = GeneralizedMoments(ang[i].k0, ang[i].k1, ang[i].k/2);
	x1 = g[i] + NGOSK;
	if (x1[0] < kmin) kmin = x1[0];
	if (x1[nk] > kmax) kmax = x1[nk];
      }
    }
	
    nk = NKINT-1.0;
    dk = (kmax - kmin)/nk;
    log_kgrid[0] = kmin;
    kgrid[0] = exp(log_kgrid[0]);
    log_kgrid[nk] = kmax;
    kgrid[nk] = exp(log_kgrid[nk]);
    for (t = 1; t < nk; t++) {
      log_kgrid[t] = log_kgrid[t-1] + dk;
      kgrid[t] = exp(log_kgrid[t]);
    }
    nk = NKINT;
    for (i = 0; i < nz; i++) {
      g1 = g[i];
      if (g1 == NULL) continue;
      x1 = g1 + NGOSK;
      InterpolateGOS(NGOSK, x1, g1, nk, log_kgrid, gos1);
      for (j = i; j < nz; j++) {
	if (ang[j].k != ang[i].k) continue;
	g2 = g[j];
	if (g2 == NULL) continue;
	x2 = g2 + NGOSK;
	InterpolateGOS(NGOSK, x2, g2, nk, log_kgrid, gos2);
	c = ang[i].coeff*ang[j].coeff;
	if (i != j) c *= 2.0;
	c *= 2.0*ck[i]*ck[j];
	if (ang[i].k == 2) {
	  r += c * (RadialMoments(1, ang[i].k0, ang[i].k1)*
		    RadialMoments(1, ang[j].k0, ang[j].k1));
	}
	c *= ang[i].k + 1.0;
	for (t = 0; t < nk; t++) {
	  gost[t] += c*gos1[t]*gos2[t];
	}
      }
    }
    free(g);
    free(ck);

    c = 0.0;
    for (t = 0; t < nk; t++) {
      if (gost[t] > c) {
	c = gost[t];
      }
    }
    if (c <= 0.0) {
      bethe[0] = -1.0;
      bethe[1] = 0.0;
      bethe[2] = 0.0;
    } else {
      r /= 3.0;
      bethe[0] = r*2.0;
      c *= EPS8;
      nk = NKINT-1;
      for (i = nk; i >= 0; i--) {
	if (gost[i] > c) break;
      }
      c = kgrid[i];
      born_egrid = 2.0*c*c/te;
      c1 = sqrt(2.0*te*born_egrid);
      c2 = sqrt(2.0*te*(born_egrid-1.0));
      kint[0] = c1 - c2;
      kint[nk] = c1 + c2;
      log_kint[0] = log(kint[0]);
      log_kint[nk] = log(kint[nk]);
      dk = (log_kint[nk] - log_kint[0])/nk;
      for (t = 1; t < nk; t++) {
	log_kint[t] = log_kint[t-1] + dk;
	kint[t] = exp(log_kint[t]);
      }
      nk = NKINT;
      InterpolateGOS(nk, log_kgrid, gost, nk, log_kint, gosint);
      born_cross = 4.0*dk*Simpson(gosint, 0, nk-1);
      if (bethe[0] > 0) bethe[1] = born_cross - bethe[0]*log(born_egrid);
      else bethe[1] = born_cross;
      bethe[2] = (born_egrid-1.0)*te;
      if (msub) {
	g2 = rq+MAXMSUB;
	for (t = 0; t < MAXMSUB; t++) g2[t] = 0.0;
	for (i = 0; i < nz; i++) {
	  for (j = i; j < nz; j++) {
	    c = ang[i].coeff * ang[j].coeff;
	    if (i != j) c *= 2.0;	  
	    nq = Min(ang[i].k, ang[j].k)/2 + 1;
	    q[0] = 0;
	    for (t = 1; t < nq; t++) {
	      q[t] = q[t-1] + 2;
	    }
	    ty = CERadialQkBornMSub(ang[i].k0, ang[i].k1, 
				    ang[j].k0, ang[j].k1, 
				    ang[i].k, ang[j].k,
				    te, bethe[2], nq, q, rq);
	    nq = Min(ang[i].k, ang[j].k);
	    kkp = (ang[i].k + ang[j].k)/2;
	    p = 0;
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
		  m = abs(m)/2;
		  g2[p] += c*rq[m]*s3j;
		  p++;
		}
	      }
	    }
	  }
	}
	p = 0;
	for (t = -j1; t <= 0; t += 2) {
	  for (h = -j2; h <= j2; h += 2) {
	    params[p] = 8.0*g2[p]/born_cross;
	    p++;
	  }
	}
      }
    }      
  } else {
    bethe[0] = -1.0;
    bethe[1] = 0.0;
    bethe[2] = 0.0;
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
	    UVIP3P(np, n_egrid, log_egrid, rqk, n_usr, log_usr, qkt);
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
	    UVIP3P(np, n_egrid, log_egrid, rqk, n_usr, log_usr, qkt);
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
  SYMMETRY *sym;
  STATE *st;
  CONFIG *cfg;
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
  double emin, emax, e, c;
  double e0, e1, te0, ei;
  double rmin, rmax, bethe[3];

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

  ei = 1E31;
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

  e0 = emin;
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
      if (e < 1.1) {
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
    if (msub) {
      r.params = (float *) malloc(sizeof(float)*nsub);
    } else if (qk_mode == QK_FIT) {
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
	r.bethe = bethe[0];
	r.born[0] = bethe[1];
	r.born[1] = bethe[2];
	r.lower = low[i];
	r.upper = up[j];
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
	for (m = 0; m < r.nsub; m++) {
	  for (ie = 0; ie < ce_hdr.n_usr; ie++) {
	    r.strength[ip] = (float) qkc[ip];
	    ip++;
	  }
	}
	WriteCERecord(f, &r);
      }
    }
    if (msub || qk_mode == QK_FIT) free(r.params);
    free(r.strength);
    DeinitFile(f, &fhdr);
    e0 = e1;
    FreeExcitationQk();
    ReinitRadial(2);
  }

  ReinitExcitation(1);

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

int FreeExcitationPk(int ie) {
  ARRAY *b;

  b = pk_array->array;
  if (b == NULL) return 0;
  if (ie < 0) {
    MultiFreeData(b, pk_array->ndim, FreeExcitationPkData);
  } else {
    b = (ARRAY *) ArrayGet(b, ie);
    if (b) {
      MultiFreeData(b, pk_array->ndim - 1, FreeExcitationPkData);
    }
  }

  return 0;
}

void FreeExcitationQkData(void *p) {
  double *dp;

  dp = *((double **) p);
  free(dp);
  *((double **) p) = NULL;
}

int FreeExcitationQk(void) {
  ARRAY *b;
  b = qk_array->array;
  if (b == NULL) return 0;
  MultiFreeData(b, qk_array->ndim, FreeExcitationQkData);
  FreeExcitationPk(-1);
  
  return 0;
}
  
int InitExcitation(void) {
  int blocks1[] = {MULTI_BLOCK4,MULTI_BLOCK4,MULTI_BLOCK4,MULTI_BLOCK4};
  int blocks2[] = {MULTI_BLOCK6,MULTI_BLOCK6,MULTI_BLOCK6,
		   MULTI_BLOCK6,MULTI_BLOCK6,MULTI_BLOCK6};
  int ndim;

  ndim = 4;
  pk_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(pk_array, sizeof(CEPK), ndim, blocks1);

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
  SetCEPWOptions(EXCLQR, EXCLMAX, EXCLCB, EXCTOL);

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
  SetCEPWOptions(EXCLQR, EXCLMAX, EXCLCB, EXCTOL);

  return 0;
}