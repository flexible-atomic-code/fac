#include "excitation.h"

#define MAXMSUB  200

static int egrid_type = -1;
static int usr_egrid_type = -1;
static int pw_type = -1;

static int interpolate_egrid = 0;
static int n_usr = 0;
static double usr_egrid[MAXNUSR];
static double log_usr[MAXNUSR];
static int n_egrid = 0;
static double egrid[MAXNE];
static double log_egrid[MAXNE];
static int n_tegrid = 0;
static double tegrid[MAXNTE];
static double log_te[MAXNTE];

static EXCIT_TIMING timing = {0, 0, 0};

static CEPW_SCRATCH pw_scratch = {1, MAXKL, 100, 5E-2, 0, 0, 10};

static MULTI *pk_array;
static MULTI *kappa0_array;
static MULTI *kappa1_array;


void uvip3p_(int *np, int *ndp, double *x, double *y, 
	     int *n, double *xi, double *yi);


CEPW_SCRATCH *GetCEPWScratch() {
  return &pw_scratch;
}

int SetCEEGridType(int utype, int etype, int ltype) {
  if (etype >= 0) egrid_type = etype;
  if (utype >= 0) usr_egrid_type = utype;
  if (ltype >= 0) pw_type = ltype;
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
    abort();
  }
  pw_scratch.max_kl = max;
  pw_scratch.kl_cb = kl_cb;
  pw_scratch.tolerence = tol;
  pw_scratch.nkl0 = 1;
  pw_scratch.kl[0] = 0;
  pw_scratch.log_kl[0] = -100.0;
  pw_scratch.nkl = 0;
  return 0;
}

int SetCEPWGrid(int ns, int *n, int *step) {
  if (pw_scratch.nkl0 <= 0) SetCEPWOptions(0, 1000, 100, 5E-2);
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
  int type, ko2, i, j, m, t, q;
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
  double **p, s;
  short **kp0, **kp1;
  double eps;

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
  max0 = Max(pw_scratch.ns, max0);
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
  eps = pw_scratch.tolerence;
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

int CERadialQk(double *rq, int ie, double te, int k0, 
		  int k1, int k2, int k3, int k) {
  int type, t, swap;
  int i, j, kl0, kl1, kl, nkappa, nkl, nkappap, nklp;
  short *kappa0, *kappa1, *kappa0p, *kappa1p, *tmp;
  double *pk, *pkp, r, s, b;
  double *pk1, *pk2, x1[MAXNTE];  
  double *qk, *y2, e1;

#ifdef PERFORM_STATISTICS
  clock_t start, stop;

  start = clock();
#endif

  qk = pw_scratch.qk;
  y2 = pw_scratch.y2;

  pk2 = NULL;
  type = CERadialPk(&nkappa, &nkl, &pk, &kappa0, &kappa1, ie,
		    k0, k1, k);
  pk1 = malloc(sizeof(double)*nkappa);
  for (i = 0; i < nkappa; i++) {
    pk1[i] = InterpolatePk(te, type, pk);
    pk += n_tegrid;
  }
  nklp = nkl;
  if (k2 != k0 || k3 != k1) {
    type = CERadialPk(&nkappap, &nklp, &pkp, &kappa0p, &kappa1p, ie,
		      k2, k3, k);
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
    qk[i] = 0.0;
  }
  t = -1;
  kl0 = -1;

  for (i = 0; i < nkappa; i++) {
    if (pw_type == 0) kl = GetLFromKappa(kappa0[i]);
    else kl = GetLFromKappa(kappa1[i]);
    if (kl != kl0) {
      t++;
      kl0 = kl;
    }
    if (k2 == k0 && k3 == k1) {
      s = pk1[i]*pk1[i];
      qk[t] += s;
    } else {
      s = 0.0;
      for (j = 0; j < nkappap; j++) {
	if (kappa0[i] == kappa0p[j] && kappa1[i] == kappa1p[j]) {
	  s = pk1[i]*pk2[j];
	  break;
	}
      }
      qk[t] += s;
    }
  }

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
 
  free(pk1);
  if (pk2) {
    free(pk2);
  }

  if (type >= CBMULTIPOLES) {
    e1 = egrid[ie];
    if (egrid_type == 0) e1 -= te;
    b = GetCoulombBetheAsymptotic(te, e1);
    s = qk[nklp]*b;
    r = r + s;
  } else if (type >= 0) {
    for (i = 0; i < n_tegrid; i++) {
      x1[i] = (GetCoulombBethe(0, i, ie, type, 1))[nklp];
    }
    b = InterpolatePk(te, -1, x1);
    s = qk[nklp]*b;
    r = r + s;       
  }

#ifdef PERFORM_STATISTICS
  stop = clock();
  timing.rad_qk += stop-start;
#endif
  *rq = r;
  return type;
}

  
int CERadialQkMSub(double *rq, int ie, double te, int k0, int k1,
		   int k2, int k3, int k, int kp, 
		   int nq, int *q) {
  int have_pk2, type1, type2, kl, np, nt;
  int i, j, kl0, klp0, kl0_2, klp0_2, kl1;
  int nkappa, nkappap, nkl, nklp;
  short  *kappa0, *kappa1, *kappa0p, *kappa1p;
  double *pk1, *pk2;
  double *pk, *pkp, *qy2;
  int km0, km1, j0, jp0, j1, kmp0, kmp1, km0_m, kmp0_m;
  int mi, mf, t, c0, cp0;
  double r, e0, e1, s, b;
  double pha0, phap0;
  double s3j1, s3j2, s3j3, s3j4;
  double dpha[MAXNTE];
  int ite, iq;
  double **qk;

#ifdef PERFORM_STATISTICS
  clock_t start, stop;
#endif     

#ifdef PERFORM_STATISTICS
  start = clock();
#endif

  qy2 = pw_scratch.y2;

  pk2 = NULL;
  have_pk2 = 0;
  type1 = CERadialPk(&nkappa, &nkl, &pk, &kappa0, &kappa1, ie, 
		     k0, k1, k);
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
    type2 = type1;
  } else {
    type2 = CERadialPk(&nkappap, &nklp, &pkp, &kappa0p, &kappa1p, ie,
		       k2, k3, kp);
    pk2 = malloc(sizeof(double)*nkappap);
    for (i = 0; i < nkappap; i++) {
      pk2[i] = InterpolatePk(te, type2, pkp);
      pkp += n_tegrid;
    }
    have_pk2 = 1;
    if (nklp < nkl) nkl = nklp;
  }
  qk = (double **) malloc(sizeof(double *)*nq); 
  for (i = 0; i < nq; i++) { 
    qk[i] = (double *) malloc(nkl*sizeof(double)); 
    for (j = 0; j < nkl; j++) { 
      qk[i][j] = 0.0; 
    }   
  } 
    
  np = 3;
  nt = 1;
  kl = -1;
  i = -1;
  if (egrid_type == 0) {
    e0 = egrid[ie];
    e1 = e0 - te;
  } else {
    e1 = egrid[ie];
  }
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
	
	if (egrid_type == 1) {
	  for (ite = 0; ite < n_tegrid; ite++) { 
	    e0 = e1 + tegrid[ite]; 
	    c0 = OrbitalIndex(0, km0_m, e0); 
	    cp0 = OrbitalIndex(0, kmp0_m, e0); 
	    pha0 = GetPhaseShift(c0); 
	    phap0 = GetPhaseShift(cp0); 
	    dpha[ite] = pha0-phap0; 
	  } 
	  r = InterpolatePk(te, -1, dpha);
	} else {
	  c0 = OrbitalIndex(0, km0_m, e0); 
	  cp0 = OrbitalIndex(0, kmp0_m, e0);
	  pha0 = GetPhaseShift(c0); 
	  phap0 = GetPhaseShift(cp0);  
	  r = pha0 - phap0;
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

  i = nkl - 1;
  if (type1 == type2) {
    if (type1 >= CBMULTIPOLES) {
      b = GetCoulombBetheAsymptotic(te, e1);
      for (iq = 0; iq < nq; iq++) {
	s = qk[iq][i]*b;
	rq[iq] += s;
      }
    } else if (type1 >= 0) {
      for (iq = 0; iq < nq; iq++) {
	for (j = 0; j < n_tegrid; j++) {
	  dpha[j] = (GetCoulombBethe(0, j, ie, type1, iq+1))[i];
	}
	b = InterpolatePk(te, -1, dpha);
	s = qk[iq][i]*b;
	rq[iq] += s;
      }
    }
  } else if (type1 >= 0) {
    b = GetCoulombBetheAsymptotic(te, e1);
    for (iq = 0; iq < nq; iq++) {
      s = qk[iq][i]*b;
      rq[iq] += s;
    }
  }

  for (i = 0; i < nq; i++) { 
    free(qk[i]); 
  } 
  free(qk); 
  free(pk1);
  if (have_pk2) {
    free(pk2);
  }


#ifdef PERFORM_STATISTICS
  stop = clock();
  timing.rad_qk += stop-start;
#endif
  
  return type1;
} 
	  
double InterpolatePk(double te, int type, double *pk) {
  double *x;
  double r;
  int np, nt;

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
  np = 3;
  nt = 1;
  uvip3p_(&np, &n_tegrid, x, pk, &nt, &te, &r);  
  return r;
}

int CollisionStrength(double *s, double *e, int lower, int upper, int msub) {
  int i, j, t, h, p, m;  
  LEVEL *lev1, *lev2;
  double te, c, r, s3j;
  ANGULAR_ZMIX *ang;
  int nz, j1, j2;
  int nq, q[MAXMSUB];
  int ie, type, tt, np;
  double rqk_tmp[MAXMSUB];
  double rq[MAXMSUB][MAXNE], *rqk;
  double *x0, *x1;

  lev1 = GetLevel(lower);
  if (lev1 == NULL) return -1;
  lev2 = GetLevel(upper);
  if (lev2 == NULL) return -1;
  te = lev2->energy - lev1->energy;
  if (te <= 0) return -1;
  *e = te;

  if (interpolate_egrid) {
    x0 = log_egrid;
    x1 = log_usr;
    if (egrid_type == 1) {
      for (i = 0; i < n_egrid; i++) {
	x0[i] = log(te + egrid[i]);
      }
    }
    if (usr_egrid_type == 1) {
      for (i = 0; i < n_usr; i++) {
	x1[i] = log(te + usr_egrid[i]);
      }
    }
  }

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
	for (ie = 0; ie < n_egrid; ie++) {
	  rq[p][ie] = 0.0;
	}
	p++;
      }
    }
  } else {
    rqk = rq[0];
    for (ie = 0; ie < n_egrid; ie++) {
      rqk[ie] = 0.0;
    }
  }


  nz = AngularZMix(&ang, lower, upper, -1, -1);
  type = -1;
  for (i = 0; i < nz; i++) {
    for (j = i; j < nz; j++) {
      c = ang[i].coeff * ang[j].coeff;
      if (i != j) c *= 2.0;
      if (!msub) {
	if (ang[i].k != ang[j].k) continue;
	for (ie = 0; ie < n_egrid; ie++) {
	  tt = CERadialQk(&r, ie, te, ang[i].k0, ang[i].k1, 
			    ang[j].k0, ang[j].k1, ang[i].k); 
	  r /= (ang[i].k + 1.0);
	  rqk[ie] += c*r;
	}
	if (type != 1 && tt >= 0) type = tt;
      } else {
	for (ie = 0; ie < n_egrid; ie++) {
	  tt = CERadialQkMSub(rqk_tmp, ie, te, ang[i].k0, ang[i].k1,
			      ang[j].k0, ang[j].k1,
			      ang[i].k, ang[j].k, nq, q);
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
	      rq[p][ie] += c*rqk_tmp[m]*s3j;
	      p++;
	    }
	  }
	}
	if (type != 1 && tt >= 0) type = tt;
      }
    }
  }
    
  if (nz > 0) {
    free(ang);
  }

  /* there is a factor of 4 coming from normalization and the 2 
     from the formula */
  np = 3;
  if (!msub) {
    if (interpolate_egrid) {
      if (type < 0) {
	for (ie = 0; ie < n_egrid; ie++) {
	  rqk[ie] = log(rqk[ie]);
	}
      }
      uvip3p_(&np, &n_egrid, x0, rqk, &n_usr, x1, s);
      if (type < 0) {
	for (ie = 0; ie < n_usr; ie++) {
	  s[ie] = exp(s[ie]);
	}
      }
      for (ie = 0; ie < n_usr; ie++) {
	s[ie] *= 8.0;
      }
    } else {
      for (ie = 0; ie < n_egrid; ie++) {
	s[ie] = rqk[ie]*8.0;
      }
    }
    return 1;
  } else {
    p = 0; 
    for (t = -j1; t <= 0; t += 2) {
      for (h = -j2; h <= j2; h += 2) {	
	if (interpolate_egrid) {
	  if (type < 0) {
	    for (ie = 0; ie < n_egrid; ie++) {
	      rq[p][ie] = log(rq[p][ie]);
	    }
	  }
	  uvip3p_(&np, &n_egrid, x0, rq[p], &n_usr, x1, s);	
	  if (type < 0) {
	    for (ie = 0; ie < n_usr; ie++) {
	      s[ie] = exp(s[ie]);
	    }
	  }
	  for (ie = 0; ie < n_usr; ie++) {
	    (*s) *= 8.0;
	    s++;
	  }
	  p++;
	} else {
	  for (ie = 0; ie < n_egrid; ie++) {
	    (*s) = rq[p][ie]*8.0;
	    s++;
	  }
	  p++;
	}
      } 
    }
    return (p);
  }
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
  double s[MAXMSUB*MAXNUSR];
  int *alev;
  LEVEL *lev1, *lev2;
  double emin, emax, e, c, b, e0;

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
    printf("No transitions can occur\n");
    return 0;
  }
  if (tegrid[0] < 0.0) {
    e = 2.0*(emax-emin)/(emax+emin);    
    if (e < 0.1) {
      SetCETEGrid(1, 0.5*(emin+emax), emax);
    } else if (e < 0.3) {
      SetCETEGrid(2, emin, emax);
    } else {
      if (m == 2) n_tegrid = 2;
      SetCETEGrid(n_tegrid, emin, emax);
    }
  }

  e = 0.5*(emin+emax);
  emin = 0.1*e;
  emax = 8.0*e;
  if (n_usr == 0) {
    n_usr = 6;
  }
  interpolate_egrid = 1;
  if (msub) {
    pw_type = 1;
    if (egrid_type < 0) egrid_type = 0;
    if (usr_egrid_type < 0) usr_egrid_type = 0;
  } else {
    if (pw_type < 0) pw_type = 0;
    if (egrid_type < 0) egrid_type = 1;
    if (usr_egrid_type < 0) usr_egrid_type = 1;
  }
    
  if (usr_egrid[0] < 0.0) {
    if (n_egrid > n_usr) {
      SetUsrCEEGridDetail(n_egrid, egrid);
      if (egrid_type == 0 && usr_egrid_type == 1) {
	for (i = 0; i < n_egrid; i++) {
	  usr_egrid[i] -= e;
	}
      } else if (egrid_type == 1 && usr_egrid_type == 0) {
	for (i = 0; i < n_egrid; i++) {
	  usr_egrid[i] += e;
	  log_usr[i] = log(usr_egrid[i]);
	}
      }
      interpolate_egrid = 0;
    } else {
      if (usr_egrid_type == 0) {
	SetUsrCEEGrid(n_usr, emin, emax, -e);
      } else {
	SetUsrCEEGrid(n_usr, emin, emax, e);
      }
    }
  }
  if (n_egrid == 0) {
    n_egrid = 6;
  }
  if (egrid[0] < 0.0) {
    if (n_usr <= 6) {
      SetCEEGridDetail(n_usr, usr_egrid);
      if (egrid_type == 0 && usr_egrid_type == 1) {
	for (i = 0; i < n_egrid; i++) {
	  egrid[i] += e;
	  log_egrid[i] = log(egrid[i]);
	}
      } else if (egrid_type == 1 && usr_egrid_type == 0) {
	for (i = 0; i < n_egrid; i++) {
	  egrid[i] -= e;
	}
      }
      interpolate_egrid = 0;
    } else {
      emin = usr_egrid[0];
      emax = usr_egrid[n_usr-1];
      if (usr_egrid_type == 0) {
	emin -= e;
	emax -= e;
      }
      if (egrid_type == 0) {
	SetCEEGrid(n_egrid, emin, emax, -e);
      } else {
	SetCEEGrid(n_egrid, emin, emax, e);
      }
    }
  }

  if (pw_scratch.nkl == 0) {
    SetCEPWGrid(0, NULL, NULL);
  }
  
  e = 0.0;
  c = GetResidualZ();
  PrepCoulombBethe(1, n_tegrid, n_egrid, c, &e, tegrid, egrid,
		   pw_scratch.nkl, pw_scratch.kl, egrid_type, 
		   pw_type, msub);
 
  fprintf(f, " TEGRID:   ");
  for (i = 0; i < n_tegrid; i++) {
    fprintf(f, "%10.4E ", tegrid[i]*HARTREE_EV);
  }
  fprintf(f, "\n");

  if (egrid_type == 0) fprintf(f, " Incident Electron ");
  else fprintf(f, " Scattered Electron ");
  fprintf(f, "EGRID:    ");
  for (i = 0; i < n_egrid; i++) {
    fprintf(f, "%10.4E ", egrid[i]*HARTREE_EV);
  }
  fprintf(f, "\n");
  if (usr_egrid_type == 0) fprintf(f, " Incident Electron UsrEGrid\n\n");
  else fprintf(f, " Scattered Electron UsrEGrid\n\n");

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
	if (usr_egrid_type == 0) {
	  e0 = usr_egrid[ie];
	} else {
	  e0 = usr_egrid[ie] + e;
	}
	b = 1.0+0.5*FINE_STRUCTURE_CONST2*e0;
	for (m = 0; m < k; m++) {
	  c = s[m*n_usr+ie];
	  fprintf(f, "%-10.3E", c);
	  c *= PI * AREA_AU20/(2*e0*b*(j1+1));
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
  n_usr = 0;
  egrid[0] = -1.0;
  usr_egrid[0] = -1.0;
  tegrid[0] = -1.0;  
}

  
  




