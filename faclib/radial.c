#include "radial.h"
#include "cf77.h"

static char *rcsid="$Id$";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

static POTENTIAL *potential;
#define Large(orb) ((orb)->wfun)
#define Small(orb) ((orb)->wfun + potential->maxrp)

static ARRAY *orbitals;
static int n_orbitals;
static int n_continua;
 
static double _dwork[MAXRP];
static double _dwork1[MAXRP];
static double _dwork2[MAXRP];
static double _dwork3[MAXRP];
static double _dwork4[MAXRP];
static double _dwork5[MAXRP];
static double _dwork6[MAXRP];
static double _dwork7[MAXRP];
static double _dwork8[MAXRP];
static double _dwork9[MAXRP];
static double _dwork10[MAXRP];
static double _dwork11[MAXRP];
static double _phase[MAXRP];
static double _dphase[MAXRP];
static double _dphasep[MAXRP];
static double _yk[MAXRP];
static double _zk[MAXRP];
static double _xk[MAXRP];

static struct {
  double stabilizer;
  double tolerance; /* tolerance for self-consistency */
  int maxiter; /* max iter. for self-consistency */
  double screened_charge; 
  int screened_kl;
  int n_screen;
  int *screened_n;
  int iprint; /* printing infomation in each iteration. */
  int iset;
} optimize_control = {OPTSTABLE, OPTTOL, OPTNITER, 
		      1.0, 1, 0, NULL, OPTPRINT, 0};
static struct {
  int kl0;
  int kl1;
} slater_cut = {100000, 1000000};

static struct {
  int se;
  int vp;
  int nms;
  int sms;
  int br;
} qed = {QEDSE, QEDVP, QEDNMS, QEDSMS, QEDBREIT};

static AVERAGE_CONFIG average_config = {0, 0, NULL, NULL, NULL};
 
static MULTI *slater_array;
static MULTI *breit_array;
static MULTI *vinti_array;
static MULTI *qed1e_array;
static MULTI *residual_array;
static MULTI *multipole_array; 
static MULTI *moments_array;
static MULTI *gos_array;
static MULTI *yk_array;

static int n_awgrid = 0;
static double awgrid[MAXNTE];
static double aw2grid[MAXNTE];

static double PhaseRDependent(double x, double eta, double b);

#ifdef PERFORM_STATISTICS
static RAD_TIMING rad_timing = {0, 0, 0, 0};
int GetRadTiming(RAD_TIMING *t) {
  memcpy(t, &rad_timing, sizeof(RAD_TIMING));
  return 0;
}
#endif

static void InitOrbitalData(void *p, int n) {
  ORBITAL *d;
  int i;
  
  d = (ORBITAL *) p;
  for (i = 0; i < n; i++) {
    d[i].wfun = NULL;
    d[i].phase = NULL;
    d[i].ilast = -1;
  }
}

static void InitYkData(void *p, int n) {
  SLATER_YK *d;
  int i;

  d = (SLATER_YK *) p;
  for (i = 0; i < n; i++) {
    d[i].npts = -1;
    d[i].yk = NULL;
  }
}

int FreeSimpleArray(MULTI *ma) {
  MultiFreeData(ma, NULL);
  return 0;
}

int FreeSlaterArray(void) {
  return FreeSimpleArray(slater_array);
}

int FreeResidualArray(void) {
  return FreeSimpleArray(residual_array);
}

static void FreeMultipole(void *p) {
  double *dp;
  dp = *((double **) p);
  free(dp);
  *((double **) p) = NULL;
}

static void FreeYkData(void *p) {
  SLATER_YK *dp;
  
  dp = (SLATER_YK *) p;
  if (dp->npts > 0) free(dp->yk);
}

int FreeMultipoleArray(void) {
  MultiFreeData(multipole_array, FreeMultipole);
  return 0;
}

int FreeMomentsArray(void) {
  MultiFreeData(moments_array, NULL);
  return 0;
}

int FreeGOSArray(void) {
  MultiFreeData(gos_array, FreeMultipole);
  return 0;
}

int FreeYkArray(void) {
  MultiFreeData(yk_array, FreeYkData);
  return 0;
}
  
double *WLarge(ORBITAL *orb) {
  return Large(orb);
}

double *WSmall(ORBITAL *orb) {
  return Small(orb);
}

POTENTIAL *RadialPotential(void) {
  return potential;
}

void SetSlaterCut(int k0, int k1) {
  if (k0 > 0) {
    slater_cut.kl0 = 2*k0;
  } else {
    slater_cut.kl0 = 1000000;
  } 
  if (k1 > 0) {
    slater_cut.kl1 = 2*k1;
  } else {
    slater_cut.kl1 = 1000000;
  }
}

int GetBoundary(double *rb, double *b, int *nmax, double *dr) {
  int ib;

  if (potential->ib > 0) {
    ib = potential->ib;
  } else {
    ib = potential->maxrp-1;
  }
  *rb = potential->rad[ib];
  *b = potential->bqp;
  *nmax = potential->nb;
  *dr = potential->rad[ib]-potential->rad[ib-1];
  return potential->ib;
}

int SetBoundary(int nmax, double p, double bqp) {
  ORBITAL *orb;
  int i, j, n, kl, kl2, kappa, k;
  double d1, d2, d;

  if (nmax == -100) {
    if (p <= 0.0) {
      printf("2nd argument must be > 0 in SetBoundary(-100,...)\n");
      return -1;
    }
    for (i = 0; i < potential->maxrp-10; i++) {
      if (potential->rad[i] >= p) break;
    }
    potential->ib1 = i;
    if (bqp > 0) {
      if (bqp >= p) {
	printf("3rd argument must be less than 2nd in SetBoundary(-100,...)\n");
	return -1;
      }
      for (i = 0; i < potential->maxrp; i++) {
	if (potential->rad[i] >= bqp) break;
      }
      if (i < potential->ib1) {
	potential->ib = i;
      }
    }
    return 0;
  }
  potential->nb = abs(nmax);
  potential->bqp = bqp;
  if (nmax == 0) {
    potential->ib = 0;
  } else if (nmax < 0) {
    d = GetResidualZ();
    d1 = potential->nb;
    if (p < 0.0) p = 2.0*d1*d1/d;
    for (i = 0; i < potential->maxrp; i++) {
      if (potential->rad[i] >= p) break;
    }
    if (i > potential->maxrp-10) {
      printf("enlarge maxrp\n");
      exit(1);
    }
    if (IsEven(i)) i++;
    potential->ib = i;
    for (n = 1; n <= potential->nb; n++) {
      for (kl = 0; kl < n; kl++) {
	kl2 = 2*kl;
	for (j = kl2 - 1; j <= kl2 + 1; j += 2) {
	  if (j < 0) continue;
	  kappa = GetKappaFromJL(j, kl2);
	  k = OrbitalIndex(n, kappa, 0);
	  orb = GetOrbitalSolved(k);
	}
      }
    }
  } else {
    if (p < 0.0) p = 1E-8;
    for (n = 1; n <= potential->nb; n++) {
      for (kl = 0; kl < n; kl++) {
	kl2 = 2*kl;
	for (j = kl2 - 1; j <= kl2 + 1; j += 2) {
	  if (j < 0) continue;
	  kappa = GetKappaFromJL(j, kl2);
	  k = OrbitalIndex(n, kappa, 0);
	  orb = GetOrbitalSolved(k);	  
	  for (i = orb->ilast-3; i >= 0; i--) {
	    d1 = Large(orb)[i];
	    d2 = Small(orb)[i];
	    d = d1*d1 + d2*d2;
	    if (d >= p) {
	      i++;
	      break;
	    }
	  }
	  if (potential->ib < i) potential->ib = i;
	}
      }
    }
    if (IsEven(i)) i++;
    if (i > potential->maxrp-10) {
      printf("enlarge maxrp\n");
      return -1;
    }
  }
  potential->ib1 = potential->ib;
  return 0;
}

int RadialOverlaps(char *fn, int kappa) {
  ORBITAL *orb1, *orb2;
  int i, j, k;
  double r;
  FILE *f;

  f = fopen(fn, "w");
  for (k = 0; k < potential->maxrp; k++) {
    _yk[k] = 1.0;
  }
  for (i = 0; i < n_orbitals; i++) {
    orb1 = GetOrbital(i);
    if (orb1->kappa != kappa) continue;
    for (j = 0; j <= i; j++) {
      orb2 = GetOrbital(j);
      if (orb2->kappa != kappa) continue;
      Integrate(_yk, orb1, orb2, 1, &r, 0);
      fprintf(f, "%2d %2d %10.3E  %2d %2d %10.3E  %12.5E\n", 
	      orb1->n, orb1->kappa, orb1->energy, 
	      orb2->n, orb2->kappa, orb2->energy, r);
    }
  }
  fclose(f);

  return 0;
}
  
void SetSE(int n) {
  qed.se = n;
}

void SetVP(int n) {
  qed.vp = n;
}

void SetBreit(int n) {
  qed.br = n;
}

void SetMS(int nms, int sms) {
  qed.nms = nms;
  qed.sms = sms;
}

int SetAWGrid(int n, double awmin, double awmax) {
  int i;
  if (awmin < 1E-3) {
    awmin = 1E-3;
    awmax = awmax + 1E-3;
  }
  awmin *= awmin;
  awmax *= awmax;
  n_awgrid = SetTEGrid(aw2grid, NULL, n, awmin, awmax);
  for (i = 0; i < n_awgrid; i++) {
    awgrid[i] = sqrt(aw2grid[i]);
  }
  return 0;
}

void SetOptimizeMaxIter(int m) {
  optimize_control.maxiter = m;
}

void SetOptimizeStabilizer(double m) {
  if (m < 0) {
    optimize_control.iset = 0;
  } else {
    optimize_control.iset = 1;
    optimize_control.stabilizer = m;
  }
}

void SetOptimizeTolerance(double c) {
  optimize_control.tolerance = c;
}

void SetOptimizePrint(int m) {
  optimize_control.iprint = m;
}

void SetOptimizeControl(double tolerance, double stabilizer, 
			int maxiter, int iprint) {
  optimize_control.maxiter = maxiter;
  optimize_control.stabilizer = stabilizer;
  optimize_control.tolerance = tolerance;
  optimize_control.iprint = iprint;  
  optimize_control.iset = 1;
}

void SetScreening(int n_screen, int *screened_n, 
		  double screened_charge, int kl) {
  optimize_control.screened_n = screened_n;
  optimize_control.screened_charge = screened_charge;
  optimize_control.n_screen = n_screen;
  optimize_control.screened_kl = kl;
}

int SetRadialGrid(int maxrp, double ratio, double asymp, double rmin) {
  if (maxrp > MAXRP) {
    printf("MAXRP must be <= %d\n", MAXRP);
    printf("to enlarge the limit, change MAXRP in global.h\n");
    return -1;
  }
  if (maxrp < 0) maxrp = DMAXRP;
  potential->maxrp = maxrp;
  if (asymp < 0 && ratio < 0) {
    asymp = GRIDASYMP;
    ratio = GRIDRATIO;
  }
  if (rmin <= 0) rmin = GRIDRMIN;
  potential->rmin = rmin;
  if (ratio == 0) potential->ratio = GRIDRATIO;
  else potential->ratio = ratio;
  if (asymp == 0) potential->asymp = GRIDASYMP;
  else potential->asymp = asymp;
  potential->flag = 0;
  return 0;
}

void AdjustScreeningParams(double *u) {
  int i;
  double c;
  
  c = 0.5*u[potential->maxrp-1];
  for (i = 0; i < potential->maxrp; i++) {
    if (u[i] > c) break;
  }
  potential->lambda = log(2.0)/potential->rad[i];
}
   
double SetPotential(AVERAGE_CONFIG *acfg, int iter) {
  int i, j, k1, k2, k, t, m, j1, j2, kl1, kl2;
  ORBITAL *orb1, *orb2;
  double large1, small1, large2, small2;
  int norbs, kmin, kmax, jmax, kmax0 = 0;
  double *u, *w, *v, w3j, a, b, c, r;

  u = potential->U;
  w = potential->W;
  v = _dwork2;

  for (j = 0; j < potential->maxrp; j++) {
    w[j] = 0.0;
  }

  norbs = 0;
  jmax = 0;
  for (i = 0; i < acfg->n_shells; i++) {
    k1 = OrbitalExists(acfg->n[i], acfg->kappa[i], 0.0);
    if (k1 < 0) continue;
    orb1 = GetOrbital(k1);
    if (orb1->wfun == NULL) continue;
    for (j = 0; j <= orb1->ilast; j++) {
      large1 = Large(orb1)[j];
      small1 = Small(orb1)[j];
      w[j] += (large1*large1 + small1*small1)*acfg->nq[i];
    }
    if (jmax < orb1->ilast) jmax = orb1->ilast;
    norbs++;
  }
  for (j = 0; j < potential->maxrp; j++) {
    u[j] = 0.0;
  }
  if (norbs && potential->N > 1+EPS3) {
    for (i = 0; i < acfg->n_shells; i++) {
      k1 = OrbitalExists(acfg->n[i], acfg->kappa[i], 0.0);
      if (k1 < 0) continue;
      orb1 = GetOrbital(k1);
      if (orb1->wfun == NULL) continue;
      GetJLFromKappa(acfg->kappa[i], &j1, &kl1);
      kmin = 0;
      kmax = 2*j1;
      kmax = Min(kmax, kmax0);
      for (k = kmin; k <= kmax; k += 2) {
	t = k/2;
	if (IsOdd(t)) continue;
	GetYk(t, _yk, orb1, orb1, k1, k1, -1);
	if (t > 0) {
	  w3j = W3j(j1, k, j1, -1, 0, -1);
	  w3j *= w3j*(j1+1.0)/j1;
	}
	for (m = 1; m <= jmax; m++) {
	  large1 = Large(orb1)[m];
	  small1 = Small(orb1)[m];
	  b = large1*large1 + small1*small1;
	  if (t == 0) {	  
	    u[m] += acfg->nq[i]*_yk[m];
	    a = _yk[m]*b*acfg->nq[i];
	    a /= w[m];
	    u[m] -= a;
	  } else {
	    a = acfg->nq[i]*(acfg->nq[i]-1.0);
	    a *= w3j*_yk[m]*b;
	    a /= w[m];
	    u[m] -= a;
	  }
	}
      }
      if (iter < 3) continue;
      for (j = 0; j < i; j++) {
	k2 = OrbitalExists(acfg->n[j], acfg->kappa[j], 0.0);
	if (k2 < 0) continue;
	orb2 = GetOrbital(k2);
	if (orb2->wfun == NULL) continue;
	GetJLFromKappa(acfg->kappa[j], &j2, &kl2);
	kmin = abs(j1 - j2);
	kmax = j1 + j2;
	kmax = Min(kmax, kmax0);
	if (IsOdd(kmin)) kmin++;
	for (k = kmin; k <= kmax; k += 2) {
	  if (IsOdd((k+kl1+kl2)/2)) continue;
	  t = k/2;
	  GetYk(t, _yk, orb1, orb2, k1, k2, -1);
	  w3j = W3j(j1, k, j2, -1, 0, -1);
	  w3j *= w3j;
	  for (m = 1; m <= jmax; m++) {
	    large1 = Large(orb1)[m];
	    large2 = Large(orb2)[m];
	    small1 = Small(orb1)[m];
	    small2 = Small(orb2)[m];
	    a = acfg->nq[i]*acfg->nq[j];
	    a *= w3j*_yk[m]*(large1*large2+small1*small2);
	    a /= w[m];
	    u[m] -= a;
	  }
	}
      }
    }

    u[0] = u[1];
   
    jmax = jmax - 8;
    for (j = jmax+1; j < potential->maxrp; j++) {
      u[j] = u[jmax];
    }
    if (potential->N > 1) {
      for (j = jmax; j > 0; j--) {
        if (fabs(u[j]-potential->N + 1.0) > EPS16) break;
      }
      potential->r_core = j+1;
    }

    if (iter < 3) {
      r = 1.0;
      for (j = 0; j < potential->maxrp; j++) {
	v[j] = u[j];
      }
    } else {	
      r = 0.0;
      k = 0;
      a = optimize_control.stabilizer;
      b = 1.0 - a;
      for (j = 0; j < potential->maxrp; j++) {
	if (u[j] + 1.0 != 1.0) {
	  r += fabs(1.0 - v[j]/u[j]);
	  k++;
	}
	u[j] = b*v[j] + a*u[j];
	v[j] = u[j];
      }
      r /= k;
    }
    AdjustScreeningParams(u);
    SetPotentialVc(potential);
    for (j = 0; j < potential->maxrp; j++) {
      a = u[j] - potential->Z[j];
      b = potential->Vc[j]*potential->rad[j];
      u[j] = a - b;
      u[j] /= potential->rad[j];
    }
    SetPotentialU(potential, 0, NULL);
  } else {
    if (potential->N < 1.0+EPS3) {
      SetPotentialVc(potential);
      SetPotentialU(potential, -1, NULL);
      return 0.0;
    }
    r = potential->Z[potential->maxrp-1];
    b = (1.0 - 1.0/potential->N);
    for (i = 0; i < acfg->n_shells; i++) {
      a = acfg->nq[i];
      c = acfg->n[i];
      c = r/(c*c);
      for (j = 0; j < potential->maxrp; j++) {
	u[j] += a*b*(1.0 - exp(-c*potential->rad[j]));
      }
    }
    AdjustScreeningParams(u);
    SetPotentialVc(potential);
    for (j = 0; j < potential->maxrp; j++) {
      a = u[j] - potential->Z[j];
      b = potential->Vc[j]*potential->rad[j];
      u[j] = a - b;
      u[j] /= potential->rad[j];
    }
    SetPotentialU(potential, 0, NULL);
    r = 1.0;
  }
  
  return r;
}

int GetPotential(char *s) {
  AVERAGE_CONFIG *acfg;
  ORBITAL *orb1;
  double large1, small1;
  int norbs, jmax;  
  FILE *f;
  int i, j, k, k1;
  double *w, *v, *ve0, *ve1;

  /* get the average configuration for the groups */
  acfg = &(average_config);

  w = potential->W;
  v = potential->dW;
  ve0 = _dphase;
  ve1 = _dphasep;

  f = fopen(s, "w");
  if (!f) return -1;
  
  fprintf(f, "# Lambda = %10.3E\n", potential->lambda);
  fprintf(f, "#      A = %10.3E\n", potential->a);
  fprintf(f, "#     ar = %10.3E\n", potential->ar);
  fprintf(f, "#     br = %10.3E\n", potential->br);
  fprintf(f, "#     rb = %10.3E\n", potential->rad[potential->ib]);
  fprintf(f, "#    rb1 = %10.3E\n", potential->rad[potential->ib1]);
  fprintf(f, "#    bqp = %10.3E\n", potential->bqp);
  fprintf(f, "#     nb = %d\n", potential->nb);
  
  for (j = 0; j < potential->maxrp; j++) {
    w[j] = 0.0;
    v[j] = -potential->Z[j]/potential->rad[j];
    ve0[j] = 0.0;
    ve1[j] = 0.0;
  }

  norbs = 0;
  jmax = 0;
  for (i = 0; i < acfg->n_shells; i++) {
    k1 = OrbitalExists(acfg->n[i], acfg->kappa[i], 0.0);
    if (k1 < 0) continue;
    orb1 = GetOrbital(k1);
    if (orb1->wfun == NULL) SolveDirac(orb1);
    for (j = 0; j <= orb1->ilast; j++) {
      large1 = Large(orb1)[j];
      small1 = Small(orb1)[j];
      w[j] += (large1*large1 + small1*small1)*acfg->nq[i];
    }
    GetYk(0, _yk, orb1, orb1, k1, k1, -1);
    for (k = 0; k < potential->maxrp; k++) {
      v[k] += _yk[k]*acfg->nq[i]/potential->rad[k];
    }
    if (jmax < orb1->ilast) jmax = orb1->ilast;
    norbs++;
  }
  
  for (k = 0; k < potential->maxrp; k++) {
    w[k] = w[k]/(potential->rad[k]*potential->rad[k]);
    w[k] = - pow(w[k], 1.0/3);
    ve1[k] += w[k]*0.4235655;
    ve0[k] = potential->Vc[k]+potential->U[k] - v[k];
  }
  ve1[0] = ve1[1];

  fprintf(f, "# Mean configuration:\n");
  for (i = 0; i < acfg->n_shells; i++) {
    fprintf(f, "# %2d %2d\t%10.3E\n", acfg->n[i], acfg->kappa[i], acfg->nq[i]);
  }
  fprintf(f, "\n\n");
  for (i = 0; i < potential->maxrp; i++) {
    fprintf(f, "%5d %14.8E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E\n",
	    i, potential->rad[i], potential->Z[i], 
	    potential->Vc[i]+potential->U[i], v[i], 
	    ve0[i], ve1[i], potential->uehling[i]);
  }

  fclose(f);  
  
  return 0;
}

double GetResidualZ(void) {
  double z;
  z = potential->Z[potential->maxrp-1];
  if (potential->N > 0) z -= potential->N - 1;
  return z;
}

double GetRMax(void) {
  return potential->rad[potential->maxrp-10];
}

int SetAverageConfig(int nshells, int *n, int *kappa, double *nq) {
  int i;
  if (nshells <= 0) return -1;
  if (average_config.n_shells > 0) {
    average_config.kappa = (int *) realloc(average_config.kappa, 
					   sizeof(int)*nshells);
    average_config.nq = (double *) realloc(average_config.nq, 
					   sizeof(double)*nshells);
    average_config.n = (int *) realloc(average_config.n, 
				       sizeof(int)*nshells);
  } else {
    average_config.kappa = (int *) malloc(sizeof(int)*nshells);
    average_config.nq = (double *) malloc(sizeof(double)*nshells);
    average_config.n = (int *) malloc(sizeof(int)*nshells);
  }
  for (i = 0; i < nshells; i++) {
    average_config.n[i] = n[i];
    average_config.kappa[i] = kappa[i];
    average_config.nq[i] = nq[i];
  }
  average_config.n_shells = nshells;
  average_config.n_cfgs = 1;
  return 0;
}
    
int OptimizeRadial(int ng, int *kg, double *weight) {
  AVERAGE_CONFIG *acfg;
  double tol;
  ORBITAL orb_old, *orb;
  int i, k, no_old;
  double a, b, z;
  int iter;
  int *frozen;

  /* get the average configuration for the groups */
  acfg = &(average_config);
  if (ng > 0) {
    if (ng > 1) {
      printf("\nWarning: more than 1 configuration groups");
      printf(" are used in OptimizeRadial.\n");
      printf("It is usually best to use the lowest lying configuration group.\n");
      printf("Make sure that you know what you are doing.\n\n");
    }
    if (acfg->n_shells > 0) {      
      acfg->n_cfgs = 0;
      acfg->n_shells = 0;
      free(acfg->n);
      free(acfg->kappa);
      free(acfg->nq);
      acfg->n = NULL;
      acfg->nq = NULL;
      acfg->kappa = NULL;
    }
    GetAverageConfig(ng, kg, weight, 
		     optimize_control.n_screen,
		     optimize_control.screened_n,
		     optimize_control.screened_charge,
		     optimize_control.screened_kl, acfg); 
  } else {
    if (acfg->n_shells <= 0) {
      printf("No average configuation exist. \n");
      printf("Specify with AvgConfig, ");
      printf("or give config groups to OptimizeRadial.\n");
      return -1;
    }
  }

  a = 0;
  for (i = 0; i < acfg->n_shells; i++) {
    if (optimize_control.iprint) 
      printf("%d %d %f\n", acfg->n[i], acfg->kappa[i], acfg->nq[i]);
    a += acfg->nq[i];
  }
  potential->N = a;  

  /* setup the radial grid if not yet */
  if (potential->flag == 0) {
    SetOrbitalRGrid(potential);
  }

  SetPotentialZ(potential, 0.0);
  z = potential->Z[potential->maxrp-1];
  if (a > 0.0) z = z - a + 1;
  potential->a = 0.0;
  potential->lambda = 0.5*z;
  if (potential->N > 1) {
    potential->r_core = potential->maxrp-5;
  } else {
    potential->r_core = 50;
  }

  if (optimize_control.iset == 0) {
    optimize_control.stabilizer = 0.25 + 0.75*(z/potential->Z[potential->maxrp-1]);
  }

  frozen = (int *) malloc(acfg->n_shells*sizeof(int));
  for (i = 0; i < acfg->n_shells; i++) {
    frozen[i] = 0;
  }

  no_old = 0;
  iter = 0;
  SetPotentialZ(potential, 0.0);
  tol = 1.0; 
  while (tol > optimize_control.tolerance) {
    if (iter > optimize_control.maxiter) break;
    a = SetPotential(acfg, iter);
    FreeYkArray();
    tol = 0.0;
    for (i = 0; i < acfg->n_shells; i++) {
      k = OrbitalExists(acfg->n[i], acfg->kappa[i], 0.0);
      if (k < 0) {
	frozen[i] = 0;
	orb_old.energy = 0.0;
	orb = GetNewOrbital();
	orb->kappa = acfg->kappa[i];
	orb->n = acfg->n[i];
	orb->energy = 1.0;
	no_old = 1;	
      } else {
	orb = GetOrbital(k);
	if (orb->wfun == NULL) {
	  frozen[i] = 0;
	  orb_old.energy = 0.0;
	  orb->energy = 1.0;
	  orb->kappa = acfg->kappa[i];
	  orb->n = acfg->n[i];
	  no_old = 1;	
	} else {
	  if (!frozen[i]) {
	    orb_old.energy = orb->energy; 
	    if (orb->wfun) free(orb->wfun);
	    no_old = 0;
	  } else {
	    continue;
	  }
	}
      }
      if (SolveDirac(orb) < 0) {
	return -1;
      }
      
      if (no_old) { 
	tol = 1.0;
	continue;
      } 
      b = fabs(1.0 - orb_old.energy/orb->energy);
      if (tol < b) tol = b;
    }
    if (optimize_control.iprint) {
      printf("%4d %13.5E %13.5E\n", iter, tol, a);
    }
    if (tol < a) tol = a;
    iter++;
  }
  
  free(frozen);
  if (potential->uehling[0] == 0.0) {
    SetPotentialUehling(potential, qed.vp);
  }
  if (iter > optimize_control.maxiter) {
    printf("Maximum iteration reached in OptimizeRadial\n");
    return 1;
  }

  return 0;
}      

static double EnergyFunc(int *n, double *x) {
  double a;

  if (x[1] < -EPS10) return 0.0;
  if (x[0] <= 0.0) return 0.0;

  potential->lambda = x[0];
  potential->a = x[1];
  SetPotentialVc(potential);
  ReinitRadial(1);
  ClearOrbitalTable(0);
  a = AverageEnergyAvgConfig(&average_config);

  return a;
}

int RefineRadial(int maxfun, int msglvl) {
  int n, ierr, mode, nfe, lw[4];
  double xtol, scale[2];
  double f0, f, x[2];
  
  if (maxfun <= 0) maxfun = 250;
  xtol = EPS3;
  n = 2;
  mode = 0;
  x[0] = potential->lambda;
  x[1] = potential->a;
  scale[0] = 0.01;
  scale[1] = 0.01;
  
  f0 = EnergyFunc(&n, x);
  if (msglvl > 0) {
    printf("%10.3E %10.3E %15.8E\n", x[0], x[1], f0);
  }
  nfe = 0;
  ierr = 0;
  SUBPLX(EnergyFunc, n, xtol, maxfun, mode, scale, x,
	 &f, &nfe, _dwork11, lw, &ierr);
  if (msglvl > 0) {
    printf("%10.3E %10.3E %15.8E %d\n", x[0], x[1], f, nfe);
  }
  f = EnergyFunc(&n, x);
  if (ierr) {
    if (f > f0) {
      printf("Error in RefineRadial: %d\n", ierr);
      return ierr;
    } else if (msglvl > 0) {
      printf("Warning in RefineRadial: %d\n", ierr);
    }
  }
  
  return 0;
}

int SolveDirac(ORBITAL *orb) {
  int err;
#ifdef PERFORM_STATISTICS
  clock_t start, stop;
  start = clock();
#endif

  err = 0;  
  potential->flag = -1;
  err = RadialSolver(orb, potential);
  if (err) { 
    printf("Error ocuured in RadialSolver, %d\n", err);
    printf("%d %d %10.3E\n", orb->n, orb->kappa, orb->energy);
    exit(1);
  }
#ifdef PERFORM_STATISTICS
  stop = clock();
  rad_timing.dirac += stop - start;
#endif

  return err;
}

int WaveFuncTable(char *s, int n, int kappa, double e) {
  int i, k;
  FILE *f;
  ORBITAL *orb;
  double z, a, ke, y;

  e /= HARTREE_EV;
  k = OrbitalIndex(n, kappa, e);
  if (k < 0) return -1;
  f = fopen(s, "w");
  if (!f) return -1;
  
  orb = GetOrbitalSolved(k);
  
  fprintf(f, "#      n = %2d\n", n);
  fprintf(f, "#  kappa = %2d\n", kappa);
  fprintf(f, "# energy = %15.8E\n", orb->energy*HARTREE_EV);
  fprintf(f, "#     vc = %15.8E\n", MeanPotential(k, k)*HARTREE_EV);
  if (n != 0) {
    fprintf(f, "\n\n");
    if (n < 0) k = potential->ib;
    else k = 0;
    for (i = k; i <= orb->ilast; i++) {
      fprintf(f, "%-4d %14.8E %13.6E %13.6E %13.6E %13.6E\n", 
	      i, potential->rad[i], 
	      (potential->Vc[i])*potential->rad[i],
	      potential->U[i] * potential->rad[i],
	      Large(orb)[i], Small(orb)[i]); 
    }
  } else {
    a = GetPhaseShift(k);
    while(a < 0) a += TWO_PI;
    a -= (int)(a/TWO_PI);
    fprintf(f, "#  phase = %15.8E\n", a);
    fprintf(f, "\n\n");
    z = GetResidualZ();
    e = orb->energy;
    a = FINE_STRUCTURE_CONST2 * e;
    ke = sqrt(2.0*e*(1.0+0.5*a));
    y = (1.0+a)*z/ke; 
    for (i = 0; i <= orb->ilast; i++) {
      fprintf(f, "%-4d %14.8E %13.6E %13.6E %13.6E %13.6E\n", 
	      i, potential->rad[i],
	      (potential->Vc[i])*potential->rad[i],
	      potential->U[i] * potential->rad[i],
	      Large(orb)[i], Small(orb)[i]); 
    }
    for (; i < potential->maxrp; i += 2) {
      a = ke * potential->rad[i];
      a = a + y*log(2.0*a);
      a = Large(orb)[i+1] - a;
      a = a - ((int)(a/(TWO_PI)))*TWO_PI;
      if (a < 0) a += TWO_PI;
      fprintf(f, "%-4d %14.8E %13.6E %13.6E %13.6E %13.6E\n",
	      i, potential->rad[i],
	      Large(orb)[i], Large(orb)[i+1], 
	      Small(orb)[i], a);
    }
  }

  fclose(f);

  return 0;
}

double PhaseRDependent(double x, double eta, double b) {
  double tau, tau2, y, y2, t, a1, a2, sb;
  
  y = 1.0/x;
  y2 = y*y;
  tau2 = 1.0 + 2.0*eta*y - b*y2;
  tau = x*sqrt(tau2);
  
  t = eta*log(x+tau+eta) + tau - eta;
  if (b > 0.0) {
    sb = sqrt(b);
    a1 = b - eta*x;
    a2 = tau*sb;
    tau2 = atan2(a1, a2);
    tau2 -= atan2(-eta, sb);
    t += sb*tau2;
  } else if (b < 0.0) {
    b = -b;
    sb = sqrt(b);
    a1 = 2.0*(b+eta*x)/(sb*b*x);
    a2 = 2.0*tau/(b*x);
    tau2 = log(a1+a2);
    tau2 -= log(2.0*(eta/sb + 1.0)/b);
    t -= sb*tau2;
  }
  
  return t;
}

double GetPhaseShift(int k) {
  ORBITAL *orb;
  double phase1, r, y, z, ke, e, a, b1;
  int i;

  orb = GetOrbitalSolved(k);
  if (orb->n > 0) return 0.0;

  if (orb->phase) return *(orb->phase);

  z = GetResidualZ();
  e = orb->energy;
  a = FINE_STRUCTURE_CONST2 * e;
  ke = sqrt(2.0*e*(1.0 + 0.5*a));
  y = (1.0 + a)*z/ke;

  i = potential->maxrp - 1;  
  phase1 = orb->wfun[i];
  r = potential->rad[i-1];  
  b1 = orb->kappa;
  b1 = b1*(b1+1.0) - FINE_STRUCTURE_CONST2*z*z;
 
  a = ke * r;
  b1 = PhaseRDependent(a, y, b1);
  phase1 = phase1 - b1;
  
  orb->phase = malloc(sizeof(double));
  *(orb->phase) = phase1;

  return phase1;  
}

int GetNumBounds(void) {
  return n_orbitals - n_continua;
}

int GetNumOrbitals(void) {
  return n_orbitals;
}

int GetNumContinua(void) {
  return n_continua;
}

int OrbitalIndex(int n, int kappa, double energy) {
  int i, j;
  ORBITAL *orb;
  int resolve_dirac;

  resolve_dirac = 0;
  for (i = 0; i < n_orbitals; i++) {
    orb = GetOrbital(i);
    if (n == 0) {
      if (orb->n == 0 &&
	  orb->kappa == kappa && 
	  orb->energy > 0.0 &&
	  fabs(orb->energy - energy) < EPS10) {
	if (orb->wfun == NULL) {
	  if (RestoreOrbital(i) == 0) return i;
	  else {
	    resolve_dirac = 1;
	    break;
	  }
	}
	return i;
      }
    } else if (orb->n == n && orb->kappa == kappa) {
      if (orb->wfun == NULL) {
	if (RestoreOrbital(i) == 0) return i;
	else {
	  resolve_dirac = 1;
	  break;
	}
      }
      return i;
    }
  }
  
  if (!resolve_dirac) {
    orb = GetNewOrbital();
  } 

  orb->n = n;
  orb->kappa = kappa;
  orb->energy = energy;
  j = SolveDirac(orb);
  if (j < 0) {
    printf("Error occured in solving Dirac eq. err = %d\n", j);
    exit(1);
  }
  
  if (n == 0 && !resolve_dirac) {
    n_continua++;
  }
  return i;
}

int OrbitalExists(int n, int kappa, double energy) {
  int i;
  ORBITAL *orb;
  for (i = 0; i < n_orbitals; i++) {
    orb = GetOrbital(i);
    if (n == 0) {
      if (orb->kappa == kappa &&
	  fabs(orb->energy - energy) < EPS10) 
	return i;
    } else if (orb->n == n && orb->kappa == kappa) {
      return i;
    }
  }
  return -1;
}

int AddOrbital(ORBITAL *orb) {

  if (orb == NULL) return -1;

  orb = (ORBITAL *) ArrayAppend(orbitals, orb, InitOrbitalData);
  if (!orb) {
    printf("Not enough memory for orbitals array\n");
    exit(1);
  }

  if (orb->n == 0) {
    n_continua++;
  }
  n_orbitals++;
  return n_orbitals - 1;
}

ORBITAL *GetOrbital(int k) {
  return (ORBITAL *) ArrayGet(orbitals, k);
}

ORBITAL *GetOrbitalSolved(int k) {
  ORBITAL *orb;
  int i;
  
  orb = (ORBITAL *) ArrayGet(orbitals, k);
  if (orb->wfun == NULL) {
    i = SolveDirac(orb);
    if (i < 0) {
      printf("Error occured in solving Dirac eq. err = %d\n", i);
      exit(1);
    }
  }
  return orb;
}

ORBITAL *GetNewOrbital(void) {
  ORBITAL *orb;

  orb = (ORBITAL *) ArrayAppend(orbitals, NULL, InitOrbitalData);
  if (!orb) {
    printf("Not enough memory for orbitals array\n");
    exit(1);
  }

  n_orbitals++;
  return orb;
}

void FreeOrbitalData(void *p) {
  ORBITAL *orb;

  orb = (ORBITAL *) p;
  if (orb->wfun) free(orb->wfun);
  if (orb->phase) free(orb->phase);
  orb->wfun = NULL;
  orb->phase = NULL;
  orb->ilast = -1;
}

int ClearOrbitalTable(int m) {
  ORBITAL *orb;
  int i;

  if (m == 0) {
    n_orbitals = 0;
    n_continua = 0;
    ArrayFree(orbitals, FreeOrbitalData);
    SetBoundary(0, 1.0, 0.0);
  } else {
    for (i = n_orbitals-1; i >= 0; i--) {
      orb = GetOrbital(i);
      if (orb->n == 0) {
	n_continua--;
      }
      if (orb->n > 0) {
	n_orbitals = i+1;
	ArrayTrim(orbitals, i+1, FreeOrbitalData);
	break;
      }
    }
  }
  return 0;
}

int SaveOrbital(int i) {
  return 0;
}

int RestoreOrbital(int i) {
  return -1;
}

int FreeOrbital(int i) {
  ORBITAL *orb;
  orb = GetOrbital(i);
  FreeOrbitalData((void *)orb);
  return 0;
}

int SaveAllContinua(int mode) {
  int i;
  ORBITAL *orb;
  for (i = 0; i < n_orbitals; i++) {
    orb = GetOrbital(i);
    if (orb->n == 0 && orb->wfun != NULL) {
      if (SaveOrbital(i) < 0) return -1;
      if (mode) {
	FreeOrbital(i);
      }
    }
  }
  return 0;
}

int SaveContinua(double e, int mode) {
  int i;
  ORBITAL *orb;
  for (i = 0; i < n_orbitals; i++) {
    orb = GetOrbital(i);
    if (orb->n == 0 && 
	orb->wfun != NULL &&
	fabs(orb->energy - e) < EPS3) {
      if (SaveOrbital(i) < 0) return -1;
      if (mode) FreeOrbital(i);
    }
  }
  return 0;
}

int FreeAllContinua(void) {
  int i;
  ORBITAL *orb;
  for (i = 0; i < n_orbitals; i++) {
    orb = GetOrbital(i);
    if (orb->n == 0 && orb->wfun != NULL) {
      FreeOrbital(i);
    }
  }
  return 0;
}

int FreeContinua(double e) {
  int i;
  ORBITAL *orb;
  for (i = 0; i < n_orbitals; i++) {
    orb = GetOrbital(i);
    if (orb->n == 0 && 
	orb->wfun != NULL &&
	fabs(orb->energy - e) < EPS3) {
      FreeOrbital(i);
    }
  }
  return 0;
}

int ConfigEnergy(int m, int mr, int ng, int *kg) {
  CONFIG_GROUP *g;
  CONFIG *cfg;
  int k, i;

  if (m == 0) {
    if (ng == 0) {
      ng = GetNumGroups();
      for (k = 0; k < ng; k++) {
	OptimizeRadial(1, &k, NULL);
	if (mr > 0) RefineRadial(mr, 0);
	g = GetGroup(k);
	for (i = 0; i < g->n_cfgs; i++) {
	  cfg = (CONFIG *) ArrayGet(&(g->cfg_list), i);
	  cfg->energy = AverageEnergyConfig(cfg);
	}
	ReinitRadial(1);
	ClearOrbitalTable(0);
      }
    } else {
      OptimizeRadial(ng, kg, NULL);
      if (mr) RefineRadial(mr, 0);
      for (k = 0; k < ng; k++) {
	g = GetGroup(kg[k]);
	for (i = 0; i < g->n_cfgs; i++) {
	  cfg = (CONFIG *) ArrayGet(&(g->cfg_list), i);
	  if (cfg->energy == 0) {
	    cfg->energy = AverageEnergyConfig(cfg);
	  }
	}
      }
      ReinitRadial(1);
      ClearOrbitalTable(0);
    }
  } else {
    ng = GetNumGroups();
    for (k = 0; k < ng; k++) {
      g = GetGroup(k);
      for (i = 0; i < g->n_cfgs; i++) {
	cfg = (CONFIG *) ArrayGet(&(g->cfg_list), i);
	if (cfg->energy != 0) {
	  cfg->delta = cfg->energy - AverageEnergyConfig(cfg);
	}
      }
    }
  }
  return 0;
}

/* calculate the total configuration average energy of a group. */
double TotalEnergyGroup(int kg) {
  CONFIG_GROUP *g;
  ARRAY *c;
  CONFIG *cfg;
  int t;
  double total_energy;

  g = GetGroup(kg);
  c = &(g->cfg_list);
  
  total_energy = 0.0;
  for (t = 0; t < g->n_cfgs; t++) {
    cfg = (CONFIG *) ArrayGet(c, t);
    total_energy += AverageEnergyConfig(cfg);
  }
  return total_energy;
}

double ZerothEnergyConfig(CONFIG *cfg) {
  int i, n, nq, kappa, k;
  double r;

  r = 0.0;
  for (i = 0; i < cfg->n_shells; i++) {
    n = (cfg->shells[i]).n;
    nq = (cfg->shells[i]).nq;
    kappa = (cfg->shells[i]).kappa;
    k = OrbitalIndex(n, kappa, 0.0);
    r += nq * (GetOrbital(k)->energy);
  }
  return r;
}

static double FKB(int ka, int kb, int k) {
  int ja, jb, ia, ib;
  double a, b;
  
  GetJLFromKappa(GetOrbital(ka)->kappa, &ja, &ia);
  GetJLFromKappa(GetOrbital(kb)->kappa, &jb, &ib);

  if (!Triangle(ia, k, ia) || !Triangle(ib, k, ib)) return 0.0;
  a = W3j(ja, k, ja, 1, 0, -1)*W3j(jb, k, jb, 1, 0, -1);
  if (fabs(a) < EPS30) return 0.0; 
  Slater(&b, ka, kb, ka, kb, k/2, 0);
  
  b *= a*(ja+1.0)*(jb+1.0);
  if (IsEven((ja+jb)/2)) b = -b;

  return b;
}

static double GKB(int ka, int kb, int k) {
  int ja, jb, ia, ib;
  double a, b;
  
  GetJLFromKappa(GetOrbital(ka)->kappa, &ja, &ia);
  GetJLFromKappa(GetOrbital(kb)->kappa, &jb, &ib);
  
  if (IsOdd((ia+k+ib)/2) || !Triangle(ia, k, ib)) return 0.0;
  a = W3j(ja, k, jb, 1, 0, -1);
  if (fabs(a) < EPS30) return 0.0;
  Slater(&b, ka, kb, kb, ka, k/2, 0);
  
  b *= a*a*(ja+1.0)*(jb+1.0);
  if (IsEven((ja+jb)/2)) b = -b;
  return b;
}  
  
static double ConfigEnergyVarianceParts0(SHELL *bra, int ia, int ib, 
					 int m2, int p) {
  int ja, jb, k, kp, k0, k1, kp0, kp1, ka, kb;
  double a, b, c, d, e;

  ja = GetJFromKappa(bra[ia].kappa);
  ka = OrbitalIndex(bra[ia].n, bra[ia].kappa, 0);
  if (p > 0) {
    jb = GetJFromKappa(bra[ib].kappa);
    kb = OrbitalIndex(bra[ib].n, bra[ib].kappa, 0);
  }
  e = 0.0;
  switch (p) {
  case 0:
    k0 = 4;
    k1 = 2*ja;
    a = 1.0/(ja*(ja+1.0));
    for (k = k0; k <= k1; k += 4) {
      for (kp = k0; kp <= k1; kp += 4) {
	b = -a + W6j(ja, ja, k, ja, ja, kp);
	if (k == kp) b += 1.0/(k+1.0);
	b *= a*FKB(ka, ka, k)*FKB(ka, ka, kp);
	e += b;
      }
    }
    break;
  case 1:
    k0 = 4;
    k1 = 2*ja;
    kp0 = 4;
    kp1 = 2*jb;
    kp1 = Min(kp1, k1);
    a = 1.0/(ja*(ja+1.0));
    for (k = k0; k <= k1; k += 4) {
      for (kp = kp0; kp <= kp1; kp += 4) {
	b = a - W6j(ja, ja, kp, ja, ja, k);
	if (k == kp) b -= 1.0/(k+1.0);
	b *= W6j(ja, ja, kp, jb, jb, m2);
	b /= 0.5*ja;
	b *= FKB(ka, ka, k)*FKB(ka, kb, kp);
	e += b;
      }
    }
    if (IsOdd((ja+jb)/2+1)) e = -e;
    break;
  case 2:
    k0 = 4;
    k1 = 2*ja;
    kp0 = abs(ja-jb);
    kp1 = ja + jb;
    a = 1.0/(ja*(ja+1.0));
    for (k = k0; k <= k1; k += 4) {
      for (kp = kp0; kp <= kp1; kp += 2) {
	b = W6j(k, kp, m2, jb, ja, ja);
	b = -b*b;
	b += W6j(jb, jb, k, ja, ja, m2)*W6j(jb, jb, k, ja, ja, kp);
	if (m2 == kp) {
	  b -= (1.0/(m2+1.0)-1.0/(jb+1.0))*a;
	} else {
	  b += a/(jb+1.0);
	}
	b /= 0.5*ja;
	b *= FKB(ka, ka, k)*GKB(ka, kb, kp);
	e += b;
      }
    }
    if (IsOdd((ja+jb)/2+1)) e = -e;
    break;
  case 3:
    k0 = 4;
    k1 = 2*ja;
    k = 2*jb;
    k1 = Min(k, k1);
    for (k = k0; k <= k1; k += 4) {
      for (kp = k0; kp <= k1; kp += 4) {
	b = 0.0;
	if (k = kp) b += 1.0/((k+1.0)*(jb+1.0));
	b -= W9j(ja, ja, k, ja, m2, jb, kp, jb, jb);
	b -= W6j(ja, jb, m2, jb, ja, k)*W6j(ja, jb, m2, jb, ja, kp)/ja;
	b /= ja;
	b *= FKB(ka, kb, k)*FKB(ka, kb, kp);
	e += b;
      }
    }
    break;
  case 4:
    k0 = abs(ja-jb);
    k1 = ja+jb;
    for (k = k0; k <= k1; k += 2) {
      for (kp = k0; kp <= k1; kp += 2) {
	b = 0.0;
	if (k == kp) b += 1.0/((k+1.0)*(jb+1.0));
	b -= W9j(ja, jb, k, jb, m2, ja, kp, ja, jb);
	c = -1.0/(jb+1.0);	
	d = c;
	if (k == m2) {
	  c += 1.0/(m2+1.0);
	}
	if (kp == m2) {	  
	  d += 1.0/(m2+1.0);
	}
	b -= c*d/ja;
	b /= ja;
	b *= GKB(ka, kb, k)*GKB(ka, kb, kp);
	e += b;
      }
    }
    break;
  case 5:
    k0 = 4;
    k1 = 2*ja;
    k = 2*jb;
    k1 = Min(k1, k);
    kp0 = abs(ja-jb);
    kp1 = ja+jb;
    for (k = k0; k <= k1; k += 4) {
      for (kp = kp0; kp <= kp1; kp += 2) {
	b = -W6j(jb, jb, k, ja, ja, kp)/(jb+1.0);
	c = W6j(k, kp, m2, ja, jb, jb)*W6j(k, kp, m2, jb, ja, ja);
	if (IsOdd((ja+jb)/2)) b += c;
	else b -= c;
	c = -1.0/(jb+1.0);
	if (kp == m2) c += 1.0/(m2+1.0);
	c *= W6j(ja, jb, m2, jb, ja, k);
	b += c/ja;
	b /= 0.5*ja;
	b *= FKB(ka, kb, k)*GKB(ka, kb, kp);
	e += b;
      }
    }
    break;
  case 6:
    k0 = 4;
    k1 = 2*ja;
    k = 2*jb;
    k1 = Min(k, k1);
    a = 1.0/((ja+1.0)*(jb+1.0));
    for (k = k0; k <= k1; k += 4) {
      b = a/(k+1.0);
      c = FKB(ka, kb, k);
      b *= c*c;
      e += b;
    }
    break;
  case 7:
    k0 = abs(ja-jb);
    k1 = ja+jb;
    a = 1.0/((ja+1.0)*(jb+1.0));
    for (k = k0; k <= k1; k += 2) {
      for (kp = k0; kp <= k1; kp += 2) {
	b = 0;
	if (k == kp) {
	  b += 1.0/(k+1.0);
	}
	b -= a;
	c = GKB(ka, kb, k);
	d = GKB(ka, kb, kp);
	e += a*b*c*d;
      }
    }
    break;
  case 8:    
    k0 = 4;
    k1 = 2*ja;
    k = 2*jb;
    k1 = Min(k, k1);
    kp0 = abs(ja-jb);
    kp1 = ja+jb;
    kp0 = k0;
    kp1 = k1;    
    a = 1.0/((ja+1.0)*(jb+1.0));
    for (k = k0; k <= k1; k += 4) {
      for (kp = kp0; kp <= kp1; kp += 2) {
	b = W6j(jb, ja, kp, ja, jb, k);
	if (fabs(b) < EPS30) continue;
	b *= 2.0*a;
	if (IsOdd(kp/2)) b = -b;
	c = FKB(ka, kb, k);
	d = GKB(ka, kb, kp);
	e += b*c*d;
      }
    }
    break;
  }

  return e;
}

static double ConfigEnergyVarianceParts1(SHELL *bra, int i, 
					 int ia, int ib, int m2, int p) {  
  int js, ja, jb, k, kp, k0, k1, kp0, kp1, ka, kb, ks;
  double a, b, e;
  
  js = GetJFromKappa(bra[i].kappa);
  ks = OrbitalIndex(bra[i].n, bra[i].kappa, 0);
  ja = GetJFromKappa(bra[ia].kappa);
  ka = OrbitalIndex(bra[ia].n, bra[ia].kappa, 0);
  jb = GetJFromKappa(bra[ib].kappa);
  kb = OrbitalIndex(bra[ib].n, bra[ib].kappa, 0);
  e = 0.0;

  switch (p) {
  case 0:
    k0 = 4;
    k1 = 2*ja;
    k = 2*jb;
    k1 = Min(k, k1);
    k = 2*js;
    k1 = Min(k, k1);
    for (k = k0; k <= k1; k += 4) {
      b = W6j(ja, ja, k, jb, jb, m2);
      if (fabs(b) < EPS30) continue;
      b *= 2.0/((k+1.0)*(js+1.0));
      b *= FKB(ks, ka, k)*FKB(ks, kb, k);
      if (IsOdd((ja+jb)/2)) b = -b;
      e += b;
    }
    break;
  case 1:
    k0 = abs(js-ja);
    k1 = js+ja;
    kp0 = abs(js-jb);
    kp1 = js+jb;
    a = 1.0/((js+1.0)*(ja+1.0)*(jb+1.0));
    for (k = k0; k <= k1; k += 2) {
      for (kp = kp0; kp <= kp1; kp += 2) {
	b = W6j(k, kp, m2, jb, ja, js);
	b = -b*b + a;
	b /= (js+1.0);
	b *= 2.0*GKB(ks, ka, k)*GKB(ks, kb, kp);
	if (IsOdd((ja+jb)/2+1)) b = -b;
	e += b;
      }
    }
    break;
  case 2:
    k0 = 4;
    k1 = 2*js;    
    k = 2*ja;
    k1 = Min(k, k1);
    kp0 = abs(js-jb);
    kp1 = js+jb;
    for (k = k0; k <= k1; k += 4) {
      for (kp = kp0; kp <= kp1; kp += 2) {
	b = W6j(jb, jb, k, ja, ja, m2);
	b *= W6j(jb, jb, k, js, js, kp);
	if (fabs(b) < EPS30) continue;
	b /= (js+1.0);
	if (IsOdd((ja+jb+kp)/2)) b = -b;
	b *= 2.0*FKB(ks, ka, k)*GKB(ks, kb, kp);
	e += b;
      }
    }
    break;
  }

  return e;
}

double ConfigEnergyVariance(int ns, SHELL *bra, int ia, int ib, int m2) {
  int i, js, p;
  double e, a, b, c;
  
  e = 0.0;
  for (i = 0; i < ns; i++) {
    js = GetJFromKappa(bra[i].kappa);
    a = bra[i].nq;
    b = js+1.0 - bra[i].nq;
    if (i == ia) {
      a -= 1.0;      
    }
    if (i == ib) {
      b -= 1.0;
    }
    if (a == 0.0 || b == 0.0) continue;
    a = a*b;
    b = 0.0;
    if (i == ia) {
      for (p = 0; p < 6; p++) {
	c = ConfigEnergyVarianceParts0(bra, ia, ib, m2, p);
	b += c;
      }
      b /= js-1.0;
    } else if (i == ib) {
      for (p = 0; p < 6; p++) {
	c = ConfigEnergyVarianceParts0(bra, ib, ia, m2, p);
	b += c;
      }
      b /= js-1.0;
    } else {
      for (p = 6; p < 9; p++) {
	c = ConfigEnergyVarianceParts0(bra, i, ia, m2, p);
	b += c;
	c = ConfigEnergyVarianceParts0(bra, i, ib, m2, p);
	b += c;
      }
      c = ConfigEnergyVarianceParts1(bra, i, ia, ib, m2, 0);
      b += c;
      c = ConfigEnergyVarianceParts1(bra, i, ia, ib, m2, 1);
      b += c;
      c = ConfigEnergyVarianceParts1(bra, i, ia, ib, m2, 2);
      b += c;
      c = ConfigEnergyVarianceParts1(bra, i, ib, ia, m2, 2);
      b += c;
      b /= js;
    }

    e += a*b;
  }
    
  if (e < 0.0) e = 0.0;
  return e;
}

double ConfigEnergyShiftCI(int nrs0, int nrs1) {
  int n0, k0, q0, n1, k1, q1;
  int j0, j1, s0, s1, kmax, k;
  double a, g, pk, w, sd, q;

  UnpackNRShell(&nrs0, &n0, &k0, &q0);
  UnpackNRShell(&nrs1, &n1, &k1, &q1);
  q = (q0-1.0)/(2.0*k0+1.0) - q1/(2.0*k1+1.0);
  if (q == 0.0) return q;
  g = 0.0;
  for (j0 = k0-1; j0 <= k0+1; j0 += 2) {
    if (j0 < 0) continue;
    s0 = OrbitalIndex(n0, GetKappaFromJL(j0, k0), 0);    
    for (j1 = k1-1; j1 <= k1+1; j1 += 2) {
      if (j1 < 0) continue;
      s1 = OrbitalIndex(n1, GetKappaFromJL(j1, k1), 0);
      w = W6j(j0, k0, 1, k1, j1, 2);
      w = w*w*(j0+1.0)*(j1+1.0)*0.5;
      pk = ((j0+1.0)*(j1+1.0))/((k0+1.0)*(k1+1.0)*4.0) - w;
      Slater(&sd, s0, s1, s0, s1, 0, 0);
      g += sd*pk;
      kmax = 2*j0;
      k = 2*j1;
      kmax = Min(k, kmax);
      for (k = 4; k <= kmax; k += 4) {	
	pk = W3j(k0, k, k0, 0, 0, 0);
	if (fabs(pk) > 0) {
	  pk *= W3j(k1, k, k1, 0, 0, 0);
	  if (fabs(pk) > 0) {
	    pk *= 0.25;
	    a = W6j(j0, 2, j1, k1, 1, k0);
	    if (fabs(a) > 0) {
	      a = a*a;
	      a *= W6j(j0, k, j0, k0, 1, k0);
	      if (fabs(a) > 0) {
		a *= W6j(j1, k, j1, k1, 1, k1);
		if (fabs(a) > 0) {
		  a *= W6j(j0, k, j0, j1, 2, j1);
		  a *= 2.0*(j0+1.0)*(j1+1.0)*(k0+1.0)*(k1+1.0);
		}
	      }
	    }
	    a += W6j(k0, k, k0, k1, 2, k1);
	    pk *= a;
	  }
	  if (fabs(pk) > 0) {
	    Slater(&sd, s0, s1, s0, s1, k/2, 0);
	    g += pk*sd*(j0+1.0)*(j1+1.0);
	  }
	}
      }
      pk = W3j(k0, 2, k1, 0, 0, 0);
      if (fabs(pk) > 0) {
	pk = pk*pk;
	a = W6j(j0, 2, j1, k1, 1, k0);
	a *= a;
	pk *= a;
	pk *= (k0+1.0)*(k1+1.0)/3.0;
	pk *= (1.0 - 0.5*(j0+1.0)*(j1+1.0)*a);
	if (fabs(pk) > 0) {
	  Slater(&sd, s0, s1, s1, s0, 1, 0);
	  g += pk*sd*(j0+1.0)*(j1+1.0);
	}
      }
    }
  }

  return q*g;
}
	  
double ConfigEnergyShift(int ns, SHELL *bra, int ia, int ib, int m2) {
  double qa, qb, a, b, c, sd, e;
  int ja, jb, k, kmin, kmax;
  int k0, k1;

  qa = bra[ia].nq;
  qb = bra[ib].nq;
  ja = GetJFromKappa(bra[ia].kappa);
  jb = GetJFromKappa(bra[ib].kappa);
  if (qa == 1 && qb == 0) e = 0.0;
  else {
    e = (qa-1.0)/ja - qb/jb;
    if (e != 0.0) {
      kmin = 4;
      kmax = 2*ja;
      k = 2*jb;
      kmax = Min(kmax, k);
      a = 0.0;
      k0 = OrbitalIndex(bra[ia].n, bra[ia].kappa, 0);
      k1 = OrbitalIndex(bra[ib].n, bra[ib].kappa, 0);
      for (k = kmin; k <= kmax; k += 4) {
	b = W6j(k, ja, ja, m2, jb, jb);	
	if (fabs(b) > EPS30) {
	  a -= b*FKB(k0, k1, k);
	}
      }
      kmin = abs(ja-jb);
      kmax = ja + jb;
      c = 1.0/((ja+1.0)*(jb+1.0));
      for (k = kmin; k <= kmax; k += 2) {
	if (k == m2) {	  
	  b = 1.0/(m2+1.0)-c;
	} else {
	  b = -c;
	}
	a += b*GKB(k0, k1, k);
      }
      if (IsEven((ja+jb)/2)) a = -a;
      e *= a;
    }
  }

  return e;
}

/* calculate the average energy of a configuration */
double AverageEnergyConfig(CONFIG *cfg) {
  int i, j, n, kappa, nq, np, kappap, nqp;
  int k, kp, kk, kl, klp, kkmin, kkmax, j2, j2p;
  double x, y, t, q, a, b, r;
 
  x = 0.0;
  for (i = 0; i < cfg->n_shells; i++) {
    n = (cfg->shells[i]).n;
    kappa = (cfg->shells[i]).kappa;
    kl = GetLFromKappa(kappa);
    j2 = GetJFromKappa(kappa);
    nq = (cfg->shells[i]).nq;
    k = OrbitalIndex(n, kappa, 0.0);
    
    if (nq > 1) {
      t = 0.0;
      for (kk = 2; kk <= j2; kk += 2) {
	Slater(&y, k, k, k, k, kk, 0);
	q = W3j(j2, 2*kk, j2, -1, 0, 1);
	t += y * q * q ;
      }
      Slater(&y, k, k, k, k, 0, 0);
      b = ((nq-1.0)/2.0) * (y - (1.0 + 1.0/j2)*t);

#if FAC_DEBUG
      fprintf(debug_log, "\nAverage Radial: %lf\n", y);
#endif

    } else {
      b = 0.0;
    }

    t = 0.0;
    for (j = 0; j < i; j++) {
      np = (cfg->shells[j]).n;
      kappap = (cfg->shells[j]).kappa;
      klp = GetLFromKappa(kappap);
      j2p = GetJFromKappa(kappap);
      nqp = (cfg->shells[j]).nq;
      kp = OrbitalIndex(np, kappap, 0.0);

      kkmin = abs(j2 - j2p);
      kkmax = (j2 + j2p);
      if (IsOdd((kkmin + kl + klp)/2)) kkmin += 2;
      a = 0.0;
      for (kk = kkmin; kk <= kkmax; kk += 4) {
	Slater(&y, k, kp, kp, k, kk/2, 0);
	q = W3j(j2, kk, j2p, -1, 0, 1);
	a += y * q * q;
#if FAC_DEBUG
	fprintf(debug_log, "exchange rank: %d, q*q: %lf, Radial: %lf\n", 
		kk/2, q*q, y);
#endif

      }
      Slater(&y, k, kp, k, kp, 0, 0);

#if FAC_DEBUG
      fprintf(debug_log, "direct: %lf\n", y);
#endif

      t += nqp * (y - a);
    }

    ResidualPotential(&y, k, k);

    r = nq * (b + t + GetOrbital(k)->energy + y);
    x += r;
  }
  return x;
}

/* calculate the average energy of an average configuration */
double AverageEnergyAvgConfig(AVERAGE_CONFIG *cfg) {
  int i, j, n, kappa, np, kappap;
  int k, kp, kk, kl, klp, kkmin, kkmax, j2, j2p;
  double x, y, t, q, a, b, r, nq, nqp;
 
  x = 0.0;
  for (i = 0; i < cfg->n_shells; i++) {
    n = cfg->n[i];
    kappa = cfg->kappa[i];
    kl = GetLFromKappa(kappa);
    j2 = GetJFromKappa(kappa);
    nq = cfg->nq[i];
    k = OrbitalIndex(n, kappa, 0.0);
    
    t = 0.0;
    for (kk = 2; kk <= j2; kk += 2) {
      Slater(&y, k, k, k, k, kk, 0);
      q = W3j(j2, 2*kk, j2, -1, 0, 1);
      t += y * q * q ;
    }
    Slater(&y, k, k, k, k, 0, 0);
    b = ((nq-1.0)/2.0) * (y - (1.0 + 1.0/j2)*t);
    
#if FAC_DEBUG
      fprintf(debug_log, "\nAverage Radial: %lf\n", y);
#endif
      
    t = 0.0;
    for (j = 0; j < i; j++) {
      np = cfg->n[j];
      kappap = cfg->kappa[j];
      klp = GetLFromKappa(kappap);
      j2p = GetJFromKappa(kappap);
      nqp = cfg->nq[j];
      kp = OrbitalIndex(np, kappap, 0.0);

      kkmin = abs(j2 - j2p);
      kkmax = (j2 + j2p);
      if (IsOdd((kkmin + kl + klp)/2)) kkmin += 2;
      a = 0.0;
      for (kk = kkmin; kk <= kkmax; kk += 4) {
	Slater(&y, k, kp, kp, k, kk/2, 0);
	q = W3j(j2, kk, j2p, -1, 0, 1);
	a += y * q * q;
#if FAC_DEBUG
	fprintf(debug_log, "exchange rank: %d, q*q: %lf, Radial: %lf\n", 
		kk/2, q*q, y);
#endif

      }
      Slater(&y, k, kp, k, kp, 0, 0);

#if FAC_DEBUG
      fprintf(debug_log, "direct: %lf\n", y);
#endif

      t += nqp * (y - a);
    }

    ResidualPotential(&y, k, k);
    a = GetOrbital(k)->energy;
    r = nq * (b + t + a + y);
    x += r;
  }
  return x;
}

/* calculate the expectation value of the residual potential:
   -Z/r - v0(r), where v0(r) is central potential used to solve 
   dirac equations. the orbital index must be valid, i.e., their 
   radial equations must have been solved. */
int ResidualPotential(double *s, int k0, int k1) {
  int i;
  ORBITAL *orb1, *orb2;
  int index[2];
  double *p, z, *p1, *p2, *q1, *q2;

  orb1 = GetOrbitalSolved(k0);
  orb2 = GetOrbitalSolved(k1);
  if (!orb1 || !orb2) return -1;
  if (orb1->wfun == NULL || orb2->wfun == NULL) {
    *s = 0.0;
    return 0;
  }

  if (k0 > k1) {
    index[0] = k1;
    index[1] = k0;
  } else {
    index[0] = k0;
    index[1] = k1;
  }
  
  p = (double *) MultiSet(residual_array, index, NULL, InitDoubleData, NULL);
  if (p && *p) {
    *s = *p;
    return 0;
  } 

  *s = 0.0;
 
  if (orb1->n < 0 || orb2->n < 0) {
    p1 = Large(orb1);
    p2 = Large(orb2);
    q1 = Small(orb1);
    q2 = Small(orb2);
    for (i = potential->ib; i <= potential->ib1; i++) {
      z = potential->U[i];
      z += potential->Vc[i];
      _yk[i] = -(potential->Z[i]/potential->rad[i]) - z;
      _yk[i] *= potential->dr_drho[i];
      _yk[i] *= p1[i]*p2[i] + q1[i]*q2[i];
    }
    *s = Simpson(_yk, potential->ib, potential->ib1);
  } else {
    for (i = 0; i < potential->maxrp; i++) {
      z = potential->U[i];
      z += potential->Vc[i];
      _yk[i] = -(potential->Z[i]/potential->rad[i]) - z;
    }
    Integrate(_yk, orb1, orb2, 1, s, -1);
  }
  *p = *s;
  return 0;
}

double MeanPotential(int k0, int k1) {
  int i;
  ORBITAL *orb1, *orb2;
  double z, *p1, *p2, *q1, *q2;

  orb1 = GetOrbitalSolved(k0);
  orb2 = GetOrbitalSolved(k1);
  if (!orb1 || !orb2) return -1;
  if (orb1->wfun == NULL || orb2->wfun == NULL) {
    return 0.0;
  }

  if (orb1->n < 0 || orb2->n < 0) {
    p1 = Large(orb1);
    p2 = Large(orb2);
    q1 = Small(orb1);
    q2 = Small(orb2);
    for (i = potential->ib; i <= potential->ib1; i++) {
      z = potential->U[i];
      z += potential->Vc[i];
      _yk[i] = z;
      _yk[i] *= potential->dr_drho[i];
      _yk[i] *= p1[i]*p2[i] + q1[i]*q2[i];
    }
    z = Simpson(_yk, potential->ib, potential->ib1);
  } else {
    for (i = 0; i < potential->maxrp; i++) {
      z = potential->U[i];
      z += potential->Vc[i];
      _yk[i] = z;
    }
    Integrate(_yk, orb1, orb2, 1, &z, -1);
  }

  return z;
}

double RadialMoments(int m, int k1, int k2) {
  int index[3];
  int npts, i0, i;
  ORBITAL *orb1, *orb2;
  double *q, r, z, *p1, *p2, *q1, *q2;
  int n1, n2;
  int kl1, kl2;
  int nh, klh;
  
  orb1 = GetOrbitalSolved(k1);
  orb2 = GetOrbitalSolved(k2);
  n1 = orb1->n;
  n2 = orb2->n;
  kl1 = orb1->kappa;
  kl2 = orb2->kappa;
  kl1 = GetLFromKappa(kl1);
  kl2 = GetLFromKappa(kl2);
  kl1 /= 2;
  kl2 /= 2;	

  GetHydrogenicNL(&nh, &klh, NULL, NULL);

  if (n1 > 0 && n2 > 0 && potential->ib <= 0) {
    if ((n1 > nh && n2 > nh) || 
	(kl1 > klh && kl2 > klh) ||
	orb1->wfun == NULL || 
	orb2->wfun == NULL) {
      if (n1 == n2 && kl1 == kl2) {
	z = GetResidualZ();
	r = HydrogenicExpectation(z, m, n1, kl1);
	if (r) {
	  return r;
	}
      } else if (m == 1) {
	z = GetResidualZ();
	if (n1 < n2) {
	  r = HydrogenicDipole(z, n1, kl1, n2, kl2);
	  return r;
	} else if (n1 < n2) {
	  r = HydrogenicDipole(z, n2, kl2, n1, kl1);
	  return r;
	}
      }
    }
  }

  if (potential->ib <= 0 && n1 == n2 && m > 0 && n1 > GetNMax()) {
    return 0.0;
  }
  if (orb1->wfun == NULL || orb2->wfun == NULL) {
    return 0.0;
  }

  if (m >= 0) {
    index[0] = 2*m;
  } else {
    index[0] = -2*m-1;
  }
    
  if (k1 < k2) {
    index[1] = k1;
    index[2] = k2;
  } else {
    index[1] = k2;
    index[2] = k1;
  }
  
  q = (double *) MultiSet(moments_array, index, NULL, InitDoubleData, NULL);
 
  if (*q) {
    return *q;
  } 

  if (n1 < 0 || n2 < 0) {
    i0 = potential->ib;
    npts = potential->ib1;
    p1 = Large(orb1);
    q1 = Small(orb1);
    p2 = Large(orb2);
    q2 = Small(orb2);
    for (i = i0; i <= npts; i++) {
      r = p1[i]*p2[i] + q1[i]*q2[i];
      r *= potential->dr_drho[i];
      _yk[i] = pow(potential->rad[i], m)*r;
    }
    r = Simpson(_yk, i0, npts);
    *q = r;
  } else {    
    npts = potential->maxrp-1;
    if (n1 != 0) npts = Min(npts, orb1->ilast);
    if (n2 != 0) npts = Min(npts, orb2->ilast);

    for (i = 0; i <= npts; i++) {
      _yk[i] = pow(potential->rad[i], m);
    }
    r = 0.0;
    Integrate(_yk, orb1, orb2, 1, &r, m);
    *q = r;
  }
  return r;
}

double MultipoleRadialNR(int m, int k1, int k2, int gauge) {
  int i, p, t;
  ORBITAL *orb1, *orb2;
  double r;
  int kappa1, kappa2;

#ifdef PERFORM_STATISTICS
  clock_t start, stop; 
  start = clock();
#endif

  orb1 = GetOrbital(k1);
  orb2 = GetOrbital(k2);
  kappa1 = orb1->kappa;
  kappa2 = orb2->kappa;

  r = 0.0;
  if (m == 1) {
    /* use the relativistic version is just simpler */
    printf("should call MultipoleRadialFR instead\n");
  } else if (m > 1) {
    t = kappa1 + kappa2;
    p = m - t;
    if (p && t) {
      r = RadialMoments(m-1, k1, k2);
      r *= p*t;
      r /= sqrt(m*(m+1.0));
      r *= -0.5 * FINE_STRUCTURE_CONST;
      for (i = 2*m - 1; i > 0; i -= 2) r /= (double) i;
      r *= ReducedCL(GetJFromKappa(kappa1), 2*m, GetJFromKappa(kappa2));
    }
  } else if (m < 0) {
    m = -m;
    if (gauge == G_BABUSHKIN) {
      r = RadialMoments(m, k1, k2);
      r *= sqrt((m+1.0)/m);
      for (i = 2*m - 1; i > 1; i -= 2) r /= (double) i;
    } else {
      /* the velocity form is not implemented yet. 
	 the following is still the length form */
      r = RadialMoments(m, k1, k2);
      r *= sqrt((m+1.0)/m);
      for (i = 2*m - 1; i > 1; i -= 2) r /= (double) i;
    }
    r *= ReducedCL(GetJFromKappa(kappa1), 2*m, GetJFromKappa(kappa2));
  }

#ifdef PERFORM_STATISTICS 
  stop = clock();
  rad_timing.radial_1e += stop - start;
#endif
  
  return r;
}

/* fully relativistic multipole operator, 
   see Grant, J. Phys. B. 1974. Vol. 7, 1458. */ 
double MultipoleRadialFR(double aw, int m, int k1, int k2, int gauge) {
  double q, ip, ipm, im, imm;
  int kappa1, kappa2;
  int am, t;
  int index[4], s;
  ORBITAL *orb1, *orb2, *orb;
  double x, a, r, rp, **p1, aw2, ef;
  int jy, n, i, j, npts;
  double rcl;

#ifdef PERFORM_STATISTICS 
  clock_t start, stop;
  start = clock();
#endif

  if (m == 0) return 0.0;
  
  if (m >= 0) {
    index[0] = 2*m;
    am = m;
  } else {
    index[0] = -2*m-1;
    am = -m;
  }

  orb1 = GetOrbitalSolved(k1);
  orb2 = GetOrbitalSolved(k2);
  if (orb1->wfun == NULL || orb2->wfun == NULL) {
    if (m == -1) {
      return MultipoleRadialNR(m, k1, k2, gauge);
    } else {
      return 0.0;
    }
  }
  
  s = 0;      
  index[1] = k1;
  index[2] = k2;
  kappa1 = orb1->kappa;
  kappa2 = orb2->kappa;
  rcl = ReducedCL(GetJFromKappa(kappa1), abs(2*m), 
		  GetJFromKappa(kappa2));
  /*
  if (orb1->energy <= orb2->energy) {
    s = 0;      
    index[1] = k1;
    index[2] = k2;
    kappa1 = orb1->kappa;
    kappa2 = orb2->kappa;
    rcl = ReducedCL(GetJFromKappa(kappa1), abs(2*m), 
		    GetJFromKappa(kappa2));
  } else {
    s = 1;     
    index[1] = k2;
    index[2] = k1;
    orb = orb1;
    orb1 = orb2;
    orb2 = orb;
    kappa1 = orb1->kappa;
    kappa2 = orb2->kappa;
    rcl = ReducedCL(GetJFromKappa(kappa2), abs(2*m), 
		    GetJFromKappa(kappa1));
  }
  */

  ef = Max(orb1->energy, orb2->energy);
  if (ef > 0.0) {
    ef *= FINE_STRUCTURE_CONST;
    if (n_awgrid > 1) {
      for (i = 0; i < n_awgrid; i++) {
	aw2grid[i] = awgrid[i] + ef;
	aw2grid[i] *= aw2grid[i];
      }
    }
  }

  if (n_awgrid > 1) {
    if (ef > 0) aw += ef;
    aw2 = aw*aw;
  }

  p1 = (double **) MultiSet(multipole_array, index, NULL, 
			    InitPointerData, FreeMultipole);
  if (*p1) {
    r = InterpolateMultipole(aw2, n_awgrid, aw2grid, *p1);
    if (s == 1 && IsOdd(am)) r = -r;
    r *= rcl;
    return r;
  }  
  *p1 = (double *) malloc(sizeof(double)*n_awgrid);
  
  npts = potential->maxrp-1;
  if (orb1->n > 0) npts = Min(npts, orb1->ilast);
  if (orb2->n > 0) npts = Min(npts, orb2->ilast);
  r = 0.0;
  jy = 1;

  for (i = 0; i < n_awgrid; i++) {
    r = 0.0;
    a = awgrid[i];
    if (ef > 0.0) a += ef;
    if (m > 0) {
      t = kappa1 + kappa2;
      if (t) {
	for (j = 0; j <= npts; j++) {
	  x = a*potential->rad[j];
	  n = m;
	  _yk[j] = BESLJN(jy, n, x);
	}
	Integrate(_yk, orb1, orb2, 4, &r, 0);
	r *= t;
	r *= (2*m + 1.0)/sqrt(m*(m+1.0));
	r /= pow(a, m);
	(*p1)[i] = r;
      }
    } else {
      if (gauge == G_COULOMB) {
	t = kappa1 - kappa2;
	q = sqrt(am/(am+1.0));
	for (j = 0; j <= npts; j++) {
	  x = a*potential->rad[j];
	  n = am+1;
	  _yk[j] = BESLJN(jy, n, x);
	  n = am-1;
	  _zk[j] = BESLJN(jy, n, x);
	}
	if (t) {
	  Integrate(_yk, orb1, orb2, 4, &ip, 0);
	  Integrate(_zk, orb1, orb2, 4, &ipm, 0);
	  r = t*ip*q - t*ipm/q;
	}
	if (k1 != k2) {
	  Integrate(_yk, orb1, orb2, 5, &im, 0);
	  Integrate(_zk, orb1, orb2, 5, &imm, 0);
	  r += (am + 1.0)*im*q + am*imm/q;
	}
	r /= pow(a,am);
	(*p1)[i] = r;
      } else if (gauge == G_BABUSHKIN) {
	t = kappa1 - kappa2;
	for (j = 0; j < npts; j++) {
	  x = a*potential->rad[j];
	  n = am+1;
	  _yk[j] = BESLJN(jy, n, x);
	  n = am;
	  _zk[j] = BESLJN(jy, n, x);
	}
	if (t) {
	  Integrate(_yk, orb1, orb2, 4, &ip, 0);
	  r = t*ip;
	}
	if (k1 != k2) {
	  Integrate(_yk, orb1, orb2, 5, &im, 0);
	} else {
	  im = 0.0;
	}
	Integrate(_zk, orb1, orb2, 1, &imm, 0);
	rp = (am + 1.0) * (imm + im);
	q = (2*am + 1.0)/sqrt(am*(am+1.0));
	q /= pow(a, am);
	r *= q;
	rp *= q;
	(*p1)[i] = r+rp;
      }
    }
  }

  r = InterpolateMultipole(aw2, n_awgrid, aw2grid, *p1);
  if (s == 1 && IsOdd(am)) r = -r;
  r *= rcl;

#ifdef PERFORM_STATISTICS 
  stop = clock();
  rad_timing.radial_1e += stop - start;
#endif
  return r;
}

double *GeneralizedMoments(int k1, int k2, int m) {
  ORBITAL *orb1, *orb2;
  int n1, i, jy, nk;
  double x, r, r0;
  double *p1, *p2, *q1, *q2;
  int index[3], t;
  double **p, k, *kg;
  double amin, amax, kmin, kmax;
  
  index[0] = m;
  if (k1 > k2) {
    index[1] = k2;
    index[2] = k1;
    orb1 = GetOrbitalSolved(k2);
    orb2 = GetOrbitalSolved(k1);
  } else {
    index[1] = k1;
    index[2] = k2;
    orb1 = GetOrbitalSolved(k1);
    orb2 = GetOrbitalSolved(k2);
  }

  p = (double **) MultiSet(gos_array, index, NULL, 
			   InitPointerData, FreeMultipole);
  if (*p) {
    return *p;
  }

  nk = NGOSK;
  *p = (double *) malloc(sizeof(double)*nk*2);
  kg = *p + nk;

  if (orb1->wfun == NULL || orb2->wfun == NULL || 
      (orb1->n <= 0 && orb2->n <= 0)) {
    for (t = 0; t < nk*2; t++) {
      (*p)[t] = 0.0;
    }
    return *p;
  }
  
  jy = 1;
  p1 = Large(orb1);
  p2 = Large(orb2);
  q1 = Small(orb1);
  q2 = Small(orb2);

  amin = sqrt(2.0*fabs(orb1->energy));
  amax = sqrt(2.0*fabs(orb2->energy));
  if (amin < amax) {
    kmin = amin;
    kmax = amax;
  } else {
    kmin = amax;
    kmax = amin;
  }
  kmin = log(0.05*kmin);
  kmax = log(50.0*kmax);
  r = (kmax - kmin)/(nk-1.0);
  kg[0] = kmin;
  for (i = 1; i < nk; i++) {
    kg[i] = kg[i-1] + r;
  }
  
  if (orb1->n > 0 && orb2->n > 0) {
    n1 = Min(orb1->ilast, orb2->ilast);
  
    for (i = 0; i <= n1; i++) {
      _phase[i] = (p1[i]*p2[i] + q1[i]*q2[i])*potential->dr_drho[i];
    }
    
    if (m == 0) {
      if (k1 == k2) r0 = 1.0;
      else if (orb1->n != orb2->n) r0 = 0.0;
      else {
	if (orb1->kappa + orb2->kappa != -1) r0 = 0.0;
	else {
	  r0 = Simpson(_phase, 0, n1);
	}
      }
    } else {
      r0 = 0.0;
    }
    
    for (t = 0; t < nk; t++) {
      k = exp(kg[t]);
      for (i = 0; i <= n1; i++) {
	x = k * potential->rad[i];
	_dphase[i] = BESLJN(jy, m, x);
	_dphase[i] *= _phase[i];
      }
      r = Simpson(_dphase, 0, n1);
      
      (*p)[t] = (r - r0)/k;
    }
  } else {
    if (orb1->n > 0) n1 = orb1->ilast;
    else n1 = orb2->ilast;
    for (t = 0; t < nk; t++) {
      k = exp(kg[t]);      
      for (i = 0; i <= n1; i++) {
	x = k * potential->rad[i];
	_yk[i] = BESLJN(jy, m, x);
      }
      Integrate(_yk, orb1, orb2, 1, &r, 0);
      (*p)[t] = r/k;
    }
  }
  return *p;
}

void PrintGeneralizedMoments(char *fn, int m, int n0, int k0, 
			     int n1, int k1, double e1) {
  FILE *f;
  int i0, i1, i;
  double *g, *x;

  i0 = OrbitalIndex(n0, k0, 0.0);
  e1 /= HARTREE_EV;
  i1 = OrbitalIndex(n1, k1, e1);
  f = fopen(fn, "w");
  if (f == NULL) {
    printf("cannot open file %s\n", fn);
    return;
  }
  g = GeneralizedMoments(i0, i1, m);
  x = g + NGOSK;
  for (i = 0; i < NGOSK; i++) {
    fprintf(f, "%15.8E %15.8E %15.8E\n", x[i], exp(x[i]), g[i]);
  }
  fclose(f);
}
  
double InterpolateMultipole(double aw2, int n, double *x, double *y) {
  double r;
  int np, nd;

  if (n == 1) {
    r = y[0];
  } else {
    np = 3;
    nd = 1;
    UVIP3P(np, n, x, y, nd, &aw2, &r);
  }

  return r;
}

int SlaterTotal(double *sd, double *se, int *j, int *ks, int k, int mode) {
  int t, kk, tt, maxn;
  int tmin, tmax;
  double e, a, d, a1, a2, am;
  int err;
  int kl0, kl1, kl2, kl3;
  int k0, k1, k2, k3;
  int js[4];
  ORBITAL *orb0, *orb1, *orb2, *orb3;

#ifdef PERFORM_STATISTICS 
  clock_t start, stop;
  start = clock();
#endif

  k0 = ks[0];
  k1 = ks[1];
  k2 = ks[2];
  k3 = ks[3];
  kk = k/2;

  maxn = 0;
  orb0 = GetOrbitalSolved(k0);
  orb1 = GetOrbitalSolved(k1);
  orb2 = GetOrbitalSolved(k2);
  orb3 = GetOrbitalSolved(k3);

  if (orb0->n <= 0) {
    maxn = -1;
  } else if (orb0->n > maxn) {
    maxn = orb0->n;
    if (orb1->n <= 0) {
      maxn = -1;
    } else if (orb1->n > maxn) {
      maxn = orb1->n;
      if (orb2->n <= 0) {
	maxn = -1;
      } else if (orb2->n > maxn) {
	maxn = orb2->n;
	if (orb3->n <= 0) {
	  maxn = -1;
	} else if (orb3->n > maxn) {
	  maxn = orb3->n;
	}
      }
    }
  }

  if (orb0->wfun == NULL || orb1->wfun == NULL ||
      orb2->wfun == NULL || orb3->wfun == NULL) {
    if (sd) *sd = 0.0;
    if (se) *se = 0.0;
    return 0;
  }

  kl0 = GetLFromKappa(orb0->kappa);
  kl1 = GetLFromKappa(orb1->kappa);
  kl2 = GetLFromKappa(orb2->kappa);
  kl3 = GetLFromKappa(orb3->kappa);

  if (orb1->n < 0 || orb3->n < 0) {
    mode = 2;
  } else {
    if (kl1 > slater_cut.kl0 && kl3 > slater_cut.kl0) {
      if (se) {
	*se = 0.0;
	se = NULL;
      }
    }
    if (kl0 > slater_cut.kl0 && kl2 > slater_cut.kl0) {
      if (se) {
	*se = 0.0;
	se = NULL;
      }
    }  
    if (kl1 > slater_cut.kl1 && kl3 > slater_cut.kl1) {
      mode = 2;
    }
    if (kl0 > slater_cut.kl1 && kl2 > slater_cut.kl1) {
      mode = 2;
    }
  }
  if (qed.br == 0 && IsOdd((kl0+kl1+kl2+kl3)/2)) {
    if (sd) *sd = 0.0;
    if (se) *se = 0.0;
    return 0;
  }

  if (j) {
    memcpy(js, j, sizeof(int)*4);
  } else {
    js[0] = 0;
    js[1] = 0;
    js[2] = 0;
    js[3] = 0;
  }

  if (js[0] <= 0) js[0] = GetJFromKappa(orb0->kappa);
  if (js[1] <= 0) js[1] = GetJFromKappa(orb1->kappa);
  if (js[2] <= 0) js[2] = GetJFromKappa(orb2->kappa);
  if (js[3] <= 0) js[3] = GetJFromKappa(orb3->kappa);  

  am = AMU * GetAtomicMass();
  if (sd) {
    d = 0.0;
    if (IsEven((kl0+kl2)/2+kk) && IsEven((kl1+kl3)/2+kk) &&
	Triangle(js[0], js[2], k) && Triangle(js[1], js[3], k)) {
      err = Slater(&d, k0, k1, k2, k3, kk, mode);
      if (kk == 1 && qed.sms && maxn > 0) {
	a1 = Vinti(k0, k2);
	a2 = Vinti(k1, k3);
	d -= a1 * a2 / am;
      }
      a1 = ReducedCL(js[0], k, js[2]);
      a2 = ReducedCL(js[1], k, js[3]); 
      d *= a1*a2;
      if (k0 == k1 && k2 == k3) d *= 0.5;
    }
    *sd = d;
    d = 0.0;
    if (qed.br < 0 || (maxn > 0 && maxn <= qed.br)) {
      if (Triangle(js[0], js[2], k) && Triangle(js[1], js[3], k)) {
	d = Breit(k0, k1, k2, k3, kk, kl0, kl1, kl2, kl3);
	a1 = ReducedCL(js[0], k, js[2]);
	a2 = ReducedCL(js[1], k, js[3]);
	d *= a1*a2;
	if (k0 == k1 && k2 == k3) d *= 0.5;
      }
    }
    *sd += d;
  }
  
  if (!se) goto EXIT;

  if (abs(mode) == 2) {
    *se = 0.0;
    goto EXIT;
  }
  *se = 0.0;
  if (k0 == k1 && (orb0->n > 0 || orb1->n > 0)) goto EXIT;
  if (k2 == k3 && (orb2->n > 0 || orb3->n > 0)) goto EXIT;
  tmin = abs(js[0] - js[3]);
  tt = abs(js[1] - js[2]);
  tmin = Max(tt, tmin);
  tmax = js[0] + js[3];
  tt = js[1] + js[2];
  tmax = Min(tt, tmax);
  tmax = Min(tmax, GetMaxRank());
  if (IsOdd(tmin)) tmin++;
  
  for (t = tmin; t <= tmax; t += 2) {
    a = W6j(js[0], js[2], k, js[1], js[3], t);
    if (fabs(a) > EPS30) {
      e = 0.0;
      if (IsEven((kl0+kl3+t)/2) && IsEven((kl1+kl2+t)/2)) {
	err = Slater(&e, k0, k1, k3, k2, t/2, mode);
	if (t == 2 && qed.sms && maxn > 0) {
	  e -= Vinti(k0, k3) * Vinti(k1, k2) / am;
	}
      }
      if (qed.br < 0 || (maxn > 0 && maxn <= qed.br)) {
	e += Breit(k0, k1, k3, k2, t/2, kl0, kl1, kl3, kl2);
      }
      if (e) {
	e *= ReducedCL(js[0], t, js[3]); 
	e *= ReducedCL(js[1], t, js[2]);
	e *= a * (k + 1.0);
	if (IsOdd(t/2 + kk)) e = -e;
	*se += e;
      }
    }
  }

 EXIT:
#ifdef PERFORM_STATISTICS 
    stop = clock();
    rad_timing.radial_slater += stop - start;
#endif
  return 0;
}

double SelfEnergyRatio(ORBITAL *orb) {
  int i, npts;
  double *p, *q, e, z;
  double *large, *small;
  double a, b;
  
  if (orb->wfun == NULL) return 1.0;

  for (i = 0; i < potential->maxrp; i++) {
    if (potential->uehling[i] > -EPS4) break;
  }
  npts = i;
  
  p = _xk;
  q = _zk;
  z = potential->Z[potential->maxrp-1];
  e = RadialDiracCoulomb(npts, p, q, potential->rad, z, 
			 orb->n, orb->kappa);
  large = Large(orb);
  small = Small(orb);  
  for (i = 0; i < npts; i++) {
    p[i] = (p[i]*p[i] + q[i]*q[i])*potential->dr_drho[i];
    p[i] *= potential->uehling[i];
    q[i] = (large[i]*large[i] + small[i]*small[i])*potential->dr_drho[i];
    q[i] *= potential->uehling[i];
  }
  a = Simpson(p, 0, npts-1);
  b = Simpson(q, 0, npts-1);

  return b/a;
}

double QED1E(int k0, int k1) {
  int i;
  ORBITAL *orb1, *orb2;
  int index[2];
  double *p, r, a;

  if (qed.nms == 0 && qed.vp == 0) {
    if (qed.se == 0 || k0 != k1) {
      return 0.0;
    }
  }

  orb1 = GetOrbitalSolved(k0);
  orb2 = GetOrbitalSolved(k1);

  if (orb1->n <= 0 || orb2->n <= 0) {
    return 0.0;
  }
  if (orb1->wfun == NULL || orb2->wfun == NULL) {
    return 0.0;
  }
  
  if (k0 > k1) {
    index[0] = k1;
    index[1] = k0;
  } else {
    index[0] = k0;
    index[1] = k1;
  }
  
  p = (double *) MultiSet(qed1e_array, index, NULL, InitDoubleData, NULL);
  if (p && *p) {
    return *p;
  }

  r = 0.0;
  
  if (qed.nms > 0) {
    for (i = 0; i < potential->maxrp; i++) {
      _yk[i] = potential->U[i] + potential->Vc[i];
    }
    a = 0.0;
    Integrate(_yk, orb1, orb2, 1, &a, -1);
    a = -a;
    if (k0 == k1) a += orb1->energy;
    a /= (AMU * GetAtomicMass());
    r += a;
    for (i = 0; i < potential->maxrp; i++) {
      _yk[i] = orb1->energy - (potential->U[i] + potential->Vc[i]);
      _yk[i] *= orb2->energy - (potential->U[i] + potential->Vc[i]);
    }
    Integrate(_yk, orb1, orb2, 1, &a, -1);
    a *= FINE_STRUCTURE_CONST2/(2.0 * AMU * GetAtomicMass());
    r += a;
  }

  if (qed.vp > 0) {
    a = 0.0;
    Integrate(potential->uehling, orb1, orb2, 1, &a, 0);
    r += a;
  }

  if (k0 == k1 && (qed.se < 0 || orb1->n <= qed.se)) {
    if (potential->ib <= 0 || orb1->n <= potential->nb) {
      a = HydrogenicSelfEnergy(potential->Z[potential->maxrp-1], 
			       orb1->n, orb1->kappa);
      if (a) {
	a *= SelfEnergyRatio(orb1);
	r += a;
      }
    }
  }
  *p = r;
  return r;
}
  
double Vinti(int k0, int k1) {
  int i;
  ORBITAL *orb1, *orb2;
  int index[2];
  double *p, *large0, *small0, *large1, *small1;
  int ka0, ka1;
  double a, b, r;

  orb1 = GetOrbitalSolved(k0);
  orb2 = GetOrbitalSolved(k1);

  if (orb1->n <= 0 || orb2->n <= 0) {
    return 0.0;
  }
  if (orb1->wfun == NULL || orb2->wfun == NULL) {
    return 0.0;
  }
  
  index[0] = k0;
  index[1] = k1;
  
  p = (double *) MultiSet(vinti_array, index, NULL, InitDoubleData, NULL);
  if (p && *p) {
    return *p;
  }

  ka0 = orb1->kappa;
  ka1 = orb2->kappa;
  large0 = Large(orb1);
  large1 = Large(orb2);
  small0 = Small(orb1);
  small1 = Small(orb2);
  a = 0.5*(ka0*(ka0+1.0) - ka1*(ka1+1.0));
  b = 0.5*(-ka0*(-ka0+1.0) + ka1*(-ka1+1.0));
  r = 0.0;

  Differential(large1, _zk, 0, potential->maxrp-1);
  for (i = 0; i < potential->maxrp; i++) {
    _yk[i] = large0[i]*_zk[i] - a*large0[i]*large1[i]/potential->rad[i];
    _yk[i] *= potential->dr_drho[i];
  }
  r += Simpson(_yk, 0, potential->maxrp-1);
  
  Differential(small1, _zk, 0, potential->maxrp-1);
  for (i = 0; i < potential->maxrp; i++) {
    _yk[i] = small0[i]*_zk[i] - b*small0[i]*small1[i]/potential->rad[i];
    _yk[i] *= potential->dr_drho[i];
  }
  r += Simpson(_yk, 0, potential->maxrp-1);
  
  *p = r;

  return r;
}

double BreitC(int n, int m, int k, int k0, int k1, int k2, int k3) {
  int ka0, ka1, ka2, ka3, kb, kp;
  double r, b, c;
  
  ka0 = GetOrbital(k0)->kappa;
  ka1 = GetOrbital(k1)->kappa;
  ka2 = GetOrbital(k2)->kappa;
  ka3 = GetOrbital(k3)->kappa;
  if (k == m) {
    r = -(ka0 + ka2)*(ka1 + ka3);
    if (r) r /= (m*(m+1.0));
  } else if (k == (m + 1)) {
    kb = ka2 - ka0;
    kp = ka3 - ka1;
    b = (m + 2.0)/(2.0*(2.0*m + 1.0));
    c = -(m - 1.0)/((2.0*m+1.0)*(2.0*m+2.0));
    switch (n) {
    case 0:
      r = (k + kb)*(b + c*kp);
      break;
    case 1:
      r = (k + kp)*(b + c*kb);
      break;
    case 2:
      r = (k - kb)*(b - c*kp);
      break;
    case 3:
      r = (k - kp)*(b - c*kb);
      break;
    case 4:
      r = -(k + kb)*(b - c*kp);
      break;
    case 5:
      r = -(k - kp)*(b + c*kb);
      break;
    case 6:
      r = -(k - kb)*(b + c*kp);
      break;
    case 7:
      r = -(k + kp)*(b - c*kb);
      break;
    default:
      r = 0;
    }
  } else {
    kb = ka2 - ka0;
    kp = ka3 - ka1;
    b = (m - 1.0)/(2.0*(2.0*m + 1.0));
    c = (m + 2.0)/(2.0*m*(2.0*m + 1.0));
    switch (n) {
    case 0:
      r = (b + c*kb)*(kp - k - 1.0);
      break;
    case 1:
      r = (b + c*kp)*(kb - k - 1.0);
      break;
    case 2:
      r = (b - c*kb)*(-kp - k - 1.0);
      break;
    case 3:
      r = (b - c*kp)*(-kb - k - 1.0);
      break;
    case 4:
      r = -(b + c*kb)*(-kp - k - 1.0);
      break;
    case 5:
      r = -(b - c*kp)*(kb - k - 1.0);
      break;
    case 6:
      r = -(b - c*kb)*(kp - k - 1.0);
      break;
    case 7:
      r = -(b + c*kp)*(-kb - k - 1.0);
      break;
    default:
      r = 0;
    }
  }

  return r;
}

double BreitS(int k0, int k1, int k2, int k3, int k) {
  ORBITAL *orb0, *orb1, *orb2, *orb3;
  int index[5], i;
  double *p, r;
  
  index[0] = k0;
  index[1] = k1;
  index[2] = k2;
  index[3] = k3;
  index[4] = k;

  p = (double *) MultiSet(breit_array, index, NULL, InitDoubleData, NULL);
  if (p && *p) {
    r = *p;
  } else {
    orb0 = GetOrbitalSolved(k0);
    orb1 = GetOrbitalSolved(k1);
    orb2 = GetOrbitalSolved(k2);
    orb3 = GetOrbitalSolved(k3);
    if (!orb0 || !orb1 || !orb2 || !orb3) return 0.0;
    
    for (i = 0; i < potential->maxrp; i++) {
      _dwork1[i] = pow(potential->rad[i], k);
    }
    
    Integrate(_dwork1, orb0, orb1, -6, _zk, 0);
    
    for (i = 0; i < potential->maxrp; i++) {
      _zk[i] /= _dwork1[i]*potential->rad[i];
    }

    Integrate(_zk, orb2, orb3, 6, &r, 0);
    *p = r;
  }

  return r;
}

double BreitI(int n, int k0, int k1, int k2, int k3, int m) {
  double r;

  switch (n) {
  case 0:
    r = BreitS(k0, k2, k1, k3, m);
    break;
  case 1:
    r = BreitS(k1, k3, k0, k2, m);
    break;
  case 2:
    r = BreitS(k2, k0, k3, k1, m);
    break;
  case 3:
    r = BreitS(k3, k1, k2, k0, m);
    break;
  case 4:
    r = BreitS(k0, k2, k3, k1, m);
    break;
  case 5:
    r = BreitS(k3, k1, k0, k2, m);
    break;
  case 6:
    r = BreitS(k2, k0, k1, k3, m);
    break;
  case 7:
    r = BreitS(k1, k3, k2, k0, m);
    break;
  default:
    r = 0.0;
  }

  return r;
}

double Breit(int k0, int k1, int k2, int k3, int k,
	     int kl0, int kl1, int kl2, int kl3) {
  int m, m0, m1, n;
  double a, c, r;
  
  m0 = k - 1;
  if (m0 < 0) m0 = 0;
  m1 = k + 1;
  r = 0.0;
  for (m = m0; m <= m1; m++) {
    if (IsEven((kl0+kl2)/2 + m) || IsEven((kl1+kl3)/2 + m)) continue;
    for (n = 0; n < 8; n++) {
      c = BreitC(n, m, k, k0, k1, k2, k3);
      a = BreitI(n, k0, k1, k2, k3, m);
      r += a*c;
    }
  }

  return r;
}

/* calculate the slater integral of rank k */
int Slater(double *s, int k0, int k1, int k2, int k3, int k, int mode) {
  int index[5];
  double *p;
  int ilast, i, npts, m;
  ORBITAL *orb0, *orb1, *orb2, *orb3;
  double norm;
#ifdef PERFORM_STATISTICS
  clock_t start, stop; 
  start = clock();
#endif

  index[0] = k0;
  index[1] = k1;
  index[2] = k2;
  index[3] = k3;
  index[4] = k;  

  if (abs(mode) < 2) {
    SortSlaterKey(index);
    p = (double *) MultiSet(slater_array, index, NULL, InitDoubleData, NULL);
  } else {
    p = NULL;
  }
  if (p && *p) {
    *s = *p;
  } else {
    orb0 = GetOrbitalSolved(k0);
    orb1 = GetOrbitalSolved(k1);
    orb2 = GetOrbitalSolved(k2);
    orb3 = GetOrbitalSolved(k3);
    *s = 0.0;
    if (!orb0 || !orb1 || !orb2 || !orb3) return -1;  

    npts = potential->maxrp;
    switch (mode) {
    case 0: /* fall through to case 1 */
    case 1: /* full relativistic with distorted free orbitals */
      GetYk(k, _yk, orb0, orb2, k0, k2, -1); 
      if (orb1->n > 0) ilast = orb1->ilast;
      else ilast = npts-1;
      if (orb3->n > 0) ilast = Min(ilast, orb3->ilast);
      for (i = 0; i <= ilast; i++) {
	_yk[i] = (_yk[i]/potential->rad[i]);
      }
      Integrate(_yk, orb1, orb3, 1, s, 0);
      break;
    
    case -1: /* quasi relativistic with distorted free orbitals */
      GetYk(k, _yk, orb0, orb2, k0, k2, -2);
      if (orb1->n > 0) ilast = orb1->ilast;
      else ilast = npts-1;
      if (orb3->n > 0) ilast = Min(ilast, orb3->ilast);
      for (i = 0; i <= ilast; i++) {
	_yk[i] /= potential->rad[i];
      }
      Integrate(_yk, orb1, orb3, 2, s, 0);

      norm  = orb0->qr_norm;
      norm *= orb1->qr_norm;
      norm *= orb2->qr_norm;
      norm *= orb3->qr_norm;

      *s *= norm;
      break;

    case 2: /* separable coulomb interaction, orb0, orb2 is inner part */
      m = k;
      *s = RadialMoments(m, k0, k2);
      if (*s != 0.0) {
	m = -m-1;
	*s *= RadialMoments(m, k1, k3);
      }
      break;

    case -2: /* separable coulomb interaction, orb1, orb3 is inner part  */
      m = k;
      *s = RadialMoments(m, k1, k3);
      if (*s != 0.0) {
	m = -m-1;
	*s *= RadialMoments(m, k0, k2);      
      }
      break;
      
    default:
      break;
    }      

    if (p) *p = *s;
  }
#ifdef PERFORM_STATISTICS 
    stop = clock();
    rad_timing.radial_2e += stop - start;
#endif
  return 0;
}


/* reorder the orbital index appears in the slater integral, so that it is
   in a form: a <= b <= d, a <= c, and if (a == b), c <= d. */ 
void SortSlaterKey(int *kd) {
  int i;

  if (kd[0] > kd[2]) {
    i = kd[0];
    kd[0] = kd[2];
    kd[2] = i;
  }

  if (kd[1] > kd[3]) {
    i = kd[1];
    kd[1] = kd[3];
    kd[3] = i;
  }
  
  if (kd[0] > kd[1]) {
    i = kd[0];
    kd[0] = kd[1];
    kd[1] = i;
    i = kd[2];
    kd[2] = kd[3];
    kd[3] = i;
  } else if (kd[0] == kd[1]) {
    if (kd[2] > kd[3]) {
      i = kd[2];
      kd[2] = kd[3];
      kd[3] = i;
    }
  }
}

void PrepSlater(int ib0, int iu0, int ib1, int iu1,
		int ib2, int iu2, int ib3, int iu3) {
  int k, kmax, kk, i, j, p, q, m, ilast;
  int j0, j1, j2, j3, k0, k1, k2, k3;
  int index[6];
  double *dp;
  ORBITAL *orb0, *orb1, *orb2, *orb3;
  int c = 0;

  kmax = GetMaxRank();
  for (kk = 0; kk <= kmax; kk += 2) {
    k = kk/2;
    for (i = ib0; i <= iu0; i++) {
      orb0 = GetOrbital(i);
      GetJLFromKappa(orb0->kappa, &j0, &k0);
      for (p = ib2; p <= iu2; p++) {
	if (p < i) continue;
	orb2 = GetOrbital(p);
	GetJLFromKappa(orb2->kappa, &j2, &k2);
	if (k0 > slater_cut.kl0 || k2 > slater_cut.kl0) continue;
	GetYk(k, _yk, orb0, orb2, i, p, -1);
	ilast = potential->maxrp-1;
	for (m = 0; m <= ilast; m++) {
	  _yk[m] /= potential->rad[m];
	}
	for (j = ib1; j <= iu1; j++) {
	  if (j < i) continue;
	  orb1 = GetOrbital(j);
	  GetJLFromKappa(orb1->kappa, &j1, &k1);
	  for (q = ib3; q <= iu3; q++) {
	    if (q < j) continue;
	    orb3 = GetOrbital(q);
	    GetJLFromKappa(orb3->kappa, &j3, &k3);
	    if (i == j && q < p) continue;
	    if (IsOdd((k0+k2)/2+k) || 
		IsOdd((k1+k3)/2+k) ||
		!Triangle(j0, j2, kk) ||
		!Triangle(j1, j3, kk)) continue;	     
	    index[0] = i;
	    index[1] = j;
	    index[2] = p;
	    index[3] = q;
	    index[4] = k;
	    index[5] = 0;
	    dp = MultiSet(slater_array, index, NULL, InitDoubleData, NULL);
	    c++;
	    if (*dp == 0) {
	      Integrate(_yk, orb1, orb3, 1, dp, 0);
	    }
	  }
	}
      }
    }
  }
  printf("PrepSlater: %d\n", c);
}
      
int GetYk0(int k, double *yk, ORBITAL *orb1, ORBITAL *orb2, int type) {
  int i, ilast, i0;
  double a, max;

  for (i = 0; i < potential->maxrp; i++) {
    _dwork1[i] = pow(potential->rad[i], k);
  }
  Integrate(_dwork1, orb1, orb2, type, _zk, 0);
  for (i = 0; i < potential->maxrp; i++) {
    _zk[i] /= _dwork1[i];
    yk[i] = _zk[i];
    _zk[i] = _dwork1[i];
  }
  if (k > 2) {
    max = 0.0;
    for (i = 0; i < potential->maxrp; i++) {
      a = fabs(yk[i]);
      if (max < a) max = a;
    }
    max *= 1E-3;
    ilast = Min(orb1->ilast, orb2->ilast);
    for (i = 0; i < ilast; i++) {
      a = Large(orb1)[i]*Large(orb2)[i]*potential->rad[i];
      if (fabs(a) > max) break;
      _dwork1[i] = 0.0;
    }
    i0 = i;
  } else i0 = 0;
  for (i = i0; i < potential->maxrp; i++) {
    _dwork1[i] = pow(potential->rad[i0]/potential->rad[i], k+1);
  }
  Integrate(_dwork1, orb1, orb2, type, _xk, 0);
  ilast = potential->maxrp - 1;    
  for (i = i0; i < potential->maxrp; i++) {
    _xk[i] = (_xk[ilast] - _xk[i])/_dwork1[i];
    yk[i] += _xk[i];
  }

  return 0;
}  

/*
** this is a better version of Yk than GetYk0.
** note that on exit, _zk consists r^k, which is used in GetYk
*/      
int GetYk1(int k, double *yk, ORBITAL *orb1, ORBITAL *orb2, int type) {
  int i, ilast;
  double r0, a;
  
  ilast = Min(orb1->ilast, orb2->ilast);
  r0 = sqrt(potential->rad[0]*potential->rad[ilast]);  
  for (i = 0; i < potential->maxrp; i++) {
    _dwork1[i] = pow(potential->rad[i]/r0, k);
  }
  Integrate(_dwork1, orb1, orb2, type, _zk, 0);
  a = pow(r0, k);
  for (i = 0; i < potential->maxrp; i++) {
    _zk[i] /= _dwork1[i];
    yk[i] = _zk[i];
    _zk[i] = _dwork1[i]*a;
  }  
  for (i = 0; i < potential->maxrp; i++) {
    _dwork1[i] = (r0/potential->rad[i])/_dwork1[i];
  }
  Integrate(_dwork1, orb1, orb2, type, _xk, -1);
  for (i = 0; i < potential->maxrp; i++) {
    yk[i] += _xk[i]/_dwork1[i];
  }
      
  return 0;
}
      
int GetYk(int k, double *yk, ORBITAL *orb1, ORBITAL *orb2, 
	  int k1, int k2, int type) {
  int i, i0, i1, n;
  double a, b, a2, b2, max, max1;
  int index[3];
  SLATER_YK *syk;

  if (k1 <= k2) {
    index[0] = k1;
    index[1] = k2;
  } else {
    index[0] = k2;
    index[1] = k1;
  }
  index[2] = k;

  syk = (SLATER_YK *) MultiSet(yk_array, index, NULL, InitYkData, FreeYkData);
  if (syk->npts < 0) {
    GetYk1(k, yk, orb1, orb2, type);
    max = 0;
    for (i = 0; i < potential->maxrp; i++) {
      _zk[i] *= yk[i];
      a = fabs(_zk[i]); 
      if (a > max) max = a;
    }
    max1 = max*EPS5;
    max = max*EPS4;
    a = _zk[potential->maxrp-1];
    for (i = potential->maxrp-2; i >= 0; i--) {
      if (fabs(_zk[i] - a) > max1) {
	break;
      }
    }
    i1 = i;
    for (i = i1; i >= 0; i--) {      
      b = fabs(a - _zk[i]);
      _zk[i] = log(b);
      if (b > max) {
	break;
      }
    }
    i0 = i;
    if (i0 == i1) {
      i0--;
      b = fabs(a - _zk[i0]);
      _zk[i0] = log(b);
    }
    syk->coeff[0] = a;    
    syk->npts = i0+1;
    syk->yk = malloc(sizeof(float)*(syk->npts));
    for (i = 0; i < syk->npts ; i++) {
      syk->yk[i] = yk[i];
    }
    n = i1 - i0 + 1;
    a = 0.0;
    b = 0.0;
    a2 = 0.0;
    b2 = 0.0;
    for (i = i0; i <= i1; i++) {      
      max = (potential->rad[i]-potential->rad[i0]);
      a += max;
      b += _zk[i];
      a2 += max*max;
      b2 += _zk[i]*max;
    }
    syk->coeff[1] = (a*b - n*b2)/(a*a - n*a2);       
    if (syk->coeff[1] >= 0) {
      i1 = i0 + (i1-i0)*0.3;
      if (i1 == i0) i1 = i0 + 1;
      for (i = i0; i <= i1; i++) {      
	max = (potential->rad[i]-potential->rad[i0]);
	a += max;
	b += _zk[i];
	a2 += max*max;
	b2 += _zk[i]*max;
      }
      syk->coeff[1] = (a*b - n*b2)/(a*a - n*a2);  
    }
    if (syk->coeff[1] >= 0) {
      syk->coeff[1] = -10.0/(potential->rad[i1]-potential->rad[i0]);
    } 
  } else {
    for (i = syk->npts-1; i < potential->maxrp; i++) {
      _dwork1[i] = pow(potential->rad[i], k);
    }
    for (i = 0; i < syk->npts; i++) {
      yk[i] = syk->yk[i];
    }
    i0 = syk->npts-1;
    a = syk->yk[i0]*_dwork1[i0];
    for (i = syk->npts; i < potential->maxrp; i++) {
      b = potential->rad[i] - potential->rad[i0];
      b = syk->coeff[1]*b;
      if (b < -20) {
	yk[i] = syk->coeff[0];
      } else {
	yk[i] = (a - syk->coeff[0])*exp(b);
	yk[i] += syk->coeff[0];
      }
      yk[i] /= _dwork1[i];
    }    
  }
      
  return 0;
}  

/* integrate a function given by f with two orbitals. */
/* type indicates the type of integral */
/* type = 1,    P1*P2 + Q1*Q2 */
/* type = 2,    P1*P2 */
/* type = 3,    Q1*Q2 */ 
/* type = 4:    P1*Q2 + Q1*P2 */
/* type = 5:    P1*Q2 - Q1*P2 */
/* type = 6:    P1*Q2 */
/* if type is positive, only the end point is returned, */
/* otherwise, the whole function is returned */
/* id indicate whether integrate inward (-1) or outward (0) */
int Integrate(double *f, ORBITAL *orb1, ORBITAL *orb2, 
	      int t, double *x, int id) {
  int i1, i2, ilast;
  int i, type;
  double *r, ext;

  if (t == 0) t = 1;
  if (t < 0) {
    r = x;
    type = -t;
  } else {
    r = _dwork;
    type = t;
  }
  for (i = 0; i < potential->maxrp; i++) {
    r[i] = 0.0;
  }

  ext = 0.0;
  ilast = Min(orb1->ilast, orb2->ilast);
  if (id >= 0) {
    IntegrateSubRegion(0, ilast, f, orb1, orb2, t, r, 0, NULL);
  } else {
    IntegrateSubRegion(0, ilast, f, orb1, orb2, t, r, -1, NULL);
  }
  i2 = ilast;
  if (orb1->ilast == ilast && orb1->n == 0) {
    i1 = ilast + 1;
    i2 = orb2->ilast;
    if (i2 > i1) {
      if (type == 6) {
	i = 7;
	if (t < 0) i = -i;
	IntegrateSubRegion(i1, i2, f, orb2, orb1, i, r, 1, NULL);
      } else {
	IntegrateSubRegion(i1, i2, f, orb1, orb2, t, r, 2, NULL);
      }
      i2--;
    }
    if (orb2->n == 0) {
      i1 = orb2->ilast + 1;
      i2 = potential->maxrp - 1;
      IntegrateSubRegion(i1, i2, f, orb1, orb2, t, r, 3, &ext);
      i2--;
    }
  } else if (orb2->ilast == ilast && orb2->n == 0) {
    i1 = ilast + 1;
    i2 = orb1->ilast;
    if (i2 > i1) {
      IntegrateSubRegion(i1, i2, f, orb1, orb2, t, r, 1, NULL);
      i2--;
    }
    if (orb1->n == 0) {
      i1 = orb1->ilast + 1;
      i2 = potential->maxrp - 1;
      IntegrateSubRegion(i1, i2, f, orb1, orb2, t, r, 3, &ext);
      i2--;
    }
  }
  if (t >= 0) {
    *x = r[i2] + ext;
  } else {
    if (id >= 0) {
      for (i = i2+1; i < potential->maxrp; i++) {
	r[i] = r[i2];
      }
    } else {
      if (i2 > ilast) {
	ext += r[i2];
	for (i = ilast + 1; i <= i2; i++) {
	  r[i] = ext - r[i];
	}
	for (i = i2+1; i < potential->maxrp; i++) {
	  r[i] = 0.0;
	}
	for (i = 0; i <= ilast; i++) {
	  r[i] = r[ilast+1] + r[i];
	}
      }
    }
  }
  
  return 0;
}

void AddEvenPoints(double *r, double *r1, int i0, int i1, int t) {
  int i;

  if (t < 0) {
    for (i = i0; i < i1; i++) {
      r[i] += r1[i];
    }
  } else {
    r[i1-1] += r1[i1-1];
  }
}

int IntegrateSubRegion(int i0, int i1, 
		       double *f, ORBITAL *orb1, ORBITAL *orb2,
		       int t, double *r, int m, double *ext) {
  int i, j, ip, i2, type;
  ORBITAL *tmp;
  double *large1, *large2, *small1, *small2;
  double *x, *y, *r1, *x1, *x2, *y1, *y2;
  double a, b, e1, e2, a2, r0;

  if (i1 <= i0) return 0;
  type = abs(t);

  x = _dwork3;
  y = _dwork4;
  x1 = _dwork5;
  x2 = _dwork6;
  y1 = _dwork7;
  y2 = _dwork8;
  r1 = _dwork9;
  i2 = i1;

  switch (m) {
  case -1: /* m = -1 same as m = 0 but integrate inward */
  case 0:  /* m = 0 */
    large1 = Large(orb1);
    large2 = Large(orb2);
    small1 = Small(orb1);
    small2 = Small(orb2);
    switch (type) {
    case 1: /* type = 1 */
      for (i = i0; i <= i1; i++) {
	x[i] = large1[i] * large2[i];
	x[i] += small1[i] * small2[i];
	x[i] *= f[i]*potential->dr_drho[i];
      }
      if (i1 == orb1->ilast && orb1->n == 0 && i < potential->maxrp) {
	if (i <= orb2->ilast) {
	  ip = i+1;
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  x[i] = large1[i]*a*large2[i];
	  x[i] += (small1[i]*b + small1[ip]*a)*small2[i];
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	} else if (orb2->n == 0) {
	  ip = i+1;
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  e1 = sin(large2[ip]);
	  e2 = cos(large2[ip]);
	  x[i] = large1[i]*a*large2[i]*e1;
	  x[i] += (small1[i]*b+small1[ip]*a)*(small2[i]*e2+small2[ip]*e1);
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	}
      } else if (i1 == orb2->ilast && orb2->n == 0 && i < potential->maxrp) {
	if (i <= orb1->ilast) {
	  ip = i+1;
	  a = sin(large2[ip]);
	  b = cos(large2[ip]);
	  x[i] = large2[i]*a*large1[i];
	  x[i] += (small2[i]*b + small2[ip]*a)*small1[i];
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	} else if (orb1->n == 0) {
	  ip = i+1;
	  a = sin(large2[ip]);
	  b = cos(large2[ip]);
	  e1 = sin(large1[ip]);
	  e2 = cos(large1[ip]);
	  x[i] = large2[i]*a*large1[i]*e1;
	  x[i] += (small2[i]*b+small2[ip]*a)*(small1[i]*e2+small1[ip]*e1);
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	}
      }
      break;
    case 2: /* type = 2 */
      for (i = i0; i <= i1; i++) {
	x[i] = large1[i] * large2[i];
	x[i] *= f[i]*potential->dr_drho[i];
      }
      if (i1 == orb1->ilast && orb1->n == 0 && i < potential->maxrp) {
	if (i <= orb2->ilast) {
	  ip = i+1;
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  x[i] = large1[i]*a*large2[i];
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	} else if (orb2->n == 0) {
	  ip = i+1;
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  e1 = sin(large2[ip]);
	  e2 = cos(large2[ip]);
	  x[i] = large1[i]*a*large2[i]*e1;
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	}
      } else if (i1 == orb2->ilast && orb2->n == 0 && i < potential->maxrp) {
	if (i <= orb1->ilast) {
	  ip = i+1;
	  a = sin(large2[ip]);
	  b = cos(large2[ip]);
	  x[i] = large2[i]*a*large1[i];
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	} else if (orb1->n == 0) {
	  ip = i+1;
	  a = sin(large2[ip]);
	  b = cos(large2[ip]);
	  e1 = sin(large1[ip]);
	  e2 = cos(large1[ip]);
	  x[i] = large2[i]*a*large1[i]*e1;
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	}
      }
      break;
    case 3: /*type = 3 */
      for (i = i0; i <= i1; i++) {
	x[i] = small1[i] * small2[i];
	x[i] *= f[i]*potential->dr_drho[i];
      }
      if (i1 == orb1->ilast && orb1->n == 0 && i < potential->maxrp) {
	if (i <= orb2->ilast) {
	  ip = i+1;
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  x[i] = (small1[i]*b + small1[ip]*a)*small2[i];
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	} else if (orb2->n == 0) {
	  ip = i+1;
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  e1 = sin(large2[ip]);
	  e2 = cos(large2[ip]);
	  x[i] = (small1[i]*b+small1[ip]*a)*(small2[i]*e2+small2[ip]*e1);
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	}
      } else if (i1 == orb2->ilast && orb2->n == 0 && i < potential->maxrp) {
	if (i <= orb1->ilast) {
	  ip = i+1;
	  a = sin(large2[ip]);
	  b = cos(large2[ip]);
	  x[i] = (small2[i]*b + small2[ip]*a)*small1[i];
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	} else if (orb1->n == 0) {
	  ip = i+1;
	  a = sin(large2[ip]);
	  b = cos(large2[ip]);
	  e1 = sin(large1[ip]);
	  e2 = cos(large1[ip]);
	  x[i] = (small2[i]*b+small2[ip]*a)*(small1[i]*e2+small1[ip]*e1);
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	}
      }
      break;
    case 4: /*type = 4 */
      for (i = i0; i <= i1; i++) {
	x[i] = large1[i] * small2[i];
	x[i] += small1[i] * large2[i];
	x[i] *= f[i]*potential->dr_drho[i];
      }
      if (i1 == orb1->ilast && orb1->n == 0 && i < potential->maxrp) {
	if (i <= orb2->ilast) {
	  ip = i+1;
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  x[i] = large1[i]*a*small2[i];
	  x[i] += (small1[i]*b + small1[ip]*a)*large2[i];
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	} else if (orb2->n == 0) {
	  ip = i+1;
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  e1 = sin(large2[ip]);
	  e2 = cos(large2[ip]);
	  x[i] = large1[i]*a*(small2[i]*e2+small2[ip]*e1);
	  x[i] += (small1[i]*b+small1[ip]*a)*large2[i]*e1;
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	}
      } else if (i1 == orb2->ilast && orb2->n == 0 && i < potential->maxrp) {
	if (i <= orb1->ilast) {
	  ip = i+1;
	  a = sin(large2[ip]);
	  b = cos(large2[ip]);
	  x[i] = (small2[i]*b + small2[ip]*a)*large1[i];
	  x[i] += large2[i]*a*small1[i];
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	} else if (orb1->n == 0) {
	  ip = i+1;
	  a = sin(large2[ip]);
	  b = cos(large2[ip]);
	  e1 = sin(large1[ip]);
	  e2 = cos(large1[ip]);
	  x[i] = (small2[i]*b+small2[ip]*a)*large1[i]*e1;
	  x[i] += large2[i]*a*(small1[i]*e2+small1[ip]*e1);
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	}
      }
      break;
    case 5: /* type = 5 */
      for (i = i0; i <= i1; i++) {
	x[i] = large1[i] * small2[i];
	x[i] -= small1[i] * large2[i];
	x[i] *= f[i]*potential->dr_drho[i];
      } 
      if (i1 == orb1->ilast && orb1->n == 0 && i < potential->maxrp) {
	if (i <= orb2->ilast) {
	  ip = i+1;
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  x[i] = large1[i]*a*small2[i];
	  x[i] -= (small1[i]*b + small1[ip]*a)*large2[i];
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	} else if (orb2->n == 0) {
	  ip = i+1;
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  e1 = sin(large2[ip]);
	  e2 = cos(large2[ip]);
	  x[i] = large1[i]*a*(small2[i]*e2+small2[ip]*e1);
	  x[i] -= (small1[i]*b+small1[ip]*a)*large2[i]*e1;
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	}
      } else if (i1 == orb2->ilast && orb2->n == 0 && i < potential->maxrp) {
	if (i <= orb1->ilast) {
	  ip = i+1;
	  a = sin(large2[ip]);
	  b = cos(large2[ip]);
	  x[i] = (small2[i]*b + small2[ip]*a)*large1[i];
	  x[i] -= large2[i]*a*small1[i];
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	} else if (orb1->n == 0) {
	  ip = i+1;
	  a = sin(large2[ip]);
	  b = cos(large2[ip]);
	  e1 = sin(large1[ip]);
	  e2 = cos(large1[ip]);
	  x[i] = (small2[i]*b+small2[ip]*a)*large1[i]*e1;
	  x[i] -= large2[i]*a*(small1[i]*e2+small1[ip]*e1);
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	}
      }
      break;
    case 6: /* type = 6 */
      for (i = i0; i <= i1; i++) {
	x[i] = large1[i] * small2[i];
	x[i] *= f[i]*potential->dr_drho[i];
      }
      if (i1 == orb1->ilast && orb1->n == 0 && i < potential->maxrp) {
	if (i <= orb2->ilast) {
	  ip = i+1;
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  x[i] = large1[i]*a*small2[i];
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	} else if (orb2->n == 0) {
	  ip = i+1;
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  e1 = sin(large2[ip]);
	  e2 = cos(large2[ip]);
	  x[i] = large1[i]*a*(small2[i]*e2+small2[ip]*e1);
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	}
      } else if (i1 == orb2->ilast && orb2->n == 0 && i < potential->maxrp) {
	if (i <= orb1->ilast) {
	  ip = i+1;
	  a = sin(large2[ip]);
	  b = cos(large2[ip]);
	  x[i] = (small2[i]*b + small2[ip]*a)*large1[i];
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	} else if (orb1->n == 0) {
	  ip = i+1;
	  a = sin(large2[ip]);
	  b = cos(large2[ip]);
	  e1 = sin(large1[ip]);
	  e2 = cos(large1[ip]);
	  x[i] = (small2[i]*b+small2[ip]*a)*large1[i]*e1;
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	}
      }
      break;
    default: /* error */
      return -1;
    }      
    NewtonCotes(r, x, i0, i2, t, m);
    break;

  case 1: /* m = 1 */
    if (type == 6) { /* type 6 needs special treatments */
      large1 = Large(orb1);
      large2 = Large(orb2);
      small1 = Small(orb1);
      small2 = Small(orb2);
      e1 = orb1->energy;
      e2 = orb2->energy;
      a2 = 0.5*FINE_STRUCTURE_CONST2;
      j = 0;
      for (i = i0; i <= i1; i+= 2) {
	ip = i+1;
	x[j] = large1[i] * small2[ip];
	x[j] *= f[i];
	y[j] = large1[i] * small2[i] * f[i];
	_phase[j] = large2[ip];
	_dphase[j] = (1.0 + a2*(e2-potential->U[i]-potential->Vc[i]))
	  /(large2[i]*large2[i]);
	j++;
      }
      if (i < potential->maxrp) {
	ip = i+1;
	if (i > orb1->ilast && orb1->n == 0) { 
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  b = small1[i]*b + small1[ip]*a;
	  a = large1[i]*a;
	  x[j] = a * small2[ip];
	  x[j] *= f[i];
	  y[j] = a * small2[i] * f[i];
	  _phase[j] = large2[ip];
	  _dphase[j] = (1.0 + a2*(e2-potential->U[i]-potential->Vc[i]))
	    /(large2[i]*large2[i]);
	  j++;
	  i2 = i;
	}
      }
      IntegrateSinCos(j, x, y, _phase, _dphase, i0, r, t, ext);
      if (IsOdd(i2)) r[i2] = r[i2-1];
      break;
    } else if (type == 7) {
      large1 = Large(orb1);
      large2 = Large(orb2);
      small1 = Small(orb1);
      small2 = Small(orb2);
      e1 = orb1->energy;
      e2 = orb2->energy;
      a2 = 0.5*FINE_STRUCTURE_CONST2;
      j = 0;
      for (i = i0; i <= i1; i+= 2) {
	ip = i+1;
	x[j] = small1[i] * large2[i];
	x[j] *= f[i];
	_phase[j] = large2[ip];
	_dphase[j] = (1.0 + a2*(e2-potential->U[i]-potential->Vc[i]))
	  /(large2[i]*large2[i]);
	j++;	
      }
      if (i < potential->maxrp) {
	ip = i+1;
	if (i > orb1->ilast && orb1->n == 0) { 
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  b = small1[i]*b + small1[ip]*a;
	  x[j] = b * large2[i];
	  x[j] *= f[i];
	  _phase[j] = large2[ip];
	  _dphase[j] = (1.0 + a2*(e2-potential->U[i]-potential->Vc[i]))
	    /(large2[i]*large2[i]);
	  j++;
	  i2 = i;
	}
      }
      IntegrateSinCos(j, x, NULL, _phase, _dphase, i0, r, t, ext);
      if (IsOdd(i2)) r[i2] = r[i2-1];
      break;
    }
  case 2: /* m = 1, 2 are essentially the same */
    /* type 6 is treated in m = 1 */
    if (m == 2) {
      tmp = orb1;
      orb1 = orb2;
      orb2 = tmp;
    }
    large1 = Large(orb1);
    large2 = Large(orb2);
    small1 = Small(orb1);
    small2 = Small(orb2);
    e1 = orb1->energy;
    e2 = orb2->energy;
    a2 = 0.5*FINE_STRUCTURE_CONST2;
    switch (type) {
    case 1: /* type = 1 */
      j = 0;
      for (i = i0; i <= i1; i+= 2) {
	ip = i+1;
	x[j] = large1[i] * large2[i];
	x[j] += small1[i] * small2[ip];
	x[j] *= f[i];
	y[j] = small1[i] * small2[i] * f[i];
	_phase[j] = large2[ip];
	_dphase[j] = (1.0 + a2*(e2-potential->U[i]-potential->Vc[i]))
	  /(large2[i]*large2[i]);
	j++;
      }
      if (i < potential->maxrp) {
	ip = i+1;
	if (i > orb1->ilast && orb1->n == 0) {
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  b = small1[i]*b + small1[ip]*a;
	  a = large1[i]*a;
	  x[j] = a * large2[i];
	  x[j] += b * small2[ip];
	  x[j] *= f[i];
	  y[j] = b * small2[i] * f[i];
	  _phase[j] = large2[ip];
	  _dphase[j] = (1.0 + a2*(e2-potential->U[i]-potential->Vc[i]))
	    /(large2[i]*large2[i]);
	  j++;
	  i2 = i;
	}
      }
      IntegrateSinCos(j, x, y, _phase, _dphase, i0, r, t, ext);
      if (IsOdd(i2)) r[i2] = r[i2-1];
      break;
    case 2: /* type = 2 */
      j = 0;
      for (i = i0; i <= i1; i+= 2) {
	ip = i+1;
	x[j] = large1[i] * large2[i];
	x[j] *= f[i]; 
	_phase[j] = large2[ip];
	_dphase[j] = (1.0 + a2*(e2-potential->U[i]-potential->Vc[i]))
	  /(large2[i]*large2[i]);
	j++;
      }
      if (i < potential->maxrp) {
	ip = i+1;
	if (i > orb1->ilast && orb1->n == 0) {
	  a = sin(large1[ip]);
	  a = large1[i]*a;
	  x[j] = a * large2[i];
	  x[j] *= f[i];
	  _phase[j] = large2[ip];
	  _dphase[j] = (1.0 + a2*(e2-potential->U[i]-potential->Vc[i]))
	    /(large2[i]*large2[i]);
	  j++;
	  i2 = i;
	}
      }
      IntegrateSinCos(j, x, NULL, _phase, _dphase, i0, r, t, ext);
      if (IsOdd(i2)) r[i2] = r[i2-1];
      break;
    case 3: /* type = 3 */
      j = 0;
      for (i = i0; i <= i1; i+= 2) {
	ip = i+1;
	x[j] = small1[i] * small2[ip];
	x[j] *= f[i];
	y[j] = small1[i] * small2[i] * f[i];
	_phase[j] = large2[ip];
	_dphase[j] = (1.0 + a2*(e2-potential->U[i]-potential->Vc[i]))
	  /(large2[i]*large2[i]);
	j++;
      }
      if (i < potential->maxrp) {
	ip = i+1;
	if (i > orb1->ilast && orb1->n == 0) {
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  b = small1[i]*b + small1[ip]*a;
	  x[j] += b * small2[ip];
	  x[j] *= f[i];
	  y[j] = b * small2[i] * f[i];
	  _phase[j] = large2[ip];
	  _dphase[j] = (1.0 + a2*(e2-potential->U[i]-potential->Vc[i]))
	    /(large2[i]*large2[i]);
	  j++;
	  i2 = i;
	}
      }
      IntegrateSinCos(j, x, y, _phase, _dphase, i0, r, t, ext);
      if (IsOdd(i2)) r[i2] = r[i2-1];
      break;  
    case 4: /* type = 4 */
      j = 0;
      for (i = i0; i <= i1; i+= 2) {
	ip = i+1;
	x[j] = large1[i] * small2[ip];
	x[j] += small1[i] * large2[i];
	x[j] *= f[i];
	y[j] = large1[i] * small2[i] * f[i];
	_phase[j] = large2[ip];
	_dphase[j] = (1.0 + a2*(e2-potential->U[i]-potential->Vc[i]))
	  /(large2[i]*large2[i]);
	j++;
      }
      if (i < potential->maxrp) {
	ip = i+1;
	if (i > orb1->ilast && orb1->n == 0) {
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  b = small1[i]*b + small1[ip]*a;
	  a = large1[i]*a;
	  x[j] = a * small2[ip];
	  x[j] += b * large2[i];
	  x[j] *= f[i];
	  y[j] = a * small2[i] * f[i];
	  _phase[j] = large2[ip];
	  _dphase[j] = (1.0 + a2*(e2-potential->U[i]-potential->Vc[i]))
	    /(large2[i]*large2[i]);
	  j++;
	  i2 = i;
	}
      }
      IntegrateSinCos(j, x, y, _phase, _dphase, i0, r, t, ext);
      if (IsOdd(i2)) r[i2] = r[i2-1];
      break;
    case 5: /* type = 5 */
      j = 0;
      if (m == 2) {
	r0 = r[i0];
	r[i0] = 0.0;
      }
      for (i = i0; i <= i1; i+= 2) {
	ip = i+1;
	x[j] = large1[i] * small2[ip];
	x[j] -= small1[i] * large2[i];
	x[j] *= f[i];
	y[j] = large1[i] * small2[i] * f[i];
	_phase[j] = large2[ip];
	_dphase[j] = (1.0 + a2*(e2-potential->U[i]-potential->Vc[i]))
	  /(large2[i]*large2[i]);
	j++;
      }
      if (i < potential->maxrp) {
	ip = i+1;
	if (i > orb1->ilast && orb1->n == 0) { 
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  b = small1[i]*b + small1[ip]*a;
	  a = large1[i]*a;
	  x[j] = a * small2[ip];
	  x[j] -= b * large2[i];
	  x[j] *= f[i];
	  y[j] = a * small2[i] * f[i];
	  _phase[j] = large2[ip];
	  _dphase[j] = (1.0 + a2*(e2-potential->U[i]-potential->Vc[i]))
	    /(large2[i]*large2[i]);
	  j++;
	  i2 = i;
	}
      }
      IntegrateSinCos(j, x, y, _phase, _dphase, i0, r, t, ext);
      if (m == 2) {
	if (t < 0) {
	  for (i = i0; i <= i2; i++) {
	    r[i] = r0 - r[i];
	  }
	} else {
	  if (IsOdd(i2)) r[i2-1] = r0 - r[i2-1];
	  else r[i2] = r0 - r[i2];
	}
      }
      if (IsOdd(i2)) r[i2] = r[i2-1];
      break;
    default:
      return -1;
    }
    break;

  case 3: /* m = 3 */
    large1 = Large(orb1);
    large2 = Large(orb2);
    small1 = Small(orb1);
    small2 = Small(orb2);
    e1 = orb1->energy;
    e2 = orb2->energy;
    a2 = 0.5*FINE_STRUCTURE_CONST2;
    r1[i0] = 0.0;
    switch (type) {
    case 1: /* type = 1 */
      j = 0;
      for (i = i0; i <= i1; i += 2) {
	ip = i+1;	
	x1[j] = small1[i] * small2[ip];
	x2[j] = small1[ip] * small2[i];	
	x[j] = x1[j] + x2[j];
	x[j] *= 0.5*f[i];
	y1[j] = small1[i] * small2[i];
	y2[j] = small1[ip] * small2[ip];
	y2[j] += large1[i] * large2[i];
	y[j] = y1[j] - y2[j];
	y[j] *= 0.5*f[i];
	_phase[j] = large1[ip] + large2[ip];
	a = (1.0 + a2*(e1-potential->U[i]-potential->Vc[i]))
	  /(large1[i]*large1[i]);
	b = (1.0 + a2*(e2-potential->U[i]-potential->Vc[i]))
	  /(large2[i]*large2[i]);
	_dphase[j] = a + b;
	_dphasep[j] = a - b;
	j++;
      }
      IntegrateSinCos(j, x, y, _phase, _dphase, i0, r, t, ext);
      j = 0;
      for (i = i0; i <= i1; i += 2) {
	ip = i+1;
	x[j] = -x1[j] + x2[j];
	x[j] *= 0.5*f[i];
	y[j] = y1[j] + y2[j];
	y[j] *= 0.5*f[i];
	_phase[j] = large1[ip] - large2[ip];
	j++;
      }
      IntegrateSinCos(j, NULL, y, _phase, _dphasep, i0, r1, t, NULL);
      AddEvenPoints(r, r1, i0, i1, t);
      break;	
    case 2: /* type = 2 */
      j = 0;
      for (i = i0; i <= i1; i += 2) {
	ip = i+1;
	y2[j] = -large1[i] * large2[i];
	y2[j] *= 0.5*f[i];
	y[j] = y2[j];
	_phase[j] = large1[ip] + large2[ip];
	a = (1.0 + a2*(e1-potential->U[i]-potential->Vc[i]))
	  /(large1[i]*large1[i]);
	b = (1.0 + a2*(e2-potential->U[i]-potential->Vc[i]))
	  /(large2[i]*large2[i]);
	_dphase[j] = a + b;
	_dphasep[j] = a - b;
	j++;
      } 
      IntegrateSinCos(j, NULL, y, _phase, _dphase, i0, r, t, ext);
      j = 0;
      for (i = i0; i <= i1; i += 2) {
	ip = i+1;
	y[j] = -y2[j];
	_phase[j] = large1[ip] - large2[ip];
	j++;
      }
      IntegrateSinCos(j, NULL, y, _phase, _dphasep, i0, r1, t, NULL);
      AddEvenPoints(r, r1, i0, i1, t);
      break;
    case 3: /* type = 3 */
      j = 0;
      for (i = i0; i <= i1; i += 2) {
	ip = i+1;
	x1[j] = small1[i] * small2[ip];
	x2[j] = small1[ip] * small2[i];
	x[j] = x1[j] + x2[j];
	x[j] *= 0.5*f[i];
	y1[j] = small1[i] * small2[i];
	y2[j] = small1[ip] * small2[ip];
	y[j] = y1[j] - y2[j];
	y[j] *= 0.5*f[i];
	_phase[j] = large1[ip] + large2[ip];
	a = (1.0 + a2*(e1-potential->U[i]-potential->Vc[i]))
	  /(large1[i]*large1[i]);
	b = (1.0 + a2*(e2-potential->U[i]-potential->Vc[i]))
	  /(large2[i]*large2[i]);
	_dphase[j] = a + b;
	_dphasep[j] = a - b;
	j++;
      }
      IntegrateSinCos(j, x, y, _phase, _dphase, i0, r, t, ext);
      j = 0;
      for (i = i0; i <= i1; i += 2) {
	ip = i+1;
	x[j] = -x1[j] + x2[j];
	x[j] *= 0.5*f[i];
	y[j] = y1[j] + y2[j];
	y[j] *= 0.5*f[i];
	_phase[j] = large1[ip] - large2[ip];
	j++;
      }
      IntegrateSinCos(j, x, y, _phase, _dphasep, i0, r1, t, NULL);
      AddEvenPoints(r, r1, i0, i1, t);
      break;
    case 4: /* type = 4 */
      j = 0;
      for (i = i0; i <= i1; i += 2) {
	ip = i+1;
	x1[j] = large1[i] * small2[i];
	x2[j] = small1[i] * large2[i];
	x[j] = x1[j] + x2[j];
	x[j] *= 0.5*f[i];
	y[j] = -large1[i] * small2[ip];
	y[j] -= small1[ip] * large2[i];
	y[j] *= 0.5*f[i];
	y2[j] = y[j];
	_phase[j] = large1[ip] + large2[ip];	
	a = (1.0 + a2*(e1-potential->U[i]-potential->Vc[i]))
	  /(large1[i]*large1[i]);
	b = (1.0 + a2*(e2-potential->U[i]-potential->Vc[i]))
	  /(large2[i]*large2[i]);
	_dphase[j] = a + b;
	_dphasep[j] = a - b;
	j++;
      }
      IntegrateSinCos(j, x, y, _phase, _dphase, i0, r, t, ext);
      j = 0;
      for (i = i0; i <= i1; i += 2) {
	ip = i+1;
	x[j] = x1[j] - x2[j];
	x[j] *= 0.5*f[i];
	y[j] = -y2[j];
	_phase[j] = large1[ip] - large2[ip];
	j++;
      }
      IntegrateSinCos(j, x, y, _phase, _dphasep, i0, r1, t, NULL);
      AddEvenPoints(r, r1, i0, i1, t);
      break;
    case 5: /* type = 5 */
      j = 0;
      for (i = i0; i <= i1; i += 2) {
	ip = i+1;
	x1[j] = large1[i] * small2[i];
	x2[j] = small1[i] * large2[i];
	x[j] = x1[j] - x2[j];
	x[j] *= 0.5*f[i];
	y[j] = -large1[i] * small2[ip];
	y[j] += small1[ip] * large2[i];
	y[j] *= 0.5*f[i];
	y2[j] = y[j];
	_phase[j] = large1[ip] + large2[ip];
	a = (1.0 + a2*(e1-potential->U[i]-potential->Vc[i]))
	  /(large1[i]*large1[i]);
	b = (1.0 + a2*(e2-potential->U[i]-potential->Vc[i]))
	  /(large2[i]*large2[i]);
	_dphase[j] = a + b;
	_dphasep[j] = a - b;
	j++;
      }
      IntegrateSinCos(j, x, y, _phase, _dphase, i0, r, t, ext);
      j = 0;
      for (i = i0; i <= i1; i += 2) {
	ip = i+1;
	x[j] = x1[j] + x2[j];
	x[j] *= 0.5*f[i];
	y[j] = -y2[j];
	_phase[j] = large1[ip] - large2[ip];
	j++;
      }
      IntegrateSinCos(j, x, y, _phase, _dphasep, i0, r1, t, NULL);
      AddEvenPoints(r, r1, i0, i1, t);
      break;
    case 6: /* type = 6 */
      j = 0;
      for (i = i0; i <= i1; i += 2) {
	ip = i+1;
	x1[j] = large1[i] * small2[i];
	x1[j] *= 0.5*f[i];
	x[j] = x1[j];
	y[j] = -large1[i] * small2[ip];
	y[j] *= 0.5*f[i];
	y2[j] = y[j];
	_phase[j] = large1[ip] + large2[ip];
	a = (1.0 + a2*(e1-potential->U[i]-potential->Vc[i]))
	  /(large1[i]*large1[i]);
	b = (1.0 + a2*(e2-potential->U[i]-potential->Vc[i]))
	  /(large2[i]*large2[i]);
	_dphase[j] = a + b;
	_dphasep[j] = a - b;
	j++;
      }
      IntegrateSinCos(j, x, y, _phase, _dphase, i0, r, t, ext);
      j = 0;
      for (i = i0; i <= i1; i += 2) {
	ip = i+1;
	x[j] = x1[j];
	y[j] = -y2[j];
	_phase[j] = large1[ip] - large2[ip];
	j++;
      }
      IntegrateSinCos(j, x, y, _phase, _dphasep, i0, r1, t, NULL);
      AddEvenPoints(r, r1, i0, i1, t);
      break;
    default:
      return -1;
    }
  default:
    return -1; 
  }

  return 0;
}

int IntegrateSinCos(int j, double *x, double *y, 
		    double *phase, double *dphase, 
		    int i0, double *r, int t, double *ext) {
  int i, k, m, n, q, s, i1;
  double si0, si1, cs0, cs1;
  double is0, is1, is2, is3;
  double ic0, ic1, ic2, ic3;
  double his0, his1, his2, his3;
  double hic0, hic1, hic2, hic3;
  double a0, a1, a2, a3;
  double d, p, h, dr;
  double *z, *u, *w, *u1, *w1;

  w = _dwork2;
  z = _dwork10;
  u = _dwork11;
  if (phase[j-1] < 0) {
    for (i = 0; i < j; i++) {
      phase[i] = -phase[i];
      dphase[i] = -dphase[i];
      if (x) x[i] = -x[i];
    }
  }
  for (i = 1, k = i0+2; i < j; i++, k += 2) {
    h = phase[i] - phase[i-1];
    z[k] = 0.0;
    if (x != NULL) z[k] += x[i]*sin(phase[i]);
    if (y != NULL) z[k] += y[i]*cos(phase[i]);
    z[k] *= potential->dr_drho[k];
    if (i < 4) {
      if (h > 0.8) break;
    } else {
      if (h > 0.4) break;
    }
  }
  if (i > 1) {
    z[i0] = 0.0;
    if (x != NULL) z[i0] += x[0]*sin(phase[0]);
    if (y != NULL) z[i0] += y[0]*cos(phase[0]);
    z[i0] *= potential->dr_drho[i0];
    if (i == j) {
      i1 = i;
    } else {
      i1 = i + 1;
    }
    u1 = u + i1;
    w1 = w + i1;
    for (m = 0, n = i0; m < i; m++, n += 2) {
      u[m] = potential->rad[n];
      w[m] = potential->rad[n+1];
      z[n+1] = 0.0;
    }
    if (i1 > i) {
      u[i] = potential->rad[n];
    }
    UVIP3P(3, i1, u, phase, i, w, w1);
    if (x) {
      UVIP3P(3, i1, u, x, i, w, u1);
      for (m = 0, n = i0+1; m < i; m++, n += 2) {
        z[n] += u1[m]*sin(w1[m]);
      }
    }
    if (y) {
      UVIP3P(3, i1, u, y, i, w, u1);
      for (m = 0, n = i0+1; m < i; m++, n += 2) {
        z[n] += u1[m]*cos(w1[m]);
      }
    }
    for (m = 0, n = i0+1; m < i; m++, n += 2) {
      z[n] *= potential->dr_drho[n];
    }
    NewtonCotes(r, z, i0, k-2, t, 0);
  }

  q = i-1;
  m = j-q;
  if (m < 2) {
    return 0;
  }

  for (n = 1; n <= 5; n++) {
    i1 = q-n;
    if (i1 < 0 || phase[i1] >= phase[i1+1]) {
      i1++;
      break;
    }
  }
  if (t < 0) {
    for (n = i1, s = i0; n < j; n++, s += 2) {
      u[n] = potential->rad[s];
      z[n] = potential->rad[s+1];
    }
    UVIP3P(3, j-i1, u+i1, phase+i1, m, z+q, w+q);
  }

  if (x != NULL) {
    for (n = i1; n < j; n++) {
      x[n] /= dphase[n];
    }
    UVIP3C(3, j-i1, phase+i1, x+i1, z+i1, z+j+i1, x+j+i1);
  }
  if (y != NULL) {
    for (n = i1; n < j; n++) {
      y[n] /= dphase[n];
    }
    UVIP3C(3, j-i1, phase+i1, y+i1, u+i1, u+j+i1, y+j+i1);
  }

  si0 = sin(phase[i-1]);
  cs0 = cos(phase[i-1]);
  for (; i < j; i++, k += 2) {
    if (t < 0) {
      dr = w[i-1] - phase[i-1];
      si1 = sin(w[i-1]);
      cs1 = cos(w[i-1]);
      his0 = -(cs1 - cs0);
      hic0 = si1 - si0;
      p = dr;
      his1 = -p * cs1 + hic0;
      hic1 = p * si1 - his0;
      p *= dr; 
      his2 = -p * cs1 + 2.0*hic1;
      hic2 = p * si1 - 2.0*his1;
      p *= dr;
      his3 = -p * cs1 + 3.0*hic2;
      hic3 = p * si1 - 3.0*his2;
      r[k-1] = r[k-2];
    }
    d = phase[i] - phase[i-1];
    si1 = sin(phase[i]);
    cs1 = cos(phase[i]);
    is0 = -(cs1 - cs0);
    ic0 = si1 - si0;
    p = d;
    is1 = -p * cs1 + ic0;
    ic1 = p * si1 - is0;
    p *= d;
    is2 = -p * cs1 + 2.0*ic1;
    ic2 = p * si1 - 2.0*is1; 
    p *= d;
    is3 = -p * cs1 + 3.0*ic2;
    ic3 = p * si1 - 3.0*is2;
    r[k] = r[k-2];
    if (x != NULL) {
      a0 = x[i-1];
      a1 = z[i-1]; 
      a2 = z[j+i-1];
      a3 = x[j+i-1];
      h = a0*is0 + a1*is1 + a2*is2 + a3*is3;
      r[k] += h;
      if (t < 0) {
        h = a0*his0 + a1*his1 + a2*his2 + a3*his3;
        r[k-1] += h;
      }
    }
    if (y != NULL) {
      a0 = y[i-1];
      a1 = u[i-1];
      a2 = u[j+i-1];
      a3 = y[j+i-1];
      h = a0*ic0 + a1*ic1 + a2*ic2 + a3*ic3;
      r[k] += h;
      if (t < 0) {
        h = a0*hic0 + a1*hic1 + a2*hic2 + a3*hic3;
        r[k-1] += h;
      }
    }
    si0 = si1;
    cs0 = cs1;
  }

  return 0;
}

void LimitArrayRadial(int m, double n) {
  int k;

  k = (int)(n*1000000);
  switch (m) {
  case 0:
    yk_array->maxelem = k;
    break;
  case 1:
    slater_array->maxelem = k;
    break;
  case 2:
    breit_array->maxelem = k;
    break;
  case 3:
    gos_array->maxelem = k;
    break;
  case 4:
    moments_array->maxelem = k;
    break;
  case 5:
    multipole_array->maxelem = k;
    break;
  case 6:
    residual_array->maxelem = k;
    break;
  default:
    printf("m > 6, nothing is done\n");
    break;
  }
}

int InitRadial(void) {
  int ndim, i;
  int blocks[5] = {MULTI_BLOCK6,MULTI_BLOCK6,MULTI_BLOCK6,
		   MULTI_BLOCK6,MULTI_BLOCK6};

  potential = malloc(sizeof(POTENTIAL));
  potential->flag = 0;
  potential->uehling[0] = 0.0;
  SetBoundary(0, 1.0, -1.0);

  n_orbitals = 0;
  n_continua = 0;
  
  orbitals = malloc(sizeof(ARRAY));
  if (!orbitals) return -1;
  if (ArrayInit(orbitals, sizeof(ORBITAL), ORBITALS_BLOCK) < 0) return -1;

  ndim = 5;
  slater_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(slater_array, sizeof(double), ndim, blocks);

  ndim = 5;
  for (i = 0; i < ndim; i++) blocks[i] = MULTI_BLOCK5;
  breit_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(breit_array, sizeof(double), ndim, blocks);
  
  ndim = 2;
  for (i = 0; i < ndim; i++) blocks[i] = MULTI_BLOCK2;
  residual_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(residual_array, sizeof(double), ndim, blocks);

  vinti_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(vinti_array, sizeof(double), ndim, blocks);

  qed1e_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(qed1e_array, sizeof(double), ndim, blocks);
  
  ndim = 3;
  for (i = 0; i < ndim; i++) blocks[i] = MULTI_BLOCK4;
  multipole_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(multipole_array, sizeof(double *), ndim, blocks);

  ndim = 3;
  for (i = 0; i < ndim; i++) blocks[i] = MULTI_BLOCK3;
  moments_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(moments_array, sizeof(double), ndim, blocks);

  gos_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(gos_array, sizeof(double *), ndim, blocks);

  yk_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(yk_array, sizeof(SLATER_YK), ndim, blocks);

  n_awgrid = 1;
  awgrid[0]= EPS3;
  
  SetRadialGrid(DMAXRP, -1.0, -1.0, -1.0);
  SetSlaterCut(-1, -1);
  return 0;
}

int ReinitRadial(int m) {
  if (m < 0) return 0;
  SetSlaterCut(-1, -1);
  ClearOrbitalTable(m);
  FreeSimpleArray(slater_array);
  FreeSimpleArray(breit_array);
  FreeSimpleArray(residual_array);
  FreeSimpleArray(qed1e_array);
  FreeSimpleArray(vinti_array);
  FreeMultipoleArray();
  FreeMomentsArray();
  FreeYkArray();
  if (m < 2) {
    FreeGOSArray();
    if (m == 0) {
      if (optimize_control.n_screen > 0) {
	free(optimize_control.screened_n);
	optimize_control.n_screen = 0;
      }
      potential->flag = 0;
      n_awgrid = 1;
      awgrid[0] = EPS3;
      SetRadialGrid(DMAXRP, -1.0, -1.0, -1.0);
      potential->uehling[0] = 0.0;
    }
  }
  return 0;
}
 
int TestIntegrate(void) {
  ORBITAL *orb1, *orb2, *orb3, *orb4;
  int k1, k2, k3, k4, k, i, i0=1500;
  double r, a, s[6];
 
  k1 = OrbitalIndex(2, -2, 0);
  k2 = OrbitalIndex(3, -1, 0);
  k3 = OrbitalIndex(0, -5, 3.30265398e3/HARTREE_EV);
  k4 = OrbitalIndex(0, -4, 2.577e3/HARTREE_EV);
  printf("# %d %d %d %d\n", k1, k2, k3, k4);
  orb1 = GetOrbital(k1);
  orb2 = GetOrbital(k2);
  orb3 = GetOrbital(k3);
  orb4 = GetOrbital(k4);

  for (k = 1; k < 7; k++) {    
    for (i = 0; i < potential->maxrp; i++) {
      _yk[i] = pow(potential->rad[i],k);
      _yk[i+i0] = 1.0/_yk[i];
    }
    Integrate(_yk+i0, orb3, orb4, -1, _zk, 0);
    Integrate(_yk+i0, orb3, orb4, -1, _zk+i0, -1);
    for (i = 0; i < potential->maxrp; i++) {
      printf("%d %4d %15.8E %15.8E %15.8E %15.8E %15.8E\n", 
	     k, i, potential->rad[i], _yk[i], _zk[i], _yk[i+i0], _zk[i+i0]);
    }
    Integrate(_yk+i0, orb3, orb4, 1, &r, 0);
    Integrate(_yk+i0, orb3, orb4, 1, &a, -1);
    printf("# %15.8E %15.8E\n", r, a);
    printf("\n\n");
  }
  
  return 0;
}
