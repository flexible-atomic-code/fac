#include "radial.h"

static char *rcsid="$Id: radial.c,v 1.67 2002/10/23 03:44:10 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

static POTENTIAL *potential;

static ARRAY *orbitals;
static int n_orbitals;
static int n_continua;
 
static double _dwork[MAX_POINTS];
static double _dwork1[MAX_POINTS];
static double _dwork2[MAX_POINTS];
static double _dwork3[MAX_POINTS];
static double _dwork4[MAX_POINTS];
static double _dwork5[MAX_POINTS];
static double _dwork6[MAX_POINTS];
static double _dwork7[MAX_POINTS];
static double _dwork8[MAX_POINTS];
static double _dwork9[MAX_POINTS];
static double _dwork10[MAX_POINTS];
static double _dwork11[MAX_POINTS];
static double _phase[MAX_POINTS];
static double _dphase[MAX_POINTS];
static double _dphasep[MAX_POINTS];
static double _yk[MAX_POINTS];
static double _zk[MAX_POINTS];
static double _xk[MAX_POINTS];

static struct {
  double stablizer;
  double tolerance; /* tolerance for self-consistency */
  int maxiter; /* max iter. for self-consistency */
  double screened_charge; 
  int screened_kl;
  int n_screen;
  int *screened_n;
  int iprint; /* printing infomation in each iteration. */
  int iset;
} optimize_control = {0.5, EPS6, 128, 1.0, 1, 0, NULL, 0, 0};

#define NPB 5
static struct {
  double b[NPB][MAX_POINTS];
  double c[NPB];
  double u[MAX_POINTS];
} pbasis;

static struct {
  int se;
  int vp;
  int nms;
  int sms;
  int br;
} qed = {5, 2, 1, 1, 5};

static AVERAGE_CONFIG average_config = {0, 0, NULL, NULL, NULL};
 
static double rgrid_min;
static double rgrid_max;    
 
static MULTI *slater_array;
static MULTI *breit_array;
static MULTI *vinti_array;
static MULTI *qed1e_array;
static MULTI *residual_array;
static MULTI *multipole_array; 
static MULTI *moments_array;
static MULTI *gos_array;

static int n_awgrid = 0;
static double awgrid[MAXNTE];
static double aw2grid[MAXNTE];

double argam_(double *x, double *y);
double besljn_(int *jy, int *n, double *x);
void uvip3p_(int *np, int *ndp, double *x, double *y, 
	     int *n, double *xi, double *yi);
double _PhaseRDependent(double x, double eta, double b);
void lmqn_(int *, int *, double *, double *, double *, double *, int *,
	   void (*sfun)(int *, double *, double *, double *),
	   int *, int *, int *, double *, double *, double *, double *);

#ifdef PERFORM_STATISTICS
static RAD_TIMING rad_timing = {0, 0, 0, 0};
int GetRadTiming(RAD_TIMING *t) {
  memcpy(t, &rad_timing, sizeof(RAD_TIMING));
  return 0;
}
#endif

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
  
void SetOptimizeControl(double tolerance, double stablizer, 
			int maxiter, int iprint) {
  optimize_control.maxiter = maxiter;
  optimize_control.stablizer = stablizer;
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

int SetRadialGrid(double rmin, double rmax) {
  if (rmin > 0.0) rgrid_min = rmin;
  if (rmax > 0.0) rgrid_max = rmax;
  potential->flag = 0;
  return 0;
}

void _AdjustScreeningParams(double *u) {
  int i;
  double c;
  
  c = 0.5*u[MAX_POINTS-1];
  for (i = 0; i < MAX_POINTS; i++) {
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

  for (j = 0; j < MAX_POINTS; j++) {
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
  for (j = 0; j < MAX_POINTS; j++) {
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
	GetYk(t, _yk, orb1, orb1, -1);
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
	  GetYk(t, _yk, orb1, orb2, -1);
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
    for (j = jmax+1; j < MAX_POINTS; j++) {
      u[j] = u[jmax];
    }
    if (potential->N > 1) {
      for (j = jmax; j > 0; j--) {
        if (fabs(u[j]-potential->N + 1.0) > EPS10) break;
      }
      potential->r_core = j+1;
    }

    if (iter < 3) {
      r = 1.0;
      for (j = 0; j < MAX_POINTS; j++) {
	v[j] = u[j];
      }
    } else {	
      r = 0.0;
      k = 0;
      a = optimize_control.stablizer;
      b = 1.0 - a;
      for (j = 0; j < MAX_POINTS; j++) {
	if (u[j] + 1.0 != 1.0) {
	  r += fabs(1.0 - v[j]/u[j]);
	  k++;
	}
	u[j] = b*v[j] + a*u[j];
	v[j] = u[j];
      }
      r /= k;
    }
    _AdjustScreeningParams(u);
    SetPotentialVc(potential);
    for (j = 0; j < MAX_POINTS; j++) {
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
    k = 0;
    a = 0.0;
    r = potential->Z[MAX_POINTS-1];
    b = 1.0 - 1.0/potential->N;
    for (i = 0; i < acfg->n_shells; i++) {
      if (acfg->n[i] != k) {
	if (k > 0) {
	  r -= a*b;
	  c = r/k;
	  for (j = 0; j < MAX_POINTS; j++) {
	    u[j] += a*b*(exp(-c*potential->rad[j])-1.0);
	  }
	  a = 0.0;
	}
      } else {
	a += acfg->nq[i];
      }
    }
    _AdjustScreeningParams(u);
    SetPotentialVc(potential);
    for (j = 0; j < MAX_POINTS; j++) {
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
  
  fprintf(f, "Lambda = %10.3E, A = %10.3E\n",
	  potential->lambda, potential->a);

  
  for (j = 0; j < MAX_POINTS; j++) {
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
    GetYk(0, _yk, orb1, orb1, -1);
    for (k = 0; k < MAX_POINTS; k++) {
      v[k] += _yk[k]*acfg->nq[i]/potential->rad[k];
    }
    if (jmax < orb1->ilast) jmax = orb1->ilast;
    norbs++;
  }
  
  for (k = 0; k < MAX_POINTS; k++) {
    w[k] = w[k]/(potential->rad[k]*potential->rad[k]);
    w[k] = - pow(w[k], 1.0/3);
    ve1[k] += w[k]*0.4235655;
    ve0[k] = potential->Vc[k]+potential->U[k] - v[k];
  }
  ve1[0] = ve1[1];

  fprintf(f, "Mean configuration:\n");
  for (i = 0; i < acfg->n_shells; i++) {
    fprintf(f, "%-2d %2d\t%-10.3E\n", acfg->n[i], acfg->kappa[i], acfg->nq[i]);
  }
  fprintf(f, "\n\n");
  for (i = 0; i < MAX_POINTS; i++) {
    fprintf(f, "%-5d %11.5E %11.5E %11.5E %11.5E %11.5E %11.5E %11.5E\n",
	    i, potential->rad[i], potential->Z[i], 
	    potential->Vc[i]+potential->U[i], v[i], 
	    ve0[i], ve1[i], potential->uehling[i]);
  }

  fclose(f);  
  
  return 0;
}

double GetResidualZ(void) {
  double z;
  z = potential->Z[MAX_POINTS-1];
  if (potential->N > 0) z -= potential->N - 1;
  return z;
}

double GetRMax(void) {
  return potential->rad[MAX_POINTS-10];
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
  int i, j, k, m, no_old;
  double a, b, z;
  int iter;
  int *frozen;

  /* get the average configuration for the groups */
  acfg = &(average_config);
  if (ng > 0) {
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
    SetOrbitalRGrid(potential, rgrid_min, rgrid_max);
  }
  SetPotentialZ(potential, 0.0);
  z = potential->Z[MAX_POINTS-1];
  if (a > 0.0) z = z - a + 1;
  potential->a = 0.0;
  potential->lambda = 0.5*z;
  if (potential->N > 1) {
    potential->r_core = MAX_POINTS-5;
  } else {
    potential->r_core = MAX_POINTS*0.75;
  }

  if (optimize_control.iset == 0) {
    optimize_control.stablizer = 0.25 + 0.75*(z/potential->Z[MAX_POINTS-1]);
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

static void TNFunc(int *n, double *x, double *f, double *g) {
  int i, j;
  int m;
  double *u, a, delta;

  u = potential->U;
  for (j = 0; j < MAX_POINTS; j++) {
    u[j] = pbasis.u[j];
  }
  for (i = 0; i < *n; i++) {
    if (x[i]) {
      for (j = 0; j < MAX_POINTS; j++) {
	u[j] += x[i]*pbasis.b[i][j];
      }
    }
  }

  SetPotentialU(potential, 0, NULL);
  ReinitRadial(1);
  ClearOrbitalTable(0);
  *f = AverageEnergyAvgConfig(&average_config);

  for (i = 0; i < *n; i++) {
    delta = 0.01*x[i];
    if (delta < EPS3) delta = EPS3;
    for (j = 0; j < MAX_POINTS; j++) {
      u[j] += delta*pbasis.b[i][j];
    }
    SetPotentialU(potential, 0, NULL);
    ReinitRadial(1);
    ClearOrbitalTable(0);
    a = AverageEnergyAvgConfig(&average_config);
    g[i] = (a - *f)/delta;
    for (j = 0; j < MAX_POINTS; j++) {
      u[j] -= delta*pbasis.b[i][j];
    }
  }

  return;
}

int RefineRadial(int maxfun, int msglvl) {
  int i, n, lw, ierr;
  int maxit;
  double eta, stepmx, accrcy, xtol;
  ORBITAL orb;
  double f0, f, g[NPB];

  if (maxfun <= 0) return 0;
  memcpy(pbasis.u, potential->U, sizeof(double)*MAX_POINTS);
  for (i = 0; i < NPB; i++) {
    pbasis.c[i] = 0.0;
    orb.n = i + 1;
    orb.kappa = -1;
    potential->flag = 1;
    ierr = RadialSolver(&orb, potential, EPS8);
    if (ierr) return ierr;
    memcpy(pbasis.b[i], orb.wfun, sizeof(double)*MAX_POINTS);    
    free(orb.wfun);
  }
  
  if (msglvl == 0) msglvl = -3;
  maxit = 60;
  eta = 0.25;
  stepmx = 0.5;
  accrcy = EPS8;
  xtol = EPS2;
  n = NPB;
  lw = MAX_POINTS;
  
  TNFunc(&n, pbasis.c, &f, g);
  f0 = f;
  lmqn_(&ierr, &n, pbasis.c, &f, g, _dwork11, &lw, TNFunc, 
	&msglvl, &maxit, &maxfun, &eta, &stepmx, &accrcy, &xtol);

  if (ierr) {
    if (f > f0) return ierr;
  }
  
  return 0;
}

int SolveDirac(ORBITAL *orb) {
  double eps;
  int err;
#ifdef PERFORM_STATISTICS
  clock_t start, stop;
  start = clock();
#endif
  
  err = 0;
  eps = 1E-8;
  
  potential->flag = -1;
  err = RadialSolver(orb, potential, eps);
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
  
  orb = GetOrbital(k);
  
  fprintf(f, "#Wave Function for n = %d, kappa = %d, energy = %12.6E\n\n",
	  n, kappa, orb->energy*HARTREE_EV);
  if (n > 0) {
    for (i = 0; i <= orb->ilast; i++) {
      fprintf(f, "%-4d %10.3E %10.3E %10.3E %10.3E %10.3E\n", 
	      i, potential->rad[i], 
	      (potential->Vc[i])*potential->rad[i],
	      potential->U[i] * potential->rad[i],
	      Large(orb)[i], Small(orb)[i]); 
    }
  } else {
    z = GetResidualZ();
    e = orb->energy;
    a = FINE_STRUCTURE_CONST2 * e;
    ke = sqrt(2.0*e*(1.0+0.5*a));
    y = (1.0+a)*z/ke; 
    for (i = 0; i <= orb->ilast; i++) {
      fprintf(f, "%-4d %10.3E %10.3E %10.3E %10.3E %10.3E\n", 
	      i, potential->rad[i],
	      (potential->Vc[i])*potential->rad[i],
	      potential->U[i] * potential->rad[i],
	      Large(orb)[i], Small(orb)[i]); 
    }
    for (; i < MAX_POINTS; i += 2) {
      a = ke * potential->rad[i];
      a = a + y*log(2.0*a);
      fprintf(f, "%-4d %10.3E %13.6E %13.6E %13.6E %13.6E\n",
	      i, potential->rad[i],
	      Large(orb)[i], Large(orb)[i+1], 
	      Small(orb)[i], a);
    }
  }

  fclose(f);

  return 0;
}

double _PhaseRDependent(double x, double eta, double b) {
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

  i = MAX_POINTS - 1;  
  phase1 = orb->wfun[i];
  r = potential->rad[i-1];  
  b1 = orb->kappa;
  b1 = b1*(b1+1.0) - FINE_STRUCTURE_CONST2*z*z;
 
  a = ke * r;
  b1 = _PhaseRDependent(a, y, b1);
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
      if (orb->kappa == kappa && 
	  orb->energy > 0.0 &&
	  fabs(orb->energy - energy) < EPS6) {
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
    orb->n = -n_continua;
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
	  fabs(orb->energy - energy) < EPS6) 
	return i;
    } else if (orb->n == n && orb->kappa == kappa) {
      return i;
    }
  }
  return -1;
}

int AddOrbital(ORBITAL *orb) {

  if (orb == NULL) return -1;

  orb = (ORBITAL *) ArrayAppend(orbitals, orb);
  if (!orb) {
    printf("Not enough memory for orbitals array\n");
    exit(1);
  }

  if (orb->n == 0) {
    n_continua++;
    orb->n = -n_continua;
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

  orb = (ORBITAL *) ArrayAppend(orbitals, NULL);
  if (!orb) {
    printf("Not enough memory for orbitals array\n");
    exit(1);
  }

  orb->wfun = NULL;
  orb->ilast = -1;

  n_orbitals++;
  return orb;
}

void _FreeOrbitalData(void *p) {
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
    ArrayFree(orbitals, _FreeOrbitalData);
  } else {
    for (i = n_orbitals-1; i >= 0; i--) {
      orb = GetOrbital(i);
      if (orb->n > 0) {
	n_continua -= n_orbitals - (i+1);
	n_orbitals = i+1;
	ArrayTrim(orbitals, i+1, _FreeOrbitalData);
	break;
      }
    }
    if (m == 2) {
      for (; i >= 0; i--) {
	FreeOrbital(i);
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
  _FreeOrbitalData((void *)orb);
  return 0;
}

int SaveAllContinua(int mode) {
  int i;
  ORBITAL *orb;
  for (i = 0; i < n_orbitals; i++) {
    orb = GetOrbital(i);
    if (orb->n <= 0 && orb->wfun != NULL) {
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
    if (orb->n <= 0 && 
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
    if (orb->n <= 0 && orb->wfun != NULL) {
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
    if (orb->n <= 0 && 
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

    } else b = 0.0;
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
  double *p, z;

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
  
  p = (double *) MultiSet(residual_array, index, NULL);
  if (p && *p) {
    *s = *p;
    return 0;
  } 

  *s = 0.0;
 
  for (i = 0; i < MAX_POINTS; i++) {
    z = potential->U[i];
    z += potential->Vc[i];
    _yk[i] = -(potential->Z[i]/potential->rad[i]) - z;
  }
  Integrate(_yk, orb1, orb2, 1, s);
  *p = *s;
  return 0;
}

double RadialMoments(int m, int k1, int k2) {
  int index[3], npts, i;
  ORBITAL *orb1, *orb2;
  double *q, r, z;
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

  if (n1 > 0 && n2 > 0) {
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

  if (n1 == n2 && m > 0 && n1 > GetNMax()) {
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
  
  q = (double *) MultiSet(moments_array, index, NULL);
 
  if (*q) {
    return *q;
  } 

  npts = MAX_POINTS-1;
  if (orb1->n > 0) npts = Min(npts, orb1->ilast);
  if (orb2->n > 0) npts = Min(npts, orb2->ilast);

  r = 0.0;
  for (i = 0; i <= npts; i++) {
    _yk[i] = pow(potential->rad[i], m);
  }
  Integrate(_yk, orb1, orb2, 1, &r);
  *q = r;

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
    /**********************************************************/ 
    /* M1 needs special treatments, because the lowest order  */
    /* non-relativistic approximation is simply the overlap   */
    /* integral, which vanishes in most cases.                */
    /**********************************************************/
    if (k1 == k2) {
      t = kappa1 + kappa2;
      p = m - t;
      if (p && t) {
	r = -0.5*FINE_STRUCTURE_CONST * (p*t) / sqrt(m*(m+1));
	r *= ReducedCL(GetJFromKappa(kappa1), 2*m, GetJFromKappa(kappa2));
      }
    } else {
      /* the M1 radial integral vanish in this approximation. 
	 instead of going to higher orders, use the relativistic 
	 version is just simpler */
      printf("should call MultipoleRadialFR instead\n");
    }
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
  ORBITAL *orb1, *orb2;
  double x, a, r, rp, **p1, **p2, aw2, ef;
  int jy, n, i, j, npts;

#ifdef PERFORM_STATISTICS 
  clock_t start, stop;
  start = clock();
#endif

  if (m == 0) return 0.0;
  
  index[0] = 0;
  if (m >= 0) {
    index[1] = 2*m;
  } else {
    index[1] = -2*m-1;
  }
 
  if (k1 <= k2) {
    s = 0;      
    index[2] = k1;
    index[3] = k2;
    orb1 = GetOrbitalSolved(k1);
    orb2 = GetOrbitalSolved(k2);
  } else {
    s = 1;     
    index[2] = k2;
    index[3] = k1;
    orb1 = GetOrbitalSolved(k2);
    orb2 = GetOrbitalSolved(k1);
  }

  if (orb1->wfun == NULL || orb2->wfun == NULL) {
    if (m == -1) {
      return MultipoleRadialNR(m, k1, k2, gauge);
    } else {
      return 0.0;
    }
  }
  
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

  p1 = (double **) MultiSet(multipole_array, index, NULL);
  p2 = NULL;
  if (m < 0 && gauge == G_BABUSHKIN) {
    index[0] = 1;
    p2 = (double **) MultiSet(multipole_array, index, NULL);
  }

  if (*p1) {
    r = InterpolateMultipole(aw2, n_awgrid, aw2grid, *p1);
    if (p2) {
      rp = InterpolateMultipole(aw2, n_awgrid, aw2grid, *p2);
      if (s == 1) rp = -rp;
      r += rp;
    } 
    if (s == 1 && gauge == G_COULOMB) r = -r;
    return r;
  }
  
  *p1 = (double *) malloc(sizeof(double)*n_awgrid);
  if (p2) {
    *p2 = (double *) malloc(sizeof(double)*n_awgrid);
  }
  
  kappa1 = orb1->kappa;
  kappa2 = orb2->kappa;
  npts = MAX_POINTS-1;
  if (orb1->n > 0) npts = Min(npts, orb1->ilast);
  if (orb2->n > 0) npts = Min(npts, orb2->ilast);
  r = 0.0;
  jy = 1;

  for (i = 0; i < n_awgrid; i++) {
    a = awgrid[i];
    if (ef > 0.0) a += ef;
    if (m > 0) {
      t = kappa1 + kappa2;
      if (t) {
	for (j = 0; j < npts; j++) {
	  x = a*potential->rad[j];
	  n = m;
	  _yk[j] = besljn_(&jy, &n, &x);
	}
	Integrate(_yk, orb1, orb2, 4, &r);
	r *= t;
	r *= (2*m + 1.0)/sqrt(m*(m+1.0));
	r /= pow(a, m);
	r *= ReducedCL(GetJFromKappa(kappa1), 2*m, GetJFromKappa(kappa2));
	(*p1)[i] = r;
      }
    } else {
      am = -m;
      if (gauge == G_COULOMB) {
	t = kappa1 - kappa2;
	q = sqrt(am/(am+1.0));
	for (j = 0; j < npts; j++) {
	  x = a*potential->rad[j];
	  n = am+1;
	  _yk[j] = besljn_(&jy, &n, &x);
	  n = am-1;
	  _zk[j] = besljn_(&jy, &n, &x);
	}
	if (t) {
	  Integrate(_yk, orb1, orb2, 4, &ip);
	  Integrate(_zk, orb1, orb2, 4, &ipm);
	  r = t*ip*q - t*ipm/q;
	}
	if (k1 != k2) {
	  Integrate(_yk, orb1, orb2, 5, &im);
	  Integrate(_zk, orb1, orb2, 5, &imm);
	  r += (am + 1.0)*im*q + am*imm/q;
	}
	r /= pow(a,am-1);
	q = ReducedCL(GetJFromKappa(kappa1), 2*am, GetJFromKappa(kappa2));
	r *= q;
	(*p1)[i] = r;
      } else if (gauge == G_BABUSHKIN) {
	t = kappa1 - kappa2;
	for (j = 0; j < npts; j++) {
	  x = a*potential->rad[j];
	  n = am+1;
	  _yk[j] = besljn_(&jy, &n, &x);
	  n = am;
	  _zk[j] = besljn_(&jy, &n, &x);
	}
	if (t) {
	  Integrate(_yk, orb1, orb2, 4, &ip);
	  r = t*ip;
	}
	if (k1 != k2) {
	  Integrate(_yk, orb1, orb2, 5, &im);
	} else {
	  im = 0.0;
	}
	Integrate(_zk, orb1, orb2, 1, &imm);
	rp = (am + 1.0) * (imm + im);
	q = (2*am + 1.0)/sqrt(am*(am+1.0));
	q /= pow(a, am);
	r *= q;
	rp *= q;
	q = ReducedCL(GetJFromKappa(kappa1), 2*am, GetJFromKappa(kappa2));
	r *= q;
	rp *= q;
	(*p1)[i] = r;
	(*p2)[i] = rp;
      }
    }
  }

  /*
  for (i = 0; i < n_awgrid; i++) {
    printf("%d %d %10.3E %10.3E %10.3E ", 
	    kappa1, kappa2, awgrid[i], aw2grid[i], (*p1)[i]);
    if (p2) printf("%10.3E ", (*p2)[i]);
    printf("\n");
  }
  printf("\n\n");
  */

  r = InterpolateMultipole(aw2, n_awgrid, aw2grid, *p1);
  if (p2) {
    rp = InterpolateMultipole(aw2, n_awgrid, aw2grid, *p2);
    if (s == 1) rp = -rp;
    r += rp;
  }

  if (s == 1 && gauge == G_COULOMB) r = -r;
  
#ifdef PERFORM_STATISTICS 
  stop = clock();
  rad_timing.radial_1e += stop - start;
#endif
  return r;
}

double *GeneralizedMoments(int nk, double *kg, int k1, int k2, int m) {
  ORBITAL *orb1, *orb2;
  int n1, i, jy;
  double x, r, r0;
  double *p1, *p2, *q1, *q2;
  int index[3], t;
  double **p, k;
  
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

  p = (double **) MultiSet(gos_array, index, NULL);
  if (*p) {
    return *p;
  }

  *p = (double *) malloc(sizeof(double)*nk);

  if (orb1->wfun == NULL || orb2->wfun == NULL) {
    for (t = 0; t < nk; t++) {
      (*p)[t] = 0.0;
    }
    return *p;
  }

  jy = 1;
  p1 = Large(orb1);
  p2 = Large(orb2);
  q1 = Small(orb1);
  q2 = Small(orb2);
  
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
    k = kg[t];
    for (i = 0; i <= n1; i++) {
      x = k * potential->rad[i];
      _dphase[i] = besljn_(&jy, &m, &x);
      _dphase[i] *= _phase[i];
    }
    r = Simpson(_dphase, 0, n1);

    (*p)[t] = r - r0;
  }
  return *p;
}
    
double InterpolateMultipole(double aw2, int n, double *x, double *y) {
  double r;
  int np, nd;

  if (n == 1) {
    r = y[0];
  } else {
    np = 3;
    nd = 1;
    uvip3p_(&np, &n, x, y, &nd, &aw2, &r);
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
	d -= Vinti(k0, k2) * Vinti(k1, k3) / am;
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
    if (fabs(a) > EPS10) {
      e = 0.0;
      if (IsEven((kl0+kl3+t)/2) && IsEven((kl1+kl2+t)/2)) {
	err = Slater(&e, k0, k1, k3, k2, t/2, mode);
      }
      if (qed.br < 0 || (maxn > 0 && maxn <= qed.br)) {
	e += Breit(k0, k1, k3, k2, t/2, kl0, kl1, kl3, kl2);
      }
      if (t == 2 && qed.sms && maxn > 0) {
	e -= Vinti(k0, k3) * Vinti(k1, k2) / am;
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
  double r0 = 1E-3;
  double *p, *q, e, z;
  double *large, *small;
  double a, b;
  
  if (orb->wfun == NULL) return 1.0;

  for (i = 0; i < MAX_POINTS; i++) {
    if (potential->rad[i] > r0) break;
  }
  npts = i;
  
  p = _xk;
  q = _zk;
  z = potential->Z[MAX_POINTS-1];
  e = RadialDiracCoulomb(npts, p, q, potential->rad, z, 
			 orb->n, orb->kappa);
  large = Large(orb);
  small = Small(orb);  
  for (i = 0; i < npts; i++) {
    p[i] = (p[i]*p[i] + q[i]*q[i])*potential->dr_drho[i];
    q[i] = (large[i]*large[i] + small[i]*small[i])*potential->dr_drho[i];
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
  
  p = (double *) MultiSet(qed1e_array, index, NULL);
  if (p && *p) {
    return *p;
  }

  r = 0.0;
  
  if (qed.nms > 0) {
    for (i = 0; i < MAX_POINTS; i++) {
      _yk[i] = potential->U[i] + potential->Vc[i];
    }
    a = 0.0;
    Integrate(_yk, orb1, orb2, 1, &a);
    a = -a;
    if (k0 == k1) a += orb1->energy;
    a /= (AMU * GetAtomicMass());
    r += a;
  }

  if (qed.vp > 0) {
    a = 0.0;
    Integrate(potential->uehling, orb1, orb2, 1, &a);
    r += a;
  }

  if (k0 == k1 && orb1->n <= qed.se) {
    a = HydrogenicSelfEnergy(potential->Z[MAX_POINTS-1], 
			     orb1->n, orb1->kappa);
    if (a) {
      a *= SelfEnergyRatio(orb1);
      r += a;
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
  
  if (k0 > k1) {
    index[0] = k1;
    index[1] = k0;
  } else {
    index[0] = k0;
    index[1] = k1;
  }
  
  p = (double *) MultiSet(vinti_array, index, NULL);
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

  Differential(large1, _zk, 0, MAX_POINTS-1);
  for (i = 0; i < MAX_POINTS; i++) {
    _yk[i] = large0[i]*_zk[i] - a*large0[i]*large1[i]/potential->rad[i];
    _yk[i] *= potential->dr_drho[i];
  }
  r += Simpson(_yk, 0, MAX_POINTS-1);
  
  Differential(small1, _zk, 0, MAX_POINTS-1);
  for (i = 0; i < MAX_POINTS; i++) {
    _yk[i] = small0[i]*_zk[i] - b*small0[i]*small1[i]/potential->rad[i];
    _yk[i] *= potential->dr_drho[i];
  }
  r += Simpson(_yk, 0, MAX_POINTS-1);
  
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

  p = (double *) MultiSet(breit_array, index, NULL);
  if (p && *p) {
    r = *p;
  } else {
    orb0 = GetOrbitalSolved(k0);
    orb1 = GetOrbitalSolved(k1);
    orb2 = GetOrbitalSolved(k2);
    orb3 = GetOrbitalSolved(k3);
    if (!orb0 || !orb1 || !orb2 || !orb3) return 0.0;
    
    for (i = 0; i < MAX_POINTS; i++) {
      _dwork1[i] = pow(potential->rad[i], k);
    }
    
    Integrate(_dwork1, orb0, orb1, -6, _zk);
    
    for (i = 0; i < MAX_POINTS; i++) {
      _zk[i] /= _dwork1[i]*potential->rad[i];
    }

    Integrate(_zk, orb2, orb3, 6, &r);
    
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
  int index[6];
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
  
  switch (mode) {
  case 0:
  case 1:
    index[5] = 0;
    break;
  case -1:
    index[5] = 1;
    break;
  case 2:
    index[5] = 2;
    break;
  case -2:
    index[5] = 3;
    break;
  default:
    printf("mode unrecognized in slater\n");
    return -1;
  }

  if (index[5] < 2) {
    SortSlaterKey(index);
    p = (double *) MultiSet(slater_array, index, NULL);
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

    npts = MAX_POINTS;
    switch (mode) {
    case 0: /* fall through to case 1 */
    case 1: /* full relativistic with distorted free orbitals */
      GetYk(k, _yk, orb0, orb2, -1); 
      if (orb1->n > 0) ilast = orb1->ilast;
      else ilast = npts-1;
      if (orb3->n > 0) ilast = Min(ilast, orb3->ilast);
      for (i = 0; i <= ilast; i++) {
	_yk[i] = (_yk[i]/potential->rad[i]);
      }
      Integrate(_yk, orb1, orb3, 1, s);
      break;
    
    case -1: /* quasi relativistic with distorted free orbitals */
      GetYk(k, _yk, orb0, orb2, -2); 
      if (orb1->n > 0) ilast = orb1->ilast;
      else ilast = npts-1;
      if (orb3->n > 0) ilast = Min(ilast, orb3->ilast);
      for (i = 0; i <= ilast; i++) {
	_yk[i] /= potential->rad[i];
      }
      Integrate(_yk, orb1, orb3, 2, s);

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

int GetYk(int k, double *yk, ORBITAL *orb1, ORBITAL *orb2, int type) {
  int i, ilast;
  int i0;
  double a, max;

  for (i = 0; i < MAX_POINTS; i++) {
    _dwork1[i] = pow(potential->rad[i], k);
  }

  Integrate(_dwork1, orb1, orb2, type, _zk);
  
  for (i = 0; i < MAX_POINTS; i++) {
    _zk[i] /= _dwork1[i];
    yk[i] = _zk[i];
  }

  if (k > 2) {
    max = 0.0;
    for (i = 0; i < MAX_POINTS; i++) {
      a = fabs(yk[i]);
      if (max < a) max = a;
    }
    max *= 1E-3;
    for (i = 0; i < Min(orb1->ilast, orb2->ilast); i++) {
      a = Large(orb1)[i]*Large(orb2)[i]*potential->rad[i];
      if (fabs(a) > max) break;
      _dwork1[i] = 0.0;
    }
    i0 = i;
  } else i0 = 0;
  for (i = i0; i < MAX_POINTS; i++) {
    _dwork1[i] = pow(potential->rad[i0]/potential->rad[i], k+1);
  }
  Integrate(_dwork1, orb1, orb2, type, _xk);
  ilast = MAX_POINTS - 1;
  
  for (i = i0; i < MAX_POINTS; i++) {
    _xk[i] = (_xk[ilast] - _xk[i])/_dwork1[i];
    yk[i] += _xk[i];
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

int Integrate(double *f, ORBITAL *orb1, ORBITAL *orb2, 
	      int t, double *x) {
  int i1, i2;
  int i, type;
  double *r;

  if (t == 0) t = 1;
  if (t < 0) {
    r = x;
    type = -t;
  } else {
    r = _dwork;
    type = t;
  }

  r[0] = 0.0;

  if (orb1->n > 0 && orb2->n > 0) {
    i2 = Min(orb1->ilast, orb2->ilast);
    IntegrateSubRegion(0, i2, f, orb1, orb2, t, r, 0);
  } else if (orb1->n > 0 && orb2->n <= 0) {
    i1 = Min(orb1->ilast, orb2->ilast);
    IntegrateSubRegion(0, i1, f, orb1, orb2, t, r, 0);
    i1 += 1;
    i2 = orb1->ilast;
    IntegrateSubRegion(i1, i2, f, orb1, orb2, t, r, 1);
  } else if (orb1->n <= 0 && orb2->n > 0) {
    i1 = Min(orb1->ilast, orb2->ilast);
    IntegrateSubRegion(0, i1, f, orb1, orb2, t, r, 0);
    i1 += 1;
    i2 = orb2->ilast;
    if (type == 6) {
      i = 7;
      if (t < 0) i = -i;
      IntegrateSubRegion(i1, i2, f, orb2, orb1, i, r, 1);
    } else {
      IntegrateSubRegion(i1, i2, f, orb1, orb2, t, r, 2);
    }
  } else {
    i1 = Min(orb1->ilast, orb2->ilast);
    IntegrateSubRegion(0, i1, f, orb1, orb2, t, r, 0);
    i1 += 1;
    if (i1 > orb1->ilast) {
      i2 = orb2->ilast;
      if (type == 6) {
	i = 7;
	if (t < 0) i = -i;
	IntegrateSubRegion(i1, i2, f, orb2, orb1, i, r, 1);
      } else {
	IntegrateSubRegion(i1, i2, f, orb1, orb2, t, r, 2);
      }
      i1 = i2+1;
      i2 = MAX_POINTS-1;
      IntegrateSubRegion(i1, i2, f, orb1, orb2, t, r, 3);
    } else {
      i2 = orb1->ilast;
      IntegrateSubRegion(i1, i2, f, orb1, orb2, t, r, 1);
      i1 = i2+1;
      i2 = MAX_POINTS-1;
      IntegrateSubRegion(i1, i2, f, orb1, orb2, t, r, 3);
    }
  }

  if (t >= 0) {
    *x = r[i2-1];
  } else {
    i2++;
    for (i = i2+1; i < MAX_POINTS; i++) {
      r[i] = r[i2];
    }
  }

  return 0;
}

int IntegrateSubRegion(int i0, int i1, 
		       double *f, ORBITAL *orb1, ORBITAL *orb2,
		       int t, double *r, int m) {
  int i, j, ip, type;
  ORBITAL *tmp;
  double *large1, *large2, *small1, *small2;
  double *x, *y, *r1, *x1, *x2, *y1, *y2;
  double a, b, e1, e2, a2;

  if (i1 <= i0) return 0;
  type = abs(t);

  x = _dwork3;
  y = _dwork4;
  x1 = _dwork5;
  x2 = _dwork6;
  y1 = _dwork7;
  y2 = _dwork8;
  r1 = _dwork9;

  switch (m) {
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
      break;
    case 2: /* type = 2 */
      for (i = i0; i <= i1; i++) {
	x[i] = large1[i] * large2[i];
	x[i] *= f[i]*potential->dr_drho[i];
      }
      break;
    case 3: /*type = 3 */
      for (i = i0; i <= i1; i++) {
	x[i] = small1[i] * small2[i];
	x[i] *= f[i]*potential->dr_drho[i];
      }
      break;
    case 4: /*type = 4 */
      for (i = i0; i <= i1; i++) {
	x[i] = large1[i] * small2[i];
	x[i] += small1[i] * large2[i];
	x[i] *= f[i]*potential->dr_drho[i];
      }
      break;
    case 5: /* type = 5 */
      for (i = i0; i <= i1; i++) {
	x[i] = large1[i] * small2[i];
	x[i] -= small1[i] * large2[i];
	x[i] *= f[i]*potential->dr_drho[i];
      } 
      break;
    case 6: /* type = 6 */
      for (i = i0; i <= i1; i++) {
	x[i] = large1[i] * small2[i];
	x[i] *= f[i]*potential->dr_drho[i];
      }
      break;
    default: /* error */
      return -1;
    }
    NewtonCotes(r, x, i0, i1, t);
    break;

  case 1:
    if (type == 6) { /* type 6 needs special treatments */
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
      if (i < MAX_POINTS) {
	ip = i+1;
	if (i > orb1->ilast) { 
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  b = small1[i]*b + small1[ip]*a;
	  a = large1[i]*a;
	} else {
	  a = large1[i];
	  b = small1[i];
	}
	x[j] = a * small2[ip];
	x[j] *= f[i];
	y[j] = a * small2[i] * f[i];
	_phase[j] = large2[ip];
	_dphase[j] = (1.0 + a2*(e2-potential->U[i]-potential->Vc[i]))
	  /(large2[i]*large2[i]);
	j++;
      }
      IntegrateSinCos(j, x, y, _phase, _dphase, i0, r, t);
      if (t < 0) {
	for (i = i0+1; i < i1; i += 2) {
	  r[i] = 0.5*(r[i-1] + r[i+1]);
	}
	if (i1 < MAX_POINTS-1) r[i1] = 0.5*(r[i1-1] + r[i1+1]);
	else r[i1] = r[i1-1];
      }
      break;
    } else if (type == 7) {
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
      if (i < MAX_POINTS) {
	ip = i+1;
	if (i > orb1->ilast) { 
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  b = small1[i]*b + small1[ip]*a;
	} else {
	  b = small1[i];
	}
	x[j] = b * large2[i];
	x[j] *= f[i];
	_phase[j] = large2[ip];
	_dphase[j] = (1.0 + a2*(e2-potential->U[i]-potential->Vc[i]))
	  /(large2[i]*large2[i]);
	j++;
      }
      IntegrateSinCos(j, x, NULL, _phase, _dphase, i0, r, t);
      if (t < 0) {
	for (i = i0+1; i < i1; i += 2) {
	  r[i] = 0.5*(r[i-1] + r[i+1]);
	}
	if (i1 < MAX_POINTS-1) r[i1] = 0.5*(r[i1-1] + r[i1+1]);
	else r[i1] = r[i1-1];
      }
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
      if (i < MAX_POINTS) {
	ip = i+1;
	if (i > orb1->ilast) {
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  b = small1[i]*b + small1[ip]*a;
	  a = large1[i]*a;
	} else {
	  a = large1[i];
	  b = small1[i];
	}
	x[j] = a * large2[i];
	x[j] += b * small2[ip];
	x[j] *= f[i];
	y[j] = b * small2[i] * f[i];
	_phase[j] = large2[ip];
	_dphase[j] = (1.0 + a2*(e2-potential->U[i]-potential->Vc[i]))
	  /(large2[i]*large2[i]);
	j++;
      }
      IntegrateSinCos(j, x, y, _phase, _dphase, i0, r, t);
      if (t < 0) {
	for (i = i0+1; i < i1; i += 2) {
	  r[i] = 0.5*(r[i-1] + r[i+1]);
	}
	if (i1 < MAX_POINTS-1) r[i1] = 0.5*(r[i1-1] + r[i1+1]);
	else r[i1] = r[i1-1];
      }
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
      if (i < MAX_POINTS) {
	ip = i+1;
	if (i > orb1->ilast) {
	  a = sin(large1[ip]);
	  a = large1[i]*a;
	} else {
	  a = large1[i];
	}
	x[j] = a * large2[i];
	x[j] *= f[i];
	_phase[j] = large2[ip];
	_dphase[j] = (1.0 + a2*(e2-potential->U[i]-potential->Vc[i]))
	  /(large2[i]*large2[i]);
	j++;
      }
      IntegrateSinCos(j, x, NULL, _phase, _dphase, i0, r, t);

      if (t < 0) {
	for (i = i0+1; i < i1; i += 2) {
	  r[i] = 0.5*(r[i-1] + r[i+1]);
	}
	if (i1 < MAX_POINTS-1) r[i1] = 0.5*(r[i1-1] + r[i1+1]);
	else r[i1] = r[i1-1];
      }
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
      if (i < MAX_POINTS) {
	ip = i+1;
	if (i > orb1->ilast) {
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  b = small1[i]*b + small1[ip]*a;
	} else {
	  b = small1[i];
	}
	x[j] += b * small2[ip];
	x[j] *= f[i];
	y[j] = b * small2[i] * f[i];
	_phase[j] = large2[ip];
	_dphase[j] = (1.0 + a2*(e2-potential->U[i]-potential->Vc[i]))
	  /(large2[i]*large2[i]);
	j++;
      }
      IntegrateSinCos(j, x, y, _phase, _dphase, i0, r, t);
      if (t < 0) {
	for (i = i0+1; i < i1; i += 2) {
	  r[i] = 0.5*(r[i-1] + r[i+1]);
	}
	if (i1 < MAX_POINTS-1) r[i1] = 0.5*(r[i1-1] + r[i1+1]);
	else r[i1] = r[i1-1];
      }
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
      if (i < MAX_POINTS) {
	ip = i+1;
	if (i > orb1->ilast) {
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  b = small1[i]*b + small1[ip]*a;
	  a = large1[i]*a;
	} else {
	  a = large1[i];
	  b = small1[i];
	}
	x[j] = a * small2[ip];
	x[j] += b * large2[i];
	x[j] *= f[i];
	y[j] = a * small2[i] * f[i];
	_phase[j] = large2[ip];
	_dphase[j] = (1.0 + a2*(e2-potential->U[i]-potential->Vc[i]))
	  /(large2[i]*large2[i]);
	j++;
      }
      IntegrateSinCos(j, x, y, _phase, _dphase, i0, r, t);
      if (t < 0) {
	for (i = i0+1; i < i1; i += 2) {
	  r[i] = 0.5*(r[i-1] + r[i+1]);
	}
	if (i1 < MAX_POINTS-1) r[i1] = 0.5*(r[i1-1] + r[i1+1]);
	else r[i1] = r[i1-1];
      }
      break;
    case 5: /* type = 5 */
      j = 0;
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
      if (i < MAX_POINTS) {
	ip = i+1;
	if (i > orb1->ilast) { 
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  b = small1[i]*b + small1[ip]*a;
	  a = large1[i]*a;
	} else {
	  a = large1[i];
	  b = small1[i];
	}
	x[j] = a * small2[ip];
	x[j] -= b * large2[i];
	x[j] *= f[i];
	y[j] = a * small2[i] * f[i];
	_phase[j] = large2[ip];
	_dphase[j] = (1.0 + a2*(e2-potential->U[i]-potential->Vc[i]))
	  /(large2[i]*large2[i]);
	j++;
      }
      IntegrateSinCos(j, x, y, _phase, _dphase, i0, r, t);
      if (m == 2) {
	if (t < 0) {
	  for (i = i0; i <= i1; i += 2) {
	    r[i] = -r[i];
	  }
	  if (i < MAX_POINTS) r[i] = -r[i];
	} else {
	  i = i1 - 1;
	  r[i] = -r[i];
	  i = i1 + 1;
	  if (i < MAX_POINTS) r[i] = -r[i];
	}
      }
      if (t < 0) {
	for (i = i0+1; i < i1; i += 2) {
	  r[i] = 0.5*(r[i-1] + r[i+1]);
	}
	if (i1 < MAX_POINTS-1) r[i1] = 0.5*(r[i1-1] + r[i1+1]);
	else r[i1] = r[i1-1];
      }
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
      IntegrateSinCos(j, x, y, _phase, _dphase, i0, r, t);
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
      IntegrateSinCos(j, x, y, _phase, _dphasep, i0, r1, t);
      if (t < 0) {
	for (i = i0; i <= i1; i += 2) {
	  r[i] += r1[i];
	}
	if (i < MAX_POINTS) r[i] += r1[i];
	for (i = i0+1; i < i1; i += 2) {
	  r[i] = 0.5*(r[i-1] + r[i+1]);
	}
	if (i1 < MAX_POINTS-1) r[i1] = 0.5*(r[i1-1] + r[i1+1]);
	else r[i1] = r[i1-1];
      } else {
	i = i1-1;
	r[i] += r1[i];
      }
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
      IntegrateSinCos(j, NULL, y, _phase, _dphase, i0, r, t);
      j = 0;
      for (i = i0; i <= i1; i += 2) {
	ip = i+1;
	y[j] = -y2[j];
	_phase[j] = large1[ip] - large2[ip];
	j++;
      }
      IntegrateSinCos(j, NULL, y, _phase, _dphasep, i0, r1, t);
      if (t < 0) {
	for (i = i0; i <= i1; i += 2) {
	  r[i] += r1[i];
	}
	if (i < MAX_POINTS) r[i] += r1[i];
	for (i = i0+1; i < i1; i += 2) {
	  r[i] = 0.5*(r[i-1] + r[i+1]);
	}
	if (i1 < MAX_POINTS-1) r[i1] = 0.5*(r[i1-1] + r[i1+1]);
	else r[i1] = r[i1-1];
      } else {
	i = i1-1;
	r[i] += r1[i];
      }
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
      IntegrateSinCos(j, x, y, _phase, _dphase, i0, r, t);
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
      IntegrateSinCos(j, x, y, _phase, _dphasep, i0, r1, t);
      if (t < 0) {
	for (i = i0; i <= i1; i += 2) {
	  r[i] += r1[i];
	}
	if (i < MAX_POINTS) r[i] += r1[i];
	for (i = i0+1; i < i1; i += 2) {
	  r[i] = 0.5*(r[i-1] + r[i+1]);
	}
	if (i1 < MAX_POINTS-1) r[i1] = 0.5*(r[i1-1] + r[i1+1]);
	else r[i1] = r[i1-1];
      } else {
	i = i1-1;
	r[i] += r1[i];
      }
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
      IntegrateSinCos(j, x, y, _phase, _dphase, i0, r, t);
      j = 0;
      for (i = i0; i <= i1; i += 2) {
	ip = i+1;
	x[j] = x1[j] - x2[j];
	x[j] *= 0.5*f[i];
	y[j] = -y2[j];
	_phase[j] = large1[ip] - large2[ip];
	j++;
      }
      IntegrateSinCos(j, x, y, _phase, _dphasep, i0, r1, t);
      if (t < 0) {
	for (i = i0; i <= i1; i += 2) {
	  r[i] += r1[i];
	}
	if (i < MAX_POINTS) r[i] += r1[i];
	for (i = i0+1; i < i1; i += 2) {
	  r[i] = 0.5*(r[i-1] + r[i+1]);
	}
	if (i1 < MAX_POINTS-1) r[i1] = 0.5*(r[i1-1] + r[i1+1]);
	else r[i1] = r[i1-1];
      } else {
	i = i1-1;
	r[i] += r1[i];
      }
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
      IntegrateSinCos(j, x, y, _phase, _dphase, i0, r, t);
      j = 0;
      for (i = i0; i <= i1; i += 2) {
	ip = i+1;
	x[j] = x1[j] + x2[j];
	x[j] *= 0.5*f[i];
	y[j] = -y2[j];
	_phase[j] = large1[ip] - large2[ip];
	j++;
      }
      IntegrateSinCos(j, x, y, _phase, _dphasep, i0, r1, t);
      if (t < 0) {
	for (i = i0; i <= i1; i += 2) {
	  r[i] += r1[i];
	}
	if (i < MAX_POINTS) r[i] += r1[i];
	for (i = i0+1; i < i1; i += 2) {
	  r[i] = 0.5*(r[i-1] + r[i+1]);
	}
	if (i1 < MAX_POINTS-1) r[i1] = 0.5*(r[i1-1] + r[i1+1]);
	else r[i1] = r[i1-1];
      } else {
	i = i1-1;
	r[i] += r1[i];
      }
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
      IntegrateSinCos(j, x, y, _phase, _dphase, i0, r, t);
      j = 0;
      for (i = i0; i <= i1; i += 2) {
	ip = i+1;
	x[j] = x1[j];
	y[j] = -y2[j];
	_phase[j] = large1[ip] - large2[ip];
	j++;
      }
      IntegrateSinCos(j, x, y, _phase, _dphasep, i0, r1, t);
      if (t < 0) {
	for (i = i0; i <= i1; i += 2) {
	  r[i] += r1[i];
	}
	if (i < MAX_POINTS) r[i] += r1[i];
	for (i = i0+1; i < i1; i += 2) {
	  r[i] = 0.5*(r[i-1] + r[i+1]);
	}
	if (i1 < MAX_POINTS-1) r[i1] = 0.5*(r[i1-1] + r[i1+1]);
	else r[i1] = r[i1-1];
      } else {
	i = i1-1;
	r[i] += r1[i];
      }
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
		    int i0, double *r, int t) {
  int i, k, m, n, q;
  double si0, si1, cs0, cs1;
  double is0, is1, is2, is3;
  double ic0, ic1, ic2, ic3;
  double a0, a1, a2, a3;
  double d, p;
  double *z, *u;
  double h, dr;

  z = _dwork10;
  u = _dwork11;
  for (i = 1, k = i0+2; i < j; i++, k += 2) {
    h = dphase[i-1]+dphase[i];
    dr = potential->rad[k] - potential->rad[k-2];
    if (h*dr > 0.1) break;
    z[i] = 0.0;
    if (x != NULL) z[i] += x[i]*sin(phase[i]);
    if (y != NULL) z[i] += y[i]*cos(phase[i]);
    z[i] *= potential->dr_drho[k];
  }
  if (i > 1) {
    z[0] = 0.0;
    if (x != NULL) z[0] += x[0]*sin(phase[0]);
    if (y != NULL) z[0] += y[0]*cos(phase[0]);
    z[0] *= potential->dr_drho[i0];
    u[0] = 0.0;
    NewtonCotes(u, z, 0, i-1, t);
    for (m = 1, n = i0+2; m < i; m++, n += 2) {
      r[n] = r[i0] + 2.0*u[m];
    }
  }
 
  q = i-1; 
  m = j-q;
  if (m < 2) {
    if (k < MAX_POINTS) {
      r[k] = r[k-2];
    }
    return 0;
  }

  if (x != NULL) {
    for (n = q; n < j; n++) {
      x[n] /= dphase[n];
    }    
    spline(phase+q, x+q, m, 1E30, 1E30, z+q);   
  } 
  if (y != NULL) {
    for (n = q; n < j; n++) {
      y[n] /= dphase[n];
    }
    spline(phase+q, y+q, m, 1E30, 1E30, u+q);      
  } 
  
  si0 = sin(phase[i-1]);
  cs0 = cos(phase[i-1]);    
  d = phase[i] - phase[i-1];
  si1 = sin(phase[i]);
  cs1 = cos(phase[i]);
  is0 = -(cs1 - cs0);
  ic0 = si1 - si0;
  is1 = -d * cs1 + ic0;
  ic1 = d * si1 - is0;
  r[k] = r[k-2];
  if (x != NULL) {
    a1 = (x[i] - x[i-1])/d;
    a0 = x[i-1];
    h = a0*is0 + a1*is1;
    r[k] += h;
  } 
  if (y != NULL) {
    a1 = (y[i] - y[i-1])/d;
    a0 = y[i-1];
    h = a0*ic0 + a1*ic1;
    r[k] += h;    
  }
  si0 = si1;
  cs0 = cs1;
  i++;
  k += 2;  

  for (; i < j; i++, k += 2) {
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
      a3 = (z[i] - z[i-1])/(6.0*d);
      a2 = z[i-1]/3.0;
      a1 = (x[i] - x[i-1])/d - (z[i] + z[i-1])*d/6.0;
      a0 = x[i-1];
      h = a0*is0 + a1*is1 + a2*is2 + a3*is3;
      r[k] += h;
    } 
    if (y != NULL) {
      a3 = (u[i] - u[i-1])/(6.0*d);
      a2 = u[i-1]/3.0;
      a1 = (y[i] - y[i-1])/d - (u[i] + u[i-1])*d/6.0;
      a0 = y[i-1];
      h = a0*ic0 + a1*ic1 + a2*ic2 + a3*ic3;
      r[k] += h;    
    }
    si0 = si1;
    cs0 = cs1;
  }
  return 0;
}

int FreeSimpleArray(MULTI *ma) {
  if (ma->array == NULL) return 0;
  MultiFreeData(ma->array, ma->ndim, NULL);
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

int FreeMultipoleArray(void) {
  if (multipole_array->array == NULL) return 0;
  MultiFreeData(multipole_array->array, multipole_array->ndim, 
		FreeMultipole);
  return 0;
}

int FreeMomentsArray(void) {
  if (moments_array->array == NULL) return 0;
  MultiFreeData(moments_array->array, moments_array->ndim, NULL);
  return 0;
}

int FreeGOSArray(void) {
  if (gos_array->array == NULL) return 0;
  MultiFreeData(gos_array->array, gos_array->ndim, FreeMultipole);
  return 0;
}

int InitRadial(void) {
  int i, ndim;
  int blocks[6] = {4, 4, 4, 4, 4, 1};

  potential = malloc(sizeof(POTENTIAL));
  potential->flag = 0;
  potential->uehling[0] = 0.0;
  n_orbitals = 0;
  n_continua = 0;
  
  orbitals = malloc(sizeof(ARRAY));
  if (!orbitals) return -1;
  if (ArrayInit(orbitals, sizeof(ORBITAL), ORBITALS_BLOCK) < 0) return -1;

  ndim = 6;
  slater_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(slater_array, sizeof(double), ndim, blocks);

  ndim = 5;
  breit_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(breit_array, sizeof(double), ndim, blocks);
  
  ndim = 2;
  residual_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(residual_array, sizeof(double), ndim, blocks);

  vinti_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(vinti_array, sizeof(double), ndim, blocks);

  qed1e_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(qed1e_array, sizeof(double), ndim, blocks);
  
  ndim = 4;
  multipole_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(multipole_array, sizeof(double *), ndim, blocks);

  ndim = 3;
  moments_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(moments_array, sizeof(double), ndim, blocks);

  gos_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(gos_array, sizeof(double *), ndim, blocks);

  n_awgrid = 1;
  awgrid[0]= EPS3;
  
  SetRadialGrid(1E-5, 5E2);

  return 0;
}

int ReinitRadial(int m) {
  int i;

  if (m < 0) return 0;
  ClearOrbitalTable(m);
  FreeSimpleArray(slater_array);
  FreeSimpleArray(breit_array);
  FreeSimpleArray(residual_array);
  FreeSimpleArray(qed1e_array);
  FreeSimpleArray(vinti_array);
  FreeMultipoleArray();
  FreeMomentsArray();
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
      SetRadialGrid(1E-5, 5E2);
      potential->uehling[0] = 0.0;
    }
  }
  return 0;
}
  
int TestIntegrate(char *s) {
  ORBITAL *orb1, *orb2, *orb3;
  int k1, k2, k3, i;
  double r;
  FILE *f;
  
  orb1 = GetOrbital(0);
  orb2 = GetOrbital(3);
  printf("%d %d %d %d \n", orb1->n, orb1->kappa, orb2->n, orb2->kappa);
  GetYk(1, _yk, orb1, orb2, -2); 
  
  for (i = 0; i < MAX_POINTS; i++) {  
    _yk[i] /= potential->rad[i];  
  }  
  i = 15;
  r = 1e4/27.2;
  k1 = OrbitalIndex(0, i-1, 6.8E2/27.2 +r);    
  k2 = OrbitalIndex(0, i, r);    
  k3 = OrbitalIndex(0, i+1, 6.8E2/27.2 +r);     
  
  orb1 = GetOrbital(k1); 
  orb2 = GetOrbital(k2);   
  orb3 = GetOrbital(k3);   
 
  Integrate(_yk, orb1, orb2, 2, &r);
  printf("%10.3E\n", r);
  Integrate(_yk, orb1, orb2, -2, _xk);  
  
  f = fopen(s, "w"); 
  for (i = 0; i < MAX_POINTS; i++) { 
    fprintf(f, "%d %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E\n",  
	    i, potential->rad[i],   
	    _yk[i], _zk[i], _xk[i], potential->Vc[i], 
	    potential->U[i], 
	    Large(orb1)[i], Large(orb2)[i],  
	    Large(orb3)[i]); 
  }

  Integrate(_yk, orb3, orb2, 2, &r);
  printf("%10.3E\n", r);
  Integrate(_yk, orb3, orb2, -2, _xk);  
  fprintf(f, "\n\n\n"); 
  for (i = 0; i < MAX_POINTS; i++) { 
    fprintf(f, "%d %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E\n",  
	    i, potential->rad[i],   
	    _yk[i], _zk[i], _xk[i], potential->Vc[i], 
	    potential->U[i], 
	    Large(orb1)[i], Large(orb2)[i],  
	    Large(orb3)[i]); 
  }


  fclose(f); 
  return 0;
}
