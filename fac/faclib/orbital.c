#include "orbital.h"

static char *rcsid="$Id: orbital.c,v 1.28 2002/04/25 16:22:27 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

/* closed Newton-Cotes formulae coeff. */
static double _CNC[5][5] = {
  {0, 0, 0, 0, 0},
  {0.5, 0.5, 0, 0, 0},
  {1.0/3.0, 4.0/3.0, 1.0/3.0, 0, 0},
  {3.0/8, 9.0/8, 9.0/8, 3.0/8, 0},
  {14.0/45, 64.0/45, 24.0/45, 64.0/45, 14.0/45,}
};

/* open Newton-Cotes formulae coeff. */
static double _ONC[9][9] = {
  {0, 0, 0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0, 0, 0},
  {0, 2.0, 0, 0, 0, 0, 0, 0, 0},
  {0, 1.5, 1.5, 0, 0, 0, 0, 0, 0},
  {0, 8./3, -4./3, 8./3, 0, 0, 0, 0, 0},
  {0, 55./24, 5./24, 5./24, 55./24, 0, 0, 0, 0},
  {0, 33./10, -42./10, 78./10, -42./10, 33./10, 0, 0, 0},
  {0, 4277./1440, -3171./1440, 3934./1440, 3934./440, 
   -3171./1440, 4277./1440, 0, 0},
  {0, 3680./945, -7632./945, 17568./945, -19672./945,
   17568./945, -7632./945, 3680./945, 0}
};

/* the following arrays provide storage space in the calculation */
static double _veff[MAX_POINTS];
static double _ACOEF[MAX_POINTS];
static double _BCOEF[MAX_POINTS];
static double _VCOEF[MAX_POINTS];
static double _dwork[MAX_POINTS];
static double _dwork1[MAX_POINTS];
static double _dwork2[MAX_POINTS];
static double _dwork3[MAX_POINTS];
static double _dwork4[MAX_POINTS];
 
static int max_iteration = 512;
static int nmax = 0;
static double wave_zero = 1E-10;

static int _SetVEffective(int kl, POTENTIAL *pot);
static int _TurningPoints(int n, double e, int *i1, int *i2, POTENTIAL *pot);
static int _Outward(double *p, double e, POTENTIAL *pot, int i1, int *k2);
static int _Inward(double *p, double e, POTENTIAL *pot, int i2);
static int _Amplitude(double *p, double e, int kl, POTENTIAL *pot, int i1);
static int _AmplitudeWKB(double *p, double e, int kl, POTENTIAL *pot, int i1);
static int _Phase(double *p, POTENTIAL *pot, int i1, double p0);
static int _DiracSmall(ORBITAL *orb, POTENTIAL *pot);

double dlogam_(double *x);
void y5n_(double *lambda, double *eta0, double *x0, 
	  double *y5, double *y5p, double *norm, int *ierr);
void cphamp_(double *lambda, double *eta0, double *x0, 
	     double *a0, double *cphase0, double *a0p, int *ierr);

int GetNMax(void) {
  return nmax;
}
 
int RadialSolver(ORBITAL *orb, POTENTIAL *pot, double tol) {
  int ierr;
  
  if (orb->n > 0) {
    if (orb->n < nmax) {
      ierr = RadialBound(orb, pot, tol);
    } else {
      ierr = RadialRydberg(orb, pot, tol);
    }
  } else {
    ierr = RadialFree(orb, pot);
  }  
  return ierr;
}

double *GetVEffective(void) { 
  return _veff;
}

static double EstimateBoundEnergy(ORBITAL *orb, POTENTIAL *pot) {
  int i;
  double a, b;
  a = orb->n;
  a = 1.0/(a*a);
  for (i = MAX_POINTS-1; i > 0; i--) {
    b = (pot->Vc[i] + pot->U[i])*pot->rad[i];
    if (b*a*pot->rad[i] < 1.0) break;
  }
  b = -0.5*b*b*a;
  return b;
}
  
int RadialBound(ORBITAL *orb, POTENTIAL *pot, double tol) {
  double z, e, e0, emin, emax;
  int i, kl, nr, nodes, nodes_old, niter;
  int i1, i2, i2p, i2m, i2p2, i2m2, ierr;
  double *p, p1, p2, qi, qo, delta, ep, norm2, fact, eps;
  
  z = (pot->Z[MAX_POINTS-1] - pot->N + 1.0);
  e = EstimateBoundEnergy(orb, pot);

  kl = orb->kappa;
  if (pot->flag < 0) {
    kl = (kl < 0)? (-kl-1):kl;
  }
  if (kl < 0 || kl >= orb->n) {
      printf("Invalid orbital angular momentum, L=%d\n", kl);
      return -1;
  }

  if (pot->flag < 0) {
    SetPotentialW(pot, e, orb->kappa); 
    p = malloc(sizeof(double)*2*MAX_POINTS);
    if (!p) return -1;
  } else { 
    p = malloc(sizeof(double)*MAX_POINTS);
    if (!p) return -1;
  }
  _SetVEffective(kl, pot);

  emin = e*2.0;
  for (i = 0; i < MAX_POINTS; i++) {
    if (emin > _veff[i]) emin = _veff[i];
  }
     
  emax = EPS10 * e;
  if (emin > emax) return 1;
  emin = 1.1*emin;
  if (e >= emax) e = emax*1.1;
  if (e <= emin) e = emin*0.9;
  ierr = 0;
  nr = orb->n - kl - 1;
  niter = -1;
  eps = 1.0;
  nodes_old = nr;
  while (1) {
    niter++;
    if (e > emax || e < emin) {
      printf("E exceeded the limits at iteration: %d\n", niter);
      printf("e = %10.3E, emin = %10.3E, emax = %10.3E\n", e, emin, emax);
      ierr = 2;
      break;
    }
    ierr = _TurningPoints(orb->n, e, &i1, &i2, pot);
    if (ierr == -1) {
      break;
    } else if (ierr == -2) {
      e *= 0.9; 
      niter--;
      continue;
    }
   
    nodes = _Outward(p, e, pot, i1, &i2);
    
    if (nodes != nr) {
      e0 = e;
      if (nodes > nr) {
	if (nodes_old > nr) {
	  eps *= 1.2;
	} else {
	  eps *= 0.8;
	}
	e = e0*(1.0+(eps*(nodes-nr))/orb->n);
      } else {
	if (nodes_old < nr) {
	  eps *= 1.2;
	} else {
	  eps *= 0.8;
	}
	e = e0*(1.0+(eps*(nodes-nr))/orb->n);
      }
      if (e >= emax) {
	e = e0 + 0.2*(emax - e0);
	eps = 1.0;
      }
      if (e <= emin) {
	e = e0 + 0.2*(emin - e0);
	eps = 1.0;
      }
      nodes_old = nodes;
      if (niter > max_iteration) {
	printf("MAX iteration reached in Bound before ");
	printf("finding the solution with correct nodes\n");
	printf("nodes expected: %d, found: %d\n", nr, nodes);
	printf("probably needs to extend the radial ");
	printf("grid to larger distances\n");
	printf("%10.3E %10.3E %10.3E %10.3E\n", e, emin, emax, eps);
	ierr = -3;
	break;
      }
      if (pot->flag < 0) {
	SetPotentialW(pot, e, orb->kappa);
	_SetVEffective(kl, pot);
      }
      continue;
    } 
    i2m = i2 - 1;
    i2p = i2 + 1;
    i2p2 = i2 + 2;
    i2m2 = i2 - 2;
    qo = (-4.0*p[i2m2-1] + 30.0*p[i2m2] - 120.0*p[i2m]
	  + 40.0*p[i2] + 60.0*p[i2p] - 6.0*p[i2p2])/120.0;
    ierr = _Inward(p, e, pot, i2);
    if (ierr) break;
    p1 = (p[i2]*_BCOEF[i2] - p[i2p]*_ACOEF[i2p])/_ACOEF[i2m];
    p2 = (p[i2m]*_BCOEF[i2m] - p[i2]*_ACOEF[i2])/_ACOEF[i2m2];
    qi = (6.0*p2 - 60.0*p1 - 40.0*p[i2] + 120.0*p[i2p]
	  - 30.0*p[i2p2] + 4.0*p[i2p2+1])/120.0;

    for (i = 0; i < MAX_POINTS; i++) {
      p[i] = p[i] * sqrt(pot->dr_drho[i]);
    }
    fact = 1.0 / sqrt(pot->dr_drho[i2]);
    qo = qo * fact;
    qi = qi * fact;
    
    norm2 = InnerProduct(MAX_POINTS, p, p, pot);
    
    delta = 0.5*p[i2]*(qo - qi)/(norm2);
    ep = e+delta;
    if (ep >= emax) {
      ep = e + 0.5*(emax - e);
    }
    if (ep <= emin) {
      ep = e + 0.5*(emin - e);
    }
    e0 = e;
    e = ep;
    ep = Min(tol, fabs(ep)*tol);
    if ((fabs(delta) < ep || fabs(e-e0) < ep)) break;
    if (niter > max_iteration) {
      printf("MAX iteration reached in Bound.\n");
      printf("probably needs to extend the radial ");
      printf("grid to larger distances.\n");
      printf("%10.3E %10.3E %10.3E %10.3E\n", e, e0, emin, delta);
      ierr = -4;
      break;
    }
    if (pot->flag < 0) {
      SetPotentialW(pot, e, orb->kappa);
      _SetVEffective(kl, pot);
    }
  }
  if (ierr) {
    free(p);
    return ierr;
  }

  qi = sqrt(norm2);
  fact = 1.0/qi;
  if ((IsOdd(nodes) && p[i2] > 0.0) || (IsEven(nodes) && p[i2] < 0.0)) {
    fact = -fact;
  }
 
  qi *= wave_zero;     
  for (i = MAX_POINTS-1; i >= 0; i--) {
    if (fabs(p[i]) > qi) break;
  }
  if (IsEven(i)) i++;
  orb->ilast = i;
       
  for (i = 0; i < MAX_POINTS; i++) {    
    p[i] *= fact;
  }
 
  orb->energy = e;
  orb->wfun = p;

  orb->qr_norm = 1.0;
  if (pot->flag < 0) {
    _DiracSmall(orb, pot);
  } 
 
  return 0;
}

int RadialRydberg(ORBITAL *orb, POTENTIAL *pot, double tol) {
  double z, e, e0, emin, emax;
  int i, kl, niter, ierr;
  double *p;
  double lambda, eta0, x0, y5, y5p, y5norm;
  int i1, i2, i2p, i2m, i2p2, i2m2, nodes;
  double qo, qi, norm2, delta, dk, zp;

  z = (pot->Z[MAX_POINTS-1] - pot->N + 1.0);
  kl = orb->kappa;
  if (pot->flag < 0) {
    kl = (kl < 0)? (-kl-1):kl;
  }
  if (kl < 0 || kl >= orb->n) {
      printf("Invalid orbital angular momentum, L=%d\n", kl);
      return -1;
  }
  emin = z/(orb->n - 0.8);
  emin = -emin*emin*0.5;
  emax = z/(orb->n + 0.2);
  emax = -emax*emax*0.5;
  eta0 = orb->n-0.05;
  dk = z/eta0;
  e = -dk*dk*0.5;
  lambda = kl;
  if (pot->flag < 0) {
    zp = FINE_STRUCTURE_CONST2*e;
    e *= (1.0 + zp)/(1.0+0.5*zp);
    dk *= sqrt(1.0+0.5*zp);
    e0 = lambda*(lambda+1.0) - FINE_STRUCTURE_CONST2*z*z;
    lambda = (sqrt(1.0+4.0*e0)-1.0)*0.5;
    zp = z*(1.0+zp);
    SetPotentialW(pot, e, orb->kappa); 
    p = malloc(sizeof(double)*2*MAX_POINTS);
    if (!p) return -1;
  } else {
    zp = z;
    p = malloc(sizeof(double)*MAX_POINTS);
    if (!p) return -1;
  }
  _SetVEffective(kl, pot);

  niter = 0;
  while (niter < max_iteration) {
    niter++;
    ierr = _TurningPoints(orb->n, e, &i1, &i2, pot);
    nodes = _Outward(p, e, pot, i1, &i2);
    x0 = dk*pot->rad[i2];
    y5n_(&lambda, &eta0, &x0, &y5, &y5p, &y5norm, &ierr);
    qi = dk*(y5p/y5);
    i2p = i2 + 1;
    i2m = i2 - 1;
    i2p2 = i2 + 2;
    i2m2 = i2 - 2;
    p[i2m2] *= sqrt(pot->dr_drho[i2m2]);
    p[i2m] *= sqrt(pot->dr_drho[i2m]);
    p[i2] *= sqrt(pot->dr_drho[i2]);
    p[i2p2] *= sqrt(pot->dr_drho[i2p2]);
    p[i2p] *= sqrt(pot->dr_drho[i2p]);
    qo = (-4.0*p[i2m2-1] + 30.0*p[i2m2] - 120.0*p[i2m]
	  + 40.0*p[i2] + 60.0*p[i2p] - 6.0*p[i2p2])/120.0;
    qo /= (p[i2]*pot->dr_drho[i2]);
    delta = qo - qi;
    e0 = eta0 - lambda;
    norm2 = dlogam_(&e0);
    e0 = eta0 + lambda + 1.0;
    norm2 += dlogam_(&e0);
    e0 = zp/(eta0*eta0);
    norm2 = -norm2 + log(e0);
    norm2 = 0.5*norm2+y5norm; 
    norm2 = exp(norm2)*y5;
    delta *= 0.5*(norm2*norm2);
    e0 = -tol*e;
    if (emax - emin < e0) break;
    if (delta > e0) {
      emin = e;
      e += delta;
      if (e > emax) e = 0.5*(emin+emax);
    } else if (delta < -e0) {
      emax = e;
      e += delta;
      if (e < emin) e = 0.5*(emin+emax);
    } else {
      break;
    }
    if (pot->flag < 0) {      
      SetPotentialW(pot, e, orb->kappa);
      _SetVEffective(kl, pot);
      zp = FINE_STRUCTURE_CONST2*e;
      dk = sqrt(-2.0*e*(1.0+0.5*zp));
      zp = z*(1.0+zp);
      eta0 = zp/dk;
    } else {
      dk = sqrt(-2.0*e);
      zp = z;
      eta0 = zp/dk;
    }
  }
  if (niter == max_iteration) {
    printf("MAX iteration reached in RadialRydberg, for N=%d, Kappa=%d\n",
	   orb->n, orb->kappa);
    printf("The wavefunction may be inaccurate\n");
  } 
  for (i = MAX_POINTS-1; i > i2+10; i--) {
    if (e > _veff[i]) break;
  }
  i2p = i;
  nodes = _Outward(p, e, pot, i1, &i2p);
  if (i2p < MAX_POINTS-1) ierr = _Inward(p, e, pot, i2p);
  for (i = 0; i < MAX_POINTS; i++) {
    p[i] *= sqrt(pot->dr_drho[i]);
  }
  qi = norm2/(p[i2]);
  for (i = 0; i < MAX_POINTS; i++) {
    p[i] *= qi;
  }

  for (i = MAX_POINTS-1; i >= 0; i--) {
    if (fabs(p[i]) > wave_zero) break;
  }
  if (IsEven(i)) i++;
  orb->ilast = i;
  orb->energy = e;
  orb->wfun = p;

  if (pot->flag < 0) {
    e0 = InnerProduct(MAX_POINTS, p, p, pot);
    orb->qr_norm = 1.0/e0;
    _DiracSmall(orb, pot);
  } else 
    orb->qr_norm = 1.0;

  return 0;
}  
  

/* note that the free states are normalized to have asymptotic 
   amplitude of 1/sqrt(k), */
int RadialFree(ORBITAL *orb, POTENTIAL *pot) {
  int i, kl, nodes;
  int i1, i2, i2p, i2m, i2p2, i2m2;
  double *p, po, qo, pm, e;
  double dfact, da, cs, si, phase0;

  e = orb->energy;
  if (e < 0.0) { 
    printf("Energy < 0 in Free\n");
    return -1;
  }
  kl = orb->kappa;
  if (pot->flag < 0) {
    if (orb->kappa == 0) {
      printf("Kappa == 0 in Free\n");
      return -1;
    }
    SetPotentialW(pot, e, kl);
    kl = (kl < 0)? (-kl-1):kl;  
    p = malloc(2*MAX_POINTS*sizeof(double));
    if (!p) return -1;
  } else { 
    if (kl < 0) {
      printf("L < 0 in Free\n");
      return -1;
    }
    p = malloc(sizeof(double)*MAX_POINTS);
    if (!p) return -1;
  }
  _SetVEffective(kl, pot);

  i = _TurningPoints(0, e, &i1, &i2, pot);
  nodes = _Outward(p, e, pot, i1, &i2);
  i2m = i2 - 1;
  i2p = i2 + 1;
  i2m2 = i2 - 2;
  i2p2 = i2 + 2;

  for (i = i2p2; i >= 0; i--) {
    p[i] *= sqrt(pot->dr_drho[i]);
  }
  dfact = 1.0 / pot->dr_drho[i2];
  qo = (-4.0*p[i2m2-1] + 30.0*p[i2m2] - 120.0*p[i2m]
	+ 40.0*p[i2] + 60.0*p[i2p] - 6.0*p[i2p2])/120.0;
  qo *= dfact;
  po = p[i2];

  pm = p[i2m];
  
  _Amplitude(p, e, kl, pot, i2m);

  da = (-24.0*p[i2m] - 130.0*p[i2] + 240.0*p[i2p] - 120.0*p[i2p2] 
	+ 40.0*p[i2p2+1] - 6.0*p[i2p2+2])/(120.0); 
  da *= dfact;
  cs = p[i2] * qo - da*po;
  si = po / p[i2];
  dfact = (si*si + cs*cs);
  dfact = 1.0/sqrt(dfact);
  phase0 = atan2(si, cs);
  if (phase0 < 0) phase0 += TWO_PI;
    
  if (IsOdd(nodes)) {
    phase0 = (phase0 < PI)?(phase0 + PI):(phase0-PI);
    dfact = -dfact;
  }
  
  p[i2m] = pm;
  for (i = 0; i < i2; i++) {
    p[i] *= dfact;
  }
    
  _Phase(p, pot, i2, phase0);

  orb->ilast = i2m;
  orb->wfun = p;
  orb->phase = NULL;

  orb->qr_norm = 1.0;
  if (pot->flag < 0) {
    _DiracSmall(orb, pot);
  } 
  return 0;
}

int _DiracSmall(ORBITAL *orb, POTENTIAL *pot) {
  int i, i1, kappa;
  double xi, e, *p, a, b;

  e = orb->energy;
  kappa = orb->kappa;
  p = orb->wfun;
  i1 = orb->ilast+1;

  for (i = 0; i < i1; i++) {
    xi = e - pot->Vc[i] - pot->U[i];
    xi = xi*FINE_STRUCTURE_CONST2*0.5;
    _dwork[i] = 1.0 + xi;
    _dwork1[i] = sqrt(_dwork[i])*p[i];
    _dwork2[i] = 1.0/(24.0*pot->dr_drho[i]);
  }
  
  for (i = 0; i < i1; i++) {
    if (fabs(_dwork1[i]) < 1E-16) {
      p[i+MAX_POINTS] = 0.0;
      continue;
    } 
    if (i == 0) {
      b = -50.0*_dwork1[0];
      b += 96.0*_dwork1[1];
      b -= 72.0*_dwork1[2];
      b += 32.0*_dwork1[3];
      b -= 6.0 *_dwork1[4];
    } else if (i == 1) {
      b = -6.0*_dwork1[0];
      b -= 20.0*_dwork1[1];
      b += 36.0*_dwork1[2];
      b -= 12.0*_dwork1[3];
      b += 2.0 *_dwork1[4];
    } else if (i == i1-1) {
      b = -50.0*_dwork1[i];
      b += 96.0*_dwork1[i-1];
      b -= 72.0*_dwork1[i-2];
      b += 32.0*_dwork1[i-3];
      b -= 6.0 *_dwork1[i-4];
      b = -b; 
    } else if (i == i1-2) {
      b = -6.0*_dwork1[i+1];
      b -= 20.0*_dwork1[i];
      b += 36.0*_dwork1[i-1];
      b -= 12.0*_dwork1[i-2];
      b += 2.0 *_dwork1[i-3];
      b = -b;
    } else {
      b = 2.0*(_dwork1[i-2] - _dwork1[i+2]);
      b += 16.0*(_dwork1[i+1] - _dwork1[i-1]);
    }
    b *= _dwork2[i];
    b += _dwork1[i]*kappa/pot->rad[i];
    b /= (2.0*_dwork[i]);
    p[i] = _dwork1[i];
    p[i+MAX_POINTS] = b*FINE_STRUCTURE_CONST;
  } 
  if (orb->n > 0) {
    for (i = i1; i < MAX_POINTS; i++) {
      xi = e - pot->Vc[i] - pot->U[i];
      xi = xi*FINE_STRUCTURE_CONST2*0.5;
      _dwork[i] = 1.0 + xi;
      p[i] = sqrt(_dwork[i])*p[i];
    }
    a = InnerProduct(i1, p+MAX_POINTS, p+MAX_POINTS, pot);
    b = InnerProduct(MAX_POINTS, p, p, pot);    
    a *= orb->qr_norm;
    b *= orb->qr_norm;
    a = sqrt(a+b);
    orb->qr_norm = a/sqrt(b);
    a = 1.0/a;
    for (i = 0; i < i1; i++) {
      p[i] *= a;
      p[i+MAX_POINTS] *= a;
    }
    for (i = i1+MAX_POINTS; i < 2*MAX_POINTS; i++) {
      p[i] = 0.0;
    }
    return 0;
  }
  
  for (i = i1; i < MAX_POINTS; i += 2) {
    xi = e - pot->Vc[i] - pot->U[i];
    xi = xi*FINE_STRUCTURE_CONST2*0.5;
    _dwork[i] = 1.0 + xi;
    _dwork1[i] = sqrt(_dwork[i])*p[i];
    _dwork2[i] = 0.25/(pot->dr_drho[i]);
  }

  for (i = i1; i < MAX_POINTS; i += 2) {
    if (i == i1) {
      b = -3.0*_dwork1[i] + 4.0*_dwork1[i+2] - _dwork1[i+4];
    } else if (i == MAX_POINTS-2) {
      b = -3.0*_dwork1[i] + 4.0*_dwork1[i-2] - _dwork1[i-4];
      b = -b;
    } else {
      b = _dwork1[i+2] - _dwork1[i-2];
    }
    b *= _dwork2[i];
    b += _dwork1[i]*kappa/pot->rad[i];
    a = _dwork1[i]/(p[i]*p[i]);
    p[i] = _dwork1[i];
    p[i+MAX_POINTS] = FINE_STRUCTURE_CONST*a/(2.0*_dwork[i]); 
    p[i+1+MAX_POINTS] = FINE_STRUCTURE_CONST*b/(2.0*_dwork[i]);
  }

  b = FINE_STRUCTURE_CONST2*orb->energy;
  orb->qr_norm = sqrt((1.0 + b)/(1.0 + 0.5*b));
  return 0;
}
  
int _AmplitudeWKB(double *p, double e, int kl, 
		  POTENTIAL *pot, int i1) {
  int i, j;
  double x, y, a, b, kl1, r, r2;

  kl1 = kl*(kl+1.0);
  for (i = i1; i < MAX_POINTS; i++) {
    x = 1.0/(2.0*(e - _veff[i]));    
    y = pow(x, 0.25);
    _dwork[i] = y;
    r = pot->rad[i];
    r2 = r*r*r;
    a = pot->dVc[i] + pot->dU[i] - kl1/r2;
    b = pot->dVc2[i] + pot->dU2[i] + 3.0*kl1/(r2*r);
    if (pot->flag < 0) {
      a += pot->dW[i];
      b += pot->dW2[i];
    }
    _dwork1[i] = 0.5*x*y*(2.5*x*a*a + b);
    _dwork2[i] = 0.0;
  }
  
  for (i = i1; i < MAX_POINTS; i++) {
    for (j = 0; j < 4; j++) {
      y = _dwork[i] + _dwork2[i];
      b = _dwork1[i];
      a = 2.0*(e - _veff[i]) + b/y;
      a = pow(a, -0.25) - _dwork[i];
      _dwork2[i] = a;
    }
  }

  for (i = i1; i < MAX_POINTS; i++) {
    p[i] = _dwork[i] + _dwork2[i];
  }

  return 0;
}
  
int _Amplitude(double *p, double e, int kl, 
	       POTENTIAL *pot, int i1) {
  int i, ierr;
  double a, b, r, z, zp, dk2, x0t;
  double a0, dk, eta0, x0, lambda, eps1, eps2;

  z = pot->Z[MAX_POINTS-1] - pot->N + 1.0;
  if (pot->flag < 0) {
    eps1 = FINE_STRUCTURE_CONST2*e;
    eps2 = FINE_STRUCTURE_CONST2*z*z;
    a = kl*(kl+1.0)-eps2;
    lambda = 0.5*(sqrt(1.0+4.0*a) - 1.0);
    dk = sqrt(2.0*e*(1.0+0.5*eps1));
    zp = z*(1.0+eps1);
  } else {
    zp = z;
    dk = sqrt(2.0*e);
    lambda = kl;
  }
  eta0 = zp/dk;
  dk2 = 1.0/sqrt(dk);
  eps1 *= e;
  if (kl > 0) {
    x0t = 0.5*(-eta0 + sqrt(eta0*eta0 + lambda*(lambda+1.0)));
  } else {
    x0t = dk*(pot->rad[i1+5]-pot->rad[i1]);
  }
  x0t += dk*pot->rad[i1];
  for (i = i1; i < MAX_POINTS; i++) {
    r = pot->rad[i];
    x0 = dk*r;
    if (x0 > x0t) break;
    ierr = 1;
    cphamp_(&lambda, &eta0, &x0, &a0, NULL, NULL, &ierr);
    a0 *= dk2;
    a = -2.0*((pot->Vc[i]+pot->U[i])+zp/r);
    if (pot->flag < 0) {
      b = -eps1;
      b -= 2.0*pot->W[i];
      b -= eps2/(r*r);
    } else {
      b = 0.0;
    }
    a0 = pow(a0, -4.0);
    a0 += a + b; 
    p[i] = pow(a0, -0.25);
  }
  if (i < MAX_POINTS) _AmplitudeWKB(p, e, kl, pot, i);

  return 0;
}

int _Phase(double *p, POTENTIAL *pot, int i1, double phase0) {
  int i;
  double fact;
  
  fact = 1.0 / 3.0;

  for (i = i1; i < MAX_POINTS; i++) {
    _dwork[i] = 1.0/p[i];
    _dwork[i] *= _dwork[i];
    _dwork[i] *= pot->dr_drho[i];
  }

  i = i1+1;
  p[i] = phase0;
  for (i = i1+3; i < MAX_POINTS; i += 2) {
    p[i] = p[i-2] + (_dwork[i-3] + 4.0*_dwork[i-2] + _dwork[i-1])*fact;
  }
  return 0;
}

int _SetVEffective(int kl, POTENTIAL *pot) {
  double kl1;
  int i;
  double r;

  kl1 = 0.5*kl*(kl+1);
 
  for (i = 0; i < MAX_POINTS; i++) {
    r = pot->rad[i];
    r *= r;
    _veff[i] = pot->Vc[i] + pot->U[i] + kl1/r;
    if (pot->flag < 0) {
      _veff[i] += pot->W[i];
    }
  }

  return 0;
}

int _TurningPoints(int n, double e, int *i1, int *i2, POTENTIAL *pot) {
  int i;
  double x, a, b;

  if (n <= 0) {
    for (i = MAX_POINTS - 5; i > 10; i--) {
      if (e < _veff[i]) break;
    }

    *i1 = i-3;
    i += 5;
    for (; i <= pot->r_core; i++) {
      x = e - _veff[i];
      b = 1.0/pot->rad[i];
      a = 20.0/(0.5*pot->ar*sqrt(b) + pot->br*b);
      x = TWO_PI/sqrt(2.0*x);
      if (x < a) break;
    }
    *i2 = i-1;
    if (*i1 <= 10) *i1 = *i2 - 5;
  } else {
    for (i = 10; i < MAX_POINTS; i++) {
      if (e > _veff[i]) break;
    }
    *i1 = i; 

    for (i = MAX_POINTS-10; i > 10; i--) {
      if (e > _veff[i]) break;
    }
    if (i <= 10) return -2;
    if (n >= nmax) {
      *i2 = i - 16;
      i = pot->r_core+16;
      i = Min(i, MAX_POINTS-10);
      if (*i2 < i) *i2 = i;
      if (*i1 > *i2 - 8) {
	*i1 = *i2 - 8;
      }
    } else {
      *i2 = i + 7;
      if (*i1 <= 10) *i1 = *i2 - 8; 
    }
  } 

  return 0;
}
     
int _Outward(double *p, double e, POTENTIAL *pot, int i1, int *k2) {
  int i, nodes, i2;
  double a, b, r, x, y, z, p0;

  i2 = *k2;
  for (i = 0; i <= i2; i++) {
    r = pot->rad[i];
    x = 2.0*(_veff[i] - e);
    x *= 4.0*r*r;
    a = pot->ar;
    b = pot->br;
    z = sqrt(r);
    y = a*z + 2.0*b;
    y = 1/(y*y);
    z = (0.75*a*a*r + 5.0*a*b*z +4.0*b*b);
    x += z*y;    
    x *= y / 12.0;
    _ACOEF[i] = 1.0 - x;
    _BCOEF[i] = 2.0 + 10.0*x;
  }

  _VCOEF[0] = 0.0;
  for (i = 1; i < i1; i++) {
    a = _BCOEF[i] - _ACOEF[i-1]*_VCOEF[i-1];
    _VCOEF[i] = _ACOEF[i+1] / a;
  }

  p[0] = 0.0;
  p[i1] = 1.0;
  nodes = 0;
  p0 = p[i1];

  for (i = i1-1; i > 0; i--) {
    p[i] = _VCOEF[i]*p[i+1];
    a = fabs(p[i]);
    if (e >= _veff[i] && a > wave_zero) {
      if ((p0 > 0 && p[i] < 0) ||
	  (p0 < 0 && p[i] > 0)) {
	nodes++;	
	p0 = p[i];
      }
    } 
    
  }
  p0 = p[i1];
  for (i = i1+1; i <= i2; i++) {
    p[i] = (_BCOEF[i-1]*p[i-1] - _ACOEF[i-2]*p[i-2])/_ACOEF[i];
    if (e < 0.0 && e >= _veff[i] && fabs(p[i]) > wave_zero) {
      if ((p0 > 0 && p[i] < 0) ||
	  (p0 < 0 && p[i] > 0)) {
	nodes++;	
	p0 = p[i];
      }
    }
  }

  i1 += 5;
  if (e > 0.0) {
    *k2 = i2-2;
    if (IsOdd(*k2)) (*k2)--;
  } else {
    if (p[i2-1] > p[i2]) {
      for (i = i2-1; i > i1; i--) {
	if (p[i-1] <= p[i]) break;
      }
    } else {
      for (i = i2-1; i > i1; i--) {
	if (p[i-1] >= p[i]) break;
      }
    }
    if (i == i2-1) *k2 = i-1;
    else *k2 = i;
  }
 
  return nodes;
}

int _Inward(double *p, double e, POTENTIAL *pot, int i2) {
  int i;
  double a, b, r, x, y, z;

  for (i = i2; i < MAX_POINTS; i++) {
    r = pot->rad[i];
    x = 2.0*(_veff[i] - e);
    x *= 4.0*r*r;
    a = pot->ar;
    b = pot->br;
    z = sqrt(r);
    y = a*z + 2.0*b;
    y = 1/(y*y);
    z = (0.75*a*a*r + 5.0*a*b*z +4.0*b*b);
    x += z*y;    
    x *= y / 12.0;
    _ACOEF[i] = 1.0 - x;
    _BCOEF[i] = 2.0 + 10.0*x;
  }
  
  i = MAX_POINTS - 1;
  p[i] = 0.0;
  _VCOEF[i] = 0.0;
  i--;
  for (; i > i2; i--) {
    a = _BCOEF[i] - _ACOEF[i+1]*_VCOEF[i+1];
    _VCOEF[i] = _ACOEF[i-1] / a;
  }
  
  for (i = i2+1; i < MAX_POINTS-1; i++) {
    p[i] = _VCOEF[i]*p[i-1];
  }
  
  return 0;
}

double InnerProduct(int n, double *p1, double *p2, POTENTIAL *pot) {
  int i;

  for (i = 0; i < n; i++) {
    _dwork[i] = p1[i]*p2[i] * pot->dr_drho[i];
  }
  _dwork4[0] = 0.0;
  NewtonCotes(_dwork4, _dwork, 0, n-1, 0);
  return _dwork4[n-2];
}

double Simpson(double *y, int ia, int ib) {
  int i;
  double a;

  a = 0.0;
  
  for (i = ia; i < ib - 1; i += 2) {
    a += (y[i] + 4.0*y[i+1] + y[i+2])/3.0;
  }
  if (i < ib) a += 0.5 * (y[i] + y[ib]);

  return a;
}

/* integration by newton-cotes formula */
int NewtonCotes(double *r, double *x, int i0, int i1, int m) {
  int i, j, n, k;
  double a;

  for (i = i0; i <= i1-4; i += 4) {
    a = 0.0;
    for (j = 0, k = i; j <= 4; j++, k++) {
      a += _CNC[4][j] * x[k];
    }
    r[i+4] = r[i] + a;
  }
  
  if (i1 < MAX_POINTS-1) {
    if (i > i0) {
      k = i - 3;
      n = i1 - i + 5;
    } else {
      k = i + 1;
      n = i1 - i + 1;
    }
    a = 0.0;
    for (j = 1; j < n; j++, k++) {
      a += _ONC[n][j] * x[k];
    }
    r[i1+1] = r[i1+1-n] + a;
  }

  if (m >= 0) {
    n = i1 - i - 1;
    if (n > 0) {
      a = 0.0;
      for (j = 0, k = i; j <= n; j++, k++) {
	a += _CNC[n][j] * x[k];
      }
      r[i1-1] = r[i] + a;
    }
    n++;
    a = 0.0;
    for (j = 0, k = i; j <= n; j++, k++) {
      a += _CNC[n][j] * x[k];
    }
    r[i1] = r[i] + a;
  } else {
    for (i = i0; i <= i1; i += 4) {
      for (n = 1; n <= Min(3, i1-i); n++) {
	a = 0.0;
	for (j = 0, k = i; j <= n; j++, k++) {
	  a += _CNC[n][j] * x[k];
	}
	r[i+n] = r[i] + a;
      }
    }
  }

  return 0;
}

int SetOrbitalRGrid(POTENTIAL *pot, double rmin, double rmax) {
  int i;  
  double z0, z, d1, d2, del;
  double a, b;

  z0 = GetAtomicNumber();
  z = z0;
  if (pot->N > 0) z = (z - pot->N + 1);
  if (pot->flag == 0) pot->flag = -1; 

  if (rmin <= 0.0) rmin = 1E-5;
  if (rmax <= 0.0) rmax = 5E+2;
  nmax = sqrt(rmax)/2.0;

  rmin /= z0;
  rmax /= z;
  
  d1 = log(rmax/rmin);
  d2 = sqrt(rmax) - sqrt(rmin);

  a = 16.0*sqrt(2.0*z)/PI;
  b = (MAX_POINTS - 1.0 - (a*d2))/d1;
  if (b < 1.0/log(1.1)) {
    printf("Not enough radial mesh points, ");
    printf("enlarge to at least %d\n", (int) (1 + a*d2 + d1/log(1.1)));
    exit(1);
  }

  d1 = b*d1;
  d2 = a*d2;
  del = (d1 + d2)/(MAX_POINTS - 1);
  pot->rad[0] = rmin;
  d1 = a*sqrt(rmin) + b*log(rmin);
  for (i = 1; i < MAX_POINTS; i++) {
    d1 += del;
    pot->rad[i] = GetRFromRho(d1, a, b, pot->rad[i-1]);
  }

  pot->ar = a;
  pot->br = b;

  for (i = 0; i < MAX_POINTS; i++) {
    d1 = a * sqrt(pot->rad[i]);
    d2 = 2.0*pot->rad[i];
    pot->dr_drho[i] = d2/(d1 + 2.0*b);
    pot->drho_dr2[i] = (d1 + 4.0*b) / (d2*d2);
  }

  return 0;
}

double GetRFromRho(double rho, double a, double b, double r0) {
  double e, d1;
  int i;

  e = 1.0;
  i = 0;
  while (fabs(e) > 1E-10) {
    if (i > 100) {
      printf("Newton iteration failed to converge in GetRFromRho\n");
      exit(1);
    }
    d1 = sqrt(r0)*a;
    e = d1 + b*log(r0) - rho;
    e /= (0.5*d1 + b);
    r0 *= (1.0 - e);
    i++;
  }

  return r0;
}

int SetPotentialZ(POTENTIAL *pot, double c) {
  int i;

  c = 1.0+c;
  for (i = 0; i < MAX_POINTS; i++) {
    pot->Z[i] = c*GetAtomicEffectiveZ(pot->rad[i]);
  }
  return 0;
}

int SetPotentialVc(POTENTIAL *pot) {
  int i;
  double n, r, r2, v, x, y, a, b, v0;

  for (i = 0; i < MAX_POINTS; i++) {
    r = pot->rad[i];
    pot->Vc[i] = - (pot->Z[i] / r);
    r2 = r*r;
    pot->dVc[i] = pot->Z[i] / r2; 
    pot->dVc2[i] = -2.0 * pot->Z[i] / (r2*r);
  }

  n = pot->N - 1;
  if (n > 0 && (pot->a > 0 || pot->lambda > 0)) {
    for (i = 0; i < MAX_POINTS; i++) {
      r = pot->rad[i];
      r2 = r*r;
      v0 = n/r;
      x = 1.0 + r*pot->a;
      a = exp(-pot->lambda * r);
      v = v0 * (1.0 - a/x);
      pot->Vc[i] += v;      
      b = (pot->lambda + pot->a/x);
      y = -v/r + (v0 - v)*b;
      pot->dVc[i] += y;
      pot->dVc2[i] -= y/r;
      pot->dVc2[i] += v/r2;
      pot->dVc2[i] -= (v0/r + y) * b;
      pot->dVc2[i] -= (v0 - v)*(pot->a*pot->a)/(x*x);
    }
  }
  return 0;
}

int SetPotentialU(POTENTIAL *pot, int n, double *u) {
  int i;
  
  if (n < 0) {
    for (i = 0; i < MAX_POINTS; i++) { 
      pot->U[i] = 0.0;
      pot->dU[i] = 0.0;
      pot->dU2[i] = 0.0;
    }
    return 0;
  }

  for (i = 0; i < n; i++) {
    pot->U[i] = u[i];    
  }

  for (i = 1; i < MAX_POINTS-1; i++) {
    pot->dU[i] = pot->U[i+1] - pot->U[i-1];
    pot->dU[i] *= 0.5;
    pot->dU[i] /= pot->dr_drho[i];
  }
  pot->dU[0] = pot->dU[1];
  pot->dU[MAX_POINTS-1] = pot->dU[MAX_POINTS-2];
  
  for (i = 1; i < MAX_POINTS-1; i++) {
    pot->dU2[i] = pot->dU[i+1] - pot->dU[i-1];
    pot->dU2[i] *= 0.5;
    pot->dU2[i] /= pot->dr_drho[i];
  }
  pot->dU2[0] = pot->dU2[1];
  pot->dU2[MAX_POINTS-1] = pot->dU2[MAX_POINTS-2];
  
  return 0;
}

int SetPotentialW (POTENTIAL *pot, double e, int kappa) {
  int i;
  double xi, r, x, y, z;

  for (i = 0; i < MAX_POINTS; i++) {
    xi = e - pot->Vc[i] - pot->U[i];
    r = xi*FINE_STRUCTURE_CONST2*0.5 + 1.0;
  
    x = pot->dU[i] + pot->dVc[i];
    y = - 2.0*kappa*x/pot->rad[i];
    x = x*x*0.75*FINE_STRUCTURE_CONST2/r;
    z = (pot->dU2[i] + pot->dVc2[i]);
    pot->W[i] = x + y + z;
    pot->W[i] /= 4.0*r;
    x = xi*xi;
    pot->W[i] = x - pot->W[i];
    pot->W[i] *= 0.5*FINE_STRUCTURE_CONST2;
    pot->W[i] = -pot->W[i];
  }

  for (i = 1; i < MAX_POINTS-1; i++) {
    pot->dW[i] = pot->W[i+1] - pot->W[i-1];   
    pot->dW[i] *= 0.5;
    pot->dW[i] /= pot->dr_drho[i];
  }
  pot->dW[0] = pot->dW[1];
  pot->dW[MAX_POINTS-1] = pot->dW[MAX_POINTS-2];
  
  for (i = 1; i < MAX_POINTS-1; i++) {
    pot->dW2[i] = pot->dW[i+1] - pot->dW[i-1];
    pot->dW2[i] *= 0.5;
    pot->dW2[i] /= pot->dr_drho[i];
  }
  pot->dW2[0] = pot->dW2[1];
  pot->dW2[MAX_POINTS-1] = pot->dW2[MAX_POINTS-2];
  
  return 0;
}
