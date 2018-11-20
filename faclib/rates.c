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

#include "rates.h"
#include "interpolation.h"
#include "cf77.h"
#include "parser.h"

static char *rcsid="$Id$";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

static int iedist = 0;
static DISTRIBUTION ele_dist[MAX_DIST];
static int ipdist = 0;
static DISTRIBUTION pho_dist[MAX_DIST];
static int ixdist = 0;
static DISTRIBUTION cxt_dist[MAX_DIST];

#define QUAD_LIMIT 64
static int _iwork[QUAD_LIMIT];
static double _dwork[4*QUAD_LIMIT];

#define N3BRI 2000
static double gamma3b = 1.0;

static double rate_epsabs = EPS8;
static double rate_epsrel = EPS3;
static int rate_iprint = 1;

static KRONOS _kronos_cx[3];
static double _cxldistmj = 1.0;

static struct {
  DISTRIBUTION *d;
  double (*Rate1E)(double, double, int, void *);
  double eth;
  int np;
  void *params;
  int i, f;
  int type;
  int xlog, elog;
  double eg[N3BRI], fg[N3BRI];
} rate_args;

#pragma omp threadprivate(rate_args, _iwork, _dwork)

#define NSEATON 19
static double log_xseaton[NSEATON];
static double xseaton[NSEATON] = {
  0.02, 0.03, 0.045, 0.07, 0.11, 0.15, 0.20, 0.30, 0.45, 
  0.70, 1.1, 1.5, 2.0, 3.0, 4.5, 7.0, 11.0, 15.0, 20.0};
static double seaton[6][NSEATON] = {
  {0.06845, 0.09148, 0.1209, 0.1615, 0.2133, 0.2552, 
   0.2987, 0.3668, 0.4413, 0.5269, 0.6146, 0.6724, 
   0.7226, 0.7862, 0.8399, 0.8865, 0.9223, 0.9408, 0.9544},
  {0.00417, 0.00538, 0.0067, 0.0080, 0.0088, 0.0086,
   0.0074, 0.0035, -0.0044, -0.0190, -0.0422, -0.0638,
   -0.0882, -0.1295, -0.1784, -0.2397, -0.3099, -0.3622, -0.4147},
  {-0.00120, -0.00175, -0.0025, -0.0038, -0.0056, -0.0074,
   -0.0095, -0.0133, -0.0187, -0.0269, -0.0388, -0.0497, 
   -0.0625, -0.0859, -0.1174, -0.1638, -0.228, -0.285, -0.350},
  {0.0877, 0.1199, 0.1627, 0.2247, 0.3090, 0.3815,
   0.4611, 0.5958, 0.7595, 0.9733, 1.231, 1.431, 
   1.632, 1.938, 2.268, 2.650, 3.059, 3.348, 3.621},
  {0.0053, 0.0073, 0.0097, 0.0130, 0.0168, 0.0195,
   0.0218, 0.0242, 0.0241, 0.0193, 0.006, -0.011,
   -0.032, -0.076, -0.138, -0.230, -0.354, -0.459, -0.570},
  {-0.0012, -0.0018, -0.0027, -0.0041, -0.0062, -0.0082, 
   -0.0106, -0.0152, -0.0216, -0.0316, -0.046, -0.060,
   -0.076, -0.106, -0.147, -0.208, -0.296, -0.375, -0.467}
};

#define NGAUNT0 11
static double bgaunt0[NGAUNT0][NGAUNT0] = {
  {2.31448E+00,-2.39863E+00, 9.79187E-01, 4.53690E-01,-1.69366E-01,
   -4.30246E-01,-3.84721E-01, 1.07579E-01, 3.78713E-01, 2.73831E-02,-1.06087E-01},
  {-5.82308E-01, 9.45299E-01, 7.50594E-01, 8.33913E-02, 3.69950E-01, 
   8.54940E-01,-2.53707E-01,-1.46431E+00,-6.03851E-01, 5.60237E-01, 3.61029E-01},
  {-7.69314E-01, 1.90311E-01,-3.88671E-01,-1.26175E+00,-1.79261E+00, 
   2.11206E+00, 6.98574E+00, 5.64434E-01,-6.05501E+00,-1.02473E+00, 1.64915E+00},
  {2.09120E-01,-6.34270E-02,-2.13061E-01,-2.26219E-01,-2.86365E+00,
   -7.69346E+00,-2.28229E+00, 1.13788E+01, 8.47077E+00,-4.33454E+00,-4.10317E+00},
  {8.47001E-01, 1.82172E-01, 3.96012E-01, 4.08453E+00, 1.01333E+01,
   -6.40153E+00,-3.14652E+01,-3.86629E+00, 2.72361E+01, 4.91736E+00,-7.45031E+00},
  {-1.21671E-01,-1.87775E-01,-5.42449E-01,-1.32846E+00, 4.46956E+00, 
   2.03609E+01, 1.11697E+01,-2.88967E+01,-2.60411E+01, 1.11268E+01, 1.20846E+01},
  {-7.67073E-01,-4.70394E-02, 9.76439E-02,-7.60194E+00,-2.17329E+01, 
   1.11517E+01, 6.35508E+01, 7.89505E+00,-5.56577E+01,-9.73728E+00, 1.54684E+01},
  {3.01370E-02, 2.44838E-01, 7.36601E-01, 2.14229E+00,-3.24693E+00,
   -2.18736E+01,-1.45136E+01, 3.07839E+01, 3.00594E+01,-1.19819E+01,-1.37950E+01},
  {4.40477E-01,-1.41621E-01,-5.24374E-01, 6.73369E+00, 2.06793E+01,
   -9.63661E+00,-5.90074E+01,-7.06822E+00, 5.22275E+01, 8.78611E+00,-1.47207E+01},
  {5.12129E-03,-9.59773E-02,-2.95736E-01,-9.46625E-01, 9.39673E-01,
   8.37046E+00, 6.04020E+00,-1.17635E+01,-1.19353E+01, 4.61643E+00, 5.46132E+00},
  {-1.12455E-01, 8.39735E-02, 2.64128E-01,-2.26435E+00,-7.24494E+00, 
   3.22187E+00, 2.04415E+01, 2.34034E+00,-1.82417E+01,-2.96404E+00, 5.19706E+00}
};

#define NGAUNT1 28
static double bgaunt1[2][NGAUNT1] = {
  {-6.00000E+00,-5.70000E+00,-5.40000E+00,-5.10000E+00,-4.80000E+00,
   -4.50000E+00,-4.20000E+00,-3.90000E+00,-3.60000E+00,-3.30000E+00,
   -3.00000E+00,-2.70000E+00,-2.40000E+00,-2.10000E+00,-1.80000E+00,
   -1.50000E+00,-1.20000E+00,-9.00000E-01,-6.00000E-01,-3.00000E-01,
   -2.22045E-16, 3.00000E-01, 6.00000E-01, 9.00000E-01, 1.20000E+00, 
   1.50000E+00, 1.80000E+00, 2.10000E+00},
  {1.98111E+00, 7.41511E-01, 1.79845E-01, 3.86527E-01, 6.29494E-01, 
   1.01035E+00, 1.12238E+00, 1.14504E+00, 1.16063E+00, 1.18338E+00, 
   1.21374E+00, 1.25252E+00, 1.29966E+00, 1.35087E+00, 1.39724E+00, 
   1.42833E+00, 1.43655E+00, 1.42038E+00, 1.38486E+00, 1.33893E+00, 
   1.29126E+00, 1.24688E+00, 1.20685E+00, 1.17155E+00, 1.14376E+00, 
   1.12047E+00, 1.05186E+00, 9.22784E-01}
};

KRONOS *KronosCX(int k) {
  if (k < 0 || k > 2) return NULL;
  return &_kronos_cx[k];
}

DISTRIBUTION *GetEleDist(int *i) {
  if (i) *i = iedist;
  return ele_dist+iedist;
}

DISTRIBUTION *GetPhoDist(int *i) {
  if (i) *i = ipdist;
  return pho_dist+ipdist;
}

DISTRIBUTION *GetCxtDist(int *i) {
  if (i) *i = ixdist;
  return cxt_dist+ixdist;
}

int SetRateAccuracy(double epsrel, double epsabs) {
  if (epsrel > 0.0) rate_epsrel = epsrel;
  if (epsabs > 0.0) rate_epsabs = epsabs;
  return 0.0;
}    

void SetGamma3B(double g) {
  gamma3b = g;
}

static void ThreeBodyDist(void) {
  int i, j, n;
  DISTRIBUTION *d;
  double emin, emax, de, y[N3BRI], x, t, *eg, c;
  
  d = ele_dist + iedist;
  if (iedist == MAX_DIST-1) {
    n = d->params[0];
    i = d->params[1];
    eg = &(d->params[3]);
    emin = eg[0];
    emax = eg[n-1];
    if (i == 1) {
      emin = exp(emin);
      emax = exp(emax);
    }
  } else {
    n = d->nparams;
    emin = d->params[n-2];
    emax = d->params[n-1];
  }
  if (emax/emin > 10.0) {
    emin = log(emin);
    emax = log(0.5*emax);
    rate_args.elog = 1;
  } else {
    emax = 0.5*emax;
    rate_args.elog = 0;
  }
  
  de = (emax - emin)/(N3BRI-1);
  rate_args.eg[0] = emin;
  rate_args.fg[0] = 0.0;
  for (i = 1; i < N3BRI; i++) {
    rate_args.eg[i] = rate_args.eg[i-1] + de;
  }
  y[0] = 0.0;
  
  c = DLOGAM(2.0*(gamma3b+1.0)) - 2.0*DLOGAM(gamma3b+1.0);
  c = exp(c);
  for (i = 1; i < N3BRI; i++) {
    x = rate_args.eg[i];
    if (rate_args.elog) x = exp(x);
    for (j = 1; j <= i; j++) {
      t = rate_args.eg[j];
      if (rate_args.elog) t = exp(t);
      y[j] = d->dist(t, d->params)/sqrt(t);
      if (rate_args.elog) y[j] *= t;
      t = 2*x - t;
      y[j] *= d->dist(t, d->params)/sqrt(t);
      t /= 2*x;
      y[j] *= c*pow(t*(1.0-t), gamma3b);
    }
    rate_args.fg[i] = Simpson(y, 0, i)*de;
  }  
}

static double RateIntegrand(double *e) {
  double a, b, x;
  double p = 1.46366E-12; /* (h^2/2m)^1.5/(4*pi) cm^3*eV^1.5 */

  if (rate_args.xlog) {
    x = exp(*e);
  } else {
    x = *e;
  }

  if (rate_args.type != -RT_CI) {
    a = rate_args.d->dist(x, rate_args.d->params);
  } else {
    if (x > rate_args.eth) {
      b = 0.5*(x - rate_args.eth);
      if (rate_args.elog) b = log(b);
      if (b < rate_args.eg[0]) a = rate_args.fg[0];
      else if (b > rate_args.eg[N3BRI-1]) a = rate_args.fg[N3BRI-1];
      else {
	UVIP3P(3, N3BRI, rate_args.eg, rate_args.fg, 1, &b, &a);
	a *= 2.0*p*sqrt(x)/(x-rate_args.eth);
	a *= VelocityFromE(x, -1.0)/VelocityFromE(x, 1.0);
      }
    } else {
      a = 0.0;
    }
  }
  b = rate_args.Rate1E(x, rate_args.eth, rate_args.np, rate_args.params);
  if (rate_args.xlog) {
    x = x*a*b;
  } else {
    x = a*b;
  }
  return x;
}

/* provide fortran access with cfortran.h */
FCALLSCFUN1(DOUBLE, RateIntegrand, RATEINTEGRAND, rateintegrand, PDOUBLE)
	    
double IntegrateRate(int idist, double eth, double bound, 
		     int np, void *params, int i0, int f0, int type, 
		     double (*Rate1E)(double, double, int, void *)) { 
  double result;
  int neval, ier, limit, lenw, last, n, ix, iy;
  double epsabs, epsrel, abserr;
  double a, b, a0, b0, r0, *eg;

  ier = 0;
  epsabs = rate_epsabs;
  epsrel = rate_epsrel;

  limit = QUAD_LIMIT;
  lenw = 4*limit;
  
  rate_args.Rate1E = Rate1E;
  if (idist == 0) rate_args.d = ele_dist + iedist;
  else if (idist == 1) rate_args.d = pho_dist + ipdist;
  else rate_args.d = cxt_dist + ixdist;
  rate_args.eth = eth;
  rate_args.np = np;
  rate_args.params = params;
  rate_args.i = i0;
  rate_args.f = f0;
  rate_args.type = type;
  rate_args.xlog = rate_args.d->xlog;
  
  if ((idist == 0 && iedist == MAX_DIST-1) ||
      (idist == 1 && ipdist == MAX_DIST-1)) {
    n = rate_args.d->params[0];
    ix = rate_args.d->params[1];
    iy = rate_args.d->params[2];
    eg = &(rate_args.d->params[3]);
    a = eg[0];
    b = eg[n-1];
    rate_args.xlog = ix;
  } else {
    n = rate_args.d->nparams;  
    b = rate_args.d->params[n-1];
    a = rate_args.d->params[n-2];    
    if (rate_args.xlog < 0) {
      if (b/a > 10) {
	rate_args.xlog = 1;
	a = log(a);
	b = log(b);
      } else {
	rate_args.xlog = 0;
      }
    }
  }
  
  if (rate_args.xlog) {
    bound = log(bound);
  }
  if (bound > a) a = bound;
  if (b <= a) return 0.0;
  if (idist == 0 && iedist == 0 && rate_args.xlog == 0) {
    a0 = rate_args.d->params[0];
    b0 = 5.0*a0;
    a0 = a;
    if (b < b0) b0 = b;
    r0 = 0.0;
    if (b0 > a0) {
      DQAGS(C_FUNCTION(RATEINTEGRAND, rateintegrand), 
	    a0, b0, epsabs, epsrel, &result, 
	    &abserr, &neval, &ier, limit, lenw, &last, _iwork, _dwork);
      r0 += result;
      if (abserr > epsabs && abserr > r0*epsrel) {
	if (ier != 0 && rate_iprint) {
	  MPrintf(-1,
		  "IntegrateRate Error0: %d %d %10.3E %10.3E %10.3E %10.3E\n", 
		  ier, neval, a0, b0, result, abserr);
	  MPrintf(-1, "%6d %6d %2d Eth = %10.3E\n", 
		  rate_args.i, rate_args.f, type, eth);
	  Abort(1);
	}
      }
      result = 0.1*r0*epsrel;
      if (epsabs < result) epsabs = result;
    } else {
      b0 = a0;
    }
    if (b > b0) {
      DQAGS(C_FUNCTION(RATEINTEGRAND, rateintegrand), 
	    b0, b, epsabs, epsrel, &result, 
	    &abserr, &neval, &ier, limit, lenw, &last, _iwork, _dwork);
      r0 += result;
      if (abserr > epsabs && abserr > r0*epsrel) {
	if (ier != 0 && rate_iprint) {
	  MPrintf(-1,
		  "IntegrateRate Error1: %d %d %10.3E %10.3E %10.3E %10.3E\n", 
		  ier, neval, b0, b, result, abserr);
	  MPrintf(-1, "%6d %6d %2d Eth = %10.3E\n", 
		  rate_args.i, rate_args.f, type, eth);
	  Abort(1);
	}
      }
    }
    if (r0 < 0.0) r0 = 0.0;
    return r0;
  } else {
    DQAGS(C_FUNCTION(RATEINTEGRAND, rateintegrand), 
	  a, b, epsabs, epsrel, &result, 
	  &abserr, &neval, &ier, limit, lenw, &last, _iwork, _dwork);
    r0 = result;
    if (abserr > epsabs && abserr > r0*epsrel) {
      if (ier != 0 && rate_iprint) {
	MPrintf(-1,
		"IntegrateRate Error2: %d %d %10.3E %10.3E %10.3E %10.3E\n", 
		ier, neval, a, b, result, abserr);
	MPrintf(-1, "%6d %6d %2d Eth = %10.3E\n", 
		rate_args.i, rate_args.f, type, eth);
	Abort(1);
      }
    }
    if (r0 < 0.0) r0 = 0.0;

    return r0;
  }
}

double IntegrateRate2(int idist, double e, int np, 
		      void *params, int i0, int f0, int type,
		      double (*Rate1E)(double, double, double, int, void *)) { 
  double a;
  
  a = Rate1E(e, e, e, np, params);
  return a;
}

double VelocityFromE(double e, double bms) {
  double k;

  if (bms < 0) {
    k = sqrt(e/fabs(bms))*5.92928E-3;
  } else if (bms == 1.0 || bms == 0.0) {
    k = e/HARTREE_EV;
    k = 2.0*k*(1.0 + 0.5*FINE_STRUCTURE_CONST2*k);
    k = FINE_STRUCTURE_CONST2*k;
    k = sqrt(k/(1.0+k));
    k /= FINE_STRUCTURE_CONST;
    k *= RBOHR*RATE_AU12*1E-6; /* in unit of 10^10 cm/s */
  } else {
    k = e/HARTREE_EV;
    k = 2.0*bms*k*(1.0 + 0.5*FINE_STRUCTURE_CONST2*k/bms);
    k = FINE_STRUCTURE_CONST2*k;
    k = sqrt(k/(bms*bms+k));
    k /= FINE_STRUCTURE_CONST;
    k *= RBOHR*RATE_AU12*1E-6; /* in unit of 10^10 cm/s */
  }
  
  return k;
}

double CXRate1E(double e1, double eth0, int np, void *p) {  
  KRONOS *cx;
  int *ip = (int *) p;
  cx = &_kronos_cx[ip[0]];
  double *x, *y, r;
  int n = 3;
  int one = 1;
  if (ip[3] == 0) {
    e1 /= cx->pmass/AMU;
  }
  double e = e1;
  if (cx->ilog & 1) {
    e = log(e);
  }
  x = cx->ep;
  if (ip[0] == 2) {
    y = cx->rcx;
  } else {
    if (ip[1] == 0) y = cx->cx0[ip[2]];
    else y = cx->cx1[ip[2]];
  }
  const double vth = 0.547;
  double v = VelocityFromE(e1, AMU);
  if (e < cx->ep[0]) e = cx->ep[0];
  if (e > cx->ep[cx->nep-1]) {
    if (v >= vth) return 0.0;
    r = y[cx->nep-1];
    double vx = cx->ep[cx->nep-1];
    if (cx->ilog & 1) {
      vx = exp(vx);
    }
    vx = VelocityFromE(vx, AMU);
    vx = vth/vx;
    if (cx->ilog & 2) {
      if (r < -300) r = 0.0;
      else {
	r = exp(r)/log(vx);
      }
    }
    r *= log(vth/v);
  } else {    
    UVIP3P(n, cx->nep, x, y, one, &e, &r);
    if (cx->ilog & 2) {
      if (r < -300) r = 0.0;
      else r = exp(r);
    }
  }
  if (ip[0] == 2) {
    r *= AREA_AU20;
  } else {
    r *= 1e4;
  }
  r *= v;
  return r;
}

int CXRate(double *dir, int *ip, int i0, int f0) {
  KRONOS *cx = &_kronos_cx[ip[0]];
  double e0 = cx->ep[0];
  if (cx->ilog & 1) e0 = exp(e0);
  *dir = IntegrateRate(2, e0, 0, 4, ip, i0, f0, RT_CX, CXRate1E);
  return 0;
}

double CERate1E(double e1, double eth0, int np, void *p) {
  double *x, *y;
  int m1, n, one;
  double *dp, a, x0, y0;
  double e0, d, c, b, b0, b1;
  double bte, bms, eth, e;

  BornFormFactorTE(&bte);
  bms = BornMass();  
  eth = (eth0 + bte*HARTREE_EV)/bms;
  e = e1/bms;

  if (e < eth) return 0.0;

  dp = (double *) p;
  m1 = np + 1;
  y = dp+2;
  x0 = log((dp[0]+e-eth)/dp[0]);
  x = y + m1;

  if (x0 <= x[np-1]) {
    n = 2;
    one = 1;
    if (fabs(bms-1.0) < EPS3 || x0 >= x[0]) {
      UVIP3P(n, np, x, y, one, &x0, &a);
      if (a < 0.0) a = 0.0;
    } else {
      a = y[0] * pow(exp(x0-x[0]), 2.5);
    }
  } else {
    x0 = (e-eth)/(dp[0]+e-eth);
    y0 = y[np-1];
    if (dp[1] > 0) {
      e0 = (x[np]*dp[0]/(1.0-x[np]) + eth)/HARTREE_EV;
      b0 = 1.0 + FINE_STRUCTURE_CONST2*e0;
      b1 = 1.0 + FINE_STRUCTURE_CONST2*(e0-eth/HARTREE_EV);
      y0 /= b0*b1;
      d = 2.0*e0*(1.0+0.5*FINE_STRUCTURE_CONST2*e0);
      c = FINE_STRUCTURE_CONST2*d;
      b = log(0.5*d*HARTREE_EV/eth0) - c/(1.0+c);
      y0 -= dp[1]*b;
      a = y[np] + (x0-1.0)*(y0-y[np])/(x[np]-1.0);
      e0 = e/HARTREE_EV;
      d = 2.0*e0*(1.0+0.5*FINE_STRUCTURE_CONST2*e0);
      c = FINE_STRUCTURE_CONST2*d;
      b = log(0.5*d*HARTREE_EV/eth0) - c/(1.0+c);  
      a += dp[1]*b;
      b0 = 1.0 + FINE_STRUCTURE_CONST2*e0;
      b1 = 1.0 + FINE_STRUCTURE_CONST2*(e0-eth/HARTREE_EV);
      a *= b0*b1;
    } else if (dp[1] + 1.0 == 1.0) {
      e0 = (x[np]*dp[0]/(1.0-x[np]) + eth)/HARTREE_EV;
      b0 = 1.0 + FINE_STRUCTURE_CONST2*e0;
      b1 = 1.0 + FINE_STRUCTURE_CONST2*(e0-eth/HARTREE_EV);
      y0 /= b0*b1;
      d = 2.0*e0*(1.0+0.5*FINE_STRUCTURE_CONST2*e0);
      c = FINE_STRUCTURE_CONST2*d;
      a = y[np] + (x0-1.0)*(y0-y[np])/(x[np]-1.0);
      e0 = e/HARTREE_EV;
      b0 = 1.0 + FINE_STRUCTURE_CONST2*e0;
      b1 = 1.0 + FINE_STRUCTURE_CONST2*(e0-eth/HARTREE_EV);
      a *= b0*b1;
    } else {
      a = y[np] + (x0-1.0)*(y0-y[np])/(x[np]-1.0);
    }
  }

  if (a <= 0.0) {
    a = 0.0;
    return a;
  }
  
  e0 = e/HARTREE_EV;
  b = 2.0*e0*(1.0+0.5*FINE_STRUCTURE_CONST2*e0);  
  a *= PI*AREA_AU20/b;
  a *= VelocityFromE(e1, bms);
  return a;
}

double DERate1E(double e1, double eth0, int np, void *p) {
  double a, x0, y0, *x, *y;
  double *dp;
  int m1, n, one;
  double e0, d, c, b, b0, b1;
  double bte, bms, eth, e;

  BornFormFactorTE(&bte);
  bms = BornMass();
  eth = (eth0 + bte*HARTREE_EV)/bms;
  e = e1/bms;

  dp = (double *) p;
  m1 = np + 1;
  x0 = log((dp[0]+e)/dp[0]);
  y = dp+2;
  x = y + m1;

  if (x0 <= x[np-1]) {
    n = 2;
    one = 1;
    if (fabs(bms-1.0) < EPS3 || x0 >= x[0]) {
      UVIP3P(n, np, x, y, one, &x0, &a);
      if (a < 0.0) a = 0.0;
    } else {
      a = y[0] * pow(exp(x0-x[0]), 2.5);
    }
  } else {
    x0 = e/(dp[0]+e);
    y0 = y[np-1]; 
    if (dp[1] > 0) {
      e0 = (x[np]*dp[0]/(1.0-x[np]) + eth)/HARTREE_EV;
      b0 = 1.0 + FINE_STRUCTURE_CONST2*e0;
      b1 = 1.0 + FINE_STRUCTURE_CONST2*(e0-eth/HARTREE_EV);
      y0 /= b0*b1;
      d = 2.0*e0*(1.0+0.5*FINE_STRUCTURE_CONST2*e0);
      c = FINE_STRUCTURE_CONST2*d;
      b = log(0.5*d*HARTREE_EV/eth0) - c/(1.0+c);  
      y0 -= dp[1]*b;
      a = y[np] + (x0-1.0)*(y0-y[np])/(x[np]-1.0);
      e0 = (e + eth)/HARTREE_EV;
      d = 2.0*e0*(1.0+0.5*FINE_STRUCTURE_CONST2*e0);
      c = FINE_STRUCTURE_CONST2*d;
      b = log(0.5*d*HARTREE_EV/eth0) - c/(1.0+c);  
      a += dp[1]*b;
      b0 = 1.0 + FINE_STRUCTURE_CONST2*e0;
      b1 = 1.0 + FINE_STRUCTURE_CONST2*(e0-eth/HARTREE_EV);
      a *= b0*b1;
    } else if (dp[1] + 1.0 == 1.0) {
      e0 = (x[np]*dp[0]/(1.0-x[np]) + eth)/HARTREE_EV;
      b0 = 1.0 + FINE_STRUCTURE_CONST2*e0;
      b1 = 1.0 + FINE_STRUCTURE_CONST2*(e0-eth/HARTREE_EV);
      y0 /= b0*b1;
      a = y[np] + (x0-1.0)*(y0-y[np])/(x[np]-1.0);
      e0 = (e + eth)/HARTREE_EV;
      b0 = 1.0 + FINE_STRUCTURE_CONST2*e0;
      b1 = 1.0 + FINE_STRUCTURE_CONST2*(e0-eth/HARTREE_EV);
      a *= b0*b1;
    } else {
      a = y[np] + (x0-1.0)*(y0-y[np])/(x[np]-1.0);
    }
  }

  if (a <= 0.0) {
    a = 0.0;
    return a;
  }

  e0 = e/HARTREE_EV;
  b = 2.0*e0*(1.0+0.5*FINE_STRUCTURE_CONST2*e0);
  a *= PI*AREA_AU20/b;
  a *= VelocityFromE(e1, bms);
  return a;
}
  
int CERate(double *dir, double *inv, int iinv, 
	   int j1, int j2, double e,
	   int m, double *params, int i0, int f0) {
  double a, e0;

  e0 = e*HARTREE_EV;	     
  a = IntegrateRate(0, e0, e0/BornMass(), m, params, 
		    i0, f0, RT_CE, CERate1E);
  *dir = a/(j1 + 1.0);
  if (iinv) {
    if (iedist == 0) {
      a *= exp(e0/ele_dist[0].params[0]);
      *inv = a/(j2 + 1.0);
    } else {
      *inv = IntegrateRate(0, e0, 0.0, m, params, 
			   i0, f0, -RT_CE, DERate1E);
      *inv /= (j2 + 1.0);
    }
  } else {
    *inv = 0.0;
  }
  return 0;
}

int TRRate(double *dir, double *inv, int iinv,
	   int j1, int j2, double e, float strength) {
  double a, b, e0;
  const double factor = 4.73313707E-2;

  a = FINE_STRUCTURE_CONST * e;
  a = 2.0*a*a*FINE_STRUCTURE_CONST;
  a *= strength;
  a *= RATE_AU;
  *dir = a/(j1 + 1.0);
  if ((*dir) < 0) printf("%10.3E %10.3E %10.3E %10.3E\n", a, e, strength, j1+1.0);
  if (iinv) {
    e0 = e*HARTREE_EV;
    b = pho_dist[ipdist].dist(e0, pho_dist[ipdist].params);
    if (b > 0) {
      *inv = factor*b*a/(e0*e0*e0*(j2+1.0));
    } else {
      *inv = 0.0;
    }
  } else {
    *inv = 0.0;
  }
  return 0;
}

double CIRate1E(double e, double eth, int np, void *p) {
  float *dp;
  double x;
  double a, b, c, f;
  double logx;

  if (e < eth) return 0.0;
  dp = (float *) p;
  x = e/eth;
  logx = log(x);
  
  a = 1.0/x;
  b = 1.0 - a;
  c = b * a;
  
  f = dp[0]*logx;
  f += dp[1]*b*b;
  f += dp[2]*c;
  f += dp[3]*a*c;
  
  c = AREA_AU20*HARTREE_EV*f/(2.0*e);
  c *= VelocityFromE(e, 0.0);

  return c;
}

int CIRate(double *dir, double *inv, int iinv, 
	   int j1, int j2, double e,
	   int m, float *params, int i0, int f0) {
  const double p = 1.65156E-12; /* 0.5*(h^2/(2pi*m*eV))^{3/2} */
  double e0, a;
  
  e0 = e*HARTREE_EV;
  a = IntegrateRate(0, e0, e0, m, params, 
		    i0, f0, RT_CI, CIRate1E);
  *dir = a/(j1 + 1.0);
  if (iinv) {
    if (iedist == 0) {
      a *= exp(e0/ele_dist[0].params[0]);
      *inv = p*pow(ele_dist[0].params[0], -1.5)*a;
      *inv /= (j2 + 1.0); 
    } else {
      a = IntegrateRate(0, e0, e0, m, params, i0, f0,
			-RT_CI, CIRate1E);
      *inv = a/(j2+1.0);
    }
  } else {
    *inv = 0.0;
  }
  return 0;
}

double RRRateHydrogenic(double t, double z, int n, double *top) {
  const double a = 5.197E-4;
  const double b = HARTREE_EV/2.0;
  double d, x, logx, y, c, c3;
  double s0, s1, s2;
  double g0, g1, g2;
  double a0[4] = {-1.0, 2.0, -6.0, 24.0};
  double a1[4] = {-8.0/3.0, 70.0/9.0, -800.0/27.0, 11440.0/81.0};
  double a2[4] = {-1.0, 32.0/9.0, -448.0/27.0, 0.0};
  double b1[4] = {4.0/3.0, -14.0/9.0, 100.0/27.0, 0.0};
  double b2[4] = {6.0/3.0, -16.0/9.0, 128.0/27.0, 0.0};
  int i, np, one;

  c = b * z*z/t;
  c3 = pow(c, -1.0/3.0);
  x = c/(n*n);
  logx = log(x);

  if (x <= 0.02) {
    s0 = x*exp(x)*(-logx-0.5772+x);
    s1 = 0.4629*x*(1+4*x) - 1.0368*pow(x, 4.0/3.0)*(1.0+15.0*x/8.0);
    s2 = -0.0672*x*(1+3.0*x) + 0.1488*pow(x, 5.0/3.0)*(1.0+9.0*x/5.0);
    if (top) {
      g0 = x*((1.0+0.5*x)*(-logx-0.5772)+(1.0+0.75*x));
      g1 = x*(0.4629*(1+2.0*x)-0.776*pow(x, 1.0/3.0)*(1.0+15.0*x/14.0));
      g2 = x*(-0.0672*(1.0+1.5*x)+0.0893*pow(x, 2.0/3.0)*(1+9.0*x/8.0));
    }
  } else if (x >= 20.0) {
    s0 = 1.0;
    s1 = 1.0;
    s2 = 1.0;
    if (top) {
      g0 = s0/x + logx + 0.5772;
      g1 = 1.0;
      g2 = 1.0;
    }
    y = 1.0/x;
    d = y;
    for (i = 0; i < 4; i++) {
      s0 += a0[i]*d;
      s1 += a1[i]*d;
      s2 += a2[i]*d;
      if (top) {
	g1 += b1[i]*d;
	g2 += b2[i]*d;
      }
      d *= y;
    }
    y = pow(x, 1.0/3.0);
    s1 *= -0.1728*y;
    s2 *= -0.0496*y*y;
    if (top) {
      g1 = 0.926 - 0.5184*y*g1;
      g2 = 0.134 - 0.0744*y*y*g2;
    }
  } else {
    i = NSEATON;
    np = 3;
    one = 1;
    UVIP3P(np, i, log_xseaton, seaton[0], one, &logx, &s0);
    UVIP3P(np, i, log_xseaton, seaton[1], one, &logx, &s1);
    UVIP3P(np, i, log_xseaton, seaton[2], one, &logx, &s2);
    if (top) {
      UVIP3P(np, i, log_xseaton, seaton[3], one, &logx, &g0);
      UVIP3P(np, i, log_xseaton, seaton[4], one, &logx, &g1);
      UVIP3P(np, i, log_xseaton, seaton[5], one, &logx, &g2);
    }
  }
  y = s0 + c3*s1 + c3*c3*s2;
  y /= n;
  d = a*z*sqrt(c);
  if (top) {
    *top = 0.5*d*(y + g0 + c3*g1 + c3*c3*g2);
  }
  y *= d;
  return y;
}
  
double RRRate1E(double e, double eth, int np, void *p) {
  int one, n;
  double *dp;
  double x0, logx0, *x, *logx, *y, *r;
  double a, b, c, f;
  
  dp = (double *) p;
  y = dp + 1;
  x = y + np;
  logx = x + np;
  r = logx + np;

  x0 = (e + eth)/eth;
  if (x0 < x[np-1]) {
    logx0 = log(x0);
    n = 3;
    one = 1;
    UVIP3P(n, np, logx, y, one, &logx0, &f);
    f = exp(f);
  } else {
    x0 = (e + r[3])/r[3];
    logx0 = log(x0);
    a = logx0*(-dp[0]+0.5*r[1]);
    b = log((1.0 + r[2])/(sqrt(x0) + r[2]))*r[1];
    if (r[0] > 0.0) {
      f = log(r[0]*((e+eth)/(e+r[3]))) + a + b;
      f = exp(f);
    } else {
      f = 0.0;
    }
  }
  c = 2.0*PI*FINE_STRUCTURE_CONST*f*AREA_AU20;
  a = FINE_STRUCTURE_CONST*(e+eth);
  a = a*a;
  a = a/(2.0*HARTREE_EV*e);
  a /= 1.0 + 0.5*FINE_STRUCTURE_CONST2*e/HARTREE_EV;
  c *= a*VelocityFromE(e, 0.0);
  
  return c;
}

double PIRate1E(double e, double eth, int np, void *p) {
  int one, n;
  double *dp;
  double x0, logx0, *x, *logx, *y, *r;
  double a, b, c, f;
  const double factor = 1.871156686E2;

  dp = (double *) p;
  y = dp + 1;
  x = y + np;
  logx = x + np;
  r = logx + np;
  
  x0 = e/eth;
  if (x0 < x[np-1]) {
    logx0 = log(x0);
    n = 3;
    one = 1;
    UVIP3P(n, np, logx, y, one, &logx0, &f);
    f = exp(f);
  } else {
    x0 = (e-eth+r[3])/r[3];
    logx0 = log(x0);
    a = logx0*(-dp[0]+0.5*r[1]);
    b = log((1.0 + r[2])/(sqrt(x0) + r[2]))*r[1];
    if (r[0] > 0) {
      f = log(r[0]*e/(e-eth+r[3])) + a + b;
      f = exp(f);
    } else {
      f = 0.0;
    }
  }
  c = 2.0*PI*FINE_STRUCTURE_CONST*f*AREA_AU20;
  c *= factor/e;

  return c;
}

double PIRateKramers(double e, double eth, int np, void *p) {
  double c, z;
  double *dp;

  dp = (double *) p;
  z = dp[0];  
  /* 
  ** the normalization of the photon distr is erg/cm3, 
  ** divide an additional e here, and convert erg to eV by 1e12/1.6
  ** cross section is in units 1E-20, with 1e10 from s.o.l, which 
  ** gives 1E2/1.6 factor.
  */
  c = 7.91E2*(sqrt(0.5*HARTREE_EV/eth)/z)*pow((eth/e), 3.0)/e;
  c *= 3.0*1E2/1.6;

  return c;
}
 
int RRRate(double *dir, double *inv, int iinv, 
	   int j1, int j2, double e,
	   int m, double *params, int i0, int f0) {
  double e0;

  e0 = e*HARTREE_EV;	
  *dir = IntegrateRate(0, e0, 0.0, m, params, 
		       i0, f0, RT_RR, RRRate1E);
  *dir /= (j1 + 1.0);
  if (iinv) {
    *inv = IntegrateRate(1, e0, e0, m, params, 
			 i0, f0, -RT_RR, PIRate1E);
    *inv /= (j2 + 1.0);
  } else {
    *inv = 0.0;
  }
  return 0;
}

int AIRate(double *dir, double *inv, int iinv, 
	   int j1, int j2, double e, float rate) {
  double a, e0;

  *dir = RATE_AU*rate;
  if (iinv) {
    e0 = e*HARTREE_EV;
    a = ele_dist[iedist].dist(e0, ele_dist[iedist].params);
    if (a > 0.0) {
      a *= 0.5*(j1 + 1.0) * PI*PI*rate/(e*(j2+1.0));
      a *= AREA_AU20*HARTREE_EV;
      a *= VelocityFromE(e0, 0.0);
    } else {
      a = 0.0;
    }
    *inv = a;
  } else {
    *inv = 0.0;
  }
  return 0;
}

double NRRFit(int z, int nele, double t) {
  double r;

  NRRFIT(z, nele, t, &r);
  return r;
}

double RRFit(int z, int nele, double t) {
  double r;

  t = t/8.617385E-5;
  RRFIT(z, nele, t, &r);
  r *= 1E10;
  return r;
}

double NDRFit(int z, int nele, double t) {
  double r;

  NDRFIT(z, nele, t, &r);
  return r;
}

double DRFit(int z, int nele, double t) {
  double r;

  if (z > 28) return 0.0;
  if (nele < 1) return 0.0;
  if (nele >= z) return 0.0;

  DRFIT(z, nele, t, &r);
  
  r *= 1E10;
  
  return r;
}

double Gaunt0NR(double g, double u) {
  int i, j;
  double a, pg[NGAUNT0], pu[NGAUNT0];
  
  g = (log10(g) + 1.5)/2.5;
  u = (log10(u) + 1.5)/2.5;
  pg[0] = 1.0;
  pu[0] = 1.0;
  for (i = 1; i < NGAUNT0; i++) {
    pg[i] = pg[i-1]*g;
    pu[i] = pu[i-1]*u;
  }
  a = 0.0;
  for (i = 0; i < NGAUNT0; i++) {
    for (j = 0; j < NGAUNT0; j++) {
      a += bgaunt0[i][j]*pg[i]*pu[j];
    }
  }
  return a;
}

double Gaunt1NR(double g) {
  double a;

  g = log10(g);
  if (g < bgaunt1[0][0]) a = bgaunt1[1][0];
  if (g > bgaunt1[0][NGAUNT1-1]) a = bgaunt1[1][NGAUNT1-1];
  UVIP3P(3, NGAUNT1, bgaunt1[0], bgaunt1[1], 1, &g, &a);
  return a;
}
  
double BremssNR(int z, double te, double e) {
  double a, g, u, c = 1.53614E-5;
  
  g = z*z*HARTREE_EV/(2.0*te);
  if (e >= 0.0) {
    u = e/te;
    a = c*Gaunt0NR(g, u)*z*z*exp(-u)/sqrt(te);
  } else {
    a = c*z*z*sqrt(te)*Gaunt1NR(g);
  }
  
  return a;
}

double EPhFit2(int z, int nele, int is) {
  double r;

  EPHFIT2(z, nele, is, &r);
  return r;
}

double PhFit2(int z, int nele, int is, double e) {
  double r;

  PHFIT2(z, nele, is, e, &r);
  r *= 1E2;
  
  return r;
}

double RBeli(int z, int nele, double t, double *a, double *dir) {
  double total;
  
  RBELI(z, nele, t, a, dir);
  total = *a + *dir;
  
  return total;
}

double CBeli(int z, int nele, double ene, 
	     double *a, double *dir, double *err) {
  double total;
  
  CBELI(z, nele, ene, a, dir, err);
  total = *a + *dir;
  
  return total;
}

double CFit(int z, int nele, double t, double *a, double *dir) {
  double total, aa, dd;
  int is;
  
  t = t/8.617385E-5;
  is = 0;
  COLFIT(z, nele, is, t, &aa, &dd);
  CFIT(z, nele, t, &dd);
  aa *= 1E10;
  dd *= 1E10;
  total = aa + dd;
  if (a) *a = aa;
  if (dir) *dir = dd;
  
  return total;
}

double ColFit(int z, int nele, int is, double t, double *a, double *dir) {
  double total, aa, dd;
  
  t = t/8.617385E-5;
  COLFIT(z, nele, is, t, &aa, &dd);
  aa *= 1E10;
  dd *= 1E10;
  total = aa + dd;
  if (a) *a = aa;
  if (dir) *dir = dd;
  
  return total;
}

double CColFit(int z, int nele, int is, double t, double *a, double *dir) {
  double total, aa, dd;
  
  CCOLFIT(z, nele, is, t, &aa, &dd);
  total = aa + dd;
  if (a) *a = aa;
  if (dir) *dir = dd;
  
  return total;
}

double EColFit(int z, int nele, int is) {
  double e;
  ECOLFIT(z, nele, is, &e);
  return e;
}

double EBeli(int z, int nele) {
  double e;
  EBELI(z, nele, &e);
  return e;
}

double Ionis(int z, int nele, double t, double *a, double *dir, int m) {
  double total, aa, dd;

  if (m == 0) {
    nele = z+1 - nele;
    IONIS(z, nele, t, &aa, &dd, &total);
    aa *= 1E10;
    dd *= 1E10;
    total *= 1E10;
  } else if (m == 1) {
    total = ColFit(z, nele, 0, t, &aa, &dd);
  } else if (m == 2) {
    total = CFit(z, nele, t, &aa, &dd);
  } else if (m == 3) {
    total = RBeli(z, nele, t, &aa, &dd);
    if (total < 0) total = ColFit(z,nele, 0, t, &aa, &dd);
  } else {
    printf("unrecognized mode in Ionis\n");
    return 0.0;
  }
  if (a) *a = aa;
  if (dir) *dir = dd;

  if (total < 1E-30) total = 1E-30;

  return total;
}

double Recomb(int z, int nele, double t, double *rr, double *dr, int m) {
  double total, r, d;
  
  if (m == 0) {
    nele = z+1 - nele;
    t = t/8.617385E-5;
    if (z == 26) {
      RECOMBFE(z, nele, t, &r, &d);
    } else {
      RECOMB(z, nele, t, &r, &d);
    }
    r *= 1E10;
    d *= 1E10;
  } else if (m == 1) {
    r = RRFit(z, nele, t);
    d = DRFit(z, nele, t);
  } else if (m == 2) {
    r = NRRFit(z, nele, t);
    if (r < 0) r = RRFit(z, nele, t);
    d = NDRFit(z, nele, t);
    if (d < 0) d = DRFit(z, nele, t);
  } else {
    printf("Unrecognized Mode in Recomb\n");
    return 0.0;
  }
  if (rr) *rr = r;
  if (dr) *dr = d;
  total = r + d;

  if (total < 1E-30) total = 1E-30;
  return total;
}

int FracAbund(int z, double t, double *a, int im, int rm) {
  int nele, k;
  double *ir, *rr, c;

  ir = (double *) malloc(sizeof(double)*(z+1));
  rr = (double *) malloc(sizeof(double)*(z+1));
  
  for (nele = 0; nele <= z; nele++) {
    if (nele != 0) {
      ir[nele] = Ionis(z, nele, t, NULL, NULL, im);
    }
    if (nele != z) {
      rr[nele] = Recomb(z, nele, t, NULL, NULL, rm);
    }
  }

  a[0] = 1.0;
  for (nele = 1; nele <= z; nele++) {
    a[nele] = a[nele-1]*(rr[nele-1]/ir[nele]);
    if (a[nele] > 1E20) {
      for (k = 0; k < nele; k++) {
	a[k] /= a[nele];
      }
      a[nele] = 1.0;
    }
  }
  c = 0.0;
  for (nele = 0; nele <= z; nele++) {
    c += a[nele];
  }
  for (nele = 0; nele <= z; nele++) {
    a[nele] /= c;
  }
  
  free(ir);
  free(rr);
  return 0;
}

double MaxAbund(int z, int nele, double *a, double eps, int im, int rm) {
  double amax, da;
  double t1, t2, t;
  double dt;
  
  t1 = 0.1*HARTREE_EV;
  t2 = z*z*HARTREE_EV;
  if (nele == 0) {
    FracAbund(z, t2, a, im, rm);
    return t2;
  }
  if (nele == z) {
    FracAbund(z, t1, a, im, rm);
    return t1;
  }
  amax = 0.0;
  t = 0.5*(t1+t2);
  dt = t*eps;
  while (t2 - t1 > dt) {
    FracAbund(z, t, a, im, rm);
    amax = a[nele];
    FracAbund(z, t+dt, a, im, rm);
    da = a[nele] - amax;
    if (da > 0.0) {
      t1 = t;
    } else {
      t2 = t;
    }
    t = 0.5*(t1+t2);
    dt = t*eps;
  }
  
  return t;  
}

/*
** two-photon rates are taken from G. W. F. Drake, PRA, 34, 2871, 1986.
*/
double TwoPhotonRate(double z, int t) {
  double a, a2, a4, z6;
  
  switch (t) {
  case 0: /* 2S_1/2 of H-like ion */
    z6 = z*z;
    z6 = z6*z6*z6;
    a = FINE_STRUCTURE_CONST*z;
    a2 = a*a;
    a4 = a2*a2;
    a = 8.22943*z6*(1.0 + 3.9448*a2 - 2.04*a4)/(1.0 + 4.6019*a2);
    break;
  case 1: /* 1s2s S_0 of He-like ion */
    a = (z - 0.806389);
    z6 = a*a;
    z6 = z6*z6*z6;
    a = FINE_STRUCTURE_CONST*a;
    a2 = a*a;
    a4 = (z+2.5);
    a4 = a4*a4;
    a = 16.458762*(z6*(1.0 + 1.539/a4) - 
		   z6*a2*(0.6571 + 2.04*a2)/(1.0 + 4.6019*a2));
    break;
  default:
    a = 0.0;
    break;
  }

  return a;
}


/*
** The following routines should be modified to 
** add more electron or photon energy distributions.
*/
static double Gaussian(double e, double *p) {
  double x;
  const double gauss_const = 0.39894228;

  if (e > p[3] || e < p[2]) return 0.0;
  x = (e - p[0])/p[1];
  x = 0.5*x*x;
  x = gauss_const * exp(-x)/p[1];

  return x;
}

static double MonoEnergy(double e, double *p) {
  double x;
  if (e < p[2]) return 0.0;
  if (e > p[3]) return 0.0;
  x = 1.0/p[1];
  return x;
}

static double SMaxwell(double e, double *p) {
  const double c = 0.282095;
  double x, y, sx, sy;

  if (e > p[3] || e < p[2]) return 0.0;
  x = e/p[0];
  y = p[1]/p[0];
  sx = sqrt(x);
  sy = sqrt(y);
  x = c*(exp(-(sx-sy)*(sx-sy))-exp(-(sx+sy)*(sx+sy)))/(sy*p[0]);
  return x;
}

static double Maxwell(double e, double *p) {
  double x;
  const double maxwell_const = 1.12837967;
  
  if (e > p[2] || e < p[1]) return 0.0;
  x = e/p[0];
  x = maxwell_const * sqrt(x) * exp(-x)/p[0];
  return x;
}

static double DoubleMaxwell(double e, double *p) {
  double x0, x1, x;
  const double maxwell_const = 1.12837967;
  
  if (e > p[4] || e < p[3]) return 0.0;
  x0 = e/p[0];
  x1 = e/p[1];
  x = ((1-p[2])*sqrt(x0)*exp(-x0)/p[0] + p[2]*sqrt(x1)*exp(-x1)/p[1]);
  x *= maxwell_const;
  return x;
}
  
static double MaxPower(double e, double *p) {
  double x, y, logy, py, g, g1, c1, c2, f;
  const double maxwell_const = 1.12837967;

  if (e < 0.0) {
    y = p[4]/p[3];
    logy = log(y);
    py = pow(y, 1.0-p[1]);
    if (p[1] == 1.0) {
      g = logy;
    } else {
      g = (1.0 - py)/(p[1] - 1.0);
    }
    if (p[1] == 2.0) {
      g1 = logy;
    } else {
      g1 = (1.0 - y*py)/(p[1] - 2.0);
    }
    y = 1.5*p[2]*p[0]/(p[3]*g1);
    c1 = maxwell_const/(1.0+y*g);
    c2 = y/(1.0+y*g);
    p[7] = c1;
    p[8] = c2;
    return 0.0;
  } 
  c1 = p[7];
  c2 = p[8];
  if (e > p[5] && e < p[6]) {
    x = e/p[0];
    f = c1*sqrt(x)*exp(-x)/p[0];
  } else {
    f = 0.0;
  }
  if (e > p[3] && e < p[4]) {
    x = e/p[3];
    f += c2*pow(x,-p[1])/p[3];
  } 

  return f;
}
  
static double PowerLaw(double e, double *p) {
  double x;

  if (e > p[2] || e < p[1]) return 0.0;
  if (p[0] == 1.0) {
    x = 1.0/(log(p[2]) - log(p[1]));
  } else {
    x = 1.0 - p[0];
    x = x/(pow(p[2], x) - pow(p[1], x));
  }
  x = x*pow(e, -p[0]);

  return x;
}

static double IncompleteGamma(double a, double x) {
  double gold=0.0, fac=1.0, b1=1.0, b0=0.0, a0=1.0;
  double a1, an, ana, anf, g, ap, del, sum, y;
  int n;

  if (x < 0.0 || a <= 0.0) return 0.0;
  y = DLOGAM(a);
  if (x < (a+1.0)) {
    ap = a;
    del = sum = 1.0/a;
    for (n = 1; n <= 100; n++) {
      ap += 1.0;
      del *= x/ap;
      sum += del;
      if (fabs(del) < fabs(sum)*EPS6) {
	break;
      }      
    }
    y = sum*exp(-x+a*log(x));
  } else {
    a1 = x;
    for (n = 1; n <= 100; n++) {
      an = n;
      ana = an-a;
      a0 = (a1 + a0*ana)*fac;
      b0 = (b1 + b0*ana)*fac;
      anf = an*fac;
      a1 = x*a0 + anf*a1;
      b1 = x*b0 + anf*b1;
      if (a1) {
	fac = 1.0/a1;
	g = b1*fac;
	if (fabs((g-gold)/g) < EPS6) {
	  break;
	}
	gold = g;
      }
    }
    y = exp(y)-g*exp(-x+a*log(x));
  }
  return y;
}
  
static double Hybrid(double e, double *p) {
  double x, y, c, x3, x4;

  if (e < 0.0) {
    x4 = p[4]/p[0];
    x3 = p[3]/p[0];
    x = e/p[0];
    c = IncompleteGamma(1.5, p[2]) - IncompleteGamma(1.5, x3);
    if (p[1] + 1.0 != 1.0) {
      c += pow(p[2], 1.5)*exp(-p[2])*(1.0-pow(x4/p[2],1.0-p[1]))/(p[1]-1.0);
    } else {
      c += pow(p[2], 1.5)*exp(-p[2])*log(x4/p[2]);
    }
    c = 1.0/c;
    p[5] = c;
    return 0.0;
  }
  c = p[5];
  x = e/p[0];
  if (e > p[4] || e < p[3]) return 0.0;
  if (x <= p[2]) {
    y = c*sqrt(x)*exp(-x)/p[0];
  } else {
    y = c*sqrt(p[2])*exp(-p[2])*pow(x/p[2], -p[1])/p[0];
  }

  return y;
}
  
static double BBody(double e, double *p) {
  double c = 21.0989; /* 8*PI*eV^4/(hc)^3 in erg/cm^3 */
  double x;

  if (e > p[2] || e < p[1]) return 0.0;
  
  x = e/p[0];
  if (x < EPS5) {    
    return c*e*e*e/x;
  } else if (x > 50.0) {
    return c*(((exp(-x)*e)*e)*e);
  } else {
    return c*e*e*e/(exp(x) - 1.0);
  }
}

static double InterpDist(double e, double *p) {
  int ng, ix, iy;
  double *eg, *dg, r;

  ng = p[0];
  ix = p[1];
  iy = p[2];
  eg = &(p[3]);
  dg = &(p[3+ng]);
  if (ix) e = log(e);
  if (e < eg[0] || e > eg[ng-1]) return 0.0;
  UVIP3P(3, ng, eg, dg, 1, &e, &r);
  if (iy) r = exp(r);
  return r;
}

int DistFromFile(char *fn, double **p) {
  int n, k, np, i, ix, iy;
  FILE *f;
  double *eg, *dg;

  f = fopen(fn, "r");
  if (f == NULL) {
    printf("cannot open file %s\n", fn);
    return -1;
  }
  np = -1;
  k = fscanf(f, "%d %d %d\n", &n, &ix, &iy);
  if (k != 3) goto ERROR;
  *p = (double *) malloc(sizeof(double)*(2*n+3));
  (*p)[0] = n;
  (*p)[1] = ix;
  (*p)[2] = iy;
  eg = (*p)+3;
  dg = (*p)+3+n;
  for (i = 0; i < n; i++) {
    k = fscanf(f, "%lf %lf\n", &(eg[i]), &(dg[i]));
    if (k != 2) goto ERROR;
  }
  if (i != n) goto ERROR;

  np = 2*n + 3;
 ERROR:
  fclose(f);
  
  return np;
}

int EleDist(char *fn, int n) {
  FILE *f;
  double y, e, de, emin, emax, *p;
  int np, i;
  DISTRIBUTION *d;

  f = fopen(fn, "w");
  if (f == NULL) {
    printf("cannot open file %s\n", fn);
    return -1;
  }

  d = ele_dist+iedist;
  p = d->params;
  np = d->nparams;
  if (iedist == MAX_DIST-1) {
    np = (np-3)/2;
    emin = p[3];
    emax = p[3+np-1];
    if (p[1] != 0) {
      emin = exp(emin);
      emax = exp(emax);
    }
  } else {
    emin = p[np-2];
    emax = p[np-1];  
  }
  de = (log(emax) - log(emin))/(n-1);
  de = exp(de);
  fprintf(f, "# electron dist=%d\n", iedist);
  if (iedist != MAX_DIST-1) {
    for (i = 0; i < np; i++) {
      fprintf(f, "# P[%d]=%10.3E\n", i, p[i]);
    } 
  }
  e = emin;
  for (i = 0; i < n; i++) {
    y = d->dist(e, p);
    fprintf(f, "%15.8E %15.8E\n", e, y);
    e *= de;
  }
  
  fclose(f);
  return 0;
}

int CxtDist(char *fn, int n) {
  FILE *f;
  double y, e, de, emin, emax, *p;
  int np, i;
  DISTRIBUTION *d;

  f = fopen(fn, "w");
  if (f == NULL) {
    printf("cannot open file %s\n", fn);
    return -1;
  }

  d = cxt_dist+ixdist;
  p = d->params;
  np = d->nparams;
  if (ixdist == MAX_DIST-1) {
    np = (np-3)/2;
    emin = p[3];
    emax = p[3+np-1];
    if (p[1] != 0) {
      emin = exp(emin);
      emax = exp(emax);
    }
  } else {
    emin = p[np-2];
    emax = p[np-1];  
  }
  de = (log(emax) - log(emin))/(n-1);
  de = exp(de);
  fprintf(f, "# cxt dist=%d\n", ixdist);
  if (ixdist != MAX_DIST-1) {
    for (i = 0; i < np; i++) {
      fprintf(f, "# P[%d]=%10.3E\n", i, p[i]);
    } 
  }
  e = emin;
  for (i = 0; i < n; i++) {
    y = d->dist(e, p);
    fprintf(f, "%15.8E %15.8E\n", e, y);
    e *= de;
  }
  
  fclose(f);
  return 0;
}

int SetEleDist(int i, int np, double *p0) {
  int k;
  double *p, c1, c2;

  if (i == -1) {
    i = MAX_DIST-1;
    ele_dist[i].nparams = np;
    if (ele_dist[i].dist != NULL) {
      free(ele_dist[i].params);
    }
    ele_dist[i].params = (double *) malloc(sizeof(double)*np);
    ele_dist[i].dist = InterpDist;
  }

  if (ele_dist[i].dist == NULL) {
    printf("Electron Dist. does not exist\n");
    return -1;
  }
  if (ele_dist[i].nparams != np) {
    printf("Num of Params for Electron Dist. %d does not match\n", i);
    return -1;
  }

  ele_dist[i].xlog = -1;
  iedist = i;
  for (k = 0; k < np; k++) {
    ele_dist[i].params[k] = p0[k];
  }
  p = ele_dist[i].params;
  switch (i) {
  case 0:
    /* Default parameters for Maxwellian */
    if (p[2] <= 0.0) {
      p[2] = 1E2*p[0];
    }
    if (p[1] <= 0.0) {
      p[1] = 1E-20*p[0];
    }
    break;
  case 1:
    /* Gaussian */
    if (p[3] <= 0.0) {
      p[3] = p[0] + 10.0*p[1];
    }
    if (p[2] <= 0.0) {
      p[2] = p[0] - 10.0*p[1];
      if (p[2] < 0) p[2] = 0.0;
    }
    ele_dist[i].xlog = 0;
    break;
  case 2:
    /* MaxPower */
    if (p[6] <= 0.0) {
      p[6] = 1E2*p[0];
    }
    if (p[5] <= 0.0) {
      p[5] = 1E-20*p[0];
    }
    MaxPower(-1.0, p);
    break;    
  case 4:
    /* shifted Maxwellian */
    if (p[3] <= 0.0) {
      p[3] = Max(p[0], p[1])*100.0;
    }
    if (p[2] <= 0.0) {
      p[2] = 1E-20;
    }
    break;
  case 5:
    /* Hybrid Maxwell, powerlaw amd maxwell continuously matched. */
    if (p[4] <= 0.0) {
      p[4] = 1E3*p[0]*p[2];
    }
    if (p[3] <= 0.0) {
      p[3] = 1E-20*p[0];
    }
    Hybrid(-1.0, p);
    break;
  case 6:
    /* double Maxwellian */
    if (p[0] < p[1]) {
      c1 = p[0];
      c2 = p[1];
    } else {
      c1 = p[1];
      c2 = p[0];
    }
    if (p[2]+1.0 == 1.0) {
      c1 = p[0];
      c2 = p[0];
    }
    if (p[2] == 1.0) {
      c1 = p[1];
      c2 = p[1];
    }
    if (p[4] <= 0.0) {
      p[4] = 1E2*c2;
    }
    if (p[3] <= 0.0) {
      p[3] = 1E-20*c1;
    }
    break;
  case 7:
    /* Mono Energy */
    p[3] = p[0] + 0.5*p[1];
    p[2] = p[0] - 0.5*p[1];
    ele_dist[i].xlog = 0;
    break;    
  default:
    break;
  }

  ThreeBodyDist();

  return 0;
}

int SetCxtDist(int i, int np, double *p0) {
  int k;
  double *p, c1, c2;

  if (i == -1) {
    i = MAX_DIST-1;
    cxt_dist[i].nparams = np;
    if (cxt_dist[i].dist != NULL) {
      free(cxt_dist[i].params);
    }
    cxt_dist[i].params = (double *) malloc(sizeof(double)*np);
    cxt_dist[i].dist = InterpDist;
  }

  if (cxt_dist[i].dist == NULL) {
    printf("Electron Dist. does not exist\n");
    return -1;
  }
  if (cxt_dist[i].nparams != np) {
    printf("Num of Params for CXT Dist. %d does not match\n", i);
    return -1;
  }

  cxt_dist[i].xlog = -1;
  ixdist = i;
  for (k = 0; k < np; k++) {
    cxt_dist[i].params[k] = p0[k];
  }
  p = cxt_dist[i].params;
  switch (i) {
  case 0:
    /* Default parameters for Maxwellian */
    if (p[2] <= 0.0) {
      p[2] = 1E2*p[0];
    }
    if (p[1] <= 0.0) {
      p[1] = 1E-20*p[0];
    }
    break;
  case 1:
    /* Gaussian */
    if (p[3] <= 0.0) {
      p[3] = p[0] + 10.0*p[1];
    }
    if (p[2] <= 0.0) {
      p[2] = p[0] - 10.0*p[1];
      if (p[2] < 0) p[2] = 0.0;
    }
    cxt_dist[i].xlog = 0;
    break;
  case 2:
    /* MaxPower */
    if (p[6] <= 0.0) {
      p[6] = 1E2*p[0];
    }
    if (p[5] <= 0.0) {
      p[5] = 1E-20*p[0];
    }
    MaxPower(-1.0, p);
    break;    
  case 4:
    /* shifted Maxwellian */
    if (p[3] <= 0.0) {
      p[3] = Max(p[0], p[1])*100.0;
    }
    if (p[2] <= 0.0) {
      p[2] = 1E-20;
    }
    break;
  case 5:
    /* Hybrid Maxwell, powerlaw amd maxwell continuously matched. */
    if (p[4] <= 0.0) {
      p[4] = 1E3*p[0]*p[2];
    }
    if (p[3] <= 0.0) {
      p[3] = 1E-20*p[0];
    }
    Hybrid(-1.0, p);
    break;
  case 6:
    /* double Maxwellian */
    if (p[0] < p[1]) {
      c1 = p[0];
      c2 = p[1];
    } else {
      c1 = p[1];
      c2 = p[0];
    }
    if (p[2]+1.0 == 1.0) {
      c1 = p[0];
      c2 = p[0];
    }
    if (p[2] == 1.0) {
      c1 = p[1];
      c2 = p[1];
    }
    if (p[4] <= 0.0) {
      p[4] = 1E2*c2;
    }
    if (p[3] <= 0.0) {
      p[3] = 1E-20*c1;
    }
    break;
  case 7:
    /* Mono Energy */
    p[3] = p[0] + 0.5*p[1];
    p[2] = p[0] - 0.5*p[1];
    cxt_dist[i].xlog = 0;
    break;        
  default:
    break;
  }

  return 0;
}

int PhoDist(char *fn, int n) {
  FILE *f;
  double y, e, de, emin, emax, *p;
  int np, i;
  DISTRIBUTION *d;

  f = fopen(fn, "w");
  if (f == NULL) {
    printf("cannot open file %s\n", fn);
    return -1;
  }

  d = pho_dist+ipdist;
  p = d->params;
  np = d->nparams;
  if (ipdist == MAX_DIST-1) {
    np = (np-3)/2;
    emin = p[3];
    emax = p[3+np-1];
    if (p[1] != 0) {
      emin = exp(emin);
      emax = exp(emax);
    }
  } else {
    emin = p[np-2];
    emax = p[np-1];
  }
  de = (log(emax) - log(emin))/(n-1);
  de = exp(de);
  fprintf(f, "# photon dist=%d\n", ipdist);
  if (ipdist != MAX_DIST-1) {
    for (i = 0; i < np; i++) {
      fprintf(f, "# P[%d]=%10.3E\n", i, p[i]);
    } 
  }
  e = emin;
  for (i = 0; i < n; i++) {
    y = d->dist(e, p);
    fprintf(f, "%15.8E %15.8E\n", e, y);
    e *= de;
  }
  
  fclose(f);
  return 0;
}

int SetPhoDist(int i, int np, double *p) {
  int k;

  if (i == -1) {
    i = MAX_DIST-1;
    pho_dist[i].nparams = np;
    if (pho_dist[i].dist != NULL) {
      free(pho_dist[i].params);
    }
    pho_dist[i].params = (double *) malloc(sizeof(double)*np);
    pho_dist[i].dist = InterpDist;
  }

  if (pho_dist[i].dist == NULL) {
    printf("Photon Dist. does not exist\n");
    return -1;
  }

  if (pho_dist[i].nparams != np) {
    printf("Num of Params for Photon Dist. %d does not match\n", i);
    return -1;
  }

  pho_dist[i].xlog = -1;
  switch (i) {
  case 0:
    /* Default parameters for BBody */
    if (p[2] <= 0.0) {
      p[2] = 1E2*p[0];
    }
    if (p[1] <= 0.0) {
      p[1] = 1E-20*p[0];
    }
    break;
  }

  ipdist = i;
  for (k = 0; k < np; k++) {
    pho_dist[i].params[k] = p[k];
  }

  return 0;
}

int InitRates(void) {
  int i;

  iedist = 0;
  ipdist = 0;
  ixdist = 0;
  
  i = 0; /* Maxwellian */
  ele_dist[i].nparams = 3;
  ele_dist[i].params = (double *) malloc(sizeof(double)*3);
  ele_dist[i].params[0] = 1E3;
  ele_dist[i].params[1] = 1E-10;
  ele_dist[i].params[2] = 1E10;
  ele_dist[i].dist = Maxwell;

  i++; /* Gaussian */
  ele_dist[i].nparams = 4;
  ele_dist[i].params = (double *) malloc(sizeof(double)*4);
  ele_dist[i].params[0] = 1.0E3;
  ele_dist[i].params[1] = 50.0;
  ele_dist[i].params[2] = 1E3-250.0;
  ele_dist[i].params[3] = 1E3+250.0;
  ele_dist[i].dist = Gaussian;

  i++; /* MaxPower */
  ele_dist[i].nparams = 7;
  ele_dist[i].params = (double *) malloc(sizeof(double)*(7+2));
  ele_dist[i].params[0] = 1.0E3;
  ele_dist[i].params[1] = 1.5;
  ele_dist[i].params[2] = 0.1;
  ele_dist[i].params[3] = 1.0E3;
  ele_dist[i].params[4] = 1E7;
  ele_dist[i].params[5] = 1E-10;
  ele_dist[i].params[6] = 1E10;
  ele_dist[i].dist = MaxPower;

  i++; /* Power */
  ele_dist[i].nparams = 3;
  ele_dist[i].params = (double *) malloc(sizeof(double)*3);
  ele_dist[i].params[0] = 1.5;
  ele_dist[i].params[1] = 1.0E3;
  ele_dist[i].params[2] = 1E7;
  ele_dist[i].dist = PowerLaw;

  i++; /* shifted Maxwellian */
  ele_dist[i].nparams = 4;
  ele_dist[i].params = (double *) malloc(sizeof(double)*4);
  ele_dist[i].params[0] = 1e1;
  ele_dist[i].params[1] = 1e3;
  ele_dist[i].params[2] = 1e-20;
  ele_dist[i].params[3] = 1e5;  
  ele_dist[i].dist = SMaxwell;

  i++; /* Maxwellian with powerlaw tail, smooth match*/
  ele_dist[i].nparams = 5;
  ele_dist[i].params = (double *) malloc(sizeof(double)*(5+1));
  ele_dist[i].params[0] = 1e2;
  ele_dist[i].params[1] = 2.0;
  ele_dist[i].params[2] = 3.0;
  ele_dist[i].params[3] = 1E-20;
  ele_dist[i].params[4] = 1E8;
  ele_dist[i].dist = Hybrid;

  i++; /* double Maxwellian */
  ele_dist[i].nparams = 5;
  ele_dist[i].params = (double *) malloc(sizeof(double)*5);
  ele_dist[i].params[0] = 1e2;
  ele_dist[i].params[1] = 1e4;
  ele_dist[i].params[2] = 0.1;
  ele_dist[i].params[3] = 1E-10;
  ele_dist[i].params[4] = 1E10;
  ele_dist[i].dist = DoubleMaxwell;
  
  i++; /* MonoEnergy */
  ele_dist[i].nparams = 4;
  ele_dist[i].params = (double *) malloc(sizeof(double)*4);
  ele_dist[i].params[0] = 1.0E3;
  ele_dist[i].params[1] = 10.0;
  ele_dist[i].params[2] = 1E3-5.0;
  ele_dist[i].params[3] = 1E3+5.0;
  ele_dist[i].dist = MonoEnergy;

  i++;

  for (; i < MAX_DIST; i++) {
    ele_dist[i].nparams = 0;
    ele_dist[i].params = NULL;
    ele_dist[i].dist = NULL;
  }

  /* photon energy distribution */
  i = 0; /* black body */
  pho_dist[i].nparams = 3;
  pho_dist[i].params = (double *) malloc(sizeof(double)*3);
  pho_dist[i].params[0] = 1E1;
  pho_dist[i].params[1] = 1E-10;
  pho_dist[i].params[2] = 1E5;
  pho_dist[i].dist = BBody;

  i++; /* Power Law */
  pho_dist[i].nparams = 3;
  pho_dist[i].params = (double *) malloc(sizeof(double)*3);
  pho_dist[i].params[0] = 2.0;
  pho_dist[i].params[1] = 1.0E1;
  pho_dist[i].params[2] = 1.0E5;
  pho_dist[i].dist = PowerLaw;
  i++;
  
  for (; i < MAX_DIST; i++) {
    pho_dist[i].nparams = 0;
    pho_dist[i].params = NULL;
    pho_dist[i].dist = NULL;
  }

  //CX distribution  
  i = 0; /* Maxwellian */
  cxt_dist[i].nparams = 3;
  cxt_dist[i].params = (double *) malloc(sizeof(double)*3);
  cxt_dist[i].params[0] = 1E3;
  cxt_dist[i].params[1] = 1E-10;
  cxt_dist[i].params[2] = 1E10;
  cxt_dist[i].dist = Maxwell;

  i++; /* Gaussian */
  cxt_dist[i].nparams = 4;
  cxt_dist[i].params = (double *) malloc(sizeof(double)*4);
  cxt_dist[i].params[0] = 1.0E3;
  cxt_dist[i].params[1] = 50.0;
  cxt_dist[i].params[2] = 1E3-250.0;
  cxt_dist[i].params[3] = 1E3+250.0;
  cxt_dist[i].dist = Gaussian;

  i++; /* MaxPower */
  cxt_dist[i].nparams = 7;
  cxt_dist[i].params = (double *) malloc(sizeof(double)*(7+2));
  cxt_dist[i].params[0] = 1.0E3;
  cxt_dist[i].params[1] = 1.5;
  cxt_dist[i].params[2] = 0.1;
  cxt_dist[i].params[3] = 1.0E3;
  cxt_dist[i].params[4] = 1E7;
  cxt_dist[i].params[5] = 1E-10;
  cxt_dist[i].params[6] = 1E10;
  cxt_dist[i].dist = MaxPower;

  i++; /* Power */
  cxt_dist[i].nparams = 3;
  cxt_dist[i].params = (double *) malloc(sizeof(double)*3);
  cxt_dist[i].params[0] = 1.5;
  cxt_dist[i].params[1] = 1.0E3;
  cxt_dist[i].params[2] = 1E7;
  cxt_dist[i].dist = PowerLaw;

  i++; /* shifted Maxwellian */
  cxt_dist[i].nparams = 4;
  cxt_dist[i].params = (double *) malloc(sizeof(double)*4);
  cxt_dist[i].params[0] = 1e1;
  cxt_dist[i].params[1] = 1e3;
  cxt_dist[i].params[2] = 1e-20;
  cxt_dist[i].params[3] = 1e5;  
  cxt_dist[i].dist = SMaxwell;

  i++; /* Maxwellian with powerlaw tail, smooth match*/
  cxt_dist[i].nparams = 5;
  cxt_dist[i].params = (double *) malloc(sizeof(double)*(5+1));
  cxt_dist[i].params[0] = 1e2;
  cxt_dist[i].params[1] = 2.0;
  cxt_dist[i].params[2] = 3.0;
  cxt_dist[i].params[3] = 1E-20;
  cxt_dist[i].params[4] = 1E8;
  cxt_dist[i].dist = Hybrid;

  i++; /* double Maxwellian */
  cxt_dist[i].nparams = 5;
  cxt_dist[i].params = (double *) malloc(sizeof(double)*5);
  cxt_dist[i].params[0] = 1e2;
  cxt_dist[i].params[1] = 1e4;
  cxt_dist[i].params[2] = 0.1;
  cxt_dist[i].params[3] = 1E-10;
  cxt_dist[i].params[4] = 1E10;
  cxt_dist[i].dist = DoubleMaxwell;

  i++; /* MonoEnergy */
  cxt_dist[i].nparams = 4;
  cxt_dist[i].params = (double *) malloc(sizeof(double)*4);
  cxt_dist[i].params[0] = 1.0E3;
  cxt_dist[i].params[1] = 10.0;
  cxt_dist[i].params[2] = 1E3-5.0;
  cxt_dist[i].params[3] = 1E3+5.0;
  cxt_dist[i].dist = MonoEnergy;

  i++;

  for (; i < MAX_DIST; i++) {
    cxt_dist[i].nparams = 0;
    cxt_dist[i].params = NULL;
    cxt_dist[i].dist = NULL;
  }

  rate_epsabs = EPS8;
  rate_epsrel = EPS3;
  rate_iprint = 1;

  for (i = 0; i < NSEATON; i++) {
    log_xseaton[i] = log(xseaton[i]);
  }

  _kronos_cx[0].nmax = 0;
  _kronos_cx[1].nmax = 0;
  return 0;
}

int LFromSym(char s) {
  char sym[] = "spdfghiklmnoqrtuvwxyz";
  int i;
  for (i = 0; i < 21; i++) {
    if (sym[i] == s) break;
  }
  return i;
}

KRONOS *InitKronos(int z, int k, int nmax, int nep, char *cxm, char *tgt) {
  KRONOS *cx = &_kronos_cx[k];
  int i, j;

  if (cx->nmax > 0 && cx->nep > 0) {
    free(cx->lnfac);
    free(cx->idn);
    free(cx->ep);
    free(cx->rcx);
    free(cx->rx);
    for (i = 0; i < cx->ncx; i++) {
      free(cx->cx0[i]);
      if (cx->cx1) free(cx->cx1[i]);
    }
    free(cx->cx0);
    if (cx->cx1) free(cx->cx1);
    cx->cx0 = NULL;
    cx->cx1 = NULL;
    cx->ep = NULL;
    cx->rcx = NULL;
    cx->rx = NULL;
    cx->idn = NULL;
  }
  cx->nmax = nmax;
  cx->z = z;
  cx->k = k;
  cx->nep = nep;
  strcpy(cx->cxm, cxm);
  strcpy(cx->tgt, tgt);
  cx->ncx = nmax*(nmax+1)/2;
  cx->idn = malloc(sizeof(int)*nmax);
  cx->lnfac = malloc(sizeof(double)*(2*nmax+1));
  cx->lnfac[0] = 0.0;
  for (i = 1; i < 2*nmax+1; i++) {
    cx->lnfac[i] = cx->lnfac[i-1]+log(i);
  }
  cx->ep = malloc(sizeof(double)*nep);
  cx->rcx = malloc(sizeof(double)*cx->ncx);
  cx->rx = malloc(sizeof(double)*cx->ncx);
  cx->cx0 = malloc(sizeof(double *)*cx->ncx);
  if (k > 0) cx->cx1 = malloc(sizeof(double *)*cx->ncx);
  else cx->cx1 = NULL;
  for (i = 0; i < cx->ncx; i++) {
    cx->rcx[i] = 0;
    cx->rx[i] = 0;
    cx->cx0[i] = malloc(sizeof(double)*nep);
    if (cx->cx1) cx->cx1[i] = malloc(sizeof(double)*nep);
    for (j = 0; j < nep; j++) {
      cx->cx0[i][j] = 0.0;
      if (cx->cx1) {
	cx->cx1[i][j] = 0.0;
      }
    }
  }
  for (i = 0; i < nmax; i++) {
    cx->idn[i] = i*(i+1)/2;
  }
  return cx;
}

int ReadKronos(char *dn, int z, int k,
	       char *prj, char *tgt, char *cxm, int md, int ilog) {
  char kfn[1024], kfn1[1024];
  char prj0[3], prj1[3];
  char tgt0[16], tgt1[16];
  char cxm0[32];
#define ncxm 7
  char cxms[ncxm][32] = {
    "rcmd", "qmocc", "mocc", "aocc", "ctmc", "mclz", "faclz"
  };
  FILE *f;
  KRONOS *cx;
  char *s, buf[8192], es[16], es1[128];
  int i;
  double dm, wkn, vr;

  if (k != 0 && k != 1) {
    printf("only cx onto bare and H-like in Kronos db: %d\n", k);
    return -1;
  }
  StrLowerUpper(prj0, prj1, prj, 3);
  StrLowerUpper(tgt0, tgt1, tgt, 16);
  sprintf(kfn, "%s/CXDatabase/elements.dat", dn);
  f = fopen(kfn, "r");
  if (f == NULL) {
    printf("cannot find KRONOS elements file: %s\n", kfn);
    return -1;
  }
  double tmass = 0;
  double pmass = 0;
  while (1) {
    if (NULL == fgets(buf, 8192, f)) break;
    sscanf(buf, "%d %s %s %lg", &i, es, es1, &dm);
    StrUpper(es1, es, 16);
    if (0 == strcmp(es1, tgt1)) {
      tmass = dm;
    } else if (i == z) {
      pmass = dm;
    }
    if (tmass > 0 && pmass > 0) break;
  }
  double rmass = pmass*tmass/(pmass+tmass);
  fclose(f);
  if (cxm && strlen(cxm) == 0) cxm = NULL;
  if (cxm) StrLower(cxm0, cxm, 32);
  int c = z-k;
  f = NULL;
  for (i = 0; i < ncxm; i++) {
    if (cxm != NULL && !strstr(cxm0, cxms[i])) continue;
    if (k == 0) {
      sprintf(kfn1, "%s/CXDatabase/Projectile_Ions/%s/Charge/%d/Targets/%s/%s%d+%s_sec_%s_nres.cs",
	      dn, prj, c, tgt, prj0, c, tgt0, cxms[i]);
      sprintf(kfn, "%s/CXDatabase/Projectile_Ions/%s/Charge/%d/Targets/%s/%s%d+%s_sec_%s.cs",
	      dn, prj, c, tgt, prj0, c, tgt0, cxms[i]);
      if (cxm != NULL && strstr(cxm0, "nres")) {
	kfn[0] = '\0';
      }
    } else {
      sprintf(kfn, "%s/CXDatabase/Projectile_Ions/%s/Charge/%d/Targets/%s/%s%d+%s_sec_%s.cs",
	      dn, prj, c, tgt, prj0, c, tgt0, cxms[i]);
      kfn1[0] = '\0';
    }
    f = NULL;
    if (kfn[0]) {
      f = fopen(kfn, "r");
    }
    if (!f && kfn1[0]) {
      strcpy(kfn, kfn1);
      f = fopen(kfn, "r");
    }
    if (f) break;
  }
  if (f == NULL) {
    printf("cannot find kronos data file: %s %s %s %d %d\n",
	   prj, tgt, cxm, z, k);
    return -1;
  }

  strcpy(cxm0, cxms[i]);
  int nep = 0; 
  int nb, ns, n, nr, nn, nmax = 0, nres = 0;
  int na[1000], ka[1000],sa[1000];
  int ncx = 0;
  while(1) {
    if (NULL == fgets(buf, 8192, f)) break;
    nb = strlen(buf);
    if (nb < 2) continue;
    StrReplace(8192, buf, '\t', ' ', ' ', ' ');
    StrTrim(buf, '\0');
    if (nb > 8) {
      if (strstr(buf, "# (eV/u)")) {
	n = StrSplit(buf, ' ');
	ncx = n-3;
	s = buf;
	for (i = 0; i < n-1; i++) {
	  ns = strlen(s);
	  if (i < 2) {
	    s += ns+1;
	    continue;
	  }
	  StrTrim(s, '\0');
	  if (strstr(s, "n=")) {
	    nres = 1;
	    nn = atoi(s+2);
	    if (nmax < nn) nmax = nn;
	    na[i] = nn;
	    ka[i] = -1;
	    sa[i] = -1;
	  } else {
	    nn = atoi(s);
	    if (nmax < nn) nmax = nn;
	    na[i] = nn;
	    char *sp = strstr(s, "[");
	    if (sp) {
	      ka[i] = atoi(sp+1);
	      sp = strstr(s, "]");
	      sa[i] = atoi(sp+2);
	    } else {
	      nr = strlen(s);
	      ka[i] = LFromSym(s[nr-5]);
	      sa[i] = atoi(s+nr-3);
	    }
	  }
	  s += ns+1;
	}
      }
    }
    if (buf[0] != '#') {
      n = StrSplit(buf, ' ');
      if (n == ncx+2) {
	nep++;
      }
    }
  }
  fseek(f, 0, SEEK_SET);
  printf("kronos file: %s\n", kfn);
  printf("kronos data: %d %d %s %s %g %g %d %d %d\n",
	 z, k, tgt0, cxm0, tmass, pmass, nmax, nep, ncx);
  cx = InitKronos(z, k, nmax, nep, cxm0, tgt0);
  cx->tmass = tmass * AMU;
  cx->pmass = pmass * AMU;
  cx->rmass = rmass * AMU;
  int j = 0, k1, irx = 0, i1;
  double a;
  while(1) {
    if (NULL == fgets(buf, 8192, f)) break;
    nb = strlen(buf);
    if (nb < 8) continue;
    irx = 0;
    if (buf[0] == '#') {
      if (buf[1] == ' ' && buf[2] == 'R' && buf[3] == 'X') {
	irx = 1;
      } else {
	irx = 0;
	continue;
      }
    }
    StrReplace(8192, buf, '\t', ' ', ' ', ' ');
    StrTrim(buf, '\0');
    n = StrSplit(buf, ' ');
    if (n != ncx+2) continue;
    s = buf;
    for (i = 0; i < n-1; i++) {
      ns = strlen(s);
      if (irx) {
	if (i >= 2) {
	  a = atof(s);
	  nn = cx->idn[na[i]-1];
	  cx->rx[nn+ka[i]] = a;
	}
	s += ns+1;
	continue;
      }
      a = atof(s);
      if (i == 0) {
	cx->ep[j] = a;
	vr = sqrt((cx->ep[j]/HARTREE_EV)*2/AMU);
      } else {
	i1 = i+1;
	nn = cx->idn[na[i1]-1];
	if (nres == 0) {
	  if (sa[i1] != 3) {
	    cx->cx0[nn+ka[i1]][j] = a;
	  } else {
	    cx->cx1[nn+ka[i1]][j] = a;
	  }
	} else {
	  for (k1 = 0; k1 < na[i1]; k1++) {
	    wkn = LandauZenerLD(cx->lnfac, na[i1], k1, c,
				vr*cx->rx[nn+k1], md, _cxldistmj);
	    cx->cx0[nn+k1][j] = a*wkn;
	  }
	}
      }
      s += ns+1;
    }
    if (!irx) j++;
  }
  fclose(f);
  if (ilog < 0) ilog = 3;
  cx->ilog = ilog;
  if (ilog & 1) {
    for (j = 0; j < cx->nep; j++) {
      cx->ep[j] = log(cx->ep[j]);
    }
  }
  double *xx, *yy;
  xx = malloc(sizeof(double)*cx->nep);
  yy = malloc(sizeof(double)*cx->nep);
  for (i = 0; i < cx->ncx; i++) {
    for (j = 0; j < cx->nep; j++) {
      if (cx->cx0[i][j] > 0) break;
    }
    nb = j;
    for (j = cx->nep-1; j >= nb; j--) {
      if (cx->cx0[i][j] > 0) break;
    }
    ns = j;
    if (ns == nb) {
      cx->cx0[i][ns] = 0;
    } else if (ns > nb) {
      n = 0;
      for (j = nb; j <= ns; j++) {
	if (cx->cx0[i][j] > 0) {
	  xx[n] = cx->ep[j];
	  yy[n] = cx->cx0[i][j];
	  n++;
	}
      }
      for (j = nb; j <= ns; j++) {
	if (!cx->cx0[i][j]) {
	  nn = 3;
	  nr = 1;
	  UVIP3P(nn, n, xx, yy, nr, &cx->ep[j], &cx->cx0[i][j]);
	}
      }
    }
    if (cx->cx1) {
      for (j = 0; j < cx->nep; j++) {
	if (cx->cx1[i][j] > 0) break;
      }
      nb = j;
      for (j = cx->nep-1; j >= nb; j--) {
	if (cx->cx1[i][j] > 0) break;
      }
      ns = j;
      if (ns == nb) {
	cx->cx1[i][ns] = 0;
      } else if (ns > nb) {
	n = 0;
	for (j = nb; j <= ns; j++) {
	  if (cx->cx1[i][j] > 0) {
	    xx[n] = cx->ep[j];
	    yy[n] = cx->cx1[i][j];
	    n++;
	  }
	}
	for (j = nb; j <= ns; j++) {
	if (!cx->cx1[i][j]) {
	  nn = 3;
	  nr = 1;
	  UVIP3P(nn, n, xx, yy, nr, &cx->ep[j], &cx->cx1[i][j]);
	}
	}
      }
    }
  }
  free(xx);
  free(yy);
  if (ilog & 2) {
    for (i = 0; i < cx->ncx; i++) {
      for (j = 0; j < cx->nep; j++) {
	if (cx->cx0[i][j] > 0) {
	  cx->cx0[i][j] = log(cx->cx0[i][j]);
	} else {
	  cx->cx0[i][j] = -500.0;
	}
	if (cx->cx1) {
	  if (cx->cx1[i][j] > 0) {
	    cx->cx1[i][j] = log(cx->cx1[i][j]);
	  } else {
	    cx->cx1[i][j] = -500.0;
	  }
	}
      }
    }  
  }
#undef ncxm
  return 0;
}

void SetOptionRates(char *s, char *sp, int ip, double dp) {
  if (0 == strcmp(s, "rates:cxldistmj")) {
    _cxldistmj = dp;
  }
}
