#include "rates.h"

static char *rcsid="$Id: rates.c,v 1.8 2002/02/04 15:48:33 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

static int iedist = 0;
static DISTRIBUTION ele_dist[MAX_DIST];
static int ipdist = 0;
static DISTRIBUTION pho_dist[MAX_DIST];

#define QUAD_LIMIT 64
static int _iwork[QUAD_LIMIT];
static double _dwork[4*QUAD_LIMIT];

#define RT_CE 1
#define RT_CI 2
#define RT_RR 3
static struct {
  DISTRIBUTION *d;
  double (*Rate1E)(double, double, int, int, void *);
  double eth;
  int np;
  int nshells;
  void *params;
  double epsabs;
  double epsrel;
  int i, f;
  int type;
} rate_args;

void dqagi_(double (*f)(double *), 
	    double *bound, int *inf, double *epsabs,
	    double *epsrel, double *result, double *abserr,
	    int *neval, int *ier, int *limit, int *lenw,
	    int *last, int *iwork, double *work);

void uvip3p_(int *np, int *ndp, double *x, double *y, 
	     int *n, double *xi, double *yi);

void ionis_(int *iz, int *ic, double *t, 
	    double *a, double *dir, double *rion);
void recomb_(int *iz, int *ic, double *t, double *rr, double *dr);
void recombfe_(int *iz, int *ic, double *t, double *rr, double *dr);


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
   

int SetEleDist(int i, int np, double *p) {
  int k;

  if (ele_dist[i].dist == NULL) {
    printf("Electron Dist. does not exist\n");
    return -1;
  }
  if (ele_dist[i].nparams != np) {
    printf("Num of Params for Electron Dist. %d does not match\n", i);
    return -1;
  }
  iedist = i;
  for (k = 0; k < np; k++) {
    ele_dist[i].params[k] = p[k];
  }

  return 0;
}

int SetPhoDist(int i, int np, double *p) {
  int k;

  if (pho_dist[i].dist == NULL) {
    printf("Photon Dist. does not exist\n");
    return -1;
  }

  if (pho_dist[i].nparams != np) {
    printf("Num of Params for Photon Dist. %d does not match\n", i);
    return -1;
  }
  ipdist = i;
  for (k = 0; k < np; k++) {
    pho_dist[i].params[k] = p[k];
  }

  return 0;
}

DISTRIBUTION *GetEleDist(int *i) {
  if (i) *i = iedist;
  return ele_dist+iedist;
}

DISTRIBUTION *GetPhoDist(int *i) {
  if (i) *i = ipdist;
  return pho_dist+iedist;
}

int SetRateAccuracy(double epsrel, double epsabs) {
  if (epsrel > 0.0) rate_args.epsrel = epsrel;
  if (epsabs > 0.0) rate_args.epsabs = epsabs;
  return 0.0;
}

static double RateIntegrand(double *e) {
  double a, b;

  a = rate_args.d->dist(*e, rate_args.d->params);
  b = rate_args.Rate1E(*e, rate_args.eth, rate_args.np,
		       rate_args.nshells, rate_args.params);
  return a*b;
}
  
double IntegrateRate(int idist, double eth, double bound, 
		     int np, int nshells, void *params, 
		     int i0, int f0, int type, 
		     double (*Rate1E)(double, double, int, int, void *)) { 
  double result;
  int neval, inf, ier, limit, lenw, last;
  double epsabs, epsrel, abserr;

  inf = 1;
  ier = 0;
  epsabs = rate_args.epsabs;
  epsrel = rate_args.epsrel;
  limit = QUAD_LIMIT;
  lenw = 4*limit;
  
  rate_args.Rate1E = Rate1E;
  if (idist == 0) rate_args.d = ele_dist + iedist;
  else rate_args.d = pho_dist + ipdist;
  rate_args.eth = eth;
  rate_args.np = np;
  rate_args.nshells = nshells;
  rate_args.params = params;
  rate_args.i = i0;
  rate_args.f = f0;
  rate_args.type = type;

  if (bound == 0.0) bound = eth*EPS3;
  dqagi_(RateIntegrand, &bound, &inf, &epsabs, &epsrel,
	 &result, &abserr, &neval, &ier, &limit, &lenw,
	 &last, _iwork, _dwork);

  if (ier == 2 || ier == 4) {
    if (abserr < epsabs) return result;
    if (abserr < fabs(epsrel*result)) return Max(0.0, result);
  }

  if (ier != 0) {
    printf("IntegrateRate Error: %d %d %10.3E %10.3E %10.3E\n", 
	   ier, neval, bound, result, abserr);
    printf("%6d %6d %2d Eth = %10.3E\n", 
	   rate_args.i, rate_args.f, type, eth);
  }

  if (result < 0.0) result = 0.0;
  return result;
}

double IntegrateRate2(int idist, double e, int np, 
		      int nshells, void *params,
		      int i0, int f0, int type,
		      double (*Rate1E)(double, double, 
				       double, int, int, void *)) { 
  double a;
  
  a = Rate1E(e, e, e, np, nshells, params);
  return a;
}

double CERate1E(double e, double eth, int np, int ns, void *p) {
  double *x, *y;
  int m1, m2, n, one;
  double *dp, a, x0;

  if (e < eth) return 0.0;
  m1 = np + 1;
  m2 = m1 + m1;
  dp = (double *) p;
  y = dp+1;
  x0 = eth/e;
  if (dp[0] >= 0) {
    x = y + m1;
  } else {
    x = y + m2;
    x0 *= x0;
  }

  n = 3;
  one = 1;
  uvip3p_(&n, &m1, x, y, &one, &x0, &a);
  if (dp[0] > 0.0) {
    a -= dp[0]*log(x0);
  }
  if (a <= 0.0) {
    a = 0.0;
    return a;
  }
  
  a *= PI*AREA_AU20*HARTREE_EV/(2.0*e);
  a *= VelocityFromE(e);
  return a;
}

double DERate1E(double e, double eth, int np, int ns, void *p) {
  double a, x0, *x, *y;
  double *dp;
  int m1, m2, n, one;

  x0 = eth/(eth+e);
  m1 = np + 1;
  m2 = m1 + m1;
  dp = (double *) p;
  y = dp+1;
  if (dp[0] >= 0) {
    x = y + m1;
  } else {
    x = y + m2;
    x0 *= x0;
  }

  n = 3;
  one = 1;
  uvip3p_(&n, &m1, x, y, &one, &x0, &a);
  if (dp[0] > 0.0) {
    a += dp[0]*log(e);
  }
  
  if (a <= 0.0) {
    a = 0.0;
    return a;
  }
  
  a *= PI*AREA_AU20*HARTREE_EV/(2.0*e);
  a *= VelocityFromE(e);
  return a;
}
  
int CERate(double *dir, double *inv, int iinv, 
	   int j1, int j2, double e,
	   int m, double *params, int i0, int f0) {
  double a, e0;

  e0 = e*HARTREE_EV;	     
  a = IntegrateRate(0, e0, e0, m, 1, params, 
		    i0, f0, RT_CE, CERate1E);
  *dir = a/(j1 + 1.0);
  if (iinv) {
    if (iedist == 0) {
      a *= exp(e0/ele_dist[0].params[0]);
      *inv = a/(j2 + 1.0);
    } else {
      *inv = IntegrateRate(0, e0, 0.0, m, 1, params, 
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
  if (iinv) {
    e0 = e*HARTREE_EV;
    b = pho_dist[ipdist].dist(e0, pho_dist[ipdist].params);
    if (b > EPS10) {
      *inv = factor*b*a/(e0*e0*e0*(j2+1.0));
    } else {
      *inv = 0.0;
    }
  } else {
    *inv = 0.0;
  }
  return 0;
}

double CIRate1E(double e, double eth, int np, int ns, void *p) {
  int i;
  float *dp;
  double x;
  double a, b, c, f;
  double c2, logc, logx;

  if (np != 6) return 0.0;

  dp = (float *) p;
  x = e/eth;
  c = 2.0/(1.0+x);
  c2 = sqrt(c);
  logc = log(c);
  logx = log(x);
  f = 0.0;
  for (i = 0; i < ns; i++) {
    a = (3.5 + dp[4])*logc;
    b = dp[1]*log((1.0 + dp[2])/(c2 + dp[2]));
    c = dp[0]*(1.0 - exp(a+b))*logx;
    c += dp[3]*(1.0 - 1.0/x - logx/(1.0+x));
    f += c*x/(x+1.0);
    dp += np;
  }

  c = AREA_AU20*HARTREE_EV*f/(2.0*e);
  c *= VelocityFromE(e);

  return c;
}

double R3BRate1E(double e1, double e2, double eth, 
		 int np, int ns, void *p) {
  return 0.0;
}

int CIRate(double *dir, double *inv, int iinv, 
	   int j1, int j2, double e,
	   int m, int nshells, float *params, int i0, int f0) {
  double e0;

  e0 = e*HARTREE_EV;
  *dir = IntegrateRate(0, e0, e0, m, nshells, params, 
		       i0, f0, RT_CI, CIRate1E);
  *dir /= (j1 + 1.0);
  if (iinv) {
    *inv = IntegrateRate2(0, e0, m, nshells, params, 
			  i0, f0, -RT_CI, R3BRate1E);
    *inv /= (j2 + 1.0); 
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
    uvip3p_(&np, &i, log_xseaton, seaton[0], &one, &logx, &s0);
    uvip3p_(&np, &i, log_xseaton, seaton[1], &one, &logx, &s1);
    uvip3p_(&np, &i, log_xseaton, seaton[2], &one, &logx, &s2);
    if (top) {
      uvip3p_(&np, &i, log_xseaton, seaton[3], &one, &logx, &g0);
      uvip3p_(&np, &i, log_xseaton, seaton[4], &one, &logx, &g1);
      uvip3p_(&np, &i, log_xseaton, seaton[5], &one, &logx, &g2);
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
  
double RRRate1E(double e, double eth, int np, int ns, void *p) {
  int i;
  float *dp;
  double x, eth0;
  double a, b, c, f;
  double x2, logx;
  
  if (np != 5) return 0.0;
  dp = (float *) p;

  f = 0.0;
  for (i = 0; i < ns; i++) {  
    eth0 = HARTREE_EV * dp[4];
    x = 1.0 + e/eth0;
    x2 = sqrt(x);
    logx = log(x);
    a = -4.5-dp[3] + 0.5*dp[1];
    b = (1.0 + dp[2])/(x2 + dp[2]);
    c = dp[1]*log(b) + a*logx;
    c = dp[0]*exp(c);
    f += c*(e+eth)/eth0;
    dp += np;
  }
  c = 2.0*PI*FINE_STRUCTURE_CONST*f*AREA_AU20;
  a = FINE_STRUCTURE_CONST*(e+eth);
  a = a*a;
  a = a/(2.0*HARTREE_EV*e);
  c *= a*VelocityFromE(e);

  return c;
}

double PIRate1E(double e, double eth, int np, int ns, void *p) {
  int i;
  float *dp;
  double x, eth0;
  double a, b, c, f;
  double x2, logx;
  const double factor = 1.871156686E2;
  
  if (np != 4) return 0.0;

  dp = (float *) p;
  f = 0.0;
  for (i = 0; i < ns; i++) {
    eth0 = HARTREE_EV * dp[4];
    x = (e-eth+eth0)/eth0;
    x2 = sqrt(x);
    logx = log(x);
    a = -4.5-dp[3] + 0.5*dp[1];
    b = (1.0 + dp[2])/(x2 + dp[2]);
    c = dp[1]*log(b) + a*logx;
    c = dp[0]*exp(c);
    f += c*e/eth0;
    dp += np;
  }
  
  c = 2.0*PI*FINE_STRUCTURE_CONST*f*AREA_AU20;
  c *= factor/e;

  return c;
}

int RRRate(double *dir, double *inv, int iinv, 
	   int j1, int j2, double e,
	   int m, int nshells, float *params, int i0, int f0) {
  double e0;

  e0 = e*HARTREE_EV;	
  *dir = IntegrateRate(0, e0, 0.0, m, nshells, params, 
		       i0, f0, RT_RR, RRRate1E);
  *dir /= (j1 + 1.0);
  if (iinv) {
    *inv = IntegrateRate(1, e0, e0, m, nshells, params, 
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
      a *= VelocityFromE(e0);
    }
    *inv = a;
  } else {
    *inv = 0.0;
  }
  return 0;
}

static double Maxwell(double e, double *p) {
  double x;
  const double maxwell_const = 1.12837967;

  x = e/p[0];
  if (x < 0.0 || x > 35) return 0.0;

  x = maxwell_const * sqrt(x) * exp(-x)/p[0];
  return x;
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

double Ionis(int z, int nele, double t, double *a, double *dir) {
  double total, aa, dd;
  
  nele = z+1 - nele;
  ionis_(&z, &nele, &t, &aa, &dd, &total);
  aa *= 1E10;
  dd *= 1E10;
  total *= 1E10;
  if (a) *a = aa;
  if (dir) *dir = dd;

  return total;
}

double Recomb(int z, int nele, double t, double *rr, double *dr) {
  double total, r, d;
  
  nele = z+1 - nele;
  t = t/8.617385E-5;
  if (z == 26) {
    recombfe_(&z, &nele, &t, &r, &d);
  } else {
    recomb_(&z, &nele, &t, &r, &d);
  }
  r *= 1E10;
  d *= 1E10;
  if (rr) *rr = r;
  if (dr) *dr = d;
  total = r + d;

  return total;
}

int FracAbund(int z, double t, double *a) {
  int nele;
  double *ir, *rr, c;

  ir = (double *) malloc(sizeof(double)*(z+1));
  rr = (double *) malloc(sizeof(double)*(z+1));
  
  for (nele = 0; nele <= z; nele++) {
    if (nele != 0) {
      ir[nele] = Ionis(z, nele, t, NULL, NULL);
    }
    if (nele != z) {
      rr[nele] = Recomb(z, nele, t, NULL, NULL);
    }
  }

  a[0] = 1.0;
  c = 1.0;
  for (nele = 1; nele <= z; nele++) {
    a[nele] = a[nele-1]*(rr[nele-1]/ir[nele]);
    c += a[nele];
  }
  for (nele = 0; nele <= z; nele++) {
    a[nele] /= c;
  }
  
  free(ir);
  free(rr);
  return 0;
}

double MaxAbund(int z, int nele, double *a, double eps) {
  double amax, da;
  double t1, t2, t;
  double dt;
  
  t1 = 0.1*HARTREE_EV;
  t2 = z*z*HARTREE_EV;
  if (nele == 0) {
    FracAbund(z, t2, a);
    return t2;
  }
  if (nele == z) {
    FracAbund(z, t1, a);
    return t1;
  }
  amax = 0.0;
  t = 0.5*(t1+t2);
  dt = t*eps;
  while (t2 - t1 > dt) {
    FracAbund(z, t, a);
    amax = a[nele];
    FracAbund(z, t+dt, a);
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
  
int InitRates(void) {
  int i;

  iedist = 0;
  ipdist = 0;
  
  i = 0;
  ele_dist[i].nparams = 1;
  ele_dist[i].params = (double *) malloc(sizeof(double));
  ele_dist[i].params[0] = 1.0E3;
  ele_dist[i].dist = Maxwell;
  i++;
  
  for (; i < MAX_DIST; i++) {
    ele_dist[i].nparams = 0;
    ele_dist[i].params = NULL;
    ele_dist[i].dist = NULL;
  }

  i = 0; 
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

  rate_args.epsabs = EPS10;
  rate_args.epsrel = EPS2;

  for (i = 0; i < NSEATON; i++) {
    log_xseaton[i] = log(xseaton[i]);
  }
  return 0;
}
