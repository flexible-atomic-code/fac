#include "rates.h"

static char *rcsid="$Id: rates.c,v 1.5 2002/01/20 06:02:56 mfgu Exp $";
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

  if (np != 5) return 0.0;

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
    f += dp[0]*(1.0 - exp(a+b));
    f += dp[3]*(1.0 - 1.0/x - logx/(1.0+x));
    p += np;
  }

  c = AREA_AU20*HARTREE_EV*f/(2.0*(e+eth));
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

double RRRate1E(double e, double eth, int np, int ns, void *p) {
  int i;
  float *dp;
  double x;
  double a, b, c, f;
  double x2, logx;
  
  if (np != 4) return 0.0;
  dp = (float *) p;

  x = 1.0 + e/eth;
  x2 = sqrt(x);
  logx = log(x);
  f = 0.0;
  for (i = 0; i < ns; i++) {  
    a = -3.5-dp[3] + 0.5*dp[1];
    b = (1.0 + dp[2])/(x2 + dp[2]);
    c = dp[1]*log(b) + a*logx;
    c = dp[0]*exp(c);
    f += c;
    p += np;
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
  double x;
  double a, b, c, f;
  double x2, logx;
  const double factor = 1.871156686E2;
  
  if (np != 4) return 0.0;

  dp = (float *) p;
  x = e/eth;
  x2 = sqrt(x);
  logx = log(x);
  f = 0.0;
  for (i = 0; i < ns; i++) {
    a = -3.5-dp[3] + 0.5*dp[1];
    b = (1.0 + dp[2])/(x2 + dp[2]);
    c = dp[1]*log(b) + a*logx;
    c = dp[0]*exp(c);
    f += c;
    p += np;
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

  return 0;
}
