#include "interpolation.h"
#include "cf77.h"

static char *rcsid="$Id: interpolation.c,v 1.14 2004/07/06 07:09:25 mfgu Exp $";
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

static struct {
  double *x;
  double *logx;
  double *y;
  double *sigma;
  void *extra;
  void (*function)(int, double *, int , double *, double *, 
		   double *, double *, int, void *);
} minpack_params;

void spline(double *x, double *y, int n, 
	    double yp1, double ypn, double *y2) {
  int i, k;
  double p, qn, sig, un, *u;

  if (n == 2) {
    y2[0] = 0.0;
    y2[n-1] = 0.0;
    return;
  } 

  u = malloc(sizeof(double)*n);

  if (yp1 > 0.99E30) {
    y2[0] = 0.0;
    u[0] = 0.0;
  } else {
    y2[0] = -0.5;
    u[0] = (3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
  }

  if (ypn > 0.99E30) {
    qn = 0.0;
    un = 0.0;
  } else {
    qn = 0.5;
    un = (3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
  }

  for (i = 1; i < n-1; i++) {
    sig = (x[i] - x[i-1])/(x[i+1]-x[i-1]);
    p = sig*y2[i-1]+2.0;
    y2[i] = (sig-1.0)/p;
    u[i] = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i] = (6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }

  y2[n-1] = (un-qn*u[n-2])/(qn*y2[n-2]+1.0);
  for (k = n-2; k >= 0; k--) {
    y2[k] = y2[k]*y2[k+1]+u[k];
  }

  free(u);

}
  
int splint(double *xa, double *ya, double *y2a, 
	   int n, double x, double *y) {
  int klo, khi, k;
  double h, b, a;
  
  *y = 0.0;
  klo = 0;
  khi = n-1;
  while (khi-klo > 1) {
    k = (khi+klo) >> 1;
    if (xa[k] > x) khi = k;
    else klo = k;
  }
  h = xa[khi] - xa[klo];
  if (h == 0.0) return -1;
  a = (xa[khi] - x) / h;
  b = (x - xa[klo]) / h;
  *y = a*ya[klo] + b*ya[khi] + ((a*a*a-a)*y2a[klo] +
				(b*b*b-b)*y2a[khi]) * (h*h)/6.0;
  return 0;
}

void splie2(double *x1a, double *x2a, double **ya, 
	    int m, int n, double **y2a) {
  int j;

  for (j = 0; j < m; j++) {
    spline(x2a, ya[j], n, 1.0E30, 1.0E30, y2a[j]);
  }
}

void splin2(double *x1a, double *x2a, double **ya, double **y2a,
	    int m, int n, double x1, double x2, double *y) {
  int j;
  double *ytmp, *yytmp;

  ytmp = malloc(sizeof(double)*m);
  yytmp = malloc(sizeof(double)*m);
  for (j = 0; j < m; j++) {
    splint(x2a, ya[j], y2a[j], n, x2, &yytmp[j]);
  }
  spline(x1a, yytmp, m, 1.0E30, 1.0E30, ytmp);
  splint(x1a, yytmp, ytmp, m, x1, y);
  free(ytmp);
  free(yytmp);
}

void PolyBasis(int n, double *y, double x, double logx) {
  int i;
  double t;
  
  t = 1.0;
  for (i = 0; i < n; i++) {
    y[i] = t;
    t *= x;
  }
}

void PolyFit(int n, double *c, int nd, double *x, double *y) {
  SVDFit(n, c, NULL, EPS6, nd, x, NULL, y, NULL, PolyBasis);
}

void SVDFit(int np, double *coeff, double *chisq, double tol,
	    int nd, double *x, double *logx, double *y, double *sig,
	    void Basis(int, double *, double, double)) {
  double *u, *v, *w, *b, *afunc;
  int i, j, k;
  double tmp, thresh, wmax, sum;
  double *dwork;
  int *iwork, lwork, infor;
  char jobz[] = "O";
 
  k = np*np;
  lwork = 5*k + 4*np;
  if (lwork < nd) lwork = nd;
  lwork += 3*k;
  lwork *= 2;
  i = (np+nd)*(np+1)+np;
  w = (double *) malloc(sizeof(double)*(i+lwork));
  iwork = (int *) malloc(sizeof(int)*8*np);

  v = w + np;
  u = v + np*np;
  b = u + np*nd;
  afunc = b + nd;
  dwork = afunc + np;
  for (i = 0; i < nd; i++) {
    if (logx) {
      Basis(np, afunc, x[i], logx[i]);
    } else {
      Basis(np, afunc, x[i], 0.0);
    }
    k = i;
    if (sig) {
      tmp = 1.0/sig[i];
      for (j = 0; j < np; j++) {
	u[k] = afunc[j]*tmp;
	k += nd;
      }
      b[i] = y[i]*tmp;
    } else {
      for (j = 0; j < np; j++) {
	u[k] = afunc[j];
	k += nd;
      }
      b[i] = y[i];
    }
  }

  DGESDD(jobz, nd, np, u, nd, w, u, nd, v, np,
	  dwork, lwork, iwork, &infor);
    
  wmax = w[0];
  thresh = tol*wmax;
  for (j = 0; j < np; j++) {
    if (w[j] < thresh) w[j] = 0.0;
  }
  
  k = 0;
  for (j = 0; j < np; j++) {
    tmp = 0.0;
    if (w[j] != 0.0) {
      for (i = 0; i < nd; i++) {
	tmp += u[k++]*b[i];
      }
      tmp /= w[j];
    }
    afunc[j] = tmp;
  }
  
  k = 0;
  for (j = 0; j < np; j++) {
    tmp = 0.0;
    for (i = 0; i < np; i++) {
      tmp += v[k++]*afunc[i];
    }
    coeff[j] = tmp;
  }
  
  if (chisq) {
    *chisq = 0.0;
    for (i = 0; i < nd; i++) {
      if (logx) {
	Basis(np, afunc, x[i], logx[i]);
      } else {
	Basis(np, afunc, x[i], 0);
      }
      sum = 0.0;
      for (j = 0; j < np; j++) {
	sum += coeff[j]*afunc[j];
      }
      tmp = y[i] - sum;
      if (sig) tmp /= sig[i];
      *chisq += tmp*tmp;
    }
  }  
  
  free(w);
  free(iwork);
}

static void MinFunc(int *m, int *n, double *x, double *fvec, 
		    double *fjac, int *ldfjac, int *iflag) {
  int ndy, i, j, k;

  if (*iflag == 1) {
    ndy = 0;
    minpack_params.function(*n, x, *m, minpack_params.x,
			    minpack_params.logx, fvec, 
			    NULL, ndy, minpack_params.extra);
    for (i = 0; i < *m; i++) {
      fvec[i] = fvec[i] - minpack_params.y[i];
      if (minpack_params.sigma) fvec[i] = fvec[i]/minpack_params.sigma[i];
    }
  } else if (*iflag == 2){
    ndy = *ldfjac;
    minpack_params.function(*n, x, *m, minpack_params.x,
			    minpack_params.logx, NULL, 
			    fjac, ndy, minpack_params.extra);
    if (minpack_params.sigma) {
      for (i = 0; i < *m; i++) {
	k = i;
	for (j = 0; j < *n; j++) {
	  fjac[k] = fjac[k]/minpack_params.sigma[i];
	  k += ndy;
	}
      }
    }
  }
}
/* provide fortran access with cfortran.h */
FCALLSCSUB7(MinFunc, MINFUNC, minfunc, PINT, PINT, 
	    DOUBLEV, DOUBLEV, DOUBLEV, PINT, PINT)
    
int NLSQFit(int np, double *p, double tol, int *ipvt,
	    double *fvec, double *fjac, int ldfjac, double *wa, int lwa,
	    int n, double *x, double *logx, double *y, double *sig,
	    void Func(int, double *, int , double *, double *, 
		      double *, double *, int, void *), 
	    void *extra) {
  int info, nprint, nfev, njev;
  int maxfev, mode;
  double factor, *diag, *qtf, *wa1, *wa2, *wa3, *wa4;
  double zero = 0.0;
  double ftol;

  minpack_params.x = x;
  minpack_params.logx = logx;
  minpack_params.y = y;
  minpack_params.sigma = sig;
  minpack_params.extra = extra;
  minpack_params.function = Func;
  
  diag = wa;
  qtf = wa+np;
  wa1 = qtf+np;
  wa2 = wa1+np;
  wa3 = wa2+np;
  wa4 = wa3+np;
  mode = 1;
  /*
  MinFunc(&n, &np, p, fvec, fjac, &ldfjac, &mode);
  CHKDER(n, np, p, fvec, fjac, ldfjac, wa3, wa4, mode, diag);
  MinFunc(&n, &np, wa3, wa4, fjac, &ldfjac, &mode);
  mode = 2;
  MinFunc(&n, &np, p, fvec, fjac, &ldfjac, &mode);
  CHKDER(n, np, p, fvec, fjac, ldfjac, wa3, wa4, mode, diag);
  for (nfev = 0; nfev < n; nfev++) {
    printf("%d %10.3E\n", nfev, diag[nfev]);
  }
  exit(1);
  */
  maxfev = 5000*np;
  mode = 1;
  nprint = 0;
  factor = 100.0;
  ftol = tol*tol*n;
  LMDER(C_FUNCTION(MINFUNC, minfunc), n, np, p, fvec, fjac, ldfjac, 
	ftol, tol, zero, maxfev, diag, mode, factor,
	nprint, &info, &nfev, &njev, ipvt, qtf, wa1, wa2, wa3, wa4);

  return info;
}


double Simpson(double *y, int ia, int ib) {
  int i;
  double a;

  a = 0.0;
  
  for (i = ia; i < ib - 1; i += 2) {
    a += y[i] + 4.0*y[i+1] + y[i+2];
  }
  a /= 3.0;
  if (i < ib) a += 0.5 * (y[i] + y[ib]);

  return a;
}

/* integration by newton-cotes formula */
int NewtonCotes(double *r, double *x, int i0, int i1, int m, int maxrp) {
  int i, j, n, k;
  double a;

  for (i = i0; i <= i1-4; i += 4) {
    a = 0.0;
    for (j = 0, k = i; j <= 4; j++, k++) {
      a += _CNC[4][j] * x[k];
    }
    r[i+4] = r[i] + a;
  }
  if (i1 < maxrp-1) {
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
