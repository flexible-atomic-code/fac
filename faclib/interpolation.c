#include "interpolation.h"

static char *rcsid="$Id: interpolation.c,v 1.6 2001/10/24 23:31:36 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

static struct {
  double *x;
  double *logx;
  double *y;
  double *sigma;
  void *extra;
  void (*function)(int, double *, int , double *, double *, 
		   double *, double *, int, void *);
} minpack_params;

void dgesdd_(char *jobz, int *m, int *n, double *a, int *lda, 
	     double *s, double *u, int *ldu, double *vt, int *ldvt,
	     double *work, int *lwork, int *iwork, int *info);
void chkder_(int *, int *, double *, double *, double *, int *, 
	     double *, double *, int *, double *);
void lmder_(void F(int *, int *, double *, double *, double *, int *, int *),
	    int *, int *, double *, double *, double *, int *, double *, 
	    double *, double *, int *, double *, int *, double *, int *,
	    int *, int *, int *, int *, double *, double *, double *,
	    double *, double *);
void lmder1_(void F(int *, int *, double *, double *, double *, int *, int *),
	     int *, int *, double *, double *, double *, int *, double *, 
	     int *, int *, double *, int *);

void spline(double *x, double *y, int n, 
	    double yp1, double ypn, double *y2) {
  int i, k;
  double p, qn, sig, un, *u;
  double a, b, c;

  u = malloc(sizeof(double)*n);
  
  if (n == 2) {
    y2[0] = 0.0;
    y2[n-1] = 0.0;
    return;
  } 

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
    Basis(np, afunc, x[i], logx[i]);
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

  dgesdd_(jobz, &nd, &np, u, &nd, w, u, &nd, v, &np, 
	  dwork, &lwork, iwork, &infor);
    
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
      Basis(np, afunc, x[i], logx[i]);
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
  chkder_(&n, &np, p, fvec, fjac, &ldfjac, wa3, wa4, &mode, diag);
  MinFunc(&n, &np, wa3, wa4, fjac, &ldfjac, &mode);
  mode = 2;
  MinFunc(&n, &np, p, fvec, fjac, &ldfjac, &mode);
  chkder_(&n, &np, p, fvec, fjac, &ldfjac, wa3, wa4, &mode, diag);
  for (nfev = 0; nfev < n; nfev++) {
    printf("%d %10.3E\n", nfev, diag[nfev]);
  }
  exit(1);
  */
  maxfev = 1000*np;
  mode = 1;
  nprint = 0;
  factor = 100.0;
  lmder_(MinFunc, &n, &np, p, fvec, fjac, &ldfjac, 
	 &tol, &tol, &zero, &maxfev, diag, &mode, &factor,
	 &nprint, &info, &nfev, &njev, ipvt, qtf, wa1, wa2, wa3, wa4);

  return info;
}

