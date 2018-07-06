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

#include "interpolation.h"
#include "cf77.h"

static char *rcsid="$Id$";
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

static double gauss_xw[2][15] = {
  {.933078120172818E-01,
   .492691740301883E+00,
   .121559541207095E+01,
   .226994952620374E+01,
   .366762272175144E+01,
   .542533662741355E+01,
   .756591622661307E+01,
   .101202285680191E+02,
   .131302824821757E+02,
   .166544077083300E+02,
   .207764788994488E+02,
   .256238942267288E+02,
   .314075191697539E+02,
   .385306833064860E+02,
   .480260855726858E+02},
  {.218234885940086E+00,
   .342210177922884E+00,
   .263027577941681E+00,
   .126425818105931E+00,
   .402068649210010E-01,
   .856387780361184E-02,
   .121243614721425E-02,
   .111674392344251E-03,
   .645992676202287E-05,
   .222631690709627E-06,
   .422743038497938E-08,
   .392189726704110E-10,
   .145651526407313E-12,
   .148302705111330E-15,
   .160059490621113E-19}
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


void spline_work(double *x, double *y, int n, 
		 double yp1, double ypn, double *y2, double *work) {
  int i, k;
  double p, qn, sig, un, *u;

  if (n == 2) {
    y2[0] = 0.0;
    y2[n-1] = 0.0;
    return;
  } 

  u = work;

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
}

void spline(double *x, double *y, int n, 
	    double yp1, double ypn, double *y2) {
  double *u;

  if (n == 2) {
    y2[0] = 0.0;
    y2[n-1] = 0.0;
    return;
  } 

  u = malloc(sizeof(double)*n);
  spline_work(x, y, n, yp1, ypn, y2, u);

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
  double *u;

  u = malloc(sizeof(double)*n);
  for (j = 0; j < m; j++) {
    spline_work(x2a, ya[j], n, 1.0E30, 1.0E30, y2a[j], u);
  }
  free(u);
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


double Simpson(double *x, int i0, int i1) {
  int i, k;
  double a, b;

  b = x[i0];
  a = 0.0;
  for (i = i0+1; i < i1; i += 2) {
    a += x[i];
  }
  b += 4.0*a;
  a = 0.0;
  k = i1-1;
  for (i = i0+2; i < k; i += 2) {
    a += x[i];
  }
  b += 2.0*a;
  if (i == i1) {
    b += x[i1];
    b /= 3.0;
  } else {
    b += x[k];
    b /= 3.0;
    b += 0.5*(x[k] + x[i1]);
  }

  return b;
}

/* integration by newton-cotes formula */
int NewtonCotes(double *r, double *x, int i0, int i1, int m, int id) {
  int i, k;
  double a, yp;

  if (id >= 0) {
    if (m >= 0) {
      r[i1] = x[i0];
      a = 0.0;
      for (i = i0+1; i < i1; i += 2) {
	a += x[i];
      }
      r[i1] += 4.0*a;
      a = 0.0;
      k = i1-1;
      for (i = i0+2; i < k; i += 2) {
	a += x[i];
      }
      r[i1] += 2.0*a;
      if (i == i1) {
	r[i1] += x[i1];
	r[i1] /= 3.0;
      } else {
	r[i1] += x[k];
	r[i1] /= 3.0;
	r[i1] += 0.5*(x[k] + x[i1]);
      }
      r[i1] += r[i0];
    } else {
      i = i0+1;      
      if (x[i0] > 0 && x[i] > 0) {
	yp = exp(0.5*(log(x[i0])+log(x[i])));
      } else if (x[i0] < 0 && x[i] < 0) {
	yp = -exp(0.5*(log(-x[i0])+log(-x[i])));
      } else {
	yp = 0.5*(x[i0]+x[i]);
      }
      r[i0+1] = r[i0] + 0.5*(_CNC[2][0]*(x[i0]+x[i0+1]) + _CNC[2][1]*yp);
      for (i = i0+2; i <= i1; i++) {
	r[i] = r[i-2] + _CNC[2][0]*(x[i-2]+x[i]) + _CNC[2][1]*x[i-1];
      }
    }
  } else {
    if (m >= 0) {
      r[i1] = x[i1];
      a = 0.0;
      for (i = i1-1; i > i0; i -= 2) {
	a += x[i];
      }
      r[i1] += 4.0*a;
      a = 0.0;
      k = i0+1;
      for (i = i1-2; i > k; i -= 2) {
	a += x[i];
      }
      r[i1] += 2.0*a;
      if (i == i0) {
	r[i1] += x[i0];
	r[i1] /= 3.0;
      } else {
	r[i1] += x[k];
	r[i1] /= 3.0;
	r[i1] += 0.5*(x[k] + x[i0]);
      }
      r[i1] += r[i0];
    } else {
      i = i1-1;
      if (x[i1] > 0 && x[i] > 0) {
	yp = exp(0.5*(log(x[i1])+log(x[i])));
      } else if (x[i1] < 0 && x[i] < 0) {
	yp = -exp(0.5*(log(-x[i1])+log(-x[i])));
      } else {
	yp = 0.5*(x[i1]+x[i]);
      }
      r[i1-1] = r[i1] + 0.5*(_CNC[2][0]*(x[i1]+x[i1-1]) + _CNC[2][1]*yp);
      for (i = i1-2; i >= i0; i--) {
	r[i] = r[i+2] + _CNC[2][0]*(x[i+2]+x[i]) + _CNC[2][1]*x[i+1];
      }
    }
  }

  return 0;
}

void PrepCEFCrossHeader(CEF_HEADER *h, double *data) {
  double *eusr, *x, bte, bms;
  int m, m1, j;

  eusr = h->egrid;
  m = h->n_egrid;
  m1 = m + 1;
  x = data+2+m1;
  BornFormFactorTE(&bte);
  bms = BornMass();
  data[0] = (h->te0*HARTREE_EV+bte)/bms;
  for (j = 0; j < m; j++) {
    x[j] = log((data[0] + eusr[j]*HARTREE_EV)/data[0]);
  }
  x[m] = eusr[m-1]/(data[0]/HARTREE_EV+eusr[m-1]);
}

void PrepCECrossHeader(CE_HEADER *h, double *data) {
  double *eusr, *x, bte, bms;
  int m, m1, j;
  
  eusr = h->usr_egrid;
  m = h->n_usr;
  m1 = m + 1;
  x = data+2+m1;
  BornFormFactorTE(&bte);
  bms = BornMass();
  data[0] = (h->te0*HARTREE_EV+bte)/bms;
  for (j = 0; j < m; j++) {
    x[j] = log((data[0] + eusr[j]*HARTREE_EV)/data[0]);
  }
  x[m] = eusr[m-1]/(data[0]/HARTREE_EV+eusr[m-1]);
}

void PrepCEFCrossRecord(CEF_RECORD *r, CEF_HEADER *h, double *data) {
  double *eusr, *x, *y, e;
  float *cs;
  int m, m1, j;
  EN_SRECORD *mem_en_table;
  int mem_en_table_size;

  mem_en_table = GetMemENFTable(&mem_en_table_size);

  eusr = h->egrid;
  m = h->n_egrid;
  m1 = m + 1;
  y = data + 2;
  x = y + m1;
  e = mem_en_table[r->upper].energy - mem_en_table[r->lower].energy;
  data[1] = r->bethe;
  
  cs = r->strength;
  for (j = 0; j < m; j++) {	
    y[j] = cs[j];
  }
  y[m] = r->born[0];
}

void PrepCECrossRecord(int k, CE_RECORD *r, CE_HEADER *h, 
		       double *data) {
  double *eusr, *x, *y, *w, e;
  float *cs;
  int m, m1, j, t;
  int j1, j2, t1, t2;
  EN_SRECORD *mem_en_table;
  int mem_en_table_size;

  mem_en_table = GetMemENTable(&mem_en_table_size);

  eusr = h->usr_egrid;
  m = h->n_usr;
  m1 = m + 1;
  y = data + 2;
  x = y + m1;
  w = x + m1;
  e = mem_en_table[r->upper].energy - mem_en_table[r->lower].energy;
  data[1] = r->bethe;
  
  cs = r->strength;
  if (k == 0) {
    if (h->msub) {
      j1 = mem_en_table[r->lower].j;
      j2 = mem_en_table[r->upper].j;
      for (j = 0; j < m; j++) {
	y[j] = 0.0;
      }
      t = 0;
      for (t1 = -j1; t1 <= 0; t1 += 2) {
	for (t2 = -j2; t2 <= j2; t2 += 2) {
	  for (j = 0; j < m; j++) {
	    y[j] += cs[t];
	    if (t1 != 0) y[j] += cs[t];
	    t++;
	  }
	}
      }
    } else {
      for (j = 0; j < m; j++) {	
	y[j] = cs[j];
      }
    }
    y[m] = r->born[0];
  }

  if (h->msub) {
    for (j = 0; j < m; j++) {
      if (y[j]) {
	w[j] = cs[k*m+j]/y[j];
      } else {
	w[j] = 1.0;
      }
    }
    if (r->bethe < 0) {
      w[j] = w[j-1];
    } else {
      w[j] = r->params[k];
    }
  }
}

double InterpolateCEFCross(double e, CEF_RECORD *r, CEF_HEADER *h, 
			   double *data) {
  double *x, *y, *w;
  int m, m1, n, one;
  double a, b, x0, y0, eth, e0, c, d, b0, b1;
  EN_SRECORD *mem_en_table;
  int mem_en_table_size;
  double bte, bms, eth1;

  mem_en_table = GetMemENFTable(&mem_en_table_size);

  eth = mem_en_table[r->upper].energy - mem_en_table[r->lower].energy;
  eth = eth * HARTREE_EV;

  a = 0.0;

  if (e < 0.0) return a;

  m = h->n_egrid;
  m1 = m + 1;
  x0 = log((data[0]+e)/data[0]);
  y = data + 2;
  x = y + m1;
  w = x + m1;
  
  BornFormFactorTE(&bte);
  bms = BornMass();
  if (x0 < x[m-1]) {
    n = 2;
    one = 1;
    if (fabs(bms-1.0) < EPS3 || x0 >= x[0]) {
      UVIP3P(n, m, x, y, one, &x0, &a);
      if (a < 0.0) a = 0.0;    
    } else {
      a = y[0] * pow(exp(x0-x[0]), 2.5);
    }
  } else {
    eth1 = (eth + bte*HARTREE_EV)/bms;
    x0 = e/(data[0] + e);
    y0 = y[m-1];
    if (data[1] > 0) {
      e0 = ((x[m]*data[0]/(1.0-x[m]))+eth1)/HARTREE_EV;
      b0 = 1.0 + FINE_STRUCTURE_CONST2*e0;
      b1 = 1.0 + FINE_STRUCTURE_CONST2*(e0-eth1/HARTREE_EV);
      y0 /= b0*b1;
      d = 2.0*e0*(1.0+0.5*FINE_STRUCTURE_CONST2*e0);
      c = FINE_STRUCTURE_CONST2*d;
      b = log(0.5*d*HARTREE_EV/eth) - c/(1.0+c);
      y0 -= data[1]*b;
      a = y[m] + (x0-1.0)*(y0-y[m])/(x[m]-1.0);
      e0 = (e + eth1)/HARTREE_EV;
      d = 2.0*e0*(1.0+0.5*FINE_STRUCTURE_CONST2*e0);
      c = FINE_STRUCTURE_CONST2*d;
      b = log(0.5*d*HARTREE_EV/eth) - c/(1.0+c);
      a += data[1]*b;
      b0 = 1.0 + FINE_STRUCTURE_CONST2*e0;
      b1 = 1.0 + FINE_STRUCTURE_CONST2*e/HARTREE_EV;
      a *= b0*b1;
    } else if (data[1]+1.0 == 1.0) {      
      e0 = ((x[m]*data[0]/(1.0-x[m]))+eth1)/HARTREE_EV;
      b0 = 1.0 + FINE_STRUCTURE_CONST2*e0;
      b1 = 1.0 + FINE_STRUCTURE_CONST2*(e0-eth1/HARTREE_EV);
      y0 /= b0*b1;
      a = y[m] + (x0-1.0)*(y0-y[m])/(x[m]-1.0);
      b0 = 1.0 + FINE_STRUCTURE_CONST2*(e+eth1)/HARTREE_EV;
      b1 = 1.0 + FINE_STRUCTURE_CONST2*e/HARTREE_EV;
      a *= b0*b1;
    } else {
      a = y[m] + (x0-1.0)*(y0-y[m])/(x[m]-1.0);
    }
  }

  return a;
}

double InterpolateCECross(double e, CE_RECORD *r, CE_HEADER *h, 
			  double *data, double *ratio) {
  double *x, *y, *w;
  int m, m1, n, one;
  double a, b, x0, y0, eth, e0, c, d, b0, b1;
  EN_SRECORD *mem_en_table;
  int mem_en_table_size;
  double bte, bms, eth1;

  mem_en_table = GetMemENTable(&mem_en_table_size);  

  eth = mem_en_table[r->upper].energy - mem_en_table[r->lower].energy;
  eth = eth * HARTREE_EV;

  a = 0.0;
  *ratio = 1.0;

  if (e < 0.0) return a;

  m = h->n_usr;
  m1 = m + 1;
  x0 = log((data[0]+e)/data[0]);
  y = data + 2;
  x = y + m1;
  w = x + m1;
  
  BornFormFactorTE(&bte);
  bms = BornMass();
  if (x0 < x[m-1]) {
    n = 2;
    one = 1;
    if (fabs(bms-1.0) < EPS3 || x0 >= x[0]) {
      UVIP3P(n, m, x, y, one, &x0, &a);
      if (a < 0.0) a = 0.0;
    } else {
      a = y[0] * pow(exp(x0-x[0]), 2.5);
    }
    if (h->msub) {
      UVIP3P(n, m, x, w, one, &x0, &b);
      if (b < 0.0) b = 0.0;
      a *= b;
      *ratio = b;
    }
  } else {
    x0 = e/(data[0] + e);
    y0 = y[m-1];
    eth1 = (eth + bte*HARTREE_EV)/bms;
    if (data[1] > 0) {
      e0 = ((x[m]*data[0]/(1.0-x[m]))+eth1)/HARTREE_EV;
      b0 = 1.0 + FINE_STRUCTURE_CONST2*e0;
      b1 = 1.0 + FINE_STRUCTURE_CONST2*(e0-eth1/HARTREE_EV);
      y0 /= b0*b1;
      d = 2.0*e0*(1.0+0.5*FINE_STRUCTURE_CONST2*e0);
      c = FINE_STRUCTURE_CONST2*d;
      b = log(0.5*d*HARTREE_EV/eth) - c/(1.0+c);
      y0 -= data[1]*b;
      a = y[m] + (x0-1.0)*(y0-y[m])/(x[m]-1.0);
      e0 = (e + eth1)/HARTREE_EV;
      d = 2.0*e0*(1.0+0.5*FINE_STRUCTURE_CONST2*e0);
      c = FINE_STRUCTURE_CONST2*d;
      b = log(0.5*d*HARTREE_EV/eth) - c/(1.0+c);
      a += data[1]*b;
      b0 = 1.0 + FINE_STRUCTURE_CONST2*e0;
      b1 = 1.0 + FINE_STRUCTURE_CONST2*e/HARTREE_EV;
      a *= b0*b1;
    } else if (data[1]+1.0 == 1.0) {
      e0 = ((x[m]*data[0]/(1.0-x[m]))+eth1)/HARTREE_EV;
      b0 = 1.0 + FINE_STRUCTURE_CONST2*e0;
      b1 = 1.0 + FINE_STRUCTURE_CONST2*(e0-eth1/HARTREE_EV);
      y0 /= b0*b1;
      a = y[m] + (x0-1.0)*(y0-y[m])/(x[m]-1.0);
      b0 = 1.0 + FINE_STRUCTURE_CONST2*(e+eth1)/HARTREE_EV;
      b1 = 1.0 + FINE_STRUCTURE_CONST2*e/HARTREE_EV;
      a *= b0*b1;
    } else {
      a = y[m] + (x0-1.0)*(y0-y[m])/(x[m]-1.0);
    }
    if (h->msub) {
      b = w[m] + (x0-1.0)*(w[m-1]-w[m])/(x[m]-1.0);
      a *= b;
      *ratio = b;
    }
  }

  return a;
}

int CEMFCross(char *ifn, char *ofn, int i0, int i1, 
	      int negy, double *egy, int mp) {
  F_HEADER fh;
  TFILE *f1;
  FILE *f2;
  int n, swp;
  CEF_HEADER h;
  CEF_RECORD r;
  CEMF_HEADER mh;
  CEMF_RECORD mr;
  int i, t, m, ith, iph;
  double data[2+(1+MAXNE)*3], e, cs, a;
  double eth, a1, cs1, k2, rp, e1, e0, b0, b1;
  double bte, bms, be;
  EN_SRECORD *mem_en_table;
  int mem_en_table_size;

  f1 = OpenFileRO(ifn, &fh, &swp); 
  if (f1 == NULL) {
    printf("File %s is not in FAC binary format\n", ifn);
    return -1;
  }  

  f2 = NULL;
  mem_en_table = GetMemENFTable(&mem_en_table_size);
  if (mem_en_table == NULL) {
    printf("Energy table has not been built in memory.\n");
    return -1;
  }
  if (strcmp(ofn, "-") == 0) {
    f2 = stdout;
  } else {
    f2 = fopen(ofn, "w");
  }
  if (f2 == NULL) {
    printf("cannot open file %s\n", ofn);
    return -1;
  }

  BornFormFactorTE(&bte);
  bms = BornMass();
  while (1) {
    n = ReadCEMFHeader(f1, &mh, swp);
    if (n == 0) break;
    CEMF2CEFHeader(&mh, &h);
    for (i = 0; i < mh.ntransitions; i++) {
      n = ReadCEMFRecord(f1, &mr, swp, &mh);
      if ((mr.lower == i0 || i0 < 0) && (mr.upper == i1 || i1 < 0)) {
	eth = mem_en_table[mr.upper].energy - mem_en_table[mr.lower].energy;
	e = eth*HARTREE_EV;
	fprintf(f2, "# %5d\t%5d\t%3d\t%5d\t%5d\t%3d\t%11.4E\t%5d\n",
		mr.lower, mem_en_table[mr.lower].p, mem_en_table[mr.lower].j,
		mr.upper, mem_en_table[mr.upper].p, mem_en_table[mr.upper].j,
		e, negy);
	be = (e + bte*HARTREE_EV)/bms;
	PrepCEFCrossHeader(&h, data);
	for (ith = 0; ith < mh.n_thetagrid; ith++) {
	  for (iph = 0; iph < mh.n_phigrid; iph++) {
	    fprintf(f2, "# %2d %2d %11.4E %11.4E\n",
		    ith, iph, mh.thetagrid[ith]*180/PI, mh.phigrid[iph]*180/PI);
	    CEMF2CEFRecord(&mr, &r, &mh, ith, iph);
	    PrepCEFCrossRecord(&r, &h, data);
	    for (t = 0; t < negy; t++) {
	      if (mp == 0) {
		e0 = egy[t];
		e1 = e0 - be;
	      } else {
		e1 = egy[t];
		e0 = e1 + be;
	      }
	      if (e1 > 0) {
		cs = InterpolateCEFCross(e1, &r, &h, data);
		a = e0/HARTREE_EV;
		b0 = 1.0 + FINE_STRUCTURE_CONST2*a;
		b1 = 1.0 + FINE_STRUCTURE_CONST2*(a-eth);
		a = a*(1.0+0.5*FINE_STRUCTURE_CONST2*a);
		a =  PI*AREA_AU20/(2.0*a);
		a *= cs;
		if (data[1] > 0.0) {	      
		  cs1 = data[1]*log(e0/e) + r.born[0];
		  k2 = e0/HARTREE_EV;
		  k2 = 2.0*k2*(1.0+0.5*FINE_STRUCTURE_CONST2*k2);
		  a1 = FINE_STRUCTURE_CONST2*k2;
		  a1 = a1/(1.0+a1);
		  a1 = data[1]*(log(0.5*k2/eth) - a1);
		  a1 += r.born[0];
		  a1 *= b0*b1;
		  k2 = cs1/a1;
		  rp = k2*(1.0+0.5*FINE_STRUCTURE_CONST2*e0/HARTREE_EV);
		  if (rp > 1.0) rp = 1.0;
		  cs1 = cs*k2;
		  a1 = a*rp;
		} else {
		  k2 = e0/HARTREE_EV;
		  rp = 1.0+0.5*FINE_STRUCTURE_CONST2*k2;
		  rp /= b0*b1;
		  cs1 = cs/(b0*b1);
		  if (rp > 1.0) rp = 1.0;
		  a1 = a*rp;
		}
	      } else {
		cs = 0.0;
		a = 0.0;
		cs1 = 0.0;
		a1 = 0.0;
	      }
	      fprintf(f2, "%11.4E %11.4E %11.4E %11.4E %11.4E\n",
		      e0, cs, a, cs1, a1);
	    }
	    fprintf(f2, "\n\n");   
	  }
	}
	     
	if (i0 >= 0 && i1 >= 0) {
	  free(mr.strength);
	  free(mr.bethe);
	  free(mr.born);
	  free(mh.tegrid);
	  free(mh.egrid);
	  free(mh.thetagrid);
	  free(mh.phigrid);
	  goto DONE;
	}
      }
      free(mr.strength);
      free(mr.bethe);
      free(mr.born);
    }
    free(mh.tegrid);
    free(mh.egrid);
    free(mh.thetagrid);
    free(mh.phigrid);
  }

 DONE:
  FCLOSE(f1);

  if (f2) {
    if (f2 != stdout) {
      fclose(f2);
    } else {
      fflush(f2);
    }
  }
  return 0;
}

int CEFCross(char *ifn, char *ofn, int i0, int i1, 
	     int negy, double *egy, int mp) {
  F_HEADER fh;
  TFILE *f1;
  FILE *f2;
  int n, swp;
  CEF_HEADER h;
  CEF_RECORD r;
  int i, t, m;
  double data[2+(1+MAXNE)*3], e, cs, a;
  double eth, a1, cs1, k2, rp, e1, e0, b0, b1;
  double bte, bms, be;
  EN_SRECORD *mem_en_table;
  int mem_en_table_size;
  
  f1 = OpenFileRO(ifn, &fh, &swp);
  if (f1 == NULL) {
    printf("cannot open file %s\n", ifn);
    return -1;
  }

  f2 = NULL;

  if (fh.type != DB_CEF || fh.nblocks == 0) {
    printf("File %s is not of DB_CE type\n", ifn);
    goto DONE;
  }

  mem_en_table = GetMemENFTable(&mem_en_table_size);
  if (mem_en_table == NULL) {
    printf("Energy table has not been built in memory.\n");
    return -1;
  }

  if (strcmp(ofn, "-") == 0) {
    f2 = stdout;
  } else {
    f2 = fopen(ofn, "w");
  }
  if (f2 == NULL) {
    printf("cannot open file %s\n", ofn);
    return -1;
  }

  BornFormFactorTE(&bte);
  bms = BornMass();
  while (1) {
    n = ReadCEFHeader(f1, &h, swp);
    if (n == 0) break;
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadCEFRecord(f1, &r, swp, &h);
      if ((r.lower == i0 || i0 < 0) && (r.upper == i1 || i1 < 0)) {
	PrepCEFCrossHeader(&h, data);
	eth = mem_en_table[r.upper].energy - mem_en_table[r.lower].energy;
	e = eth*HARTREE_EV;
	fprintf(f2, "# %5d\t%5d\t%3d\t%5d\t%5d\t%3d\t%11.4E\t%5d\n",
		r.lower, mem_en_table[r.lower].p, mem_en_table[r.lower].j,
		r.upper, mem_en_table[r.upper].p, mem_en_table[r.upper].j,
		e, negy);
	be = (e + bte*HARTREE_EV)/bms;
	PrepCEFCrossRecord(&r, &h, data);
	for (t = 0; t < negy; t++) {
	  if (mp == 0) {
	    e0 = egy[t];
	    e1 = e0 - be;
	  } else {
	    e1 = egy[t];
	    e0 = e1 + be;
	  }
	  if (e1 > 0) {	
	    cs = InterpolateCEFCross(e1, &r, &h, data);
	    a = e0/HARTREE_EV;
	    b0 = 1.0 + FINE_STRUCTURE_CONST2*a;
	    b1 = 1.0 + FINE_STRUCTURE_CONST2*(a-eth);
	    a = a*(1.0+0.5*FINE_STRUCTURE_CONST2*a);
	    a = PI*AREA_AU20/(2.0*a);
	    a *= cs;
	    if (data[1] > 0.0) {	      
	      cs1 = data[1]*log(e0/e) + r.born[0];
	      k2 = e0/HARTREE_EV;
	      k2 = 2.0*k2*(1.0+0.5*FINE_STRUCTURE_CONST2*k2);
	      a1 = FINE_STRUCTURE_CONST2*k2;
	      a1 = a1/(1.0+a1);
	      a1 = data[1]*(log(0.5*k2/eth) - a1);
	      a1 += r.born[0];
	      a1 *= b0*b1;
	      k2 = cs1/a1;
	      a1 = e0/HARTREE_EV;
	      rp = k2*(1.0+0.5*FINE_STRUCTURE_CONST2*a1);
	      if (rp > 1.0) rp = 1.0;
	      cs1 = cs*k2;
	      a1 = a*rp;
	    } else {
	      k2 = e0/HARTREE_EV;
	      rp = 1.0+0.5*FINE_STRUCTURE_CONST2*k2;
	      rp /= b0*b1;
	      cs1 = cs/(b0*b1);
	      if (rp > 1.0) rp = 1.0;
	      a1 = a*rp;
	    }
	  } else {
	    cs = 0.0;
	    a = 0.0;
	    cs1 = 0.0;
	    a1 = 0.0;
	  }
	  fprintf(f2, "%11.4E %11.4E %11.4E %11.4E %11.4E\n",
		  e0, cs, a, cs1, a1);
	}
	fprintf(f2, "\n\n");
      
	if (i0 >= 0 && i1 >= 0) {
	  free(r.strength);
	  free(h.tegrid);
	  free(h.egrid);
	  goto DONE;
	}
      }
      free(r.strength);
    }
    free(h.tegrid);
    free(h.egrid);
  }

 DONE:
  FCLOSE(f1);

  if (f2) {
    if (f2 != stdout) {
      fclose(f2);
    } else {
      fflush(f2);
    }
  }
  return 0;
}

int CECross(char *ifn, char *ofn, int i0, int i1, 
	    int negy, double *egy, int mp) {
  F_HEADER fh;
  TFILE *f1;
  FILE *f2;
  int n, swp;
  CE_HEADER h;
  CE_RECORD r;
  int i, t, m, k;
  double data[2+(1+MAXNUSR)*3], e, cs, a, ratio;
  double eth, a1, cs1, k2, rp, e1, e0, b0, b1;
  double bte, bms, be;
  EN_SRECORD *mem_en_table;
  int mem_en_table_size;
  
  f1 = OpenFileRO(ifn, &fh, &swp);
  if (f1 == NULL) {
    printf("cannot open file %s\n", ifn);
    return -1;
  }
  
  f2 = NULL;

  if (fh.type != DB_CE || fh.nblocks == 0) {
    printf("File %s is not of DB_CE type\n", ifn);
    goto DONE;
  }

  mem_en_table = GetMemENTable(&mem_en_table_size);
  if (mem_en_table == NULL) {
    printf("Energy table has not been built in memory.\n");
    return -1;
  }

  if (strcmp(ofn, "-") == 0) {
    f2 = stdout;
  } else {
    f2 = fopen(ofn, "w");
  }
  if (f2 == NULL) {
    printf("cannot open file %s\n", ofn);
    return -1;
  }
    
  BornFormFactorTE(&bte);
  bms = BornMass();
  while (1) {
    n = ReadCEHeader(f1, &h, swp);
    if (n == 0) break;
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadCERecord(f1, &r, swp, &h);
      if ((r.lower == i0 || i0 < 0) && (r.upper == i1 || i1 < 0)) {
	eth = mem_en_table[r.upper].energy - mem_en_table[r.lower].energy;
	e = eth*HARTREE_EV;
	fprintf(f2, "# %5d\t%2d\t%5d\t%2d\t%11.4E\t%5d\t%d\n",
		r.lower, mem_en_table[r.lower].j,
		r.upper, mem_en_table[r.upper].j,
		e, negy, r.nsub);
	be = (e + bte*HARTREE_EV)/bms;	
	PrepCECrossHeader(&h, data);
	for (k = 0; k < r.nsub; k++) {
	  PrepCECrossRecord(k, &r, &h, data);
	  for (t = 0; t < negy; t++) {
	    if (mp == 0) {
	      e0 = egy[t];
	      e1 = e0 - be;
	    } else {
	      e1 = egy[t];
	      e0 = e1 + be;
	    }
	    if (e1 > 0) {
	      cs = InterpolateCECross(e1, &r, &h, data, &ratio);
	      a = e0/HARTREE_EV;
	      b0 = 1.0 + FINE_STRUCTURE_CONST2*a;
	      b1 = 1.0 + FINE_STRUCTURE_CONST2*(a-eth);
	      a = a*(1.0+0.5*FINE_STRUCTURE_CONST2*a);
	      a = PI*AREA_AU20/(2.0*a);
	      if (!h.msub) a /= (mem_en_table[r.lower].j+1.0);
	      a *= cs;
	      if (data[1] > 0.0) {	      
		cs1 = data[1]*log(e0/e) + r.born[0];
		k2 = e0/HARTREE_EV;
		k2 = 2.0*k2*(1.0+0.5*FINE_STRUCTURE_CONST2*k2);
		a1 = FINE_STRUCTURE_CONST2*k2;
		a1 = a1/(1.0+a1);
		a1 = data[1]*(log(0.5*k2/eth) - a1);
		a1 += r.born[0];
		a1 *= b0*b1;
		k2 = cs1/a1;
		a1 = e0/HARTREE_EV;
		rp = k2*(1.0+0.5*FINE_STRUCTURE_CONST2*a1);
		if (rp > 1.0) rp = 1.0;
		cs1 = cs*k2;
		a1 = a*rp;
	      } else {
		k2 = e0/HARTREE_EV;
		rp = 1.0+0.5*FINE_STRUCTURE_CONST2*k2;
		rp /= (b0*b1);
		cs1 = cs/(b0*b1);
		if (rp > 1.0) rp = 1.0;
		a1 = a*rp;
	      }
	    } else {
	      cs = 0.0;
	      a = 0.0;
	      cs1 = 0.0;
	      a1 = 0.0;
	      ratio = 0.0;
	    }
	    fprintf(f2, "%11.4E %11.4E %11.4E %11.4E %11.4E %11.4E\n",
		    e0, cs, a, cs1, a1, ratio);
	  }
	  fprintf(f2, "\n\n");
	}
	if (i0 >= 0 && i1 >= 0) {
	  if (h.msub || h.qk_mode == QK_FIT) free(r.params);
	  free(r.strength);
	  free(h.tegrid);
	  free(h.egrid);
	  free(h.usr_egrid);
	  goto DONE;
	}
      }
      if (h.msub || h.qk_mode == QK_FIT) free(r.params);
      free(r.strength);
    }
    free(h.tegrid);
    free(h.egrid);
    free(h.usr_egrid);
  }

 DONE:
  FCLOSE(f1);

  if (f2) {
    if (f2 != stdout) {
      fclose(f2);
    } else {
      fflush(f2);
    }
  }
  return 0;
}

int CEMFMaxwell(char *ifn, char *ofn, int i0, int i1, 
		int nt, double *temp) {
  F_HEADER fh;
  TFILE *f1;
  FILE *f2;
  int n, swp;
  CEF_HEADER h;
  CEF_RECORD r;
  CEMF_HEADER mh;
  CEMF_RECORD mr;
  int i, t, m, p, ith, iph;
  double data[2+(1+MAXNE)*4], e, cs, a, c;
  double *xg = gauss_xw[0];  
  double *wg = gauss_xw[1];
  EN_SRECORD *mem_en_table;
  int mem_en_table_size;
  
  f1 = OpenFileRO(ifn, &fh, &swp);
  if (f1 == NULL) {
    printf("cannot open file %s\n", ifn);
    return -1;
  }
  f2 = NULL;
  mem_en_table = GetMemENFTable(&mem_en_table_size);
  if (mem_en_table == NULL) {
    printf("Energy table has not been built in memory.\n");
    return -1;
  }

  if (strcmp(ofn, "-") == 0) {
    f2 = stdout;
  } else {
    f2 = fopen(ofn, "w");
  }
  if (f2 == NULL) {
    printf("cannot open file %s\n", ofn);
    return -1;
  }

  while (1) {
    n = ReadCEMFHeader(f1, &mh, swp);
    if (n == 0) break;
    CEMF2CEFHeader(&mh, &h);
    for (i = 0; i < mh.ntransitions; i++) {
      n = ReadCEMFRecord(f1, &mr, swp, &mh);
      if ((mr.lower == i0 || i0 < 0) && (mr.upper == i1 || i1 < 0)) {	
	PrepCEFCrossHeader(&h, data);
	e = mem_en_table[mr.upper].energy - mem_en_table[mr.lower].energy;
	e *= HARTREE_EV;
	fprintf(f2, "# %5d\t%5d\t%3d\t%5d\t%5d\t%3d\t%11.4E\t%5d\n",
		mr.lower, mem_en_table[mr.lower].p, mem_en_table[mr.lower].j,
		mr.upper, mem_en_table[mr.upper].p, mem_en_table[mr.upper].j,
		e, nt);
	for (ith = 0; ith < mh.n_thetagrid; ith++) {
	  for (iph = 0; iph < mh.n_phigrid; iph++) {
	    fprintf(f2, "# %2d %2d %11.4E %11.4E\n",
		    ith, iph, mh.thetagrid[ith]*180.0/PI, mh.phigrid[iph]*180.0/PI);
	    CEMF2CEFRecord(&mr, &r, &mh, ith, iph);
	    PrepCEFCrossRecord(&r, &h, data);
	    for (t = 0; t < nt; t++) {
	      cs = 0.0;
	      for (p = 0; p < 15; p++) {
		a = temp[t]*xg[p];
		c = InterpolateCEFCross(a, &r, &h, data);
		cs += wg[p]*c;
	      }
	      a = 217.16*sqrt(HARTREE_EV/(2.0*temp[t]));
	      a *= cs*exp(-e/temp[t]);
	      fprintf(f2, "%11.4E\t%11.4E\t%11.4E\n", 
		      temp[t], cs, a);
	    }
	    fprintf(f2, "\n\n");
	  }
	}
	if (i0 >= 0 && i1 >= 0) {
	  free(mr.strength);
	  free(mr.bethe);
	  free(mr.born);
	  free(mh.tegrid);
	  free(mh.egrid);
	  free(mh.thetagrid);
	  free(mh.phigrid);
	  goto DONE;
	}    
      }
      free(mr.strength);
      free(mr.bethe);
      free(mr.born);
    }
    free(mh.tegrid);
    free(mh.egrid);
    free(mh.thetagrid);
    free(mh.phigrid);
  }

 DONE:
  FCLOSE(f1);

  if (f2) {
    if (f2 != stdout) {
      fclose(f2);
    } else {
      fflush(f2);
    }
  }

  return 0;
}

int CEFMaxwell(char *ifn, char *ofn, int i0, int i1, 
	       int nt, double *temp) {
  F_HEADER fh;
  TFILE *f1;
  FILE *f2;
  int n, swp;
  CEF_HEADER h;
  CEF_RECORD r;
  int i, t, m, p;
  double data[2+(1+MAXNE)*4], e, cs, a, c;
  double *xg = gauss_xw[0];  
  double *wg = gauss_xw[1];
  EN_SRECORD *mem_en_table;
  int mem_en_table_size;
  
  f1 = OpenFileRO(ifn, &fh, &swp);
  if (f1 == NULL) {
    printf("cannot open file %s\n", ifn);
    return -1;
  }
   
  f2 = NULL;

  if (fh.type != DB_CEF || fh.nblocks == 0) {
    printf("File %s is not of DB_CE type\n", ifn);
    goto DONE;
  }
   
  mem_en_table = GetMemENFTable(&mem_en_table_size);
  
  if (mem_en_table == NULL) {
    printf("Energy table has not been built in memory.\n");
    return -1;
  }

  if (strcmp(ofn, "-") == 0) {
    f2 = stdout;
  } else {
    f2 = fopen(ofn, "w");
  }
  if (f2 == NULL) {
    printf("cannot open file %s\n", ofn);
    return -1;
  }

  while (1) {
    n = ReadCEFHeader(f1, &h, swp);
    if (n == 0) break;
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadCEFRecord(f1, &r, swp, &h);
      if ((r.lower == i0 || i0 < 0) && (r.upper == i1 || i1 < 0)) {	
	PrepCEFCrossHeader(&h, data);
	e = mem_en_table[r.upper].energy - mem_en_table[r.lower].energy;
	e *= HARTREE_EV;
	fprintf(f2, "# %5d\t%5d\t%3d\t%5d\t%5d\t%3d\t%11.4E\t%5d\n",
		r.lower, mem_en_table[r.lower].p, mem_en_table[r.lower].j,
		r.upper, mem_en_table[r.upper].p, mem_en_table[r.upper].j,
		e, nt);
	PrepCEFCrossRecord(&r, &h, data);
	for (t = 0; t < nt; t++) {
	  cs = 0.0;
	  for (p = 0; p < 15; p++) {
	    a = temp[t]*xg[p];
	    c = InterpolateCEFCross(a, &r, &h, data);
	    cs += wg[p]*c;
	  }
	  a = 217.16*sqrt(HARTREE_EV/(2.0*temp[t]));
	  a *= cs*exp(-e/temp[t]);
	  fprintf(f2, "%11.4E\t%11.4E\t%11.4E\n", 
		  temp[t], cs, a);
	}
	fprintf(f2, "\n\n");
	if (i0 >= 0 && i1 >= 0) {
	  free(r.strength);
	  free(h.tegrid);
	  free(h.egrid);
	  goto DONE;
	}    
      }
      free(r.strength);
    }
    free(h.tegrid);
    free(h.egrid);
  }

 DONE:
  FCLOSE(f1);

  if (f2) {
    if (f2 != stdout) {
      fclose(f2);
    } else {
      fflush(f2);
    }
  }

  return 0;
}

int CEMaxwell(char *ifn, char *ofn, int i0, int i1, 
	      int nt, double *temp) {
  F_HEADER fh;
  TFILE *f1;
  FILE *f2;
  int n, swp;
  CE_HEADER h;
  CE_RECORD r;
  int i, t, m, k, p;
  double data[2+(1+MAXNUSR)*4], e, cs, a, c, ratio;
  double *xg = gauss_xw[0];  
  double *wg = gauss_xw[1]; 
  EN_SRECORD *mem_en_table;
  int mem_en_table_size;
  
  f1 = OpenFileRO(ifn, &fh, &swp);
  if (f1 == NULL) {
    printf("cannot open file %s\n", ifn);
    return -1;
  }
   
  f2 = NULL;
  if (n == 0) {
    printf("File %s is not in FAC binary format\n", ifn);
    goto DONE;
  }
  
  if (fh.type != DB_CE || fh.nblocks == 0) {
    printf("File %s is not of DB_CE type\n", ifn);
    goto DONE;
  }
   
  mem_en_table = GetMemENTable(&mem_en_table_size);
  
  if (mem_en_table == NULL) {
    printf("Energy table has not been built in memory.\n");
    return -1;
  }

  if (strcmp(ofn, "-") == 0) {
    f2 = stdout;
  } else {
    f2 = fopen(ofn, "w");
  }
  if (f2 == NULL) {
    printf("cannot open file %s\n", ofn);
    return -1;
  }

  while (1) {
    n = ReadCEHeader(f1, &h, swp);
    if (n == 0) break;
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadCERecord(f1, &r, swp, &h);
      if ((r.lower == i0 || i0 < 0) && (r.upper == i1 || i1 < 0)) {
	e = mem_en_table[r.upper].energy - mem_en_table[r.lower].energy;
	e *= HARTREE_EV;
	fprintf(f2, "# %5d\t%2d\t%5d\t%2d\t%11.4E\t%5d\t%d\n",
		r.lower, mem_en_table[r.lower].j,
		r.upper, mem_en_table[r.upper].j,
		e, nt, r.nsub);
	PrepCECrossHeader(&h, data);
	for (k = 0; k < r.nsub; k++) {
	  PrepCECrossRecord(k, &r, &h, data);
	  for (t = 0; t < nt; t++) {
	    cs = 0.0;
	    for (p = 0; p < 15; p++) {
	      a = temp[t]*xg[p];
	      c = InterpolateCECross(a, &r, &h, data, &ratio);
	      cs += wg[p]*c;
	    }
	    a = 217.16*sqrt(HARTREE_EV/(2.0*temp[t]));
	    a *= cs*exp(-e/temp[t]);
	    if (!h.msub) a /= (mem_en_table[r.lower].j+1.0);
	    fprintf(f2, "%11.4E\t%11.4E\t%11.4E\n", 
		    temp[t], cs, a);
	  }
	  fprintf(f2, "\n\n");
	}
	if (i0 >= 0 && i1 >= 0) {
	  if (h.msub || h.qk_mode == QK_FIT) free(r.params);
	  free(r.strength);
	  free(h.tegrid);
	  free(h.egrid);
	  free(h.usr_egrid);
	  goto DONE;
	}
      }
      if (h.msub || h.qk_mode == QK_FIT) free(r.params);
      free(r.strength);
    }
    free(h.tegrid);
    free(h.egrid);
    free(h.usr_egrid);
  }

 DONE:
  FCLOSE(f1);

  if (f2) {
    if (f2 != stdout) {
      fclose(f2);
    } else {
      fflush(f2);
    }
  }
  return 0;
}

double InterpolateCICross(double e1, double eth, CI_RECORD *r, CI_HEADER *h) {
  double y[MAXNE], x, a, b, tc;
  int i;

  for (i = 0; i < h->n_usr; i++) {
    y[i] = r->strength[i];
  }
  if (e1 < 0) return 0.0; 
  if (e1 < h->usr_egrid[0] || e1 > h->usr_egrid[h->n_usr-1]) {
    x = 1.0 + e1/eth;
    a = 1.0/x;
    b = 1.0 - a;
    tc = r->params[0]*log(x) + r->params[1]*b*b;
    tc += r->params[2]*a*b + r->params[3]*a*a*b;
    return tc;
  } else {
    UVIP3P(2, h->n_usr, h->usr_egrid, y, 1, &e1, &tc);
    return tc;
  }
}
  
double InterpolateCIMCross(double e1, double eth, CIM_RECORD *r, CIM_HEADER *h,
			   int q) {
  double y[MAXNE], z[MAXNE], x, tc;
  int i, j;

  for (i = 0; i < h->n_usr; i++) {
    j = q*h->n_usr + i;
    y[i] = r->strength[j];
    z[i] = log(1.0 + h->usr_egrid[i]/eth);
  }
  if (e1 < 0) return 0.0; 
  x = log(1.0 + e1/eth);
  if (e1 < h->usr_egrid[0] || e1 > h->usr_egrid[h->n_usr-1]) {
    UVIP3P(1, h->n_usr, z, y, 1, &x, &tc);
  } else {
    UVIP3P(2, h->n_usr, z, y, 1, &x, &tc);
  }
  return tc;
}

int TotalCICross(char *ifn, char *ofn, int ilev, 
		 int negy, double *egy, int imin, int imax) {
  F_HEADER fh;
  TFILE *f1;
  FILE *f2;
  int n, swp;
  CI_HEADER h;
  CI_RECORD r;
  int i, t, nb, m;
  double *c, tc, a, e, bte, bms, be;
  EN_SRECORD *mem_en_table;
  int mem_en_table_size;

  mem_en_table = GetMemENTable(&mem_en_table_size);
  
  if (mem_en_table == NULL) {
    printf("Energy table has not been built in memory.\n");
    return -1;
  }

  nb = 0;
  
  f1 = OpenFileRO(ifn, &fh, &swp);
  if (f1 == NULL) {
    printf("cannot open file %s\n", ifn);
    return -1;
  }

  if (strcmp(ofn, "-") == 0) {
    f2 = stdout;
  } else {
    f2 = fopen(ofn, "w");
  }
  if (f2 == NULL) {
    printf("cannot open file %s\n", ofn);
    return -1;
  }
  
  if (fh.type != DB_CI || fh.nblocks == 0) {
    printf("File %s is not of DB_CI type\n", ifn);
    goto DONE;
  }
  
  c = (double *) malloc(sizeof(double)*negy);
  for (i = 0; i < negy; i++) {
    c[i] = 0.0;
    egy[i] /= HARTREE_EV;
  }

  if (imin < 0) imin = 0;
  if (imax < 0) imax = mem_en_table_size - 1;

  BornFormFactorTE(&bte);
  bms = BornMass();
  while (1) {
    n = ReadCIHeader(f1, &h, swp);
    if (n == 0) break;
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadCIRecord(f1, &r, swp, &h);
      if (n == 0) break;
      if (r.b != ilev) continue;
      if (r.f < imin || r.f > imax) continue;
      e = mem_en_table[r.f].energy - mem_en_table[r.b].energy; 
      be = (e + bte)/bms;
      for (t = 0; t < negy; t++) {
	if (egy[t] < be) continue;
	tc = InterpolateCICross(egy[t]-be, e, &r, &h);
	a = egy[t]*(1.0 + FINE_STRUCTURE_CONST2*egy[t]);
	tc *= AREA_AU20/(2.0*a*(mem_en_table[r.b].j + 1.0));
	c[t] += tc;
      }
      free(r.params);
      free(r.strength);
    }

    free(h.tegrid);
    free(h.egrid);
    free(h.usr_egrid);
    
    nb++;
  }

  fprintf(f2, "#Energy (eV)   CI Cross (10^-20 cm2)\n");
  for (t = 0; t < negy; t++) {
    fprintf(f2, " %11.4E    %15.8E\n", egy[t]*HARTREE_EV, c[t]);
  }

  free(c);  

 DONE:
  FCLOSE(f1);

  if (f2 != stdout) {
    fclose(f2);
  } else {
    fflush(f2);
  }

  return nb;
} 

int CICross(char *ifn, char *ofn, int i0, int i1,
	    int negy, double *egy, int mp) {
  F_HEADER fh;
  TFILE *f1;
  FILE *f2;
  int n, swp;
  CI_HEADER h;
  CI_RECORD r;
  int i, t, nb, m;
  double tc, a, b, e, e0, bte, bms, be;
  EN_SRECORD *mem_en_table;
  int mem_en_table_size;

  mem_en_table = GetMemENTable(&mem_en_table_size);
  
  if (mem_en_table == NULL) {
    printf("Energy table has not been built in memory.\n");
    return -1;
  }

  nb = 0;
  
  f1 = OpenFileRO(ifn, &fh, &swp);
  if (f1 == NULL) {
    printf("cannot open file %s\n", ifn);
    return -1;
  }

  if (strcmp(ofn, "-") == 0) {
    f2 = stdout;
  } else {
    f2 = fopen(ofn, "w");
  }
  if (f2 == NULL) {
    printf("cannot open file %s\n", ofn);
    return -1;
  }

  if (fh.type != DB_CI || fh.nblocks == 0) {
    printf("File %s is not of DB_CI type\n", ifn);
    goto DONE;
  }
  
  for (i = 0; i < negy; i++) {
    egy[i] /= HARTREE_EV;
  }
  BornFormFactorTE(&bte);
  bms = BornMass();
  while (1) {
    n = ReadCIHeader(f1, &h, swp);
    if (n == 0) break;
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadCIRecord(f1, &r, swp, &h);
      if (n == 0) break;      
      if ((r.b == i0 || i0 < 0) && (r.f == i1 || i1 < 0)) {
	e = mem_en_table[r.f].energy - mem_en_table[r.b].energy;
	fprintf(f2, "# %5d\t%2d\t%5d\t%2d\t%11.4E\t%5d\n",
		r.b, mem_en_table[r.b].j,
		r.f, mem_en_table[r.f].j,
		e*HARTREE_EV, negy);
	be = (e + bte)/bms;
	for (t = 0; t < negy; t++) {
	  if (mp == 0) {
	    e0 = egy[t];
	  } else {
	    e0 = egy[t] + be;
	  }
	  if (e0 < be) {
	    tc = 0.0; 
	    b = 0.0;
	  } else {
	    tc = InterpolateCICross(e0-be, e, &r, &h);
	    b = tc;
	    a = e0*(1.0 + FINE_STRUCTURE_CONST2*e0);
	    tc *= AREA_AU20/(2.0*a*(mem_en_table[r.b].j + 1.0));
	  }
	  fprintf(f2, "%11.4E %11.4E %11.4E\n",
		  e0*HARTREE_EV, b, tc);
	}
	fprintf(f2, "\n\n");
      
	if (i0 >= 0 && i1 >= 0) {
	  free(r.params);
	  free(r.strength);
	  free(h.tegrid);
	  free(h.egrid);
	  free(h.usr_egrid);
	  goto DONE;
	}
      }
      free(r.params);
      free(r.strength);
    }

    free(h.tegrid);
    free(h.egrid);
    free(h.usr_egrid);
    
    nb++;
  }

 DONE:
  FCLOSE(f1);

  if (f2 != stdout) {
    fclose(f2);
  } else {
    fflush(f2);
  }

  return nb;
} 

int CIMaxwell(char *ifn, char *ofn, int i0, int i1,
	      int negy, double *egy) {
  F_HEADER fh;
  TFILE *f1;
  FILE *f2;
  int n, swp;
  CI_HEADER h;
  CI_RECORD r;
  int i, t, nb, m, p;
  double tc, e, e0, cs;
  double *xg = gauss_xw[0];  
  double *wg = gauss_xw[1]; 
  EN_SRECORD *mem_en_table;
  int mem_en_table_size;

  mem_en_table = GetMemENTable(&mem_en_table_size);
  
  if (mem_en_table == NULL) {
    printf("Energy table has not been built in memory.\n");
    return -1;
  }

  nb = 0;
  
  f1 = OpenFileRO(ifn, &fh, &swp);
  if (f1 == NULL) {
    printf("cannot open file %s\n", ifn);
    return -1;
  }

  if (strcmp(ofn, "-") == 0) {
    f2 = stdout;
  } else {
    f2 = fopen(ofn, "w");
  }
  if (f2 == NULL) {
    printf("cannot open file %s\n", ofn);
    return -1;
  }

  if (fh.type != DB_CI || fh.nblocks == 0) {
    printf("File %s is not of DB_CI type\n", ifn);
    goto DONE;
  }
  
  for (i = 0; i < negy; i++) {
    egy[i] /= HARTREE_EV;
  }
  while (1) {
    n = ReadCIHeader(f1, &h, swp);
    if (n == 0) break;
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadCIRecord(f1, &r, swp, &h);
      if (n == 0) break;      
      if ((r.b == i0 || i0 < 0) && (r.f == i1 || i1 < 0)) {
	e = mem_en_table[r.f].energy - mem_en_table[r.b].energy;
	fprintf(f2, "# %5d\t%2d\t%5d\t%2d\t%11.4E\t%5d\n",
		r.b, mem_en_table[r.b].j,
		r.f, mem_en_table[r.f].j,
		e*HARTREE_EV, negy);
	for (t = 0; t < negy; t++) {
	  cs = 0.0;
	  for (p = 0; p < 15; p++) {
	    e0 =  egy[t]*xg[p];
	    tc = InterpolateCICross(e0, e, &r, &h);
	    cs += wg[p]*tc;
	  }
	  tc = (217.16/PI)*sqrt(1.0/(2.0*egy[t]));
	  tc *= cs*exp(-e/egy[t]);
	  tc /= (mem_en_table[r.b].j + 1.0);
	  fprintf(f2, "%11.4E %11.4E %11.4E\n",
		  egy[t]*HARTREE_EV, cs, tc);
	}
	fprintf(f2, "\n\n");
      
	if (i0 >= 0 && i1 >= 0) {
	  free(r.params);
	  free(r.strength);
	  free(h.tegrid);
	  free(h.egrid);
	  free(h.usr_egrid);
	  goto DONE;
	}
      }
      free(r.params);
      free(r.strength);
    }

    free(h.tegrid);
    free(h.egrid);
    free(h.usr_egrid);
    
    nb++;
  }

 DONE:
  FCLOSE(f1);

  if (f2 != stdout) {
    fclose(f2);
  } else {
    fflush(f2);
  }

  return nb;
} 

int CIMCross(char *ifn, char *ofn, int i0, int i1,
	     int negy, double *egy, int mp) {
  F_HEADER fh;
  TFILE *f1;
  FILE *f2;
  int n, swp;
  CIM_HEADER h;
  CIM_RECORD r;
  int i, t, nb, m, k;
  double tc, a, b, e, e0, bte, bms, be;
  EN_SRECORD *mem_en_table;
  int mem_en_table_size;

  mem_en_table = GetMemENTable(&mem_en_table_size);
  
  if (mem_en_table == NULL) {
    printf("Energy table has not been built in memory.\n");
    return -1;
  }

  nb = 0;
  
  f1 = OpenFileRO(ifn, &fh, &swp);
  if (f1 == NULL) {
    printf("cannot open file %s\n", ifn);
    return -1;
  }

  if (strcmp(ofn, "-") == 0) {
    f2 = stdout;
  } else {
    f2 = fopen(ofn, "w");
  }
  if (f2 == NULL) {
    printf("cannot open file %s\n", ofn);
    return -1;
  }
  
  if (fh.type != DB_CIM || fh.nblocks == 0) {
    printf("File %s is not of DB_CIM type\n", ifn);
    goto DONE;
  }
  
  for (i = 0; i < negy; i++) {
    egy[i] /= HARTREE_EV;
  }
  BornFormFactorTE(&bte);
  bms = BornMass();
  while (1) {
    n = ReadCIMHeader(f1, &h, swp);
    if (n == 0) break;
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadCIMRecord(f1, &r, swp, &h);
      if (n == 0) break;      
      if ((r.b == i0 || i0 < 0) && (r.f == i1 || i1 < 0)) {
	e = mem_en_table[r.f].energy - mem_en_table[r.b].energy;
	fprintf(f2, "# %5d\t%2d\t%5d\t%2d\t%11.4E\t%5d\t%2d\n",
		r.b, mem_en_table[r.b].j,
		r.f, mem_en_table[r.f].j,
		e*HARTREE_EV, negy, r.nsub);
	be = (e + bte)/bms;
	for (k = 0; k < r.nsub; k++) {
	  for (t = 0; t < negy; t++) {
	    if (mp == 0) {
	      e0 = egy[t];
	    } else {
	      e0 = egy[t] + be;
	    }
	    if (e0 < be) {
	      tc = 0.0; 
	      b = 0.0;
	    } else {
	      tc = InterpolateCIMCross(e0-be, e, &r, &h, k);
	      b = tc;
	      a = e0*(1.0 + FINE_STRUCTURE_CONST2*e0);
	      tc *= AREA_AU20/(2.0*a);
	    }
	    fprintf(f2, "%11.4E %11.4E %11.4E\n",
		    e0*HARTREE_EV, b, tc);
	  }
	  fprintf(f2, "\n\n");
	}
	if (i0 >= 0 && i1 >= 0) {
	  free(r.strength);
	  free(h.egrid);
	  free(h.usr_egrid);
	  goto DONE;
	}
      }
      free(r.strength);      
    }
    free(h.egrid);
    free(h.usr_egrid);
    
    nb++;
  }

 DONE:
  FCLOSE(f1);

  if (f2 != stdout) {
    fclose(f2);
  } else {
    fflush(f2);
  }

  return nb;
} 

int TotalPICross(char *ifn, char *ofn, int ilev, 
		 int negy, double *egy, int imin, int imax) {
  F_HEADER fh;
  TFILE *f1;
  FILE *f2;
  int n, swp;
  RR_HEADER h;
  RR_RECORD r;
  int i, t, nb, m;
  float e, eph, ee, phi;
  double *xusr, *dstrength, *c, tc, emax;
  double x, y;
  int np=3, one=1, nele;
  EN_SRECORD *mem_en_table;
  int mem_en_table_size;

  mem_en_table = GetMemENTable(&mem_en_table_size);
  
  if (mem_en_table == NULL) {
    printf("Energy table has not been built in memory.\n");
    return -1;
  }

  nb = 0;

  f1 = OpenFileRO(ifn, &fh, &swp);
  if (f1 == NULL) {
    printf("cannot open file %s\n", ifn);
    return -1;
  }

  if (strcmp(ofn, "-") == 0) {
    f2 = stdout;
  } else {
    f2 = fopen(ofn, "w");
  }
  if (f2 == NULL) {
    printf("cannot open file %s\n", ofn);
    return -1;
  }

  if (fh.type != DB_RR || fh.nblocks == 0) {
    printf("File %s is not of DB_RR type\n", ifn);
    goto DONE;
  }
  
  c = (double *) malloc(sizeof(double)*negy);
  for (i = 0; i < negy; i++) {
    c[i] = 0.0;
    egy[i] /= HARTREE_EV;
  }
  
  if (imin < 0) imin = 0;
  if (imax < 0) imax = mem_en_table_size - 1;

  while(1) {
    n = ReadRRHeader(f1, &h, swp);
    if (n == 0) break;
    nele = h.nele;
    xusr = (double *) malloc(sizeof(double)*h.n_usr); 
    dstrength = (double *) malloc(sizeof(double)*h.n_usr);
    emax = h.usr_egrid[h.n_usr-1];
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadRRRecord(f1, &r, swp, &h);
      if (n == 0) break;
      if (r.b != ilev) continue;
      if (r.f < imin || r.f > imax) continue;
      e = mem_en_table[r.f].energy - mem_en_table[r.b].energy;

      for (t = 0; t < h.n_usr; t++) {
	dstrength[t] = log(r.strength[t]);
	xusr[t] = log(1.0 + h.usr_egrid[t]/e);
      }
      
      for (t = 0; t < negy; t++) {
	eph = egy[t];
	ee = eph - e;
	if (ee <= 0.0) continue;
	if (h.qk_mode != QK_FIT || ee <= emax) {
	  x = log(eph/e);
	  UVIP3P(np, h.n_usr, xusr, dstrength, one, &x, &tc);
	  tc = exp(tc);
	} else {
	  x = (ee + r.params[3])/r.params[3];
	  y = (1 + r.params[2])/(sqrt(x) + r.params[2]);
	  tc = (-3.5 - r.kl + 0.5*r.params[1])*log(x) + r.params[1]*log(y);
	  if (r.params[0] > 0.0) {
	    tc = tc + log(r.params[0]*(eph/(ee+r.params[3])));
	    tc = exp(tc);
	  } else {
	    tc = 0.0;
	  }
	}
	phi = 2.0*PI*FINE_STRUCTURE_CONST*tc*AREA_AU20;	
	c[t] += phi/(mem_en_table[r.b].j + 1.0);
      }
      if (h.qk_mode == QK_FIT) free(r.params);
      free(r.strength);
    }

    free(dstrength);
    free(h.tegrid);
    free(h.egrid);
    free(h.usr_egrid);      
    free(xusr);
    
    nb++;
  }

  fprintf(f2, "# Energy (eV)   PI Cross (10^-20 cm2)\n");
  for (t = 0; t < negy; t++) {
    fprintf(f2, " %12.5E    %15.8E\n", egy[t]*HARTREE_EV, c[t]);
  }

  free(c);  

 DONE:
  FCLOSE(f1);

  if (f2 != stdout) {
    fclose(f2);
  } else {
    fflush(f2);
  }

  return nb;
}
  
int RRCross(char *ifn, char *ofn, int i0, int i1,
	    int negy, double *egy, int mp) {
  F_HEADER fh;
  TFILE *f1;
  FILE *f2;
  int n, swp;
  RR_HEADER h;
  RR_RECORD r;
  int i, t, nb, m;
  float e, eph, ee, phi, rr;
  double *xusr, *dstrength, tc, emax;
  double x, y;
  int np=3, one=1, nele;
  EN_SRECORD *mem_en_table;
  int mem_en_table_size;

  mem_en_table = GetMemENTable(&mem_en_table_size);
  
  if (mem_en_table == NULL) {
    printf("Energy table has not been built in memory.\n");
    return -1;
  }

  nb = 0;

  f1 = OpenFileRO(ifn, &fh, &swp);
  if (f1 == NULL) {
    printf("cannot open file %s\n", ifn);
    return -1;
  }

  if (strcmp(ofn, "-") == 0) {
    f2 = stdout;
  } else {
    f2 = fopen(ofn, "w");
  }
  if (f2 == NULL) {
    printf("cannot open file %s\n", ofn);
    return -1;
  }

  if (fh.type != DB_RR || fh.nblocks == 0) {
    printf("File %s is not of DB_RR type\n", ifn);
    goto DONE;
  }
  
  for (i = 0; i < negy; i++) {
    egy[i] /= HARTREE_EV;
  }
  
  while(1) {
    n = ReadRRHeader(f1, &h, swp);
    if (n == 0) break;
    nele = h.nele;
    xusr = (double *) malloc(sizeof(double)*h.n_usr); 
    dstrength = (double *) malloc(sizeof(double)*h.n_usr);
    emax = h.usr_egrid[h.n_usr-1];
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadRRRecord(f1, &r, swp, &h);
      if (n == 0) break;
      if ((r.b == i0 || i0 < 0) && (r.f == i1 || i1 < 0)) {
	e = mem_en_table[r.f].energy - mem_en_table[r.b].energy;
	fprintf(f2, "# %5d\t%2d\t%5d\t%2d\t%11.4E\t%5d\n",
		r.b, mem_en_table[r.b].j,
		r.f, mem_en_table[r.f].j,
		e*HARTREE_EV, negy);	
	for (t = 0; t < h.n_usr; t++) {
	  dstrength[t] = log(r.strength[t]);
	  xusr[t] = log(1.0 + h.usr_egrid[t]/e);
	}	
	for (t = 0; t < negy; t++) {
	  if (mp == 0) {
	    eph = egy[t];
	    ee = eph - e;
	  } else {
	    eph = egy[t] + e;
	    ee = egy[t];
	  }
	  if (ee <= 0.0) {
	    tc = 0.0;
	    phi = 0.0;
	    rr = 0.0;
	  } else {
	    if (h.qk_mode != QK_FIT || ee <= emax) {
	      x = log(eph/e);
	      UVIP3P(np, h.n_usr, xusr, dstrength, one, &x, &tc);
	      tc = exp(tc);	    
	    } else {
	      x = (ee + r.params[3])/r.params[3];
	      y = (1 + r.params[2])/(sqrt(x) + r.params[2]);
	      tc = (-3.5 - r.kl + 0.5*r.params[1])*log(x) + r.params[1]*log(y);
	      if (r.params[0] > 0.0) {
		tc = tc + log(r.params[0]*(eph/(ee+r.params[3])));
		tc = exp(tc);
	      } else {
		tc = 0.0;
	      }	  
	    }
	    phi = 2.0*PI*FINE_STRUCTURE_CONST*tc*AREA_AU20;
	    rr = phi * pow(FINE_STRUCTURE_CONST*eph, 2) / (2.0*ee);
	    rr /= (mem_en_table[r.f].j + 1.0);	
	    phi /= (mem_en_table[r.b].j + 1.0);
	  }
	  fprintf(f2, "%11.4E %11.4E %11.4E %11.4E %11.4E\n",
		  ee*HARTREE_EV, eph*HARTREE_EV, tc, rr, phi);
	}
	fprintf(f2, "\n\n");

	if (i0 >= 0 && i1 >= 0) {
	  if (h.qk_mode == QK_FIT) free(r.params);
	  free(r.strength);
	  free(dstrength);
	  free(h.tegrid);
	  free(h.egrid);
	  free(h.usr_egrid);      
	  free(xusr);
	  goto DONE;
	}
      }
      if (h.qk_mode == QK_FIT) free(r.params);
      free(r.strength);
    }

    free(dstrength);
    free(h.tegrid);
    free(h.egrid);
    free(h.usr_egrid);      
    free(xusr);
    
    nb++;
  }

 DONE:
  FCLOSE(f1);

  if (f2 != stdout) {
    fclose(f2);
  } else {
    fflush(f2);
  }

  return nb;
}
  
int RRMaxwell(char *ifn, char *ofn, int i0, int i1,
	      int negy, double *egy) {
  F_HEADER fh;
  TFILE *f1;
  FILE *f2;
  int n, swp;
  RR_HEADER h;
  RR_RECORD r;
  int i, t, nb, m, p;
  float e, eph, ee;
  double *xusr, *dstrength, tc, cs, rr, emax;
  double x, y;
  int np=3, one=1, nele;
  double *xg = gauss_xw[0];  
  double *wg = gauss_xw[1]; 
  EN_SRECORD *mem_en_table;
  int mem_en_table_size;

  mem_en_table = GetMemENTable(&mem_en_table_size);
  
  if (mem_en_table == NULL) {
    printf("Energy table has not been built in memory.\n");
    return -1;
  }

  nb = 0;

  f1 = OpenFileRO(ifn, &fh, &swp);
  if (f1 == NULL) {
    printf("cannot open file %s\n", ifn);
    return -1;
  }

  if (strcmp(ofn, "-") == 0) {
    f2 = stdout;
  } else {
    f2 = fopen(ofn, "w");
  }
  if (f2 == NULL) {
    printf("cannot open file %s\n", ofn);
    return -1;
  }

  if (fh.type != DB_RR || fh.nblocks == 0) {
    printf("File %s is not of DB_RR type\n", ifn);
    goto DONE;
  }
  
  for (i = 0; i < negy; i++) {
    egy[i] /= HARTREE_EV;
  }
  while(1) {
    n = ReadRRHeader(f1, &h, swp);
    if (n == 0) break;
    nele = h.nele;
    xusr = (double *) malloc(sizeof(double)*h.n_usr); 
    dstrength = (double *) malloc(sizeof(double)*h.n_usr);
    emax = h.usr_egrid[h.n_usr-1];
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadRRRecord(f1, &r, swp, &h);
      if (n == 0) break;
      if ((r.b == i0 || i0 < 0) && (r.f == i1 || i1 < 0)) {
	e = mem_en_table[r.f].energy - mem_en_table[r.b].energy;
	fprintf(f2, "# %5d\t%2d\t%5d\t%2d\t%11.4E\t%5d\n",
		r.b, mem_en_table[r.b].j,
		r.f, mem_en_table[r.f].j,
		e*HARTREE_EV, negy);	
	for (t = 0; t < h.n_usr; t++) {
	  dstrength[t] = log(r.strength[t]);
	  xusr[t] = log(1.0 + h.usr_egrid[t]/e);
	}	
	for (t = 0; t < negy; t++) {
	  cs = 0.0;
	  for (p = 0; p < 15; p++) {
	    ee = egy[t]*xg[p];
	    eph = ee + e;	    
	    if (h.qk_mode != QK_FIT || ee <= emax) {
	      x = log(eph/e);
	      UVIP3P(np, h.n_usr, xusr, dstrength, one, &x, &tc);
	      tc = exp(tc);	    
	    } else {
	      x = (ee + r.params[3])/r.params[3];
	      y = (1 + r.params[2])/(sqrt(x) + r.params[2]);
	      tc = (-3.5 - r.kl + 0.5*r.params[1])*log(x) + r.params[1]*log(y);
	      if (r.params[0] > 0.0) {
		tc = tc + log(r.params[0]*(eph/(ee+r.params[3])));
		tc = exp(tc);
	      } else {
		tc = 0.0;
	      }	  
	    }	    
	    tc *= 2.0*FINE_STRUCTURE_CONST;
	    tc *= pow(FINE_STRUCTURE_CONST*eph, 2);
	    cs += tc * wg[p];
	  }
	  rr = 217.16 * sqrt(1.0/(2.0*egy[t]));
	  rr *= cs;
	  rr /= (mem_en_table[r.f].j + 1.0);	
	  fprintf(f2, "%11.4E %11.4E %11.4E\n",
		  egy[t]*HARTREE_EV, tc, rr);
	}
	fprintf(f2, "\n\n");
      
	if (i0 >= 0 && i1 >= 0) {
	  if (h.qk_mode == QK_FIT) free(r.params);
	  free(r.strength);
	  free(dstrength);
	  free(h.tegrid);
	  free(h.egrid);
	  free(h.usr_egrid);      
	  free(xusr);
	  goto DONE;
	}
      }
      if (h.qk_mode == QK_FIT) free(r.params);
      free(r.strength);
    }

    free(dstrength);
    free(h.tegrid);
    free(h.egrid);
    free(h.usr_egrid);      
    free(xusr);
    
    nb++;
  }

 DONE:
  FCLOSE(f1);

  if (f2 != stdout) {
    fclose(f2);
  } else {
    fflush(f2);
  }

  return nb;
}

double RRCrossHn(double z, double e, int n) {
  double x, z2;
  double y;
  double f = 7.499E-6;
  
  z2 = z*z;
  x = n;
  e = 2.0*e;
  y = f*z2*z2/(x*e*(z2 + e*x*x));

  return y;
}

int TotalRRCross(char *ifn, char *ofn, int ilev, 
		 int negy, double *egy, int n0, int n1, int nmax,
		 int imin, int imax) {
  F_HEADER fh;
  TFILE *f1;
  FILE *f2;
  int n, swp;
  RR_HEADER h;
  RR_RECORD r;
  int i, t, nb, m;
  float e, eph, ee, phi, rr;
  double *xusr, *dstrength, *c, tc, emax;
  double x, y;
  int np=3, one=1, nele;
  EN_SRECORD *mem_en_table;
  int mem_en_table_size;

  mem_en_table = GetMemENTable(&mem_en_table_size);

  if (mem_en_table == NULL) {
    printf("Energy table has not been built in memory.\n");
    return -1;
  }

  nb = 0;

  f1 = OpenFileRO(ifn, &fh, &swp);
  if (f1 == NULL) {
    printf("cannot open file %s\n", ifn);
    return -1;
  }

  if (strcmp(ofn, "-") == 0) {
    f2 = stdout;
  } else {
    f2 = fopen(ofn, "w");
  }
  if (f2 == NULL) {
    printf("cannot open file %s\n", ofn);
    return -1;
  }

  if (fh.type != DB_RR || fh.nblocks == 0) {
    printf("File %s is not of DB_RR type\n", ifn);
    goto DONE;
  }
  
  c = (double *) malloc(sizeof(double)*negy);
  for (i = 0; i < negy; i++) {
    c[i] = 0.0;
    egy[i] /= HARTREE_EV;
  }
  
  if (imin < 0) imin = 0;
  if (imax < 0) imax = mem_en_table_size - 1;

  while(1) {
    n = ReadRRHeader(f1, &h, swp);
    if (n == 0) break;
    nele = h.nele;
    xusr = (double *) malloc(sizeof(double)*h.n_usr); 
    dstrength = (double *) malloc(sizeof(double)*h.n_usr);
    emax = h.usr_egrid[h.n_usr-1];
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadRRRecord(f1, &r, swp, &h);
      if (n == 0) break;
      if (r.f != ilev) continue;
      if (r.b < imin || r.b > imax) continue;
      e = mem_en_table[r.f].energy - mem_en_table[r.b].energy;
      
      for (t = 0; t < h.n_usr; t++) {
	dstrength[t] = log(r.strength[t]);
	xusr[t] = log(1.0 + h.usr_egrid[t]/e);
      }
      
      for (t = 0; t < negy; t++) {
	ee = egy[t];
	eph = ee + e;
	if (h.qk_mode != QK_FIT || ee <= emax) {
	  x = log(eph/e);
	  UVIP3P(np, h.n_usr, xusr, dstrength, one, &x, &tc);
	  tc = exp(tc);
	} else {
	  x = (ee + r.params[3])/r.params[3];
	  y = (1 + r.params[2])/(sqrt(x) + r.params[2]);
	  tc = (-3.5 - r.kl + 0.5*r.params[1])*log(x) + r.params[1]*log(y);
	  if (r.params[0] > 0.0) {
	    tc = tc + log(r.params[0]*(eph/(ee+r.params[3])));
	    tc = exp(tc);
	  } else {
	    tc = 0.0;
	  }
	}
	phi = 2.0*PI*FINE_STRUCTURE_CONST*tc*AREA_AU20;
	ee *= 1.0 + 0.5*FINE_STRUCTURE_CONST2*ee;
	rr = phi * pow(FINE_STRUCTURE_CONST*eph, 2) / (2.0*ee);
	rr /= (mem_en_table[r.f].j + 1.0);
	c[t] += rr;
      }
      if (h.qk_mode == QK_FIT) free(r.params);
      free(r.strength);
    }  

    free(dstrength);
    free(h.tegrid);
    free(h.egrid);
    free(h.usr_egrid);      
    free(xusr);
    
    nb++;
  }
  
  x = fh.atom - nele + 1.0;
  for (n = n0+1; n < n1; n++) {
    for (t = 0; t < negy; t++) {
      c[t] += AREA_AU20*RRCrossHn(x, egy[t], n);
    }
  }
  for (n = n1+1; n <= nmax; n++) {
    for (t = 0; t < negy; t++) {
      c[t] += AREA_AU20*RRCrossHn(x, egy[t], n);
    }
  }
  
  fprintf(f2, "#Energy (eV)   RR Cross (10^-20 cm2)\n");
  for (t = 0; t < negy; t++) {
    fprintf(f2, " %11.4E    %15.8E\n", egy[t]*HARTREE_EV, c[t]);
  }

  free(c);

 DONE:
  FCLOSE(f1);
  if (f2 != stdout) {
    fclose(f2);
  } else {
    fflush(f2);
  }

  return nb;
}

double voigt(double alpha, double v) {
  double a[8], b[8], c[8];
  double v2, v3, fac1, fac2;
  double p1,p2,p3,p4,p5,p6,p7;
  double o1,o2,o3,o4,o5,o6,o7;
  double q1,q2;
  double r1,r2;
  double H, n, w, vw, vb;
  int i;

  a[1]=122.607931777104326;
  a[2]=214.382388694706425;
  a[3]=181.928533092181549;
  a[4]=93.155580458134410;
  a[5]=30.180142196210589;
  a[6]=5.912626209773153;
  a[7]=0.564189583562615;

  b[1]=122.607931773875350;
  b[2]=352.730625110963558;
  b[3]=457.334478783897737;
  b[4]=348.703917719495792;
  b[5]=170.354001821091472;
  b[6]=53.992906912940207;
  b[7]=10.479857114260399;

  c[1]=0.5641641;
  c[2]=0.8718681;
  c[3]=1.474395;
  c[4]=-19.57862;
  c[5]=802.4513;
  c[6]=-4850.316;
  c[7]=8031.468;
 
  if (alpha <= 1e-3 && fabs(v) > 2.5) {
    v2 = v*v;
    v3 = 1.0;
    fac1 = c[1];
    fac2 = c[1] * (v2 - 1.0);     
    for (i=1; i<=7; i++) {
      v3     = v3 * v2;
      fac1 = fac1 + c[i] / v3;
      fac2 = fac2 + c[i] / v3 * (v2 - i);
    }
    H = exp(-v2) * (1. + fac2*alpha*alpha * (1. - 2.*v2)) + fac1 * (alpha/v2);
  } else {
    p1 = alpha;
    o1 = -v;  
    p2 = (p1 * alpha + o1 * v);
    o2 = (o1 * alpha - p1 * v);
    p3 = (p2 * alpha + o2 * v);
    o3 = (o2 * alpha - p2 * v);
    p4 = (p3 * alpha + o3 * v);
    o4 = (o3 * alpha - p3 * v);
    p5 = (p4 * alpha + o4 * v);
    o5 = (o4 * alpha - p4 * v);
    p6 = (p5 * alpha + o5 * v);
    o6 = (o5 * alpha - p5 * v);
    p7 = (p6 * alpha + o6 * v);
    o7 = (o6 * alpha - p6 * v);
    
    q1 = a[1] + p1 * a[2] + p2 * a[3] + p3 * a[4] +
      p4 * a[5] + p5 * a[6] + p6 * a[7];
    r1 =        o1 * a[2] + o2 * a[3] + o3 * a[4] +
      o4 * a[5] + o5 * a[6] + o6 * a[7];
    q2 = b[1] + p1 * b[2] + p2 * b[3] + p3 * b[4] +
      p4 * b[5] + p5 * b[6] + p6 * b[7] + p7;
    r2 =        o1 * b[2] + o2 * b[3] + o3 * b[4] +
      o4 * b[5] + o5 * b[6] + o6 * b[7] + o7;

    H = (q1 * q2 + r1 * r2) / (q2 * q2 + r2 * r2);
  }

  return H/sqrt(PI);
}

int InterpCross(char *ifn, char *ofn, int i0, int i1, 
		int negy, double *egy, int mp) {  
  F_HEADER fh;
  TFILE *f1;
  int n, swp;
  
  f1 = OpenFileRO(ifn, &fh, &swp);
  if (f1 == NULL) {
    printf("cannot open file %s\n", ifn);
    return -1;
  }
  FCLOSE(f1);

  switch (fh.type) {
  case DB_CE:
    CECross(ifn, ofn, i0, i1, negy, egy, mp);
    break;
  case DB_CEF:
    CEFCross(ifn, ofn, i0, i1, negy, egy, mp);
    break;
  case DB_CEMF:
    CEMFCross(ifn, ofn, i0, i1, negy, egy, mp);
    break;
  case DB_CI:
    CICross(ifn, ofn, i0, i1, negy, egy, mp);
    break;
  case DB_CIM:
    CIMCross(ifn, ofn, i0, i1, negy, egy, mp);
    break;
  case DB_RR:
    RRCross(ifn, ofn, i0, i1, negy, egy, mp);
    break;
  default:
    printf("Improper FAC binary file for cross sections\n");
  }

  return 0;
}

int MaxwellRate(char *ifn, char *ofn, int i0, int i1, 
		int negy, double *egy) {  
  F_HEADER fh;
  TFILE *f1;
  int n, swp;
  
  f1 = OpenFileRO(ifn, &fh, &swp);
  if (f1 == NULL) {
    printf("cannot open file %s\n", ifn);
    return -1;
  }
  FCLOSE(f1);
  if (n == 0) {
    printf("File %s is not in FAC binary format\n", ifn);
    return -1;
  }  

  switch (fh.type) {
  case DB_CE:
    CEMaxwell(ifn, ofn, i0, i1, negy, egy);
    break;
  case DB_CEF:
    CEFMaxwell(ifn, ofn, i0, i1, negy, egy);
    break;
  case DB_CEMF:
    CEMFMaxwell(ifn, ofn, i0, i1, negy, egy);
    break;
  case DB_CI:
    CIMaxwell(ifn, ofn, i0, i1, negy, egy);
    break;
  case DB_RR:
    RRMaxwell(ifn, ofn, i0, i1, negy, egy);
    break;
  default:
    printf("Improper FAC binary file for rate coefficients\n");
  }

  return 0;
}

double ModifyEntry(double r, MOD_RECORD *m) {
  switch (m->op) {
  case 0:
    return m->c;
  case 1:
    return r + m->c;
  case 2:
    return r*m->c;
  default:
    return r;
  }
}

int CompareModRecord(const void *p1, const void *p2) {
  MOD_RECORD *m1, *m2;
  
  m1 = (MOD_RECORD *) p1;
  m2 = (MOD_RECORD *) p2;
  
  if (m1->m < m2->m) return -1;
  else if (m1->m > m2->m) return 1;
  else {
    if (m1->i0 < m2->i0) return -1;
    else if (m1->i0 > m2->i0) return 1;
    else {
      if (m1->i1 < m2->i1) return -1;
      else if (m1->i1 > m2->i1) return 1;
      else return 0;
    }
  }
}

int CompareModRecordI0(const void *p1, const void *p2) {
  MOD_RECORD *m1, *m2;
  
  m1 = (MOD_RECORD *) p1;
  m2 = (MOD_RECORD *) p2;
  
  if (m1->i0 < m2->i0) return -1;
  else if (m1->i0 > m2->i0) return 1;
  else return 0;
}

int CompareModRecordI1(const void *p1, const void *p2) {
  MOD_RECORD *m1, *m2;
  
  m1 = (MOD_RECORD *) p1;
  m2 = (MOD_RECORD *) p2;
  
  if (m1->i1 < m2->i1) return -1;
  else if (m1->i1 > m2->i1) return 1;
  else return 0;
}
	
void ModifyEN(int nc, MOD_RECORD *mr, int nr, MOD_RECORD *pr, 
	      TFILE *f0, TFILE *f1, F_HEADER *fh, int swp) {
  int n, i, k;
  double e0, e1;
  EN_HEADER h;
  EN_RECORD r;
  MOD_RECORD mr0;
  
  if (nr > 0) {
    e1 = 1e30;
    for (i = 0; i < nr; i++) {
      if (e1 > pr->c) e1 = pr->c;
    }
  } else {
    e1 = 0.0;
  }
  e0 = 1E30;
  while (1) {
    n = ReadENHeader(f0, &h, swp);
    if (n == 0) break;
    InitFile(f1, fh, &h);
    for (i = 0; i < h.nlevels; i++) {
      n = ReadENRecord(f0, &r, swp);
      if (e0 > 0) e0 = r.energy;
      if (n == 0) break;
      mr0.m = 0;
      mr0.i0 = 0;
      mr0.i1 = r.ilev;
      if (nr == 0) {
	k = Bisect(&mr0, nc, sizeof(MOD_RECORD), mr, CompareModRecordI1);
	if (k >= 0) {
	  r.energy = ModifyEntry(r.energy, mr+k);
	}
      } else {
	k = Bisect(&mr0, nc, sizeof(MOD_RECORD), mr, CompareModRecordI1);
	if (k >= 0) {
	  mr0.i0 = 0;
	  mr0.i1 = mr[k].i0;
	  k = Bisect(&mr0, nr, sizeof(MOD_RECORD), pr, CompareModRecordI1);
	  if (k >= 0) {
	    r.energy = e0 + pr[k].c-e1;
	  }
	}
      }
      WriteENRecord(f1, &r);
    }
    DeinitFile(f1, fh);
  }
}

void ModifyTR(int nc, MOD_RECORD *mr, int nr, MOD_RECORD *pr, 
	      TFILE *f0, TFILE *f1, F_HEADER *fh, int swp) {
  int n, i, k;
  TR_HEADER h;
  TR_RECORD r, *r0;
  TR_EXTRA rx;
  MOD_RECORD mr0;
  double a;
  
  while (1) {
    n = ReadTRHeader(f0, &h, swp);
    if (n == 0) break;
    InitFile(f1, fh, &h);
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadTRRecord(f0, &r, &rx, swp);
      if (n == 0) break;
      if (nr == 0) {
	mr0.m = h.multipole;
	mr0.i0 = r.lower;
	mr0.i1 = r.upper;
	k = Bisect(&mr0, nc, sizeof(MOD_RECORD), mr, CompareModRecord);
	if (k >= 0) {
	  a = r.strength;
	  r.strength = ModifyEntry(a, mr+k);
	}
      } else {
	mr0.i1 = r.lower;
	k = Bisect(&mr0, nc, sizeof(MOD_RECORD), mr, CompareModRecordI1);
	if (k >= 0) {
	  mr0.i0 = mr[k].i0;
	  mr0.i1 = r.upper;
	  k = Bisect(&mr0, nc, sizeof(MOD_RECORD), mr, CompareModRecordI1);
	  if (k >= 0) {
	    mr0.i1 = mr[k].i0;
	    mr0.m = h.multipole;
	    k = Bisect(&mr0, nr, sizeof(MOD_RECORD), pr, CompareModRecord);
	    if (k >= 0) {
	      r0 = (TR_RECORD *) pr[k].r;
	      r.strength = r0->strength;
	    }
	  }
	}
      }
      WriteTRRecord(f1, &r, &rx);
    }
    DeinitFile(f1, fh);
  }
}

void ModifyCE(int nc, MOD_RECORD *mr, int nr, MOD_RECORD *pr, 
	      TFILE *f0, TFILE *f1, F_HEADER *fh, int swp) {
  int n, i, k, j, t, p;
  CE_HEADER h, *h0;
  CE_RECORD r, *r0;  
  MOD_RECORD mr0;
  double a, data[2+(1+MAXNUSR)*3];

  if (GetMemENTable(&n) == NULL) {
    printf("must build the in-memory enegy table for CE modification\n");
    return;
  }

  while (1) {
    n = ReadCEHeader(f0, &h, swp);
    if (n == 0) break;
    InitFile(f1, fh, &h);
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadCERecord(f0, &r, swp, &h);
      if (n == 0) break;
      if (nr == 0) {
	if (r.bethe > 0) {
	  mr0.m = -1;
	} else if (r.bethe + 1.0 == 1.0) {
	  mr0.m = -2;
	}
	mr0.i0 = r.lower;
	mr0.i1 = r.upper;
	k = Bisect(&mr0, nc, sizeof(MOD_RECORD), mr, CompareModRecord);
	if (k >= 0) {
	  if (r.bethe > 0) {
	    a = r.bethe;
	    r.bethe = ModifyEntry(a, mr+k);
	  }
	  a = r.born[0];
	  r.born[0] = ModifyEntry(a, mr+k);
	  p = 0;
	  for (j = 0; j < r.nsub; j++) {
	    if (h.msub) {
	      a = r.params[j];
	      r.params[j] = ModifyEntry(a, mr+k);
	    }
	    for (t = 0; t < h.n_usr; t++) {
	      a = r.strength[p];
	      r.strength[p] = ModifyEntry(a, mr+k);
	      p++;
	    }
	  }
	}
      } else {
	mr0.i1 = r.lower;
	k = Bisect(&mr0, nc, sizeof(MOD_RECORD), mr, CompareModRecordI1);
	if (k >= 0) {
	  mr0.i0 = mr[k].i0;
	  mr0.i1 = r.upper;
	  k = Bisect(&mr0, nc, sizeof(MOD_RECORD), mr, CompareModRecordI1);
	  if (k >= 0) {
	    mr0.i1 = mr[k].i0;
	    mr0.m = 0;
	    k = Bisect(&mr0, nr, sizeof(MOD_RECORD), pr, CompareModRecord);
	    if (k >= 0) {
	      r0 = (CE_RECORD *) pr[k].r;
	      h0 = (CE_HEADER *) pr[k].h;
	      r.bethe = r0->bethe;
	      r.born[0] = r0->born[0];
	      r.born[1] = r0->born[1];
	      p = 0;
	      PrepCECrossHeader(h0, data);	      
	      for (j = 0; j < r.nsub; j++) {
		if (h.msub) {
		  r.params[j] = r0->params[j];
		}
		PrepCECrossRecord(j, r0, h0, data);
		for (t = 0; t < h.n_usr; t++) {
		  r.strength[p] = InterpolateCECross(h.usr_egrid[t]*HARTREE_EV,
						     r0, h0, data, &a);
		  p++;
		}
	      }
	    }
	  }
	}
      }
      WriteCERecord(f1, &r);
      if (h.msub || h.qk_mode == QK_FIT) free(r.params);
      free(r.strength);
    }
    DeinitFile(f1, fh);
    free(h.tegrid);
    free(h.egrid);
    free(h.usr_egrid);
  }
}

int ReadModRecord(char *fn, MOD_RECORD **mr, int nr) {
  FILE *f;
#define NBUF 512
  char buf[NBUF];
  int nc, n, m, i0, i1, op;
  double c;

  f = fopen(fn, "r");
  if (f == NULL) {
    printf ("cannot open file %s", fn);
    return 0;
  }
  
  nc = 0;
  while (1) {
    if (NULL == fgets(buf, NBUF, f)) break;
    if (nr > 0) {
      n = sscanf(buf, "%d %d", &i0, &i1);
      if (n != 2) continue;
    } else {
      n = sscanf(buf, "%d %d %d %lf %d", &m, &i0, &i1, &c, &op);
      if (n != 5) continue;
    }
    nc++;
  }

  if (nc == 0) {
    fclose(f);
    return 0;
  }

  (*mr) = malloc(sizeof(MOD_RECORD)*nc);
  fseek(f, 0, SEEK_SET);
  nc = 0;
  while (1) {
    if (NULL == fgets(buf, NBUF, f)) break;
    if (nr > 0) {
      n = sscanf(buf, "%d %d", &i0, &i1);
      if (n != 2) continue;
      (*mr)[nc].m = 0;
      (*mr)[nc].i0 = i0;
      (*mr)[nc].i1 = i1;
      (*mr)[nc].c = 0.0;
      (*mr)[nc].op = 0;
    } else {
      n = sscanf(buf, "%d %d %d %lf %d", &m, &i0, &i1, &c, &op);
      if (n != 5) continue;
      (*mr)[nc].m = m;
      (*mr)[nc].i0 = i0;
      (*mr)[nc].i1 = i1;
      (*mr)[nc].c = c;
      (*mr)[nc].op = op;
    }
    nc++;
  }
  
  fclose(f);

  if (nr == 0) {
    qsort(*mr, nc, sizeof(MOD_RECORD), CompareModRecord);
  } else {
    qsort(*mr, nc, sizeof(MOD_RECORD), CompareModRecordI1);
  }

  return nc;
#undef NBUF
}
	
int ReadENTable(char *fn, int *nh, EN_HEADER **h, 
		int *nr, EN_RECORD **r, MOD_RECORD **mr) {
  F_HEADER fh;
  EN_HEADER h0;
  int n, swp, i, j, m;
  TFILE *f;

  f = OpenFileRO(fn, &fh, &swp);
  if (f == NULL) return -1;

  if (fh.type != DB_EN) {
    FCLOSE(f);
    return -1;
  }
  
  *nh = fh.nblocks;
  *nr = 0;
  while (1) {
    n = ReadENHeader(f, &h0, swp);
    if (n == 0) break;
    *nr += h0.nlevels;
    FSEEK(f, h0.length, SEEK_CUR);
  }
  if (*nh == 0 || *nr == 0) {
    FCLOSE(f);
    return -1;
  }
  *h = malloc(sizeof(EN_HEADER)*(*nh));
  *r = malloc(sizeof(EN_RECORD)*(*nr));
  if (mr) {
    *mr = malloc(sizeof(MOD_RECORD)*(*nr));
  }
  FSEEK(f, 0, SEEK_SET);
  n = ReadFHeader(f, &fh, &swp);
  i = 0;
  j = 0;
  while (1) {
    n = ReadENHeader(f, (*h)+i, swp);
    if (n == 0) break;
    for (m = 0; m < (*h)[i].nlevels; m++) {
      n = ReadENRecord(f, (*r)+j, swp);      
      if (n == 0) break;
      if (mr) (*mr)[j].h = (*h) + i;
      j++;
    }
    i++;
  }
  
  if (mr) {
    for (i = 0; i < (*nr); i++) {
      (*mr)[i].m = 0;
      (*mr)[i].i0 = 0;
      (*mr)[i].i1 = (*r)[i].ilev;
      (*mr)[i].c = (*r)[i].energy;
      (*mr)[i].r = (*r) + i;
    }
  }

  FCLOSE(f);
  return 0;
}
	
int ReadTRTable(char *fn, int *nh, TR_HEADER **h, 
		int *nr, TR_RECORD **r, MOD_RECORD **mr) {
  F_HEADER fh;
  TR_HEADER h0;
  TR_EXTRA rx;
  int n, swp, i, j, m;
  TFILE *f;
  
  f = OpenFileRO(fn, &fh, &swp);
  if (f == NULL) return -1;

  if (fh.type != DB_TR) {
    FCLOSE(f);
    return -1;
  }
  
  *nh = fh.nblocks;
  *nr = 0;
  while (1) {
    n = ReadTRHeader(f, &h0, swp);
    if (n == 0) break;
    *nr += h0.ntransitions;
    FSEEK(f, h0.length, SEEK_CUR);
  }
  if (*nh == 0 || *nr == 0) {
    FCLOSE(f);
    return -1;
  }
  *h = malloc(sizeof(TR_HEADER)*(*nh));
  *r = malloc(sizeof(TR_RECORD)*(*nr));
  if (mr) {
    *mr = malloc(sizeof(MOD_RECORD)*(*nr));
  }
  FSEEK(f, 0, SEEK_SET);
  n = ReadFHeader(f, &fh, &swp);
  i = 0;
  j = 0;
  while (1) {
    n = ReadTRHeader(f, (*h)+i, swp);
    if (n == 0) break;
    for (m = 0; m < (*h)[i].ntransitions; m++) {
      n = ReadTRRecord(f, (*r)+j, &rx, swp);
      if (n == 0) break;
      if (mr) {
	(*mr)[j].m = (*h)[i].multipole;
	(*mr)[j].h = (*h) + i;
      }
      j++;
    }
    i++;
  }
  
  if (mr) {
    for (i = 0; i < (*nr); i++) {
      (*mr)[i].i0 = (*r)[i].lower;
      (*mr)[i].i1 = (*r)[i].upper;
      (*mr)[i].r = (*r) + i;
    }
  }

  FCLOSE(f);
  return 0;
}
	
int ReadCETable(char *fn, int *nh, CE_HEADER **h, 
		int *nr, CE_RECORD **r, MOD_RECORD **mr) {
  F_HEADER fh;
  CE_HEADER h0;
  int n, swp, i, j, m;
  TFILE *f;

  f = OpenFileRO(fn, &fh, &swp);
  if (f == NULL) return -1;

  if (fh.type != DB_CE) {
    FCLOSE(f);
    return -1;
  }
  
  *nh = fh.nblocks;
  *nr = 0;
  while (1) {
    n = ReadCEHeader(f, &h0, swp);
    if (n == 0) break;
    *nr += h0.ntransitions;
    free(h0.tegrid);
    free(h0.egrid);
    free(h0.usr_egrid);
    FSEEK(f, h0.length, SEEK_CUR);
  }
  if (*nh == 0 || *nr == 0) {
    FCLOSE(f);
    return -1;
  }
  *h = malloc(sizeof(CE_HEADER)*(*nh));
  *r = malloc(sizeof(CE_RECORD)*(*nr));
  if (mr) {
    *mr = malloc(sizeof(MOD_RECORD)*(*nr));
  }
  FSEEK(f, 0, SEEK_SET);
  n = ReadFHeader(f, &fh, &swp);
  i = 0;
  j = 0;
  while (1) {
    n = ReadCEHeader(f, (*h)+i, swp);
    if (n == 0) break;
    for (m = 0; m < (*h)[i].ntransitions; m++) {
      n = ReadCERecord(f, (*r)+j, swp, (*h)+i);      
      if (n == 0) break;
      if (mr) (*mr)[j].h = (*h)+i;
      j++;
    }
    i++;
  }
  
  if (mr) {
    for (i = 0; i < (*nr); i++) {
      (*mr)[i].m = 0;
      (*mr)[i].i0 = (*r)[i].lower;
      (*mr)[i].i1 = (*r)[i].upper;
      (*mr)[i].r = (*r) + i;
    }
  }

  FCLOSE(f);
  return 0;
}

void ModifyTable(char *fn, char *fn0, char *fn1, char *fnm) {
  int nc, swp, n, nr, nh, i;
  TFILE *f0, *f1;
  MOD_RECORD *mr, *mr0;
  F_HEADER fh;  
  EN_HEADER *eh;
  EN_RECORD *er;
  TR_HEADER *th;
  TR_RECORD *tr;
  CE_HEADER *ch;  
  CE_RECORD *cr;
  void *hp, *rp;
  
  f0 = OpenFileRO(fn0, &fh, &swp);
  if (f0 == NULL) {
    printf("cannot open file %s\n", fn0);
    return;
  }
  
  if (fnm) {
    switch (fh.type) {
    case DB_EN:
      n = ReadENTable(fnm, &nh, &eh, &nr, &er, &mr0);
      hp = eh;
      rp = er;
      break;
    case DB_TR:
      n = ReadTRTable(fnm, &nh, &th, &nr, &tr, &mr0);
      hp = th;
      rp = tr;
      break;
    case DB_CE:
      n = ReadCETable(fnm, &nh, &ch, &nr, &cr, &mr0);
      hp = ch;
      rp = cr;
      break;
    default:
      printf("file %s not of right type %d\n", fnm, fh.type);
      FCLOSE(f0);
      return;
      break;
    }
    qsort(mr0, nr, sizeof(MOD_RECORD), CompareModRecord);
  } else {
    nr = 0; 
    nh = 0;
  }

  nc = ReadModRecord(fn, &mr, nr);
  if (nc <= 0) {
    printf("no records in file %s\n", fn);
    FCLOSE(f0);
    return;
  }
  
  f1 = OpenFile(fn1, &fh);
  
  switch (fh.type) {
  case DB_EN:
    if (nr == 0) {
      for (n = 0; n < nc; n++) {
	mr[n].c /= HARTREE_EV;
      }
    }
    ModifyEN(nc, mr, nr, mr0, f0, f1, &fh, swp);
    break;
  case DB_TR:
    ModifyTR(nc, mr, nr, mr0, f0, f1, &fh, swp);
    break;
  case DB_CE:
    ModifyCE(nc, mr, nr, mr0, f0, f1, &fh, swp);
    for (i = 0; i < nh; i++) {
      free(ch[i].tegrid);
      free(ch[i].egrid);
      free(ch[i].usr_egrid);
    }
    for (i = 0; i < nr; i++) {
      if (cr[i].params) free(cr[i].params);
      free(cr[i].strength);
    }
    break;
  default:
    printf("cannot modify table type %d\n", fh.type);
    break;
  }
    
  if (nh) {
    free(hp);
  }
  if (nr) {
    free(rp);
    free(mr0);
  }
  free(mr);
  CloseFile(f1, &fh);
  FCLOSE(f0);
}

void F77Flush(void) {
  fflush(stdout);
}
FCALLSCSUB0(F77Flush, F77FLUSH, f77flush);
