#include "interpolation.h"
#include "cf77.h"

static char *rcsid="$Id: interpolation.c,v 1.15 2005/04/07 21:27:57 mfgu Exp $";
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

void PrepCECrossHeader(CE_HEADER *h, double *data) {
  double *eusr, *x;
  int m, m1, j;

  eusr = h->usr_egrid;
  m = h->n_usr;
  m1 = m + 1;
  x = data+2+m1;
  data[0] = h->te0*HARTREE_EV;
  for (j = 0; j < m; j++) {
    x[j] = log((h->te0 + eusr[j])/h->te0);
  }
  x[m] = eusr[m-1]/(h->te0+eusr[m-1]);
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

double InterpolateCECross(double e, CE_RECORD *r, CE_HEADER *h, 
			  double *data, double *ratio) {
  double *x, *y, *w;
  int m, m1, n, one;
  double a, b, x0, y0, eth, e0, c, d;
  EN_SRECORD *mem_en_table;
  int mem_en_table_size;

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
  
  if (x0 < x[m-1]) {
    n = 2;
    one = 1;
    UVIP3P(n, m, x, y, one, &x0, &a);
    if (h->msub) {
      UVIP3P(n, m, x, w, one, &x0, &b);
      if (b < 0.0) b = 0.0;
      a *= b;
      *ratio = b;
    }
  } else {
    x0 = e/(data[0] + e);
    y0 = y[m-1];
    if (data[1] > 0) {
      e0 = ((x[m]*data[0]/(1.0-x[m]))+eth)/HARTREE_EV;
      d = 2.0*e0*(1.0+0.5*FINE_STRUCTURE_CONST2*e0);
      c = FINE_STRUCTURE_CONST2*d;
      b = log(0.5*d*HARTREE_EV/eth) - c/(1.0+c);
      y0 /= 1.0 + c;
      y0 -= data[1]*b;
      a = y[m] + (x0-1.0)*(y0-y[m])/(x[m]-1.0);
      e0 = (e + eth)/HARTREE_EV;
      d = 2.0*e0*(1.0+0.5*FINE_STRUCTURE_CONST2*e0);
      c = FINE_STRUCTURE_CONST2*d;
      b = log(0.5*d*HARTREE_EV/eth) - c/(1.0+c);
      a += data[1]*b;
      a *= 1.0 + c;
    } else if (data[1]+1.0 == 1.0) {
      e0 = ((x[m]*data[0]/(1.0-x[m]))+eth)/HARTREE_EV;
      d = 2.0*e0*(1.0+0.5*FINE_STRUCTURE_CONST2*e0);
      c = FINE_STRUCTURE_CONST2*d;
      y0 /= 1.0 + c;
      a = y[m] + (x0-1.0)*(y0-y[m])/(x[m]-1.0);
      e0 = (e + eth)/HARTREE_EV;
      d = 2.0*e0*(1.0+0.5*FINE_STRUCTURE_CONST2*e0);
      c = FINE_STRUCTURE_CONST2*d;
      a *= 1.0 + c;
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

int CECross(char *ifn, char *ofn, int i0, int i1, 
	    int negy, double *egy, int mp) {
  F_HEADER fh;
  FILE *f1, *f2;
  int n, swp;
  CE_HEADER h;
  CE_RECORD r;
  int i, t, m, k;
  double data[2+(1+MAXNUSR)*3], e, cs, a, ratio;
  double eth, a1, cs1, k2, rp, e1, e0;
  EN_SRECORD *mem_en_table;
  int mem_en_table_size;

  mem_en_table = GetMemENTable(&mem_en_table_size);
  
  if (mem_en_table == NULL) {
    printf("Energy table has not been built in memory.\n");
    return -1;
  }
  
  f1 = fopen(ifn, "r");
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
    
  n = ReadFHeader(f1, &fh, &swp);
  if (n == 0) {
    printf("File %s is not in FAC binary format\n", ifn);
    goto DONE;
  }
  
  if (fh.type != DB_CE || fh.nblocks == 0) {
    printf("File %s is not of DB_CE type\n", ifn);
    goto DONE;
  }
   
  while (1) {
    n = ReadCEHeader(f1, &h, swp);
    if (n == 0) break;
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadCERecord(f1, &r, swp, &h);
      if ((r.lower == i0 || i0 < 0) && (r.upper == i1 || i1 < 0)) {
	PrepCECrossHeader(&h, data);
	eth = mem_en_table[r.upper].energy - mem_en_table[r.lower].energy;
	e = eth*HARTREE_EV;
	fprintf(f2, "#%5d\t%2d\t%5d\t%2d\t%11.4E\t%5d\t%d\n",
		r.lower, mem_en_table[r.lower].j,
		r.upper, mem_en_table[r.upper].j,
		e, negy, r.nsub);
	for (k = 0; k < r.nsub; k++) {
	  PrepCECrossRecord(k, &r, &h, data);
	  for (t = 0; t < negy; t++) {
	    if (mp == 0) {
	      e0 = egy[t];
	      e1 = e0 - e;
	    } else {
	      e1 = egy[t];
	      e0 = e1 + e;
	    }
	    if (e1 > 0) {
	      cs = InterpolateCECross(e1, &r, &h, data, &ratio);
	      a = e0/HARTREE_EV;
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
		a1 *= 1.0 + FINE_STRUCTURE_CONST2*k2;
		k2 = cs1/a1;
		rp = k2*(1.0+0.5*FINE_STRUCTURE_CONST2*e0/HARTREE_EV);
		if (rp > 1.0) {
		  k2 /= rp;
		  rp = 1.0;
		}
		cs1 = cs*k2;
		a1 = a*rp;
	      } else {
		k2 = (egy[t]+e)/HARTREE_EV;
		k2 = 2.0*k2*(1.0+0.5*FINE_STRUCTURE_CONST2*k2);
		k2 = 1.0 + FINE_STRUCTURE_CONST2*k2;
		cs1 = cs/k2;
		a1 = a*(1.0+0.5*FINE_STRUCTURE_CONST2*e0/HARTREE_EV)/k2;
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
  fclose(f1);

  if (f2 != stdout) {
    fclose(f2);
  } else {
    fflush(f2);
  }

  return 0;
}

int CEMaxwell(char *ifn, char *ofn, int i0, int i1, 
	      int nt, double *temp) {
  F_HEADER fh;
  FILE *f1, *f2;
  int n, swp;
  CE_HEADER h;
  CE_RECORD r;
  int i, t, m, k, p;
  double data[2+(1+MAXNUSR)*4], e, cs, a, b, c, ratio;
  double xg[15] = {.933078120172818E-01,
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
		   .480260855726858E+02};
  double wg[15] = {.218234885940086E+00,
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
		   .160059490621113E-19};
  EN_SRECORD *mem_en_table;
  int mem_en_table_size;

  mem_en_table = GetMemENTable(&mem_en_table_size);
  
  if (mem_en_table == NULL) {
    printf("Energy table has not been built in memory.\n");
    return -1;
  }
  
  f1 = fopen(ifn, "r");
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
    
  n = ReadFHeader(f1, &fh, &swp);
  if (n == 0) {
    printf("File %s is not in FAC binary format\n", ifn);
    goto DONE;
  }
  
  if (fh.type != DB_CE || fh.nblocks == 0) {
    printf("File %s is not of DB_CE type\n", ifn);
    goto DONE;
  }
   
  while (1) {
    n = ReadCEHeader(f1, &h, swp);
    if (n == 0) break;
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadCERecord(f1, &r, swp, &h);
      if ((r.lower == i0 || i0 < 0) && (r.upper == i1 || i1 < 0)) {
	PrepCECrossHeader(&h, data);
	e = mem_en_table[r.upper].energy - mem_en_table[r.lower].energy;
	e *= HARTREE_EV;
	fprintf(f2, "#%5d\t%2d\t%5d\t%2d\t%11.4E\t%5d\t%d\n",
		r.lower, mem_en_table[r.lower].j,
		r.upper, mem_en_table[r.upper].j,
		e, nt, r.nsub);
	for (k = 0; k < r.nsub; k++) {
	  PrepCECrossRecord(k, &r, &h, data);
	  for (t = 0; t < nt; t++) {
	    cs = 0.0;
	    for (p = 0; p < 15; p++) {
	      a = temp[t]*xg[p];
	      b = a/HARTREE_EV;
	      a = InterpolateCECross(a, &r, &h, data, &ratio);
	      c = 2.0*b*(1.0+0.5*FINE_STRUCTURE_CONST2*b);
	      c = 1.0+FINE_STRUCTURE_CONST2*c;
	      c = c*(1.0+0.5*FINE_STRUCTURE_CONST2*b);
	      c = sqrt(c);
	      cs += wg[p]*a/c;
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
  fclose(f1);

  if (f2 != stdout) {
    fclose(f2);
  } else {
    fflush(f2);
  }

  return 0;
}

int TotalCICross(char *ifn, char *ofn, int ilev, 
		 int negy, double *egy, int imin, int imax) {
  F_HEADER fh;
  FILE *f1, *f2;
  int n, swp;
  CI_HEADER h;
  CI_RECORD r;
  int i, t, nb, m;
  double *c, tc, a, b, x, e;
  EN_SRECORD *mem_en_table;
  int mem_en_table_size;

  mem_en_table = GetMemENTable(&mem_en_table_size);
  
  if (mem_en_table == NULL) {
    printf("Energy table has not been built in memory.\n");
    return -1;
  }

  nb = 0;
  
  f1 = fopen(ifn, "r");
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
  
  n = ReadFHeader(f1, &fh, &swp);
  if (n == 0) {
    printf("File %s is not in FAC binary format\n", ifn);
    goto DONE;
  }

  if (fh.type != DB_CI || fh.nblocks == 0) {
    printf("File %s is not of DB_CI type\n", ifn);
    goto DONE;
  }
  
  c = (double *) malloc(sizeof(double)*negy);
  for (i = 0; i < negy; i++) {
    c[i] = 0.0;
  }

  if (imin < 0) imin = 0;
  if (imax < 0) imax = mem_en_table_size - 1;

  while (1) {
    n = ReadCIHeader(f1, &h, swp);
    if (n == 0) break;
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadCIRecord(f1, &r, swp, &h);
      if (n == 0) break;
      if (r.b != ilev) continue;
      if (r.f < imin || r.f > imax) continue;
      e = mem_en_table[r.f].energy - mem_en_table[r.b].energy; 
      
      for (t = 0; t < negy; t++) {
	if (egy[t] < e) continue;
	x = egy[t]/e;
	a = 1.0/x;
	b = 1.0 - a;
	tc = r.params[0]*log(x) + r.params[1]*b*b;
	tc += r.params[2]*a*b + r.params[3]*a*a*b;
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
  fclose(f1);

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
  FILE *f1, *f2;
  int n, swp;
  RR_HEADER h;
  RR_RECORD r;
  int i, t, nb, m;
  float e, eph, ee, phi;
  double *xusr, *dstrength, *c, tc, emax;
  double x, y, a;
  int np=3, one=1, nele;
  EN_SRECORD *mem_en_table;
  int mem_en_table_size;

  mem_en_table = GetMemENTable(&mem_en_table_size);
  
  if (mem_en_table == NULL) {
    printf("Energy table has not been built in memory.\n");
    return -1;
  }

  nb = 0;

  f1 = fopen(ifn, "r");
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

  n = ReadFHeader(f1, &fh, &swp);
  if (n == 0) {
    printf("File %s is not in FAC binary format\n", ifn);
    goto DONE;
  }

  if (fh.type != DB_RR || fh.nblocks == 0) {
    printf("File %s is not of DB_RR type\n", ifn);
    goto DONE;
  }
  
  c = (double *) malloc(sizeof(double)*negy);
  for (i = 0; i < negy; i++) {
    c[i] = 0.0;
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
	a = FINE_STRUCTURE_CONST2*ee;
	a = (1.0+a)/(1+0.5*a);
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
	phi = a*2.0*PI*FINE_STRUCTURE_CONST*tc*AREA_AU20;	
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
  fclose(f1);

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
  FILE *f1, *f2;
  int n, swp;
  RR_HEADER h;
  RR_RECORD r;
  int i, t, nb, m;
  float e, eph, ee, phi, rr;
  double *xusr, *dstrength, *c, tc, emax;
  double x, y, a;
  int np=3, one=1, nele;
  EN_SRECORD *mem_en_table;
  int mem_en_table_size;

  mem_en_table = GetMemENTable(&mem_en_table_size);

  if (mem_en_table == NULL) {
    printf("Energy table has not been built in memory.\n");
    return -1;
  }

  nb = 0;

  f1 = fopen(ifn, "r");
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

  n = ReadFHeader(f1, &fh, &swp);
  if (n == 0) {
    printf("File %s is not in FAC binary format\n", ifn);
    goto DONE;
  }

  if (fh.type != DB_RR || fh.nblocks == 0) {
    printf("File %s is not of DB_RR type\n", ifn);
    goto DONE;
  }
  
  c = (double *) malloc(sizeof(double)*negy);
  for (i = 0; i < negy; i++) {
    c[i] = 0.0;
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
	a = FINE_STRUCTURE_CONST2*ee;
	a = (1.0+a)/(1+0.5*a);
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
	rr = a * phi * pow(FINE_STRUCTURE_CONST*eph, 2) / (2.0*ee);
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
  fclose(f1);
  if (f2 != stdout) {
    fclose(f2);
  } else {
    fflush(f2);
  }

  return nb;
}
