#include "coulomb.h"

static char *rcsid="$Id: coulomb.c,v 1.5 2001/09/14 13:16:59 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

static int _ncb = 0;
static int _cbindex[CBMULTIPOLES];
static double *_cb[MAXNE][MAXNTE][MAXNE][MAXNCB];
static double *_dwork = NULL;
static int _nm_min = 100;
static int _nm_max = 5000;
static int _nm_factor = 100;
static int _nm = 0;

int AddPW(int *nkl0, double *kl, double *logkl, 
	  int maxkl, int n, int step) {
  int i;
  for (i = *nkl0; i < n+(*nkl0); i++) {
    if (i >= MAXNKL) {
      printf("Maximum partial wave grid points reached: "); 
      printf("%d > %d in constructing grid\n",  i, MAXNKL);
      abort();
    }
    kl[i] = kl[i-1] + step;
    logkl[i] = log(kl[i]);
    if ((int)(kl[i]) > maxkl) break;
  }
  (*nkl0) = i;
  return 0;
}

int SetPWGrid(int *nkl0, double *kl, double *logkl, 
	      int maxkl, int *ns, int *n, int *step) {
  int i, m, k, j;

  if ((*ns) > 0) {
    for (i = 0; i < (*ns); i++) {
      AddPW(nkl0, kl, logkl, maxkl, n[i], step[i]);
    }
    k = step[(*ns)-1]*2;
    j = 2;
  } else {
    (*ns) = -(*ns);
    if ((*ns) == 0) (*ns) = 8;
    AddPW(nkl0, kl, logkl, maxkl, (*ns), 1);
    k = 2;
    j = 2;
  }   

  m = kl[(*nkl0)-1];
  while (m+k <= maxkl) {
    AddPW(nkl0, kl, logkl, maxkl, j, k);
    m = kl[(*nkl0)-1];
    if (k < 50) k *= 2;
    else k += 50;
  }
  kl[(*nkl0)] = maxkl+1;
  return (*nkl0);
}

int SetTEGridDetail(double *te, double *logte, int n, double *x) {
  int i;
  
  if (n > MAXNTE) {
    printf("Max # of grid points reached \n");
    return -1;
  }

  for (i = 0; i < n; i++) {
    te[i] = x[i];
    logte[i] = log(te[i]);
  }
  return n;
}

int SetTEGrid(double *te, double *logte, int n, double emin, double emax) {
  int i;
  double del;

  if (n < 1) {
    te[0] = -1.0;
    return 0;
  }

  if (emin < 0.0) {
    te[0] = emin;
    return n;
  }

  if (n > MAXNTE) {
    printf("Max # of grid points reached \n");
    return -1;
  }

  if (n == 1) {
    te[0] = emin;
    logte[0] = log(emin);
    return n;
  }

  if (n == 2) {
    te[0] = emin;
    te[1] = emax;
    logte[0] = log(emin);
    logte[1] = log(emax);
    return n;
  }

  if (emax < emin) {
    printf("emin must > 0 and emax < emin\n");
    return -1;
  }
  
  del = emax - emin;
  del /= n-1.0;
  te[0] = emin;
  logte[0] = log(emin);
  for (i = 1; i < n; i++) {
    te[i] = te[i-1] + del;
    logte[i] = log(te[i]);
  }
  
  return n;
}
  
int SetEGridDetail(double *e, double *log_e, int n, double *xg) {
  int i;
  
  for (i = 0; i < n; i++) {
    e[i] = xg[i];
    log_e[i] = log(e[i]);
  }

  return n;
}

int SetEGrid(double *e, double *log_e, 
	     int n, double emin, double emax, double eth) {
  double del, et;
  int i;

  if (n < 1) {
    e[0] = -1.0;
    return 0;
  }
  if (emin < 0.0) {
    e[0] = emin;
    return 0;
  }

  if (emax < emin) {
    printf("emin must > 0 and emax < emin\n");
    return -1;
  }

  et = fabs(eth);
  if (et > 1E-30) {
    emin += et;
    emax += et;
  }
  
  e[0] = emin;
  log_e[0] = log(emin);
  e[n-1] = emax;
  log_e[n-1] = log(emax);
  del = (log_e[n-1] - log_e[0])/(n-1.0);
  del = exp(del);
  for (i = 1; i < n-1; i++) {
    e[i] = e[i-1]*del;
    log_e[i] = log(e[i]);
  }

  if (eth > 1E-30) {
    for (i = 0; i < n; i++) {
      e[i] -= eth;
    }
  }
  return n;
}

double CoulombPhaseShift(double z, double e, int kappa) {
  double phase, r, y, ke, a, b1, b2;

  a = FINE_STRUCTURE_CONST2 * e;
  ke = sqrt(2.0*e*(1.0 + 0.5*a));
  a += 1.0;
  y = a*z/ke;

  r = kappa;
  b1 = y/(fabs(r)*a);
  r = sqrt(r*r - FINE_STRUCTURE_CONST2*z*z);
  b2 = y/r;

  if (kappa < 0) {
    phase = 0.5*(atan(b1) - atan(b2));
  } else {
    phase = -0.5*(atan(b1) + atan(b2) + PI);
  }

  phase -= argam_(&r, &y);
  phase += (1.0 - r)*0.5*PI;

  return phase;
}

double *GetCoulombBethe(int ie2, int ite, int ie1, int m, int k) {
  return _cb[ie2][ite][ie1][_cbindex[m]+k];
}

double GetCoulombBetheAsymptotic(double te, double e1) {
  double r;

  r = e1/(e1+te);
  r = r/(1.0-r);
  
  return r;
}

int PrepCBIndex(int mode) {
  int i;
  _cbindex[0] = 0;
  for (i = 1; i < CBMULTIPOLES; i++) {
    if (mode) _cbindex[i] = _cbindex[i-1] + i+1;
    else _cbindex[i] = _cbindex[i-1] + 2;
  }
  if (mode) _ncb = _cbindex[i-1] + i+1;
  else _ncb = _cbindex[i-1] + 2;

  return 0;
}

/* calculates the Coulomb Bethe contributions from L to Inf. */
int CoulombBetheTail(int n, double *w3, int nkl, double *kl, double *tcb) {
  double r, a, b, tn;
  int i, j, k, p;

  r = w3[n-1];
  tn = r/(1.0-r);
  for (i = n-2; i > 0 ; i--) {
    if (fabs(1.0 - w3[i]/r) > 1E-3) break;
  }	
  for (j = nkl - 1; j > 0; j--) {
    k = kl[j];
    if (k >= i) {
      tcb[j] = tn;
    } else {
      a = 1.0;
      b = 0.0;
      for (p = k; p < i; p++) {
	a *= w3[p];
	b += a;
      }
      b += a*tn;
      tcb[j] = b;
      tn = b;
    }
  }
}

int PrepCoulombBethe(int ne2, int nte, int ne1, double z,
		     double *e2, double *te, double *e1,
		     int nkl, double *kl, 
		     int etype, int ltype, int mode) {
  double xi, ee0, ee1, z2, a, b, c, d0, d1, d2;
  int i, j, k, n, ie1, ie2, ite, i1p, i1m, nm;
  double *w0, *w1, *w2, *w3, *w4, *tcb, k0, k1, r, eta, eta2;

  if (mode) ltype = 1;

  if (ne2 > MAXNE || ne1 > MAXNE || nte > MAXNTE) {
    printf("Array multipoles not large enough in CoulombMultipoles\n");
    abort();
  }
  PrepCBIndex(mode);
  _nm = _nm_factor*(e1[ne1-1]/(te[0]+e2[0]));
  if (_nm > _nm_max) _nm = _nm_max;
  _nm = Max(_nm, kl[nkl-1] + _nm_min);

  free(_dwork);
  _dwork = malloc(sizeof(double)*_nm*4);
  for (ie2 = 0; ie2 < MAXNE; ie2++) {
    for (ite = 0; ite < MAXNTE; ite++) {
      for (ie1 = 0; ie1 < MAXNE; ie1++) {
	for (i = 0; i < _ncb; i++) {
	  free(_cb[ie2][ite][ie1][i]);
	  _cb[ie2][ite][ie1][i] = NULL;
	}
      }
    }
  }

  w0 = _dwork;
  w1 = w0 + _nm;
  w2 = w1 + _nm;
  w3 = w2 + _nm;
  nm = _nm;
  n = nm - 1;
  
  z2 = z*z;
  if (etype == 0) {
    w4 = w1;
  } else {
    w4 = w0;
  }
  for (ie1 = 0; ie1 < ne1; ie1++) {
    if (etype == 0) {
      ee0 = e1[ie1];
      k0 = sqrt(2.0*ee0); 
      for (i = 0; i < nm; i++) {
	w0[i] = sqrt(z2 + k0*k0*i*i);
      }
    } else {
      ee1 = e1[ie1];
      k1 = sqrt(2.0*ee1); 
      for (i = 0; i < nm; i++) {
	w1[i] = sqrt(z2 + k1*k1*i*i);
      }
    }
    for (ie2 = 0; ie2 < ne2; ie2++) {
      for (ite = 0; ite < nte; ite++) {
	for (k = 0; k < _ncb; k++) {
	  _cb[ie2][ite][ie1][k] = malloc(sizeof(double)*nkl);
	}
	if (etype == 0) {
	  ee1 = ee0 - e2[ie2] - te[ite];  
	  k1 = sqrt(2.0*ee1);
	  for (i = 0; i < nm; i++) {
	    w1[i] = sqrt(z2 + k1*k1*i*i);	  
	  }
	} else {
	  ee0 = ee1 + e2[ie2] + te[ite];  
	  k0 = sqrt(2.0*ee0);
	  for (i = 0; i < nm; i++) {
	    w0[i] = sqrt(z2 + k0*k0*i*i);	  
	  }
	}
	/* monopole radial integrals using recursion */
	i1p = n;
	i = n-1;
	i1m = n-2;
	c = w0[i1p]*w1[i1p];
	b = (1.0/(2*i) + 1.0) * 2.0 * (z2 + i*(i+1.0)*(ee0 + ee1));
	a = (1.0 + 1.0/i) * w0[i]*w1[i];
	d0 = 0.5*k0*k1/(ee0-ee1);
	d1 = d0/(2.0*i+1.0);
	d1 = (w1[i] - w0[i]*d1)/(w0[i] - w1[i]*d1);
	d2 = d0/(2.0*i1p+1.0);
	d2 = (w1[i1p] - w0[i1p]*d2)/(w0[i1p] - w1[i1p]*d2);
	w2[i] = (b-sqrt(b*b-4.0*a*c*(d2/d1)))/(2.0*c);
	
	for (;;) {
	  w2[i1m] = a/(b-c*w2[i]);
	  if (i == 1) break;
	  i1p = i;
	  i = i1m;
	  i1m--;	  
	  c = w0[i1p]*w1[i1p];
	  b = (1.0 + 1.0/(2*i)) * 2.0 * (z2 + i*(i+1.0)*(ee0 + ee1));
	  a = (1.0 + 1.0/i) * w0[i]*w1[i];
	}
	
	/* Coulomb Bethe contribution for monopole */
	for (i = 0; i < n; i++) {
	  w3[i] = w2[i]*w2[i] * (i+1.0)/i;
	}
	tcb = GetCoulombBethe(ie2, ite, ie1, 0, 0);
	for (i = 1; i < nkl; i++) {
	  k = kl[i];
	  tcb[i] = w3[k];
	}
	tcb = GetCoulombBethe(ie2, ite, ie1, 0, 1);
	CoulombBetheTail(n, w3, nkl, kl, tcb);

	/* Coulomb Bethe for dipoles */
	i1m = 0;
	i = 1;
	i1p = 2;
	b = w2[0];
	for (; i < n; ) {
	  a = w2[i];
	  d1 = ((double)i)/i1p;
	  if (ltype == 0) {
	    d2 = (w0[i1p] - w1[i1p]*a) / (w0[i]/b - w1[i]);
	    w3[i] = d1 * d2;
	    d2 = (w1[i1p] - w0[i1p]*a) / (w0[i]/b - w1[i]);
	    w2[i] = d1 * d2;
	  } else {
	    d2 = (w1[i1p] - w0[i1p]*a) / (w1[i]/b - w0[i]);
	    w3[i] = d1 * d2;
	    d2 = (w0[i1p] - w1[i1p]*a) / (w1[i]/b - w0[i]);
	    w2[i] = d1 * d2;
	  }
	  b = a;
	  i1m = i;
	  i = i1p;
	  i1p++;
	}		
	for (i = 1; i < n-1; i++) {
	  d1 = (i+1.0)/i;
	  d2 = (i+2.0)/(i+1.0);
	  d2 = (1.0 + d2*w2[i+1]*w2[i+1])/(1.0 + d1*w2[i]*w2[i]);
	  w3[i] = d1*d2*w3[i]*w3[i];
	}
	w3[i] = w3[i-1];
	tcb = GetCoulombBethe(ie2, ite, ie1, 1, 0); 
	for (i = 1; i < nkl; i++) {
	  k = kl[i];
	  tcb[i] = w3[k];
	} 
	if (mode == 0) {
	  tcb = GetCoulombBethe(ie2, ite, ie1, 1, 1);
	  CoulombBetheTail(n, w3, nkl, kl, tcb);
	} else {
	  d2 = 0.0;
	  c = 1.0;
	  eta = z/k0;
	  eta2 = eta*eta;
	  for (i = 1; i < n; i++) {
	    a = w2[i]*w2[i];
	    b = 1.0 + 1.0/i;
	    d1 = 1.0 + b*a;
	    d0 = 1.0 + a;
	    a = i*i + eta2;
	    b = (i+1.0)*(i+1.0) + eta2;
	    a = i*(i+1.0)/sqrt(a*b);
	    b = 1.0 - eta2/(i*(i+1.0));
	    d0 -= 2.0*w2[i]*a*b;
	    w2[i] = d0/d1;
	    if (fabs(1.0 - c/d2) > 1E-3) {
	      d2 = c;
	      k = 2*i;
	      c = AngularMSub(k, k-2, k-2, -2);
	      c = c/i;
	    } else {
	      d2 = c;
	      c = 0.5 + (c-0.5)*(i/(i+1.0));
	    }
	    w2[i] *= c;
	  }
	  
	  for (i = 1; i < n-1; i++) {
	    w4[i] = w3[i] * w2[i+1]/w2[i];
	  }
	  tcb = GetCoulombBethe(ie2, ite, ie1, 1, 1);	
	  CoulombBetheTail(n, w4, nkl, kl, tcb);	
	  for (i = 1; i < n-1; i++) {
	    w4[i] = w3[i] * (1.0-w2[i+1])/(1.0-w2[i]);
	  }
	  tcb = GetCoulombBethe(ie2, ite, ie1, 1, 2);	
	  CoulombBetheTail(n, w4, nkl, kl, tcb);	
	}
      }
    }
  } 
  return 0;
}
 
double AngularMSub(int lf, int li1, int li2, int q) {
  int mi, mf, ji1, ji2, jf;
  double a, b, c, d, as, bs, cs, ds;

  a = sqrt((li1+1.0)*(li2+1.0));
  as = 0.0;
  for (jf = lf - 1; jf <= lf + 1; jf += 2) {
    bs = 0.0;
    for (ji1 = li1 - 1; ji1 <= li1 + 1; ji1 += 2) {
      b = sqrt(ji1+1.0);
      b *= ReducedCL(ji1, 2, jf);
      cs = 0.0;
      for (ji2 = li2 - 1; ji2 <= li2 + 1; ji2 += 2) {
	c = sqrt(ji2+1.0);
	c *= ReducedCL(ji2, 2, jf);
	ds = 0.0;
	for (mi = -1; mi <= 1; mi += 2) {
	  d = W3j(ji1, 1, li1, -mi, mi, 0);
	  d *= W3j(ji2, 1, li2, -mi, mi, 0);
	  mf = q + mi;
	  d *= W3j(ji1, 2, jf, -mi, -q, mf);
	  d *= W3j(ji2, 2, jf, -mi, -q, mf);
	  ds += d;
	}
	cs += c*ds;
      }
      bs += b*cs;
    }
    as += bs;
  }
  a *= as;
  return a;
}

int TestCoulomb(char *s) {
#define M 80

  double kl[M], z;
  int i, j;
  FILE *f;
  int nte, ne1, ne2, ie1, ite;
  double e2[] = {0.0};
  double e1[] = {600.0, 1100.0};
  double te[] = {50.0, 100.0};
  
  ne2 = 1;
  ne1 = 2;
  nte = 2;
  z = 10;

  f = fopen(s, "w");
  for (i = 0; i < 50; i++) {
    kl[i] = i;
  }
  j = 2;
  for (i = 50; i < M; i++) {
    kl[i] = kl[i-1] + j;
    j += 1;
  }
  
  PrepCoulombBethe(ne2, nte, ne1, z, e2, te, e1, M, kl, 0, 0, 0);
  for (ite = 0; ite < nte; ite++) {
    for (ie1 = 0; ie1 < ne1; ie1++) {
      fprintf(f, "## TE = %10.3E, E1 = %10.3E\n", te[ite], e1[ie1]);
      for (i = 0; i < M; i++) {
	fprintf(f, "%d %d ", i, (int)kl[i]);
	for (j = 0; j < _ncb; j++) {
	  fprintf(f, "%12.5E ", _cb[0][ite][ie1][j][i]);
	}
	fprintf(f, "\n");
      }
      fprintf(f, "\n\n");
    }
  }

  return 0;
#undef M
}
    
int InitCoulomb() {
  int i, ie1, ie2, ite;

  for (ie2 = 0; ie2 < MAXNE; ie2++) {
    for (ite = 0; ite < MAXNTE; ite++) {
      for (ie1 = 0; ie1 < MAXNE; ie1++) {
	for (i = 0; i < MAXNCB; i++) {
	  _cb[ie2][ite][ie1][i] = NULL;
	}
      }
    }
  }

  return 0;
}
