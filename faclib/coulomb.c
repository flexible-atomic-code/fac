#include "coulomb.h"

/*************************************************************
  Implementation for module "coulomb". 
  This module calculates quatities related to the H-like ions.

  Author: M. F. Gu, mfgu@space.mit.edu
**************************************************************/

static char *rcsid="$Id: coulomb.c,v 1.21 2003/01/13 02:57:41 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

static int n_hydrogenic;
static int kl_hydrogenic;
static int n_hydrogenic_max;
static int kl_hydrogenic_max;
static ARRAY *dipole_array;

static int _ncb = 0;
static int _cbindex[CBMULTIPOLES];
static double *_cb[MAXNE][MAXNTE][MAXNE][MAXNCB];
static double *_dwork = NULL;
static int _nm_min = 100;
static int _nm_max = 50000;
static int _nm_factor = 100;
static int _nm = 0;


double argam_(double *x, double *y);
void acofz1_(double *, double *, int *, int *, double *, 
	     double *, int *, int *);
void uvip3p_(int *np, int *ndp, double *x, double *y, 
	     int *n, double *xi, double *yi);

void SetHydrogenicNL(int n, int kl, int nm, int klm) {
  if (n > 0) n_hydrogenic = n;
  else n_hydrogenic = NHYDROGEN;
  if (kl >= 0) kl_hydrogenic = kl;
  else kl_hydrogenic = LHYDROGEN;
  if (nm > 0) n_hydrogenic_max = nm;
  else n_hydrogenic_max = NHYDROGENMAX;
  if (klm >= 0) kl_hydrogenic_max = klm;
  else kl_hydrogenic_max = LHYDROGENMAX;
  if (n_hydrogenic_max < n_hydrogenic) 
    n_hydrogenic_max = n_hydrogenic;
  if (kl_hydrogenic_max < kl_hydrogenic) 
    kl_hydrogenic_max = kl_hydrogenic;
}

void GetHydrogenicNL(int *n, int *kl, int *nm, int *klm) {
  if (n) *n = n_hydrogenic;
  if (kl) *kl = kl_hydrogenic;
  if (nm) *nm = n_hydrogenic_max;
  if (klm) *klm = kl_hydrogenic_max;
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

double HydrogenicDipole(double z, int n0, int kl0, int n1, int kl1) {
  double anc, am;
  double z0 = 1.0;
  int nmax = 512;
  double ac[1024];
  static iopt = 2;
  double **qk, *t;
  int n, i;
  
  if (n1 > 512) return 0.0;
  if (n0 >= n1) {
    return 0.0;
  } 
  if (kl1 != kl0 + 1 && kl1 != kl0 - 1) {
    return 0.0;
  }
  qk = (double **) ArraySet(dipole_array, n1, NULL);
  if (*qk == NULL) {
    *qk = (double *) malloc(sizeof(double)*n1*(n1-1));
    t = *qk;
    am = 100.0;
    if (iopt == 2) {
      acofz1_(&z0, &am, &nmax, &n0, ac, &anc, &n1, &iopt);
    }
    iopt = 1;
    for (n = 1; n < n1; n++) {
      acofz1_(&z0, &am, &n1, &n, ac, &anc, &n1, &iopt);
      for (i = 0; i < n; i++) {
	*(t++) = ac[i];
      }
      for (i = 0; i < n; i++) {
	*(t++) = ac[i+n1];
      }
    }
  }

  t = (*qk) + n0*(n0-1);
  if (kl1 == kl0 - 1) {
    return t[kl1]/z;
  } else {
    return t[n0 + kl0]/z;
  }
}

double TRRateHydrogenic(double z, int n0, int kl0, int n1, int kl1, int s) {
  double e, c, g, al;
  double factor, e2, r;

  c = HydrogenicDipole(z, n0, kl0, n1, kl1);
  if (c == 0) return c;
  e2 = z*z*(1.0/(n0*n0) - 1.0/(n1*n1));
  al = (double) Max(kl0, kl1);
  if (s == 0) {
    factor = 2.6775015E9*e2*e2*e2;
    r = (al/(2.0*kl1 + 1.0))*c*c*factor;
  } else if (s == 1) {
    factor = 2.6775015E9*e2;
    r = (4.0*al/(2.0*kl1 + 1.0))*c*c*factor;
    c = 1.0/(2.0*pow(FINE_STRUCTURE_CONST, 3.0));
    r /= RATE_AU;
    g = 2.0*(2.0*kl1+1.0);
    r *= g*c;
  } else {
    r = c;
  }
  return r;
}

double HydrogenicExpectation(double z, int m, int n, int kl) {
  double r, n2, k, e;
  
  k = kl*(kl+1.0);
  n2 = n*n;
  e = -0.5*z*z/n2;
  r = 0.0;
  
  switch (m) {
  case -3:
    r = 2*z*z/(n2*n*(2.0*kl+1.0));
    r *= z/k; 
    break;
  case -2:
    r = 2*z*z/(n2*n*(2.0*kl+1.0));
    break;
  case -1:
    r = z/n2;
    break;
  case 0:
    r = 1.0;
    break;
  case 1:
    r = (3.0*n2 - k)/(2.0*z);
    break;
  case 2:
    r = (n2/(2*z*z))*(5*n2 + 1.0 - 3.0*k);
    break;
  case 3:
    r = (n2/(8*z*z*z))*(5*n2*(7*n2+5)+4*k*(k-10*n2-2));
    break;
  default:
    break;
  }

  return r;
}

double HydrogenicSelfEnergy(double z, int n, int k) {
  double zd[12] = {1.0, 10.0, 20.0, 30.0, 40.0, 50.0, 
		   60.0, 70.0, 80.0, 90.0, 100.0, 110};
  double sd[16][12] = {{10.3168,4.6540,3.2460,2.5519,
			2.1351,1.8644,1.6838,1.5675,
			1.5032,1.4880,1.5317,1.6614},
		       {10.5468,  4.8930,   3.5063,   2.8391, 
			2.4550,  2.2244,   2.0948,   2.0435,
			2.0650,  2.1690,   2.3870,   2.7980},
		       {-0.1264, -0.1145,  -0.0922,  -0.0641,
			-0.0308,  0.0082,   0.0549,   0.1129,
			0.1884,  0.2934,   0.4530,   0.7250 },
		       {0.1235,  0.1303,   0.1436,   0.1604,
			0.1794,  0.1999,   0.2215,   0.2440,
			0.2671,  0.2906,   0.3141,   0.3367},
		       {10.5000,  4.9524,   3.5633,   2.8940,
			2.5083,  2.2757,   2.1431,   2.0874,
			2.1018,  2.1935,   2.3897,   2.7609},
		       {-0.1300, -0.1021,  -0.0760,  -0.0430,
			0.0041,  0.0414,   0.0956,   0.1623,
			0.2483,  0.3660,   0.5408,   0.8322},
		       {0.1250,  0.1421,   0.1572,   0.1761,
			0.1977,  0.2214,   0.2470,   0.2745,
			0.3038,  0.3350,   0.3679,   0.4020},
		       {-0.0440, -0.0428,  -0.0420,  -0.0410,
			-0.0396, -0.0378,  -0.0353,  -0.0321,
			-0.0279, -0.0225,  -0.0154,  -0.0062},
		       {10.5000,  4.9749,   3.5834,   2.9110,
			2.5215,  2.2842,   2.1455,   2.0814,
			2.0840,  2.1582,   2.3262,   2.6484},
		       {-0.1200, -0.0963,  -0.0690,  -0.0344,
			0.0064,  0.0538,   0.1098,   0.1780,
			0.2649,  0.3819,   0.5525,   0.8311},
		       {0.1250,  0.1477,   0.1630,   0.1827,
			0.2052,  0.2299,   0.2568,   0.2858,
			0.3170,  0.3507,   0.3868,   0.4247},
		       {-0.0410, -0.0403,  -0.0399,  -0.0387,
			-0.0371, -0.0348,  -0.0317,  -0.0276,
			-0.0222, -0.0149,  -0.0053,   0.0074},
		       {10.5000,  4.9858,   3.5923,   2.9173,
			2.5246,  2.2833,   2.1395,   2.0686,
			2.0619,  2.1225,   2.2696,   2.5566},
		       {-0.1200, -0.0933,  -0.0652,  -0.0299,
			0.0116,  0.0597,   0.1161,   0.1843,
			0.2703,  0.3848,   0.5497,   0.8150},
		       {0.1300,  0.1502,   0.1662,   0.1861,
			0.2089,  0.2341,   0.2614,   0.2910,
			0.3229,  0.3574,   0.3946,   0.4338},
		       {-0.0405, -0.0396,  -0.0387,  -0.0374,
			-0.0356, -0.0331,  -0.0297,  -0.0252,
			-0.0190, -0.0108,   0.0001,   0.0145}};
  int nd[5] = {0, 1, 4, 8, 12};
  int id, np = 3, nx = 12, m = 1;
  double c;
  double r;

  if (n <= 5) {
    id = nd[n-1];
  } else {
    id = nd[4];
  }
  switch (k) {
  case -1:
    break;
  case 1:
    id += 1;
    break;
  case -2:
    id += 2;
    break;
  case 2:
    id += 3;
    break;
  default:
    return 0.0;
  }

  uvip3p_(&np, &nx, zd, sd[id], &m, &z, &r);
  
  c = FINE_STRUCTURE_CONST*z;
  c = c*c*c*c;
  c /= PI*n*n*n*FINE_STRUCTURE_CONST;
  r *= c;
  
  return r;
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

  i = n-1;
  r = w3[i];
  tn = r/(1.0-r);
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
      i = k;
    }
  }
  return 0;
}

int PrepCoulombBethe(int ne2, int nte, int ne1, double z,
		     double *e2, double *te, double *e1,
		     int nkl, double *kl, 
		     int etype, int ltype, int mode) {
  double ee0, ee1, z2, a, b, c, d0, d1, d2;
  int i, k, n, ie1, ie2, ite, i1p, i1m, nm;
  double *w0, *w1, *w2, *w3, *w4, *tcb, k0, k1, eta, eta2;

  if (mode) ltype = 1;

  if (ne2 > MAXNE || ne1 > MAXNE || nte > MAXNTE) {
    printf("Array multipoles not large enough in CoulombMultipoles\n");
    exit(1);
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
  double e1[] = {6.0986E1,5.296E2};
  double te[] = {1.68686E-1,8.434314E-1};
  
  ne2 = 1;
  ne1 = 2;
  nte = 2;
  z = 17;

  f = fopen(s, "w");
  for (i = 0; i < 50; i++) {
    kl[i] = i;
  }
  j = 2;
  for (i = 50; i < M; i++) {
    kl[i] = kl[i-1] + j;
    j += 1;
  }
  
  PrepCoulombBethe(ne2, nte, ne1, z, e2, te, e1, M, kl, 1, 0, 0);
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
    
int InitCoulomb(void) {
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

  SetHydrogenicNL(-1, -1, -1, -1);

  dipole_array = (ARRAY *) malloc(sizeof(ARRAY));
  ArrayInit(dipole_array, sizeof(double *), 10);

  return 0;
}
