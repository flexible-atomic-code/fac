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

#include "coulomb.h"
#include "cf77.h"

/*************************************************************
  Implementation for module "coulomb". 
  This module calculates quatities related to the H-like ions.

  Author: M. F. Gu, mfgu@stanford.edu
**************************************************************/

static char *rcsid="$Id: coulomb.c,v 1.31 2006/08/04 07:43:53 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

static int n_hydrogenic;
static int kl_hydrogenic;
static int n_hydrogenic_max;
static int kl_hydrogenic_max;
static ARRAY *dipole_array;

static int _cbindex[CBMULT][CBMULT+1];
static double *_cb[MAXNE][MAXNTE][MAXNE][MAXNCB];

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

double HydrogenicDipole(double z, int n0, int kl0, int n1, int kl1) {
  double anc, am;
  double z0 = 1.0;
  int nmax = 512;
  double ac[1024];
  static int iopt = 2;
  double **qk, *t;
  int n, i;
  
  if (n1 > 512) return 0.0;
  if (n0 >= n1) {
    return 0.0;
  } 
  if (kl1 != kl0 + 1 && kl1 != kl0 - 1) {
    return 0.0;
  }
  qk = (double **) ArraySet(dipole_array, n1, NULL, InitPointerData);
  if (*qk == NULL) {
    *qk = (double *) malloc(sizeof(double)*n1*(n1-1));
    t = *qk;
    am = 100.0;
    if (iopt == 2) {
      ACOFZ1(z0, am, nmax, n0, ac, &anc, n1, iopt);
    }
    iopt = 1;
    for (n = 1; n < n1; n++) {
      ACOFZ1(z0, am, n1, n, ac, &anc, n1, iopt);
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
  double c, g, al;
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

  UVIP3P(np, nx, zd, sd[id], m, &z, &r);
  
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

  phase -= ARGAM(r, y);
  phase += (1.0 - r)*0.5*PI;

  return phase;
}

double *GetCoulombBethe(int ie2, int ite, int ie1, int t, int q) {
  return _cb[ie2][ite][ie1][_cbindex[t-1][q]];
}

double GetCoulombBetheAsymptotic(double te, double e1) {
  double r;

  r = e1/(e1+te);
  r = r/(1.0-r);
  
  return r;
}

void PrepCBIndex(void) {
  int t, q, i;

  i = 0;
  for (t = 1; t <= CBMULT; t++) {
    for (q = 0; q <= t; q++) {
      _cbindex[t-1][q] = i;
      i++;
    }
  }
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

double PartialOmegaMSub(int k1, double *r, int k, int m, double z, double e) {
  int k0, k0p, m0, j0, j0p, j1, i, ip, k1m;
  int km, m2, k0m, k0pm, kappa0, kappa0p;
  double p0, p0p, a, w1, w2, w3, w4, w5, w6, x;

  km = k*2;
  m2 = m*2;
  x = 0.0;
  for (i = 0; i <= k; i++) {
    k0 = k1 + 2*i - k;
    k0m = k0*2;    
    for (j0 = k0m - 1; j0 <= k0m + 1; j0 += 2) {
      kappa0 = GetKappaFromJL(j0, k0m);
      p0 = CoulombPhaseShift(z, e, kappa0);
      for (ip = 0; ip <= k; ip++) {
	k0p = k1 + 2*ip - k;
	k0pm = k0p*2;
	for (j0p = k0pm - 1; j0p <= k0pm + 1; j0p += 2) {
	  kappa0p = GetKappaFromJL(j0p, k0pm);
	  p0p = CoulombPhaseShift(z, e, kappa0p);	  
	  a = sqrt((k0m+1.0)*(k0pm+1.0)*(j0+1.0)*(j0p+1.0))*cos(p0p-p0);
	  k1m = k1*2;
	  for (j1 = k1m - 1; j1 <= k1m + 1; j1 += 2) {
	    w5 = ReducedCL(j0, km, j1);
	    w6 = ReducedCL(j0p, km, j1);
	    for (m0 = -1; m0 <= 1; m0 += 2) {
	      w1 = W3j(j0, 1, k0m, -m0, m0, 0);
	      w2 = W3j(j0p, 1, k0pm, -m0, m0, 0);
	      w3 = W3j(j0, km, j1, -m0, -m2, m0+m2);
	      w4 = W3j(j0p, km, j1, -m0, -m2, m0+m2);
	      x += w1*w2*w3*w4*w5*w6*a*r[i]*r[ip];
	    }
	  }
	}
      }
    }
  }

  return x;
}

double PartialOmega(int k0, double *r, int k) {
  int i, k1, g, g2, i0, i1, i2, i3;
  double a, b, c, x;

  x = 0.0;
  for (i = 0; i <= k; i++) {
    k1 = k0 + i*2 - k;
    a = 2.0*(2.0*k0+1.0)*(2.0*k1+1.0);
    g2 = k0 + k1 + k;
    g = g2/2;
    i0 = g2 - 2*k0;
    i1 = g2 - 2*k1;
    i2 = g2 - 2*k;
    if (i2 < 0) i2 = 0;
    i3 = g2 + 1;
    b = LnFactorial(i0)+LnFactorial(i1)+LnFactorial(i2)-LnFactorial(i3);
    i0 = g - k0;
    i1 = g - k1;
    i2 = g - k;
    if (i2 < 0) i2 = 0;
    c = LnFactorial(g)-LnFactorial(i0)-LnFactorial(i1)-LnFactorial(i2);
    b += 2.0*c;
    b = exp(b);
    x += a*b*r[i]*r[i];
  }

  return x;
}

int PrepCoulombBethe(int ne2, int nte, int ne1, double z,
		     double *e2, double *te, double *e1,
		     int nkl, double *kl, int mode) {
  int i, ie1, ie2, ite, i1, k, m, ik, ierr, q0, q1;
  int t1, ik2, mq;
  double ee1, ee0, k0, k1, k0s, k1s, a, b, r0, r1;
  double w2, w3, *wr, *ws, *wt, *w, *tcb, *r[CBMULT];
  
  if (ne2 > MAXNE || ne1 > MAXNE || nte > MAXNTE) {
    printf("Array multipoles not large enough in CoulombMultipoles\n");
    exit(1);
  }
  
  for (i = 0; i < CBMULT; i++) {
    r[i] = malloc(sizeof(double)*(i+2)*(CBLMAX+1));
  }
  for (ie2 = 0; ie2 < MAXNE; ie2++) {
    for (ite = 0; ite < MAXNTE; ite++) {
      for (ie1 = 0; ie1 < MAXNE; ie1++) {
	for (i = 0; i < MAXNCB; i++) {
	  free(_cb[ie2][ite][ie1][i]);
	  _cb[ie2][ite][ie1][i] = malloc(sizeof(double)*nkl);
	}
      }
    }
  }
  
  w = malloc(sizeof(double)*(CBLMAX+1));
  wt = malloc(sizeof(double)*(CBLMAX+1));
  ws = malloc(sizeof(double)*(CBLMAX+1));
  wr = malloc(sizeof(double)*nkl);
  for (ie1 = 0; ie1 < ne1; ie1++) {
    ee1 = e1[ie1];
    k1s = 2.0*ee1;
    k1 = sqrt(k1s);
    for (ie2 = 0; ie2 < ne2; ie2++) {
      for (ite = 0; ite < nte; ite++) {
	ee0 = ee1 + e2[ie2] + te[ite];
	k0s = 2.0*ee0;
	k0 = sqrt(k0s);
	a = ee1/ee0;
	b = -log(10.0)/log(a);
	q0 = 3.0*b;
	q1 = 7.0*b;
	if (q0 > CBLMIN) q0 = CBLMIN;
	if (q1 > CBLMAX) q1 = CBLMAX;
	if (a < 0.6 || q1 <= q0) {
	  a = a/(1-a);
	  for (i = 0; i < MAXNCB; i++) {
	    for (m = 0; m < nkl; m++) {
	      _cb[ie2][ite][ie1][i][m] = a;
	    }
	  }
	  continue;
	}
	if (q0 > 1) {
	  r0 = 0.5*(sqrt(z*z+2.0*q0*(q0+1.0)*ee0)-z)/(2.0*ee0);
	} else {
	  r0 = 0.1/z;
	}
	r1 = -1.0;
	for (i = 0; i < CBMULT; i++) {
	  i1 = i + 1;
	  if (mode == 0) {
	    CMULTIP(r0, r1, z, k0, k1, i1, q0, q1, 
		    r[i], 1, &ierr);
	  } else {
	    CMULTIP(r0, r1, z, k1, k0, i1, q0, q1, 
		    r[i], 1, &ierr);
	  }
	  for (k = q0; k <= q1; k++) {
	    ik = k - q0;
	    ik2 = ik*(i+2);
	    wt[k] = PartialOmega(k, r[i]+ik2, i1);
	    if (k > q0) {
	      w[k-1] = wt[k]/wt[k-1];
	      if (wt[k]/wt[q0] < EPS6) {
		k++;
		break;
	      }
	    }
	    if (i1 == 1) {
	      if (r[i][ik2] < 0 || r[i][ik2+1] < 0) 
		break;
	    }
	  }
	  if (i1 == 1) {
	    q1 = k-1;
	    if (q1 <= q0) {
	      a = ee1/(ee0 - ee1);
	      for (i = 0; i < MAXNCB; i++) {
		for (m = 0; m < nkl; m++) {
		  _cb[ie2][ite][ie1][i][m] = a;
		}
	      }
	      goto NEXTE;
	    }
	    for (k = q0; k < q1; k++) {
	      ik = k - q0;
	      ik2 = ik*2;
	      w3 = r[i][ik2+2]/r[i][ik2];
	      w2 = r[i][ik2+1]/r[i][ik2];
	      w3 *= w3;
	      w2 *= w2;
	      if (mode == 0) {
		a = z*z + (k+1.0)*(k+1.0)*k0s;
		a *= w3-w2;
	      } else {
		a = z*z + (k+1.0)*(k+1.0)*k1s;
		a *= w2 - w3;
	      }
	      b = k + (k+1.0)*w2;
	      b *= (k+1.0)*2.0*(te[ite]+e2[ie2]);
	      ws[k] = a/b;
	    }
	    if (mode == 0) {
	      tcb = GetCoulombBethe(ie2, ite, ie1, i1, 0);
	      for (m = 0; m < nkl; m++) {
		k = (int) kl[m];
		if (k >= q1) {
		  tcb[m] = ws[q1-1];
		} else if (k < q0) {
		  tcb[m] = ws[q0];
		} else {
		  tcb[m] = ws[k];
		}
	      }
	    }
	    k = q1;
	  }
	  for (; k <= CBLMAX; k++) {
	    w[k-1] = w[k-2];
	  }
	  w[k-1] = w[k-2];
	  for (k = 0; k < q0; k++) {
	    w[k] = w[q0];
	  }
	  if (mode == 0) {
	    if (i1 > 1) {
	      tcb = GetCoulombBethe(ie2, ite, ie1, i1, 0);
	      CoulombBetheTail(CBLMAX+1, w, nkl, kl, tcb);
	    }
	  } else if (i1 == 1) {
	    CoulombBetheTail(CBLMAX+1, w, nkl, kl, wr);
	    for (m = 0; m < nkl; m++) {
	      k = (int) kl[m];
	      if (k < q0) {
		wr[m] = ws[q0]/wr[m];
	      } else if (k > q1) {
		wr[m] = ws[q1-1]/wr[m];
	      } else {
		wr[m] = ws[k]/wr[m];
	      }
	    }
	  }
	}
	if (mode) {
	  for (t1 = 0; t1 < CBMULT; t1++) {
	    for (mq = 0; mq <= t1+1; mq++) {
	      for (k = q0; k <= q1; k++) {		  
		ik = k - q0;
		ik2 = ik*(t1+2);
		wt[k] = PartialOmegaMSub(k, r[t1]+ik2, t1+1, mq, z, ee0);
		if (k > q0) w[k-1] = wt[k]/wt[k-1];
	      }		
	      for (; k <= CBLMAX; k++) {
		w[k-1] = w[k-2];
	      }
	      w[k-1] = w[k-2];
	      for (k = 0; k < q0; k++) {
		w[k] = w[q0];
	      }
	      tcb = GetCoulombBethe(ie2, ite, ie1, t1+1, mq);
	      CoulombBetheTail(CBLMAX+1, w, nkl, kl, tcb);
	      if (t1 == 0) {
		for (m = 0; m < nkl; m++) {
		  tcb[m] *= wr[m];
		}
	      }
	    }
	  }
	}
      NEXTE:
	continue;
      }
    }
  }
	
  for (i = 0; i < CBMULT; i++) {
    free(r[i]);
  }
  free(w);
  free(wt);
  free(ws);
  free(wr);

  return 0;
}

int CoulombBethe(char *s, double z, double te, double e1) {
#define M 200
  double kl[M];
  int i, j, k, m;
  FILE *f;
  int nte, ne1, ne2, ie1, ite;
  double e2, *tcb;
  
  ne2 = 1;
  ne1 = 1;
  nte = 1;
  e2 = 0.0;
  f = fopen(s, "w");
  fprintf(f, "## TE = %10.3E, E1 = %10.3E\n", te, e1);

  te /= HARTREE_EV;
  e1 /= HARTREE_EV;
  for (i = 0; i < M; i++) {
    kl[i] = i;
  }

  PrepCoulombBethe(ne2, nte, ne1, z, &e2, &te, &e1, M, kl, 1);
  for (i = 1; i < M; i++) {
    fprintf(f, "%3d %3d ", i, (int)kl[i]);
    for (j = 0; j < CBMULT; j++) {
      for (m = 0; m <= j+1; m++) {
	tcb = GetCoulombBethe(0, 0, 0, j+1, m);
	fprintf(f, "%12.5E ", tcb[i]);
      }
    }
    fprintf(f, "\n");
  }

  fclose(f);

  return 0;
#undef M
}

int CoulombMultip(char *fn, double z, double te, double e1,
		   int k, int q0, int q1, int m) {
  FILE *f;
  double *r, r0, r1, k0, k1;
  int i, j, ierr;

  k1 = sqrt(2.0*e1/HARTREE_EV);
  k0 = sqrt(2.0*(e1+te)/HARTREE_EV);

  r0 = 0.1/z;
  r1 = -1.0;
  r = (double *) malloc(sizeof(double)*(k+1)*(q1-q0+1));
  CMULTIP(r0, r1, z, k0, k1, k, q0, q1, r, m, &ierr);
  if (ierr != 0) {
    printf("error in CMULTIP: %d\n", ierr);
    free(r);
    return -1;
  }
  f = fopen(fn, "w");
  fprintf(f, "# %10.3E %10.3E %10.3E %15.8E %15.8E %2d\n",
	  z, te, e1, k0, k1, k);
  j = 0;
  for (m = q0; m <= q1; m++) {
    fprintf(f, "%4d", m);
    for (i = -k; i <= k; i += 2) {
      fprintf(f, " %17.10E", r[j]);
      j++;
    }
    fprintf(f, "\n");
  }
  fclose(f);
  free(r);
  
  return 0;
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
  PrepCBIndex();
  SetHydrogenicNL(-1, -1, -1, -1);

  dipole_array = (ARRAY *) malloc(sizeof(ARRAY));
  ArrayInit(dipole_array, sizeof(double *), DIPOLE_BLOCK);

  return 0;
}
