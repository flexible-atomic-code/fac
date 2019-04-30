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

#include "orbital.h"
#include "coulomb.h"
#include "cf77.h"
#include "mpiutil.h"
#include "fftsg.h"

static char *rcsid="$Id$";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

/* the following arrays provide storage space in the calculation */
static double *_veff, *ABAND, *_dwork, *_dwork1, *_dwork2, *_dwork3;
static int *_ipiv;

#pragma omp threadprivate(_veff, ABAND, _ipiv, _dwork, _dwork1, _dwork2, _dwork3)

static int _relativistic_fermi = 1;
static double _fermi_abserr = 1e-7;
static double _fermi_relerr = 1e-5;
static double _fermi_stablizer = 0.5;
#define NFRM1 128
#define NFRM2 (NFRM1+2)
static double _fermi_rm1[4][NFRM2];
static char _fermi_rmf[256] = "";
static char _sp_ofn[256] = "";
static double _sp_vmax = 1e30;
static double _sp_rmax = 1.0;
static double _sp_yeps = 1e-7;
static double _sp_neps = 1e-3;
static double _sp_z0 = 1.0;
static double _sp_z1 = 1.0;
static int _sp_nzs = 0;
static double *_sp_zs = NULL;
static double *_sp_zw = NULL;
static int _sp_mode = 2;
static int _sp_print = 0;
static double _icf_tol = EPS3;
static int _icf_maxiter = 1024;
static int _icf_nfft = 0;
static double _icf_rmax = 50.0;
static double _icf_kb = 2.0;
static double _icf_kd = 0.5;
static double _icf_stablizer = 0.5;
static int _icf_ozc = 0;
static int _icf_spmi = 0;
static double _icf_tmax = 5.0;
static double _icf_tidx = 3.0;
static double *_icf_dw = NULL;
static char _icf_ofn[256] = "";
static int max_iteration = 512;
static double wave_zero = 1E-10;

static double _wk_q[16] = {2.09400223235,
			   1.42425723499,
			   0.76419630941,
			   0.36124746165,
			   0.15118743043,
			   0.05247168246,
			   0.01797758276,
			   -0.04761316284,
			   0.19713924639,
			   -0.62421312791,
			   1.29259217859,
			   -1.89507086023,
			   1.89089237003,
			   -1.24611625668,
			   0.48793792307,
			   -0.08858438257};
static double _wk_p[36] = {2.094001968517539e-2,
			   -6.328682056509706e-2,
			   1.656675047469226e-1,
			   6.218254052794603e-2,
			   8.561450880149488e-1,
			   -1.324099530525646,
			   1.022894990055742e-1,
			   1.818708880809046e-1,
			   5.002561255243687e-4,
			   2.167778174657628e-1,
			   5.453835857080859e-1,
			   3.137663571113079e-1,
			   -1.619774880547481e-1,
			   -6.749903602516955e-2,
			   -2.935144138820853e-4,
			   9.755800734714044e-2,
			   -1.994016822313688e-1,
			   -4.777623414220549e-2,
			   7.503472520308676e-3,
			   3.440742164594281e-5,
			   2.8294209897e-3,
			   1.4195518237e-4,
			   1.3489833978e-5,
			   2.0134008556e-7,
			   8.6029184536e-6,
			   -1.7870165724e-5,
			   2.3269607856e-5,
			   -1.7798463760e-5,
			   8.7274744996e-6,
			   -2.9235606710e-6,
			   6.8929360802e-7,
			   -1.1491699126e-7,
			   1.3302835807e-8,
			   -1.0190988753e-9,
			   4.6508041566e-11,
			   -9.5785924366e-13};
static int _wk_s[4][4] = {{1, 1, 0, 0},
			  {0, -3, 0, 10},
			  {-1, 2, 0, 10},
			  {-13, -6, 10, 10}};
static double _wk_a[4][4][8] =
  {{{2.9675,-4.4288,5.9339e1,-5.682e2,4.6131e3,-2.6282e4,8.9448e4,-1.3412e5},
    {7.8894e-1,-2.9946e1,8.4881e2,-1.4888e4,1.6003e5,-1.0259e6,3.5980e6,-5.3067e6},
    {-3.0201e-4,-9.0441e-1,5.8388e1,-1.5899e3,2.880e4,-1.8036e5,7.3716e5,-1.2228e6},
    {4.9725e-4,2.4826,-1.0745e2,2.7128e3,-3.8249e4,3.01e5,-1.2368e6,2.0678e6}},
   {{-6.64673e-1,5.70540e-1,-9.43955e-1,1.18405,-8.91689e-1,3.64707e-1,-6.17902e-2,0.0},
    {4.68647e-5,-9.83836e-4,8.28746e-3,5.94372e-2,-5.31431e-2,7.52653e-3,-4.54081e-4,0.0},
    {2.56839e-2,-3.21021e-3,1.01163e-2,-1.62093e-2,1.17762e-2,-4.69814e-3,7.76469e-4,0.0},
    {0,0,0,0,0,0,0,0}},
   {{1.93743e-2,-6.43399e-1,2.56565e-1,-1.09235e-1,4.14226e-2,-7.28744e-3,4.54914e-4,0.0},
    {1.64381e-1,-3.05149e-1,2.31560e-1,-9.49948e-2,2.22450e-2,-2.80168e-3,1.47382e-4,0.0},
    {2.53250e-2,3.56651e-3,-8.14980e-3,6.21168e-3,-3.29527e-3,7.60885e-4,-6.03993e-5,0.0},
    {0,0,0,0,0,0,0,0}},
   {{-2.4734758e9,3.2270784e9,-1.7207915e9,4.7635740e8,-7.1414759e7,5.4305859e6,-1.6479928e5,0.0},
    {-1.5115262e4,1.4945586e4,-5.8361389e3,1.1476935e3,-1.2098290e2,6.5439613,-1.4296421e-1,0.0},
    {0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0}}};
static double _wk_fd[5] = {225.0, 1323.0, 4725.0, 190575.0, 12297285.0};
static double _wk_f[16][5] =
  {{2.0, 0, 0, 0, 0},
   {59.0, 0, 0, 0, 0},
   {1977.0, 20, 0, 0, 0},
   {1258680.0, 34144, 0, 0, 0},
   {1960420032.0, 93618070.0, 96096.0, 0, 0},
   {0.3812769,0.0266795,0.000086,0,0},
   {0.4946161,0.0458093,0.0002918,6e-8,0},
   {0.6150374,0.0707397,0.0007216,5e-7,0},
   {0.7412290,0.1014082,0.0014705,23e-7,1e-10},
   {0.8721898,0.1376368,0.0026288,69e-7,2e-9},
   {1.0071432,0.1791838,0.0042776,167e-7,8e-9},
   {1.1454776,0.2257770,0.0064866,346e-7,3e-8},
   {1.2867042,0.2771342,0.0093134,643e-7,9e-8},
   {1.4304276,0.3329753,0.0128043,1097e-7,2e-7},
   {1.5763244,0.3930296,0.0169953,1750e-7,5e-7},
   {1.7241274,0.4570401,0.0219132,2644e-7,9e-7}};
    
static int SetVEffective(int kl, int kv, POTENTIAL *pot);
static int TurningPoints(int n, double e, POTENTIAL *pot);
static int CountNodes(double *p, int i1, int i2);
static int IntegrateRadial(double *p, double e, POTENTIAL *pot, 
			   int i1, double p1, int i2, double p2, int q);
static double Amplitude(double *p, double e, int kl, POTENTIAL *pot, int i1);
static int Phase(double *p, POTENTIAL *pot, int i1, double p0);

int SPMode() {
  return _sp_mode;
}

void SetOrbitalWorkSpace(double *p, int n) {
  _veff = p;
  p += n;
  ABAND = p;
  p += n*4;
  _dwork = p;
  p += n;
  _dwork1 = p;
  p += n;
  _dwork2 = p;
  p += n;
  _dwork3 = p;
  p += n;
  _ipiv = (int *) p;
}
 
double EneTol(double e) {
  e = fabs(e);
  double d0 = e*ENERELERR;
  if (d0 > ENEABSERR) {
    d0 = e*ENERELERR1;
    d0 = Max(ENEABSERR, d0);
  }
  return d0;
}

double EnergyDefect(double z, double dn, int j, int kl, double dkl) {
  double dki, e0, e1;
  int ki0, ki1, ji0, ji1, ka0, ka1, x;
  
  dki = dkl*dn;
  ki0 = 2*((int)dki);
  ki1 = ki0+2;
  if (j > kl) {
    ji0 = ki0 + 1;
    ji1 = ki1 + 1;
  } else {
    if (ki0 == 0) ji0 = 1;
    else ji0 = ki0 - 1;
    if (ki1 == 0) ji1 = 1;
    else ji1 = ki1 - 1;
  }
  ka0 = GetKappaFromJL(ji0, ki0);
  ka1 = GetKappaFromJL(ji1, ki1);
  e0 = EnergyH(z, dn, ka0);
  e1 = EnergyH(z, dn, ka1);
  x = (2*dki - ki0)/2.0;
  return e0*(1-x) + e1*x;
}

double QuantumDefect(double z, int n, int ka, double e) {
  int j, kl, i;
  double dkl, dn, eh, dn0, dn1, x;
  
  GetJLFromKappa(ka, &j, &kl);
  kl /= 2;
  dn = (double) n;
  dkl = kl/dn;
  eh = EnergyH(z, dn, ka);
  if (eh == e) return 0.0;
  dn0 = dn;
  dn1 = dn;
  if (eh < e) {
    i = 0;
    while (1) {
      i++;
      if (i > 1000) return 0.0;
      dn1 *= 1.1;
      eh = EnergyDefect(z, dn1, j, kl, dkl);
      if (eh > e) break;
    }
  } else {
    i = 0;
    while (1) {
      i++;
      if (i > 1000) return 0.0;
      dn0 /= 1.1;
      eh = EnergyDefect(z, dn0, j, kl, dkl);
      if (eh < e) break;
    }
  }
  i = 0;
  while (fabs(dn1-dn0)/fabs(dn1+dn0) > 1e-8) {
    i++;
    if (i > 1000) return 0.0;
    x = 0.5*(dn0+dn1);
    eh = EnergyDefect(z, x, j, kl, dkl);
    if (eh < e) {
      dn0 = x;
    } else if (eh > e) {
      dn1 = x;
    } else {
      break;
    }
  }
  x = 0.5*(dn0+dn1);
  return n-x;
}
  
double EnergyH(double z, double n, int ka) {
  double a, np;

  ka = abs(ka);
  a = FINE_STRUCTURE_CONST*z;
  
  np = n + sqrt(ka*ka - a*a) - ka;
  a = a/np;
  a = (1.0/sqrt(1.0 + a*a)) - 1.0;
  a /= FINE_STRUCTURE_CONST2;
  
  return a;
}
    
double RadialDiracCoulomb(int npts, double *p, double *q, double *r,
			  double z, int n, int kappa) {
  int i, iordr1, iordr2, k, nr;
  double a, alfa, an1, an2, argr, b, bn, bign, bignmk, nrfac;
  double eps, fac, facn, fden, ff, fg, fk, f1, f2, gamma;
  double ovlfac, rgamm1, rgamm2, rho, rhon, twogp1, zalfa;
  double *ta, *tb;
  double energy, t;
  
  ta = _dwork1;
  tb = _dwork2;
  
  alfa = FINE_STRUCTURE_CONST;
  k = abs(kappa);
  nr = n - k;
  fk = (double) k;
  zalfa = z*alfa;
  gamma = sqrt(fk*fk - zalfa*zalfa);
  twogp1 = gamma*2.0 + 1.0;
  bign = sqrt(n*n - 2.0*(n - fk)*(fk - gamma));
  t = zalfa/(gamma + n - fk);
  eps = 1.0/sqrt(1.0 + t*t);

  energy = -(1.0 - eps)/FINE_STRUCTURE_CONST2;
  
  nrfac = 0.0;
  for (i = 1; i <= nr; i++) nrfac += log(i);
  
  argr = twogp1 + nr;
  rgamm1 = DLOGAM(argr);
  argr = twogp1;
  rgamm2 = DLOGAM(argr);

  /*
  fac = - sqrt(rgamm1) / (rgamm2*sqrt((double)nrfac)) *
    sqrt(z/(2.0*bign*bign*(bign-kappa)));
  */
  fac = -exp(0.5*rgamm1 - rgamm2 - 0.5*nrfac)*sqrt(z/(2.0*bign*bign*(bign-kappa)));

  if (kappa > 0) fac = -fac;

  fg = fac*sqrt(1.0 + eps);
  ff = fac*sqrt(1.0 - eps);
  
  if (nr == 0) {
    iordr1 = 0;
    iordr2 = 0;
  } else {
    iordr1 = nr - 1;
    iordr2 = nr;
  }

  fac = 1.0;
  facn = 1.0;
  a = -nr;
  an1 = a + 1.0;
  an2 = a;
  b = twogp1;
  bn = b;
  
  k = 0;
  
  while (1) {
    k = k + 1;
    fden = 1.0/(facn*bn);
    if (k <= iordr1) {
      ta[k] = an1 * fden;
    }

    if (k > iordr2) break;
    tb[k] = an2*fden;
    a += 1.0;
    an1 *= a + 1.0;
    an2 *= a;
    b += 1.0;
    bn *= b;
    fac += 1.0;
    facn *= fac;
  }

  p[0] = 0.0;
  q[0] = 0.0;
  fac = 2.0*z/bign;
  bignmk = bign - kappa;
  for (i = 0; i < npts; i++) {
    rho = fac * r[i];
    rhon = rho;
    k = 0;
    f1 = 1.0;
    f2 = 1.0;
    while (1) {
      k = k+1;
      if (k <= iordr1) {
	f1 += ta[k]*rhon;
      }
      if (k > iordr2) break;
      f2 += tb[k]*rhon;
      rhon *= rho;
    }
    f1 *= nr;
    f2 *= bignmk;
    ovlfac = exp(-0.5*rho)*pow(rho, gamma);
    p[i] = fg*ovlfac * (f1-f2);
    q[i] = ff*ovlfac * (f1+f2);
  }

  return energy;
}

void InitOrbitalData(void *p, int n) {
  ORBITAL *d;
  int i;
  
  d = (ORBITAL *) p;
  for (i = 0; i < n; i++) {
    d[i].wfun = NULL;
    d[i].phase = NULL;
    d[i].ilast = -1;
    d[i].pdx = 0;
    d[i].bqp0 = 0;
    d[i].bqp1 = 0;
    d[i].se = 1e31;
    d[i].ose = 1e31;
    d[i].qed = 0.0;
    d[i].kv = -1000000000;
    d[i].horb = NULL;
    d[i].rorb = NULL;
    d[i].isol = 0;
  }
}

int RadialSolver(ORBITAL *orb, POTENTIAL *pot) {
  int ierr;
  int nm, km, k;
  double z;
  
  orb->rfn = 0;
  if (orb->n > 0) {
    if (orb->n == 1000000) {
      ierr = RadialFreeInner(orb, pot);
    } else {
      if (pot->ib > 0 && orb->n > pot->nb) {
	ierr = RadialBasis(orb, pot);
      } else {
	GetHydrogenicNL(NULL, NULL, &nm, &km);
	k = GetLFromKappa(orb->kappa);
	k /= 2;
	z = GetResidualZ();
	if (orb->n > nm || k > km) {
	  if (z < 1) z = 1;
	  orb->energy = EnergyH(z, (double)(orb->n), orb->kappa);
	  orb->ilast = -1;
	  orb->wfun = NULL;
	  orb->isol = 1;
	  orb->kv = 0;
	  orb->dn = 0.0;
	  if (pot->pse) orb->kv = IdxVT(orb->kappa);
	  return 0;
	}
	if (orb->n < pot->nmax) {
	  ierr = RadialBound(orb, pot);
	} else {
	  ierr = RadialRydberg(orb, pot);
	}
	orb->dn = QuantumDefect(z, orb->n, orb->kappa, orb->energy);
      }
      for (k = 0; k < pot->maxrp; k++) {
	if (pot->rad[k] > pot->atom->rms0) break;
      }
      for (; k < pot->maxrp; k++) {
	if (orb->wfun[k] > 0 && orb->wfun[k+1] <= 0) break;
	if (orb->wfun[k] < 0 && orb->wfun[k+1] >= 0) break;
      }
      if (k > 0) {
	orb->rfn = pot->rad[k-1];
      }
    }
  } else if (orb->n == 0) {
    ierr = RadialFree(orb, pot);
  } else {
    ierr = RadialBasisOuter(orb, pot);
  }
  if (orb->wfun != NULL && orb->ilast > 0) orb->isol = 1;
  if (orb->wfun) {
    if (isnan(orb->wfun[0]) || isnan(orb->wfun[pot->maxrp])) {
      MPrintf(-1, "wfun nan: %d %d %12.5E %12.5E %12.5E\n",
	      orb->n, orb->kappa, orb->energy,
	      orb->wfun[0], orb->wfun[pot->maxrp]);
    }
  }
  return ierr;
}

double *GetVEffective(void) { 
  return _veff;
}

int LastMaximum(double *p, int i1, int i2) {
  int i;

  if (p[i2] > p[i2-1]) {
    for (i = i2-1; i > i1; i--) {
      if (p[i-1] > p[i]) break;
    }
  } else if (p[i2] < p[i2-1]) {
    for (i = i2-1; i > i1; i--) {
      if (p[i-1] < p[i]) break;
    }
  } else {
    i = i2;
  }
  if (i == i1) {    
    i = i2;
  }
  return i;
} 

void DrLargeSmall(ORBITAL *orb, POTENTIAL *pot, double *pr, double *qr) {
  int i, kv;
  double a = FINE_STRUCTURE_CONST;
  double *p = orb->wfun;
  double *q = orb->wfun+pot->maxrp;
  double b, c;
  
  kv = orb->kv;
  if (kv < 0 || kv > NKSEP) kv = 0;

  for (i = 0; i <= orb->ilast; i++) {
    b = a*(orb->energy - pot->VT[kv][i]);
    c = orb->kappa/pot->rad[i];
    pr[i] = (b+2.0/a)*q[i];
    pr[i] -= p[i]*c;
    qr[i] = -b*p[i];
    qr[i] += q[i]*c;
  }
  for (i = orb->ilast+1; i < pot->maxrp; i++) {
    pr[i] = 0;
    qr[i] = 0;
  }
}

void Differential(double *p, double *dp, int i1, int i2, double *drdrho) {
  double c = 1.0/24.0;
  double b;
  int i;

  b = -50.0*p[i1];
  b += 96.0*p[i1+1];
  b -= 72.0*p[i1+2];
  b += 32.0*p[i1+3];
  b -= 6.0 *p[i1+4];
  dp[i1] = b*c;

  b = -6.0*p[i1];
  b -= 20.0*p[i1+1];
  b += 36.0*p[i1+2];
  b -= 12.0*p[i1+3];
  b += 2.0 *p[i1+4];
  dp[i1+1] = b*c;

  b = -50.0*p[i2];
  b += 96.0*p[i2-1];
  b -= 72.0*p[i2-2];
  b += 32.0*p[i2-3];
  b -= 6.0 *p[i2-4];
  dp[i2] = -b*c; 

  b = -6.0*p[i2];
  b -= 20.0*p[i2-1];
  b += 36.0*p[i2-2];
  b -= 12.0*p[i2-3];
  b += 2.0 *p[i2-4];
  dp[i2-1] = -b*c;

  for (i = i1 + 2; i < i2-1; i++) {
    b = 2.0*(p[i-2] - p[i+2]);
    b += 16.0*(p[i+1] - p[i-1]);
    dp[i] = b*c;
  }
  for (i = i1; i <= i2; i++) {
    dp[i] /= drdrho[i];
  }
}    

double DpDr(int kappa, int k, int i, double e, POTENTIAL *pot,
	    double b, int m, double *bqp) {
  double x2, dx, dr, pr, z0, kl;
  
  x2 = 1 + 0.5*FINE_STRUCTURE_CONST2*(e - pot->VT[k][i]);
  if (m == 0) {
    b = (b + kappa)*FINE_STRUCTURE_CONST;
    b /= (2.0*pot->rad[i]);
    if (bqp) *bqp = b;
    b = 2*x2*b/FINE_STRUCTURE_CONST - kappa/pot->rad[i];
  } else if (m == 1) {
    if (bqp) *bqp = b;
    b = 2*x2*b/FINE_STRUCTURE_CONST - kappa/pot->rad[i];
  } else if (m == -1) {
    if (pot->atom->z1 > 0) {
      z0 = pot->atom->z1;
      kl = GetLFromKappa(kappa);
      if (kappa < 0) {
	b = (e + z0)*FINE_STRUCTURE_CONST;
	b *= -pot->rad[i]/(kl+3.0);
      } else {
	b =(e + z0 + 2/FINE_STRUCTURE_CONST2)*FINE_STRUCTURE_CONST;
	b *= pot->rad[i]/(kl+1.0);
	b = 1.0/b;
      }   
      if (bqp) {
	*bqp = b;
      }
      b = (1 + kl/2)/pot->rad[i];
    } else {
      z0 = pot->Z[0];
      double g = sqrt(kappa*kappa - z0*z0*FINE_STRUCTURE_CONST2);
      b = z0*FINE_STRUCTURE_CONST/(kappa - g);   
      if (bqp) {
	*bqp = b;
      }
      b = g/pot->rad[i];
    } 
  } else if (m == -2) {
    b = -sqrt(1.0/FINE_STRUCTURE_CONST2 - e*e*FINE_STRUCTURE_CONST2);
    if (bqp) {
      *bqp = (b + kappa/pot->rad[i])*FINE_STRUCTURE_CONST/(2*x2);
    }
  }
  b *= pot->dr_drho[i];
  dx = -0.25*FINE_STRUCTURE_CONST2*pot->dVT[k][i]/x2;
  dx *= pot->dr_drho[i];
  dr = 0.5*pot->dr_drho[i]/pot->rad[i];
  dr *= (1-pot->qr*pot->qr*pot->dr_drho[i]*pot->ar*pow(pot->rad[i], pot->qr-1));
  pr = b - dx - dr;
  return pr;
}
  
int RadialBasisOuter(ORBITAL *orb, POTENTIAL *pot) {  
  double e, de, delta, emin, emax, dr, ke;
  double *p, p1, p2, qo, qi, bqp, bqp1;
  int n, i, k, kl, nr, nodes, niter;
  int i2, i2m1, i2m2, i2p1, i2p2, i1;
  int ib0, ib1;

  int kv = 0;
  if (pot->pse) kv = IdxVT(orb->kappa);
  orb->kv = kv;
  kl = orb->kappa;
  kl = (kl < 0)? (-kl-1):kl;

  n = -orb->n;
  nr = n - kl - 1;
  if (kl < 0 || kl >= n) {
    printf("Invalid orbital angular momentum, L=%d, %d %d\n", 
	   kl, orb->n, orb->kappa);
    return -1;
  }
  
  p = malloc(sizeof(double)*2*pot->maxrp);
  if (!p) return -1;
  
  ib0 = pot->ib-1;
  ib1 = pot->ib1+1;
  niter = 0;
  
  SetPotentialW(pot, 0.0, orb->kappa, kv);
  SetVEffective(kl, kv, pot);
  emin = 1E30;
  for (i = pot->ib; i <= pot->ib1; i++) {
    if (_veff[i] < emin) emin = _veff[i];
  }
  dr = pot->rad[pot->ib1]-pot->rad[pot->ib];
  ke = TWO_PI*(2*nr+5)/dr;
  emax = 0.5*ke*ke;
  ke = TWO_PI*(nr+1)/dr;
  e = 0.5*ke*ke;
  de = 0.5*e;
  while (niter < max_iteration) {
    niter++;
    SetPotentialW(pot, e, orb->kappa, kv);
    SetVEffective(kl, kv, pot);
    bqp = DpDr(orb->kappa, kv, pot->ib, e, pot, pot->bqp, 0, NULL);
    bqp1 = DpDr(orb->kappa, kv, pot->ib1, e, pot, pot->bqp, 0, NULL);
    i2 = TurningPoints(orb->n, e, pot);
    nodes = IntegrateRadial(p, e, pot, ib0, bqp, i2, 1.0, 2);
    if (nodes > 1) {
      i2 = LastMaximum(p, pot->ib+20, i2);
      nodes = CountNodes(p, ib0, i2);
    }
    nodes += IntegrateRadial(p, e, pot, i2, 1.0, ib1, bqp1, 1);
    if (nodes == nr) break;
    else if (nodes > nr) {
      emax = e;
      e = 0.5*(emin + e);
    } else {
      emin = e;
      e = 0.5*(emax + e);
    }
  }
  if (niter == max_iteration) {
    printf("Max iteration before finding correct nodes in RadialBasisOuter %d %d\n",
	   nodes, nr);
    free(p);
    return -2;
  }
  
  niter = 0;
  de = 0.25*ke*ke;
  while (de > ENEABSERR) {
    de *= 0.25;
    if (niter == max_iteration) break;
    while (nodes == nr) {
      niter++;
      if (niter == max_iteration) break;
      e += de;
      SetPotentialW(pot, e, orb->kappa, kv);
      SetVEffective(kl, kv, pot);
      bqp = DpDr(orb->kappa, kv, pot->ib, e, pot, pot->bqp, 0, NULL);
      bqp1 = DpDr(orb->kappa, kv, pot->ib1, e, pot, pot->bqp, 0, NULL);
      i2 = TurningPoints(orb->n, e, pot);
      nodes = IntegrateRadial(p, e, pot, ib0, bqp, i2, 1.0, 2);
      if (nodes > 1) {
	i2 = LastMaximum(p, pot->ib+20, i2);
	nodes = CountNodes(p, ib0, i2);
      }
      nodes += IntegrateRadial(p, e, pot, i2, 1.0, ib1, bqp1, 1);
    }
    e -= de;
    nodes = nr;
  }
  if (niter == max_iteration) {
    printf("Max iteration before finding solution in RadialBasisOuter %d %d\n",
	   nodes, nr);
    free(p);
    return -3;
  }
  niter = 0;
  emax = e;
  if (nr > 0) {
    de = 0.25*ke*ke;
    while (de > ENEABSERR) {
      de *= 0.25;
      if (niter == max_iteration) break;
      while (nodes == nr) {
	niter++;
	if (niter == max_iteration) break;
	e -= de;
	SetPotentialW(pot, e, orb->kappa, kv);
	SetVEffective(kl, kv, pot);
	bqp = DpDr(orb->kappa, kv, pot->ib, e, pot, pot->bqp, 0, NULL);
	bqp1 = DpDr(orb->kappa, kv, pot->ib1, e, pot, pot->bqp, 0, NULL);
	i2 = TurningPoints(orb->n, e, pot);
	nodes = IntegrateRadial(p, e, pot, ib0, bqp, i2, 1.0, 2);
	if (nodes > 1) {
	  i2 = LastMaximum(p, pot->ib+20, i2);
	  nodes = CountNodes(p, ib0, i2);
	}
	nodes += IntegrateRadial(p, e, pot, i2, 1.0, ib1, bqp1, 1);
      }
      e += de;
      nodes = nr;
    }  
    if (niter == max_iteration) {
      printf("Max iteration before finding solution in RadialBasisOuter %d %d\n",
	     nodes, nr);
      free(p);
      return -4;
    }
    emin = e;
  }
  e = 0.5*(emin+emax);
  niter = 0;
  while (niter < max_iteration) {
    niter++;
    SetPotentialW(pot, e, orb->kappa, kv);
    SetVEffective(kl, kv, pot); 
    bqp = DpDr(orb->kappa, kv, pot->ib, e, pot, pot->bqp, 0, NULL);
    bqp1 = DpDr(orb->kappa, kv, pot->ib1, e, pot, pot->bqp, 0, NULL);
    i2 = TurningPoints(orb->n, e, pot);
    i2p2 = i2 + 2;
    nodes = IntegrateRadial(p, e, pot, ib0, bqp, i2p2, 1.0, 2);
    for (i = ib0; i <= i2p2; i++) {
      p[i] *= pot->dr_drho2[i];
    }
    if (nodes > 1) {
      i2 = LastMaximum(p, pot->ib+20, i2);
    }
    i2m1 = i2 - 1;
    i2m2 = i2 - 2;
    i2p1 = i2 + 1;    
    i2p2 = i2 + 2;
    p1 = p[i2m2]/pot->dr_drho2[i2m2];
    qo = (-4.0*p[i2m2-1] + 30.0*p[i2m2] - 120.0*p[i2m1]
	  + 40.0*p[i2] + 60.0*p[i2p1] - 6.0*p[i2p2])/120.0;
    qo /= p[i2];
    k = IntegrateRadial(p, e, pot, i2m2, p1, ib1, bqp1, 1);
    for (i = i2m2; i <= ib1; i++) {
      p[i] *= pot->dr_drho2[i];
    }
    qi = (6.0*p[i2m2] - 60.0*p[i2m1] - 40.0*p[i2] + 120.0*p[i2p1]
	  - 30.0*p[i2p2] + 4.0*p[i2p2+1])/120.0;
    qi /= p[i2];
    delta = qo - qi;
    /*
    printf("%3d %2d %2d %3d %15.8E %15.8E %15.8E %12.5E %12.5E %12.5E\n",
	   niter, nodes, k, i2, e, emin, emax, qo, qi, delta);
    */
    if (delta < -EPS8) {
      emax = e;
      e = 0.5*(emin + e);
    } else if (delta > EPS8) {
      emin = e;
      e = 0.5*(emax + e);
    } else {
      break;
    }
  }
  
  if (niter == max_iteration) {
    printf("Max iteration before finding solution in RadialBasisOuter %d %d\n",
	   nodes, nr);
    free(p);
    return -5;
  }
  
  if (p[pot->ib] < 0) {
    for (i = ib0; i <= ib1; i++) {
      p[i] = -p[i];
    }
  }
  orb->ilast = pot->ib1;
  orb->energy = e;
  orb->wfun = p;
  orb->qr_norm = 1.0;
  if (pot->flag == -1) {
    DiracSmall(orb, pot, ib1, kv);
  }

  return 0;
}

int RadialBasis(ORBITAL *orb, POTENTIAL *pot) {
  double z, z0, e, de, ep, delta, emin, emax;
  double *p, norm2, fact, p1, p2, qo, qi, bqp, bqp0;
  int i, k, kl, nr, nodes, niter;
  int i2, i2m1, i2m2, i2p1, i2p2, i1;

  int kv = 0;
  if (pot->pse) kv = IdxVT(orb->kappa);
  orb->kv = kv;
  kl = orb->kappa;
  kl = (kl < 0)? (-kl-1):kl;

  if (kl < 0 || kl >= orb->n) {
    printf("Invalid orbital angular momentum, L=%d, %d %d\n", 
	   kl, orb->n, orb->kappa);
    return -1;
  }

  p = malloc(sizeof(double)*2*pot->maxrp);
  if (!p) {
    MPrintf(-1, "cannot alloc memory RadialBasis: %d %d\n", orb->n, orb->kappa);
    return -1;
  }

  nr = orb->n - kl - 1;

  niter = 0;
  z0 = GetAtomicNumber();
  z = GetResidualZ();
  //if (z < 0.1) z = 0.1;

  double emin0 = 1.5*EnergyH(z0, orb->n, orb->kappa);
  if (kl > 0) {
    SetPotentialW(pot, 0.0, orb->kappa, kv);
    SetVEffective(kl, kv, pot);
    emin = 1E30;
    for (i = 0; i < pot->maxrp; i++) {
      if (_veff[i] < emin) emin = _veff[i];
    }
    emin += ENEABSERR;
    if (emin < emin0) emin = emin0;
  } else {
    emin = EnergyH(z, orb->n, orb->kappa);
  }

  e = 1.1*fabs(emin);
  while (niter < max_iteration) {
    niter++;
    SetPotentialW(pot, e, orb->kappa, kv);
    SetVEffective(kl, kv, pot);
    i2 = TurningPoints(orb->n, e, pot);
    if (i2 < 0) {
      nodes = -1;
    } else {
      nodes = IntegrateRadial(p, e, pot, 0, 0.0, i2, 1.0, 0);
    }
    if (nodes > nr) break;
    e = e*2.0;
  }
  if (niter == max_iteration) {
    printf("Max iteration before finding correct nodes in RadialBasis %d %d %d %d %g %g\n",
	   nodes, nr, orb->n, kl, e, emin);
    free(p);
    return -2;
  }
  emax = e;

  if (kl == 0) {
    e = emin;      
    while (niter < max_iteration) {
      niter++;
      SetPotentialW(pot, e, orb->kappa, kv);
      SetVEffective(kl, kv, pot);
      i2 = TurningPoints(orb->n, e, pot);
      if (i2 < 0) {
	nodes = -1;
      } else {
	nodes = IntegrateRadial(p, e, pot, 0, 0.0, i2, 1.0, 0);
      }
      if (nodes < nr) break;
      e = e*2.0;
      if (nr == 0 && e < emin0) break;
    }
    emin = e;
    if (niter == max_iteration) {
      printf("Max iteration before finding correct nodes in RadialBasis %d %d\n",
	     nodes, nr);
      free(p);
      return -3;
    }
  }
  while (niter < max_iteration) {
    niter++;
    SetPotentialW(pot, e, orb->kappa, kv);
    SetVEffective(kl, kv, pot);
    i2 = TurningPoints(orb->n, e, pot);
    nodes = IntegrateRadial(p, e, pot, 0, 0.0, i2, 1.0, 0);
    if (nodes == nr) break;
    else if (nodes > nr) {
      emax = e;
      e = 0.5*(emin + e);
    } else {
      emin = e;
      e = 0.5*(emax + e);
    }
  }
  if (niter == max_iteration) {
    printf("Max iteration before finding correct nodes in RadialBasis %d %d\n",
	   nodes, nr);
    free(p);
    return -4;
  }

  de = 0.5*(emax - emin);
  while (de > ENEABSERR) {
    de *= 0.25;
    while (nodes == nr) {
      e += de;
      SetPotentialW(pot, e, orb->kappa, kv);
      SetVEffective(kl, kv, pot);
      i2 = TurningPoints(orb->n, e, pot);
      nodes = IntegrateRadial(p, e, pot, 0, 0.0, i2, 1.0, 0);
    }
    e -= de;
    nodes = nr;
  }
  niter = 0;
  if (fabs(pot->bqp) < 1E10) {
    bqp = (pot->bqp + orb->kappa)*FINE_STRUCTURE_CONST;
    bqp /= (2.0*pot->rad[pot->ib]);    
    emax = e;
    de = emax - emin;
    while (de > ENEABSERR) {
      de *= 0.25;
      while (nodes == nr) {
	e -= de;
	SetPotentialW(pot, e, orb->kappa, kv);
	SetVEffective(kl, kv, pot);
	i2 = TurningPoints(orb->n, e, pot);
	if (i2 < 0) {
	  nodes = -1;
	  break;
	}
	nodes = IntegrateRadial(p, e, pot, 0, 0.0, i2, 1.0, 0);
      }
      e += de;
      nodes = nr;
    }
    emin = e;
    e = 0.5*(emin+emax);
    while (niter < max_iteration) {
      niter++;
      SetPotentialW(pot, e, orb->kappa, kv);
      SetVEffective(kl, kv, pot); 
      i2 = TurningPoints(orb->n, e, pot); 
      if (i2 == pot->ib) {
	i2m1 = i2 - 1;
	i2m2 = i2 - 2;
	i2p1 = i2 + 1;
	i2p2 = i2 + 2;      
	bqp0 = DpDr(orb->kappa, kv, 0, e, pot, 0, -1, &orb->bqp0);
	nodes = IntegrateRadial(p, e, pot, 0, bqp0, i2p2, 1.0, 2);
	for (k = i2m2-1; k <= i2p2; k++) {
	  p[k] *= pot->dr_drho2[k];
	  p1 = e - pot->VT[kv][k];
	  p1 = sqrt(1.0 + FINE_STRUCTURE_CONST2*p1*0.5);
	  p[k] *= p1;
	}
	qo = (-4.0*p[i2m2-1] + 30.0*p[i2m2] - 120.0*p[i2m1]
	      + 40.0*p[i2] + 60.0*p[i2p1] - 6.0*p[i2p2])/120.0;
	p2 = qo/(p[i2]*pot->dr_drho[i2]) + orb->kappa/pot->rad[i2];
	p1 = e - pot->VT[kv][i2];
	p1 = 1.0 + 0.5*FINE_STRUCTURE_CONST2*p1;
	p2 *= 0.5*FINE_STRUCTURE_CONST/p1;
	qo = p2 - bqp;
	if (qo > EPS8) {
	  emin = e;
	  e = 0.5*(emax + e);
	} else if (qo < -EPS8) {
	  emax = e;
	  e = 0.5*(emin + e);
	} else {
	  break;
	}
      } else {
	i1 = pot->ib;
	for (k = i1-1; k <= i1+1; k++) {
	  p[k] = 0.5*FINE_STRUCTURE_CONST2*(e-pot->VT[kv][k]);
	  p[k] += 1.0;
	  if (k == i1) {
	    p2 = 2.0*p[k]*bqp/FINE_STRUCTURE_CONST;
	    p2 -= orb->kappa/pot->rad[i1];
	    p2 *= pot->dr_drho[i1];
	  }
	  p[k] = sqrt(p[k]);
	}
	p1 = 0.5*(p[i1+1] - p[i1-1])/p[i1];
	p2 -= p1;
	p1 = 0.5*(pot->dr_drho2[i1+1] - pot->dr_drho2[i1-1]);
	p1 /= pot->dr_drho2[i1];
	p2 -= p1;
	i2p2 = i2 + 2;      
	bqp0 = DpDr(orb->kappa, kv, 0, e, pot, 0, -1, &orb->bqp0);
	nodes = IntegrateRadial(p, e, pot, 0, bqp0, i2p2, 1.0, 2);
	for (i = 0; i <= i2p2; i++) {
	  p[i] *= pot->dr_drho2[i];
	}
	if (i2 > i1-10) {
	  i2 = i1 - 10;
	  i2 = LastMaximum(p, 0, i2);
	}
	i2m1 = i2 - 1;
	i2m2 = i2 - 2;
	i2p1 = i2 + 1;
	i2p2 = i2 + 2;
	p1 = p[i2m2]/pot->dr_drho2[i2m2];
	qo = (-4.0*p[i2m2-1] + 30.0*p[i2m2] - 120.0*p[i2m1]
	      + 40.0*p[i2] + 60.0*p[i2p1] - 6.0*p[i2p2])/120.0;
	qo /= p[i2];
	IntegrateRadial(p, e, pot, i2m2, p1, i1+1, p2, 1);
	for (i = i2m2; i <= pot->ib+1; i++) {
	  p[i] *= pot->dr_drho2[i];
	}
	qi = (6.0*p[i2m2] - 60.0*p[i2m1] - 40.0*p[i2] + 120.0*p[i2p1]
	      - 30.0*p[i2p2] + 4.0*p[i2p2+1])/120.0;
	qi /= p[i2];
	delta = qo - qi;      
	if (delta < -EPS8) {
	  emax = e;
	  e = 0.5*(emin + e);
	} else if (delta > EPS8) {
	  emin = e;
	  e = 0.5*(emax + e);
	} else {
	  break;
	}
      }
    }

    if (niter == max_iteration) {
      printf("Max iteration before finding solution in RadialBasis %d %d\n",
	     nodes, nr);
      free(p);
      return -7;
    }

    if (i2 == pot->ib) {      
      bqp0 = DpDr(orb->kappa, kv, 0, e, pot, 0, -1, &orb->bqp0);
      nodes = IntegrateRadial(p, e, pot, 0, bqp0, i2p2, 1.0, 2);
      for (i = 0; i <= i2p2; i++) {
	p[i] = p[i] * pot->dr_drho2[i];
      }      
    } else {
      i2p2 = pot->ib + 1;
    }
    if (IsOdd(nodes)) {
      for (i = 0; i <= i2p2; i++) {
	p[i] = -p[i];
      }
    }
    orb->ilast = pot->ib;
    orb->energy = e;
    orb->wfun = p;
    orb->qr_norm = 1.0;
    
    if (pot->flag == -1) {
      DiracSmall(orb, pot, i2p2, kv);
    }

    return 0;
  }

  niter = 0;
  if (i2 < pot->ib) {
    while (niter < max_iteration) {
      niter++;
      SetPotentialW(pot, e, orb->kappa, kv);
      SetVEffective(kl, kv, pot);
      i2 = TurningPoints(orb->n, e, pot);
      i2p2 = i2 + 2;      
      bqp0 = DpDr(orb->kappa, kv, 0, e, pot, 0, -1, &orb->bqp0);
      nodes = IntegrateRadial(p, e, pot, 0, bqp0, i2p2, 1.0, 2);
      for (i = 0; i <= i2p2; i++) {
	p[i] = p[i] * pot->dr_drho2[i];
      }
      i2 = LastMaximum(p, 0, i2);
      i2m1 = i2 - 1;
      i2m2 = i2 - 2;
      i2p1 = i2 + 1;
      i2p2 = i2 + 2;
      p1 = p[i2m2]/pot->dr_drho2[i2m2];
      p2 = p[i2];    
      qo = (-4.0*p[i2m2-1] + 30.0*p[i2m2] - 120.0*p[i2m1]
	    + 40.0*p[i2] + 60.0*p[i2p1] - 6.0*p[i2p2])/120.0;
      qo /= p2*pot->dr_drho[i2];
      norm2 = InnerProduct(0, i2, p, p, pot);
      nodes += IntegrateRadial(p, e, pot, i2m2, p1, pot->ib, 0.0, 0);
      for (i = i2m2; i < pot->maxrp; i++) {
	p[i] = p[i] * pot->dr_drho2[i];
      }
      qi = (6.0*p[i2m2] - 60.0*p[i2m1] - 40.0*p[i2] + 120.0*p[i2p1]
	    - 30.0*p[i2p2] + 4.0*p[i2p2+1])/120.0;
      qi /= p[i2]*pot->dr_drho[i2];
      fact = p2/p[i2];
      norm2 += fact*fact*InnerProduct(i2, pot->ib, p, p, pot);
      
      delta = 0.5*p2*p2*(qo - qi)/norm2;
      double e0 = e;
      e = e + delta;
      ep = EneTol(e);
      if (fabs(delta) < ep) break;
      if (niter > 20) {
	ep = niter-20;
	ep = 1.0 - 0.75*Min(100,ep)/100;
	e = e*ep + e0*(1-ep);
      }
    }
    nodes = 0;
    i = pot->ib-1;
    p1 = p[i];
    for (; i > 0; i--) {
      if (pot->rad[i] < pot->atom->rms0) break;
      p2 = fabs(p[i]);
      if (p2 > wave_zero) {
	if ((p1 > 0 && p[i] < 0) ||
	    (p1 < 0 && p[i] > 0)) {
	  nodes++;
	  p1 = p[i];
	}
      }
    }
    if (nodes != nr) {
      printf("RadialBasis: No. nodes changed in iteration\n");
      free(p);
      return -6;
    }    
    if (niter == max_iteration) {
      printf("Max iteration before finding solution in RadialBasis %d %d\n",
	     nodes, nr);
      free(p);
      return -7;
    }
    i2 = pot->ib;
    i2p2 = i2;
  } else {
    de *= 0.05;
    p2 = 0.0;
    while (niter < max_iteration) {
      niter++;
      SetPotentialW(pot, e, orb->kappa, kv);
      SetVEffective(kl, kv, pot);      
      bqp0 = DpDr(orb->kappa, kv, 0, e, pot, 0, -1, &orb->bqp0);
      nodes = IntegrateRadial(p, e, pot, 0, bqp0, i2, 1.0, 2);
      for (i = 0; i <= i2; i++) {
	p[i] = p[i] * pot->dr_drho2[i];
      }
      norm2 = InnerProduct(0, i2, p, p, pot);
      if (norm2 < p2) break;
      p2 = norm2;
      e += de;
    }
    if (niter == max_iteration) {
      printf("Max iteration before finding solution in RadialBasis %d %d\n",
	     nodes, nr);
      free(p);
      return -8;
    }
    e -= de;
    SetPotentialW(pot, e, orb->kappa, kv);
    SetVEffective(kl, kv, pot);
    i2 = pot->ib;
    i2p2 = i2 + 2;      
    bqp0 = DpDr(orb->kappa, kv, 0, e, pot, 0, -1, &orb->bqp0);
    IntegrateRadial(p, e, pot, 0, bqp0, i2p2, 1.0, 2);
    for (i = 0; i <= i2p2; i++) {
      p[i] = p[i] * pot->dr_drho2[i];
    }
    norm2 = InnerProduct(0, i2, p, p, pot);
  }

  ep = sqrt(norm2);
  fact = 1.0/ep;
  if (IsOdd(nodes)) {
    fact = -fact;
  }
  for (i = 0; i <= i2p2; i++) {
    p[i] *= fact;
  }
  for (i = i2p2+1; i < pot->maxrp; i++) {
    p[i] = 0.0;
  }
  
  orb->ilast = pot->ib;
  orb->energy = e;
  orb->wfun = p;
  orb->qr_norm = 1.0;
  if (pot->flag == -1) {
    DiracSmall(orb, pot, i2p2, kv);
  }
  return 0;
}
 
int RadialBound(ORBITAL *orb, POTENTIAL *pot) {
  double z, z0, e, de, ep, delta, emin, emax;
  double *p, norm2, fact, p1, p2, qo, qi, bqp;
  int i, kl, nr, nodes, niter;
  int i2, i2m1, i2m2, i2p1, i2p2;
  
  int kv = 0;
  if (pot->pse) kv = IdxVT(orb->kappa);
  orb->kv = kv;
  kl = orb->kappa;
  kl = (kl < 0)? (-kl-1):kl;

  if (kl < 0 || kl >= orb->n) {
    printf("Invalid orbital angular momentum, L=%d, %d %d\n", 
	   kl, orb->n, orb->kappa);
    return -1;
  }
  
  p = malloc(sizeof(double)*2*pot->maxrp);
  if (!p) return -1;

  nr = orb->n - kl - 1;
  z0 = pot->atom->atomic_number;
  z = (z0 - pot->N1);
  //if (z < 1) z = 1.0;
  double emin0 = 1.5*EnergyH(z0, orb->n, orb->kappa);
  if (kl > 0) {   
    SetPotentialW(pot, 0, orb->kappa, kv);
    SetVEffective(kl, kv, pot);
    emin = 1E30;
    for (i = 0; i < pot->maxrp; i++) {
      if (_veff[i] < emin) emin = _veff[i];
    }
    emin += ENEABSERR;
    if (emin < emin0) emin = emin0;
  } else {
    emin = EnergyH(z, orb->n, orb->kappa);
  }

  emax = EnergyH(Max(z,1.0), orb->n*10, orb->kappa);
  e = 0.5*emin;
  niter = 0;
  while (niter < max_iteration) {
    niter++;
    SetPotentialW(pot, e, orb->kappa, kv);
    SetVEffective(kl, kv, pot);
    i2 = TurningPoints(orb->n, e, pot);
    if (i2 > 0) {
      nodes = IntegrateRadial(p, e, pot, 0, 0.0, i2, 1.0, 0);
      if (nodes > nr) break;
      else if (nodes == nr && e >= emax) break;
    }
    e *= 0.5;
  }
  if (niter == max_iteration) {
    printf("Max iteration before finding correct nodes in RadialBound a: %d %d %d %d %d %g\n",
	   nodes, nr, i2, orb->n, orb->kappa, e);
    free(p);
    return -2;
  }
  emax = e;
  niter = 0;
  if (kl == 0) {
    e = emin;   
    while (niter < max_iteration) {
      niter++;
      SetPotentialW(pot, e, orb->kappa, kv);
      SetVEffective(kl, kv, pot);
      i2 = TurningPoints(orb->n, e, pot);
      if (i2 < 0) {
	nodes = -1;
      } else {
	nodes = IntegrateRadial(p, e, pot, 0, 0.0, i2, 1.0, 0);
      }
      if (nodes < nr) break;
      if (nr == 0 && e < emin0) break;
      e *= 1.25;
    }
    emin = e;
    if (niter == max_iteration) {
      printf("Max iteration before finding correct nodes in RadialBound b: %d %d\n",
	     nodes, nr);
      free(p);
      return -3;
    }
  }
  while (niter < max_iteration) {
    niter++;
    SetPotentialW(pot, e, orb->kappa, kv);
    SetVEffective(kl, kv, pot);
    i2 = TurningPoints(orb->n, e, pot);
    if (i2 < 0) {
      nodes = -1;
    } else {
      bqp = DpDr(orb->kappa, kv, 0, e, pot, 0, -1, &orb->bqp0);
      nodes = IntegrateRadial(p, e, pot, 0, bqp, i2, 1.0, 2);
    }
    if (nodes == nr) break;
    else if (nodes > nr) {
      emax = e;
      e = 0.5*(emin + e);
    } else {
      emin = e;
      e = 0.5*(emax + e);
    }
  }
  if (niter == max_iteration) {
    printf("Max iteration before finding correct nodes in RadialBound c: %d %d %d %d %g %g %g\n",
	   orb->n, orb->kappa, nodes, nr, e, emin, emax);
    free(p);
    return -4;
  }
  
  de = 0.5*(emax - emin);
  while (de > ENEABSERR) {
    de *= 0.25;
    while (nodes == nr) {
      e += de;
      SetPotentialW(pot, e, orb->kappa, kv);
      SetVEffective(kl, kv, pot);
      i2 = TurningPoints(orb->n, e, pot);
      bqp = DpDr(orb->kappa, kv, 0, e, pot, 0, -1, &orb->bqp0);
      nodes = IntegrateRadial(p, e, pot, 0, bqp, i2, 1.0, 2);
    }
    e -= de;
    nodes = nr;
  }

  niter = 0;
  de = 0.0;
  while (niter < max_iteration) {
    niter++;
    SetPotentialW(pot, e, orb->kappa, kv);
    SetVEffective(kl, kv, pot);
    i2 = TurningPoints(orb->n, e, pot);
    i2p2 = i2 + 2;
    bqp = DpDr(orb->kappa, kv, 0, e, pot, 0, -1, &orb->bqp0);
    nodes = IntegrateRadial(p, e, pot, 0, bqp, i2p2, 1.0, 2); 
    for (i = 0; i <= i2p2; i++) {
      p[i] = p[i] * pot->dr_drho2[i];
    }
    i2 = LastMaximum(p, 0, i2);
    //orb->im = i2;
    i2m1 = i2 - 1;
    i2m2 = i2 - 2;
    i2p1 = i2 + 1;
    i2p2 = i2 + 2;
    p1 = p[i2m2]/pot->dr_drho2[i2m2];
    p2 = p[i2];    
    qo = (-4.0*p[i2m2-1] + 30.0*p[i2m2] - 120.0*p[i2m1]
	  + 40.0*p[i2] + 60.0*p[i2p1] - 6.0*p[i2p2])/120.0;
    qo /= p2*pot->dr_drho[i2];
    norm2 = InnerProduct(0, i2, p, p, pot);
    //bqp = DpDr(orb->kappa, kv, pot->maxrp-1, e, pot, 0, -2, &orb->bqp1);
    IntegrateRadial(p, e, pot, i2m2, p1, pot->maxrp-1, 0, 0);
    for (i = i2m2; i < pot->maxrp; i++) {
      p[i] = p[i] * pot->dr_drho2[i];
    }
    qi = (6.0*p[i2m2] - 60.0*p[i2m1] - 40.0*p[i2] + 120.0*p[i2p1]
	  - 30.0*p[i2p2] + 4.0*p[i2p2+1])/120.0;
    qi /= p[i2]*pot->dr_drho[i2];
    fact = p2/p[i2];
    norm2 += fact*fact*InnerProduct(i2, pot->maxrp-1, p, p, pot);

    delta = 0.5*p2*p2*(qo - qi)/norm2;
    double e0 = e;
    e = e + delta;
    ep = EneTol(e);
    if (fabs(delta) < ep) break;
    if (niter > 10) {
      ep = niter-10;
      ep = 1.0 - 0.75*Min(100,ep)/100;
    }
    e = e*ep + e0*(1-ep);
  }
  if (niter == max_iteration) {
    printf("Max iteration reached in RadialBound: %d %d\n", orb->n, orb->kappa);
    Abort(1);
  }
  ep = sqrt(norm2);
  fact = 1.0/ep;
  if (IsOdd(nodes)) {
    fact = -fact;
  }
  ep *= wave_zero;     
  for (i = pot->maxrp-1; i >= 0; i--) {
    if (fabs(p[i]) > ep) break;
  }
  if (IsEven(i)) i++;
  orb->ilast = i;
  for (i = 0; i < pot->maxrp; i++) {    
    p[i] *= fact;
  }

  nodes = 0;
  for (i = orb->ilast; i >= 0; i--) {
    if (fabs(p[i]) > EPS5) break;
  }
  p1 = p[i];
  for (; i > 0; i--) {
    if (pot->rad[i] < pot->atom->rms0) break;
    p2 = fabs(p[i]);
    if (p2 > wave_zero) {
      if ((p1 > 0 && p[i] < 0) ||
	  (p1 < 0 && p[i] > 0)) {
	nodes++;
	p1 = p[i];
      }
    }
  }
  if (nodes != nr) {
    printf("RadialBound: No. nodes changed in iteration %d %d\n", nodes, nr);
    free(p);
    return -6;
  }

  orb->energy = e;
  orb->wfun = p;
  orb->qr_norm = 1.0;
  if (pot->flag == -1) {
    DiracSmall(orb, pot, -1, kv);
  }

  return 0;
}

int RadialRydberg(ORBITAL *orb, POTENTIAL *pot) {
#define ME 20
  double z, e, e0;
  int i, kl, niter, ierr;
  double x0, pp, qq, ppi, qqi;
  int i2, i2p, i2m, i2p2, i2m2, nodes, nr;
  double qo, qi, norm2, delta, zp, *p;
  double ep, p1, p2, fact, bqp;
  double en[ME], dq[ME], dn, zero=0.0;
  int j, np, nme, one=1;

  int kv = 0;
  if (pot->pse) kv = IdxVT(orb->kappa);
  orb->kv = kv;
  z = GetResidualZ();
  if (z < 1) z = 1;
  kl = orb->kappa;
  kl = (kl < 0)? (-kl-1):kl;
  if (kl < 0 || kl >= orb->n) {
    printf("Invalid orbital angular momentum, L=%d, %d %d\n", 
	   kl, orb->n, orb->kappa);
    return -1;
  }
  e = EnergyH(z, orb->n, orb->kappa);
  p = malloc(sizeof(double)*2*pot->maxrp);
  if (!p) return -1;
  nr = orb->n - kl - 1;

  SetPotentialW(pot, e, orb->kappa, kv);
  SetVEffective(kl, kv, pot);
  i2 = TurningPoints(orb->n, e, pot);
  if (i2 < pot->maxrp - 1.5*pot->asymp) {
    niter = 0;
    dn = orb->n;
    dn = z*z/(dn*dn*dn);
    while (1) {
      SetPotentialW(pot, e, orb->kappa, kv);
      SetVEffective(kl, kv, pot);
      i2 = TurningPoints(orb->n, e, pot);
      if (i2 < 0) {
	printf("The orbital angular momentum = %d too high\n", kl);
	return -2;
      }
      i2p2 = i2 + 2;
      nodes = IntegrateRadial(p, e, pot, 0, 0, i2p2, 1.0, 0);
      if (nodes != nr) {
	e += 0.5*(nr-nodes)*dn;
      } else {
	break;
      }
    }
    dn *= 0.01;
    while (nodes == nr) {
      e += dn;
      SetPotentialW(pot, e, orb->kappa, kv);
      SetVEffective(kl, kv, pot);
      i2 = TurningPoints(orb->n, e, pot);
      if (i2 < 0) {
	printf("The orbital angular momentum = %d too high\n", kl);
	return -2;
      }
      i2p2 = i2 + 2;
      nodes = IntegrateRadial(p, e, pot, 0, 0.0, i2p2, 1.0, 0);
    }
    e -= 1.25*dn;
    while (niter < max_iteration) {
      niter++;
      SetPotentialW(pot, e, orb->kappa, kv);
      SetVEffective(kl, kv, pot);
      i2 = TurningPoints(orb->n, e, pot);
      if (i2 < 0) {
	printf("The orbital angular momentum = %d too high\n", kl);
	return -2;
      }
      i2p2 = i2 + 2;
      bqp = DpDr(orb->kappa, kv, 0, e, pot, 0, -1, &orb->bqp0);
      nodes = IntegrateRadial(p, e, pot, 0, bqp, i2p2, 1.0, 2);
      for (i = 0; i <= i2p2; i++) {
	p[i] = p[i] * pot->dr_drho2[i];
      }
      i2 = LastMaximum(p, 0, i2);
      i2p = i2 + 1;
      i2m = i2 - 1;
      i2p2 = i2 + 2;
      i2m2 = i2 - 2;
      p1 = p[i2m2]/pot->dr_drho2[i2m2];
      p2 = p[i2];
      qo = (-4.0*p[i2m2-1] + 30.0*p[i2m2] - 120.0*p[i2m]
	    + 40.0*p[i2] + 60.0*p[i2p] - 6.0*p[i2p2])/120.0;
      qo /= p2*pot->dr_drho[i2];
      norm2 = InnerProduct(0, i2, p, p, pot);
      
      IntegrateRadial(p, e, pot, i2m2, p1, pot->maxrp-1, 0.0, 0);
      for (i = i2m2; i < pot->maxrp; i++) {
	p[i] = p[i] * pot->dr_drho2[i];
      }
      qi = (6.0*p[i2m2] - 60.0*p[i2m] - 40.0*p[i2] + 120.0*p[i2p]
	    - 30.0*p[i2p2] + 4.0*p[i2p2+1])/120.0;
      qi /= p[i2]*pot->dr_drho[i2];
      fact = p2/p[i2];
      norm2 += fact*fact*InnerProduct(i2, pot->maxrp-1, p, p, pot);
      
      delta = 0.5*p2*p2*(qo - qi)/norm2;
      e0 = e;
      e = e + delta;
      ep = EneTol(e);
      if (fabs(delta) < ep) {
	fact = 1.0/sqrt(norm2);
	if (IsOdd(nodes)) fact = -fact;
	for (i = 0; i < pot->maxrp; i++) {
	  p[i] *= fact;
	}
	break;
      }
      if (niter > 20) {
	ep = niter-20;
	ep = 1.0 - 0.75*Min(100,ep)/100;
	e = e*ep + e0*(1-ep);
      }
    }
    if (niter == max_iteration) {
      printf("Max iteration reached in RadialRydberg\n");
      free(p);
      return -3;
    }    
    pp = 0.0;
  } else {
    if (pot->N <= 1) j = 3;
    else j = 1;
    i2p2 = pot->maxrp - j*pot->asymp - 5;
    for (i2 = pot->r_core; i2 < i2p2; i2++) {
      if (_veff[i2+1] > _veff[i2]) break;
    }
    i2 += j*pot->asymp;    
    i2p2 = i2 + 2;
    if (i2p2 >= pot->maxrp) {
      i2p2 = pot->maxrp-1;
      i2 = i2p2-2;
    }
    nodes = IntegrateRadial(p, e, pot, 0, 0.0, i2p2, 1.0, 0);
    for (i = 0; i <= i2p2; i++) {
      p[i] = p[i] * pot->dr_drho2[i];
    }
    i2 = LastMaximum(p, pot->r_core, i2);
    i2p = i2 + 1;
    i2m = i2 - 1;
    i2p2 = i2 + 2;
    i2m2 = i2 - 2;
    dn = 1.2/(ME-1.0);
    en[0] = orb->n - 0.95;
    for (j = 1; j < ME; j++) {
      en[j] = en[j-1] + dn;
    }
    for (j = 0; j < ME; j++) {
      e = EnergyH(z, en[j], orb->kappa);
      en[j] = e;
      SetPotentialW(pot, e, orb->kappa, kv);
      SetVEffective(kl, kv, pot);
      bqp = DpDr(orb->kappa, kv, 0, e, pot, 0, -1, &orb->bqp0);
      nodes = IntegrateRadial(p, e, pot, 0, bqp, i2p2, 1.0, 2);
      p2 = p[i2];
      qo = (-4.0*p[i2m2-1] + 30.0*p[i2m2] - 120.0*p[i2m]
	    + 40.0*p[i2] + 60.0*p[i2p] - 6.0*p[i2p2])/120.0;
      qo /= p2;
      zp = FINE_STRUCTURE_CONST2*e;
      x0 = pot->rad[i2];
      ierr = 1;
      DCOUL(z, e, orb->kappa, x0, &pp, &qq, &ppi, &qqi, &ierr);
      p1 = qq/pp;
      qi = DpDr(orb->kappa, kv, i2, e, pot, p1, 1, NULL);
      delta = qo-qi;
      dq[j] = delta;
    }
    for (j = 0; j < ME-1; j++) {
      if (dq[j] > 0 && dq[j+1] < dq[j]) break;
    }
    i = j;
    for (; j < ME-1; j++) {
      if (dq[j+1] >= dq[j]) break;
    }
    nme = j - i + 1;
    for (np = i; np <= j; np++) {
      dq[np] = -dq[np];
    }
    np = 3;
    UVIP3P(np, nme, &(dq[i]), &(en[i]), one, &zero, &e);
    i2p2 = pot->maxrp-1;
    bqp = DpDr(orb->kappa, kv, 0, e, pot, 0, -1, &orb->bqp0);
    nodes = IntegrateRadial(p, e, pot, 0, bqp, i2p2, 1.0, 2);
    for (i = 0; i <= i2p2; i++) {
      p[i] = p[i] * pot->dr_drho2[i];
    }
    i2 = LastMaximum(p, pot->r_core, i2);
    zp = FINE_STRUCTURE_CONST2*e;
    x0 = pot->rad[i2];
    ierr = 1;
    DCOUL(z, e, orb->kappa, x0, &pp, &qq, &ppi, &qqi, &ierr);
    norm2 = pp;
    fact = norm2/p[i2];
    if (IsOdd(nodes)) {
      fact = -fact;
    }
    for (i = 0; i <= i2p2; i++) {
      p[i] *= fact;
    }
  }

  i = pot->maxrp-1;
  orb->ilast = i;
  orb->energy = e;
  orb->wfun = p;
  
  e0 = InnerProduct(0, pot->maxrp-1, p, p, pot);
  orb->qr_norm = 1.0/e0;

  if (pot->flag == -1) {
    DiracSmall(orb, pot, -1, kv);
    if (pp != 0) {
      fact = fabs(pp/p[i2]);
      for (i = 0; i < pot->maxrp-1; i++) {
	orb->wfun[i] *= fact;
	orb->wfun[i+pot->maxrp] *= fact;
      }
    }
  }

  return 0;
#undef ME
}
  
int RadialFreeInner(ORBITAL *orb, POTENTIAL *pot) {
  int i, kl, nodes;
  int i2;
  double *p, e, bqp0;

  int kv = 0;
  if (pot->pse) kv = IdxVT(orb->kappa);
  orb->kv = kv;
  e = orb->energy;
  kl = orb->kappa;
  if (orb->kappa == 0) {
    printf("Kappa == 0 in RadialFreeInner\n");
    return -1;
  }
  if (pot->ib1 <= 0 && pot->ib <= 0) {
    printf("Boundary not set\n");
    return -2;
  }
  SetPotentialW(pot, e, kl, kv);
  kl = (kl < 0)? (-kl-1):kl;  
  p = malloc(2*pot->maxrp*sizeof(double));
  if (!p) return -1;
  SetVEffective(kl, kv, pot);

  i2 = pot->ib1;
  i2 = Max(i2,pot->ib) + 2;      
  bqp0 = DpDr(orb->kappa, kv, 0, e, pot, 0, -1, &orb->bqp0);
  nodes = IntegrateRadial(p, e, pot, 0, bqp0, i2, 1.0, 2);
  for (i = i2; i >= 0; i--) {
    p[i] *= pot->dr_drho2[i];
  }
  orb->ilast = i2;
  orb->wfun = p;
  orb->phase = NULL;
  
  orb->qr_norm = 1.0;
  
  if (pot->flag == -1) {
    DiracSmall(orb, pot, i2, kv);
  }

  return 0;
}

/* note that the free states are normalized to have asymptotic 
   amplitude of 1/sqrt(k), */
int RadialFree(ORBITAL *orb, POTENTIAL *pot) {
  int i, kl, nodes;
  int i2, i2p, i2m, i2p2, i2m2;
  double *p, po, qo, e, po1, bqp0;
  double dfact, da, cs, si, phase0;

  int kv = 0;
  if (pot->pse) kv = IdxVT(orb->kappa);
  orb->kv = kv;
  e = orb->energy;
  if (e < 0.0) { 
    printf("Energy < 0 in Free\n");
    return -1;
  }
  kl = orb->kappa;
  if (orb->kappa == 0) {
    printf("Kappa == 0 in Free\n");
    return -1;
  }
  SetPotentialW(pot, e, kl, kv);
  kl = (kl < 0)? (-kl-1):kl;  
  p = malloc(2*pot->maxrp*sizeof(double));
  if (!p) return -1;
  SetVEffective(kl, kv, pot);
  
  i2 = TurningPoints(0, e, pot);
  i2m = i2 - 1;
  i2p = i2 + 1;
  i2m2 = i2 - 2;
  i2p2 = i2 + 2;
      
  bqp0 = DpDr(orb->kappa, kv, 0, e, pot, 0, -1, &orb->bqp0);
  nodes = IntegrateRadial(p, e, pot, 0, bqp0, i2p2, 1.0, 2);
  for (i = i2p2; i >= 0; i--) {
    p[i] *= pot->dr_drho2[i];
  }
  dfact = 1.0 / pot->dr_drho[i2];
  qo = (-4.0*p[i2m2-1] + 30.0*p[i2m2] - 120.0*p[i2m]
	+ 40.0*p[i2] + 60.0*p[i2p] - 6.0*p[i2p2])/120.0;
  qo *= dfact;
  po = p[i2];
  po1 = p[i2p];

  da = Amplitude(p, e, orb->kappa, pot, i2);
  cs = p[i2] * qo - da*po;
  si = po / p[i2];
  dfact = (si*si + cs*cs);
  dfact = 1.0/sqrt(dfact);
  phase0 = atan2(si, cs);
  if (phase0 < 0) phase0 += TWO_PI;
    
  if (IsOdd(nodes)) {
    phase0 = (phase0 < PI)?(phase0 + PI):(phase0-PI);
    dfact = -dfact;
  }
  Phase(p, pot, i2, phase0);
  
  p[i2] = po;
  p[i2p] = po1;
  for (i = 0; i < i2p2; i++) {
    p[i] *= dfact;
  }
    
  orb->ilast = i2p;
  orb->wfun = p;
  orb->phase = NULL;

  orb->qr_norm = 1.0;
  
  if (pot->flag == -1) {
    DiracSmall(orb, pot, -1, kv);
  }

  return 0;
}

/*
** The storage of amplitude and phase in the wfun array is arranged:
** large component: large[i]=amplitude_i, large[i+1]=phase_i,
** small component: small[i]=a_cos_i, small[i+1]=a_sin_i,
** so: P[i] = large[i]*sin(large[i+1])
** and Q[i] = small[i]*cos(large[i+1])+small[i+1]*sin(large[i+1])
*/
int DiracSmall(ORBITAL *orb, POTENTIAL *pot, int i2, int kv) {
  int i, i0, i1, ia, kappa, imax, irn, is;
  double xi, e, *p, *q, a, b, c;

  if (orb->n >= 0) i0 = 0;
  else i0 = pot->ib-1;
  e = orb->energy;
  kappa = orb->kappa;
  p = orb->wfun;
  q = p + pot->maxrp;
  if (i2 < 0) i2 = orb->ilast;
  i1 = orb->ilast+1;

  if (i2 > i0) {
    is = i0;
    if (i0 == 0) {
      a = 0;
      for (is = i0; is < i1; is++) {
	b = fabs(p[is]);
	if (b > a) a = b;
      }
      for (is = i0; is < i1; is++) {
	if (fabs(p[is])>1e-30*a &&
	    ((p[is] > 0 && p[is+1] > 0) ||
	     (p[is] < 0 && p[is+1] < 0))) break;
      }
    }
    if (i0 == 0 && orb->kappa < 0) {
      imax = -1;
      irn = -1;
      for (i = i0; i <= i2; i++) {
	if (kv >= 0) {
	  xi = e - pot->VT[kv][i];
	} else {
	  xi = e;
	}
	_dwork3[i] = xi;
	_dwork[i] = 1 + xi*FINE_STRUCTURE_CONST2*0.5;
	if (orb->isol != 3) p[i] = sqrt(_dwork[i])*p[i];
	if (irn < 0 && pot->rad[i] > pot->atom->rms0) {
	  irn = i-1;
	}
	if (irn > 0 && i >= irn && i >= is && imax < 0 &&
	    ((p[i-1] > wave_zero && p[i] < p[i-1]) ||
	     (p[i-1] < -wave_zero && p[i] > p[i-1]))) {
	  imax = i-1;
	}
      }
      if (imax < 0) {
	imax = i2-2;
      } else {
	imax = Max(imax, irn);
	imax = Max(imax, is+10);
	imax = Min(imax, i2-2);
      }
      orb->pdx = -kappa;
      if (pot->atom->z1 <= 0) {
	orb->pdx = sqrt(orb->pdx*orb->pdx -
			FINE_STRUCTURE_CONST2*pot->atom->atomic_number);
      }
      for (i = is; i <= imax; i++) {
	a = log(p[i+1]/p[i])/log(pot->rad[i+1]/pot->rad[i]);
	if (isnan(a)) continue;
	if (fabs(p[i]) > wave_zero ||
	    fabs(a-orb->pdx) < 1e-4*orb->pdx ||
	    a < orb->pdx) break;
      }
      is = i;
      if (p[is] > 0) {
	b = log(p[is])-a*log(pot->rad[is]);
      } else {
	b = log(-p[is])-a*log(pot->rad[is]);
      }
      for (i = i0; i < is; i++) {
	p[i] = b + a*log(pot->rad[i]);
      }
      orb->pdx = a;
      double a1 = a;
      if (pot->atom->z1 > 0) {
	a1 = 1+a;
      }
      for (i = i0; i <= imax; i++) {
	_dwork1[i] = a1*log(pot->rad[i]);
      }
      for (i = i0; i <= imax; i++) {
	xi = _dwork3[i];
	if (i < is) {
	  if (p[is] > 0) {
	    _dwork2[i] = -exp(p[i]-_dwork1[i])*FINE_STRUCTURE_CONST;
	  } else {
	    _dwork2[i] = exp(p[i]-_dwork1[i])*FINE_STRUCTURE_CONST;
	  }
	} else {
	  if (p[i] > 0) {
	    _dwork2[i] = -exp(log(p[i]*FINE_STRUCTURE_CONST)-_dwork1[i]);
	  } else if (p[i] < 0) {
	    _dwork2[i] = exp(log(-p[i]*FINE_STRUCTURE_CONST)-_dwork1[i]);
	  } else {
	    _dwork2[i] = 0;
	  }
	}
	_dwork2[i] *= xi*pot->dr_drho[i];
	_dwork3[i] = (a1-kappa)/pot->rad[i];
	_dwork3[i] *= pot->dr_drho[i];
      }
      if (i0 < is) {
	q[i0] = exp(p[i0]-_dwork1[i0])*orb->bqp0;
	if (p[is] < 0) q[i0] = -q[i0];
      } else {
	q[i0] = (p[i0]/exp(_dwork1[i0]))*orb->bqp0;
      }
      //printf("%d %d %d %d %g %g %g %g %g %g\n", i0, is, imax, i1, a1, orb->bqp0, (pot->rad[is]/pot->rad[i0]), _dwork1[is], q[is], q[is]*exp(_dwork1[is]-log(p[is])));
      for (i = i0+1; i <= imax; i++) {
	q[i] = 0.5*(_dwork2[i-1]+_dwork2[i]) + (1-0.5*_dwork3[i-1])*q[i-1];
	q[i] /= 1 + 0.5*_dwork3[i];	
      }
      for (i = i0; i <= imax; i++) {
	q[i] *= exp(_dwork1[i]);
      }
      for (i = i0; i < is; i++) {
	p[i] = exp(p[i]);
	if (p[is] < 0) p[i] = -p[i];
      }
      if (imax < i2) {
	a = q[imax];
	ia = Max(i0, imax-2);
	for (i = ia; i <= i2; i++) {
	  _dwork2[i] = 1.0/(24.0*pot->dr_drho[i]);
	  _dwork1[i] = p[i];
	}
	for (i = imax; i < i1; i++) {
	  if (i == ia) {
	    b = -50.0*_dwork1[imax];
	    b += 96.0*_dwork1[imax+1];
	    b -= 72.0*_dwork1[imax+2];
	    b += 32.0*_dwork1[imax+3];
	    b -= 6.0 *_dwork1[imax+4];
	  } else if (i == ia+1) {
	    b = -6.0*_dwork1[imax];
	    b -= 20.0*_dwork1[imax+1];
	    b += 36.0*_dwork1[imax+2];
	    b -= 12.0*_dwork1[imax+3];
	    b += 2.0 *_dwork1[imax+4];
	  } else if (i == i2) {
	    b = -50.0*_dwork1[i];
	    b += 96.0*_dwork1[i-1];
	    b -= 72.0*_dwork1[i-2];
	    b += 32.0*_dwork1[i-3];
	    b -= 6.0 *_dwork1[i-4];
	    b = -b; 
	  } else if (i == i2-1) {
	    b = -6.0*_dwork1[i+1];
	    b -= 20.0*_dwork1[i];
	    b += 36.0*_dwork1[i-1];
	    b -= 12.0*_dwork1[i-2];
	    b += 2.0 *_dwork1[i-3];
	    b = -b;
	  } else {
	    b = 2.0*(_dwork1[i-2] - _dwork1[i+2]);
	    b += 16.0*(_dwork1[i+1] - _dwork1[i-1]);
	  }
	  b *= _dwork2[i];
	  b = b + _dwork1[i]*kappa/pot->rad[i];
	  b /= (2.0*_dwork[i]);
	  q[i] = b*FINE_STRUCTURE_CONST;
	}
	//printf("qm: %d %d %d %d %g %g %g\n", orb->n, orb->kappa, is, imax, pot->rad[imax], a, q[imax]);
	a = q[imax]/a;
	for (i = i0; i < imax; i++) {
	  q[i] *= a;
	}
      }
    } else {
      imax = -1;
      irn = -1;
      for (i = i0; i <= i2; i++) {
	if (kv >= 0) {
	  xi = e - pot->VT[kv][i];
	} else {
	  xi = e;
	}
	xi = xi*FINE_STRUCTURE_CONST2*0.5;
	_dwork[i] = 1.0 + xi;
	if (orb->isol != 3) p[i] = sqrt(_dwork[i])*p[i];
	_dwork2[i] = 1.0/(24.0*pot->dr_drho[i]);
	if (irn < 0 && pot->rad[i] > pot->atom->rms0) {
	  irn = i-1;
	}
	if (irn > 0 && i >= irn && i >= is && imax < 0 &&
	    ((p[i-1] > wave_zero && p[i] < p[i-1]) ||
	     (p[i-1] < -wave_zero && p[i] > p[i-1]))) {
	  imax = i-1;
	}
      }
      if (imax < 0) {
	imax = i2-2;
      } else {
	imax = Max(imax, irn);
	imax = Max(imax, is+10);
	imax = Min(imax, i2-2);
      }
      if (i0 == 0) {
	if (pot->atom->z1 > 0) {
	  orb->pdx = kappa < 0?-kappa:kappa+1;
	} else {
	  orb->pdx = sqrt(kappa*kappa -
			  FINE_STRUCTURE_CONST2*pot->atom->atomic_number);
	}
	for (i = is; i <= imax; i++) {
	  a = log(p[i+1]/p[i])/log(pot->rad[i+1]/pot->rad[i]);
	  if (fabs(p[i]) > wave_zero ||
	      fabs(a-orb->pdx) < 1e-4*orb->pdx ||
	      a < orb->pdx) break;      
	}
	is = i;
	if (p[is] > 0) {
	  b = log(p[is]) - a*log(pot->rad[is]);
	} else {
	  b = log(-p[is]) - a*log(pot->rad[is]);
	}
	for (i = i0; i < is; i++) {
	  p[i] = b + a*log(pot->rad[i]);
	}
	orb->pdx = a;
	for (i = i0; i <= i2; i++) {
	  _dwork3[i] = a*log(pot->rad[i]/(pot->rad[i]+pot->rad[imax]));
	  if (i < is) {
	    if (p[is] > 0) {
	      _dwork1[i] = exp(p[i] - _dwork3[i]);
	    } else {
	      _dwork1[i] = -exp(p[i] - _dwork3[i]);
	    }
	  } else {
	    if (p[i] > 0) {
	      _dwork1[i] = exp(log(p[i])-_dwork3[i]);
	    } else if (p[i] < 0) {
	      _dwork1[i] = -exp(log(-p[i])-_dwork3[i]);
	    } else {
	      _dwork1[i] = 0;
	    }
	  }
	}
	for (i = i0; i < is; i++) {
	  p[i] = exp(p[i]);
	  if (p[is] < 0) p[i] = -p[i];
	}
      } else {
	is = i0;
	for (i = i0; i <= i2; i++) {
	  _dwork3[i] = 0.0;
	  _dwork1[i] = p[i];
	}
      }
      double afs = log(FINE_STRUCTURE_CONST);
      for (i = i0; i < i1; i++) {
	if (i == i0) {
	  b = -50.0*_dwork1[i0];
	  b += 96.0*_dwork1[i0+1];
	  b -= 72.0*_dwork1[i0+2];
	  b += 32.0*_dwork1[i0+3];
	  b -= 6.0 *_dwork1[i0+4];
	} else if (i == i0+1) {
	  b = -6.0*_dwork1[i0];
	  b -= 20.0*_dwork1[i0+1];
	  b += 36.0*_dwork1[i0+2];
	  b -= 12.0*_dwork1[i0+3];
	  b += 2.0 *_dwork1[i0+4];
	} else if (i == i2) {
	  b = -50.0*_dwork1[i];
	  b += 96.0*_dwork1[i-1];
	  b -= 72.0*_dwork1[i-2];
	  b += 32.0*_dwork1[i-3];
	  b -= 6.0 *_dwork1[i-4];
	  b = -b; 
	} else if (i == i2-1) {
	  b = -6.0*_dwork1[i+1];
	  b -= 20.0*_dwork1[i];
	  b += 36.0*_dwork1[i-1];
	  b -= 12.0*_dwork1[i-2];
	  b += 2.0 *_dwork1[i-3];
	  b = -b;
	} else {
	  b = 2.0*(_dwork1[i-2] - _dwork1[i+2]);
	  b += 16.0*(_dwork1[i+1] - _dwork1[i-1]);
	}
	b *= _dwork2[i];
	c = a*pot->rad[imax]/(pot->rad[imax]+pot->rad[i]) + kappa;
	b = (b + _dwork1[i]*c/pot->rad[i]);
	b /= (2.0*_dwork[i]);
	if (b > 0) {
	  q[i] = exp(log(b)+afs+_dwork3[i]);
	} else if (b < 0) {
	  q[i] = -exp(log(-b)+afs+_dwork3[i]);
	} else {
	  q[i] = 0;
	}
      }
    }
  }
  if (orb->n != 0) {
    for (i = i1; i < pot->maxrp; i++) {
      p[i] = 0;
      p[i+pot->maxrp] = 0.0;
    }
    a = InnerProduct(i0, i1-1, p+pot->maxrp, p+pot->maxrp, pot);
    b = InnerProduct(i0, i1-1, p, p, pot);
    a *= orb->qr_norm;
    b *= orb->qr_norm;
    a = sqrt(a+b);
    orb->qr_norm = a/sqrt(b);
    a = 1.0/a;
    for (i = i0; i < i1; i++) {
      p[i] *= a;
      p[i+pot->maxrp] *= a;
    }
    for (i = 0; i < i0; i++) {
      p[i] = 0.0;
      p[i+pot->maxrp] = 0.0;
    }
    for (i = i1+pot->maxrp; i < 2*pot->maxrp; i++) {
      p[i] = 0.0;
    }
    if (i0 > 0) {
      p[i0-1] = 0.0;
      p[i0+pot->maxrp-1] = 0.0;
    }
    return 0;
  }
  
  for (i = i1; i < pot->maxrp; i += 2) {
    if (kv >= 0) {
      xi = e - pot->VT[kv][i];
    } else {
      xi = e;
    }
    xi = xi*FINE_STRUCTURE_CONST2*0.5;
    _dwork[i] = 1.0 + xi;
    _dwork1[i] = sqrt(_dwork[i])*p[i];
    _dwork2[i] = 0.25/(pot->dr_drho[i]);
  }

  for (i = i1; i < pot->maxrp; i += 2) {
    if (i == i1) {
      b = -3.0*_dwork1[i] + 4.0*_dwork1[i+2] - _dwork1[i+4];
    } else if (i == pot->maxrp-2) {
      b = -3.0*_dwork1[i] + 4.0*_dwork1[i-2] - _dwork1[i-4];
      b = -b;
    } else {
      b = _dwork1[i+2] - _dwork1[i-2];
    }
    b *= _dwork2[i];
    b += _dwork1[i]*kappa/pot->rad[i];
    a = _dwork1[i]/(p[i]*p[i]);
    p[i] = _dwork1[i];
    q[i] = FINE_STRUCTURE_CONST*a/(2.0*_dwork[i]);
    q[i+1] = FINE_STRUCTURE_CONST*b/(2.0*_dwork[i]);
  }

  b = FINE_STRUCTURE_CONST2*orb->energy;
  orb->qr_norm = sqrt((1.0 + b)/(1.0 + 0.5*b));
  return 0;
}

void DerivODE(int *neq, double *t, double *y, double *ydot) {
  double w0, w;
  double t0, s, e;

  t0 = y[2];
  w0 = y[3];
  s = y[4];
  e = y[5];

  w = w0 + (*t - t0)*s;
  w = 2.0*(e - w/(*t));
  
  ydot[0] = y[1];
  ydot[1] = 1.0/(y[0]*y[0]*y[0]) - y[0]*w;
}

/* provide fortran access with cfortran.h */
FCALLSCSUB4(DerivODE, DERIVODE, derivode, PINT, PDOUBLE, DOUBLEV, DOUBLEV)
  
double Amplitude(double *p, double e, int ka, POTENTIAL *pot, int i0) {
  int i, n, kl1;
  double a, b, xi, r2, r3;
  double z, dk, r0, r1, r, w, v1;
  
  int neq, itol, itask, istate, iopt, lrw, iwork[22], liw, mf;
  double y[6], rtol, atol[2], *rwork;

  n = pot->maxrp-1;
  z = GetResidualZ();
  if (z < 1) z = 1;
  kl1 = ka*(ka+1);
  r1 = pot->rad[n];
  _dwork[0] = r1;
  _dwork1[0] = _veff[n];
  dk = sqrt(2.0*e*(1.0+0.5*FINE_STRUCTURE_CONST2*e));
  dk = EPS5*e*dk;
  for (i = 1; i < pot->maxrp; i++) {
    r = _dwork[i-1]*1.05;
    _dwork[i] = r;
    r2 = r*r;
    r3 = r2*r;
    a = -z/r;
    b = e - a;
    xi = sqrt(1.0 + 0.5*FINE_STRUCTURE_CONST2*b);
    xi = xi*xi;
    b = b*b;
    v1 = z/r2;
    w = (-2.0*v1/r + 0.75*FINE_STRUCTURE_CONST2*v1*v1/xi - 2*ka*v1/r);
    w /= 4.0*xi;
    _dwork1[i] = a + 0.5*kl1/r2 - 0.5*FINE_STRUCTURE_CONST2*(b - w);
    a = z/r2;
    b = kl1/r3;
    if (a < dk && b < dk) break;
  }
  r0 = r;
  w = 2.0*(e - _dwork1[i]);
  y[0] = pow(w, -0.25);
  a = FINE_STRUCTURE_CONST2*e;
  b = FINE_STRUCTURE_CONST*z;
  y[1] = 0.5*(y[0]/w)*(z*(1+a)/r2-(kl1-b*b)/r3);

  rwork = _dwork2;
  lrw = pot->maxrp;
  liw = 22;
  neq = 2;
  itol = 2;
  rtol = EPS4;
  atol[0] = 0.0;
  atol[1] = EPS6;
  itask = 1;
  istate = 1;
  iopt = 0;
  mf = 10;

  i--;
  for (; i >= 0; i--) {
    r = _dwork[i];
    y[2] = r0;
    y[3] = _dwork1[i+1]*r0;
    y[4] = (_dwork1[i]*r - y[3])/(r - r0);
    y[5] = e;
    while (r0 != r) {
      LSODE(C_FUNCTION(DERIVODE, derivode), neq, y, &r0, r, itol, rtol, atol, 
	    itask, &istate, iopt, rwork, lrw, iwork, liw, NULL, mf);
      if (istate == -1) istate = 2;
      else if (istate < 0) {
	printf("LSODE Error %d\n", istate);
	exit(1);
      }
    }
  }

  p[n] = y[0];
  for (i = n-1; i >= i0; i--) {
    r = pot->rad[i];
    y[2] = r0;
    y[3] = _veff[i+1]*r0;
    y[4] = (_veff[i]*r - y[3])/(r - r0);
    y[5] = e;
    while (r0 != r) {
      LSODE(C_FUNCTION(DERIVODE, derivode), neq, y, &r0, r, itol, rtol, atol,
	    itask, &istate, iopt, rwork, lrw, iwork, liw, NULL, mf);
      if (istate == -1) istate = 2;
      else if (istate < 0) {
	printf("LSODE Error %d\n", istate);
	exit(1);
      }
    }
    p[i] = y[0];
    /*
    if (isnan(p[i])) {
      printf("nan amplitude: %d %g %g\n", i, y[0], y[1]);
      Abort(1);
    }
    */
  }

  return y[1];
}    

int Phase(double *p, POTENTIAL *pot, int i1, double phase0) {
  int i;
  double fact;
  
  fact = 1.0 / 3.0;

  for (i = i1; i < pot->maxrp; i++) {
    _dwork[i] = 1.0/p[i];
    _dwork[i] *= _dwork[i];
    _dwork[i] *= pot->dr_drho[i];
  }

  i = i1+1;
  p[i] = phase0;
  for (i = i1+3; i < pot->maxrp; i += 2) {
    p[i] = p[i-2] + (_dwork[i-3] + 4.0*_dwork[i-2] + _dwork[i-1])*fact;
    /*
    if (isnan(p[i])) {
      printf("nan phase: %d %g\n", i, p[i-2]);
    }
    */
  }
  return 0;
}

int SetVEffective(int kl, int kv, POTENTIAL *pot) {
  double kl1;
  int i;
  double r;

  kl1 = 0.5*kl*(kl+1);
 
  for (i = 0; i < pot->maxrp; i++) {
    r = pot->rad[i];
    r *= r;
    _veff[i] = pot->VT[kv][i] + kl1/r;
    _veff[i] += pot->W[i];
  }

  return 0;
}

static int TurningPoints(int n, double e, POTENTIAL *pot) {
  int i, i2;
  double x, a, b;

  if (n == 0) {
    for (i = 10; i < pot->maxrp-5; i++) {
      x = e - _veff[i];
      if (x <= 0) continue;
      b = 1.0/pot->rad[i];
      a = 20.0/(pot->qr*pot->ar*pow(b,1-pot->qr) + pot->br*b);
      x = TWO_PI/sqrt(2.0*x);
      if (x < a) break;
    }
    i2 = i-2;
    if (IsOdd(i2)) (i2)--;
  } else if (n > 0) {
    for (i = pot->maxrp-1; i > 10; i--) {
      if (e - _veff[i] > wave_zero) {
	i++;
	break;
      }
    }
    if (i <= 10) return -2;
    i2 = pot->maxrp-5;
    i2 = Min(i2, i);
    if (pot->ib && n > pot->nb) {
      if (i2 > pot->ib) {
	i2 = pot->ib;
      }
    }
  } else {
    for (i = pot->ib; i < pot->ib1; i++) {
      if (e - _veff[i] > wave_zero) {
	if (i > pot->ib) i--;
	break;
      }
    }
    i2 = i+20;
    for (i = pot->ib1-20; i > i2; i--) {
      if (e - _veff[i] > wave_zero) {
	i++;
	break;
      }
    }
    i2 = i;
  }
  return i2;
}

static int CountNodes(double *p, int i1, int i2) {
  int n, i;
  double a, p0;

  n = 0;
  i = i2;
  p0 = p[i];
  for (; i >= i1; i--) {
    a = fabs(p[i]);
    if (a > wave_zero) {
      if ((p0 > 0 && p[i] < 0) ||
	  (p0 < 0 && p[i] > 0)) {
	n++;
	p0 = p[i];
      }
    }
  }

  return n;
}

static int IntegrateRadial(double *p, double e, POTENTIAL *pot,
			   int i1, double p1, int i2, double p2, int q) {
  double a, b, r, x, y, z, p0, a1, a2;
  int kl=1, ku=1, nrhs=1;
  int i, info, n, m, j, k;
  int *ipiv = _ipiv;
  
  m = i2 - i1 - 1;
  if (m < 0) return 0;
  if (m == 0) {
    if (q == 0) {
      p[i1] = p1;
      p[i2] = p2;
    } else if (q == 1) {
      p[i1] = p1;
      p[i2] = (1.0+p2)*p1;
    } else if (q == 2) {
      p[i1] = (1.0-p1)*p2;
      p[i2] = p2;
    }
    return 0;
  }

  for (i = i1+1; i < i2; i++) {
    p[i] = 0.0;
  }

  j = 1;
  n = 2*kl + ku + 1;
  k = kl + ku;
  for (i = i1+1; i < i2; i++, j++, k += n) {
    y = pot->dr_drho[i]*pot->dr_drho[i];
    x = (2.0*(_veff[i] - e) + pot->vtr[i])*y/12.0;
    a = 1.0 - x;
    b = -2.0*(1.0 + 5.0*x);
    ABAND[k-1] = a;
    ABAND[k] = b;
    ABAND[k+1] = a;
  }

  i = i2;
  y = pot->dr_drho[i]*pot->dr_drho[i];
  x = (2.0*(_veff[i] - e) + pot->vtr[i])*y/12.0;
  a2 = 1.0 - x;
  if (q == 1) {
    k -= n;
    ABAND[k] += 2*p2*a2;
    k -= n;
    ABAND[k+1] += a2;
  } else {
    p[i2] = p2;
    p[i2-1] += -a2*p2;
  }
    
  i = i1;
  y = pot->dr_drho[i]*pot->dr_drho[i];
  x = (2.0*(_veff[i] - e) + pot->vtr[i])*y/12.0;
  a1 = 1.0 - x;
  if (q == 2) {
    k = kl + ku;
    ABAND[k] -= 2*p1*a1;
    k += n;
    ABAND[k-1] += a1;
  } else {
    p[i1] = p1;
    p[i1+1] += -a1*p1;
  }

  DGBSV(m, kl, ku, nrhs, ABAND, n, ipiv, p+i1+1, m, &info);
  if (info) {
    printf("Error in Integrating the radial equation: %d\n", info);
    exit(1);
  }

  if (q == 1) {
    p[i2] = 2*p2*p[i2-1] + p[i2-2];
  }
  if (q == 2) {
    p[i1] = -2*p1*p[i1+1] + p[i1+2];
  }

  n = 0;
  i = i2;
  p0 = p[i];
  for (; i >= i1; i--) {
    if (pot->rad[i] < pot->atom->rms0) break;
    a = fabs(p[i]);
    if (a > wave_zero) {
      if ((p0 > 0 && p[i] < 0) ||
	  (p0 < 0 && p[i] > 0)) {
	n++;
	p0 = p[i];
      }
    }
  }
    
  return n;
}

double InnerProduct(int i1, int n, double *p1, double *p2, POTENTIAL *pot) {
  int k;

  for (k = i1; k <= n; k++) {
    _dwork[k] = p1[k]*p2[k] * pot->dr_drho[k];
  }
  return Simpson(_dwork, i1, n);
}

void MaxRGrid(POTENTIAL *pot, double gasymp, double gratio, double z,
	      double rmin, int maxrp, double *ap, double *cp, double *rmp) {
  double a, c, d1, d2, rmax, r1;
  if (gasymp > 0 && gratio > 0) {
    a = gasymp*sqrt(2.0*z)/PI;
    c = 1.0/log(gratio);    
    d2 = maxrp-10.0 + a*pow(rmin, pot->qr) + c*log(rmin);
    rmax = d2/a;
    rmax = pow(rmax, 1.0/pot->qr);
    d1 = 1.0;
    while (d1 > EPS5) {
      r1 = d2 - c*log(rmax);
      r1 = r1/a;
      r1 = pow(r1, 1.0/pot->qr);
      d1 = fabs(r1/rmax-1.0);
      rmax = r1;
    }
  } else if (gratio > 0) {
    rmax = -gasymp;
    if (rmax > 1e10 && pot->rb > rmin) {
      rmax = pot->rb*1.001;
    }
    c = 1.0/log(gratio);
    a = maxrp-15.0 + c*(log(rmin)-log(rmax));
    a /= pow(rmax, pot->qr) - pow(rmin, pot->qr);
  } else if (gasymp > 0) {
    rmax = -gratio;
    if (rmax > 1e10 && pot->rb > rmin) {
      rmax = pot->rb*1.001;
    }
    a = gasymp*sqrt(2.0*z)/PI;
    c = maxrp-15.0 + a*(pow(rmin, pot->qr)-pow(rmax, pot->qr));
    c /= log(rmax) - log(rmin);
  }
  *ap = a;
  *cp = c;
  *rmp = rmax;
}

int SetOrbitalRGrid(POTENTIAL *pot) {
  int i, maxrp;  
  double z0, z, d1, d2, del, gratio, gasymp;
  double a, b, c, q, rn, r1, rmin, rmax;

  gratio = pot->ratio;
  gasymp = pot->asymp;
  z0 = GetAtomicNumber();
  rn = GetAtomicR();
  z = z0 - pot->N1;
  if (pot->flag == 0) pot->flag = -1; 
  rmin = pot->rmin/z0;
  if (rn > 0) {
    a = rn*GRIDRMINN0;
    if (rmin < a) rmin = a;
    a = rn*GRIDRMINN1;
    if (rmin > a) rmin = a;
  }

  pot->xps = 1.0;
  pot->gps = 0.0;
  if (pot->mps >= 0) {
    if (pot->zps <= 0) {
      pot->zps = (z0-pot->N)+_sp_z0 + (_sp_z1-_sp_z0)*(1-pot->N/z0);
    }
    if (pot->ups < 0) {
      pot->ups = pot->zps;
    }
    pot->dps = sqrt(pot->tps/(FOUR_PI*pot->nps*(1+pot->ups)));
    if (pot->nps > 0) {
      pot->rps = pow(3*pot->zps/(FOUR_PI*pot->nps),ONETHIRD);
    } else {
      pot->rps = pot->dps;
    }
    pot->gps = (pot->zps*pot->ups)/(pot->tps*pot->rps);
    if (pot->mps == 2) {
      pot->aps = FermiDegeneracy(pot->nps, pot->tps, &pot->fps);
      PrepFermiRM1(pot->aps, pot->fps, pot->tps);
    }    
  }  
  MaxRGrid(pot, gasymp, gratio, z, rmin, pot->maxrp, &a, &c, &rmax);
  maxrp = pot->maxrp;
  d1 = log(rmax/rmin);
  d2 = pow(rmax, pot->qr) - pow(rmin, pot->qr);
  b = (maxrp - 1.0 - (a*d2))/d1;
  i = (int)(1+a*d2+c*d1);
  if (i > pot->maxrp) {
    printf("Not enough radial mesh points: %d %g %g %g %g %g %g, ",
	   maxrp, gasymp, gratio, rmin, rmax, b, c);
    printf("enlarge to at least %d\n", i);
    exit(1);
  }
  if (pot->rps > 0) {
    rmax = pot->rps;
    d1 = log(rmax/rmin);
    d2 = pow(rmax, pot->qr) - pow(rmin, pot->qr);
    maxrp = (int)(1+a*d2 + b*d1);
    b = (maxrp - 1.0 - (a*d2))/d1;
  }
  d1 = b*d1;
  d2 = a*d2;
  del = (d1 + d2)/(maxrp - 1);
  if (fabs(del-1) > EPS8) {
    printf("inconsistent grid transformation: %g %g %g\n", a, b, del);
    Abort(1);
  }
  pot->rad[0] = rmin;
  d1 = a*pow(rmin, pot->qr) + b*log(rmin);
  pot->rho[0] = d1;
  for (i = 1; i < pot->maxrp; i++) {
    d1 += del;
    pot->rho[i] = d1;
    pot->rad[i] = GetRFromRho(d1, a, b, pot->qr, pot->rad[i-1]);
  }
  rmax = pot->rad[pot->maxrp-1];
  if (pot->rb > rmin) {
    pot->nmax = (int)(10+sqrt(rmax*z)/2.0);
  } else {
    pot->nmax = (int)(sqrt(rmax*z)/2.0);
  }
  pot->ar = a;
  pot->br = b;
  double tp2, tp3;
  q = pot->qr*(pot->qr-1);
  c = q*(pot->qr-2);
  for (i = 0; i < pot->maxrp; i++) {
    d1 = a * pow(pot->rad[i], pot->qr);
    d2 = pot->rad[i]/pot->qr;
    pot->dr_drho[i] = d2/(d1 + b/pot->qr);
    pot->dr_drho2[i] = sqrt(pot->dr_drho[i]);
    r1 = pot->rad[i]*pot->rad[i];
    tp2 = (q*d1-b)/r1;
    tp3 = (c*d1+2*b)/(r1*pot->rad[i]);
    pot->vtr[i] = pot->dr_drho[i]*(0.5*tp3-0.75*pot->dr_drho[i]*tp2*tp2);
  }

  RGMQED(&a, &b);
  for (i = 0; i < pot->maxrp; i++) {
    pot->mqrho[i] = b*log(pot->rad[i]) + a*pot->rad[i];
  }
  
  if (pot->rps > 0) {
    for (i = 0; i < pot->maxrp; i++) {
      if (pot->rad[i] > pot->rps) {
	pot->ips = i-1;
	break;
      }
    }
    if (pot->ips <= 0) {
      pot->rps = 0;
    }
  }
  return 0;
}

double GetRFromRho(double rho, double a, double b, double q, double r0) {
  double e, d1;
  int i;

  e = 1.0;
  i = 0;
  while (fabs(e) > 1E-10) {
    if (i > 100) {
      printf("Newton iteration failed to converge in GetRFromRho\n");
      exit(1);
    }
    d1 = pow(r0, q)*a;
    e = d1 + b*log(r0) - rho;
    e /= (q*d1 + b);
    r0 *= (1.0 - e);
    i++;
  }

  return r0;
}

int SetPotentialZ(POTENTIAL *pot) {
  int i;

  for (i = 0; i < pot->maxrp; i++) {
    pot->qdist[i] = GetAtomicChargeDist(pot->rad[i]);
    pot->Z[i] = GetAtomicEffectiveZ(pot->rad[i]);
  }
  for (i = 0; i < pot->atom->nep; i++) {
    SetPotentialExtraZ(pot, i);
  }
  if (pot->atom->nepr > 0) {
    SetPotentialExtraZ(pot, pot->atom->nep);
  }
  if (pot->atom->rn > 0 || pot->atom->nep > 0 || pot->atom->nepr > 0) {
    Differential(pot->Z, pot->dZ, 0, pot->maxrp-1, pot->dr_drho);
    Differential(pot->dZ, pot->dZ2, 0, pot->maxrp-1, pot->dr_drho);
  } else {
    for (i = 0; i < pot->maxrp; i++) {
      pot->dZ[i] = 0;
      pot->dZ2[i] = 0;
    }
  }  
  SetPotentialVP(pot);
  SetPotentialSE(pot);
  if (pot->sps == 0) SetPotentialPS(pot, NULL, NULL, 1);
  return 0;
}

int SetPotentialVc(POTENTIAL *pot) {
  int i;
  double n, r, r2, v, x, y, y2, a, b, v0;

  for (i = 0; i < pot->maxrp; i++) {
    r = pot->rad[i];
    a = pot->Z[i];
    b = pot->dZ[i];
    y = pot->dZ2[i];
    pot->Vc[i] = - (a / r);
    r2 = r*r;
    pot->dVc[i] = a/r2 - b/r;
    pot->dVc2[i] = 2.0*(-a/r + b)/r2 - y/r;
  }
  
  n = pot->N1;
  if (n > 0 && (pot->a > 0 || pot->lambda > 0)) {
    for (i = 0; i < pot->maxrp; i++) {
      r = pot->rad[i];
      r2 = r*r;
      v0 = n/r;
      x = 1.0 + r*pot->a;
      a = exp(-pot->lambda * r);
      v = v0 * (1.0 - a/x);
      pot->Vc[i] += v;
      b = (pot->lambda + pot->a/x);
      y = -v/r + (v0 - v)*b;
      pot->dVc[i] += y;
      y2 = -y/r + v/r2 - (v0/r+y)*b - (v0-v)*(pot->a*pot->a)/(x*x);      
      pot->dVc2[i] += y2;
    }
  }
  return 0;
}

int SetPotentialVT(POTENTIAL *pot) {
  int i, k, k1;
  double r, r2, y, a, b;
  for (i = 0; i < pot->maxrp; i++) {
    pot->VT[0][i] = pot->Vc[i] + pot->U[i];
    pot->dVT[0][i] = pot->dVc[i] + pot->dU[i];
    pot->dVT2[0][i] = pot->dVc2[i] + pot->dU2[i];
  }
  if (pot->pvp && pot->mvp) {
    for (i = 0; i < pot->maxrp; i++) {
      r = pot->rad[i];
      r2 = r*r;
      a = pot->ZVP[i];
      b = pot->dZVP[i];
      y = pot->dZVP2[i];
      pot->VT[0][i] -= a/r;
      pot->dVT[0][i] += a/r2 - b/r;
      pot->dVT2[0][i] += 2.0*(-a/r +b)/r2 - y/r;
    }
  }
  if (pot->zps > 0 && pot->nps > 0) {
    for (i = 0; i < pot->maxrp; i++) {
      r = pot->rad[i];
      r2 = r*r;
      a = pot->ZPS[i];
      b = pot->dZPS[i];
      y = pot->dZPS2[i];
      pot->VT[0][i] -= a/r;
      pot->dVT[0][i] += a/r2 - b/r;
      pot->dVT2[0][i] += 2.0*(-a/r +b)/r2 - y/r;
    }
  } 
  if (pot->pse && pot->nse) {
    for (k = 0; k < NKSEP; k++) {
      k1 = k+1;
      for (i = 0; i < pot->maxrp; i++) {
	r = pot->rad[i];
	r2 = r*r;
	a = pot->ZSE[k][i];
	b = pot->dZSE[k][i];
	y = pot->dZSE2[k][i];
	pot->VT[k1][i] = pot->VT[0][i] - a/r;
	pot->dVT[k1][i] = pot->dVT[0][i] + a/r2 - b/r;
	pot->dVT2[k1][i] = pot->dVT2[0][i] + 2.0*(-a/r +b)/r2 - y/r;
      }
    }
  } else {
    for (k = 0; k < NKSEP; k++) {
      k1 = k+1;
      for (i = 0; i < pot->maxrp; i++) {
	pot->VT[k1][i] = pot->VT[0][i];
	pot->dVT[k1][i] = pot->dVT[0][i];
	pot->dVT2[k1][i] = pot->dVT2[0][i];
      }
    }
  }
  return 0;
}

int SetPotentialU(POTENTIAL *pot, int n, double *u) {
  int i;
  double a;
  
  if (n < 0) {
    for (i = 0; i < pot->maxrp; i++) { 
      pot->U[i] = 0.0;
      pot->dU[i] = 0.0;
      pot->dU2[i] = 0.0;
    }
    if (pot->vxf > 0) {
      for (i = 0; i < pot->maxrp; i++) {
	pot->U[i] -= pot->VXF[i]/pot->rad[i];
      }
    }  
    return 0;
  }

  if (u) {
    for (i = 0; i < n; i++) {
      pot->U[i] = u[i];    
    }
  }

  if (pot->vxf > 0) {
    for (i = 0; i < pot->maxrp; i++) {
      pot->U[i] -= pot->VXF[i]/pot->rad[i];
    }
  }
  for (i = pot->maxrp-1; i >= 0; i--) {
    a = pot->rad[i]*pot->U[i];
    if (fabs(a) > EPS10) break;
    else pot->U[i] = 0.0;
  }
  
  Differential(pot->U, pot->dU, 0, pot->maxrp-1, pot->dr_drho);
  Differential(pot->dU, pot->dU2, 0, pot->maxrp-1, pot->dr_drho);
  return 0;
}

int IdxVT(int kappa) {
  int k;
  k = 2*(abs(kappa)-1)+(kappa>0);
  if (k >= NKSEP) k = -1;
  return k+1;
}

int SetPotentialW(POTENTIAL *pot, double e, int kappa, int k) {
  int i;
  double xi, r, r2, x, y, z;

  for (i = 0; i < pot->maxrp; i++) {
    xi = e - pot->VT[k][i];
    r = xi*FINE_STRUCTURE_CONST2*0.5 + 1.0;  
    x = pot->dVT[k][i];
    y = - 2.0*kappa*x/pot->rad[i];
    x = x*x*0.75*FINE_STRUCTURE_CONST2/r;
    z = pot->dVT2[k][i];
    pot->W[i] = x + y + z;
    pot->W[i] /= 4.0*r;
    x = xi*xi;
    pot->W[i] = x - pot->W[i];
    pot->W[i] *= 0.5*FINE_STRUCTURE_CONST2;
    pot->W[i] = -pot->W[i];
  }

  Differential(pot->W, pot->dW, 0, pot->maxrp-1, pot->dr_drho);
  Differential(pot->dW, pot->dW2, 0, pot->maxrp-1, pot->dr_drho);
  return 0;
}

static double Poly(double x, double n, double c[]) {
  int i;
  double t, r;

  t = 1.0;
  r = 0.0;
  for (i = 0; i < n; i++) {
    r += c[i]*t;
    t *= x;
  }
  
  return r;
}
  
static double UehlingK0(double x) {
  double a[8] = {+0.88357293375,
		 -0.28259817381,
		 -0.58904879578,
		 +0.12500133434,
		 -0.032729913852,
		 +0.0082888574511,
		 -0.0000103277658,
		 +0.0000636436689};
  double b[3] = {-319.999594328,
		 +2.53900995981,
		 1.0};
  double c[3] = {-319.999594333,
		 +2.53901020662,
		 0.0};
  double d[5] = {5.018065179,
		 71.51891262,
		 211.6209929,
		 31.40327478,
		 -1.0};
  double e[5] = {2.669207401,
		 51.72549669,
		 296.9809720,
		 536.4324164,
		 153.5335924};

  double r, x2;

  if (x <= 1.0) {
    r = Poly(x, 8, a);
    if (x) {
      x2 = x*x;
      r += x*log(x)*Poly(x2, 3, b)/Poly(x2, 3, c);
    }
  } else {
    r = exp(-x)/pow(x, 1.5);
    x2 = 1.0/x;
    r *= Poly(x2, 5, d)/Poly(x2, 5, e);
  }

  return r;
}

static double UehlingK1(double x) {
  double a[8] = {-0.71740181754,
		 1.1780972274,
		 -0.37499963087,
		 0.13089675530,
		 -0.038258286439,
		 -0.0000242972873,
		 -0.0003592014867,
		 -0.0000171700907};
  double b[3] = {-64.0514843293,
		 0.711722714285,
		 1.0};
  double c[3] = {64.0514843287,
		 -0.711722686403,
		 0.0008042207748};
  double d[5] = {217.2386409,
		 1643.364528,
		 2122.244512,
		 1.0};
  double e[5] = {115.5589983,
		 1292.191441,
		 3831.198012,
		 2904.410075,
		 0.0};

  double r, x2;

  if (x <= 1.0) {
    r = Poly(x, 8, a);
    x2 = x*x;
    r += log(x)*Poly(x2, 3, b)/Poly(x2, 3, c);
  } else {
    r = exp(-x)/pow(x, 1.5);
    x2 = 1.0/x;
    r *= Poly(x2, 5, d)/Poly(x2, 5, e);
  }

  return r;
}

double UehlingL0(double x) {
  double f[6] = {1.990159,
		 -2.397605,
		 1.046471,
		 -0.367066,
		 0.063740,
		 -0.037058};
  double g[3] = {0.751198,
		 0.128889,
		 0.020886};
  double h[2] = {-0.444444,
		 -0.003472};
  
  double logx, x2, r;
  
  if (x <= 2) {
    logx = log(x);
    r = Poly(x, 6, f);
    if (x) {
      x2 = x*x;
      r += x*Poly(x2, 3, g)*logx;
      x2 = x2*x2;
      r += x*Poly(x2, 2, h)*logx*logx;
    }
  } else {
    r = 0.0;
  }
  
  return r;
}

double UehlingL1(double x) {
  double f[6] = {1.646407,
		 -2.092942,
		 0.962310,
		 -0.254960,
		 0.164404,
		 0.0};
  double g[3] = {0.137691,
		 -0.416667,
		 -0.097486};
  double h[2] = {0.444444,
		 0.017361};
  
  double logx, x2, r;
  
  if (x <= 2) {
    logx = log(x);
    r = Poly(x, 6, f);
    x2 = x*x;
    r += Poly(x2, 3, g)*logx;
    x2 = x2*x2;
    r += Poly(x2, 2, h)*logx*logx;
  } else {
    r = 0.0;
  }
  
  return r;
}

int SetPotentialSE(POTENTIAL *pot) {
  int i, j;
  for (i = 0; i < NKSEP; i++) {
    LOCSEP(i, pot->maxrp, pot->mqrho, pot->ZSE[i]);
    Differential(pot->ZSE[i], pot->dZSE[i], 0, pot->maxrp-1, pot->dr_drho);
    Differential(pot->dZSE[i], pot->dZSE2[i], 0, pot->maxrp-1, pot->dr_drho);
    for (j = 0; j < pot->maxrp; j++) {
      pot->ZSE[i][j] = -pot->ZSE[i][j];
      pot->dZSE[i][j] = -pot->dZSE[i][j];
      pot->dZSE2[i][j] = -pot->dZSE2[i][j];
    }
    RADFND(i, &pot->rfn[i]);
    pot->nfn[i] = 0;
  }
  return 0;
}

int SetPotentialVP(POTENTIAL *pot) {
  int i, j, k, p, m, n, mvp, vp;
  double a, b, z0, r0, r;
  double v3, za, za2, za3, a2, lam, r2, r3, r5, x, y, f0, fs, fm;

  mvp = pot->mvp%100;
  vp = mvp%10;
  mvp = mvp/10;
  r0 = 3.86159E-3/RBOHR;
  z0 = pot->atom->atomic_number;
  a = -2.0*z0*FINE_STRUCTURE_CONST/(3.0*PI);
  b = -z0*FINE_STRUCTURE_CONST2/(PI*PI);  
  for (i = 0; i < pot->maxrp; i++) {
    r = pot->rad[i]*2.0/r0;
    pot->ZVP[i] = -a*UehlingK1(r);
  }
  if (vp > 1) {
    for (i = 0; i < pot->maxrp; i++) {
      pot->ZVP[i] -= b*UehlingL1(r);
    }
  }
  if (vp > 2) {
    za = z0*FINE_STRUCTURE_CONST;
    za2 = za*za;
    za3 = za2*za;
    a = FINE_STRUCTURE_CONST;
    a2 = FINE_STRUCTURE_CONST2;
    lam = 1-sqrt(1-za2);
    for (i = 0; i < pot->maxrp; i++) {
      r = pot->rad[i]/a;
      v3 = 0.0;
      r2 = r*r;
      r3 = r2*r;
      r5 = r2*r3;
      if (r >= 10) {
	y = 1.0/r5;
	fm = 1.0;
	for (m = 0; m < 16; m++) {
	  x = 1.0;
	  f0 = 0.0;
	  for (n = 0; n < 5; n++) {
	    f0 += _wk_f[m][n]*x;
	    x *= za2;
	  }
	  if (m < 5) {
	    f0 /= _wk_fd[m];
	  } else {
	    f0 *= fm*fm;
	  }
	  v3 += f0*y;
	  y /= r2;
	  fm *= m+1.0;
	}
	v3 /= PI;
      } else {
	if (r < 4.25) {
	  x = 1.0/r;
	  y = log(r);
	  for (n = 0; n < 9; n++) {
	    v3 += _wk_p[n]*x;
	    x *= r;
	  }
	  x = r*r*y;
	  for (n = 0; n < 6; n++) {
	    v3 += _wk_p[n+9]*x;
	    x *= r;
	  }
	  x = r3*y*y;
	  for (n = 0; n < 5; n++) {
	    v3 += _wk_p[n+15]*x;
	    x *= r;
	  }
	} else {
	  v3 = 0.0;
	  y = 100/r2;
	  x = 1.0;
	  for (n = 0; n < 16; n++) {
	    v3 += _wk_p[n+20]*x;
	    x *= y;
	  }
	  v3 /= r5;
	}
	if (r < 0.15) {
	  k = 0;
	  x = lam/r;
	  for (n = 1; n < 16; n++) {
	    v3 += 1e-2*_wk_q[n]*x;
	    x *= lam;
	  }
	} else if (r >= 0.15 && r < 1.5) {
	  k = 1;
	} else if (r >= 1.5 && r < 4) {
	  k = 2;
	} else if (r >= 4) {
	  k = 3;
	}
	y = lam;
	fs = 0;
	for (n = 0; n < 4; n++) {
	  if (_wk_s[k][n] == 10) continue;
	  x = pow(r, _wk_s[k][n]);
	  f0 = 0;
	  for (m = 0; m < 8; m++) {
	    f0 += _wk_a[k][n][m]*x;
	    x *= r;
	  }
	  fs += f0*y;
	  y *= lam;
	}
	v3 /= (1+fs);
      }
      pot->ZVP[i] -= (za3*v3/a)*pot->rad[i];
    }
  }
  
  m = pot->maxrp;
  if (pot->atom->rn > 0 && mvp == 0) {
    r5 = pot->atom->rn*5.0;
    r3 = r5*3;
    for (i = 0; i < m; i++) {
      _dwork[i] = pot->ZVP[i]*pot->dr_drho[i];
      _dwork1[i] = pot->rho[i];
    }
    _dwork2[m-1] = 0;
    NewtonCotesIP(_dwork2, _dwork, 0, m-1, -1, -1);
    n = 3;
    for (k = 0; k < m; k++) {
      if (pot->qdist[k] < 1e-15*pot->qdist[0]) break;
    }
    for (i = 0; i < m; i++) {
      for (j = 0; j < k; j++) {
	r = fabs(pot->rad[i] - pot->rad[j]);
	if (r < pot->rad[0]) r = pot->rad[0];
	pot->dZVP[j] = pot->ar*pow(r, pot->qr) + pot->br*log(r);
      }
      UVIP3P(n, m, _dwork1, _dwork2, k, pot->dZVP, _dwork3);
      for (j = 0; j < k; j++) {	
	r = pot->rad[i] + pot->rad[j];
	pot->dZVP[j] = pot->ar*pow(r, pot->qr) + pot->br*log(r);
      }
      UVIP3P(n, m, _dwork1, _dwork2, k, pot->dZVP, pot->dZVP2);
      for (j = 0; j < k; j++) {
	_dwork3[j] -= pot->dZVP2[j];
	_dwork3[j] *= pot->qdist[j]*pot->rad[j];
	_dwork3[j] *= pot->dr_drho[j];
      }
      pot->dW2[i] = Simpson(_dwork3, 0, k-1)/(2*pot->atom->atomic_number);
      if (pot->rad[i] > r5 &&
	  fabs(pot->dW2[i]-pot->ZVP[i]) <= 1e-5*fabs(pot->ZVP[i])) {
	break;
      }
      if (pot->rad[i] > r3) {
	break;
      }
    }
    p = i;
    for (i = 0; i < p; i++) {
      pot->ZVP[i] = pot->dW2[i];
    }
  }
  Differential(pot->ZVP, pot->dZVP, 0, m-1, pot->dr_drho);
  Differential(pot->dZVP, pot->dZVP2, 0, m-1, pot->dr_drho);
  return 0;
}

int SetPotentialExtraZ(POTENTIAL *pot, int iep) {
  int i, m, j, k, n, p;
  double r, r5, r3;
  
  m = pot->maxrp;
  for (i = 0; i < m; i++) {
    r = pot->rad[i];
    pot->dW[i] = GetExtraZ(r, iep);
  }
  int conv = 0;
  if (iep >= pot->atom->nep) {
    if (pot->atom->nepr > 0 && pot->atom->cepr > 0) {
      conv = 1;
    }
  } else if (pot->atom->epm[iep] >= 100 && pot->atom->rn > 0) {
    conv = 1;
  }
  if (conv) {
    r5 = pot->atom->rn*5.0;
    r3 = r5*3;
    for (i = 0; i < m; i++) {
      _dwork[i] = pot->dW[i]*pot->dr_drho[i];
      _dwork1[i] = pot->rho[i];
    }
    _dwork2[m-1] = 0;
    NewtonCotesIP(_dwork2, _dwork, 0, m-1, -1, -1);
    n = 3;
    for (k = 0; k < m; k++) {
      if (pot->qdist[k] < 1e-15*pot->qdist[0]) break;
    }
    for (i = 0; i < m; i++) {
      for (j = 0; j < k; j++) {
	r = fabs(pot->rad[i] - pot->rad[j]);
	if (r < pot->rad[0]) r = pot->rad[0];
	pot->dZ[j] = pot->ar*pow(r, pot->qr) + pot->br*log(r);
      }
      UVIP3P(n, m, _dwork1, _dwork2, k, pot->dZ, _dwork3);
      for (j = 0; j < k; j++) {	
	r = pot->rad[i] + pot->rad[j];
	pot->dZ[j] = pot->ar*pow(r, pot->qr) + pot->br*log(r);
      }
      UVIP3P(n, m, _dwork1, _dwork2, k, pot->dZ, pot->dZ2);
      for (j = 0; j < k; j++) {
	_dwork3[j] -= pot->dZ2[j];
	_dwork3[j] *= pot->qdist[j]*pot->rad[j];
	_dwork3[j] *= pot->dr_drho[j];
      }
      pot->dW2[i] = Simpson(_dwork3, 0, k-1)/(2*pot->atom->atomic_number);
      if (pot->rad[i] > r5 &&
	  fabs(pot->dW2[i]-pot->dW[i]) <= 1e-5*fabs(pot->dW[i])) {
	break;
      }
      if (pot->rad[i] > r3) {
	break;
      }
    }
    p = i;   
    for (i = 0; i < p; i++) {
      pot->dW[i] = pot->dW2[i];
    }
  }
  for (i = 0; i < m; i++) {
    pot->Z[i] += pot->dW[i];
  }
  return 0;
}

int SetPotentialPS(POTENTIAL *pot, double *vt, double *wb, int init) {
  int i;
  
  if (pot->mps < 0 || init) {
    for (i = 0; i < pot->maxrp; i++) {
      pot->NPS[i] = 0;
      pot->EPS[i] = 0;
      pot->ZPS[i] = 0;
      pot->dZPS[i] = 0;
      pot->dZPS2[i] = 0;
    }
    return 0;
  }
  SetPotentialIPS(pot, vt, wb);
  return 0;
}

void SetPotentialIPS(POTENTIAL *pot, double *vt, double *wb) {
  int i;
  double n0, r0, x;
  if (pot->mps == 1) {//debye screening
    if (vt) return;
    for (i = 0; i < pot->maxrp; i++) {
      pot->NPS[i] = 0;
      pot->EPS[i] = 0;
      r0 = pot->rad[i]/pot->dps;
      x = exp(-r0);
      pot->ZPS[i] = pot->zps*(x-1.0);
      pot->dZPS[i] = -pot->zps*x/pot->dps;
      pot->dZPS2[i] = pot->zps*x/(pot->dps*pot->dps);
    }
    return;
  }

  if (vt && pot->vxf > 0) {
    for (i = 0; i < pot->maxrp; i++) {
      x = pot->VXF[i];
      vt[i] += x/pot->rad[i];      
    }
  }

  if (pot->mps == 2) {//stewart&pyatt
    if (vt) {
      if (_sp_mode == 0 ||
	  _sp_mode == 2 ||
	  _sp_mode == 4) return;
      if (_sp_mode != 3) {
	for (i = 0; i < pot->maxrp; i++) {
	  x = pot->rad[i]/pot->dps;
	  ABAND[i] = (x*fabs(vt[i]))/pot->tps;
	}
      }
    }
    int nicf = _icf_nfft;
    _icf_nfft = 0;
    i = 0;
    double x0 = 0.75;
    double xps = 2.0;
    while (xps > 1) {
      xps = StewartPyatt(pot, vt, wb, x0);
      x0 *= 0.5;
      i++;
      if (i > 256) {
	printf("maxiter StewartPyatt0: %g %g\n", x0, xps);
	Abort(1);
      }
    }
    double x1 = 1.5;
    xps = 0.5;
    i = 0;
    while (xps < 1) {
      xps = StewartPyatt(pot, vt, wb, x1);
      x1 *= 2.0;
      i++;
      if (i > 256) {
	printf("maxiter StewartPyatt1: %g %g\n", x1, xps);
	Abort(1);
      }
    }
    i = 0;
    while (x1-x0 > EPS5) {
      x = 0.5*(x0+x1);
      xps = StewartPyatt(pot, vt, wb, x);
      if (xps < 1) {
	x0 = x;
      } else if (xps > 1) {
	x1 = x;
      } else {
	break;
      }
    }
    x = 0.5*(x0+x1);
    pot->xps = x;
    _icf_nfft = nicf;
    double *ve0 = ABAND+pot->maxrp*2;
    double *ve1 = ABAND+pot->maxrp*3;
    for (i = 0; i < pot->maxrp; i++) {
      ve0[i] = _veff[i];
      ve1[i] = 0.0;
    }
    if (_icf_nfft > 0) {
      int iter = 0;
      while(1) {
	IonCF(pot, _veff, _icf_rmax, _icf_nfft, _icf_dw);	
	xps = StewartPyatt(pot, vt, wb, x);
	if (fabs(xps-1) < EPS3) break;
	x0 = x;
	x1 = x;
	i = 0;
	if (xps < 1) {
	  while (xps < 1) {
	    xps = StewartPyatt(pot, vt, wb, x1);
	    x1 *= 2.0;
	    i++;
	    if (i > 256) {
	      printf("maxiter StewartPyatt1: %g %g\n", x1, xps);
	      Abort(1);
	    }
	  }
	} else {
	  while (xps > 1) {
	    xps = StewartPyatt(pot, vt, wb, x0);
	    x0 *= 0.5;
	    i++;
	    if (i > 256) {
	      printf("maxiter StewartPyatt0: %g %g\n", x0, xps);
	      Abort(1);
	    }
	  }
	}
	while (x1-x0 > EPS5) {
	  x = 0.5*(x0+x1);
	  xps = StewartPyatt(pot, vt, wb, x);
	  if (xps < 1) {
	    x0 = x;
	  } else if (xps > 1) {
	    x1 = x;
	  } else {
	    break;
	  }
	}
	x = 0.5*(x0+x1);
	iter++;	
	if (iter > _icf_spmi) {
	  //printf("maxiter ioncf: %d %g\n", iter, xps);
	  break;
	}
      }
      pot->xps = x;	
      double *rw = _icf_dw + 3*_icf_nfft;
      double *yw = _icf_dw + _icf_nfft;
      UVIP3P(3, _icf_nfft, rw, yw, pot->maxrp, pot->rho, pot->ICF);
      yw = _icf_dw;
      UVIP3P(3, _icf_nfft, rw, yw, pot->maxrp, pot->rho, ve1);
      double dr = pot->rps*_icf_rmax/_icf_nfft;
      double dk = PI/(_icf_nfft*dr);
      if (_icf_ofn[0]) {
	FILE *f = fopen(_icf_ofn, "w");
	if (f != NULL) {
	  for (i = 0; i < _icf_nfft; i++) {
	    fprintf(f, "%5d %15.8E %15.8E %12.5E %12.5E %12.5E\n",
		    i, i*dr, i*dk, _icf_dw[i],
		    _icf_dw[_icf_nfft+i],
		    _icf_dw[_icf_nfft*2+i]);
	  }
	  fclose(f);
	} else {
	  printf("cannot open file %s\n", _icf_ofn);
	}
      }
      yw = _icf_dw + _icf_nfft*4;
      double nz = pot->nps/pot->ups;
      double kf = pow(3*PI*PI*nz, ONETHIRD);
      double ku = (1.0/pot->rps)/kf;
      ku *= ku;
      double ga = 3*pot->gps;
      for (i = 1; i < _icf_nfft; i++) {
	x0 = i*dk/kf;
	x1 = x0*x0;
	yw[i] = (_icf_dw[_icf_nfft*2+i])/(x1+ku);
      }
      double sf1 = Simpson(yw, 1, _icf_nfft-1)*(dk/kf);
      double sf0 = ga/sqrt(3.68317+ga);
      pot->qps = (pot->ups)*2*sf0*sf1/(PI*pot->rps);
      pot->sf0 = sf0;
      pot->sf1 = sf1;
    }
    if (_sp_ofn[0]) {
      FILE *f = fopen(_sp_ofn, "w");
      if (f != NULL) {
	POTENTIAL *potential = pot;
	fprintf(f, "#    mps = %d\n", potential->mps);
	fprintf(f, "#    vxf = %d\n", potential->vxf);
	fprintf(f, "#    zps = %12.5E\n", potential->zps);
	fprintf(f, "#    nps = %12.5E\n", potential->nps);
	fprintf(f, "#    tps = %12.5E\n", potential->tps);
	fprintf(f, "#    ups = %12.5E\n", potential->ups);
	fprintf(f, "#    rps = %12.5E\n", potential->rps);
	fprintf(f, "#    dps = %12.5E\n", potential->dps);
	fprintf(f, "#    aps = %12.5E\n", potential->aps);
	fprintf(f, "#    fps = %12.5E\n", potential->fps);
	fprintf(f, "#    xps = %12.5E\n", potential->xps);
	fprintf(f, "#    jps = %12.5E\n", potential->jps);
	fprintf(f, "#    gps = %12.5E\n", potential->gps);
	fprintf(f, "#    qps = %12.5E\n", potential->qps);
	fprintf(f, "#    sf0 = %12.5E\n", potential->sf0);
	fprintf(f, "#    sf1 = %12.5E\n", potential->sf1);
	fprintf(f, "#    ips = %d\n", potential->ips);
	fprintf(f, "#    sps = %d\n", potential->sps);
	fprintf(f, "#   nmax = %d\n", potential->nmax);
	fprintf(f, "#  maxrp = %d\n", potential->maxrp);
	fprintf(f, "\n\n");
	for (i = 0; i < potential->maxrp; i++) {
	  fprintf(f, "%5d %14.8E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E\n",
		  i, potential->rad[i], potential->Z[i],
		  potential->Z[i]-GetAtomicEffectiveZ(potential->rad[i]),
		  potential->Vc[i]*potential->rad[i],
		  potential->U[i]*potential->rad[i],
		  ve0[i], _veff[i], ve1[i], wb?wb[i]:0.0,
		  potential->NPS[i], potential->EPS[i],
		  potential->VXF[i], potential->ICF[i],
		  potential->ZPS[i]);
	}    
	fclose(f);    
      } else {
	printf("cannot open file %s\n", _sp_ofn);
      }
    }
    Differential(pot->ZPS, pot->dZPS, 0, pot->maxrp-1, pot->dr_drho);
    Differential(pot->dZPS, pot->dZPS2, 0, pot->maxrp-1, pot->dr_drho);
    return;
  }
  
  if (pot->ips == 0) {
    for (i = 0; i < pot->maxrp; i++) {
      pot->NPS[i] = 0;
      pot->EPS[i] = 0;
      pot->ZPS[i] = 0;
      pot->dZPS[i] = 0;
      pot->dZPS2[i] = 0;
    }
    return;
  }
  n0 = pot->nps;
  r0 = pot->rps;
  if (pot->tps <= 0 || vt == NULL) {
    pot->ups = 0.0;
    for (i = 0; i < pot->maxrp; i++) {
      if (i <= pot->ips) {
	x = pot->rad[i]/r0;
	pot->NPS[i] = n0;
	pot->EPS[i] = n0;
	pot->ZPS[i] = -0.5*x*pot->zps*(3-x*x);
	pot->dZPS[i] = 1.5*pot->zps*(x*x-1.0)/r0;
	pot->dZPS2[i] = 3.0*pot->zps*x/(r0*r0);
      } else {
	pot->NPS[i] = 0.0;
	pot->EPS[i] = 0.0;
	pot->ZPS[i] = -pot->zps;
	pot->dZPS[i] = 0.0;
	pot->dZPS2[i] = 0.0;
      }
    }
    return;
  }
  for (i = 0; i < 128; i++) {
    n0 = pot->ups;
    FreeElectronDensity(pot, vt);
    if (fabs(1-exp(pot->ups-n0)) < _fermi_abserr) break;
    x = _fermi_stablizer;
    pot->ups = x*n0 + (1-x)*pot->ups;
  }
  for (i = 0; i <= pot->ips; i++) {
    x = pot->rad[i];
    _dwork[i] = x*x*pot->NPS[i]*pot->dr_drho[i];
  }
  _dwork1[0] = (pow(pot->rad[0],3)/3.0)*pot->NPS[0];
  NewtonCotesIP(_dwork1, _dwork, 0, pot->ips, -1, 0);
  for (i = 0; i <= pot->ips; i++) {
    x = pot->rad[i];
    _dwork[i] = x*pot->NPS[i]*pot->dr_drho[i];
  }
  _dwork2[pot->ips] = 0.0;
  NewtonCotesIP(_dwork2, _dwork, 0, pot->ips, -1, -1);
  for (i = 0; i <= pot->ips; i++) {
    pot->ZPS[i] = -FOUR_PI*(_dwork1[i] + pot->rad[i]*_dwork2[i]);
  }
  for (i = pot->ips+1; i < pot->maxrp; i++) {
    pot->ZPS[i] = pot->ZPS[pot->ips];
  }
  pot->jps = FOUR_PI*_dwork2[0];
  Differential(pot->ZPS, pot->dZPS, 0, pot->maxrp-1, pot->dr_drho);
  Differential(pot->dZPS, pot->dZPS2, 0, pot->maxrp-1, pot->dr_drho);
}

void FreeElectronDensity(POTENTIAL *pot, double *vt) {
  int i, k;
  double a = 1.0/(PI*PI);
  if (_relativistic_fermi == 2) {
    double c2 = 1.0/FINE_STRUCTURE_CONST2;
    double a2 = FINE_STRUCTURE_CONST2;
    int nk = Min(pot->maxrp, 1001);
    double v, k0, k1, k2, dk, y;
    for (i = 0; i <= pot->ips; i++) {
      v = vt[i];
      k0 = (-2*v+a2*v*v);
      k1 = 50.0*pot->tps+(c2-v);
      k1 = k1*k1*a2 - c2;
      if (k1 < k0) k1 = 2*k0;
      dk = (k1 - k0)/(nk-1);
      _dwork[0] = k0;
      for (k = 1; k < nk; k++) {
	_dwork[k] = _dwork[k-1] + dk;
      }
      for (k = 0; k < nk; k++) {
	k2 = _dwork[k];
	y = (sqrt(k2*c2 + c2*c2)-c2+v)/pot->tps;
	y = exp(pot->ups-y);
	_dwork1[k] = sqrt(k2)*y/(y + 1);
      }
      pot->NPS[i] = a*0.5*dk*Simpson(_dwork1, 0, nk-1);
      pot->EPS[i] = pot->NPS[i];
      _dwork2[i] = pot->NPS[i]*pot->rad[i]*pot->rad[i]*FOUR_PI*pot->dr_drho[i];
    }
  } else {
    double g, y;
    a *= sqrt(2*pot->tps)*pot->tps;
    if (_relativistic_fermi) g = pot->tps*FINE_STRUCTURE_CONST2;
    else g = 0;
    for (i = 0; i <= pot->ips; i++) {
      y = (-vt[i])/pot->tps;
      pot->NPS[i] = a*FermiIntegral(y+pot->ups, y, g);
      pot->EPS[i] = pot->NPS[i];
      _dwork2[i] = pot->NPS[i]*pot->rad[i]*pot->rad[i]*FOUR_PI*pot->dr_drho[i];
    }
  }
  _dwork1[0] = pot->NPS[0]*FOUR_PI*pow(pot->rad[0],3)/3.0;
  NewtonCotesIP(_dwork1, _dwork2, 0, pot->ips, -1, 0);
  a = pot->zps/_dwork1[pot->ips];
  pot->ups += log(a);
  for (i = 0; i <= pot->ips; i++) {
    pot->NPS[i] *= a;
    pot->EPS[i] = pot->NPS[i];
  }
}

static double _fermi_ax = 0.0;
static double _fermi_ag = 0.0;
static double FermiIntegrand(double *tp) {
  double et, x, g, t, gt, r;
  t = *tp;
  x = _fermi_ax;
  g = _fermi_ag;
  if (g <= 0) {
    gt = t;
  } else {
    gt = g*t;
    if (gt < 1e-6) {
      gt = t - 0.5*gt*t;
    } else {    
      gt = (sqrt(1+2*gt)-1.0)/g;
    }
  }
  if (gt > x) {
    et = exp(x-gt);
    r = sqrt(t)*(et/(1.0+et));
  } else {
    et = exp(gt-x);
    r = sqrt(t)/(1.0+et);
  }
  return r;
}
FCALLSCFUN1(DOUBLE, FermiIntegrand, FERMIINTEGRAND, fermiintegrand, PDOUBLE)
double FermiIntegral(double x, double y, double g) {
  double t0, t1;
  double r, s, a;
  int n, ier, lim, lenw, last, iter;
  int iwork[16];
  double dwork[16];

  lim = 4;
  lenw = 16;  
  r = 0.0;
  _fermi_ax = x;
  _fermi_ag = g;
  if (y > 0 && g > 0) y += 0.5*g*y*y;
  t0 = y;
  iter = 0;
  while (1) {
    double dt = Max(0.05*t0, 1.0);
    t1 = t0 + dt;
    DQAGS(C_FUNCTION(FERMIINTEGRAND, fermiintegrand),
	  t0, t1, _fermi_abserr, _fermi_relerr, &s, &a, &n, &ier,
	  lim, lenw, &last, iwork, dwork);
    r += s;
    t0 = t1;
    if (iter > 10 && s < 1e-3*r) {
      if (g <= 0) {
	a = t0;
      } else {
	a = g*t0;
	if (a < 1e-6) {
	  a = t0 - 0.5*a*t0;
	} else {
	  a = (sqrt(1+2*a)-1.0)/g;
	}
      }
      if (a-x > 35 && a > 10.0) {
	s = exp(x-a)*sqrt(a);
	if (s < _fermi_relerr*r) {
	  break;
	}
      }
    }
    iter++;
  }
  return r;
}

double FermiDegeneracy(double ne, double te, double *yi) {
  double a1, y1, g, d;
  int i;
  if (_relativistic_fermi) {
    g = te*FINE_STRUCTURE_CONST2;
  } else {
    g = 0.0;
  }
  double x = ne/(0.143289792*te*sqrt(te));
  double a0 = log(1.128379167*x);
  double y0 = FermiIntegral(a0, 0.0, g);
  if (y0 > x) {
    a1 = a0;
    i = 0;
    while(y0 > x) {
      d = 0.5*fabs(a0);
      d = Max(d, 0.5);
      a0 -= d;
      y0 = FermiIntegral(a0, 0.0, g);
      i++;
      if (i > 256) {
	printf("FermiDegeneracy maxiter0: %d %g %g %g\n",
	       i, a0, y0, x);
	Abort(1);
      }
    }
  } else {
    a1 = a0;
    y1 = y0;
    i = 0;
    while (y1 < x) {
      d = 0.5*fabs(a1);
      d = Max(d, 0.5);
      a1 += d;
      y1 = FermiIntegral(a1, 0.0, g);
      i++;
      if (i > 256) {
	printf("FermiDegeneracy maxiter1: %d %g %g %g\n",
	       i, a1, y1, x);
	Abort(1);
      }
    }    
  }
  double a, y;
  i = 0;
  while (a1-a0 > EPS6*fabs(a0+a1)) {
    a = 0.5*(a0+a1);
    y = FermiIntegral(a, 0.0, g);
    if (y < x) {
      a0 = a;
    } else if (y > x) {
      a1 = a;
    } else {
      break;
    }
    i++;
    if (i > 256) {
      printf("FermiDegeneracy maxiter2: %d %g %g %g %g\n",
	     i, a0, a1, y, x);
      Abort(1);
    }
  }
  a = 0.5*(a0+a1);
  *yi = FermiIntegral(a, 0.0, g);
  return a;
}

double ExpM1(double x) {
  if (fabs(x) < 1e-5) {
    return x+0.5*x*x;
  } 
  return exp(x)-1;
}

double FermiAsymRM1(double a, double y, int m) {
  const double a0 = 0.3535533905933;
  const double a1 = 0.7522527780637;
  const double a2 = 1.2337005501362;
  const double a3 = 1.1283791670955;
  double ay = fabs(y);
  double r;
  if (ay < 0.1) {
    if (m == 0) {
      r = ExpM1(y)+a0*exp(2*y+a)*ExpM1(-y);
    } else {
      r = y + a1*y*sqrt(y) + 0.5*y*y;
    }
    return r;
  } else {
    if (y < 0) return -1.0;
    if (m == 0) {
      double ya = y+a;
      return a1*exp(-a)*ya*sqrt(ya)*(1+a2/(ya*ya))-1.0;
    } else {
      return a3*sqrt(y)*(1+0.5/y-0.25/(y*y))-1.0;
    }
  }
}

void FermiRM1(int n, double *y, double *r, double *rn,
	      double a, double fa, double g, int m) {
  int i, im;
  double y0;
  im = -1;
  for (i = 0; i < n; i++) {
    y0 = exp(y[i]);
    if (m == 0) {
      r[i] = FermiIntegral(y0+a, 0.0, g)/fa-1.0;
      if (im >= 0) {
	rn[i] = 0.0;
      } else {
	rn[i] = FermiIntegral(-y0+a, 0.0, g)/fa;
      }
      if (im < 0 && rn[i] < 1e-8) {
	im = i;
      }
      rn[i] -= 1;
    } else {
      r[i] = FermiIntegral(y0+a, y0, g)/fa-1.0;
    }
  }
  double y0a=1.0, y0b=1.0;
  y0 = exp(y[0]);
  if (m > 0) {
    y0a = FermiAsymRM1(a, y0, 0);
  }
  y0 += 0.5*g*y0*y0;
  if (m > 0) {
    y0b = FermiAsymRM1(a, y0, 0);
  }
  r[n] = r[0]/FermiAsymRM1(a, y0, m);
  if (m > 0) r[n] *= y0a/y0b;
  else {
    y0a = FermiAsymRM1(a, -y0, 0);
    rn[n] = rn[0]/y0a;
  }
  y0 = exp(y[n-1]);
  if (m > 0) {
    y0a = FermiAsymRM1(a, y0, 0);
  }
  y0 += 0.5*g*y0*y0;
  if (m > 0) {
    y0b = FermiAsymRM1(a, y0, 0);
  }
  r[n+1] = r[n-1]/FermiAsymRM1(a, y0, m);
  if (m > 0) r[n+1] *= y0a/y0b;
}

double InterpFermiRM1(double y0, int m) {
  double y, r, a, g;
  int i;
  if (y0 < 0) {
    m = 0;
    i = 3;
  } else {
    if (m == 0) {
      i = 1;
    } else {
      i = 2;
    }
  }
  a = _fermi_rm1[0][NFRM1];
  g = _fermi_rm1[0][NFRM1+1];
  double ya = fabs(y0);
  y = log(ya);
  if (y < _fermi_rm1[0][0]) {
    r = FermiAsymRM1(a, y0, m);
    r *= _fermi_rm1[i][NFRM1];
  } else if (y > _fermi_rm1[0][NFRM1-1]) {
    if (y0 < 0) return -1;
    if (y0 > 1e10) {
      y0 = 1e10;
    }
    y = y0 + 0.5*g*y0*y0;
    r = FermiAsymRM1(a, y, m);
    if (m > 0) {
      double a1 = FermiAsymRM1(a, y, 0);
      double a0 = FermiAsymRM1(a, y0, 0);
      r *= (a1/a0)*_fermi_rm1[i][NFRM1+1];
    } else {
      r *= _fermi_rm1[i][NFRM1+1];
    }
  } else {
    UVIP3P(3, NFRM1, _fermi_rm1[0], _fermi_rm1[i], 1, &y, &r);
  }
  return r;
}

void PrepFermiRM1(double a, double fa, double t) {
  double g;
  if (_relativistic_fermi) {
    g = t*FINE_STRUCTURE_CONST2;
  } else {
    g = 0.0;
  }
  double aa = fabs(a);
  double y0 = log(Min(1.0,aa)*1e-4);
  double y1 = log(Max(1.0,aa)*1e4);
  double dy = (y1-y0)/(NFRM1-1.0);
  int i;
  _fermi_rm1[0][0] = y0;
  for (i = 1; i < NFRM1; i++) {
    _fermi_rm1[0][i] = _fermi_rm1[0][i-1] + dy;
  }
  _fermi_rm1[0][i] = a;
  _fermi_rm1[0][i+1] = g;
  _fermi_rm1[3][i+1] = fa;
  FermiRM1(NFRM1, _fermi_rm1[0], _fermi_rm1[1], _fermi_rm1[3], a, fa, g, 0);
  FermiRM1(NFRM1, _fermi_rm1[0], _fermi_rm1[2], NULL, a, fa, g, 1);
  if (_fermi_rmf[0]) {
    FILE *f = fopen(_fermi_rmf, "w");
    if (f == NULL) {
      printf("cannot open file %s\n", _fermi_rmf);
    } else {
      for (i = 0; i < NFRM2; i++) {
	if (i < NFRM1) {
	  y0 = exp(_fermi_rm1[0][i]);
	} else {
	  y0 = 0;
	}
	fprintf(f, "%3d %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E\n",
		i, _fermi_rm1[0][i], y0, _fermi_rm1[1][i], _fermi_rm1[2][i],
		_fermi_rm1[3][i], a, g, fa);
      }
      fclose(f);
    }
  }
}

double IonSF0(double k, double g, double sk) {
  double sm = 1.211 + 1.079e-2*g;
  double km = 4.152 + 7.51e-4*g;
  double w = -0.978+4.84e-2*g;
  double kl = 1.25;
  double dk = 0.25;
  double wk = 0.5*(1+tanh((k-kl)/dk));
  double k2 = k*k;
  sk = sk*wk + (k2/(k2+3*g))*(1-wk);
  return sk;
}

double OZBridge(double x, double g, double k2) {
  double lg = log(g);
  double lg2 = lg*lg;
  double b0 = 0.258 - 0.0612*lg + 0.0123*lg2 - 1/g;
  double b1 = 0.0269 + 0.0318*lg + 0.00814*lg2;
  double c1 = 0.498 - 0.28*lg + 0.0294*lg2;
  double c2 = -0.412 + 0.219*lg - 0.0251*lg2;
  double c3 = 0.0988 - 0.0534*lg + 0.00682*lg2;
  double x2 = x*x;
  double x4 = x2*x2;
  double x6 = x2*x4;
  double x8 = x4*x4;
  double r = (-b0 + c1*x4 + c2*x6 + c3*x8)*exp(-b1*x2/b0);
  r *= exp(-0.25*k2)*g;
  return r;
}
  
void IonCF(POTENTIAL *pot, double *y, double rmax, int n, double *w) {
  double *t, *t0, *c, *u, *b, *dw1, *dw2, nz, dr, dk, r, k, a1, a2, tol, ut;
  double tmax, tidx;
  int i, iter, nw1, nw2, *iw;

  nw1 = n/2;
  nw2 = (n*5)/8;
  
  t = w;
  c = t + n;
  u = c + n;
  b = u + n;
  t0 = b + n;
  dw1 = t0 + n;
  dw2 = dw1 + nw1;
  iw = (int *)(dw2 + nw2);
  iter = 0;
  tmax = Max(pot->rps, pot->dps);  
  rmax = tmax*rmax;
  tmax *= _icf_tmax;
  tidx = _icf_tidx;
  dr = rmax/n;
  dk = PI/(n*dr);
  iw[0] = 0;
  for (i = 1; i < n; i++) {
    b[i] = i*dr;
    b[i] = pot->ar*pow(b[i], pot->qr) + pot->br*log(b[i]);
  }
  b[0] = b[1]-(b[2]-b[1]);
  nz = pot->nps/pot->zps;
  double qd2 = FOUR_PI*pot->ups*pot->zps*nz/pot->tps;
  double qd = sqrt(qd2);
  for (i = 0; i < pot->maxrp; i++) {
    _dwork[i] = log(y[i]*pot->ups);
  }  
  UVIP3P(3, pot->maxrp, pot->rho, _dwork, n, b, u);  
  for (i = 0; i < n; i++) {
    if (u[i] < -250) u[i] = 0;
    else u[i] = exp(u[i]);
    if (i == 0) {
      r = 0.01*dr;
    } else {
      r = i*dr;
    }
    u[i] -= OZBridge(r/pot->rps, pot->gps, qd2*pot->rps*pot->rps);
    t[i] = 0.0;
  }
  /*
  for (i = 0; i < n; i++) {
    if (i == 0) {
      r = 0.01*dr;
    } else {
      r = i*dr;
    }
    u[i] = (pot->zps*pot->ups/r)*exp(-qd*r)/pot->tps;
    u[i] -= OZBridge(r/pot->rps, pot->gps, qd2*pot->rps*pot->rps);
    t[i] = 0;
  }
  */
  a1 = FOUR_PI*dr;
  a2 = dk/(2*PI*PI);
  double tt = 0, tt0 = 0;
  double wk, k2, sk;
  while (1) {
    c[0] = 0.0;
    for (i = 1; i < n; i++) {
      t0[i] = t[i];
      ut = u[i]-t[i];
      if (_icf_ozc == 1 || ut >= 0) {
	c[i] = ExpM1(-ut)-t[i];
      } else {
	c[i] = -u[i];
      }
      c[i] *= i*dr;
    }
    dfst(n, c, dw1, iw, dw2);    
    for (i = 1; i < n; i++) {      
      k = i*dk;
      c[i] *= a1/k;
      t[i] = c[i]/(1-nz*c[i]) - c[i];
      t[i] *= k;
    }
    dfst(n, t, dw1, iw, dw2);
    tt0 = tt;
    tt = 0.0;
    for (i = 1; i < n; i++) {
      r = i*dr;
      t[i] *= a2/r;
      tt += t[i]*t[i];
    }
    tol = fabs(tt-tt0)/Max(tt,tt0);
    if (_sp_print == 2 || _sp_print == 3) {
      printf("icf: %d %d %g %g %g %g %g %g %g\n", iter, n, nz, dr, qd*pot->rps, tt0, tt, tol, rmax);
    }
    if (tol < _icf_tol) break;
    for (i = 1; i < n; i++) {
      t[i] = t0[i]*_icf_stablizer + t[i]*(1-_icf_stablizer);
    }
    iter++;
    if (iter > _icf_maxiter) {
      printf("maxiter reached in IonCF: %d %g\n", iter, tt-tt0);
      Abort(1);
    }    
  }
  t[0] = t[1];
  for (i = 0; i < n; i++) {
    sk = 1/(1-nz*c[i]);
    k = i*dk;
    k2 = k*k;
    wk = 0.5*(1+tanh((k*pot->rps-_icf_kb)/_icf_kd));
    sk = wk*sk + (1-wk)*k2/(k2+qd2);
    c[i] = u[i]-t[i];
    u[i] = sk;
    r = i*dr/tmax;
    r = pow(r, tidx);
    t[i] *= exp(-r);
  }
  c[0] = c[1];
  u[0] = u[1];
}

double StewartPyattIntegrand(double x, double a, double fa, double y,
			     double y0, double g, double z, double nb) {
  double ys, zs, r;
  int m;
  if (y0 > 0 || nb >= 0) {
    m = 1;
  } else {
    m = 0;
  }
  r = InterpFermiRM1(y, m);
  double r0 = 0.0;
  if (y0 <= 0 && nb >= 0) {
    r0 = InterpFermiRM1(y, 0);
    r += nb;
  }
  //printf("r01: %g %g %g %g %g\n", x, y, r0, r, nb);
  if (_sp_nzs <= 0) {
    if (_icf_nfft > 0) {
      double *rw = _icf_dw + _icf_nfft*3;
      if (x >= rw[1] && x <= rw[_icf_nfft-1]) {
	double *yw = _icf_dw;// + _icf_nfft;
	UVIP3P(3, _icf_nfft, rw, yw, 1, &x, &zs);
	ys = z*y - zs;
      } else {
	ys = z*y;
      }
    } else {
      ys = z*y;
    }
    if (z > 0) {
      if (ys >= -1) {
	r -= ExpM1(-ys);
      } else {
	r -= ExpM1(1.0)-(ys+1.0);
      }
    } else {
      if (x < 1) {
	r += 1.0;
      }
    }
  } else {
    int i;
    ys = 0.0;
    zs = 0.0;
    for (i = 0; i < _sp_nzs; i++) {
      zs += _sp_zw[i]*_sp_zs[i];
      ys += _sp_zw[i]*_sp_zs[i]*ExpM1(-_sp_zs[i]*y);
    }
    ys /= zs;
    r -= ys;
  }
  r /= (1+z);
  return r;
}

void DerivSPY(int *neq, double *t, double *y, double *yd) {
  double x = exp(*t);
  yd[0] = x*y[1];
  double rho = 0.0;  
  if (_icf_nfft > 0) {
    rho = x*y[9];
    rho = y[6]*pow(rho, y[8]) + y[7]*log(rho);
  } else {
    rho = x/y[10];
  }
  double s = StewartPyattIntegrand(rho, y[2], y[3], y[0]/x, 0.0,
				   y[4], y[5], y[11]);
  yd[1] = x*x*s;
}

FCALLSCSUB4(DerivSPY, DERIVSPY, derivspy, PINT, PDOUBLE, DOUBLEV, DOUBLEV)
double StewartPyatt(POTENTIAL *pot, double *vt, double *wb, double xps) {
  int i, m, ib;
  double a = pot->aps;
  double fa = pot->fps;
  double x, x0, x1, x2, y, y0, k, j, c0, c, z1, g, kc, t0, t;
  double *vb = ABAND;
  
  if (vt) {
    if (_sp_mode == 0 ||
	_sp_mode == 2 ||
	//_sp_mode == 3 ||
	_sp_mode == 4) return 1.0;
  }
  kc = 1.0;
  if (_relativistic_fermi) {
    g = pot->tps*FINE_STRUCTURE_CONST2;
  } else {
    g = 0.0;
  }
  m = pot->maxrp-1;
  double ys[12];
  ys[0] = 0;
  k = pot->zps/(pot->dps*pot->tps);
  z1 = 1.0 + pot->ups;
  x = 3*z1;
  j = (pow(x*k+1.0, 0.66666666667)-1)*0.5/z1;
  x1 = sqrt(2*z1*j+1)-1.0;
  c0 = exp(x1)*(j-x1*x1/(2*z1));
  if (xps > 0) {
    c = c0*xps;
  }
  x0 = pot->rps/pot->dps;
  x2 = _sp_rmax*x1;
  x = Max(1.0, x0);
  x2 = Max(x2, x);
  y0 = -1e30;
  if (xps > 0) {
    if (_sp_mode >= 2) {
      double *rwork = _dwork3;
      int iwork[22];
      double rtol, atol[2];
      int neq, itol, itask, istate, iopt, mf, lrw, liw;
      lrw = pot->maxrp;
      liw = 22;
      neq = 2;
      itol = 2;
      rtol = EPS3;
      atol[0] = EPS6;
      atol[1] = EPS6;
      itask = 1;
      istate = 1;
      iopt = 0;
      mf = 10;
      for (i = pot->maxrp-1; i > 10; i--) {      
	x = pot->rad[i]/pot->dps;
	y0 = z1*c0*exp(-x)/x;
	_veff[i] = c*exp(-x)/x;
	//if (pot->rad[i] < x2) break;
	if (y0 > _sp_yeps) {
	  break;
	}
      }
      ys[2] = a;
      ys[3] = fa;
      ys[4] = g;
      ys[5] = pot->ups;
      ys[6] = pot->ar;
      ys[7] = pot->br;
      ys[8] = pot->qr;
      ys[9] = pot->dps;
      ys[10] = pot->rps/pot->dps;
      //printf("x2: %d %g %g %g %g %g %g %g\n", i, pot->rad[i], x, x1, x2, y0, ys[10], k);
      ys[0] = c*exp(-x);
      ys[1] = -c*exp(-x);
      x0 = x;
      t0 = log(x0);
      i--;
      ib = 0;
      for (; i >= 0; i--) {
	x = pot->rad[i]/pot->dps;
	t = log(x);
	if (wb) {
	  ys[11] = (wb[i]+wb[i+1])*0.5/pot->nps;
	} else {
	  ys[11] = -1;
	}
	while (t0 > t) {
	  LSODE(C_FUNCTION(DERIVSPY, derivspy), neq, ys, &t0, t,
		itol, rtol, atol, itask, &istate, iopt, rwork, lrw, iwork,
		liw, NULL, mf);
	  if (istate == -1) istate = 2;
	  else if (istate < 0) {
	    printf("PSODE error in StewartPyatt: %d %d %g",
		   istate, vt==NULL, rwork[13]);
	    Abort(1);
	  }
	}
	x0 = x;
	_veff[i] = ys[0]/x;
	if (_veff[i] > _sp_vmax) {
	  if (vt) {
	    kc = vb[i]/ys[0];
	  } else {
	    kc = 1.0;
	  }
	  ib = i;
	  break;
	}
	if (_sp_mode >= 4 && vt) {
	  y = _veff[i];
	  x2 = exp(-pot->zps*y);
	  if (x2 < _sp_neps) {
	    kc = vb[i]/ys[0];
	    ib = i;
	    //printf("ib: %d %g %g %g %g %g %g\n", i, pot->rad[i], x, y, pot->NPS[i], pot->EPS[i], kc);
	    break;
	  }
	}
      }
      if (ib > 0) {
	if (_sp_mode < 6) {
	  for (i--; i >= 0; i--) {
	    x = pot->rad[i]/pot->dps;
	    if (vt) {
	      _veff[i] = vb[i]/x/kc;
	    } else {
	      _veff[i] = _veff[ib]*pot->rad[ib]/pot->rad[i];
	    }
	  }
	} else {
	  for (i--; i >= 0; i--) {
	    x = pot->rad[i]/pot->dps;
	    if (vt) {
	      _veff[i] = vb[i]/x;
	    } else {
	      _veff[i] = _veff[ib]*pot->rad[ib]/pot->rad[i];
	    }
	  }
	  for (i = ib; i < pot->maxrp; i++) {
	    x = pot->rad[i]/pot->dps;
	    _veff[i] *= kc;
	  }
	}
      }
    } else {
      for (i = 0; i < pot->maxrp; i++) {
	x = pot->rad[i]/pot->dps;
	if (vt == NULL) {
	  if (x >= x1) {
	    y = c*exp(-x)/x;
	  } else {
	    y = ((k/x) - j + x*x/(6.0*z1))*kc;
	  }
	} else {
	  y = fabs(vt[i])/pot->tps;
	  if (x2 > 0 && x >= x2) {
	    if (y0 < 0) {
	      y0 = c*exp(-x)/x;	
	      y0 = y/y0;
	    }
	    y = c*y0*exp(-x)/x;
	  }
	}
	_veff[i] = y;
      }
    }
  }
  /*
  if (_icf_nfft > 0) {
    IonCF(pot, _veff, _icf_rmax, _icf_nfft, _icf_dw);
  }
  */
  for (i = 0; i <= m; i++) {
    x = pot->rad[i]/pot->dps;
    y = _veff[i];
    if (pot->ups <= 0) {
      t = pot->rad[i]/pot->rps;
    } else {
      t = pot->rho[i];
    }
    _dwork1[i] = StewartPyattIntegrand(t, a, fa, y, y, g, pot->ups, -1);
    pot->NPS[i] = _dwork1[i];
    pot->EPS[i] = (1+InterpFermiRM1(y, 1))/(1+pot->ups);
    //pot->EPS[i] = StewartPyattIntegrand(t, a, fa, y, 0, g, pot->ups);
    _dwork1[i] *= pot->dr_drho[i]*x/pot->dps;
    _dwork2[i] = _dwork1[i]*x;
  }
  _dwork[m] = 0.0;
  NewtonCotesIP(_dwork, _dwork1, 0, m, -1, -1);
  x = pot->rad[0]/pot->dps;
  _dwork3[0] = _dwork1[0]*x*x*x/3.0;
  NewtonCotesIP(_dwork3, _dwork2, 0, m, -1, 0);
  double cs = pot->zps/_dwork3[m];
  double ds = cs/(FOUR_PI*pow(pot->dps,3));
  double xs = _dwork3[m]/k;
  if (_sp_print == 1 || _sp_print == 3) {
    printf("cs: %g %g %g %g %g %g %g %g %g %g %g %g\n", pot->zps, _dwork3[m], c, k, xps, xs, cs, x1, x2, kc, ys[0], _veff[0]);
  }
  for (i = 0; i <= m; i++) {
    x = pot->rad[i]/pot->dps;
    y = -cs*(_dwork3[i] + x*_dwork[i]);
    pot->ZPS[i] = y;
    pot->NPS[i] *= ds;
    pot->EPS[i] *= ds;
  }
  pot->jps = _dwork[0]*pot->tps;
  return xs;
}

void SetLatticePotential(POTENTIAL *pot, double *vt) {
  int n = 3;
  int m = 500;
  int i, j, ix, iy, iz;
  double r, x, y, z, v;
  double dx, dy, dz, dr, rho, u1, u2, vi;

  for (i = 0; i < pot->maxrp; i++) {
    _dwork[i] = vt[i]*pot->rad[i];
  }
  for (i = 0; i < pot->maxrp; i++) {
    v = 0.0;
    for (j = 0; j < m; j++) {
      r = pot->rad[i];
      if (i < pot->maxrp-1) {
	u1 = drand48();
	r += u1*(pot->rad[i+1]-pot->rad[i]);
      }
      u1 = 2*drand48()-1;
      u2 = TWO_PI*drand48();
      x = r*u1*cos(u2);
      y = r*u1*sin(u2);
      z = r*sqrt(1-u1*u1);
      for (ix = -n; ix <= n; ix++) {
	for (iy = -n; iy <= n; iy++) {
	  for (iz = -n; iz <= n; iz++) {
	    if (ix == 0 && iy == 0 && iz == 0) continue;
	    dx = x-ix*pot->rps;
	    dy = y-iy*pot->rps;
	    dz = z-iz*pot->rps;
	    dx *= dx;
	    dy *= dy;
	    dz *= dz;
	    dr = sqrt(dx+dy+dz);
	    rho = pot->ar*pow(dr, pot->qr) + pot->br*log(dr);
	    UVIP3P(3, pot->maxrp, pot->rho, _dwork, 1, &rho, &vi);
	    v += vi/dr;
	  }
	}
      }
    }
    v /= m;
    printf("zps: %d %g %g %g %g\n", i, pot->rad[i], vt[i]*pot->rad[i], pot->ZPS[i], v*pot->rad[i]);
    pot->ZPS[i] -= v*pot->rad[i];
  }
  Differential(pot->ZPS, pot->dZPS, 0, pot->maxrp-1, pot->dr_drho);
  Differential(pot->dZPS, pot->dZPS2, 0, pot->maxrp-1, pot->dr_drho);
}

void SetOptionOrbital(char *s, char *sp, int ip, double dp) {
  if (0 == strcmp(s, "orbital:relativistic_fermi")) {
    _relativistic_fermi = ip;
    return;
  }
  if (0 == strcmp(s, "orbital:fermi_abserr")) {
    _fermi_abserr = dp;
    return;
  }
  if (0 == strcmp(s, "orbital:fermi_relerr")) {
    _fermi_relerr = dp;
    return;
  }
  if (0 == strcmp(s, "orbital:fermi_stablizer")) {
    _fermi_stablizer = dp;
    return;
  }
  if (0 == strcmp(s, "orbital:fermi_rmf")) {
    strncpy(_fermi_rmf, sp, 256);
    return;
  }
  if (0 == strcmp(s, "orbital:sp_ofn")) {
    strncpy(_sp_ofn, sp, 256);
    return;
  }
  if (0 == strcmp(s, "orbital:sp_rmax")) {
    _sp_rmax = dp;
    return;
  }
  if (0 == strcmp(s, "orbital:sp_yeps")) {
    _sp_yeps = dp;
    return;
  }
  if (0 == strcmp(s, "orbital:sp_neps")) {
    _sp_neps = dp;
    return;
  }
  if (0 == strcmp(s, "orbital:sp_z0")) {
    _sp_z0 = dp;
    return;
  }
  if (0 == strcmp(s, "orbital:sp_z1")) {
    _sp_z1 = dp;
    return;
  }
  if (0 == strcmp(s, "orbital:sp_mode")) {
    _sp_mode = ip;
    return;
  }
  if (0 == strcmp(s, "orbital:sp_print")) {
    _sp_print = ip;
    return;
  }
  if (0 == strcmp(s, "orbital:icf_tol")) {
    _icf_tol = dp;
    return;
  }
  if (0 == strcmp(s, "orbital:icf_maxiter")) {
    _icf_maxiter = ip;
    return;
  }  
  if (0 == strcmp(s, "orbital:icf_rmax")) {
    _icf_rmax = dp;
    return;
  }  
  if (0 == strcmp(s, "orbital:icf_kd")) {
    _icf_kd = dp;
    return;
  }
  if (0 == strcmp(s, "orbital:icf_kb")) {
    _icf_kb = dp;
    return;
  }
  if (0 == strcmp(s, "orbital:icf_tmax")) {
    _icf_tmax = dp;
    return;
  }
  if (0 == strcmp(s, "orbital:icf_tidx")) {
    _icf_tidx = dp;
    return;
  }
  if (0 == strcmp(s, "orbital:icf_stablizer")) {
    _icf_stablizer = dp;
    return;
  }
  if (0 == strcmp(s, "orbital:icf_spmi")) {
    _icf_spmi = ip;
    return;
  }
  if (0 == strcmp(s, "orbital:icf_ozc")) {
    _icf_ozc = ip;
    return;
  }
  if (0 == strcmp(s, "orbital:icf_nfft")) {
    if (_icf_nfft > 0) {
      free(_icf_dw);
      _icf_nfft = 0;
      _icf_dw = NULL;
    }
    if (ip >= 0) {
      if (ip > 5) ip = 5;
      _icf_nfft = 1<<(10+ip);
      _icf_dw = malloc(sizeof(double)*_icf_nfft*7);
      int *iw;
      iw = (int *)(_icf_dw + 5*_icf_nfft + (_icf_nfft/8)*9);
      iw[0] = 0;
    }
    return;
  }  
  if (0 == strcmp(s, "orbital:icf_ofn")) {
    strncpy(_icf_ofn, sp, 256);
    return;
  }
}

double SetSPZW(int n, double *zw) {
  if (_sp_nzs > 0) {
    free(_sp_zs);
    free(_sp_zw);
    _sp_zs = NULL;
    _sp_zw = NULL;
    _sp_nzs = 0;
  }
  if (n <= 0) return -1;
  int n2 = (n+1)/2;
  _sp_nzs = n2;
  _sp_zs = malloc(sizeof(double)*n2);
  _sp_zw = malloc(sizeof(double)*n2);
  int i, k;
  double w = 0, z1 = 0, z2 = 0;
  for (i = 0; i < n2; i++) {
    k = 2*i;
    _sp_zs[i] = zw[k];
    if (k+1 < n) {
      _sp_zw[i] = zw[k+1];
    } else {
      _sp_zw[i] = 1.0;
    }
    w += _sp_zw[i];
    z1 += _sp_zs[i]*_sp_zw[i];
    z2 += _sp_zs[i]*_sp_zs[i]*_sp_zw[i];
  }
  for (i = 0; i < n2; i++) {
    _sp_zw[i] /= w;
  }
  z1 /= w;
  z2 /= w;
  
  return z2/z1;
}
