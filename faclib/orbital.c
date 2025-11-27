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
static double *_veff, *_veff0, *ABAND, *_dwork, *_dwork1, *_dwork2, *_dwork3;
static int *_ipiv;

#pragma omp threadprivate(_veff, _veff0, ABAND, _ipiv, _dwork, _dwork1, _dwork2, _dwork3)

static double _relativistic_xfermi = 0.025;
static int _relativistic_fermi0 = -1;
static int _relativistic_fermi = 0;
static double _fermi_abserr = 1e-7;
static double _fermi_relerr = 1e-5;
#define NFRM1 256
#define NFRM2 (NFRM1+2)
#define NFRM3 201
#define NFRM4 651
static double _fermi_nr[8][NFRM4];
static int _nfermi_nr[8];
static double _fermi_a0 = 1.0;
static double _fermi_a1 = 10.0;
static double _fermi_y0 = 1e-4;
static double _fermi_y1 = 1e5;
static double _fermi_tmin = 0.01;
static double _fermi_tmax = 100.0;
static double _fermi_rm1[4][NFRM2];
static double _fermi_ymx;
static double _fermi_rmx;
static char _fermi_rmf[256] = "";
static char _sp_ofn[256] = "";
static double _sp_rmax = 1.0;
static double _sp_yeps = 1e-5;
static double _sp_neps = 1e-3;
static int _sp_nzs = 0;
static double *_sp_zs = NULL;
static double *_sp_zw = NULL;
static double _sp_zu = 0.0;
// _sp_mode = 1, use SP approx. solution for y
// _sp_mode = 2, integrate SP ode, only on the 1st optimize loop
// _sp_mode = 3, integrate SP ode, on every optimize loop
static int _sp_mode = 3;
static int _debye_mode = 0; 
static int _sp_print = 0;
static int _sp_piter0 = 0;
static int _sp_piter1 = 0;
static int _ionsph_miniter = 10;
static int _ionsph_maxiter = 1024;
static double _ionsph_ifermi = 1.0;
// _ionsph_bmode=0, mu not updated during SP solution iteration
// _ionsph_bmode=1, mu updated, just as regular ionsph model
static int _ionsph_bmode = 0;
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
static int _on_error = 0;
static double _zcoll = 0;
static double _mcoll = 0;
static int _pwa = 0;
static double _emin_amp = 0.05;
static double _sturm_rmx = 10.0;
static double _sc_bqp = 1E31;
static double _sc_rbf = 1.0;
static double _sc_rsf = 1.0;
static double _sc_ewf = 1.0;
static int _sc_ewm = 0;
static double _sc_ewr = 0.0;
static int _debug = 0;
static double _matchtol = 1e-3;
static int _veff_corr = 1;
static int _rydnorm = 1;

#define MINNKF 55
#define MAXNKF 155
static double _kwork[MAXNKF];

static double _enerelerr = ENERELERR;
static double _eneabserr = ENEABSERR;
static double _enerelerr1 = ENERELERR1;

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
static int SetVEffectiveZMC(int kl, int kv, POTENTIAL *pot);
static int TurningPoints(int n, double e, int kappa, POTENTIAL *pot, int isol);
static int IntegrateRadial(double *p, double e, POTENTIAL *pot, 
			   int i1, double p1, int i2, double p2, int q);
static double Amplitude(double *p, double e, int kl, POTENTIAL *pot, int i1);
static int Phase(double *p, POTENTIAL *pot, int i1, double p0);

int IonSphBMode() {
  return _ionsph_bmode;
}

int OrbPWA() {
  return _pwa;
}

double SPZU() {
  return _sp_zu;
}

double SCEWF() {
  return _sc_ewf;
}

int SCEWM() {
  return _sc_ewm;
}

double SCRSF() {
  return _sc_rsf;
}

double SCRBF() {
  return _sc_rbf;
}

double SCBQP() {
  return _sc_bqp;
}

double ZColl() {
  return _zcoll;
}

double MColl() {
  return _mcoll;
}

void SetZColl(double z) {
  _zcoll = z;
}

void SetMColl(double m) {
  if (m > 0) {
    _mcoll = m*AMU;
    //SetBornMass(m);
  } else if (m < 0) {
    _mcoll = -m;
    SetBornMass(m);
  } else {
    //_mcoll = 0.0;
  }
}

int OnErrorOrb() {
  return _on_error;
}

void SetOnErrorOrb(int e) {
  _on_error = e;
}

int SPMode() {
  return _sp_mode;
}

void SetOrbitalWorkSpace(double *p, int n) {
  _veff0 = p;
  p += n;
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
  e = Max(e, _eneabserr);
  double d0 = e*_enerelerr;
  if (d0 > _eneabserr) {
    d0 = e*_enerelerr1;
    d0 = Max(_eneabserr, d0);
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
  orb->isol = 0;
  if (orb->n > 0) {
    if (orb->n == 1000000) {
      ierr = RadialFreeInner(orb, pot);
    } else {
      if (pot->ib > 0 && orb->n > pot->nb) {
	if (pot->sturm_idx > 0) {
	  ierr = RadialSturm(orb, pot);
	  if (orb->isol < 0) orb->isol = 0;
	} else {
	  ierr = RadialBasis(orb, pot);
	}
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
	if (orb->n <= pot->nmax) {
	  ierr = RadialBound(orb, pot);
	} else {
	  ierr = RadialRydberg(orb, pot);
	}
	if (!ierr) {
	  orb->dn = QuantumDefect(z, orb->n, orb->kappa, orb->energy);
	}
      }
      if (!ierr) {
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
    }
  } else if (orb->n == 0) {
    if (orb->energy <= 0) {
      printf("continuum invalid energy: %d %d %.10e\n", orb->n, orb->kappa, orb->energy);
      orb->wfun = NULL;
      return -97;
    }
    ierr = RadialFree(orb, pot);
  } else {
    ierr = RadialBasisOuter(orb, pot);
  }
  if (orb->wfun != NULL && orb->ilast > 0) orb->isol = 1;
  else if (ierr == 0) {
    printf("orb ilast not set: %d %d %g %d\n",
	   orb->n, orb->kappa, orb->energy, orb->ilast);
    return -98;
  }
  if (orb->wfun && _on_error >= 0) {
    if (isnan(orb->wfun[0]) || isnan(orb->wfun[pot->maxrp])) {
      MPrintf(-1, "wfun nan: %d %d %12.5E %12.5E %12.5E\n",
	      orb->n, orb->kappa, orb->energy,
	      orb->wfun[0], orb->wfun[pot->maxrp]);    
      return -99;
    }
  }
  return ierr;
}

double *GetVEffective(void) { 
  return _veff;
}

int FirstMaximum(double *p, int i1, int i2, POTENTIAL *pot) {
  int i, im;
  double fm, fi;

  fm = 0.0;
  im = -1;
  for (i = i1; i <= i2; i++) {
    if (pot->rad[i] < pot->atom->rms0) continue;
    fi = fabs(p[i]);
    if (fm < fi) {
      fm = fi;
      im = i;
    } else {
      break;
    }
  }
  return im;
}

int LastMaximum(double *p, int i1, int i2, POTENTIAL *pot) {
  int i, im;
  double fm, fi;

  fm = 0.0;
  im = -1;
  for (i = i1; i <= i2; i++) {
    if (pot->rad[i] < pot->atom->rms0) continue;
    fi = fabs(p[i]);
    if (fm < fi) {
      fm = fi;
      im = i;
    }
  }
  return im;
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

void Differential(double *p, double *dp, int i1, int i2,
		  POTENTIAL *pot) {
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
    dp[i] /= pot->dr_drho[i];
  }
}    

double DpDrPWA(int kappa, int i, double e, POTENTIAL *pot,
	       double b, int m, double *bqp) {
  double x2, dr, pr, kl;
  
  x2 = 1 + 0.5*FINE_STRUCTURE_CONST2*e;
  if (m == 0) {
    b = (b + kappa/pot->rad[i])*FINE_STRUCTURE_CONST*0.5;
    if (bqp) *bqp = b;
    b = 2*x2*b/FINE_STRUCTURE_CONST - kappa/pot->rad[i];
  } else if (m == 1) {
    if (bqp) *bqp = b;
    b = 2*x2*b/FINE_STRUCTURE_CONST - kappa/pot->rad[i];
  } else if (m == -1) {
    kl = GetLFromKappa(kappa);
    if (kappa < 0) {
      b = e*FINE_STRUCTURE_CONST;
      b *= -pot->rad[i]/(kl+3.0);
    } else {
      b =(e + 2/FINE_STRUCTURE_CONST2)*FINE_STRUCTURE_CONST;
      b *= pot->rad[i]/(kl+1.0);
      b = 1.0/b;
    }   
    if (bqp) {
      *bqp = b;
    }
    b = (1 + kl/2)/pot->rad[i];    
  } else if (m == -2) {
    b = -sqrt(1.0/FINE_STRUCTURE_CONST2 - e*e*FINE_STRUCTURE_CONST2);
    if (bqp) {
      *bqp = (b + kappa/pot->rad[i])*FINE_STRUCTURE_CONST/(2*x2);
    }
  }
  b *= pot->dr_drho[i];
  dr = 0.5*pot->dr_drho[i]/pot->rad[i];
  dr *= (1-pot->qr*pot->qr*pot->dr_drho[i]*pot->ar*pow(pot->rad[i], pot->qr-1));
  pr = b - dr;
  
  return pr;
}

double DpDr(int kappa, int k, int i, double e, POTENTIAL *pot,
	    double b, int m, double *bqp) {
  double x2, dx, dr, pr, z0, kl;

  if (k < 0) {
    return DpDrPWA(kappa, i, e, pot, b, m, bqp);
  }
  x2 = 1 + 0.5*FINE_STRUCTURE_CONST2*(e - pot->VT[k][i]);
  if (m == 0) {
    b = (b + kappa/pot->rad[i])*FINE_STRUCTURE_CONST*0.5;
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
  double e, ep, de, delta, emin, emax, dr, ke;
  double *p, p0, p1, p2, qo, qi, bqp, bqp1, norm2, fact;
  int n, i, k, kl, nr, nodes, niter;
  int i2, i2o, i2m1, i2m2, i2p1, i2p2, i1;
  int ib0, ib1;

  int kv = 0;
  if (pot->pse) kv = IdxVT(orb->kappa);
  orb->kv = kv;
  kl = orb->kappa;
  kl = (kl < 0)? (-kl-1):kl;

  n = -orb->n;
  nr = n - kl - 1;
  if (kl < 0 || kl >= n) {
    if (_on_error >= 0) {
      printf("Invalid orbital angular momentum, L=%d, %d %d\n", 
	     kl, orb->n, orb->kappa);
    }
    return -1;
  }
  
  p = malloc(sizeof(double)*2*pot->maxrp);
  if (!p) return -1;
  orb->wfun = p;
  
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
  ke = TWO_PI*(5*nr+5)/dr;
  emax = 0.5*ke*ke;
  e = emin + 0.25*emax;
  emax = Max(e, emax);
  ke = TWO_PI*(nr+1)/dr;
  e = 0.5*ke*ke;
  de = 0.5*e;
  while (niter < max_iteration) {
    niter++;
    SetPotentialW(pot, e, orb->kappa, kv);
    SetVEffective(kl, kv, pot);
    bqp = DpDr(orb->kappa, kv, pot->ib, e, pot, pot->bqp, 0, NULL);
    bqp1 = DpDr(orb->kappa, kv, pot->ib1, e, pot, pot->bqp, 0, NULL);
    i2 = TurningPoints(orb->n, e, orb->kappa, pot, 0);
    nodes = IntegrateRadial(p, e, pot, ib0, bqp, i2, 1.0, 2);
    if (nodes > 1) {
      i2 = LastMaximum(p, pot->ib+20, i2, pot);
      nodes = CountNodes(p, pot, ib0, i2);
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
    if (_on_error >= 0) {
      printf("Max iteration before finding correct nodes in RadialBasisOuter1 %d %d %d %d %g %g %g\n",
	     nodes, nr, orb->n, orb->kappa, emin, emax, e);
    }
    //free(p);
    return -2;
  }
  
  niter = 0;
  de = emax-emin;
  ep = _eneabserr;
  while (niter < max_iteration) {
    while (nodes == nr && niter < max_iteration) {
      niter++;
      e += de;
      SetPotentialW(pot, e, orb->kappa, kv);
      SetVEffective(kl, kv, pot);
      bqp = DpDr(orb->kappa, kv, pot->ib, e, pot, pot->bqp, 0, NULL);
      bqp1 = DpDr(orb->kappa, kv, pot->ib1, e, pot, pot->bqp, 0, NULL);
      i2 = TurningPoints(orb->n, e, orb->kappa, pot, 0);
      nodes = IntegrateRadial(p, e, pot, ib0, bqp, i2, 1.0, 2);
      if (nodes > 1) {
	i2 = LastMaximum(p, pot->ib+20, i2, pot);
	nodes = CountNodes(p, pot, ib0, i2);
      }
      nodes += IntegrateRadial(p, e, pot, i2, 1.0, ib1, bqp1, 1);
    }
    if (nodes-nr == 1 && de < ep) break;
    e -= de;
    nodes = nr;
    de *= 0.25;
  }
  if (niter == max_iteration) {
    if (_on_error >= 0) {
      printf("Max iteration before finding solution in RadialBasisOuter2 %d %d %d %d %d %g %g\n", niter, orb->n, orb->kappa, nodes, nr, e, de);
    }
    //free(p);
    return -3;
  }
  emax = e;
  
  if (nr > 0) {
    e -= de;
    nodes = nr;
    de = emax-emin;
    niter = 0;
    while (niter < max_iteration) {
      while (nodes == nr && niter < max_iteration) {
	niter++;
	e -= de;
	SetPotentialW(pot, e, orb->kappa, kv);
	SetVEffective(kl, kv, pot);
	bqp = DpDr(orb->kappa, kv, pot->ib, e, pot, pot->bqp, 0, NULL);
	bqp1 = DpDr(orb->kappa, kv, pot->ib1, e, pot, pot->bqp, 0, NULL);
	i2 = TurningPoints(orb->n, e, orb->kappa, pot, 0);
	nodes = IntegrateRadial(p, e, pot, ib0, bqp, i2, 1.0, 2);
	if (nodes > 1) {
	  i2 = LastMaximum(p, pot->ib+20, i2, pot);
	  nodes = CountNodes(p, pot, ib0, i2);
	}
	nodes += IntegrateRadial(p, e, pot, i2, 1.0, ib1, bqp1, 1);
	//printf("nd1: %d %d %d %d %d %d %15.8E %15.8E %15.8E\n", niter, orb->n, orb->kappa, i2, nodes, nr, e, de, ep);
      }
      if (nodes-nr == -1 && de < ep) break;
      e += de;
      nodes = nr;
      de *= 0.25;
    }  
    if (niter == max_iteration) {
      if (_on_error >= 0) {
	printf("Max iteration before finding solution in RadialBasisOuter3 %d %d %d %d %d %g %g\n", niter, orb->n, orb->kappa, nodes, nr, e, de);
      }
      //free(p);
      return -4;
    }
    emin = e;
  }

  niter = 0;
  while (niter < max_iteration) {
    niter++;
    e = 0.5*(emin+emax);
    SetPotentialW(pot, e, orb->kappa, kv);
    SetVEffective(kl, kv, pot); 
    bqp = DpDr(orb->kappa, kv, pot->ib, e, pot, pot->bqp, 0, NULL);
    bqp1 = DpDr(orb->kappa, kv, pot->ib1, e, pot, pot->bqp, 0, NULL);
    i2o = TurningPoints(orb->n, e, orb->kappa, pot, 0);
    i2p2 = i2o + 2;
    nodes = IntegrateRadial(p, e, pot, ib0, bqp, i2p2, 1.0, 2);
    if (nodes > 1) {
      i2 = LastMaximum(p, pot->ib+20, i2o, pot);
    } else {
      i2 = i2o;
    }
    i2m1 = i2 - 1;
    i2m2 = i2 - 2;
    i2p1 = i2 + 1;    
    i2p2 = i2 + 2;
    p1 = p[i2m2];
    p0 = p[i2m1];
    p2 = p[i2];
    qo = (-4.0*p[i2m2-1] + 30.0*p[i2m2] - 120.0*p[i2m1]
	  + 40.0*p[i2] + 60.0*p[i2p1] - 6.0*p[i2p2])/120.0;
    qo /= p[i2];
    k = IntegrateRadial(p, e, pot, i2m2, 1.0, ib1, bqp1, 1);
    qi = (6.0*p[i2m2] - 60.0*p[i2m1] - 40.0*p[i2] + 120.0*p[i2p1]
	  - 30.0*p[i2p2] + 4.0*p[i2p2+1])/120.0;
    qi /= p[i2];
    delta = qo - qi;
    fact = p[i2]/p2;
    p[i2m1] = p0;
    p[i2m2] = p1;
    for (i = 0; i < i2; i++) p[i] *= fact;
    nodes = CountNodes(p, pot, ib0, ib1); 
    //printf("rb: %d %d %d %d %d %d %d %d %d %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E\n", niter, orb->n, orb->kappa, i2, i2o, ib0, ib1, nodes, nr, pot->rad[i2], e, _veff[i2], delta, emin, emax, qo, qi, bqp, bqp1); 
    if (nodes > nr) {
      emax = e;
      continue;
    } else if (nodes < nr) {
      emin = e;
      continue;
    }
    if (delta < 0) {
      emax = e;
    } else if (delta > 0) {
      emin = e;
    } else {
      break;
    }
    if (fabs(emax-emin) < EneTol(e)*0.1) break;
  }
  
  if (niter == max_iteration) {
    if (_on_error >= 0) {
      printf("Max iteration before finding solution in RadialBasisOuter %d %d\n",
	     nodes, nr);
    }
    //free(p);
    return -5;
  }
  
  for (i = ib0; i <= ib1; i++) {
    p[i] *= pot->dr_drho2[i];
  }
  if (p[pot->ib] < 0) {
    for (i = ib0; i <= ib1; i++) {
      p[i] = -p[i];
    }
  }
  norm2 = InnerProduct(ib0, ib1, p, p, pot);
  if (fabs(delta) > _matchtol) {
    if (_on_error >= 0) {
     printf("Discontinuous Match in RadialBasisOuter: %d %d %g %g\n",
	    orb->n, orb->kappa, delta, norm2);
    }
    return -6;
  }
  orb->ilast = pot->ib1;
  orb->energy = e;
  orb->qr_norm = 1.0;
  if (pot->flag == -1) {
    DiracSmall(orb, pot, ib1, kv);
  }

  return 0;
}

int RadialBasisBQP(ORBITAL *orb, POTENTIAL *pot, double pbqp) {
  double z, z0, e, de, ep, delta, emin, emax, emin0, pp, qq, uu, vv;
  double *p, norm2, fact, p0, p1, p2, qo, qi, bqp, bqp0, bqp1;
  int i, k, kl, nr, nodes, niter, mb;
  int i2, i2o, i2m1, i2m2, i2p1, i2p2;

  int kv = 0;
  if (pot->pse) kv = IdxVT(orb->kappa);
  orb->kv = kv;
  kl = orb->kappa;
  kl = (kl < 0)? (-kl-1):kl;

  if (kl < 0 || kl >= orb->n) {
    if (_on_error >= 0) {
      MPrintf(-1, "Invalid orbital angular momentum, L=%d n=%d ka=%d\n", 
	      kl, orb->n, orb->kappa);
    }
    return -1;
  }

  nr = orb->n - kl - 1;

  int ib = pot->ib;
  double rb;
  if (orb->isol == -1) {
    ib = pot->maxrp-4;
  }

  rb = pot->rad[ib];
  
  if (pbqp > 1e10) {
    mb = 0;
  } else if (pbqp > -1e10+EPS5) {
    bqp = (pbqp + orb->kappa/pot->rad[ib])*FINE_STRUCTURE_CONST*0.5;    
    mb = 1;
  } else {
    mb = 1;
    if (pbqp > -1e12+EPS5) {
      if (orb->kappa == -1) {
	bqp = 0.0;
      } else if (orb->kappa == 1) {
	bqp = FINE_STRUCTURE_CONST/pot->rad[ib];
      } else {
	p1 = FINE_STRUCTURE_CONST*(1-orb->kappa); 
	p2 = FINE_STRUCTURE_CONST*(1+orb->kappa);
	qi = pot->rad[ib];
	qo = qi*qi;
	bqp = (1/p1)*(qi-sqrt(qo-p1*p2));
      }
    } else {
      if (pbqp > -1e15+EPS5) {
	p1 = ((pbqp/(-1e12))+1.0)/2;
      } else {
	p1 = (-(pbqp/(-1e15))+3.0)/2;
      }
      p2 = p1 - orb->kappa;
      qi = p1 + orb->kappa;
      qo = FINE_STRUCTURE_CONST/pot->rad[ib];
      if (fabs(p2) < 1e-3) {
	bqp = 0.5*qo*qi;
      } else {
	qi *= p2*qo*qo;
	if (qi > 1) {
	  if (_on_error >= 0) {
	    MPrintf(-1, "Invalid boundary condition in RadialBasis: %d %d %g %g %g\n",
		    orb->n, orb->kappa, pbqp, p1, qi);
	  }
	  return -1;
	}
	qi = 1-sqrt(1-qi);
	bqp = (qi/p2)/qo;
      }
    }  
  }
  if (orb->isol == 0) {
    p = malloc(sizeof(double)*2*pot->maxrp);
    orb->wfun = p;
  } else {
    p = orb->wfun;
  }
  if (!p) {
    if (_on_error >= 0) {
      MPrintf(-1, "cannot alloc memory RadialBasis: %d %d\n",
	      orb->n, orb->kappa);
    }
    return -1;
  }

  niter = 0;
  z0 = GetAtomicNumber();
  z = GetResidualZ();
  z = 0.5*(z0+z);

  emin0 = EnergyH(z0, 1, -1)*1.1;
  ep = EnergyH(z, orb->n, orb->kappa);
  e = fabs(ep);
  while (niter < max_iteration) {
    niter++;
    SetPotentialW(pot, e, orb->kappa, kv);
    SetVEffective(kl, kv, pot);
    i2 = TurningPoints(orb->n, e, orb->kappa, pot, orb->isol);
    if (i2 >= 5) {
      bqp0 = DpDr(orb->kappa, kv, 0, e, pot, 0, -1, &orb->bqp0);    
      nodes = IntegrateRadial(p, e, pot, 0, bqp0, i2, 1.0, 2);
      if (nodes > nr) break;
    }
    e = e*2.0;
  }
  if (niter == max_iteration) {
    if (_on_error >= 0) {
      printf("Max iteration before finding correct nodes in RadialBasis0 %d %d %d %d %d %d %g %g\n",
	     nodes, nr, orb->n, kl, i2, orb->isol, e, emin);
    }
    return -2;
  }
  emax = e;

  e = ep;
  niter = 0;
  while (niter < max_iteration) {
    niter++;
    SetPotentialW(pot, e, orb->kappa, kv);
    SetVEffective(kl, kv, pot);
    i2 = TurningPoints(orb->n, e, orb->kappa, pot, orb->isol);
    if (i2 < 5) {
      nodes = -1;
    } else {
      bqp0 = DpDr(orb->kappa, kv, 0, e, pot, 0, -1, &orb->bqp0);
      nodes = IntegrateRadial(p, e, pot, 0, bqp0, i2, 1.0, 2);
    }
    if (nodes < nr) break;    
    e = e*2.0;
    if (e <= emin0) break;
  }
  emin = e;
  if (niter == max_iteration) {
    if (_on_error >= 0) {
      printf("Max iteration before finding correct nodes in RadialBasis1 %d %d %d\n",
	     nodes, nr, i2);
    }
    return -3;
  }
  
  e = 0.5*(emin+emax);
  niter = 0;
  while (niter < max_iteration) {
    niter++;
    SetPotentialW(pot, e, orb->kappa, kv);
    SetVEffective(kl, kv, pot);
    i2 = TurningPoints(orb->n, e, orb->kappa, pot, orb->isol);
    if (i2 < 5) {
      nodes = -1;
    } else {
      bqp0 = DpDr(orb->kappa, kv, 0, e, pot, 0, -1, &orb->bqp0);
      nodes = IntegrateRadial(p, e, pot, 0, bqp0, i2, 1.0, 2);
    }
    //printf("nda: %d %d %d %d %d %d %d %15.8E %15.8E %15.8E %15.8E\n", niter, orb->n, orb->kappa, i2, nodes, nr, orb->isol, e, emin, emax, pot->sturm_ene[kv]);
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
    if (_on_error >= 0) {
      printf("Max iteration before finding correct nodes in RadialBasis2 %d %d %d\n",
	     nodes, nr, i2);
    }
    return -4;
  }
  
  int minib = 16;
  de = emax-emin;
  ep = _eneabserr;
  niter = 0;
  while (niter < max_iteration) {
    de *= 0.25;
    nodes = nr;
    while (nodes <= nr+1-mb && niter < max_iteration) {
      niter++;
      e += de;
      SetPotentialW(pot, e, orb->kappa, kv);
      SetVEffective(kl, kv, pot);
      i2 = TurningPoints(orb->n, e, orb->kappa, pot, orb->isol);
      bqp0 = DpDr(orb->kappa, kv, 0, e, pot, 0, -1, &orb->bqp0);
      nodes = IntegrateRadial(p, e, pot, 0, bqp0, i2, 1.0, 2);
      //printf("nd0: %d %d %d %d %d %d %15.8E %15.8E %15.8E\n", niter, orb->n, orb->kappa, i2, nodes, nr, e, de, ep);
    }
    if (nodes - nr == 2-mb && de < ep) break;
    e -= de;
    nodes = nr;    
  }
  if (niter == max_iteration) {
    if (_on_error >= 0) {
      printf("Max iteration before finding correct nodes in RadialBasis5 %d %d %d %g %g %g\n",
	     nodes, nr, i2, de, emin, emax);
    }
    return -5;
  }
  emax = e;
  
  if (nr > 0) {
    e -= de;
    nodes = nr;
    de = emax-emin;
    niter = 0;    
    while (niter < max_iteration) {
      de *= 0.25;
      while (nodes >= nr && niter < max_iteration) {
	niter++;
	e -= de;
	SetPotentialW(pot, e, orb->kappa, kv);
	SetVEffective(kl, kv, pot);
	i2 = TurningPoints(orb->n, e, orb->kappa, pot, orb->isol);
	if (i2 < 5) {
	  nodes = -2;
	} else {
	  bqp0 = DpDr(orb->kappa, kv, 0, e, pot, 0, -1, &orb->bqp0);
	  nodes = IntegrateRadial(p, e, pot, 0, bqp0, i2, 1.0, 2);
	}
	//printf("nd1: %d %d %d %d %d %d %15.8E %15.8E %15.8E\n", niter, orb->n, orb->kappa, i2, nodes, nr, e, de, ep);
      }
      if (nodes - nr == -1 && de < ep) break;
      e += de;
      nodes = nr;
    }
    if (niter == max_iteration) {
      if (_on_error >= 0) {
	printf("Max iteration before finding correct nodes in RadialBasis6 %d %d %d\n",
	       nodes, nr, i2);
      }
      return -6;
    }    
    emin = e;
  }

  bqp0 = 0.0;
  bqp1 = 0.0;
  niter = 0;
  i2 = 0;
  i2o = 0;
  while (niter < max_iteration) {
    niter++;
    e = 0.5*(emin+emax);
    SetPotentialW(pot, e, orb->kappa, kv);
    SetVEffective(kl, kv, pot);
    i2o = TurningPoints(orb->n, e, orb->kappa, pot, orb->isol);
    if (i2o < 5) {
      emin = e;
      continue;
    }
    bqp0 = DpDr(orb->kappa, kv, 0, e, pot, 0, -1, &orb->bqp0);
    nodes = IntegrateRadial(p, e, pot, 0, bqp0, i2o, 1.0, 2);

    if (nodes > nr+1) {
      emax = e;
      continue;
    } else if (nodes < nr) {
      emin = e;
      continue;
    }
    i2 = i2o;
    if (i2 > ib-minib) i2 = ib-minib;    
    i2 = LastMaximum(p, 0, i2, pot);
    if (i2 > i2o-2) {
      i2 = i2o-2;
    }
    i2m1 = i2 - 1;
    i2m2 = i2 - 2;
    i2p1 = i2 + 1;
    i2p2 = i2 + 2;
    p1 = p[i2m2];
    p0 = p[i2m1];
    p2 = p[i2];    
    qo = (-4.0*p[i2m2-1] + 30.0*p[i2m2] - 120.0*p[i2m1]
	  + 40.0*p[i2] + 60.0*p[i2p1] - 6.0*p[i2p2])/120.0;
    qo /= p2;
    if (mb == 0) {
      IntegrateRadial(p, e, pot, i2m2, 1.0, ib, 0.0, 0);
    } else {
      bqp1 = DpDr(orb->kappa, kv, ib, e, pot, bqp, 1, &orb->bqp1);
      IntegrateRadial(p, e, pot, i2m2, 1.0, ib, bqp1, 1);
    }
    qi = (6.0*p[i2m2] - 60.0*p[i2m1] - 40.0*p[i2] + 120.0*p[i2p1]
	  - 30.0*p[i2p2] + 4.0*p[i2p2+1])/120.0;
    qi /= p[i2];
    delta = qo - qi;      
    fact = p[i2]/p2;
    p[i2m1] = p0;
    p[i2m2] = p1;
    for (i = 0; i < i2; i++) p[i] *= fact;
    nodes = CountNodes(p, pot, 0, ib);
    //printf("rb: %d %d %d %d %d %d %d %d %d %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E\n", niter, orb->n, orb->kappa, mb, i2, i2o, ib, nodes, nr, pot->rad[i2], e, _veff[i2], delta, emin, emax, qo, qi, pbqp, bqp, bqp0, bqp1);    
    
    if (nodes > nr) {
      emax = e;
      continue;
    } else if (nodes < nr) {
      emin = e;
      continue;
    }
    if (delta < 0) {
      emax = e;
    } else if (delta > 0) {
      emin = e;
    } else {
      break;
    }
    if (fabs(emax-emin) < EneTol(e)*0.1) break;
  }

  orb->energy = e;
  if (niter == max_iteration) {
    if (_on_error >= 0) {
      printf("Max iteration before finding solution in RadialBasis %d %d %d %d\n",
	     nodes, nr, i2, i2o);
    }
    return -7;
  }
  for (i = 0; i <= ib; i++) {
    p[i] *= pot->dr_drho2[i];
  }
  norm2 = InnerProduct(0, ib, p, p, pot);
  if (fabs(delta) > _matchtol) {
    if (_on_error >= 0) {
     printf("Discontinuous Match in RadialBasis: %d %d %d %d %d %g %g %g\n",
	    niter, orb->n, orb->kappa, i2, i2o, e, delta, norm2);
    }
    return -8;
  }
  fact = 1.0/sqrt(norm2);
  for (i = 1; i <= ib; i++) {
    if (fabs(p[i]*fact) < wave_zero) continue;
    if (p[i-1] > 0 && p[i] < p[i-1]) break;
    else if (p[i-1] < 0 && p[i] > p[i-1]) {
      fact = -fact;
      break;
    }
  }
  if (i > ib && p[ib] < 0) fact = -fact;
  for (i = 0; i <= ib; i++) {
    p[i] *= fact;
  }
  for (i = ib+1; i < pot->maxrp; i++) {
    p[i] = 0.0;
  }  
  orb->ilast = ib;
  orb->qr_norm = 1.0;
  nodes = CountNodes(p, pot, 0, ib);
  if (nodes != nr) {
    if (_on_error >= 0) {
      printf("RadialBasis: No. nodes changed in iteration: %d %d %d %d %d %d %g\n",
	     niter, i2, nodes, nr, orb->n, orb->kappa, orb->energy);
    }
    return -9;
  }    
  if (pot->flag == -1) {
    DiracSmall(orb, pot, ib, kv);
  }
  return 0;
}

int RadialBasis(ORBITAL *orb, POTENTIAL *pot) {
  if (orb->isol == -1) {
    return RadialBasisBQP(orb, pot, 1E30);
  }
  if (pot->bqp > -1E12+EPS5) {
    return RadialBasisBQP(orb, pot, pot->bqp);
  }
  //the mode with bqp<=-1E12 may fail, we revert to 0 boundary when that happens
  int on_error = _on_error;
  _on_error = -1;
  int r = RadialBasisBQP(orb, pot, pot->bqp);
  if (r < -1) {
    free(orb->wfun);
    r = RadialBasisBQP(orb, pot, 1E30);
  }
  _on_error = on_error;
  return r;
}

int RadialSturm(ORBITAL *orb, POTENTIAL *pot) {
  int ierr, k, niter;
  double e0, etol, de;

  e0 = orb->energy;
  ierr = RadialBasis(orb, pot);
  if (ierr < 0) return ierr;
  k = orb->kv;
  orb->isol = -1;
  etol = EneTol(e0)*0.1;
  pot->sturm_ene[k] = 0.0;
  for (niter = 0; niter <= max_iteration; niter++) {
    de = orb->energy-e0;
    //printf("RadialSturm0: %d %d %d %d %12.5E %12.5E %12.5E\n", niter, orb->n, orb->kappa, k, e0, orb->energy, de);
    if (fabs(de) <= etol) break;
    AddPotentialSturm(pot, k, de);
    ierr = RadialBasis(orb, pot);
    if (ierr < 0) return ierr;
  }
  AddPotentialSturm(pot, k, -pot->sturm_ene[k]);
  if (niter > max_iteration) return -niter;
  return 0;
}

int RadialBound(ORBITAL *orb, POTENTIAL *pot) {
  double z, z0, e, de, ep, delta, emin, emax;
  double *p, norm2, fact, p0, p1, p2, qo, qi, bqp;
  double qoa[5], qia[5];
  int i, kl, nr, nodes, niter;
  int i2, i2o, i2m1, i2m2, i2p1, i2p2;
  
  int kv = 0;
  if (pot->pse) kv = IdxVT(orb->kappa);
  orb->kv = kv;
  kl = orb->kappa;
  kl = (kl < 0)? (-kl-1):kl;

  if (kl < 0 || kl >= orb->n) {
    if (_on_error >= 0) {
      printf("Invalid orbital angular momentum, L=%d, %d %d\n", 
	     kl, orb->n, orb->kappa);
    }
    return -1;
  }

  p = malloc(sizeof(double)*2*pot->maxrp);
  if (!p) return -1;
  orb->wfun = p;
  
  nr = orb->n - kl - 1;
  z0 = pot->atom->atomic_number;
  z = GetResidualZ();
  z = 0.5*(z0+z);

  double emin0 = 1.1*EnergyH(z0, 1, -1);
  ep = EnergyH(z, orb->n, orb->kappa);
  e = ep;
  de = fabs(ep)*0.5;
  delta = fabs(emin0)*0.05;
  de = Min(de, delta);
  niter = 0;
  while (niter < max_iteration) {
    niter++;
    SetPotentialW(pot, e, orb->kappa, kv);
    SetVEffective(kl, kv, pot);
    i2 = TurningPoints(orb->n, e, orb->kappa, pot, 0);
    if (i2 > 0) {
      nodes = IntegrateRadial(p, e, pot, 0, 0.0, i2, 1.0, 0);
      if (nodes > nr) break;
    }
    //printf("ndm0: %d %d %d %d %d %d %15.8E %15.8E\n", niter, orb->n, orb->kappa, i2, nodes, nr, e, de);
    e += de;
    delta = 0.5*fabs(e);
    de = Min(de, delta);
  }
  if (niter == max_iteration) {
    if (_on_error >= 0) {
      printf("Max iteration before finding correct nodes in RadialBound a: %d %d %d %d %d %g\n",
	     nodes, nr, i2, orb->n, orb->kappa, e);
    }
    //free(p);
    return -2;
  }
  emax = e;

  e = ep;
  niter = 0;
  while (niter < max_iteration) {
    niter++;
    SetPotentialW(pot, e, orb->kappa, kv);
    SetVEffective(kl, kv, pot);
    i2 = TurningPoints(orb->n, e, orb->kappa, pot, 0);
    if (i2 < 0) {
      nodes = -1;
    } else {
      bqp = DpDr(orb->kappa, kv, 0, e, pot, 0, -1, &orb->bqp0);
      nodes = IntegrateRadial(p, e, pot, 0, bqp, i2, 1.0, 2);
    }
    if (nodes < nr) break;
    e = e*2.0;
    if (e <= emin0) break;
  }
  emin = e;
  if (niter == max_iteration) {
    if (_on_error >= 0) {
      printf("Max iteration before finding correct nodes in RadialBound c: %d %d %d %d %g %g %g\n",
	     orb->n, orb->kappa, nodes, nr, e, emin, emax);
    }
    //free(p);
    return -3;
  }
  
  niter = 0;
  while (niter < max_iteration) {
    niter++;
    e = 0.5*(emin+emax);
    SetPotentialW(pot, e, orb->kappa, kv);
    SetVEffective(kl, kv, pot);
    i2 = TurningPoints(orb->n, e, orb->kappa, pot, orb->isol);
    if (i2 < 0) {
      nodes = -1;
    } else {
      bqp = DpDr(orb->kappa, kv, 0, e, pot, 0, -1, &orb->bqp0);
      nodes = IntegrateRadial(p, e, pot, 0, bqp, i2, 1.0, 2);
    }
    //printf("nda: %d %d %d %d %d %d %d %15.8E %15.8E %15.8E %15.8E\n", niter, orb->n, orb->kappa, i2, nodes, nr, orb->isol, e, emin, emax, pot->sturm_ene[kv]);
    if (nodes == nr) break;
    else if (nodes > nr) {
      emax = e;
    } else {
      emin = e;
    }
  }
  if (niter == max_iteration) {
    if (_on_error >= 0) {
      printf("Max iteration before finding correct nodes in RadialBound d: %d %d\n",
	     nodes, nr);
    }
    return -4;
  }
  
  de = emax-emin;
  ep = _eneabserr;
  niter = 0;
  while (niter < max_iteration) {
    de *= 0.25;
    nodes = nr;
    while (nodes == nr && niter < max_iteration) {
      niter++;
      e += de;
      SetPotentialW(pot, e, orb->kappa, kv);
      SetVEffective(kl, kv, pot);
      i2 = TurningPoints(orb->n, e, orb->kappa, pot, orb->isol);
      bqp = DpDr(orb->kappa, kv, 0, e, pot, 0, -1, &orb->bqp0);
      nodes = IntegrateRadial(p, e, pot, 0, bqp, i2, 1.0, 2);
      //printf("nd0: %d %d %d %d %d %d %15.8E %15.8E %15.8E %15.8E\n", niter, orb->n, orb->kappa, i2, nodes, nr, e, de, emin, emax);
    }
    if (nodes - nr == 1 && de < ep) break;
    e -= de;
    nodes = nr;    
  }
  if (niter == max_iteration) {
    if (_on_error >= 0) {
      printf("Max iteration before finding correct nodes in RadialBound5 %d %d %d\n",
	     nodes, nr, i2);
    }
    return -5;
  }
  emax = e;
  
  if (nr > 0) {
    e -= de;
    nodes = nr;
    de = emax-emin;
    niter = 0;    
    while (niter < max_iteration) {
      de *= 0.25;
      while (nodes == nr && niter < max_iteration) {
	niter++;
	e -= de;
	SetPotentialW(pot, e, orb->kappa, kv);
	SetVEffective(kl, kv, pot);
	i2 = TurningPoints(orb->n, e, orb->kappa, pot, orb->isol);
	bqp = DpDr(orb->kappa, kv, 0, e, pot, 0, -1, &orb->bqp0);
	nodes = IntegrateRadial(p, e, pot, 0, bqp, i2, 1.0, 2);
	//printf("nd1: %d %d %d %d %d %d %15.8E %15.8E %15.8E\n", niter, orb->n, orb->kappa, i2, nodes, nr, e, de, bqp);
      }
      if (nodes - nr == -1 && de < ep) break;
      e += de;
      nodes = nr;
    }
    if (niter == max_iteration) {
      if (_on_error >= 0) {
	printf("Max iteration before finding correct nodes in RadialBasis6 %d %d %d\n",
	       nodes, nr, i2);
      }
      return -6;
    }    
    emin = e;
  }

  niter = 0;
  delta = 1e10;
  int minib = 16;
  while (niter < max_iteration) {
    niter++;
    e = 0.5*(emin + emax);
    SetPotentialW(pot, e, orb->kappa, kv);
    SetVEffective(kl, kv, pot);
    i2o = TurningPoints(orb->n, e, orb->kappa, pot, 0);
    if (i2o < 5) {
      emin = e;
      continue;
    }
    bqp = DpDr(orb->kappa, kv, 0, e, pot, 0, -1, &orb->bqp0);
    nodes = IntegrateRadial(p, e, pot, 0, bqp, i2o, 1.0, 2);  
    if (nodes > nr) {
      emax = e;
      continue;
    } else if (nodes < nr) {
      emin = e;
      continue;
    }
    i2 = i2o;
    if (i2 > pot->maxrp-minib) i2 = pot->maxrp-minib;
    i2 = LastMaximum(p, 0, i2, pot);
    if (i2 > i2o-2) i2 = i2o-2;
    i2m1 = i2 - 1;
    i2m2 = i2 - 2;
    i2p1 = i2 + 1;
    i2p2 = i2 + 2;
    p1 = p[i2m2];
    p0 = p[i2m1];
    p2 = p[i2];    
    qo = (-4.0*p[i2m2-1] + 30.0*p[i2m2] - 120.0*p[i2m1]
	  + 40.0*p[i2] + 60.0*p[i2p1] - 6.0*p[i2p2])/120.0;
    qo /= p2;
    IntegrateRadial(p, e, pot, i2m2, 1.0, pot->maxrp-1, 0.0, 0);
    qi = (6.0*p[i2m2] - 60.0*p[i2m1] - 40.0*p[i2] + 120.0*p[i2p1]
	  - 30.0*p[i2p2] + 4.0*p[i2p2+1])/120.0;
    qi /= p[i2];
    delta = qo - qi;
    fact = p[i2]/p2;
    p[i2m1] = p0;
    p[i2m2] = p1;
    for (i = 0; i < i2; i++) p[i] *= fact;
    nodes = CountNodes(p, pot, 0, pot->maxrp-1);

    //printf("rb: %d %d %d %d %d %d %d %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E\n", niter, orb->n, orb->kappa, i2, i2o, nodes, nr, pot->rad[i2], e, de, delta, emin, emax, qo, qi); 
  
    if (nodes > nr) {
      emax = e;
      continue;
    } else if (nodes < nr) {
      emin = e;
      continue;
    }
    
    if (delta < 0) {
      emax = e;
    } else if (delta > 0) {
      emin = e;
    } else {
      break;
    }
    if (fabs(emax-emin) < EneTol(e)*0.1) break;
  }

  orb->energy = e;
  if (niter == max_iteration) {
    if (_on_error >= 0) {
      printf("Max iteration reached in RadialBound: %d %d\n",
	     orb->n, orb->kappa);
    }
    return -7;
  }
	     
  for (i = 0; i < pot->maxrp; i++) {
    p[i] *= pot->dr_drho2[i];
  }
  norm2 = InnerProduct(0, pot->maxrp-1, p, p, pot);
  if (fabs(delta) > _matchtol) {
    if (_on_error >= 0) {
      printf("Discontinuous Match in RadialBound: %d %d %g %g\n",
	     orb->n, orb->kappa, delta, norm2);
    }
    return -8;
  }
  if (isinf(norm2)) {
    if (_on_error >= 0) {
      printf("RadialBound: inf norm: %d %d %d %d %d %d %g %g %g %g %g\n", niter, orb->n, orb->kappa, orb->ilast, i2o, i2, orb->energy, norm2, delta, p2, p[i2]);
    }
    return -9;
  }
  if (orb->energy >= 0) {
    if (_on_error >= 0) {
      printf("RadialBound: positive energy for bound %d %d %g\n", orb->n, orb->kappa, orb->energy);
    }
    return -10;
  }
  
  fact = 1.0/sqrt(norm2);
  for (i = 1; i < pot->maxrp; i++) {
    if (fabs(p[i]*fact) < wave_zero) continue;
    if (p[i-1] > 0 && p[i] < p[i-1]) break;
    else if (p[i-1] < 0 && p[i] > p[i-1]) {
      fact = -fact;
      break;
    }
  }
  for (i = 0; i < pot->maxrp; i++) {    
    p[i] *= fact;
  }

  for (i = pot->maxrp-1; i >= 0; i--) {
    if (fabs(p[i]) > wave_zero) break;
  }
  if (IsEven(i)) i++;
  orb->ilast = i;
  
  nodes = CountNodes(p, pot, 0, orb->ilast);
  if (nodes != nr) {
    if (_on_error >= 0) {
      printf("RadialBound: No. nodes changed in iteration: %d %d %d %d %d %d %g\n",
	     niter, i2, nodes, nr, orb->n, orb->kappa, orb->energy);
    }
    return -10;
  }    

  orb->qr_norm = 1.0;
  if (pot->flag == -1) {
    DiracSmall(orb, pot, -1, kv);
  }

  return 0;
}

int RadialRydberg(ORBITAL *orb, POTENTIAL *pot) {
  double z, e, e0;
  int i, j, kl, niter, ierr;
  double x0, pp, qq, ppi, qqi;
  int i2, i2p, i2m, i2p2, i2m2, nodes, nr;
  double qo, qi, norm2, delta, zp, *p;
  double dn, dn0, dd, ep, p1, p2, fact, bqp;

  int kv = 0;
  if (pot->pse) kv = IdxVT(orb->kappa);
  orb->kv = kv;
  i = pot->maxrp-1;
  z = -(pot->VT[kv][i]*pot->rad[i]);
  if (z < 1) z = 1;
  kl = orb->kappa;
  kl = (kl < 0)? (-kl-1):kl;
  if (kl < 0 || kl >= orb->n) {
    if (_on_error >= 0) {
      printf("Invalid orbital angular momentum, L=%d, %d %d\n", 
	     kl, orb->n, orb->kappa);
    }
    return -1;
  }
  e = EnergyH(z, orb->n, orb->kappa);
  SetPotentialW(pot, e, orb->kappa, kv);
  SetVEffective(kl, kv, pot);
  i2 = TurningPoints(orb->n, e, orb->kappa, pot, 0);
  if (i2 < pot->maxrp - 3*pot->asymp) {
    return RadialBound(orb, pot);
  }
  p = malloc(sizeof(double)*2*pot->maxrp);
  if (!p) return -1;
  orb->wfun = p;
  nr = orb->n - kl - 1;
  if (pot->N <= 1) j = 3;
  else j = 1;
  i2p2 = pot->maxrp - j*pot->asymp - 5;
  qo = e-_veff[i2p2];
  for (i2 = i2p2-2; i2 > pot->r_core; i2--) {
    p1 = pot->VT[kv][i2]*pot->rad[i2];
    if (fabs(p1+z) > EPS3) break;
    qi = e-_veff[i2];
    if (qi-qo > 10*qo) break;
    if (_veff[i2-1] >= _veff[i2]) break;
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
  i2 = LastMaximum(p, pot->r_core, i2, pot);
  i2p = i2 + 1;
  i2m = i2 - 1;
  i2p2 = i2 + 2;
  i2m2 = i2 - 2;

  dn = 0.05;
  dn0 = orb->n+dn;
  dd = 0.0;
  j = 0;
  while (1) {
    e = EnergyH(z, dn0, orb->kappa);
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
    //printf("dn: %d %d %d %d %g %15.8E %15.8E %g %g %g %g %g %g %g\n", orb->n, orb->kappa, i2, ierr, z, e, dn0, dn, qq, pp, qo, qi, delta, dd);
    if (delta > 0) {
      if (j > 0) {
	if (dn < EPS5) {
	  if (delta > dd) {
	    p1 = delta/(delta-dd);	  
	    dn0 = dn0*(1-p1) + (dn0+dn)*p1;
	  } else {
	    dn0 += 0.5*dn;
	  }
	  break;
	} else {
	  dn /= 2;
	  dn0 += dn;
	  continue;	
	}
      }
    } else {
      j = 1;
    }
    dd = delta;
    dn0 -= dn;
  }
  e = EnergyH(z, dn0, orb->kappa);
  i2p2 = pot->maxrp-1;
  bqp = DpDr(orb->kappa, kv, 0, e, pot, 0, -1, &orb->bqp0);
  nodes = IntegrateRadial(p, e, pot, 0, bqp, i2p2, 1.0, 2);
  for (i = 0; i <= i2p2; i++) {
    p[i] = p[i] * pot->dr_drho2[i];
  }
  if (_rydnorm <= 0) {
    i2 = LastMaximum(p, pot->r_core, i2, pot);
    zp = FINE_STRUCTURE_CONST2*e;
    x0 = pot->rad[i2];
    ierr = 1;
    DCOUL(z, e, orb->kappa, x0, &pp, &qq, &ppi, &qqi, &ierr);
    norm2 = pp;
    fact = fabs(norm2/p[i2]);
    i2 = FirstMaximum(p, 0, i2, pot);
    if (p[i2] < 0) {
      fact = -fact;
    }
    for (i = 0; i <= i2p2; i++) {
      p[i] *= fact;
    }
  }
  i = pot->maxrp-1;
  orb->ilast = i;
  orb->energy = e;
      
  e0 = InnerProduct(0, pot->maxrp-1, p, p, pot);
  orb->qr_norm = 1.0/e0;

  if (pot->flag == -1) {
    DiracSmall(orb, pot, -1, kv);
    if (_rydnorm > 0) {
      for (i = 0; i <= pot->ifm; i++) {
	pp = orb->wfun[i];
	qq = orb->wfun[i+pot->maxrp];
	pot->dW2[i] = (pp*pp + qq*qq)*pot->dr_drho[i];
      }
      p2 = pot->nmax;
      p2 /= orb->n;
      p2 = pow(p2, 3-p2*(3.0-pot->xfm));
      p2 /= Simpson(pot->dW2, 0, pot->ifm);
      fact = sqrt(pot->cfm*p2);
      i2 = FirstMaximum(p, 0, i2, pot);
      if (p[i2] < 0) fact = -fact;
      for (i = 0; i < pot->maxrp; i++) {
	orb->wfun[i] *= fact;
	orb->wfun[i+pot->maxrp] *= fact;
      }
    }
  }

  return 0;
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
    if (_on_error >= 0) {
      printf("Kappa == 0 in RadialFreeInner\n");
    }
    return -1;
  }
  if (pot->ib1 <= 0 && pot->ib <= 0) {
    if (_on_error >= 0) {
      printf("Boundary not set\n");
    }
    return -2;
  }
  SetPotentialW(pot, e, kl, kv);
  kl = (kl < 0)? (-kl-1):kl;  
  p = malloc(2*pot->maxrp*sizeof(double));
  if (!p) return -1;
  orb->wfun = p;
  SetVEffective(kl, kv, pot);

  i2 = pot->ib1;
  i2 = Max(i2,pot->ib) + 2;      
  bqp0 = DpDr(orb->kappa, kv, 0, e, pot, 0, -1, &orb->bqp0);
  nodes = IntegrateRadial(p, e, pot, 0, bqp0, i2, 1.0, 2);
  for (i = i2; i >= 0; i--) {
    p[i] *= pot->dr_drho2[i];
  }
  orb->ilast = i2;
  //orb->wfun = p;
  orb->phase = NULL;
  
  orb->qr_norm = 1.0;
  
  if (pot->flag == -1) {
    DiracSmall(orb, pot, i2, kv);
  }

  return 0;
}

/* note that the free states are normalized to have asymptotic 
   amplitude of 1/sqrt(k), */
// Solve the Dirac equation for a free (continuum) electron 
// wave function(energy > 0)
int RadialFree(ORBITAL *orb, POTENTIAL *pot) {
  int i, kl, nodes;
  int i2, i2p, i2m, i2p2, i2m2;
  double *p, po, qo, e, po1, bqp0;
  double dfact, da, cs, si, phase0;

  if (_zcoll || _mcoll) return RadialFreeZMC(orb, pot);
  
  int kv = 0;
  if (_pwa > 0) {
    kv = -1;
  } else {
    if (pot->pse) kv = IdxVT(orb->kappa);
    orb->kv = kv;
  }
  e = orb->energy;
  //Catch erroneous energy inputs
  if (e < 0.0) {
    if (_on_error >= 0) {
      printf("Energy < 0 in Free\n");
    }
    return -1;
  }
  //Catch erroneous angular momentum quantum number inputs
  kl = orb->kappa;
  if (orb->kappa == 0) {
    if (_on_error >= 0) {
      printf("Kappa == 0 in Free\n");
    }
    return -1;
  }
  // Initialize effective W potential with correct constants 
  // etc
  SetPotentialW(pot, e, kl, kv);
  kl = (kl < 0)? (-kl-1):kl;  
  p = malloc(2*pot->maxrp*sizeof(double));
  if (!p) return -1;
  orb->wfun = p;
  // Calculates the effective potential for the
  // transformed radial Dirac equation Eq. (8) of
  // Gu (2002), Veff(r) = U + k(k+1)/r^2/2
  // effective potential is stored in pot
  SetVEffective(kl, kv, pot);
  
  // Determine cut-off points between region I
  // (Numerov solution to radial Dirac equation)
  // and region II (Hullac style Phase-Amplitude
  // solution)
  i2 = TurningPoints(0, e, orb->kappa, pot, 0);
  i2m = i2 - 1;
  i2p = i2 + 1;
  i2m2 = i2 - 2;
  i2p2 = i2 + 2;
      
  bqp0 = DpDr(orb->kappa, kv, 0, e, pot, 0, -1, &orb->bqp0);
  // Solution to the second order radial equation via the Numerov
  // method, between 0 and i2p2, solution stored in p
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
  // Solution to the Radial Dirac equation using the phase-
  // amplitude method of Bar-Shalom et al. Computer Physics
  // Communications 83 (1996) 21-32
  // First calculate amplitude solution for all r>i2
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
  // Next calculate phase solution for all r>i2.
  // Amplitude and phase will be stored in alternating
  // (ie odd and even) entries of the array p
  Phase(p, pot, i2, phase0);
  
  p[i2] = po;
  p[i2p] = po1;
  // Adjust amplitude of region I solution to match
  // that in region II
  for (i = 0; i < i2p2; i++) {
    p[i] *= dfact;
  }
    
  orb->ilast = i2p;
  //orb->wfun = p;
  orb->phase = NULL;

  orb->qr_norm = 1.0;
  
  if (pot->flag == -1) {
    DiracSmall(orb, pot, -1, kv);
  }

  return 0;
}

void DiracODE(int *neq, double *t, double *y, double *yd) {
  double v, r, a, a2, kr;

  a = FINE_STRUCTURE_CONST;
  a2 = 2/FINE_STRUCTURE_CONST2;
  if (y[7] > 0) a2 *= y[7];
  r = *t;
  v = (y[2] + y[3]*(r - y[4]))/r;
  kr = y[6]/r;
  yd[0] = a*(y[5] - v + a2)*y[1] - kr*y[0];
  yd[1] = a*(-y[5] + v)*y[0] + kr*y[1];
}
  
int RadialFreeZMC(ORBITAL *orb, POTENTIAL *pot) {
  double a, a2, b, z0, e, r, r0;
  int i, j, kl, i2, i2p, i2m, i2p2, i2m2, nodes;
  double *p, *q, po, qo, po1;
  double dfact, da, cs, si, phase0;
  int neq, itol, itask, istate, iopt, lrw, iwork[22], liw, mf;
  double y[8], rtol, atol[2], *rwork;
  double zc = _zcoll;
  double mc = _mcoll;
  int kv = 0;
  if (pot->pse) kv = IdxVT(orb->kappa);
  orb->kv = kv;
  e = orb->energy;
  if (e < 0.0) {
    if (_on_error >= 0) {
      printf("Energy < 0 in Free\n");
    }
    return -1;
  }
  kl = orb->kappa;
  if (orb->kappa == 0) {
    if (_on_error >= 0) {
      printf("Kappa == 0 in Free\n");
    }
    return -1;
  }
  SetPotentialWZMC(pot, e, kl, kv);
  kl = (kl < 0)? (-kl-1):kl;  
  p = malloc(2*pot->maxrp*sizeof(double));
  if (!p) return -1;
  orb->wfun = p;
  SetVEffectiveZMC(kl, kv, pot);
  q = p + pot->maxrp;
  
  i2 = TurningPoints(0, e, orb->kappa, pot, 0);
  i2m = i2 - 1;
  i2p = i2 + 1;
  i2m2 = i2 - 2;
  i2p2 = i2 + 2;
  
  a2 = FINE_STRUCTURE_CONST2;
  if (mc > 0) a2 /= mc;
  if (pot->atom->z1 > 0) {
    z0 = pot->atom->z1;
    if (zc) z0 *= zc;
    if (orb->kappa < 0) {
      b = (e + z0)*FINE_STRUCTURE_CONST;
      b *= -pot->rad[i]/(2*kl+3.0);
    } else {
      b =(e + z0 + 2/a2)*FINE_STRUCTURE_CONST;
      b *= pot->rad[i]/(2*kl+1.0);
      b = 1.0/b;
    }
  } else {
    z0 = pot->Z[0];
    if (zc) z0 *= zc;
    double g = sqrt(orb->kappa*orb->kappa - z0*z0*FINE_STRUCTURE_CONST2);
    b = z0*FINE_STRUCTURE_CONST/(orb->kappa - g);   
  }    
  
  rwork = _dwork2;
  lrw = pot->maxrp;
  liw = 22;
  neq = 2;
  itol = 2;
  rtol = EPS6;
  atol[0] = 0.0;
  atol[1] = EPS8;
  itask = 4;
  istate = 1;
  iopt = 0;
  mf = 10;
  p[0] = 1.0;
  q[0] = b;
  r0 = pot->rad[0];
  for (i = 1; i <= i2p2; i++) {
    r = pot->rad[i];
    rwork[0] = r;
    y[0] = p[i-1];
    y[1] = q[i-1];
    y[2] = pot->VT[kv][i-1]*r0;
    y[3] = (pot->VT[kv][i]*r-y[2])/(r-r0);
    if (zc) {
      y[2] *= zc;
      y[3] *= zc;
    }
    y[4] = r0;
    y[5] = e;
    y[6] = orb->kappa;
    y[7] = mc;
    while (r0 < r) {
      LSODE(DiracODE, &neq, y, &r0, r, itol, rtol, atol,
	    itask, &istate, iopt, rwork, lrw, iwork, liw, NULL, mf);
      if (istate == -1) istate = 2;
      else if (istate < 0) {
	printf("RadialFreeZMC LSODE Error: %d\n", istate);
	Abort(1);
      }
    }
    p[i] = y[0];
    q[i] = y[1];
    
    a = sqrt(p[i]*p[i] + q[i]*q[i]);
    if (a > 1E10) {
      for (j = 0; j <= i; j++) {
	if (fabs(p[j]) < 1e-250 || fabs(q[j]) < 1e-250) {
	  p[j] = 0;
	  q[j] = 0;
	} else {
	  break;
	}
      }      
      for (j = 0; j <= i; j++) {
	p[j] /= a;
	q[j] /= a;
      }
      istate = 1;
    }

    r0 = r;
  }

  nodes = CountNodes(p, pot, 0, i2p2);
  
  for (i = i2m2-1; i <= i2p2; i++) {
    _dwork3[i] = e - zc*pot->VT[kv][i];
    ABAND[i] = sqrt(1 + _dwork3[i]*a2*0.5);
    _dwork1[i] = p[i]/ABAND[i];
  }
  dfact = 1.0/pot->dr_drho[i2];
  qo = (-4.0*_dwork1[i2m2-1] + 30.0*_dwork1[i2m2] - 120.0*_dwork1[i2m]
	+ 40.0*_dwork1[i2] + 60.0*_dwork1[i2p] - 6.0*_dwork1[i2p2])/120.0;
  qo *= dfact;
  po = _dwork1[i2];
  po1 = _dwork1[i2p];

  da = Amplitude(p, e, orb->kappa, pot, i2);  
  cs = p[i2] * qo - da*po;
  si = po / p[i2];
  dfact = (si*si + cs*cs);
  dfact = 1.0/sqrt(dfact);
  phase0 = atan2(si, cs);
  if (phase0 < 0) phase0 += TWO_PI;
  Phase(p, pot, i2, phase0);
  
  p[i2] = po*ABAND[i2];
  p[i2p] = po1*ABAND[i2p];  
  for (i = 0; i < i2p2; i++) {
    p[i] *= dfact;
    q[i] *= dfact;
  }

  orb->ilast = i2p;
  //orb->wfun = p;
  orb->phase = NULL;
  orb->qr_norm = 1.0;
  if (pot->flag == -1) {
    DiracSmallZMC(orb, pot, i2p2, kv);
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
  if (orb->n > 0 && pot->mps >= 3) {
    if (orb->ilast > pot->ips) orb->ilast = pot->ips;
  }
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

int DiracSmallZMC(ORBITAL *orb, POTENTIAL *pot, int i1, int kv) {
  double *p, *q, xi, e, a, b, a2, zc;
  int i;

  e = orb->energy;
  p = orb->wfun;
  q = p + pot->maxrp;
  a2 = FINE_STRUCTURE_CONST2;
  zc = _zcoll;
  if (!zc) zc = 1.0;
  if (_mcoll > 0) a2 /= _mcoll;
  for (i = i1; i < pot->maxrp; i += 2) {
    if (kv >= 0) {
      xi = e - zc*pot->VT[kv][i];
    } else {
      xi = e;
    }
    xi = xi*a2*0.5;
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
    b += _dwork1[i]*orb->kappa/pot->rad[i];
    a = _dwork1[i]/(p[i]*p[i]);
    p[i] = _dwork1[i];
    q[i] = FINE_STRUCTURE_CONST*a/(2.0*_dwork[i]);
    q[i+1] = FINE_STRUCTURE_CONST*b/(2.0*_dwork[i]);
    if (_mcoll > 0) {
      q[i] /= _mcoll;
      q[i+1] /= _mcoll;
    }
  }

  b = a2*orb->energy;
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
  if (_mcoll > 0) w *= _mcoll;
  ydot[0] = y[1];
  ydot[1] = 1.0/(y[0]*y[0]*y[0]) - y[0]*w;
}

double Amplitude(double *p, double e, int ka, POTENTIAL *pot, int i0) {
  // Variables
  // p:   Electron wave function
  // e:   Electron energy
  // ka:
  // pot: Electron potential
  // i0:  Array index at which Phase-Amplitude solution to radial
  //      Schrodinger equation begins
  int i, ii, n;
  double a, b, xi, r2, r3, kl1;
  double z, dk, r0, r1, r, w, v1;
  double zc, mc, me, a2;
  int neq, itol, itask, istate, iopt, lrw, iwork[22], liw, mf;
  double y[6], rtol, atol[2], *rwork;

  zc = _zcoll;
  mc = _mcoll;
  n = pot->maxrp-1;
  if (_pwa > 0) {
    z = 0.0;
  } else {
    z = GetResidualZ()-pot->ZPS[n];
  }
  me = e;
  a2 = FINE_STRUCTURE_CONST2;
  if (zc) z *= zc;
  kl1 = ka*(ka+1.0);
  double mkl1 = kl1;
  double az = fabs(z);
  double emin = _emin_amp*e;
  if (mc) {
    me *= mc;
    a2 /= mc;
    az *= mc;
    mkl1 /= mc;
  }
  r1 = pot->rad[n];
  _dwork[0] = r1;
  _dwork1[0] = _veff0[n];
  dk = sqrt(2.0*me*(1.0+0.5*a2*e));
  dk = EPS5*e*dk;
  for (i = 1; i < pot->maxrp; i++) {
    r = _dwork[i-1]*1.05;
    _dwork[i] = r;
    r2 = r*r;
    r3 = r2*r;
    a = -z/r;
    b = e - a;
    xi = 1.0 + 0.5*a2*b;
    //xi = xi*xi;
    b = b*b;
    v1 = z/r2;
    w = (-2.0*v1/r + 0.75*a2*v1*v1/xi - 2*ka*v1/r);
    if (mc > 0) {
      w /= mc;
    }
    w /= 4.0*xi;
    _dwork1[i] = a + 0.5*mkl1/r2 - 0.5*a2*(b - w);
    a = e-_dwork1[i];
    if (a < emin) _dwork1[i] = e-emin;
    a = az/r2;
    b = mkl1/r3;
    if (a < dk && b < dk) break;
  }
  r0 = r;
  w = 2*(e - _dwork1[i]);
  if (mc > 0) w *= mc;
  y[0] = pow(w, -0.25);
  a = a2*e;
  b = a2*z*z;
  y[1] = 0.5*(y[0]/w)*(z*(1+a)/r2-(mkl1-b)/r3);

  rwork = _dwork2;
  lrw = pot->maxrp;
  liw = 22;
  neq = 2;
  itol = 2;
  rtol = EPS6;
  atol[0] = 0.0;
  atol[1] = EPS8;
  itask = 4;
  istate = 1;
  iopt = 0;
  mf = 10;
  
  i--;
  ii = i;
  for (; i >= 0; i--) {
    r = _dwork[i];
    rwork[0] = r;
    y[2] = r0;
    y[3] = _dwork1[i+1]*r0;    
    y[4] = (_dwork1[i]*r - y[3])/(r - r0);
    y[5] = e;
    while (r0 != r) {
      LSODE(DerivODE, &neq, y, &r0, r, itol, rtol, atol, 
	    itask, &istate, iopt, rwork, lrw, iwork, liw, NULL, mf);
      if (istate == -1) istate = 2;
      else if (istate < 0) {
	printf("Amplitude0 LSODE Error: %d %d %g %d %d %d %g\n",
	       istate, ka, e, i, ii, i0, r);
	Abort(1);
      }
    }
  }

  itask = 4;
  istate = 1;
  p[n] = y[0];
  for (i = n-1; i >= i0; i--) {
    r = pot->rad[i];
    rwork[0] = r;
    y[2] = r0;
    a = e-_veff0[i+1];
    if (a < emin) {
      y[3] = (e-emin)*r0;
    } else {
      y[3] = _veff0[i+1]*r0;
    }
    a = e-_veff0[i];
    if (a < emin) {
      y[4] = ((e-emin)*r - y[3])/(r-r0);
    } else {
      y[4] = (_veff0[i]*r - y[3])/(r - r0);
    }
    y[5] = e;
    while (r0 != r) {
      LSODE(DerivODE, &neq, y, &r0, r, itol, rtol, atol,
	    itask, &istate, iopt, rwork, lrw, iwork, liw, NULL, mf);
      if (istate == -1) istate = 2;
      else if (istate < 0) {
	printf("Amplitude1 LSODE Error %d\n", istate);
	Abort(1);
      }
    }
    v1 = z;
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
  // Numerical solution of Eq. (19) from A. Bar-Shalom et al./
  // Computer Physics Communications 93 (1996) 21-32
  // Input variables
  // p:      Electron wave function, even values to be overwritten
  //         by phase solution
  // pot:    Potential to solve Dirac equation for
  // i1:     Index of point at which the Phase-Amplitude solution
  //         begins
  // phase0: Value of the phase component, ie. Boundary condition
  //         at i1
  
  int i;

  for (i = i1; i < pot->maxrp; i++) {
    _dwork[i] = 1.0/p[i];
    _dwork[i] *= _dwork[i];
    _dwork[i] *= pot->dr_drho[i];
  }

  i = i1+1;
  //Start with boundary condition
  p[i] = phase0;
  //Integrate from i1, every even value is overwritten
  for (i = i1+3; i < pot->maxrp; i += 2) {
    p[i] = p[i-2] + (_dwork[i-3] + 4.0*_dwork[i-2] + _dwork[i-1])*ONETHIRD;
    /*
    if (isnan(p[i])) {
      printf("nan phase: %d %g\n", i, p[i-2]);
    }
    */
  }
  return 0;
}

int SetVEffectiveZMC(int kl, int kv, POTENTIAL *pot) {
  double kl1;
  int i;
  double r;

  if (!_zcoll) return SetVEffective(kl, kv, pot);
  kl1 = 0.5*kl*(kl+1);
  if (_mcoll > 0) kl1 /= _mcoll;
  int m5 = pot->maxrp-5;
  if (pot->mps >= 3 || _veff_corr == 0) {
    m5 = pot->maxrp;
  }
  for (i = 0; i < pot->maxrp; i++) {
    r = pot->rad[i];
    r *= r;
    _veff0[i] = _zcoll*pot->VT[kv][i];
    if (i < m5) _veff0[i] += kl1/r;
  }
  
  if (kl1 > 0 && _veff_corr > 1) {
    for (i = m5-1; i > 0; i--) {
      if (pot->rad[i] <= pot->atom->rms0) break;
      if (_veff0[i-1] < _veff0[i]) {
	i = Max(i,pot->ips);
	if (i < m5-1) {
	  i++;
	  r = _veff0[i];
	  for (; i < m5; i++) _veff0[i] = r;
	}
	break;
      }
    }
  }
  
  for (i = 0; i < pot->maxrp; i++) {
    _veff0[i] += pot->W[i];
    _veff[i] = _veff0[i] + pot->vtr[i];
  }

  return 0;
}

int SetVEffective(int kl, int kv, POTENTIAL *pot) {
  // Calculates the effective potential for the
  // transformed radial Dirac equation Eq. (8) of
  // Gu (2002), Veff(r) = U + k(k+1)/r^2/2
  double kl1;
  int i;
  double r;

  kl1 = 0.5*kl*(kl+1);
  if (kv < 0) {
    for (i = 0; i < pot->maxrp; i++) {
      r = pot->rad[i];
      r *= r;
      _veff0[i] = kl1/r;
      _veff[i] = _veff0[i] + pot->vtr[i];
    }
    return 0;
  }
  int m5 = pot->maxrp-5;
  if (pot->mps >= 3 || _veff_corr == 0) {
    m5 = pot->maxrp;
  }
  for (i = 0; i < pot->maxrp; i++) {
    r = pot->rad[i];
    r *= r;
    _veff0[i] = pot->VT[kv][i];
    if (i < m5) _veff0[i] += kl1/r;
  }
  
  if (kl1 > 0 && _veff_corr > 1) {
    for (i = m5-1; i > 0; i--) {
      if (pot->rad[i] <= pot->atom->rms0) break;
      if (_veff0[i-1] < _veff0[i]) {
	i = Max(i,pot->ips);
	if (i < m5-1) {
	  i++;
	  r = _veff0[i];
	  for (; i < m5; i++) _veff0[i] = r;
	}
	break;
      }
    }
  }
  
  for (i = 0; i < pot->maxrp; i++) {
    _veff0[i] += pot->W[i];
    _veff[i] = _veff0[i] + pot->vtr[i];
  }

  return 0;
}

static int TurningPoints(int n, double e, int kappa,
			 POTENTIAL *pot, int isol) {
  // For continuum wave function solutions this function finds 
  // the boundary between Region I (numerov integration) and 
  // Region II (Hullac style Amplitude-phase) solutions to the 
  // radial Dirac equation 
  int i, i2, ip;
  double x, a, b, xp;

  if (n == 0) { //Case of a continuum wave function
    // Start from maximum radial coordinate in grid and 
    // iterate downward (cut off at grid-point 20)
    // Calculate energy below potential and break when this
    // goes to 0 or below
    for (i = pot->maxrp-6; i >= 20; i--) {
      x = e - _veff[i];
      if (x <= 0) break; 
    }    
    ip = i+1;
    for (i = ip+1; i < pot->maxrp-5; i++) {
      x = e - _veff[i];
      b = 1.0/pot->rad[i];
      a = 20.0/(pot->qr*pot->ar*pow(b,1-pot->qr) + pot->br*b);
      x = 2.0*x;
      if (_mcoll > 0) x *= _mcoll;
      x = TWO_PI/sqrt(x);
      //printf("tp: %d %d %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E\n", i, ip, pot->rad[i], e, _veff[i], a, x, pot->rad[i]-pot->rad[ip]);
      if (x < a) break;
      if (pot->rad[i]-pot->rad[ip] > x) break;
    }
    i2 = i-2;
    if (IsOdd(i2)) (i2)--;
  } else if (n > 0) {    
    ip = pot->maxrp-1;
    if (pot->ib && n > pot->nb && isol == 0) ip = pot->ib;
    for (i = ip-1; i >= 0; i--) {
      x = e-_veff[i];
      if (x > 0) {
	i++;
	break;
      }
    }
    i2 = i;
    /*
    if (kappa != -1) {
      for (i = i2-1; i >= 0; i--) {
	x = e-_veff[i];
	if (x < 0) break;
      }
      for (; i >= 0; i--) {
	x = e-_veff[i];
	if (x > 0) {
	  break;
	}
      }
      if (i > 0 && x > 0 && pot->rad[i] > pot->atom->rms0) {
	i2 = i+1;
      }
    }
    */
    //printf("tp: %d %d %d %d %d %d %g %g %g %g\n", n, kappa, i, pot->ib, pot->nb, i2, e, _veff[i], _veff[i+1], _veff[pot->maxrp-1]);
  } else {
    for (i = pot->ib; i < pot->ib1; i++) {
      if (e - _veff[i] > 0) {
	if (i > pot->ib) i--;
	break;
      }
    }
    i2 = i+20;
    for (i = pot->ib1-20; i > i2; i--) {
      if (e - _veff[i] > 0) {
	i++;
	break;
      }
    }
    i2 = i;
  }
  return i2;
}

static int IntegrateRadial(double *p, double e, POTENTIAL *pot,
			   int i1, double p1, int i2, double p2, int q) {
  double a, b, r, x, y, z, a1, a2;
  int kl=1, ku=1, nrhs=1;
  int i, info, n, m, j, k;
  int *ipiv = _ipiv;
  
  m = i2 - i1 - 1;
  if (m < 0) return 0;

  // Apply boundary conditions at i1 and i2
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

  // Initialize p to zero for i1<i<i2
  for (i = i1+1; i < i2; i++) {
    p[i] = 0.0;
  }

  // Numerov method to solve radial second order ODE
  // This is achieved by setting up a banded matrix
  // Containing the coefficients of the Numerov solution
  // and then using the Lapack function DGBSV to solve the 
  // resulting matrix equation

  j = 1;
  n = 2*kl + ku + 1;
  k = kl + ku;
  for (i = i1+1; i < i2; i++, j++, k += n) {
    y = pot->dr_drho[i]*pot->dr_drho[i];
    x = (2.0*(_veff[i] - e))*y/12.0;
    a = 1.0 - x;
    b = -2.0*(1.0 + 5.0*x);
    ABAND[k-1] = a;
    ABAND[k] = b;
    ABAND[k+1] = a;
  }

  i = i2;
  y = pot->dr_drho[i]*pot->dr_drho[i];
  x = (2.0*(_veff[i] - e))*y/12.0;
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
  x = (2.0*(_veff[i] - e))*y/12.0;
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

  // Call to LAPACK function DGBSV to solve matrix equation
  DGBSV(m, kl, ku, nrhs, ABAND, n, ipiv, p+i1+1, m, &info);
  if (info) {
    printf("Error in Integrating the radial equation: %d\n", info);
    Abort(1);
  }

  if (q == 1) {
    p[i2] = 2*p2*p[i2-1] + p[i2-2];
  }
  if (q == 2) {
    p[i1] = -2*p1*p[i1+1] + p[i1+2];
  }

  return CountNodes(p, pot, i1, i2);
}

int CountNodes(double *p, POTENTIAL *pot, int i1, int i2) {
  int n, i;
  double p0, a, b, ma;

  ma = 1/sqrt(InnerProduct(i1, i2, p, p, pot));
  n = 0;
  i = i2;
  for (; i >= i1; i--) {
    if (pot->rad[i] < pot->atom->rms0) return n;
    a = fabs(p[i])*ma;
    if (a > wave_zero) break;
  }
  p0 = p[i]*ma;
  for (; i >= i1; i--) {    
    if (pot->rad[i] < pot->atom->rms0) return n;
    b = p[i]*ma;
    a = fabs(b);
    if (a > wave_zero) {
      if ((p0 > 0 && b < 0) ||
	  (p0 < 0 && b > 0)) {
	n++;
	p0 = b;
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
      if (pot->sturm_idx > 0) {
	rmax = pot->rb*(_sturm_rmx/pot->sturm_idx);
      } else {
	rmax = pot->rb*1.001;
      }
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

void InitializePS(POTENTIAL *pot) {
  double a, z0, x0, x1;
  z0 = GetAtomicNumber();
  pot->xps = 1.0;
  pot->gps = 0.0;
  pot->aps = 0.0;
  pot->bps = 0.0;
  pot->nqf = 0.0;
  pot->iqf = 0;
  if (pot->mps >= 0) {
    if (_relativistic_fermi0 < 0) {
      a = pot->tps*FINE_STRUCTURE_CONST2;
      if (a > _relativistic_xfermi) {
	_relativistic_fermi = 1;
      } else {
	_relativistic_fermi = 0;
      }
    } else {
      _relativistic_fermi = _relativistic_fermi0;
    } 
    if (pot->zps <= 0 || (pot->mps == 0 && 1+pot->ups != 1)) {
      pot->zps = (z0-pot->NC);
      if (pot->ihx < 0 && pot->NC > 0.999) pot->zps -= pot->ihx;
    }
    if (pot->ups < 0) {      
      pot->ups = pot->zps;
    }
    if (pot->mps == 1 && _debye_mode >= 2) {
      pot->dps = pow(1.0/(FOUR_PI*pot->nps), 0.25);
    } else {
      pot->dps = sqrt(pot->tps/(FOUR_PI*pot->nps*(1+pot->ups)));
    }
    pot->ewd = 0.0;
    if (pot->nps > 0) {
      if (pot->mps >= 3) {
	pot->rps = pow(3/(FOUR_PI*pot->nps),ONETHIRD);
	pot->nbt = z0 - pot->NC;
	pot->nbs = 0.0;
      } else {
	pot->rps = pow(3*pot->zps/(FOUR_PI*pot->nps),ONETHIRD);
	if (pot->kps <= 1) {
	  pot->nbt = -1.0;
	} else if (pot->kps == 2) {
	  pot->nbt = z0-pot->zps-pot->NC;	  
	} else if (pot->kps == 3) {
	  pot->nbt = 0;
	}
      }
      if (fabs(1-_sc_rsf)>1e-10) pot->rps *= _sc_rsf;
      if (pot->mps >= 3) {
	if (_sc_ewr > 0) {
	  a = TWO_PI/_sc_ewr;
	} else {
	  a = TWO_PI/pot->rps;
	}
	pot->ewd = (_sc_ewf/2.355)*0.5*a*a;
      } else {
	pot->ewd = 0.0;
      }
    } else {
      pot->rps = pot->dps;
    }
    pot->gps = (pot->zps*pot->ups)/(pot->tps*pot->rps);
    if (pot->mps == 2) {
      pot->aps = FermiDegeneracy(pot->nps, pot->tps, &pot->fps);
      PrepFermiRM1(pot->aps, pot->fps, pot->tps);     
    } else if (pot->mps == 0 && pot->tps > 0) {
      pot->aps = FermiDegeneracy(pot->nps, pot->tps, &pot->fps);
      if (pot->ups > 0 && _ionsph_bmode == 0) {
	PrepFermiRM1(pot->aps, pot->fps, pot->tps);     
      }
    }
  }  
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
  if (z < 1) z = 1.0;
  if (pot->flag == 0) pot->flag = -1; 
  rmin = pot->rmin/z0;
  if (rn > 0) {
    a = rn*GRIDRMINN0;
    if (rmin < a) rmin = a;
    a = rn*GRIDRMINN1;
    if (rmin > a) rmin = a;
  }
  InitializePS(pot);
  maxrp = pot->maxrp;
  if (pot->mps >= 3) {
    pot->qr = 1.0;
    gasymp = -pot->rps;
    if (_sc_rbf > 1) {
      gasymp *= _sc_rbf+1e-3;
    }
    maxrp -= 10;
  }  
  MaxRGrid(pot, gasymp, gratio, z, rmin, maxrp, &a, &c, &rmax);
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
  /*
  if (pot->rps > 0 && pot->mps == 0 && pot->ups+1 == 1) {
    rmax = pot->rps;
    d1 = log(rmax/rmin);
    d2 = pow(rmax, pot->qr) - pow(rmin, pot->qr);
    maxrp = (int)(1+a*d2 + b*d1);
    b = (maxrp - 1.0 - (a*d2))/d1;
  }
  */
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
    pot->vtr[i] = 0.5*pot->dr_drho[i]*(0.5*tp3-0.75*pot->dr_drho[i]*tp2*tp2);
  }

  UVIP3C(3, pot->maxrp, pot->rho, pot->rad,
	 pot->a1r, pot->a2r, pot->a3r);
  UVIP3C(3, pot->maxrp, pot->rho, pot->dr_drho,
	 pot->a1dr, pot->a2dr, pot->a3dr);
  
  RGMQED(&a, &b);
  for (i = 0; i < pot->maxrp; i++) {
    pot->mqrho[i] = b*log(pot->rad[i]) + a*pot->rad[i];
  }
  
  if (pot->rps > 0) {
    for (i = 1; i < pot->maxrp; i++) {
      if (pot->rad[i] > pot->rps+0.5*(pot->rad[i]-pot->rad[i-1])) {
	pot->ips = i-1;
	break;
      }
    }
    if (pot->ips <= 0) {
      pot->rps = 0;
    } else {
      if (pot->mps >= 3 ||
	  ((pot->mps == 0 || pot->mps == 2) && pot->kps > 0)) {
	pot->ib = pot->ips;
	pot->ib1 = pot->ips;
	pot->rb = pot->rad[pot->ib];
	pot->nb = 0;
	if (_sc_rbf > 1) {
	  for (i = pot->ips; i < pot->maxrp-1; i++) {
	    if (pot->rad[i] > _sc_rbf*pot->rad[pot->ib]) break;
	  }
	  pot->ib = i;
	  pot->ib1 = i;
	  pot->rb = pot->rad[pot->ib];	  
	}
	pot->bqp = _sc_bqp;
      }
    }
  }
  if (pot->mps > 0 && pot->mps < 3) {
    if (pot->rad[pot->maxrp-1] < 5*Max(pot->dps,pot->rps)) {
      printf("rmax < 5x ionsphere or debye length: %d %g %g %g\n",
	     pot->maxrp, pot->rad[pot->maxrp-1], pot->dps, pot->rps);
    }
  }
  return 0;
}

double GetRFromRho(double rho, double a, double b, double q, double r) {
  double x, r0, r1;
  int i;

  x = a*pow(r, q) + b*log(r);
  if (x > rho) {
    r1 = r;
    r /= 2;
    for (i = 0; i < 128; i++) {
      x = a*pow(r, q) + b*log(r);
      if (x < rho) break;
      r /= 2;
    }
    if (i == 128) {
      printf("maxiter reached in GetRFromRho0: %g %g %g\n", r, x, rho);
      Abort(1);
    }
    r0 = r;
  } else if (x < rho) {
    r0 = r;
    r *= 2;
    for (i = 0; i < 128; i++) {
      x = a*pow(r, q) + b*log(r);
      if (x > rho) break;
      r *= 2;
    }
    if (i == 128) {
      printf("maxiter reached in GetRFromRho1: %g %g %g\n", r, x, rho);
      Abort(1);
    }
    r1 = r;
  } else {
    return r;
  }
  for (i = 0; i < 256; i++) {
    r = 0.5*(r0+r1);
    x = a*pow(r, q) + b*log(r);
    if (x < rho-EPS10) {
      r0 = r;
    } else if (x > rho+EPS10) {
      r1 = r;
    } else {
      break;
    }    
    if (r1-r0 < EPS6*Min(r,EPS6)) break;
  }
  if (i == 256) {
    printf("maxiter reached in GetRFromRho: %g %g %g %g %g %g %g\n",
	   r0, r1, r, x, rho, r1-r0, x-rho);
    Abort(1);
  }
  return r;
}

double DrRho(POTENTIAL *p, int i, double x) {
  double r;

  r = x - p->rho[i];
  r = p->dr_drho[i] + r*(p->a1dr[i] + r*(p->a2dr[i] + r*p->a3dr[i]));
  return r;
}

double RadRho(POTENTIAL *p, int i, double x) {
  double r;

  r = x - p->rho[i];
  r = p->rad[i] + r*(p->a1r[i] + r*(p->a2r[i] + r*p->a3r[i]));
  return r;
}  
  
void FitDiff(double *y0, double *y1, double *y2, int i1, int i2,
	     POTENTIAL *pot) {
  int i, n, k, nd;
  double c[3], z1m, z2m, z1, z2, r, r2, r3;

  if (i2 < i1+10) return;
  n = 10;
  nd = 0;
  z1m = 0;
  z2m = 0;
  while (n <= i2-i1) {
    SVDFit(3, c, NULL, 1E-15, n, pot->rad+i1, NULL, y1+i1, NULL, PolyBasis);
    k = i1+n-1;
    r = pot->rad[k];
    r2 = r*r;
    z1 = c[0] + c[1]*r + c[2]*r2;
    z2 = c[1] + 2*c[2]*r;
    if (fabs(z1) > z1m) z1m = fabs(z1);
    if (fabs(z2) > z2m) z2m = fabs(z2);
    if (fabs(z1-y1[k]) < 1e-3*z1m &&
	fabs(z2-y2[k]) < 1e-1*z2m) {
      nd++;
    } else {
      nd = 0;
    }
    if (nd > 1) {
      break;
    }
    n = (int)(n*1.25);
  }
  if (nd > 1) {
    for (i = i1; i < i1+n; i++) {
      r = pot->rad[i];
      y1[i] = c[0] + c[1]*r + c[2]*r*r;
      y2[i] = c[1] + 2*c[2]*r;
    }
  }
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
    Differential(pot->Z, pot->dZ, 0, pot->maxrp-1, pot);
    Differential(pot->dZ, pot->dZ2, 0, pot->maxrp-1, pot);
    FitDiff(pot->Z, pot->dZ, pot->dZ2, 0, pot->maxrp-1, pot);
  } else {
    for (i = 0; i < pot->maxrp; i++) {
      pot->dZ[i] = 0;
      pot->dZ2[i] = 0;
    }
  }  
  SetPotentialVP(pot);
  SetPotentialSE(pot);

  if (pot->sps == 0) {
    SetPotentialPS(pot, -1);
  }
  
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
    if (pot->a) {
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
	y2 = -y/r + (v/r2 - (v0/r+y)*b) - (v0-v)*(pot->a*pot->a)/(x*x);
	pot->dVc2[i] += y2;
	//printf("vc1: %d %g %g %g %g %g %g\n", i, pot->Vc[i], pot->dVc[i], pot->dVc2[i], v, y, y2);
      }
    } else {
      for (i = 0; i < pot->maxrp; i++) {
	r = pot->rad[i];
	r2 = r*r;
	a = pot->lambda*r;
	if (a < 1e-3) {
	  double a2 = a*a;
	  double a3 = a2*a;
	  double nlam = n*pot->lambda;
	  double nlam2 = nlam*pot->lambda;
	  v = nlam*(1-0.5*a+a2/6.0-a3/24.0);
	  y = nlam2*(-0.5+a/3.0-a2/8.0);
	  y2 = nlam2*pot->lambda*(ONETHIRD-a*0.25);
	} else {
	  b = exp(-a);
	  v0 = n/r;
	  v = (1-b)*v0;
	  y = -v/r + v0*b*pot->lambda;
	  double nlamb = n*pot->lambda*b;
	  y2 = (-y - nlamb*pot->lambda)/r;
	  y2 += (v - nlamb)/r2;
	}
	pot->Vc[i] += v;
	pot->dVc[i] += y;
	pot->dVc2[i] += y2;
	//printf("vc1: %d %g %g %g %g %g %g %g\n", i, pot->Vc[i], pot->dVc[i], pot->dVc2[i], v, y, y2, a);
      }
    }
  }
  return 0;
}

int SetPotentialSturm(POTENTIAL *pot) {
  int i;
  double d;
  for (i = 0; i < pot->maxrp; i++) {
    d = pot->rad[i]/pot->rb;
    d = pow(d, pot->sturm_idx);
    if (d < 1e-5) {
      pot->ZPS[i] = -(1-0.5*d+d*d/6.0);
    } else {
      pot->ZPS[i] = -(1-exp(-d))/d;
    }
  }
  Differential(pot->ZPS, pot->dZPS, 0, pot->maxrp-1, pot);
  Differential(pot->dZPS, pot->dZPS2, 0, pot->maxrp-1, pot);
  return 0;
}

int AddPotentialSturm(POTENTIAL *pot, int k, double e) {
  int i;
  //printf("aps: %d %12.5E %12.5E\n", k, e, pot->sturm_ene[k]);
  pot->sturm_ene[k] += e;
  for (i = 0; i < pot->maxrp; i++) {
    pot->VT[k][i] += e*pot->ZPS[i];
    pot->dVT[k][i] += e*pot->dZPS[i];
    pot->dVT2[k][i] += e*pot->dZPS2[i];
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
      a = -pot->ZPS[i];
      b = -pot->dZPS[i];
      y = -pot->dZPS2[i];
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
    return 0;
  }

  if (u) {
    for (i = 0; i < n; i++) {
      pot->U[i] = u[i];    
    }
  }

  for (i = pot->maxrp-1; i >= 0; i--) {
    a = pot->rad[i]*pot->U[i];
    if (fabs(a) > EPS10) break;
    else pot->U[i] = 0.0;
  }
  
  Differential(pot->U, pot->dU, 0, pot->maxrp-1, pot);
  Differential(pot->dU, pot->dU2, 0, pot->maxrp-1, pot);
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

  if (k < 0) {
    for (i = 0; i < pot->maxrp; i++) {
      pot->W[i] = 0;
      pot->dW[i] = 0;
      pot->dW2[i] = 0;
    }
    return 0;
  }
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

  Differential(pot->W, pot->dW, 0, pot->maxrp-1, pot);
  Differential(pot->dW, pot->dW2, 0, pot->maxrp-1, pot);
  return 0;
}

int SetPotentialWZMC(POTENTIAL *pot, double e, int kappa, int k) {
  int i;
  double xi, r, r2, x, y, z, a2;

  if (!_zcoll) return SetPotentialW(pot, e, kappa, k);

  a2 = FINE_STRUCTURE_CONST2;
  if (_mcoll > 0) a2 /= _mcoll;
  for (i = 0; i < pot->maxrp; i++) {
    xi = e - _zcoll*pot->VT[k][i];
    r = xi*a2*0.5 + 1.0;  
    x = _zcoll*pot->dVT[k][i];
    y = - 2.0*kappa*x/pot->rad[i];
    x = x*x*0.75*a2/r;
    z = _zcoll*pot->dVT2[k][i];
    pot->W[i] = x + y + z;
    if (_mcoll > 0) pot->W[i] /= _mcoll;
    pot->W[i] /= 4.0*r;
    x = xi*xi;
    pot->W[i] = x - pot->W[i];
    pot->W[i] *= 0.5*a2;
    pot->W[i] = -pot->W[i];
  }

  Differential(pot->W, pot->dW, 0, pot->maxrp-1, pot);
  Differential(pot->dW, pot->dW2, 0, pot->maxrp-1, pot);
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
    Differential(pot->ZSE[i], pot->dZSE[i], 0, pot->maxrp-1, pot);
    Differential(pot->dZSE[i], pot->dZSE2[i], 0, pot->maxrp-1, pot);
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
  double a, b, z0, r0, r, rr;
  double v3, za, za2, za3, a2, lam, r2, r3, r5, x, y, f0, fs, fm;

  mvp = pot->mvp%100;
  vp = mvp%10;
  mvp = mvp/10;
  r0 = 3.86159E-3/RBOHR;
  rr = RBOHR/_RBOHR;
  z0 = pot->atom->atomic_number;
  a = -2.0*z0*FINE_STRUCTURE_CONST/(3.0*PI);
  b = -z0*FINE_STRUCTURE_CONST2/(PI*PI);  
  for (i = 0; i < pot->maxrp; i++) {
    r = pot->rad[i]*2.0/r0;
    pot->ZVP[i] = -a*UehlingK1(r);
  }
  if (vp > 1) {
    for (i = 0; i < pot->maxrp; i++) {
      r = pot->rad[i]*2.0/r0;
      pot->ZVP[i] -= b*UehlingL1(r);
    }
  }
  if (LEPTON_TYPE > 0) {
    r0 = 3.86159E-3;
    for (i = 0; i < pot->maxrp; i++) {
      r = pot->rad[i]*2.0/r0;
      pot->ZVP[i] -= a*UehlingK1(r);
    }
    if (vp > 1) {
      for (i = 0; i < pot->maxrp; i++) {
	r = pot->rad[i]*2.0/r0;
	pot->ZVP[i] -= b*UehlingL1(r);
      }
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
      r = (pot->rad[i]/a)*rr;
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
      pot->ZVP[i] -= (rr*za3*v3/a)*pot->rad[i];
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
  Differential(pot->ZVP, pot->dZVP, 0, m-1, pot);
  Differential(pot->dZVP, pot->dZVP2, 0, m-1, pot);
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

int SetPotentialPS(POTENTIAL *pot, int iter) {
  int i;
  if (pot->mps < 0 || pot->mps > 2 || iter < 0) {
    for (i = 0; i < pot->maxrp; i++) {
      pot->EPS[i] = 0;
      pot->VXF[i] = 0;
      pot->ZPS[i] = 0;
      pot->dZPS[i] = 0;
      pot->dZPS2[i] = 0;
    }
    return 0;
  }
  SetPotentialIPS(pot, iter);
  return 0;
}

void SetPotentialIPS(POTENTIAL *pot, int iter) {
  int i, ii, maxit;
  double n0, r0, x;
  if (pot->mps == 1) {//debye screening
    if (iter > 0) return;
    for (i = 0; i < pot->maxrp; i++) {
      pot->EPS[i] = 0;
      r0 = pot->rad[i]/pot->dps;
      x = exp(-r0);
      double xd = x/pot->dps;
      double xdd = xd/pot->dps;
      if (_debye_mode == 0) {
	pot->ZPS[i] = pot->zps*(x-1.0);
	pot->dZPS[i] = -pot->zps*xd;
	pot->dZPS2[i] = pot->zps*xdd;
      } else if (_debye_mode == 1) {
	pot->ZPS[i] = pot->Z[i]*(x-1.0);
	double xd = x/pot->dps;
	pot->dZPS[i] = -pot->Z[i]*xd + pot->dZ[i]*(x-1.0);
	pot->dZPS2[i] = pot->Z[i]*xdd - 2*pot->dZ[i]*xd + pot->dZ2[i]*(x-1.0);
      } else if (_debye_mode == 2) {
	double y = x*cos(r0);
	pot->ZPS[i] = pot->zps*(y-1.0);	
      } else if (_debye_mode == 3) {
	double y = x*cos(r0);
	pot->ZPS[i] = pot->Z[i]*(y-1.0);
      }
    }
    if (_debye_mode >= 2) {
      Differential(pot->ZPS, pot->dZPS, 0, pot->maxrp-1, pot);
      Differential(pot->dZPS, pot->dZPS2, 0, pot->maxrp-1, pot);      
    }
    return;
  }
  
  if (pot->ips == 0) {
    for (i = 0; i < pot->maxrp; i++) {
      pot->EPS[i] = 0;
      pot->ZPS[i] = 0;
      pot->dZPS[i] = 0;
      pot->dZPS2[i] = 0;
    }
    return;
  }

  if (pot->tps <= 0) {
    r0 = pot->rps;
    n0 = pot->nps;
    for (i = 0; i < pot->maxrp; i++) {      
      if (pot->rad[i] <= r0) {
	x = pot->rad[i]/r0;
	pot->EPS[i] = FOUR_PI*pot->rad[i]*pot->rad[i]*n0;
	pot->ZPS[i] = 0.5*x*pot->zps*(3-x*x);
	pot->dZPS[i] = -1.5*pot->zps*(x*x-1.0)/r0;
	pot->dZPS2[i] = -3.0*pot->zps*x/(r0*r0);
      } else {
	pot->EPS[i] = 0.0;
	pot->ZPS[i] = pot->zps;
	pot->dZPS[i] = 0.0;
	pot->dZPS2[i] = 0.0;
      }
    }
    return;
  }

  if (pot->ups > 0) {
    if (iter <= 0) return;
    if (iter > 1 && (_sp_mode == 0 || _sp_mode == 2)) return;
  }
  
  double x0, x1, dx, tol;
  maxit = OptimizeMaxIter()+1;
  tol = EPS7;
  i = 0;
  ii = 0;
  dx = 0;
  x0 = 0.0;
  x1 = 0.0;
  x = pot->aps;
  r0 = Max(1,pot->zps)*tol;
  pot->xps = 1.0;
  n0 = FreeElectronDensity(pot, 0.0, x, pot->zps, -1, iter);
  if (pot->ups > 0 && _ionsph_bmode == 0) {
    if (_sp_mode < 0) {
      goto END;
    }
    if (fabs(n0-pot->zps) < r0) goto END;
    dx = 2.0;
    if (n0 < pot->zps) {
      x0 = pot->xps;
      x1 = x0;
      while (n0 < pot->zps) {
	x0 = x1;
	x1 *= dx;
	pot->xps = x1;
	n0 = FreeElectronDensity(pot, 0.0, x, pot->zps, -1, ii*maxit+iter);
	if (ii > _ionsph_maxiter-10 || (_sp_print > 0 && ii%_sp_print == 0)) {
	  printf("ion sphere iter0: %d %d %g %g %g %g %g\n",
		 iter, ii, pot->zps, n0, pot->xps, x0, x1);
	}
	ii++;
	if (fabs(n0-pot->zps) < r0) goto END;
	dx *= 2;
	if (ii >= _ionsph_maxiter) {
	  printf("max ion sphere iter0: %d %d %g %g %g %g %g\n",
		 iter, ii, pot->zps, n0, pot->xps, x0, x1);
	  break;
	}
      }
    } else {
      x1 = pot->xps;
      x0 = x1;
      while (n0 > pot->zps) {
	x1 = x0;
	x0 /= dx;
	pot->xps = x0;
	n0 = FreeElectronDensity(pot, 0.0, x, pot->zps, -1, ii*maxit+iter);
	if (ii > _ionsph_maxiter-10 || (_sp_print > 0 && ii%_sp_print == 0)) {
	  printf("ion sphere iter1: %d %d %g %g %g %g %g\n",
		 iter, ii, pot->zps, n0, pot->xps, x0, x1);
	}
	ii++;
	if (fabs(n0-pot->zps) < r0) goto END;
	dx *= 2;
	if (ii >= _ionsph_maxiter) {
	  printf("max ion sphere iter1: %d %d %g %g %g %g %g\n",
		 iter, ii, pot->zps, n0, pot->xps, x0, x1);
	  break;
	}
      }
    }
    i = ii;
    pot->xps = sqrt(x0*x1);
    while(fabs(n0-pot->zps) > r0) {
      i++;
      n0 = FreeElectronDensity(pot, 0.0, x, pot->zps, -1, i*maxit+iter);
      if (i > _ionsph_maxiter-10 || (_sp_print > 0 && i%_sp_print == 0)) {
	printf("ion sphere iter2: %d %d %d %g %g %g %g %g %g %g %g %g %g\n",
	       iter, i, ii, pot->zps, n0, dx, x0, x1, x,
	       pot->xps, pot->aps, pot->rps, pot->dps);
      }
      if (n0 < pot->zps) {
	x0 = pot->xps;
      } else if (n0 > pot->zps) {
	x1 = pot->xps;
      } else {
	break;
      }
      if (x1/x0 > 2) {
	pot->xps = sqrt(x0*x1);
      } else {
	pot->xps = 0.5*(x0+x1);
      }
      if (fabs(x1-x0)/pot->xps < EPS3*tol && n0 <= pot->zps) break;
      if (i >= _ionsph_maxiter) {
	printf("max ion sphere iter2: %d %d %g %g %g %g %g\n",
	       iter, i, pot->zps, n0, pot->xps, x0, x1);
	break;
      }
    }
    goto END;
  }
  dx = 2*fabs(pot->aps-x);
  if (dx < 0.1) dx = 0.1;
  if (fabs(n0-pot->zps) < r0) goto END;
  if (n0 < pot->zps) {
    x0 = x;
    x1 = x;
    while (n0 < pot->zps) {
      x0 = x1;
      x1 = pot->aps + dx;
      n0 = FreeElectronDensity(pot, 0.0, x1, pot->zps, -1, ii*maxit+iter);      
      ii++;
      if (fabs(n0-pot->zps) < r0) goto END;
      if (n0 < pot->zps) dx *= 2;
      if (ii >= _ionsph_maxiter) {
	printf("max ion sphere iter3: %d %d %g %g %g %g %g\n",
	       iter, ii, pot->zps, n0, dx, x0, x1);
	break;
      }
    }
  } else {
    x1 = x;
    x0 = x;
    while (n0 > pot->zps) {
      x1 = x0;
      x0 = pot->aps - dx;
      n0 = FreeElectronDensity(pot, 0.0, x0, pot->zps, -1, ii*maxit+iter);
      ii++;
      if (fabs(n0-pot->zps) < r0) goto END;
      if (n0 > pot->zps) dx *= 2;
      if (ii >= _ionsph_maxiter) {
	printf("max ion sphere iter4: %d %d %g %g %g %g %g\n",
	       iter, ii, pot->zps, n0, dx, x0, x1);
	break;
      }
    }
  }
  i = ii;
  x = 0.5*(x0+x1);
  while (fabs(n0-pot->zps) > r0) {
    i++;    
    n0 = FreeElectronDensity(pot, 0.0, x, pot->zps, -1, i*maxit+iter);
    if (i > _ionsph_maxiter-10 || (_sp_print > 0 && i%_sp_print == 0)) {
      printf("ion sphere iter5: %d %d %d %g %g %g %g %g %g %g %g %g %g\n",
	     iter, i, ii, pot->zps, n0, dx, x0, x1, x,
	     pot->xps, pot->aps, pot->rps, pot->dps);
    }
    if (n0 < pot->zps) {
      x0 = x;
    } else if (n0 > pot->zps) {
      x1 = x;
    } else {
      break;
    }
    x = 0.5*(x0+x1);
    if (i >= _ionsph_maxiter) {
      printf("max ion sphere iter5: %d %d %g %g %g\n",
	     iter, i, pot->zps, n0, pot->aps);
      break;
    }
  }
 END:
  if (_sp_print < 0) {
      printf("ion sphere iter6: %d %d %d %g %g %g %g %g %g %g %g %g %g\n",
	     iter, i, ii, pot->zps, n0, dx, x0, x1, x,
	     pot->xps, pot->aps, pot->rps, pot->dps);
  }
}

double XCPotential(double tx, double n, int md) {
  double p, r, rs, ti, tis, t, t2, t3, t4, at, bt, ct, dt, et;
  double nx[101], nx0, nx1, dn;
  int i;
  
  const double a[6] = {0.610887, 3.04363, -0.09227, 1.7035, 8.31051, 5.1105};
  const double b[5] = {0.283997, 48.932154, 0.370919, 61.095357, 0.871837};
  const double c[3] = {0.870089, 0.193077, 2.414644};
  const double d[5] = {0.579824, 94.537454, 97.839603, 59.939999, 24.388037};
  const double e[5] = {0.212036, 16.731249, 28.485792, 34.028876, 17.235515};
  const double f[5] = {2.8343, -0.2151, 5.2759, 3.9431, 7.9138};  

  if (n < 1E-35) return 0.0;
  switch (md) {
  case 0:
    return pow(FOUR_PI*n, ONETHIRD);
  case 1:
    t = 2*tx/pow(3*PI*PI*n, TWOTHIRD);
    p = pow(FOUR_PI*n, ONETHIRD);
    if (t < 1e-10) {
      return p;
    } else if (t > 1e10) {
      return (p/t)*(f[2]/f[4]);
    } else {
      ti = 1/ti;
      t2 = t*t;
      t3 = t2*t;
      t4 = t3*t;
      p *= tanh(ti);
      p *= 1 + f[0]*t2 + f[1]*t3 + f[2]*t4;
      p /= 1 + f[3]*t2 + f[4]*t4;
      return p;
    }
  case 11:
  case 21:
  case 31:
    nx0 = log(0.01*n);
    nx1 = log(n);
    dn = (nx1 - nx0)/100.;
    r = exp(nx0);
    p = 0.75*r*XCPotential(tx, r, 1);
    for (i = 0; i < 101; i++) {
      nx[i] = exp(nx0+dn*i);
      nx[i] = XCPotential(tx, nx[i], 1)*nx[i];
    }
    p += Simpson(nx, 0, 100)*dn;
    p *= 1.33333333333/n;
    if (md == 11) {
      r = pow(3/(FOUR_PI*n), ONETHIRD);
      return (r/3.14788045)*p;
    } else if (md == 21) {
      return p/3.14788045;
    } else {    
      return p;
    }
  case 2:    
    r = pow(3/(FOUR_PI*n), ONETHIRD);
    t = n*EPS3;
    t2 = n - t;
    t3 = n + t;
    p = (XCPotential(tx, t3, 12) - XCPotential(tx, t2, 12))/(2*t);
    p = p*n/r + 1.3333333333333*XCPotential(tx, n, 22);
    p /= 0.42356543;
    return p;
  case 12:
  case 22:
  case 32:
    if (tx > 0) {
      t = 2*tx/pow(3*PI*PI*n, TWOTHIRD);
    } else {
      t = 0.0;
    }
    r = pow(3/(FOUR_PI*n), ONETHIRD);
    if (t < 1e-10) {
      at = a[0]*0.75;
      bt = b[0];
      dt = d[0];
      et = e[0];
      ct = c[0]*et;
    } else if (t > 1e10) {
      ti = 1/t;
      tis = sqrt(ti);
      at = a[0]*ti*a[3]/a[5];
      bt = tis*b[2]/b[4];
      dt = tis*d[2]/d[4];
      et = ti*e[2]/e[4];
      ct = (c[0]+c[1])*et;
    } else {
      ti = 1/t;
      tis = tanh(sqrt(ti));
      ti = tanh(ti);
      t2 = t*t;
      t4 = t2*t2;
      t3 = t2*t;
      at = a[0]*ti*(0.75 + a[1]*t2 + a[2]*t3 + a[3]*t4);
      at /= (1.0 + a[4]*t2 + a[5]*t4);
      bt = tis*(b[0] + b[1]*t2 + b[2]*t4);
      bt /= (1 + b[3]*t2 + b[4]*t4);
      dt = tis*(d[0] + d[1]*t2 + d[2]*t4);
      dt /= (1 + d[3]*t2 + d[4]*t4);
      et = ti*(e[0] + e[1]*t2 + e[2]*t4);
      et /= (1 + e[3]*t2 + e[4]*t4);
      ct = (c[0] + c[1]*exp(-c[2]/t))*et;
    }
    rs = sqrt(r);  
    p = (at + bt*rs + ct*r)/(1 + dt*rs + et*r);
    if (md == 12) {
      return p;
    } else if (md == 22) {
      return p/r;
    } else {
      return (3.14788045/r)*p;
    }
  default:
    return 0.0;
  }
}
    
int DensityToSZ(POTENTIAL *pot, double *d, double *z,
		double *zx, double *jps, int md) {
  int i, im, maxrp;
  double nx, rx;

  if (pot->ips > 0 && pot->ups+1 == 1) {
    maxrp = pot->ips+1;
  } else {
    maxrp = pot->maxrp;
  }
  for (im = maxrp-1; im >= 0; im--) {
    if (d[im] || zx[im]) break;
  }
  if (im <= 0) {
    for (i = 0; i < pot->maxrp; i++) {
      z[i] = 0.0;
      zx[i] = 0.0;
    }
    *jps = 0.0;
    return 0;
  }
  im = maxrp-1;
  for (i = 0; i <= im; i++) {
    _dwork[i] = d[i]*pot->dr_drho[i];
  }
  _dwork1[0] = d[0]*pot->rad[0]/3.0;
  NewtonCotesIP(_dwork1, _dwork, 0, im, -1, 0);
  for (i = 0; i <= im; i++) {
    _dwork[i] = (d[i]/pot->rad[i])*pot->dr_drho[i];
  }
  _dwork2[im] = 0.0;
  NewtonCotesIP(_dwork2, _dwork, 0, im, -1, -1);
  for (i = 0; i <= im; i++) {
    z[i] = _dwork1[i] + pot->rad[i]*_dwork2[i];
    if (pot->ahx) {
      rx = pot->rad[i];       
      nx = zx[i]/(FOUR_PI*rx*rx);
      if (nx > 0) {
	zx[i] = pot->ahx*XCPotential(pot->tps, nx, pot->vxm)*rx;
      } else {
	zx[i] = 0.0;
      }
    } else {
      zx[i] = 0.0;
    }
  }
  for (i = im+1; i < pot->maxrp; i++) {
    z[i] = z[im];
    zx[i] = zx[im];
  }

  *jps = _dwork2[0];
  return im;
}

double BoundFactor(double e, double eth, double de, int mps) {
  double x, y;

  if (_sc_ewm < 0 || mps < 3) {
    return 0.0;
  }
  
  if (de <= 0) {
    if (eth-e > 1e-7) return 1.0;
    return 0.0;
  }

  if (_sc_ewm > 0) {
    if (e < eth) return 1.0;
    x = (e - eth)/de;
    y = 1 - ERF(x);
    y = Min(1.0, y);
    y = Max(0.0, y);
  } else {
    x = (e - eth)/de;
    y = 0.5*(1 - ERF(x));
    y = Min(1.0, y);
    y = Max(0.0, y);
  }
  return y;  
}

double FreeElectronDensity(POTENTIAL *pot, double e0, double u,
			   double zn, int md, int iter) {
  int i;
  if (pot->ups > 0 && pot->dps > 0) {
    for (i = 0; i < pot->maxrp; i++) {
      pot->EPS[i] = 0.0;
      pot->ICF[i] = 0.0;
    }
    double *aps = &pot->aps;
    if (pot->mps == 0 && pot->ups > 0 && _ionsph_bmode == 0) {
      aps = NULL;
    }
    return FreeElectronIntegral(pot, 0, pot->ips, pot->maxrp-1,
				pot->EPS, pot->ICF,
				e0, pot->tps, pot->eth, pot->ewd,
				u, zn, md, iter, pot->ups, pot->fps,
				pot->rps, pot->dps, aps);
  } else {
    for (i = 0; i < pot->maxrp; i++) {
      pot->EPS[i] = 0.0;
    }
    return FreeElectronIntegral(pot, 0, pot->ips, pot->ips,
				pot->EPS, NULL,
				e0, pot->tps, pot->eth, pot->ewd,
				u, zn, md, iter, 0.0, 0.0,
				0.0, 0.0, &pot->aps);
  }
}

double StewartPyattIntegrand(double x, double a, double fa, double y,
			     double y0, double g, double z, double nb) {
  double ys, zs, r, r0;
  int m;

  if (y0 > 0) {
    m = 1;
  } else {
    m = 0;
  }
  r0 = InterpFermiRM1(y, m);
  if (nb > 0) {
    r0 += nb;
  }
  r = 0.0;
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
    if (z >= 0) {
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
  r = (r0+r)/(1+z);
  return r;
}

void DerivSPY(int *neq, double *t, double *y, double *yd) {
  double x = *t;
  yd[0] = y[1];
  double rho = 0.0;  
  if (_icf_nfft > 0) {
    rho = x*y[9];
    rho = y[6]*pow(rho, y[8]) + y[7]*log(rho);
  } else {
    rho = x/y[10];
  }
  double z;
  z = (x - y[14])/(y[15] - y[14]);
  z = (1-z)*y[12] + z*y[13];  
  double y0 = (y[0]-z)/x;
  double y1 = 1.0;
  double s = StewartPyattIntegrand(rho, y[2], y[3], y0, y1,
				   y[4], y[5], y[11]);
  yd[1] = x*s;
}

double FreeElectronIntegral(POTENTIAL *pot, int i0, int i1, int i2,
			    double *eps, double *icf,
			    double eth, double tps,
			    double eref, double ewd,
			    double u, double zn, int md, int miter,
			    double ups, double fa,
			    double rps, double dps, double *aps) {
  int i, k, ii, mps, iter, ips, maxit;
  double a = 4.0/PI, g, e0, x, y, ye, r2, y0, y1, xk, xj, xs, x1;
  double x0, xm, xd, yr, y2, ym, u1, u2, u3, u6, x2, x3, b;
  double *rad, *drdx, *vt;

  maxit = OptimizeMaxIter() + 1;
  iter = miter%maxit;
  ips = miter/maxit;
  mps = pot->mps;
  rad = pot->rad;
  drdx = pot->dr_drho;
  vt = pot->VT[0];
  if (md < 0) {
    md = _relativistic_fermi;
  }
  if (md >= 2) {
    double c2 = 1.0/FINE_STRUCTURE_CONST2;
    double a2 = FINE_STRUCTURE_CONST2;
    int nk = MAXNKF;
    double v, vr, k0, k1, k2, dk, r2, ek, ek0, ekr;
    double e0i, e1i, uv, emax, e0v, g2, g3, g32, ig32;
    double eth10 = 0.0;
    if (eth > 0) eth10 = 10*eth;
    double ew5, ews, em, ep, eb, ts, ef, ef0, ef1, bf, es[8];
    
    if (md < 10) {
      g = md-2.0;
    } else {
      g = 0.0;
    }
    ts = tps*0.1;
    ew5 = ewd*5;
    ews = ewd*0.25;
    g2 = 0.5*g;
    g3 = g+3;
    g32 = 0.5*(g+3);
    ig32 = 1/g32;
    a /= g3;
    for (i = i0; i <= i1; i++) {
      r2 = rad[i]*rad[i];
      v = vt[i];
      vr = v/tps;
      uv = u-vr;
      emax = tps*(Max(uv,vr)+100.0);
      emax = Max(eth10, emax);
      em = (uv-3.0)*tps;
      ep = (uv+3.0)*tps;
      eb = (uv+10.0)*tps;
      ef = eref - v;
      if (_sc_ewm <= 0) {
	es[0] = ef - ew5;
      } else {
	es[0] = ef;
      }
      es[1] = ef + ew5;
      ef0 = es[0];
      ef1 = es[1];
      es[2] = em;
      es[3] = ep;
      es[4] = eb;
      es[5] = (uv+20.0)*tps;
      es[6] = (uv+50.0)*tps;
      es[7] = emax;
      qsort(es, 8, sizeof(double), CompareDouble);
      e0 = eth;
      if (e0 < v) e0 = v;
      e0v = e0-v;
      eps[i] = 0.0;
      for (ii = 0; ii < 8; ii++) {
	if (e0v < es[ii]) break;
      }
      e0i = e0v;
      for (; ii < 8; ii++) {
	e1i = es[ii];
	if (e1i-e0i < 1e-10*tps) continue;
	if (e0i < ef1 && e1i > ef0) {
	  nk = 1+(e1i-e0i)/Min(ews,ts);	  
	} else {
	  nk = 1+(e1i-e0i)/ts;
	}
	if (nk < MINNKF) nk = MINNKF;
	if (nk > MAXNKF) nk = MAXNKF;
	if (_relativistic_fermi) {
	  k0 = 2*e0i + a2*e0i*e0i;
	  k1 = 2*e1i + a2*e1i*e1i;
	} else {
	  k0 = 2*e0i;
	  k1 = 2*e1i;
	}
	k0 = pow(k0, g32);
	k1 = pow(k1, g32);
	dk = (k1 - k0)/(nk-1);
	for (k = 0; k < nk; k++) {
	  k2 = k0 + k*dk;
	  k2 = pow(k2, ig32);
	  if (k2 < 1e-6*c2 || !_relativistic_fermi) {
	    ek0 = 0.5*k2;
	    ekr = 0.5;
	  } else {
	    ek0 = sqrt(k2*c2 + c2*c2)-c2;
	    ekr = ek0/k2;
	  }
	  ek = ek0 + v;
	  y = ek/tps - u;
	  if (y > 0) {
	    y = exp(-y);
	    y /= 1+y;
	  } else {
	    y = 1/(1+exp(y));
	  }
	  bf = BoundFactor(ek, eref, ewd, mps);
	  _kwork[k] = 1.0-bf;
	  if (md < 10) {
	    _kwork[k] *= y;
	    if (md > 2) {
	      _kwork[k] *= pow(ekr, g2);
	    }
	  } else {
	    ye = 0.0;
	    if (y > 0) {
	      ye += y*log(y);
	    }
	    if (y < 1) {
	      ye += (1-y)*log(1-y);
	    }
	    _kwork[k] *= -ye;
	  }
	}
	y = dk*Simpson(_kwork, 0, nk-1);
	eps[i] += y;
	if (y < EPS8*eps[i]) break;
	e0i = e1i;
      }
      eps[i] *= a*r2;
      if (eps[i] < 1E-99) eps[i] = 0.0;
    }
    if (ups > 0 && icf) {
      if (_relativistic_fermi) {
	g = tps*FINE_STRUCTURE_CONST2;
      } else {
	g = 0.0;
      }
      a = (4/PI)*sqrt(2*tps)*tps;
    }
  } else {    
    a *= sqrt(2*tps)*tps;
    if (md) g = tps*FINE_STRUCTURE_CONST2;
    else g = 0;
    int ifermi = 0;
    y0 = 0.0;
    if (mps == 0 && ups > 0 && _ionsph_bmode == 0) {
      ifermi = 1;
      y0 = fa;
    } else {
      if (_ionsph_ifermi > 0 && i2-i0+1 > _ionsph_ifermi*NFRM1) {
	y0 = FermiIntegral(u, 0.0, g);
	if (y0 > 1E-30) {
	  PrepFermiRM1(u, y0, tps);
	} else {
	  y0 = 0.0;
	}
	ifermi = 1;
      }
    }
    if (ups > 0 && icf && dps > 0) {      
      u1 = ups+1.0;
      if (ifermi) {
	b = ups+_fermi_rm1[2][NFRM1];
	b = sqrt(b/u1);	
      } else {
	b = 1.0;
      }
      u2 = 2*u1;
      u3 = 3*u1;
      u6 = 6*u1;
      x = b*rps/dps;
      x2 = x*x;
      x3 = x2*x;
      if (x < 1e-5) {
	x1 = x3/(3*b);
	xj = (2*x1/b + x1*x1)/u2;
      } else if (x > 1e10) {
	x1 = (x - 1.0)/b;
      } else {
	x1 = (pow(1.0 + x3, ONETHIRD)-1)/b;
      }
      xj = (2*x1/b +x1*x1)/u2;
      xk = x3/u3/(b*b*b);
      ye = (xj - x1*x1/u2)/b;
      if (_ionsph_ifermi &&
	  (_sp_mode%10 == 2 || _sp_mode%10 == 3)
	  && iter > 0) {
	double *rwork = _dwork3;
	int iwork[22];
	double rtol, atol[2], ysp[16], t, t0, kc, dt, dx, hdx;
	int neq, itol, itask, istate, iopt, mf, lrw, liw, ib;
	lrw = pot->maxrp;
	liw = 22;
	neq = 2;
	itol = 2;
	rtol = EPS3;
	atol[0] = EPS6;
	atol[1] = EPS6;
	itask = 4;
	istate = 1;
	iopt = 0;
	mf = 22;
	ysp[2] = u;
	ysp[3] = y0;
	ysp[4] = g;
	ysp[5] = pot->ups;
	ysp[6] = pot->ar;
	ysp[7] = pot->br;
	ysp[8] = pot->qr;
	ysp[9] = pot->dps;
	ysp[10] = pot->rps/pot->dps;
	dt = pot->dps*pot->tps;
	if (_sp_mode >= 10) {
	  t0 = -pot->zps/dt;
	  for (i = pot->maxrp-1; i >= 0; i--) {
	    x0 = pot->rad[i]/pot->dps;	    
	    t = ye*exp(b*(x1-x0));
	    _veff[i] = (t0+t*pot->xps)/x0;
	    if (fabs(t/t0) >= EPS2) break;
	  }
	  ysp[0] = t0 + t*pot->xps;
	  ysp[1] = -b*t*pot->xps;
	} else {
	  i = 0;
	  x0 = pot->rad[i]/pot->dps;
	  ysp[1] = -xj/pot->xps;
	  ysp[0] = x0*ysp[1];
	}
	t0 = x0;	
	ysp[12] = (vt[i]*pot->rad[i]-pot->ZPS[i])/dt;
	ysp[14] = t0;
	_veff[i] = ysp[0]/x0;
	_dwork1[i] = ysp[1];
	if (_sp_mode >= 10) {
	  for (i--; i >= 0; i--) {
	    x = pot->rad[i]/pot->dps;
	    dx = x - x0;
	    hdx = 0.5*dx;
	    ysp[11] = 0.0;
	    ysp[13] = (vt[i]*pot->rad[i]-pot->ZPS[i])/dt;
	    ysp[15] = x;
	    t = x0;
	    rwork[0] = ysp[0];
	    rwork[1] = ysp[1];
	    DerivSPY(&neq, &t, ysp, rwork+10);
	    t = x0 + hdx;
	    ysp[0] = rwork[0] + rwork[10]*hdx;
	    ysp[1] = rwork[1] + rwork[11]*hdx;
	    DerivSPY(&neq, &t, ysp, rwork+12);
	    ysp[0] = rwork[0] + rwork[12]*hdx;
	    ysp[1] = rwork[1] + rwork[13]*hdx;
	    DerivSPY(&neq, &t, ysp, rwork+14);
	    ysp[0] = rwork[0] + rwork[14]*dx;
	    ysp[1] = rwork[1] + rwork[15]*dx;
	    t = x;
	    DerivSPY(&neq, &t, ysp, rwork+16);
	    ysp[0] = rwork[0] + (hdx/3.)*(rwork[10]+
					  2*rwork[12]+
					  2*rwork[14]+
					  rwork[16]);
	    ysp[1] = rwork[1] + (hdx/3.)*(rwork[11]+
					  2*rwork[13]+
					  2*rwork[15]+
					  rwork[17]);
	    x0 = x;
	    _veff[i] = ysp[0]/x;
	    ysp[12] = ysp[13];
	    ysp[14] = ysp[15];
	  }
	} else {
	  t0 = -pot->zps/dt;
	  rtol = ysp[0] - t0;
	  for (i++; i < pot->maxrp; i++) {
	    x = pot->rad[i]/pot->dps;
	    dx = x - x0;
	    hdx = 0.5*dx;
	    ysp[11] = 0.0;
	    ysp[13] = (vt[i]*pot->rad[i]-pot->ZPS[i])/dt;
	    ysp[15] = x;
	    t = x0;
	    rwork[0] = ysp[0];
	    rwork[1] = ysp[1];
	    DerivSPY(&neq, &t, ysp, rwork+10);
	    _dwork2[i-1] = rwork[11];
	    t = x0 + hdx;
	    ysp[0] = rwork[0] + rwork[10]*hdx;
	    ysp[1] = rwork[1] + rwork[11]*hdx;
	    DerivSPY(&neq, &t, ysp, rwork+12);
	    ysp[0] = rwork[0] + rwork[12]*hdx;
	    ysp[1] = rwork[1] + rwork[13]*hdx;
	    DerivSPY(&neq, &t, ysp, rwork+14);
	    ysp[0] = rwork[0] + rwork[14]*dx;
	    ysp[1] = rwork[1] + rwork[15]*dx;
	    t = x;
	    DerivSPY(&neq, &t, ysp, rwork+16);
	    ysp[0] = rwork[0] + (hdx/3.)*(rwork[10]+
					  2*rwork[12]+
					  2*rwork[14]+
					  rwork[16]);
	    ysp[1] = rwork[1] + (hdx/3.)*(rwork[11]+
					  2*rwork[13]+
					  2*rwork[15]+
					  rwork[17]);
	    if (fabs((ysp[13]-ysp[12])/ysp[13]) < EPS3) {
	      if (x > x1) {
		t = ye*exp(b*(x1-x));
		if (fabs(t/t0) < EPS2) break;
		t = ysp[0]-t0;
		if (t < rtol) {
		  rtol = t;
		} else {
		  break;
		}
	      }
	      t = (ysp[0]-ysp[13])/x;
	      if (t < EPS4/pot->ups) break;
	    }
	    ysp[12] = ysp[13];
	    ysp[14] = ysp[15];
	    _veff[i] = ysp[0]/x;
	    _dwork1[i] = ysp[1];
	    x0 = x;
	  }
	  k = i-1;
	  if (i < pot->maxrp) {
	    ysp[0] = _veff[k]*x0;
	    t = ysp[0]-t0;
	    for (; i < pot->maxrp; i++) {
	      x = pot->rad[i]/pot->dps;
	      _veff[i] = (t0 + t*exp(b*(x0-x)))/x;
	    }
	  }
	  if (_sp_piter0 > 0 || _sp_piter1 > 0) {
	    if ((_sp_piter0 == 0 || iter == _sp_piter0) &&
		(_sp_piter1 == 0 || ips == _sp_piter1)) {
	      for (i = 0; i < pot->maxrp; i++) {
		x = pot->rad[i]/pot->dps;
		t = (vt[i]*pot->rad[i]-pot->ZPS[i])/dt;
		dx = ExpM1(-pot->ups*(_veff[i]-t/x));
		hdx = pot->NPS[i]/(FOUR_PI*pot->rad[i]*pot->rad[i]*pot->nps);
		MPrintf(0, "spy: %4d %4d %4d %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E\n", iter, ips, i, pot->rad[i], b, ye, x1, x, t, _veff[i]*x, _veff[i]*x-t, _dwork1[i], _dwork2[i], hdx);
	      }
	    }
	  }
	}
      }    
     
      for (i = i0; i <= i2; i++) {
	if (ifermi && 1+y0==1) {
	  eps[i] = 0.0;
	  icf[i] = 0.0;
	  continue;
	}
	r2 = rad[i]*rad[i];
	x = rad[i]/dps;
	ym = xk/x;
	if (x <= x1) {
	  y1 = ym - xj + x*x/u6;
	} else {
	  y1 = (ye/x)*exp(b*(x1-x));
	}
	double vc = (vt[i]-pot->ZPS[i]/rad[i])/tps;
	double xf = pot->xps;
	if (_sp_mode == 1) {
	  y1 = xf * ((y1-ym) - vc);
	} else if (_sp_mode == 2 || _sp_mode == 12) {
	  y1 = _veff[i] - vc;
	} else if (_sp_mode == 3 || _sp_mode == 13) {
	  y1 = _veff[i] - vc;
	} else if (_sp_mode == 4) {
	  if (i <= pot->ips) {
	    y1 = xf * (-vt[i]/tps);
	  } else {
	    k = pot->ips;
	    x2 = rad[k]/dps;
	    xs = xk/x2;
	    if (x2 <= x1) {
	      y2 = xs - xj + x2*x2/u6;
	    } else {
	      y2 = (ye/x2)*exp(b*(x1-x2));
	    }
	    y1 = (y1-ym) - vc;
	    y2 = (y2-xs) - vc;
	    y1 = xf * y1 * (-vt[k])/(tps*y2);
	  }
	} else if (_sp_mode == 5) {
	  y1 = xf * (-vt[i]/tps);
	}

	if (ifermi) {
	  xs = InterpFermiRM1(y1, 1);	  
	  x2 = ExpM1(-ups*y1);
	  if (xs > x2) {
	    eps[i] = a*r2*y0*(xs - x2);
	  } else {
	    eps[i] = 0.0;
	  }
	  icf[i] = 1 - (x2+1)/(xs+1.0);
	  if (icf[i] < 0) icf[i] = 0.0;
	  _dwork1[i] = -ExpM1(-ups*y1);
	} else {
	  eps[i] = a*r2*FermiIntegral(y1+u, y1, g);
	  FERMID(0.5, u, EPS10, &y0, &ii);
	  FERMID(0.5, y1+u, EPS10, &y, &ii);
	  icf[i] = (1 - (y0/y)*exp(-ups*y1));
	  if (icf[i] < 1E-99) icf[i] = 0.0;
	  eps[i] *= icf[i];
	  _dwork1[i] = -ExpM1(-ups*y1);
	}
	if (eps[i] < 1E-99) eps[i] = 0.0;
      }
    } else {      
      for (i = i0; i <= i2; i++) {
	if (ifermi && 1+y0 == 1) {
	  eps[i] = 0.0;
	  continue;
	}
	r2 = rad[i]*rad[i];
	y = -vt[i]/tps;
	if (ifermi) {
	  ye = InterpFermiRM1(y, 1) + 1;
	  eps[i] = a*r2*ye*y0;
	} else {
	  eps[i] = a*r2*FermiIntegral(y+u, y, g);
	}
	if (eps[i] < 1E-99) eps[i] = 0.0;
	_dwork1[i] = -ExpM1(-ups*y);
      }
    }
    if (drdx && i2 > i0) {
      for (i = i0; i <= i2; i++) {
	_dwork1[i] *= rad[i]*rad[i]*drdx[i];
      }
      pot->cps = 3*Simpson(_dwork1, i0, i2);
      if (i0 == 0) {
	pot->cps += rad[0]*rad[0]*rad[0];
      }
      pot->cps = pow(pot->cps, ONETHIRD);
    }
  }
  if (drdx && i2 > i0) {
    if (i0 == 0) {
      _dwork1[i0] = eps[i0]*rad[i0]/3.0;
    } else {
      _dwork1[i0] = 0.0;
    }
    for (i = 0; i <= i2; i++) {
      _dwork2[i] = eps[i]*drdx[i];
    }    
    NewtonCotesIP(_dwork1, _dwork2, i0, i2, -1, 0);
    y0 = _dwork1[i2];
    if (zn > 0) {
      a = zn/y0;
      if (mps == 0 && aps) {
	*aps = u + log(a);
      }
      for (i = i0; i <= i2; i++) {
	eps[i] *= a;
      }
    }
    return y0;
  } else {
    return eps[i0];
  }
}

static double _fermi_ax = 0.0;
static double _fermi_ag = 0.0;
double FermiIntegrand(double *tp) {
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
    x -= gt;
    if (x < -100.) return 0.0;
    et = exp(x);
    r = sqrt(t)*(et/(1.0+et));
  } else {
    x = gt-x;
    if (x < -100.) et = 0.0;
    else et = exp(x);
    r = sqrt(t)/(1.0+et);
  }
  return r;
}

double FermiIntegral(double x, double y, double g) {
  double t0, t1;
  double r, s, a, ra;
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
    double dt = Max(_fermi_tmin*t0, 1.0);
    dt = Min(_fermi_tmax, dt);
    t1 = t0 + dt;
    DQAGS(FermiIntegrand,
	  t0, t1, _fermi_abserr, _fermi_relerr, &s, &a, &n, &ier,
	  lim, lenw, &last, iwork, dwork);
    r += s;
    t0 = t1;
    ra = Max(r, _fermi_abserr);
    if (iter > 10 && s < 1e-3*ra) {
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
	if (s < _fermi_relerr*ra) {
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
      if (i > 1024) {
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
      if (i > 1024) {
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
    if (i > 1024) {
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
      r = ExpM1(y);//+a0*exp(2*y+a)*ExpM1(-y);
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
  double y0, r0;
  im = -1;
  for (i = 0; i < n; i++) {
    y0 = exp(y[i]);
    if (m == 0) {
      r0 = FermiIntegral(y0+a, 0.0, g)/fa;
      r[i] = r0 - 1.0;
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
      r0 = FermiIntegral(y0+a, y0, g)/fa;
      r[i] = r0-1.0;
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
    } else if (m == 1) {
      i = 2;
    }
  }
  a = _fermi_rm1[0][NFRM1];
  g = _fermi_rm1[0][NFRM1+1];
  double ya = fabs(y0);
  if (1+ya == 1) {
    return 0.0;
  }
  y = log(ya);
  if (y < _fermi_rm1[0][0]) {
    r = FermiAsymRM1(a, y0, m);
    r *= _fermi_rm1[i][NFRM1];
  } else if (y > _fermi_rm1[0][NFRM1-1]) {
    if (y0 < 0) return -1;
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

void PrepFermiNR() {
  int i, ie;
  double du, ru, u, y;
  double umin = -15.0, umid = 5.0;
  du = (umid-umin)/(NFRM3-1);
  ru = 1 + du/umid;
  u = umin;
  for (i = 0; i < NFRM4; i++) {
    if (i < NFRM3) {
      _fermi_nr[0][i] = u;
    } else {
      _fermi_nr[0][i] = log(u);
    }
    FERMID(0.5, u, EPS10, &y, &ie);
    _fermi_nr[1][i] = log(y);
    FERMID(1.0, u, EPS10, &y, &ie);
    _fermi_nr[2][i] = log(y);
    FERMID(1.5, u, EPS10, &y, &ie);
    _fermi_nr[3][i] = log(y);
    _fermi_nr[4][i] = 2.5*_fermi_nr[1][i]-1.5*_fermi_nr[3][i];
    _fermi_nr[5][i] = 4*_fermi_nr[1][i]-3*_fermi_nr[2][i];
    _fermi_nr[6][i] = _fermi_nr[3][i]-_fermi_nr[1][i];
    _fermi_nr[7][i] = _fermi_nr[2][i]-_fermi_nr[1][i];
    /*
    printf("fermi: %d %g %g %g %g %g %g %g %g %g\n", i, u,
	   _fermi_nr[0][i],
	   _fermi_nr[1][i],
	   _fermi_nr[2][i],
	   _fermi_nr[3][i],
	   _fermi_nr[4][i],
	   _fermi_nr[5][i],
	   _fermi_nr[6][i],
	   _fermi_nr[7][i]);
    */
    if (i < NFRM3-1) {
      u += du;
    } else if (i == NFRM3-1) {
      _fermi_nr[0][i+1] = log(u);
      for (ie = 1; ie < 8; ie++) {
	_fermi_nr[ie][i+1] = _fermi_nr[ie][i];
      }
      i++;
      u *= ru;
    } else {
      u *= ru;
    }
  }
  _nfermi_nr[0] = NFRM3;
  for (ie = 1; ie < 8; ie++) {
    _nfermi_nr[ie] = NFRM4-1;
    for (i = NFRM4-1; i > NFRM3; i--) {
      if (_fermi_nr[ie][i-1] < _fermi_nr[ie][i]) {
	_nfermi_nr[ie] = i;
	break;
      }
    }
  }
}

double FreeEta0(double ne, double te) {
  double x, u, *au, *fd1;
  int one, three, nu;

  au = _fermi_nr[0];
  fd1 = _fermi_nr[1];
  x = log(ne/(0.126987272*te*sqrt(te)));
  if (x < fd1[0]) {
    return log(x);
  }
  if (x > fd1[_nfermi_nr[1]]) {
    return pow(1.3293403882*x, TWOTHIRD);
  }
  one = 1;
  three = 3;
  if (x <= fd1[NFRM3-1]) {
    nu = NFRM3;
    UVIP3P(three, nu, fd1, au, one, &x, &u);
  } else {
    nu = _nfermi_nr[1]-NFRM3;
    UVIP3P(three, nu, fd1+NFRM3, au+NFRM3, one, &x, &u);
    u = exp(u);
  }
  return u;
}

double FreeEta1(double ne, double ke) {
  double x, u, *au, *fd1;
  int one, three, nu;

  au = _fermi_nr[0];
  fd1 = _fermi_nr[4];
  x = log(14.46694050558*ne/(ke*sqrt(ke)));
  if (x < fd1[0]) {
    return log(x);
  }
  if (x > fd1[_nfermi_nr[4]]) {
    return exp(au[NFRM4-1]);
  }
  one = 1;
  three = 3;
  if (x <= fd1[NFRM3-1]) {
    nu = NFRM3;
    UVIP3P(three, nu, fd1, au, one, &x, &u);
  } else {
    nu = _nfermi_nr[4]-NFRM3;
    UVIP3P(three, nu, fd1+NFRM3, au+NFRM3, one, &x, &u);
    u = exp(u);
  }
  return u;
}

double FreeEta2(double ne, double ve) {
  double x, u, *au, *fd1;
  int one, three, nu;

  au = _fermi_nr[0];
  fd1 = _fermi_nr[5];
  x = log(32.0*ne/(ve*ve*ve));
  if (x < fd1[0]) {
    return log(x);
  }
  if (x > fd1[_nfermi_nr[5]]) {
    return au[NFRM4-1];
  }
  one = 1;
  three = 3;
  if (x <= fd1[NFRM3-1]) {
    nu = NFRM3;
    UVIP3P(three, nu, fd1, au, one, &x, &u);
  } else {
    nu = _nfermi_nr[5]-NFRM3;
    UVIP3P(three, nu, fd1+NFRM3, au+NFRM3, one, &x, &u);
    u = exp(u);
  }
  return u;
}

double FreeKe(double u, double te) {
  double *au, *fd1, y;
  int one, three, nu;

  au = _fermi_nr[0];
  fd1 = _fermi_nr[6];
  
  if (u < au[0]) {
    return 1.5*te;
  }
  if (u > exp(au[NFRM4-1])) {
    return 0.6*u*te;
  }
  one = 1;
  three = 3;
  if (u <= au[NFRM3-1]) {
    nu = NFRM3;
    UVIP3P(three, nu, au, fd1, one, &u, &y);
  } else {
    nu = _nfermi_nr[6]-NFRM3;
    u = log(u);
    UVIP3P(three, nu, au+NFRM3, fd1+NFRM3, one, &u, &y);
  }
  return 1.5*exp(y)*te;
}

double FreeVe(double u, double te) {
  double *au, *fd1, y;
  int one, three, nu;

  au = _fermi_nr[0];
  fd1 = _fermi_nr[7];

  if (u < au[0]) {
    return 1.5957691216*sqrt(te);
  }
  if (u > exp(au[NFRM4-1])) {
    return 1.0606601718*sqrt(te*u);
  }
  one = 1;
  three = 3;
  if (u <= au[NFRM3-1]) {
    nu = NFRM3;
    UVIP3P(three, nu, au, fd1, one, &u, &y);
  } else {
    u = log(u);
    nu = _nfermi_nr[7]-NFRM3;
    UVIP3P(three, nu, au+NFRM3, fd1+NFRM3, one, &u, &y);
  }
  return 1.5957691216*exp(y)*sqrt(te);
}

double FreeTe(double u, double ne) {
  double *au, *fd1, y, t;
  int one, three, nu;

  au = _fermi_nr[0];
  fd1 = _fermi_nr[1];
  
  if (u < au[0]) {
    t = TWOTHIRD*(log(ne/2)-u);
    if (t > 100) t = 100;
    return TWO_PI*exp(t);
  }

  if (u > exp(au[NFRM4-1])) {
    t = 4.78539*pow(ne,TWOTHIRD)/u;
    return t;
  }

  one = 1;
  three = 3;
  if (u <= au[NFRM3-1]) {
    nu = NFRM3;
    UVIP3P(three, nu, au, fd1, one, &u, &y);
  } else {
    nu = _nfermi_nr[1]-NFRM3;
    u = log(u);
    UVIP3P(three, nu, au+NFRM3, fd1+NFRM3, one, &u, &y);
  }
  t = TWO_PI*pow(0.5*ne/exp(y), TWOTHIRD);
  return t;
}

double InterpFermiNR(int m, double u) {
  if (u < _fermi_nr[0][0]) {
    if (m < 4) {
      return exp(u);
    } else {
      return 1.0;
    }
  }
  if (u > exp(_fermi_nr[0][NFRM4-1])) {
    switch (m) {
    case 1:
      return 0.7522527780636751*u*sqrt(u);
    case 2:
      return 0.5*u*u;
    case 3:
      return 0.30090111122547003*u*u*sqrt(u);
    case 4:
      return 2.9735401935879517;
    case 5:
      return 2.561799803697627;
    case 6:
      return 0.4*u;
    case 7:
      return 0.6646701940895685*sqrt(u);
    default:
      return -2.0;
    }
  }

  double r;
  int one, three;
  one = 1;
  three = 3;
  if (u <= _fermi_nr[0][NFRM3-1]) {
    UVIP3P(three, NFRM3, _fermi_nr[0], _fermi_nr[m], one, &u, &r);
  } else {
    u = log(u);
    UVIP3P(three, _nfermi_nr[m]-NFRM3, _fermi_nr[0]+NFRM3,
	   _fermi_nr[m]+NFRM3, one, &u, &r);
  }
  return exp(r);
}

void PrepFermiRM1(double a, double fa, double t) {
  double g;
  if (_relativistic_fermi) {
    g = t*FINE_STRUCTURE_CONST2;
  } else {
    g = 0.0;
  }
  double aa = fabs(a);
  double y0 = log(Min(_fermi_a0, aa)*_fermi_y0);
  double y1 = log(Max(_fermi_a1, aa)*_fermi_y1);
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
  y1 = exp(y1);
  _fermi_ymx = y1;
  double y2 = y1 + 0.5*g*y1*y1;
  _fermi_rmx = FermiAsymRM1(a, y2, 0)/FermiAsymRM1(a, y1, 0);  
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

void SetOptionOrbital(char *s, char *sp, int ip, double dp) {
  if (0 == strcmp(s, "orbital:relativistic_fermi")) {
    _relativistic_fermi0 = ip;
    return;
  }
  if (0 == strcmp(s, "orbital:relativistic_xfermi")) {
    _relativistic_xfermi = dp;
    return;
  }
  if (0 == strcmp(s, "orbital:fermi_a0")) {
    _fermi_a0 = dp;
    return;
  }
  if (0 == strcmp(s, "orbital:fermi_a1")) {
    _fermi_a1 = dp;
    return;
  }
  if (0 == strcmp(s, "orbital:fermi_y0")) {
    _fermi_y0 = dp;
    return;
  }
  if (0 == strcmp(s, "orbital:fermi_y1")) {
    _fermi_y1 = dp;
    return;
  }
  if (0 == strcmp(s, "orbital:fermi_tmin")) {
    _fermi_tmin = dp;
    return;
  }
  if (0 == strcmp(s, "orbital:fermi_tmax")) {
    _fermi_tmax = dp;
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
  if (0 == strcmp(s, "orbital:fermi_rmf")) {
    strncpy(_fermi_rmf, sp, 256);
    return;
  }
  if (0 == strcmp(s, "orbital:ionsph_ifermi")) {
    _ionsph_ifermi = dp;
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
  if (0 == strcmp(s, "orbital:sp_mode")) {
    _sp_mode = ip;
    return;
  }
  if (0 == strcmp(s, "orbital:debye_mode")) {
    _debye_mode = ip;
    return;
  }
  if (0 == strcmp(s, "orbital:sp_print")) {
    _sp_print = ip%10;
    _sp_piter0 = ip/10;
    _sp_piter1 = _sp_piter0%10000;
    _sp_piter0 = _sp_piter0/10000;
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
  if (0 == strcmp(s, "orbital:ionsph_miniter")) {
    _ionsph_miniter = ip;
    return;
  }  
  if (0 == strcmp(s, "orbital:ionsph_maxiter")) {
    _ionsph_maxiter = ip;
    return;
  }  
  if (0 == strcmp(s, "orbital:ionsph_bmode")) {
    _ionsph_bmode = ip;
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
  if (0 == strcmp(s, "orbital:on_error")) {
    _on_error = ip;
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
  if (0 == strcmp(s, "orbital:zcoll")) {
    SetZColl(dp);
    return;
  }
  if (0 == strcmp(s, "orbital:mcoll")) {
    SetMColl(dp);
    return;
  }
  if (0 == strcmp(s, "orbital:emin_amp")) {
    _emin_amp = dp;
    return;
  }
  if (0 == strcmp(s, "orbital:sturm_rmx")) {
    _sturm_rmx = dp;
    return;
  }
  if (0 == strcmp(s, "orbital:maxiter")) {
    max_iteration = ip;
    return;
  }
  if (0 == strcmp(s, "orbital:enerelerr")) {
    _enerelerr = dp;
    return;
  }
  if (0 == strcmp(s, "orbital:enerelerr1")) {
    _enerelerr1 = dp;
    return;
  }
  if (0 == strcmp(s, "orbital:eneabserr")) {
    _eneabserr = dp;
    return;
  }
  if (0 == strcmp(s, "orbital:sc_bqp")) {
    _sc_bqp = dp;
    return;
  }
  if (0 == strcmp(s, "orbital:sc_rbf")) {
    _sc_rbf = dp;
    return;
  }
  if (0 == strcmp(s, "orbital:sc_rsf")) {
    _sc_rsf = dp;
    return;
  }
  if (0 == strcmp(s, "orbital:sc_ewr")) {
    _sc_ewr = dp;
    return;
  }
  if (0 == strcmp(s, "orbital:sc_ewf")) {
    _sc_ewf = dp;
    return;
  }
  if (0 == strcmp(s, "orbital:sc_ewm")) {
    _sc_ewm = ip;
    return;
  }
  if (0 == strcmp(s, "orbital:debug")) {
    _debug = ip;
    return;
  }
  if (0 == strcmp(s, "orbital:matchtol")) {
    _matchtol = dp;
    return;
  }
  if (0 == strcmp(s, "orbital:veff_corr")) {
    _veff_corr = ip;
    return;
  }
  if (0 == strcmp(s, "orbital:pwa")) {
    _pwa = ip;
    return;
  }
  if (0 == strcmp(s, "orbital:rydnorm")) {
    _rydnorm = ip;
    return;
  }
}

double SetSPZW(int n, double *zw) {
  if (_sp_nzs > 0) {
    free(_sp_zs);
    free(_sp_zw);
    _sp_zs = NULL;
    _sp_zw = NULL;
    _sp_zu = 0.0;
    _sp_nzs = 0;
  }
  if (n <= 0) return -1;
  int n2 = n/2;
  _sp_nzs = n2;
  _sp_zs = malloc(sizeof(double)*n2);
  _sp_zw = malloc(sizeof(double)*n2);
  int i, k;
  double w = 0, z1 = 0, z2 = 0;
  for (i = 0; i < n2; i++) {
    k = 2*i;
    _sp_zs[i] = zw[k];
    _sp_zw[i] = zw[k+1];
    w += _sp_zw[i];
    z1 += _sp_zs[i]*_sp_zw[i];
    z2 += _sp_zs[i]*_sp_zs[i]*_sp_zw[i];
  }
  for (i = 0; i < n2; i++) {
    _sp_zw[i] /= w;
  }
  z1 /= w;
  z2 /= w;
  if (z1) {
    _sp_zu = z2/z1;
  } else {
    _sp_zu = 0.0;
  }
}
