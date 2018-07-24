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

#include "radial.h"
#include "mpiutil.h"
#include "cf77.h"
#include "structure.h"
#include <errno.h>

static char *rcsid="$Id$";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

typedef struct _FLTARY_ {
  short npts;
  float *yk;
} FLTARY;

static POTENTIAL *potential;
static POTENTIAL *hpotential;
static POTENTIAL *rpotential;

#define Large(orb) ((orb)->wfun)
#define Small(orb) ((orb)->wfun + potential->maxrp)

static ARRAY *orbitals;
static int n_orbitals;
static int n_continua;
 
static double _dwork[MAXRP];
static double _dwork1[MAXRP];
static double _dwork2[MAXRP];
static double _dwork3[MAXRP];
static double _dwork4[MAXRP];
static double _dwork5[MAXRP];
static double _dwork6[MAXRP];
static double _dwork7[MAXRP];
static double _dwork8[MAXRP];
static double _dwork9[MAXRP];
static double _dwork10[MAXRP];
static double _dwork11[MAXRP];
static double _dwork12[MAXRP];
static double _dwork13[MAXRP];
static double _dwork14[MAXRP];
static double _dwork15[MAXRP];
static double _dwork16[MAXRP];
static double _phase[MAXRP];
static double _dphase[MAXRP];
static double _dphasep[MAXRP];
static double _yk[MAXRP];
static double _zk[MAXRP];
static double _xk[MAXRP];

#pragma omp threadprivate(potential,hpotential,rpotential,_dwork,_dwork1,_dwork2,_dwork3,_dwork4,_dwork5,_dwork6,_dwork7,_dwork8,_dwork9,_dwork10,_dwork11,_dwork12,_dwork13,_dwork14,_dwork15,_dwork16,_phase,_dphase,_dphasep,_yk,_zk,_xk)

static struct {
  double stabilizer;
  double tolerance; /* tolerance for self-consistency */
  int maxiter; /* max iter. for self-consistency */
  double screened_charge; 
  int screened_kl;
  int n_screen;
  int *screened_n;
  int iprint; /* printing infomation in each iteration. */
  int iset;
  int mce;
} optimize_control = {OPTSTABLE, OPTTOL, OPTNITER, 
		      1.0, 1, 0, NULL, OPTPRINT, 0, 0};

static struct {
  int kl0;
  int kl1;
} slater_cut = {1000000, 1000000};

static struct {
  int se;
  int mse;
  int sse;
  int pse;
  int vp;
  int nms;
  int sms;
  int br;
  int mbr;
  int nbr;
  double xbr;
  int minbr;
  double ose0, ose1, ase;
  double cse0, cse1, ise;
} qed = {QEDSE, QEDMSE, 0, 0, QEDVP, QEDNMS, QEDSMS,
	 QEDBREIT, QEDMBREIT, QEDNBREIT, 0.01, 0,
	 1.0, 1.0, 3.0,
	 0.07, 0.05, 1.5};

static int _orthogonalize_mode = 1;
static int _korbmap = KORBMAP;
static int _norbmap0 = NORBMAP0;
static int _norbmap1 = NORBMAP1;
static int _norbmap2 = NORBMAP2;
static ORBMAP *_orbmap = NULL;

static AVERAGE_CONFIG average_config = {0, 0, NULL, NULL, NULL, 0, NULL, NULL};

static MULTI *slater_array;
static MULTI *xbreit_array[5];
static MULTI *wbreit_array;
static MULTI *breit_array;
static MULTI *vinti_array;
static MULTI *qed1e_array;
static MULTI *residual_array;
static MULTI *multipole_array; 
static MULTI *moments_array;
static MULTI *gos_array;
static MULTI *yk_array;

static int n_awgrid = 0;
static double awgrid[MAXNTE];

static double PhaseRDependent(double x, double eta, double b);

#ifdef PERFORM_STATISTICS
static RAD_TIMING rad_timing = {0, 0, 0, 0};
int GetRadTiming(RAD_TIMING *t) {
  memcpy(t, &rad_timing, sizeof(RAD_TIMING));
  return 0;
}
#endif

double *WorkSpace(int i) {
  switch (i) {
  case 0:
    return _dwork;
  case 1:
    return _dwork1;
  case 2:
    return _dwork2;
  case 3:
    return _dwork3;
  case 4:
    return _dwork4;
  case 5:
    return _dwork5;
  case 6:
    return _dwork6;
  case 7:
    return _dwork7;
  case 8:
    return _dwork8;
  case 9:
    return _dwork9;
  case 10:
    return _dwork10;
  case 11:
    return _dwork11;
  case 12:
    return _dwork12;
  case 13:
    return _dwork13;
  case 14:
    return _dwork14;
  case 15:
    return _dwork15;
  case 16:
    return _dwork16;
  default:
    return NULL;
  }  
}

void SetOrbMap(int k, int n0, int n1, int n2) {
  if (k <= 0) k = KORBMAP;
  if (n0 <= 0) n0 = NORBMAP0;
  if (n1 <= 0) n1 = NORBMAP1;
  if (n2 <= 0) n2 = NORBMAP2;
  int i, j;
  if (_orbmap != NULL) {
    for (i = 0; i < _korbmap; i++) {
      free(_orbmap[i].opn);
      free(_orbmap[i].onn);
      free(_orbmap[i].ozn);
    }
    free(_orbmap);
  }

  _korbmap = k;
  _norbmap0 = n0;
  _norbmap1 = n1;
  _norbmap2 = n2;  
  _orbmap = malloc(sizeof(ORBMAP)*_korbmap);  
  for (i = 0; i < _korbmap; i++) {
    _orbmap[i].nzn = 0;
    _orbmap[i].opn = malloc(sizeof(ORBITAL *)*_norbmap0);
    _orbmap[i].onn = malloc(sizeof(ORBITAL *)*_norbmap1);
    _orbmap[i].ozn = malloc(sizeof(ORBITAL *)*_norbmap2);
    for (j = 0; j < _norbmap0; j++) {
      _orbmap[i].opn[j] = NULL;
    }
    for (j = 0; j < _norbmap1; j++) {
      _orbmap[i].onn[j] = NULL;
    }
    for (j = 0; j < _norbmap2; j++) {
      _orbmap[i].ozn[j] = NULL;
    }
  }
}

int RestorePotential(char *fn, POTENTIAL *p) {
  BFILE *f;
  int n, i, k;

  f = BFileOpen(fn, "r", -1);
  if (f == NULL) {
    MPrintf(0, "cannot open potential file: %s\n", fn);
    return -1;
  }

  n = BFileRead(&p->mode, sizeof(int), 1, f);
  n = BFileRead(&p->flag, sizeof(int), 1, f);
  n = BFileRead(&p->r_core, sizeof(int), 1, f);
  n = BFileRead(&p->nmax, sizeof(int), 1, f);
  n = BFileRead(&p->maxrp, sizeof(int), 1, f);
  n = BFileRead(&p->hxs, sizeof(double), 1, f);
  n = BFileRead(&p->ahx, sizeof(double), 1, f);
  n = BFileRead(&p->ihx, sizeof(double), 1, f);
  n = BFileRead(&p->rhx, sizeof(double), 1, f);
  n = BFileRead(&p->dhx, sizeof(double), 1, f);
  n = BFileRead(&p->chx, sizeof(double), 1, f);
  n = BFileRead(&p->hx0, sizeof(double), 1, f);
  n = BFileRead(&p->hx1, sizeof(double), 1, f);
  n = BFileRead(&p->ratio, sizeof(double), 1, f);
  n = BFileRead(&p->asymp, sizeof(double), 1, f);
  n = BFileRead(&p->rmin, sizeof(double), 1, f);
  n = BFileRead(&p->N, sizeof(double), 1, f);
  n = BFileRead(&p->lambda, sizeof(double), 1, f);
  n = BFileRead(&p->a, sizeof(double), 1, f);
  n = BFileRead(&p->ar, sizeof(double), 1, f);
  n = BFileRead(&p->br, sizeof(double), 1, f);
  n = BFileRead(&p->ib, sizeof(int), 1, f);
  n = BFileRead(&p->nb, sizeof(int), 1, f);
  n = BFileRead(&p->ib1, sizeof(int), 1, f);
  n = BFileRead(&p->bqp, sizeof(double), 1, f);
  n = BFileRead(&p->rb, sizeof(double), 1, f);
  n = BFileRead(p->Z, sizeof(double), p->maxrp, f);
  n = BFileRead(p->dZ, sizeof(double), p->maxrp, f);
  n = BFileRead(p->dZ2, sizeof(double), p->maxrp, f);
  n = BFileRead(p->rad, sizeof(double), p->maxrp, f);
  n = BFileRead(p->mqrho, sizeof(double), p->maxrp, f);
  n = BFileRead(p->dr_drho, sizeof(double), p->maxrp, f);
  n = BFileRead(p->dr_drho2, sizeof(double), p->maxrp, f);
  n = BFileRead(p->vtr, sizeof(double), p->maxrp, f);
  n = BFileRead(p->Vc, sizeof(double), p->maxrp, f);
  n = BFileRead(p->dVc, sizeof(double), p->maxrp, f);
  n = BFileRead(p->dVc2, sizeof(double), p->maxrp, f);
  n = BFileRead(p->qdist, sizeof(double), p->maxrp, f);
  n = BFileRead(p->U, sizeof(double), p->maxrp, f);
  n = BFileRead(p->dU, sizeof(double), p->maxrp, f);
  n = BFileRead(p->dU2, sizeof(double), p->maxrp, f);
  n = BFileRead(p->W, sizeof(double), p->maxrp, f);
  n = BFileRead(p->dW, sizeof(double), p->maxrp, f);
  n = BFileRead(p->dW2, sizeof(double), p->maxrp, f);
  for (i = p->maxrp; i < MAXRP; i++) {
    p->Z[i] = 0;
    p->dZ[i] = 0;
    p->dZ2[i] = 0;
    p->rad[i] = 0;
    p->dr_drho[i] = 0;
    p->dr_drho2[i] = 0;
    p->Vc[i] = 0;
    p->dVc[i] = 0;
    p->dVc2[i] = 0;
    p->qdist[i] = 0;
    p->qdist[i] = 0;
    p->U[i] = 0;
    p->dU[i] = 0;
    p->dU2[i] = 0;
    p->W[i] = 0;
    p->dW[i] = 0;
    p->dW2[i] = 0;
    p->ZVP[i] = 0;
    p->dZVP[i] = 0;
    p->dZVP2[i] = 0;
    for (k = 0; k < NKSEP; k++) {
      p->ZSE[k][i] = 0;
      p->dZSE[k][i] = 0;
      p->dZSE2[k][i] = 0;
    }
  }
  p->atom = GetAtomicNucleus();
  BFileClose(f);
  ReinitRadial(1);
  SetReferencePotential(hpotential, p, 1);
  SetReferencePotential(rpotential, p, 0);
  CopyPotentialOMP(0);
  return 0;
}

void SetReferencePotential(POTENTIAL *h, POTENTIAL *p, int hlike) {
  int i, k;
  memcpy(h, p, sizeof(POTENTIAL));  
  if (hlike) {
    h->N = 1;
    h->a = 0;
    h->lambda = 0;
    h->hlike = 1;
  } else {
    h->hlike = 0;
  }
  h->nse = 0;
  h->mse = 0;
  h->pse = 0;
  h->mvp = 0;
  h->pvp = 0;
  for (i = 0; i < p->maxrp; i++) {
    if (hlike) {
      h->Vc[i] = 0;
      h->dVc[i] = 0;
      h->dVc2[i] = 0;
      h->qdist[i] = 0;
      h->U[i] = 0;
      h->dU[i] = 0;
      h->dU2[i] = 0;
      h->W[i] = 0;
      h->dW[i] = 0;
      h->dW2[i] = 0;
      if (h->atom->epm >= 0) {
	h->Z[i] = GetAtomicEffectiveZ(h->rad[i]);
      }
    }
    h->ZVP[i] = 0;
    h->dZVP[i] = 0;
    h->dZVP2[i] = 0;
    for (k = 0; k < NKSEP; k++) {
      h->ZSE[k][i] = 0;
      h->dZSE[k][i] = 0;
      h->dZSE2[k][i] = 0;
    }
  }
  if (hlike) {
    if (h->atom->epm >= 0) {
      Differential(h->Z, h->dZ, 0, h->maxrp-1, h->dr_drho);
      Differential(h->dZ, h->dZ2, 0, h->maxrp-1, h->dr_drho);
    }
    SetPotentialVc(h);
  }
  SetPotentialVT(h);
}

int SavePotential(char *fn, POTENTIAL *p) {
  FILE *f;
  int n, i;

  if (MyRankMPI() != 0) return 0;
  
  f = fopen(fn, "w");
  if (f == NULL) {
    MPrintf(0, "cannot open potential file: %s\n", fn);
    return -1;
  }

  n = fwrite(&p->mode, sizeof(int), 1, f);
  n = fwrite(&p->flag, sizeof(int), 1, f);
  n = fwrite(&p->r_core, sizeof(int), 1, f);
  n = fwrite(&p->nmax, sizeof(int), 1, f);
  n = fwrite(&p->maxrp, sizeof(int), 1, f);
  n = fwrite(&p->hxs, sizeof(double), 1, f);
  n = fwrite(&p->ahx, sizeof(double), 1, f);
  n = fwrite(&p->ihx, sizeof(double), 1, f);
  n = fwrite(&p->rhx, sizeof(double), 1, f);
  n = fwrite(&p->dhx, sizeof(double), 1, f);
  n = fwrite(&p->chx, sizeof(double), 1, f);
  n = fwrite(&p->hx0, sizeof(double), 1, f);
  n = fwrite(&p->hx1, sizeof(double), 1, f);
  n = fwrite(&p->ratio, sizeof(double), 1, f);
  n = fwrite(&p->asymp, sizeof(double), 1, f);
  n = fwrite(&p->rmin, sizeof(double), 1, f);
  n = fwrite(&p->N, sizeof(double), 1, f);
  n = fwrite(&p->lambda, sizeof(double), 1, f);
  n = fwrite(&p->a, sizeof(double), 1, f);
  n = fwrite(&p->ar, sizeof(double), 1, f);
  n = fwrite(&p->br, sizeof(double), 1, f);
  n = fwrite(&p->ib, sizeof(int), 1, f);
  n = fwrite(&p->nb, sizeof(int), 1, f);
  n = fwrite(&p->ib1, sizeof(int), 1, f);
  n = fwrite(&p->bqp, sizeof(double), 1, f);
  n = fwrite(&p->rb, sizeof(double), 1, f);
  n = fwrite(p->Z, sizeof(double), p->maxrp, f);
  n = fwrite(p->dZ, sizeof(double), p->maxrp, f);
  n = fwrite(p->dZ2, sizeof(double), p->maxrp, f);
  n = fwrite(p->rad, sizeof(double), p->maxrp, f);
  n = fwrite(p->mqrho, sizeof(double), p->maxrp, f);
  n = fwrite(p->dr_drho, sizeof(double), p->maxrp, f);
  n = fwrite(p->dr_drho2, sizeof(double), p->maxrp, f);
  n = fwrite(p->vtr, sizeof(double), p->maxrp, f);
  n = fwrite(p->Vc, sizeof(double), p->maxrp, f);
  n = fwrite(p->dVc, sizeof(double), p->maxrp, f);
  n = fwrite(p->dVc2, sizeof(double), p->maxrp, f);
  n = fwrite(p->qdist, sizeof(double), p->maxrp, f);
  n = fwrite(p->U, sizeof(double), p->maxrp, f);
  n = fwrite(p->dU, sizeof(double), p->maxrp, f);
  n = fwrite(p->dU2, sizeof(double), p->maxrp, f);
  n = fwrite(p->W, sizeof(double), p->maxrp, f);
  n = fwrite(p->dW, sizeof(double), p->maxrp, f);
  n = fwrite(p->dW2, sizeof(double), p->maxrp, f);
  fclose(f);
  return 0;
}

int ModifyPotential(char *fn, POTENTIAL *p) {
  BFILE *f;
  int n, i, k, np;
  char buf[BUFLN];
  char *c;

  f = BFileOpen(fn, "r", -1);
  if (f == NULL) {
    MPrintf(0, "cannot open potential file: %s\n", fn);
    return -1;
  }
  i = 0;
  while (1) {
    if (NULL == BFileGetLine(buf, BUFLN, f)) break;
    buf[BUFLN-1] = '\0';
    c = buf;
    while (*c && (*c == ' ' || *c == '\t')) c++;
    if (*c == '\0' || *c == '#') continue;
    if (i >= p->maxrp) {
      MPrintf(0, "potential file exceeds max grid points: %d %d\n",
	      i, p->maxrp);
      break;
    }
    k = sscanf(buf, "%lg %lg", _dwork12+i, _dwork13+i);
    if (k != 2) continue;    
    i++;
  }
  BFileClose(f);
  n = i;  
  for (i = 0; i < n; i++) {
    _dwork13[i] *= _dwork12[i];
    _dwork12[i] = log(_dwork12[i]);
  }
  for (i = 0; i < p->maxrp; i++) {
    p->W[i] = log(p->rad[i]);
  }
  np = 3;
  UVIP3P(np, n, _dwork12, _dwork13, p->maxrp, p->W, p->dW);
  for (i = 0; i < p->maxrp; i++) {
    p->U[i] = p->dW[i]/p->rad[i] - p->Vc[i];
  }
  SetPotentialU(p, p->maxrp, NULL);
  ReinitRadial(1);
  return 0;
}

static void InitFltAryData(void *p, int n) {
  FLTARY *d;
  int i;

  d = (FLTARY *) p;
  for (i = 0; i < n; i++) {
    d[i].npts = -1;
    d[i].yk = NULL;
  }
}

int FreeSimpleArray(MULTI *ma) {
  MultiFreeData(ma, NULL);
  return 0;
}

int FreeSlaterArray(void) {
  FreeSimpleArray(slater_array);
  return 0;
}

int FreeResidualArray(void) {
  return FreeSimpleArray(residual_array);
}

static void FreeMultipole(void *p) {
  double *dp;
  dp = *((double **) p);
  free(dp);
  *((double **) p) = NULL;
}

static void FreeFltAryData(void *p) {
  FLTARY *dp;
  
  dp = (FLTARY *) p;
  if (dp->npts > 0) {
    free(dp->yk);    
    dp->yk = NULL;
    dp->npts = -1;
  }
}

int FreeMultipoleArray(void) {
  MultiFreeData(multipole_array, FreeMultipole);
  return 0;
}

int FreeMomentsArray(void) {
  MultiFreeData(moments_array, NULL);
  return 0;
}

int FreeGOSArray(void) {
  MultiFreeData(gos_array, FreeMultipole);
  return 0;
}

int FreeVintiArray(void) {
  MultiFreeData(vinti_array, FreeMultipole);
  return 0;
}

int FreeBreitArray(void) {
  int i;
  FreeSimpleArray(breit_array);
  FreeSimpleArray(wbreit_array);
  for (i = 0; i < 5; i++) {
    MultiFreeData(xbreit_array[i], FreeFltAryData);
  }
  return 0;
}

int FreeYkArray(void) {
  MultiFreeData(yk_array, FreeFltAryData);
  return 0;
}
  
double *WLarge(ORBITAL *orb) {
  return Large(orb);
}

double *WSmall(ORBITAL *orb) {
  return Small(orb);
}

POTENTIAL *RadialPotential(void) {
  return potential;
}

void SetSlaterCut(int k0, int k1) {
  if (k0 > 0) {
    slater_cut.kl0 = 2*k0;
  } else {
    slater_cut.kl0 = 1000000;
  } 
  if (k1 > 0) {
    slater_cut.kl1 = 2*k1;
  } else {
    slater_cut.kl1 = 1000000;
  }
}

void SolvePseudo(int kmin, int kmax, int nb, int nmax, int nd, double xdf) {
  int k, n, i, k2, j, na, ka, kb, kn, maxrp2;
  double *p, *q, pn, qn, a, e0;
  ORBITAL *orb0, *orb;
  
  if (nb <= 0) nb = potential->nb;
  if (nb <= 0) {
    printf("invalid nb in SolvePseudoOrbitals: %d\n", nb);
    return;
  }
  potential->npseudo = nb;
  potential->mpseudo = nmax;
  if (nd == -1) nd = nb;
  else if (nd == -2) nd = nmax;
  else if (nd == 0) nd = 1;
  potential->dpseudo = nd;
  maxrp2 = 2*potential->maxrp;
  ResetWidMPI();
#pragma omp parallel default(shared) private(k, k2, na, j, ka, kb, orb0, n, kn, orb, p, q, a, e0, pn, qn)
  {
  for (k = kmin; k <= kmax; k++) {
    k2 = 2*k;
    na = k + nd;
    if (na < nb) na = nb;
    if (na > nmax) na = nmax;
    for (j = k2-1; j <= k2+1; j += 2) {
      if (j < 0) continue;
      int skip = SkipMPI();
      if (skip) continue;
      double wt0 = WallTime();
      ka = GetKappaFromJL(j, k2);
      for (n = k+1; n <= na; n++) {
	kb = OrbitalIndex(n, ka, 0);
      }
      orb0 = GetOrbitalSolved(kb);
      double wt1 = WallTime();
      for (n = na+1; n <= nmax; n++) {
	kn = OrbitalIndex(n, ka, 0);
	orb = GetOrbital(kn);
	if (orb->isol <= 0) {
	  orb->wfun = malloc(sizeof(double)*maxrp2);	  
	}
	if (orb->isol != 3) {
	  orb->isol = 3;
	  memcpy(orb->wfun, orb0->wfun, sizeof(double)*maxrp2);
	  p = Large(orb);
	  q = Small(orb);
	  orb->bqp0 = orb0->bqp0;
	  orb->bqp1 = orb0->bqp1;
	  //orb->im = orb0->im;
	  orb->qr_norm = orb0->qr_norm;
	  orb->ilast = orb0->ilast;
	  orb->kv = orb0->kv;
	  for (i = 0; i <= orb->ilast; i++) {
	    a = potential->rad[i]/potential->rad[orb->ilast];
	    p[i] *= sin(TWO_PI*a);
	  }
	  e0 = orb0->energy;
	  orb->energy = e0;
	  DiracSmall(orb, potential, -1, orb->kv);
	  pn = InnerProduct(0, orb->ilast, p, p, potential);
	  qn = InnerProduct(0, orb->ilast, q, q, potential);
	  a = 1.0/sqrt(pn+qn);
	  for (i = 0; i <= orb->ilast; i++) {
	    p[i] *= a;
	    q[i] *= a;
	  }	  
	  Orthogonalize(orb);
	}
	orb0 = orb;
      }
      double wt2 = WallTime();
      if (xdf >= 0) SolveDFKappa(ka, nmax, xdf);
      double wt3 = WallTime();
      MPrintf(-1, "SolvePseudo: %3d %3d %2d %g %11.4E %11.4E %11.4E\n",
	      na, nmax, ka, xdf, wt1-wt0, wt2-wt1, wt3-wt2);
    }
  }
  }
}

void SolveDFKappa(int ka, int nmax, double xdf) {
  int j, k, k2, n, i, isym, m, s, t, nn, maxrp2;
  int ig, ic;
  double *mix;
  ORBITAL **orb;
  double **wf, r, r0, ec;
  HAMILTON *h;  
  AVERAGE_CONFIG *acfg;
  CONFIG *c;
  CONFIG_GROUP *g;

  acfg = &average_config;
  
  GetJLFromKappa(ka, &j, &k2);
  k = k2/2;
  maxrp2 = potential->maxrp*2;
  nn = nmax-k;
  orb = malloc(sizeof(ORBITAL *)*nn);
  wf = malloc(sizeof(double *)*nn);
  isym = j*2;
  if (IsOdd(k)) isym++;
  h = GetHamilton(isym);
  AllocHamMem(h, nn, nn);
  i = 0;
  for (n = k+1; n <= nmax; n++,i++) {
    m = OrbitalIndex(n, ka, 0);
    orb[i] = GetOrbitalSolved(m);
    wf[i] = malloc(sizeof(double)*maxrp2);
    memcpy(wf[i], orb[i]->wfun, sizeof(double)*maxrp2);
    h->basis[i] = n;
  }
  for (j = 0; j < nn; j++) {
    t = j*(j+1)/2;
    for (i = 0; i <= j; i++) {
      r0 = ZerothHamilton(orb[i], orb[j]);
      r = 0;
      if (xdf >= 0) {
	for (ig = 0; ig < acfg->ng; ig++) {
	  g = GetGroup(acfg->kg[ig]);
	  double tec = 0;
	  for (ic = 0; ic < g->n_cfgs; ic++) {
	    c = GetConfigFromGroup(acfg->kg[ig], ic);
	    ec = ConfigHamilton(c, orb[i], orb[j], xdf);
	    if (i == j) {
	      ec += AverageEnergyConfig(c);
	    }
	    tec += ec;
	  }
	  if (g->n_cfgs > 0) {
	    r += acfg->weight[ig]*tec/g->n_cfgs;
	  }
	}
      }
      h->hamilton[i+t] = r0 + r;
      //printf("dh: %d %d %12.5E %12.5E %12.5E\n", orb[i]->n, orb[j]->n, r0, r, h->hamilton[i+t]);
    }
  }
  
  if (0 > DiagnolizeHamilton(h)) {
    MPrintf(-1, "Diag DF Ham error: %d %d\n", ka, nmax);
    Abort(1);
  }
  
  mix = h->mixing + h->dim;
  for (i = 0; i < h->dim; i++) {
    orb[i]->energy = h->mixing[i];
    for (s = 0; s < maxrp2; s++) {
      orb[i]->wfun[s] = 0;
    }

    if (mix[i] < 0) {
      for (t = 0; t < h->dim; t++) {
	mix[t] = -mix[t];
      }
    }

    for (t = 0; t < h->dim; t++) {
      for (s = 0; s < maxrp2; s++) {
	orb[i]->wfun[s] += wf[t][s]*mix[t];
      }
    }
    mix += h->dim;
  }
  for (i = 0; i < h->dim; i++) {
    free(wf[i]);
  }
  free(wf);
  free(orb);
  AllocHamMem(h, -1, -1);
  AllocHamMem(h, 0, 0);
  ReinitRadial(2);
}

void SetPotentialMode(int m, double h, double ih, double h0, double h1) {
  potential->mode = m;
  if (h > 1e10) {
    potential->hxs = POTHXS;
  } else {
    potential->hxs = h;
  }
  if (ih > 1e10) {
    potential->ihx = POTIHX;
  } else {
    potential->ihx = ih;
  }
  if (h0 >= 0) {
    potential->hx0 = h0;
  } else {
    potential->hx0 = POTHX0;
  }
  if (h1 >= 0) {
    potential->hx1 = h1;
  } else {
    potential->hx1 = POTHX1;
  }
}

void PrintQED() {
  MPrintf(0, "SE: %d %d %d %d %d\n", qed.se, qed.mse,
	  potential->pse, qed.sse, qed.pse);
  MPrintf(0, "MSE: %g %g %g %g %g %g\n",
	  qed.ose0, qed.ose1, qed.ase, qed.cse0, qed.cse1, qed.ise);
  MPrintf(0, "VP: %d %d\n", qed.vp, potential->pvp);
  MPrintf(0, "MS: %d %d\n", qed.nms, qed.sms);
  MPrintf(0, "BR: %d %d %d %g %d\n",
	  qed.br, qed.mbr, qed.nbr, qed.xbr, qed.minbr);
}

int GetBoundary(double *rb, double *b, int *nmax, double *dr) {
  int ib;

  if (potential->ib > 0) {
    ib = potential->ib;
  } else {
    ib = potential->maxrp-1;
  }
  *rb = potential->rad[ib];
  *b = potential->bqp;
  *nmax = potential->nb;
  *dr = potential->rad[ib]-potential->rad[ib-1];
  return potential->ib;
}

int SetBoundaryMaster(int nmax, double p, double bqp) {
  ORBITAL *orb;
  int i, j, n, kl, kl2, kappa, k;
  double d1, d2, d;

  if (nmax == -100) {
    if (p <= 0.0) {
      printf("2nd argument must be > 0 in SetBoundary(-100,...)\n");
      return -1;
    }
    for (i = 0; i < potential->maxrp-10; i++) {
      if (potential->rad[i] >= p) break;
    }
    potential->ib1 = i;
    if (bqp > 0) {
      if (bqp >= p) {
	printf("3rd argument must be less than 2nd in SetBoundary(-100,...)\n");
	return -1;
      }
      for (i = 0; i < potential->maxrp; i++) {
	if (potential->rad[i] >= bqp) break;
      }
      if (i < potential->ib1) {
	potential->ib = i;
      }
    }
    return 0;
  }
  potential->nb = abs(nmax);
  potential->bqp = bqp;
  if (nmax == 0) {
    potential->ib = 0;
  } else if (nmax < 0) {
    d = GetResidualZ();
    d1 = potential->nb;
    if (p < 0.0) p = 2.0*d1*d1/d;
    for (i = 0; i < potential->maxrp; i++) {
      if (potential->rad[i] >= p) break;
    }
    if (i > potential->maxrp-10) {
      printf("enlarge maxrp\n");
      exit(1);
    }
    if (IsEven(i)) i++;
    potential->ib = i;
    for (n = 1; n <= potential->nb; n++) {
      for (kl = 0; kl < n; kl++) {
	kl2 = 2*kl;
	for (j = kl2 - 1; j <= kl2 + 1; j += 2) {
	  if (j < 0) continue;
	  kappa = GetKappaFromJL(j, kl2);
	  k = OrbitalIndex(n, kappa, 0);
	  orb = GetOrbitalSolved(k);
	}
      }
    }
  } else {
    if (p < 0.0) p = 1E-8;
    for (n = 1; n <= potential->nb; n++) {
      for (kl = 0; kl < n; kl++) {
	kl2 = 2*kl;
	for (j = kl2 - 1; j <= kl2 + 1; j += 2) {
	  if (j < 0) continue;
	  kappa = GetKappaFromJL(j, kl2);
	  k = OrbitalIndex(n, kappa, 0);
	  orb = GetOrbitalSolved(k);	  
	  for (i = orb->ilast-3; i >= 0; i--) {
	    d1 = Large(orb)[i];
	    d2 = Small(orb)[i];
	    d = d1*d1 + d2*d2;
	    if (d >= p) {
	      i++;
	      break;
	    }
	  }
	  if (potential->ib < i) potential->ib = i;
	}
      }
    }
    if (IsEven(i)) i++;
    if (i > potential->maxrp-10) {
      printf("enlarge maxrp\n");
      return -1;
    }
  }
  potential->ib1 = potential->ib;
  if (potential->ib > 0 && potential->ib < potential->maxrp) {
    potential->rb = potential->rad[potential->ib];
  }
  return 0;
}

int SetBoundary(int nmax, double p, double bqp) {
  int r = SetBoundaryMaster(nmax, p, bqp);
  CopyPotentialOMP(0);
  return r;
}

int RadialOverlaps(char *fn, int kappa) {
  ORBITAL *orb1, *orb2;
  int i, j, k;
  double r, b1, b2, a1, a2;
  FILE *f;

  f = fopen(fn, "w");
  for (k = 0; k < potential->maxrp; k++) {
    _dwork10[k] = 1.0;
  }
  for (i = 0; i < n_orbitals; i++) {
    orb1 = GetOrbital(i);
    if (orb1->kappa != kappa) continue;
    a1 = ZerothHamilton(orb1, orb1);
    for (j = 0; j < n_orbitals; j++) {
      orb2 = GetOrbital(j);
      if (orb2->kappa != kappa) continue;
      a2 = ZerothHamilton(orb2, orb2);
      Integrate(_dwork10, orb1, orb2, 1, &r, 0);
      b1 = ZerothHamilton(orb1, orb2);
      b2 = ZerothHamilton(orb2, orb1);
      fprintf(f, "%2d %2d %12.5E %12.5E  %2d %2d %12.5E %12.5E %12.5E %12.5E %12.5E\n",
	      orb1->n, orb1->kappa, orb1->energy, a1, 
	      orb2->n, orb2->kappa, orb2->energy, a2, r, b1, b2);
    }
  }
  fclose(f);

  return 0;
}
  
void SetSE(int n, int m, int s, int p) {
  qed.se = n;
  if (m >= 0) qed.mse = m%100;
  if (s >= 0) qed.sse = s;
  if (p >= 0) qed.pse = p;
  if (m >= 1000) {
    potential->hpvs = 1;
    m = m%1000;
  } else {
    potential->hpvs = 0;
  }
  if (m >= 140) {
    potential->pse = 1;
  } else {
    potential->pse = 0;
  }
  potential->mse = qed.mse;
  potential->nse = qed.se;
}

void SetModSE(double ose0, double ose1, double ase,
	      double cse0, double cse1, double ise) {
  qed.ose0 = ose0;
  qed.ose1 = ose1;
  qed.ase = ase;
  if (cse0 > 0) {
    qed.cse0 = cse0;
  }
  if (cse1 > 0) {
    qed.cse1 = cse1;
  }
  if (ise > 0) {
    qed.ise = ise;
  }
}

void SetVP(int n) {
  qed.vp = n;
  potential->mvp = qed.vp%100;
  potential->pvp = qed.vp > 100;
}

void SetBreit(int n, int m, int n0, double x0, int n1) {
  qed.br = n;
  if (m >= 0) qed.mbr = m;
  if (n0 >= 0) qed.nbr = n0;
  if (x0 > 0) qed.xbr = x0;
  if (n1 >= 0) qed.minbr = n1;
}

void SetMS(int nms, int sms) {
  qed.nms = nms;
  qed.sms = sms;
}

int SetAWGrid(int n, double awmin, double awmax) {
  int i;
  if (awmin < 1E-3) {
    awmin = 1E-3;
    awmax = awmax + 1E-3;
  }
  n_awgrid = SetTEGrid(awgrid, NULL, n, awmin, awmax);

  return 0;
}

int GetAWGrid(double **a) {
  *a = awgrid;
  return n_awgrid;
}

void SetOptimizeMaxIter(int m) {
  optimize_control.maxiter = m;
}

void SetOptimizeStabilizer(double m) {
  if (m < 0) {
    optimize_control.iset = 0;
  } else {
    optimize_control.iset = 1;
    optimize_control.stabilizer = m;
  }
}

void SetOptimizeTolerance(double c) {
  optimize_control.tolerance = c;
}

void SetOptimizePrint(int m) {
  optimize_control.iprint = m;
}

void SetConfigEnergyMode(int m) {
  optimize_control.mce = m;
}

int ConfigEnergyMode() {
  return optimize_control.mce;
}

void SetOptimizeControl(double tolerance, double stabilizer, 
			int maxiter, int iprint) {
  optimize_control.maxiter = maxiter;
  optimize_control.stabilizer = stabilizer;
  optimize_control.tolerance = tolerance;
  optimize_control.iprint = iprint;  
  optimize_control.iset = 1;
}

void SetScreening(int n_screen, int *screened_n, 
		  double screened_charge, int kl) {
  optimize_control.screened_n = screened_n;
  optimize_control.screened_charge = screened_charge;
  optimize_control.n_screen = n_screen;
  optimize_control.screened_kl = kl;
}

int SetRadialGrid(int maxrp, double ratio, double asymp,
		  double rmin, double qr) {
  if (maxrp > MAXRP) {
    printf("MAXRP must be <= %d\n", MAXRP);
    printf("to enlarge the limit, change MAXRP in global.h\n");
    return -1;
  }
  if (maxrp < 0) maxrp = DMAXRP;
  potential->maxrp = maxrp;
  if (asymp < 0 && ratio < 0) {
    asymp = GRIDASYMP;
    ratio = GRIDRATIO;
  }
  if (rmin <= 0) rmin = GRIDRMIN;
  potential->rmin = rmin;
  if (ratio == 0) potential->ratio = GRIDRATIO;
  else potential->ratio = ratio;
  if (asymp == 0) potential->asymp = GRIDASYMP;
  else potential->asymp = asymp;
  if (qr <= 0) potential->qr = GRIDQR;
  else potential->qr = qr;
  potential->flag = 0;
  return 0;  
}

void AdjustScreeningParams(double *u) {
  int i;
  double c;
  
  c = 0.5*u[potential->maxrp-1];
  for (i = 0; i < potential->maxrp; i++) {
    if (u[i] > c) break;
  }
  potential->lambda = log(2.0)/potential->rad[i];
}

int PotentialHX(AVERAGE_CONFIG *acfg, double *u) {
  int md, md1, jmax, j, i, m, jm;
  double n1, a, b, d, d0, d1;
  double *u0, *ue, *ue1;
  
  if (potential->N < 1+EPS3) return -1; 
  if (acfg->n_shells <= 0) return 0;
  md = potential->mode % 10;
  md1 = potential->mode / 10;
  ue1 = _dwork14;
  ue = _dwork15;
  u0 = _dwork16;
  for (i = 0; i < potential->maxrp; i++) {
    ue1[i] = 0;
    ue[i] = 0;
    u0[i] = 0;
    u[i] = 0;
  }
  if (md == 0) {
    jmax = PotentialHX1(acfg, -1);
    for (i = 0; i < potential->maxrp; i++) {
      ue1[i] = _phase[i];
      ue[i] = _dphase[i];
      u0[i] = _dphasep[i];
    }
  } else if (md == 1) {
    jmax = 0;
    for (i = 0; i < acfg->n_shells; i++) {
      j = PotentialHX1(acfg, i);
      if (j > jmax) jmax = j;
      for (m = 0; m < potential->maxrp; m++) {
	ue1[m] += _phase[m];
	ue[m] += _dphase[m];
	u0[m] += _dphasep[m];
      }
    }
    for (m = 0; m < potential->maxrp; m++) {
      ue1[m] /= acfg->n_shells;
      ue[m] /= acfg->n_shells;
      u0[m] /= acfg->n_shells;
    }
  }
  if (jmax <= 0) return jmax;
  for (m = 0; m < potential->maxrp; m++) {
    u[m] = u0[m] - ue[m];
  }
  jm = 0;
  for (j = 0; j < potential->maxrp; j++) {
    if (u[j] > potential->Z[j]) {
      u[j] = potential->Z[j];
      jm = j;
    }
  }

  n1 = potential->N-1;
  for (jm = jmax; jm >= 10; jm--) {
    if (n1 > u[jm] && u[jm] > u[jm-1]) {
      break;
    }
  }
  d0 = log(potential->rad[jm-1]);
  d1 = log(potential->rad[jm]);
  a = log(n1 - u[jm-1]);
  b = log(n1 - u[jm]);
  d = (b-a)/(d1-d0);    
  for (j = jm+1; j < potential->maxrp; j++) {
    u[j] = d*(log(potential->rad[j]/potential->rad[jm])) + b;
    u[j] = n1 - exp(u[j]);
  }

  for (m = jmax; m > 50; m--) {
    if (fabs(u[m]-n1) > EPS6) break;
  }
  potential->r_core = m+1;
  return jmax;
}

int PotentialHX1(AVERAGE_CONFIG *acfg, int ik) {
  int i, j, k, kk, kk0, kk1, k1, k2, j1, j2;
  int ic, jmax, jmaxk, m, jm;
  ORBITAL *orb1, *orb2;
  double large, small, a, b, c, d, d0, d1, c0, c1, fk, gk;
  CONFIG_GROUP *gc;
  CONFIG *cfg;
  SHELL *s1, *s2;
  double *ue1, *ue, *u, *w;
  
  if (potential->N < 1+EPS3) return -1; 
  ue1 = _phase;
  ue = _dphase;
  u = _dphasep;
  w = potential->W;
  for (m = 0; m < potential->maxrp; m++) {
    w[m] = 0.0;
    ue1[m] = 0.0;
    ue[m] = 0.0;
    u[m] = 0.0;
  }
  if (acfg->n_shells <= 0) return 0;
  if (ik >= acfg->n_shells) ik = -1;
  jmax = -1;
  orb2 = NULL;
  for (i = 0; i < acfg->n_shells; i++) {
    k1 = OrbitalExists(acfg->n[i], acfg->kappa[i], 0.0);
    if (k1 < 0) continue;
    orb1 = GetOrbital(k1);
    if (orb1->wfun == NULL) continue;
    if (ik == i) orb2 = orb1;
    for (m = 0; m <= orb1->ilast; m++) {
      large = Large(orb1)[m];
      small = Small(orb1)[m];
      w[m] += acfg->nq[i]*(large*large + small*small);
    }
    if (jmax < orb1->ilast) jmax = orb1->ilast;    
  }
  if (jmax < 0) return jmax;    
  for (i = 0; i < acfg->n_shells; i++) {
    k1 = OrbitalExists(acfg->n[i], acfg->kappa[i], 0.0);
    if (k1 < 0) continue;
    orb1 = GetOrbital(k1);
    GetYk(0, _yk, orb1, orb1, k1, k1, -1);
    for (m = 0; m <= jmax; m++) {
      if (i == ik) {
	u[m] += (acfg->nq[i]-1.0)*_yk[m];
      } else {
	u[m] += acfg->nq[i]*_yk[m];
      }
    }
  }    

  potential->rhx = 0;
  potential->dhx = 0;
  potential->ahx = 0;
  potential->chx = 0;
  double n1 = potential->N-1;
  if (potential->hxs) {
    potential->rhx = 0;
    j = -1;
    c = potential->Z[potential->maxrp-1];
    c = potential->hx0 + potential->hx1*(c-potential->N)/c;
    potential->ahx = potential->hxs*c;
    double nr = 0;
    if (orb2 != NULL) {
      nr = (double)(orb2->n - GetLFromKappa(orb2->kappa)/2);
    }
    for (m = 0; m <= jmax; m++) {
      a = w[m]*potential->rad[m];
      if (orb2 != NULL) {
	large = Large(orb2)[m];
	small = Small(orb2)[m];
	d = potential->rad[m];
	b = (large*large+small*small)*d;
	if (acfg->nq[ik] < 2) {
	  a += (2-acfg->nq[ik])*b;
	}
	a = pow(a, ONETHIRD)-pow(2*b,ONETHIRD);
      } else {
	a = pow(a, ONETHIRD);
      }
      ue1[m] = potential->ahx*a;
      if (j < 0 && m > 0 && ue1[m-1] > 0 && u[m] - ue1[m] > n1) {
	potential->rhx = potential->rad[m];
	j = m;
      }
    }
    a = fabs(potential->ihx);
    if (ik < 0 && a > 0.01) {
      potential->dhx = potential->rad[j];
      for (m = j; m > 0; m--) {
	if (ue1[m-1] < ue1[m]) break;
      }
      potential->rhx = potential->rad[m];
      b = potential->dhx - potential->rhx;
      c = 0.25*potential->dhx;
      if (b < c) {
	potential->rhx = potential->dhx - c;
	potential->dhx = c;
      } else {
	potential->dhx = b;
      }      
      c = 1.0 - n1/potential->N;
      for (m = 0; m <= jmax; m++) {
	d = (potential->rad[m]-potential->rhx)*a/potential->dhx;
	d = 1.0/(1.0 + exp(-d));
	d0 = u[m]*c;
	d1 = ue1[m];
	ue[m] = d0*d + d1*(1-d);
	if (d > 0.5 && u[m]-ue[m] > n1) {
	  ue[m] = u[m]-n1;
	}
      }
    } else {
      for (m = 0; m <= jmax; m++) {
	ue[m] = ue1[m];
      }
    }
  }

  for (m = jmax+1; m < potential->maxrp; m++) {
    u[m] = u[jmax];
    ue1[m] = 0.0;
    if (ue[jmax] > ue1[jmax]) ue[m] = 1.0;
    else ue[m] = 0.0;
  }
  return jmax;
}

double SetPotential(AVERAGE_CONFIG *acfg, int iter) {
  int jmax, i, j, k;
  double *u, *v, a, b, c, r, rn;

  u = potential->U;
  v = _dwork2;

  jmax = PotentialHX(acfg, u);
  rn = GetAtomicR()*2;
  if (jmax > 0) {
    if (iter < 3) {
      r = 1.0;
      for (j = 0; j < potential->maxrp; j++) {
	v[j] = u[j];
      }
    } else {	
      r = 0.0;
      k = 0;
      a = optimize_control.stabilizer;
      b = fabs(potential->hxs);
      if (b > 0.75) {
	a *= pow(0.75/b,2);
      }
      b = 1.0 - a;
      for (j = 0; j < potential->maxrp; j++) {
	if (u[j] + 1.0 != 1.0 && potential->rad[j] > rn) {
	  r += fabs(1.0 - v[j]/u[j]);
	  k++;
	}
	u[j] = b*v[j] + a*u[j];
	v[j] = u[j];
      }
      r /= k;
    }
    AdjustScreeningParams(u);
    SetPotentialVc(potential);
    for (j = 0; j < potential->maxrp; j++) {
      a = u[j] - potential->Z[j];
      b = potential->Vc[j]*potential->rad[j];
      u[j] = a - b;
      u[j] /= potential->rad[j];
    }
    SetPotentialU(potential, 0, NULL);
    SetPotentialVT(potential);
  } else {
    if (potential->N < 1.0+EPS3) {
      SetPotentialVc(potential);
      SetPotentialU(potential, -1, NULL);
      SetPotentialVT(potential);
      return 0.0;
    }
    r = potential->Z[potential->maxrp-1];
    b = (1.0 - 1.0/potential->N);
    for (i = 0; i < acfg->n_shells; i++) {
      a = acfg->nq[i];
      c = acfg->n[i];
      c = r/(c*c);
      for (j = 0; j < potential->maxrp; j++) {
	if (potential->rad[j] < 0.1*potential->atom->rn) {
	  u[j] = 0.0;
	} else {
	  u[j] += a*b*(1.0 - exp(-c*potential->rad[j]));
	}
      }
    }
    AdjustScreeningParams(u);
    SetPotentialVc(potential);
    for (j = 0; j < potential->maxrp; j++) {
      a = u[j] - potential->Z[j];
      b = potential->Vc[j]*potential->rad[j];
      u[j] = a - b;
      u[j] /= potential->rad[j];
    }
    SetPotentialU(potential, 0, NULL);
    SetPotentialVT(potential);
    r = 1.0;
  }
  
  return r;
}

int GetPotential(char *s) {
  AVERAGE_CONFIG *acfg;
  ORBITAL *orb1, *orb2;
  double large1, small1, large2, small2;
  int norbs, jmax, kmin, kmax;  
  FILE *f;
  int i, j;
  double rb, rb1, rc;

  /* get the average configuration for the groups */
  acfg = &(average_config);

  f = fopen(s, "w");
  if (!f) return -1;
  
  fprintf(f, "# Lambda = %12.5E\n", potential->lambda);
  fprintf(f, "#      A = %12.5E\n", potential->a);
  fprintf(f, "#     ar = %12.5E\n", potential->ar);
  fprintf(f, "#     br = %12.5E\n", potential->br);
  fprintf(f, "#     qr = %12.5E\n", potential->qr);
  rc = potential->r_core > 0?potential->rad[potential->r_core]:0;
  fprintf(f, "#     rc = %12.5E\n", rc);
  rb = potential->ib>0?potential->rad[potential->ib]:0;
  rb1 = potential->ib1>0?potential->rad[potential->ib1]:0;
  fprintf(f, "#     rb = %12.5E\n", rb);
  fprintf(f, "#    rb1 = %12.5E\n", rb1);
  fprintf(f, "#    bqp = %12.5E\n", potential->bqp);
  fprintf(f, "#     nb = %d\n", potential->nb);
  fprintf(f, "#   mode = %d\n", potential->mode);
  fprintf(f, "#    HXS = %12.5E\n", potential->hxs);
  fprintf(f, "#    AHX = %12.5E\n", potential->ahx);
  fprintf(f, "#    IHX = %12.5E\n", potential->ihx);
  fprintf(f, "#    RHX = %12.5E\n", potential->rhx);
  fprintf(f, "#    DHX = %12.5E\n", potential->dhx);
  fprintf(f, "#    CHX = %12.5E\n", potential->chx);
  fprintf(f, "#    HX0 = %12.5E\n", potential->hx0);
  fprintf(f, "#    HX1 = %12.5E\n", potential->hx1);
  fprintf(f, "#   nmax = %d\n", potential->nmax);
  fprintf(f, "#  maxrp = %d\n", potential->maxrp);
  fprintf(f, "# Mean configuration: %d\n", acfg->n_shells);
  for (i = 0; i < acfg->n_shells; i++) {
    fprintf(f, "# %2d %2d\t%10.3E\n", acfg->n[i], acfg->kappa[i], acfg->nq[i]);
  }
  fprintf(f, "\n\n");
  for (i = 0; i < potential->maxrp; i++) {
    fprintf(f, "%5d %14.8E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E\n",
	    i, potential->rad[i], potential->Z[i],
	    potential->Z[i]-GetAtomicEffectiveZ(potential->rad[i]),
	    potential->Vc[i]*potential->rad[i],
	    potential->U[i]*potential->rad[i],
	    _dwork14[i], _dwork15[i], _dwork16[i],
	    potential->ZSE[0][i],
	    potential->ZSE[1][i],
	    potential->ZSE[2][i],
	    potential->ZSE[3][i],
	    potential->ZSE[4][i],
	    potential->ZVP[i]);
  }    
  fclose(f);    
  return 0;
}

double GetResidualZ(void) {
  double z;
  z = potential->Z[potential->maxrp-1];
  if (potential->N > 0) z -= potential->N - 1;
  return z;
}

double GetRMax(void) {
  return potential->rad[potential->maxrp-10];
}

int SetAverageConfig(int nshells, int *n, int *kappa, double *nq) {
  int i;
  if (nshells <= 0) return -1;
  if (average_config.n_shells > 0) {
    average_config.kappa = (int *) realloc(average_config.kappa, 
					   sizeof(int)*nshells);
    average_config.nq = (double *) realloc(average_config.nq, 
					   sizeof(double)*nshells);
    average_config.n = (int *) realloc(average_config.n, 
				       sizeof(int)*nshells);
  } else {
    average_config.kappa = (int *) malloc(sizeof(int)*nshells);
    average_config.nq = (double *) malloc(sizeof(double)*nshells);
    average_config.n = (int *) malloc(sizeof(int)*nshells);
  }
  for (i = 0; i < nshells; i++) {
    average_config.n[i] = n[i];
    average_config.kappa[i] = kappa[i];
    average_config.nq[i] = nq[i];
  }
  average_config.n_shells = nshells;
  average_config.n_cfgs = 1;
  return 0;
}

void FreezeOrbital(char *s, int m) {
  char s0[1024], *p;
  int i, j, k, n, nc;
  CONFIG *cfg;
  ORBITAL *orb;

  if (m >= 0) _orthogonalize_mode = m;
  strncpy(s0, s, 1023);
  s = s0;
  while (*s == ' ') s++;
  if (*s == '\0') return;
  n = StrSplit(s, ' ');
  p = s;
  for (i = 0; i < n; i++) {
    while (*p == ' ') p++;
    nc = GetConfigFromString(&cfg, p);
    for (j = 0; j < nc; j++) {
      if (cfg[j].n_shells != 1) {
	printf("incorrect freeze orbital format: %d %s\n", i, p);
	continue;
      }
      k = OrbitalExists((cfg[j].shells)[0].n, (cfg[j].shells)[0].kappa, 0);
      if (k < 0) continue;
      orb = GetOrbital(k);
      if (orb->isol == 1) {
	orb->isol = 2;
	potential->nfrozen++;
      }
    }
    if (nc > 0) free(cfg);
    while (*p) p++;
    p++;
  }
}

int OptimizeLoop(AVERAGE_CONFIG *acfg) {
  double tol, atol, tol0, atol0, tol1, a, b, ahx, hxs0;
  ORBITAL orb_old, *orb;
  int i, k, iter, no_old;
  
  no_old = 0;
  iter = 0;
  tol = 1.0;
  atol = 1e1;
  tol0 = optimize_control.tolerance*ENERELERR;
  tol1 = optimize_control.tolerance*ENERELERR1;
  atol0 = optimize_control.tolerance*ENEABSERR;
  ahx = 1.0;
  hxs0 = potential->hxs;
  if (fabs(hxs0)<1e-5) ahx = 0.0;
  while (((tol > tol0 || atol > atol0) && (tol > tol1)) || ahx > 1e-5) {
    if (iter > optimize_control.maxiter) break;
    if (fabs(hxs0) > 1e-5) {
      ahx = iter>5?0.0:exp(-(1+iter)*0.75);
      if (ahx < 0.01) ahx = 0.0;
      potential->hxs = hxs0*(1-ahx);
    }
    a = SetPotential(acfg, iter);
    FreeYkArray();
    tol = 0.0;
    atol = 0.0;
    for (i = 0; i < acfg->n_shells; i++) {
      k = OrbitalExists(acfg->n[i], acfg->kappa[i], 0.0);
      if (k < 0) {
	orb_old.energy = 0.0;
	orb = GetNewOrbital(acfg->n[i], acfg->kappa[i], 1.0);
	orb->energy = 1.0;
	no_old = 1;	
      } else {
	orb = GetOrbital(k);
	if (orb->isol == 0 || orb->wfun == NULL) {
	  orb_old.energy = 0.0;
	  orb->energy = 1.0;
	  orb->kappa = acfg->kappa[i];
	  orb->n = acfg->n[i];
	  no_old = 1;	
	} else if (orb->isol == 1) {
	  orb_old.energy = orb->energy; 
	  free(orb->wfun);
	  orb->wfun = NULL;
	  orb->isol = 0;
	  no_old = 0;
	} else {
	  continue;
	}
      }

      if (SolveDirac(orb) < 0) {
	return -1;
      }
      
      if (no_old) { 
	tol = 1.0;
	atol = 1e1;
	continue;
      } 
      b = fabs(1.0 - orb_old.energy/orb->energy);
      if (tol < b) tol = b;
      b = fabs(orb_old.energy - orb->energy);
      if (atol < b) atol = b;
    }
    if (tol < a) tol = a;
    if (optimize_control.iprint) {
      printf("optimize loop: %4d %13.5E %13.5E %13.5E %13.5E %13.5E %13.5E %13.5E\n", iter, tol, atol, a, tol0, tol1, atol0, ahx);
    }
    iter++;
  }

  return iter;
}

void CopyPotentialOMP(int init) {
  SetReferencePotential(hpotential, potential, 1);
  SetReferencePotential(rpotential, potential, 0);
#if USE_MPI == 2
  if (!MPIReady()) {
    InitializeMPI(0, 0);
    return;
  }
  POTENTIAL pot;
  memcpy(&pot, potential, sizeof(POTENTIAL));
#pragma omp parallel shared(pot)
  {
    if (init && MyRankMPI() != 0) {
      potential = malloc(sizeof(POTENTIAL));
    }
    memcpy(potential, &pot, sizeof(POTENTIAL));    
  }
  memcpy(&pot, hpotential, sizeof(POTENTIAL));
#pragma omp parallel shared(pot)
  {
    if (init && MyRankMPI() != 0) {
      hpotential = malloc(sizeof(POTENTIAL));
    }
    memcpy(hpotential, &pot, sizeof(POTENTIAL));    
  }
  memcpy(&pot, rpotential, sizeof(POTENTIAL));
#pragma omp parallel shared(pot)
  {
    if (init && MyRankMPI() != 0) {
      rpotential = malloc(sizeof(POTENTIAL));
    }
    memcpy(rpotential, &pot, sizeof(POTENTIAL));    
  }
#endif
}

#define NXS 7
int OptimizeRadial(int ng, int *kg, int ic, double *weight, int ife) {
  AVERAGE_CONFIG *acfg;
  double a, b, c, z, emin, smin, hxs[NXS], ehx[NXS], mse;
  int iter, i, j;

  mse = qed.se;
  qed.se = -1000000;
  /* get the average configuration for the groups */
  acfg = &(average_config);
  if (ng > 0) {
    /*
    if (ng > 1) {
      printf("\nWarning: more than 1 configuration groups");
      printf(" are used in OptimizeRadial.\n");
      printf("It is usually best to use the lowest lying configuration group.\n");
      printf("Make sure that you know what you are doing.\n\n");
    }
    */
    if (acfg->n_shells > 0) {      
      acfg->n_cfgs = 0;
      acfg->n_shells = 0;
      free(acfg->n);
      free(acfg->kappa);
      free(acfg->nq);
      if (acfg->ng > 0) {
	free(acfg->kg);
	free(acfg->weight);
	acfg->ng = 0;
      }
      acfg->n = NULL;
      acfg->nq = NULL;
      acfg->kappa = NULL;
      acfg->kg = NULL;
      acfg->weight = NULL;
    }
    GetAverageConfig(ng, kg, ic, weight,
		     optimize_control.n_screen,
		     optimize_control.screened_n,
		     optimize_control.screened_charge,
		     optimize_control.screened_kl, acfg); 
  } else {
    if (acfg->n_shells <= 0) {
      printf("No average configuation exist. \n");
      printf("Specify with AvgConfig, ");
      printf("or give config groups to OptimizeRadial.\n");
      qed.se = mse;
      return -1;
    }
  }
  
  a = 0.0;
  for (i = 0; i < acfg->n_shells; i++) {
    if (optimize_control.iprint) {
      MPrintf(-1, "avgcfg: %d %d %d %g\n",
	      i, acfg->n[i], acfg->kappa[i], acfg->nq[i]);
    }
    a += acfg->nq[i];
  }
  if (a > potential->atom->atomic_number) {
    for (i = 0; i < acfg->n_shells; i++) {
      acfg->nq[i] *= potential->atom->atomic_number/a;
      if (optimize_control.iprint) {
	MPrintf(-1, "avgcfg: %d %d %d %g\n",
		i, acfg->n[i], acfg->kappa[i], acfg->nq[i]);
      }
    }
    a = potential->atom->atomic_number;
  }
  potential->N = a;  

  /* setup the radial grid if not yet */
  if (potential->flag == 0) {
    SetOrbitalRGrid(potential);
  }
  
  int nmax = potential->nmax-1;
  if (potential->nb > 0 && nmax < potential->nb) nmax = potential->nb;
  for (i = 0; i < acfg->n_shells; i++) {
    if (acfg->n[i] > nmax) {
      printf("too large n in avgcfg: %d %d %d %d %d %g\n",
	     ife, nmax, i, acfg->n[i], acfg->kappa[i], acfg->nq[i]);
      if (ife) {
	return -1;
      }
      j = GetLFromKappa(acfg->kappa[i]);
      acfg->n[i] = potential->nmax;
      if (j/2 >= acfg->n[i]) {
	acfg->kappa[i] = -1;
      }
    }
  }
  SetPotentialZ(potential);
  SetReferencePotential(hpotential, potential, 1);
  SetReferencePotential(rpotential, potential, 0);
  z = potential->Z[potential->maxrp-1];
  if (a > 0.0) z = z - a + 1;
  potential->a = 0.0;
  potential->lambda = 0.5*z;
  if (potential->N > 1) {
    potential->r_core = potential->maxrp-5;
  } else {
    potential->r_core = 50;
  }

  if (optimize_control.iset == 0) {
    optimize_control.stabilizer = 0.25 + 0.5*(z/potential->Z[potential->maxrp-1]);
  }

  if (potential->mode/10 == 2) {
    a = 0.7/(NXS-1.0);
    hxs[0] = 0.5;
    for (i = 1; i < NXS; i++) {
      hxs[i] = hxs[i-1] + a;
    }      
    for (i = 0; i < NXS; i++) {
      ReinitRadial(1);
      potential->hxs = hxs[i];
      iter = OptimizeLoop(acfg);
      if (iter > optimize_control.maxiter) {
	printf("Maximum iteration reached in OptimizeRadial %d %d\n", i, iter);
	return -1;
      }
      if (ng > 0) {
	ehx[i] = 0.0;
	for (j = 0; j < ng; j++) {
	  ehx[i] += TotalEnergyGroup(kg[j])*acfg->weight[j];
	}
      } else {
	ehx[i] = AverageEnergyAvgConfig(acfg);
      }
      if (optimize_control.iprint) {
	printf("hxs iter: %d %d %g %g %18.10E\n", ng, i, hxs[i], potential->ahx, ehx[i]);
      }
    }
    
    b = 0.001;
    a = hxs[0];
    emin = 1e10;
    smin = 0.0;
    while (a <= hxs[NXS-1]) {
      UVIP3P(2, NXS, hxs, ehx, 1, &a, &c);
      if (c < emin) {
	emin = c;
	smin = a;
      }
      a += b;
    }
    
    ReinitRadial(1);
    potential->hxs = smin;
    
    if (optimize_control.iprint) {
      printf("hxs: %g %18.10E\n", smin, emin);
    }
    iter = OptimizeLoop(acfg);
  } else {
    iter = OptimizeLoop(acfg);
  }

  qed.se = mse;
  CopyPotentialOMP(0);
  return iter;
}      
#undef NXS

static double **_refine_wfb = NULL;
static int _refine_msglvl = 0;
static double EnergyFunc(int *n, double *x) {
  double a;
  int i, k;
  ORBITAL *orb;
  
  for (i = 0; i < potential->maxrp; i++) {
    potential->U[i] = _dphasep[i];
    for (k = 0; k < *n; k++) {      
      potential->U[i] += x[k]*_refine_wfb[k][i];
    }
  }  
  //SetPotentialVc(potential);
  SetPotentialU(potential, 0, NULL);
  SetPotentialVT(potential);
  ReinitRadial(1);
  ClearOrbitalTable(0);
  if (average_config.ng > 0) {
    a = 0.0;
    for (k = 0; k < average_config.ng; k++) {
      a += TotalEnergyGroup(average_config.kg[k])*average_config.weight[k];
    }
  } else {
    a = AverageEnergyAvgConfig(&average_config);
  }
  if (_refine_msglvl > 1) {
    printf("refine step: ");
    for (i = 0; i < *n; i++) {
      printf("%g ", x[i]);
    }
    printf("--- %18.10E\n", a);
  }
  return a;
}

int RefineRadial(int maxfun, int msglvl) {
  int n, i, k, ierr, mode, nfe, se, *lw;
  double xtol;
  double f0, f, *x, *scale, *dw;
  AVERAGE_CONFIG *a;
  ORBITAL *orb;

  se = qed.se;
  qed.se = -1000000;
  if (maxfun <= 0) maxfun = 10000;
  _refine_msglvl = msglvl;
  xtol = EPS3;
  n = 0;
  k = 0;
  a = &average_config;  
  lw = malloc(sizeof(int)*a->n_shells*2);
  for (i = 0; i < a->n_shells; i++) {
    if (a->n[i] != k) {
      k = a->n[i];
      lw[n] = k;
      n++;
    }
  }
  _refine_wfb = malloc(sizeof(double *)*n);
  for (i = 0; i < n; i++) {
    _refine_wfb[i] = malloc(sizeof(double)*potential->maxrp);
    k = OrbitalIndex(lw[i], -1, 0);
    orb = GetOrbitalSolved(k);
    for (k = 0; k < potential->maxrp; k++) {
      _refine_wfb[i][k] = orb->wfun[k];
    }
  }
  mode = 0;
  x = malloc(sizeof(double)*n);
  scale = malloc(sizeof(double)*n);
  dw = malloc(sizeof(double)*(n*2+n*(n+4)+1));
  for (i = 0; i < n; i++) {
    x[i] = 0.0;
    scale[i] = 1e-5;
  }

  for (i = 0; i < potential->maxrp; i++) {
    _dphasep[i] = potential->U[i];
  }
  nfe = 0;
  ierr = 0;
  SUBPLX(EnergyFunc, n, xtol, maxfun, mode, scale, x,
	 &f, &nfe, dw, lw, &ierr);
  if (msglvl > 0) {
    printf("refine final: ");
    for (i = 0; i < n; i++) {
      printf("%g ", x[i]);
    }
    printf("--- %18.10E %d %d\n", f, ierr, nfe);
  }
  free(x);
  free(scale);
  free(lw);
  free(dw);
  for (i = 0; i < n; i++) {
    free(_refine_wfb[i]);
  }
  free(_refine_wfb);
  qed.se = se;
  return 0;
}

void Orthogonalize(ORBITAL *orb) {
  int i, k;

  for (i = 0; i < potential->maxrp; i++) {
    _yk[i] = 1;
  }
  k = ((abs(orb->kappa)-1)*2)+(orb->kappa>0);
  ORBMAP *om = &_orbmap[k];  
  ORBITAL *orb0;
  double a, b, *p, *q, *p0, *q0;
  double qn, qn0;
  Integrate(_yk, orb, orb, 1, &qn0, 0);
  p = Large(orb);
  q = Small(orb);
  b = 1.0;
  for (k = 0; k < _norbmap0; k++) {
    orb0 = om->opn[k];
    if (orb0 == NULL || orb->n == orb0->n) continue;
    if (orb0->isol <= 0) continue;
    if (_orthogonalize_mode > 1 && orb0->isol < 2) continue;
    if (orb->isol == 3 && orb0->n > orb->n) continue;
    Integrate(_yk, orb, orb0, 1, &a, 0);
    p0 = Large(orb0);
    q0 = Small(orb0);
    for (i = 0; i <= orb->ilast; i++) {
      p[i] -= p0[i]*a;
      q[i] -= q0[i]*a;
    }
  }
  Integrate(_yk, orb, orb, 1, &qn, 0);
  b = 1.0/sqrt(qn);
  for (i = 0; i <= orb->ilast; i++) {
    p[i] *= b;
    q[i] *= b;
  }
  orb->energy = ZerothHamilton(orb, orb);
}

int SolveDirac(ORBITAL *orb) {
  int err;
#ifdef PERFORM_STATISTICS
  clock_t start, stop;
  start = clock();
#endif

  err = 0;
  if (potential->npseudo > 0 &&
      orb->n < 1000000 &&
      orb->n > potential->npseudo) {
    int k = GetLFromKappa(orb->kappa)/2;
    if (orb->n > k+potential->dpseudo) {
      return 0;
    }
  }
  potential->flag = -1;
  err = RadialSolver(orb, potential);
  if (err) { 
    printf("Error ocuured in RadialSolver, %d\n", err);
    printf("%d %d %10.3E\n", orb->n, orb->kappa, orb->energy);
    exit(1);
  }
  if (_orthogonalize_mode > 0) {
    if (potential->nfrozen > 0 && orb->n > 0) {
      Orthogonalize(orb);
    }
  }
#ifdef PERFORM_STATISTICS
  stop = clock();
  rad_timing.dirac += stop - start;
#endif

  return err;
}

int WaveFuncTable(char *s, int n, int kappa, double e) {
  int i, k;
  FILE *f;
  ORBITAL *orb;
  double z, a, ke, y;

  e /= HARTREE_EV;
  k = OrbitalIndex(n, kappa, e);
  if (k < 0) return -1;
  f = fopen(s, "w");
  if (!f) return -1;
  
  orb = GetOrbitalSolved(k);    
  fprintf(f, "#      n = %2d\n", n);
  fprintf(f, "#  kappa = %2d\n", kappa);
  fprintf(f, "# energy = %15.8E\n", orb->energy*HARTREE_EV);
  fprintf(f, "#SelfEne = %15.8E\n", orb->se*HARTREE_EV);
  fprintf(f, "#     vc = %15.8E\n", MeanPotential(k, k)*HARTREE_EV);
  fprintf(f, "#  ilast = %4d\n", orb->ilast);
  //fprintf(f, "#     im = %4d\n", orb->im);
  fprintf(f, "#    rfn = %15.8E\n", orb->rfn);
  fprintf(f, "#    pdx = %15.8E\n", orb->pdx);
  fprintf(f, "#   bqp0 = %15.8E\n", orb->bqp0);
  fprintf(f, "#   bqp1 = %15.8E\n", orb->bqp1);
  fprintf(f, "#    idx = %d\n", k);
  ORBITAL *horb = orb->horb;
  ORBITAL *rorb = orb->rorb;
  if (n != 0) {
    fprintf(f, "\n\n");
    if (n < 0) k = potential->ib;
    else k = 0;
    for (i = k; i <= orb->ilast; i++) {
      fprintf(f, "%-4d %14.8E %13.6E %13.6E %13.6E %13.6E", 
	      i, potential->rad[i], 
	      (potential->Vc[i])*potential->rad[i],
	      potential->U[i] * potential->rad[i],
	      Large(orb)[i], Small(orb)[i]); 
      if (rorb != NULL && rorb->wfun != NULL) {
	fprintf(f, " %13.6E %13.6E", Large(rorb)[i], Small(rorb)[i]);
      } else {
	fprintf(f, " %13.6E %13.6E", 0.0, 0.0);
      }
      if (horb != NULL && horb->wfun != NULL) {
	fprintf(f, " %13.6E %13.6E\n", Large(horb)[i], Small(horb)[i]);
      } else {
	fprintf(f, " %13.6E %13.6E\n", 0.0, 0.0);
      }
    }
  } else {
    a = GetPhaseShift(k);
    while(a < 0) a += TWO_PI;
    a -= (int)(a/TWO_PI);
    fprintf(f, "#  phase = %15.8E\n", a);
    fprintf(f, "\n\n");
    z = GetResidualZ();
    e = orb->energy;
    a = FINE_STRUCTURE_CONST2 * e;
    ke = sqrt(2.0*e*(1.0+0.5*a));
    y = (1.0+a)*z/ke; 
    for (i = 0; i <= orb->ilast; i++) {
      fprintf(f, "%-4d %14.8E %13.6E %13.6E %13.6E %13.6E", 
	      i, potential->rad[i],
	      (potential->Vc[i])*potential->rad[i],
	      potential->U[i] * potential->rad[i],
	      Large(orb)[i], Small(orb)[i]);
      if (horb != NULL && horb->wfun != NULL) {
	fprintf(f, " %13.6E %13.6E\n", Large(horb)[i], Small(horb)[i]);
      } else {
	fprintf(f, " %13.6E %13.6E\n", 0.0, 0.0);
      }
    }
    for (; i < potential->maxrp; i += 2) {
      a = ke * potential->rad[i];
      a = a + y*log(2.0*a);
      a = Large(orb)[i+1] - a;
      a = a - ((int)(a/(TWO_PI)))*TWO_PI;
      if (a < 0) a += TWO_PI;
      fprintf(f, "%-4d %14.8E %13.6E %13.6E %13.6E %13.6E %13.6E %13.6E\n",
	      i, potential->rad[i],
	      Large(orb)[i], Large(orb)[i+1], 
	      Small(orb)[i], a, 0.0, 0.0);
    }
  }

  fclose(f);

  return 0;
}

double PhaseRDependent(double x, double eta, double b) {
  double tau, tau2, y, y2, t, a1, a2, sb;
  
  y = 1.0/x;
  y2 = y*y;
  tau2 = 1.0 + 2.0*eta*y - b*y2;
  tau = x*sqrt(tau2);
  
  t = eta*log(x+tau+eta) + tau - eta;
  if (b > 0.0) {
    sb = sqrt(b);
    a1 = b - eta*x;
    a2 = tau*sb;
    tau2 = atan2(a1, a2);
    tau2 -= atan2(-eta, sb);
    t += sb*tau2;
  } else if (b < 0.0) {
    b = -b;
    sb = sqrt(b);
    a1 = 2.0*(b+eta*x)/(sb*b*x);
    a2 = 2.0*tau/(b*x);
    tau2 = log(a1+a2);
    tau2 -= log(2.0*(eta/sb + 1.0)/b);
    t -= sb*tau2;
  }
  
  return t;
}

double GetPhaseShift(int k) {
  ORBITAL *orb;
  double phase1, r, y, z, ke, e, a, b1;
  int i;

  orb = GetOrbitalSolved(k);
  if (orb->n > 0) return 0.0;

  if (orb->phase) return *(orb->phase);

  z = GetResidualZ();
  e = orb->energy;
  a = FINE_STRUCTURE_CONST2 * e;
  ke = sqrt(2.0*e*(1.0 + 0.5*a));
  y = (1.0 + a)*z/ke;

  i = potential->maxrp - 1;  
  phase1 = orb->wfun[i];
  r = potential->rad[i-1];  
  b1 = orb->kappa;
  b1 = b1*(b1+1.0) - FINE_STRUCTURE_CONST2*z*z;
 
  a = ke * r;
  b1 = PhaseRDependent(a, y, b1);
  phase1 = phase1 - b1;
  
  orb->phase = malloc(sizeof(double));
  *(orb->phase) = phase1;

  return phase1;  
}

int GetNumBounds(void) {
  return n_orbitals - n_continua;
}

int GetNumOrbitals(void) {
  return n_orbitals;
}

int GetNumContinua(void) {
  return n_continua;
}

int OrbitalIndex(int n, int kappa, double energy) {
  ORBITAL *orb = NULL;
  int k = ((abs(kappa)-1)*2)+(kappa>0);
  if (k >= _korbmap) {
    printf("too large kappa, enlarge korbmap: %d >= %d\n",
	   k, _korbmap);
    Abort(1);
  }
  ORBMAP *om = &_orbmap[k];
  if (n > 0) {
    k = n-1;
    orb = om->opn[k];
  } else if (n < 0) {
    k = -n-1;
    orb = om->onn[k];
  } else {
    for (k = 0; k < om->nzn; k++) {
      if (fabs(energy-om->ozn[k]->energy) < EPS10) {
	orb = om->ozn[k];
	break;
      }
    }
  }
  if (orb == NULL) {
    if (orbitals->lock) {
      SetLock(orbitals->lock);
      if (n > 0) {
	orb = om->opn[k];
      } else if (n < 0) {
	orb = om->onn[k];
      } else {
	for (; k < om->nzn; k++) {
	  if (fabs(energy-om->ozn[k]->energy) < EPS10) {
	    orb = om->ozn[k];
	    break;
	  }
	}
      }
    }
    if (orb == NULL) {
      orb = GetNewOrbitalNoLock(n, kappa, energy);
      k = SolveDirac(orb);
      if (k < 0) {
	MPrintf(-1, "Error occured in solving Dirac eq. err = %d\n", k);
	Abort(1);
      }
    }
    if (orbitals->lock) ReleaseLock(orbitals->lock);
  } else if (orb->isol == 0) {
    if (orbitals->lock) SetLock(orbitals->lock);
    if (orb->isol == 0) {
      k = SolveDirac(orb);
      if (k < 0) {
	MPrintf(-1, "Error occured in solving Dirac eq. err = %d\n", k);
	Abort(1);
      }
    }
    if (orbitals->lock) ReleaseLock(orbitals->lock);
  }
  if (potential->npseudo > 0 && orb->n > potential->npseudo) {
    k = GetLFromKappa(orb->kappa)/2;
    if (orb->n > k+potential->dpseudo) {
      return orb->idx;
    }
  }
  if (!orb->isol) {
    printf("isol0a: %d %d %d %g\n", orb->idx, orb->n, orb->kappa, orb->energy);
    Abort(1);
  }
  return orb->idx;
}

int OrbitalIndexNoLock0(int n, int kappa, double energy) {
  int i, j;
  ORBITAL *orb;
  int resolve_dirac;

  resolve_dirac = 0;
  for (i = 0; i < n_orbitals; i++) {
    orb = GetOrbital(i);
    if (n == 0) {
      if (orb->n == 0 &&
	  orb->kappa == kappa && 
	  orb->energy > 0.0 &&
	  fabs(orb->energy - energy) < EPS10) {
	if (orb->isol == 0) {
	  resolve_dirac = 1;
	  break;
	}
	return i;
      }
    } else if (orb->n == n && orb->kappa == kappa) {
      if (orb->isol == 0) {
	resolve_dirac = 1;
	break;
      }
      return i;
    }
  }    
  if (!resolve_dirac) {
    orb = GetNewOrbitalNoLock(n, kappa, energy);
  } 
  j = SolveDirac(orb);
  if (j < 0) {
    MPrintf(-1, "Error occured in solving Dirac eq. err = %d\n", j);
    Abort(1);
  }
#pragma omp flush
  return i;
}

int OrbitalIndex0(int n, int kappa, double energy) {
  int i;
  if (orbitals->lock) SetLock(orbitals->lock);
  i = OrbitalIndexNoLock0(n, kappa, energy);
  if (orbitals->lock) ReleaseLock(orbitals->lock);
  return i;
}

int OrbitalExistsNoLock(int n, int kappa, double energy) {
  ORBITAL *orb = NULL;
  int k = ((abs(kappa)-1)<<1)+(kappa>0);
  if (k >= _korbmap) {
    printf("too large kappa, enlarge korbmap: %d >= %d\n",
	   k, _korbmap);
    Abort(1);
  }
  ORBMAP *om = &_orbmap[k];
  if (n > 0) {
    k = n-1;
    orb = om->opn[k];
  } else if (n < 0) {
    k = -n-1;
    orb = om->onn[k];
  } else {
    int i;
    for (i = 0; i < om->nzn; i++) {
      if (fabs(energy-om->ozn[i]->energy) < EPS10) {
	orb = om->ozn[i];
	break;
      }
    }
  }
  if (orb != NULL) return orb->idx;  
  return -1;
}

int OrbitalExistsNoLock0(int n, int kappa, double energy) {
  int i;
  ORBITAL *orb;
  for (i = 0; i < n_orbitals; i++) {
    orb = GetOrbital(i);
    if (n == 0) {
      if (orb->kappa == kappa &&
	  fabs(orb->energy - energy) < EPS10) {
	return i;
      }
    } else if (orb->n == n && orb->kappa == kappa) {
      return i;
    }
  }
  return -1;
}

int OrbitalExists(int n, int kappa, double energy) {
  int i;
  i = OrbitalExistsNoLock(n, kappa, energy);
  if (i >= 0) return i;
  if (orbitals->lock) {
    SetLock(orbitals->lock);
    i = OrbitalExistsNoLock(n, kappa, energy);
    ReleaseLock(orbitals->lock);
  }
  return i;
}

void AddOrbMap(ORBITAL *orb) {
  int k = ((abs(orb->kappa)-1)<<1)+(orb->kappa>0);
  if (k >= _korbmap) {
    printf("too large kappa, enlarge korbmap: %d >= %d\n",
	   k, _korbmap);
    Abort(1);
  }
  ORBMAP *om = &_orbmap[k];  
  if (orb->n > 0) {
    k = orb->n-1;
    if (k >= _norbmap0) {
      printf("too many bound orbitals, enlarge norbmap0: %d >= %d\n",
	     k, _norbmap0);
      Abort(1);
    }
    om->opn[k] = orb;
  } else if (orb->n < 0) {
    k = -orb->n-1;
    if (k >= _norbmap1) {
      printf("too many basis orbitals, enlarge norbmap1: %d >= %d\n",
	     k, _norbmap1);
      Abort(1);
    }
    om->onn[k] = orb;
  } else {
    if (om->nzn >= _norbmap2) {
      printf("too many free orbitals, enlarge norbmap2: %d >= %d\n",
	     om->nzn, _norbmap2);
      Abort(1);
    }
    om->ozn[om->nzn] = orb;
    om->nzn++;
  }
}

void RemoveOrbMap(int m) {
  int k, i;
  for (k = 0; k < _korbmap; k++) {
    ORBMAP *om = &_orbmap[k];
    if (m == 0) {
      for (i = 0; i < _norbmap0; i++) {
	om->opn[i] = NULL;
      }
      for (i = 0; i < _norbmap1; i++) {
	om->onn[i] = NULL;
      }
    }
    for (i = 0; i < om->nzn; i++) {
      om->ozn[i] = NULL;
    }
    om->nzn = 0;
  }
}

/*
int AddOrbital(ORBITAL *orb) {

  if (orb == NULL) return -1;
  if (orbitals->lock) SetLock(orbitals->lock);
  orb = (ORBITAL *) ArrayAppend(orbitals, orb, InitOrbitalData);
  if (!orb) {
    printf("Not enough memory for orbitals array\n");
    Abort(1);
  }
  orb->idx = n_orbitals;
  if (orb->n == 0) {
    n_continua++;
  }
  n_orbitals++;
  AddOrbMap(orb);
  if (orbitals->lock) ReleaseLock(orbitals->lock);
#pragma omp flush
  return n_orbitals - 1;
}
*/

ORBITAL *GetOrbital(int k) {
  return (ORBITAL *) ArrayGet(orbitals, k);
}

ORBITAL *GetOrbitalSolved(int k) {
  ORBITAL *orb;
  int i;
  
  orb = (ORBITAL *) ArrayGet(orbitals, k);
  if (orb != NULL && orb->isol) return orb;  
  if (orbitals->lock) {
    SetLock(orbitals->lock);
    orb = (ORBITAL *) ArrayGet(orbitals, k);
  }
  if (orb->isol == 0) {
    MPrintf(-1, "isol0b: %d %d %d %d %g %x\n",
	    k, orbitals->dim, orb->n, orb->kappa, orb->energy, orbitals->lock);
    Abort(1);
    i = SolveDirac(orb);
    if (i < 0) {
      printf("Error occured in solving Dirac eq. err = %d\n", i);
      Abort(1);
    }
  }
  if (orbitals->lock) ReleaseLock(orbitals->lock);
  return orb;
}

ORBITAL *GetOrbitalSolvedNoLock(int k) {
  ORBITAL *orb;
  int i;
  
  orb = (ORBITAL *) ArrayGet(orbitals, k);
  if (orb != NULL && orb->isol) return orb;
  
  orb = (ORBITAL *) ArrayGet(orbitals, k);
  if (orb->isol == 0) {
    i = SolveDirac(orb);
    if (i < 0) {
      printf("Error occured in solving Dirac eq. err = %d\n", i);
      Abort(1);
    }
  }
  return orb;
}

ORBITAL *GetNewOrbitalNoLock(int n, int kappa, double e) {
  ORBITAL *orb;

  orb = (ORBITAL *) ArrayAppend(orbitals, NULL, InitOrbitalData);
  if (!orb) {
    printf("Not enough memory for orbitals array\n");
    Abort(1);
  }
  orb->idx = n_orbitals;
  orb->n = n;
  orb->kappa = kappa;
  orb->energy = e;
  n_orbitals++;
  if (n == 0) {
    n_continua++;
  }
  AddOrbMap(orb);
#pragma omp flush
  return orb;
}

ORBITAL *GetNewOrbital(int n, int kappa, double e) {
  ORBITAL *orb;

  if (orbitals->lock) SetLock(orbitals->lock);
  orb = GetNewOrbitalNoLock(n, kappa, e);
  if (orbitals->lock) ReleaseLock(orbitals->lock);
  return orb;
}

void FreeOrbitalData(void *p) {
  ORBITAL *orb;

  orb = (ORBITAL *) p;
  //RemoveOrbMap(orb);
  if (orb->wfun) free(orb->wfun);
  if (orb->phase) free(orb->phase);
  orb->wfun = NULL;
  orb->phase = NULL;
  orb->isol = 0;
  orb->ilast = -1;
  //orb->im = -1;
  if (orb->horb) {
    FreeOrbitalData(orb->horb);
    orb->horb = NULL;
  }
}

int ClearOrbitalTable(int m) {
  ORBITAL *orb;
  int i;

  if (orbitals->lock) SetLock(orbitals->lock);
  if (m == 0) {
    n_orbitals = 0;
    n_continua = 0;
    ArrayFree(orbitals, FreeOrbitalData);
    RemoveOrbMap(0);
  } else {
    for (i = n_orbitals-1; i >= 0; i--) {
      orb = GetOrbital(i);
      if (orb->n == 0) {
	n_continua--;
      }
      if (orb->n > 0) {
	n_orbitals = i+1;
	ArrayTrim(orbitals, i+1, FreeOrbitalData);
	RemoveOrbMap(1);
	break;
      }
    }
  }
  if (orbitals->lock) ReleaseLock(orbitals->lock);
  return 0;
}

int ConfigEnergy(int m, int mr, int ng, int *kg) {
  CONFIG_GROUP *g;
  CONFIG *cfg;
  int k, kk, i, md, md1, ic;
  double e0;

  if (optimize_control.mce < 0) return 0;
  if (kg == NULL) ng = 0;
  else if (ng <= 0) kg = NULL;
  if (ng == 0) {
    ng = GetNumGroups();
  }  
  md = optimize_control.mce%10;  
  md1 = 10*(optimize_control.mce/10);
  if (m == 0) {
    for (k = 0; k < ng; k++) {
      if (kg != NULL) kk = kg[k];
      else kk = k;
      g = GetGroup(kk);
      int nmax = potential->nmax-1;
      if (potential->nb > 0 && nmax < potential->nb) nmax = potential->nb;
      if (nmax > 0 && g->nmax > nmax) {
	for (i = 0; i < g->n_cfgs; i++) {
	  cfg = (CONFIG *) ArrayGet(&(g->cfg_list), i);
	  cfg->energy = 0;
	  cfg->delta = 0;
	}
	continue;
      }
      if (md == 0) {
	if (OptimizeRadial(1, &kk, -1, NULL, 1) < 0) {
	  ReinitRadial(1);
	  ClearOrbitalTable(0);
	  continue;
	}
	if (mr > 0) RefineRadial(mr, 0);
      }
      for (i = 0; i < g->n_cfgs; i++) {
	if (md > 0) {
	  if (OptimizeRadial(1, &kk, i, NULL, 1) < 0) {
	    ReinitRadial(1);
	    ClearOrbitalTable(0);
	    continue;
	  }
	  if (mr > 0) RefineRadial(mr, 0);
	}
	cfg = (CONFIG *) ArrayGet(&(g->cfg_list), i);
	nmax = potential->nmax-1;
	if (potential->nb > 0 && nmax < potential->nb) nmax = potential->nb;
	if (nmax > 0 && g->nmax > nmax) {
	  cfg->energy = 0;
	  cfg->delta = 0;
	  ReinitRadial(1);
	  ClearOrbitalTable(0);
	  continue;
	}
	e0 = AverageEnergyConfigMode(cfg, md1);
	cfg->energy = e0;
	ReinitRadial(1);
	ClearOrbitalTable(0);
      }
    }
  } else if (md1 < 20) {
    for (k = 0; k < ng; k++) {
      if (kg != NULL) kk = kg[k];
      else kk = k;
      g = GetGroup(kk);
      for (i = 0; i < g->n_cfgs; i++) {
	cfg = (CONFIG *) ArrayGet(&(g->cfg_list), i);
	if (cfg->energy != 0) {
	  e0 = AverageEnergyConfigMode(cfg, md1);
	  cfg->delta = cfg->energy - e0;	
	  if (optimize_control.iprint) {
	    MPrintf(-1, "ConfigEnergy: %d %d %d %d %g %g %g\n", 
		    md, md1, kk, i, cfg->energy, e0, cfg->delta);
	  }
	}
      }
    }
  }
  return 0;
}

double TotalEnergyGroup(int kg) {
  return TotalEnergyGroupMode(kg, 0);
}

/* calculate the total configuration average energy of a group. */
double TotalEnergyGroupMode(int kg, int md) {
  CONFIG_GROUP *g;
  ARRAY *c;
  CONFIG *cfg;
  int t;
  double total_energy, a;

  g = GetGroup(kg);
  c = &(g->cfg_list);
  
  total_energy = 0.0;
  g->sweight = 0.0;
  for (t = 0; t < g->n_cfgs; t++) {
    cfg = (CONFIG *) ArrayGet(c, t);
    a = AverageEnergyConfigMode(cfg, md);
    total_energy += a*cfg->sweight;
    g->sweight += cfg->sweight;
  }
  total_energy /= g->sweight;
  return total_energy;
}

double ZerothEnergyConfig(CONFIG *cfg) {
  int i, n, nq, kappa, k;
  double r, e;

  r = 0.0;
  for (i = 0; i < cfg->n_shells; i++) {
    n = (cfg->shells[i]).n;
    nq = (cfg->shells[i]).nq;
    kappa = (cfg->shells[i]).kappa;
    k = OrbitalIndex(n, kappa, 0.0);
    e = GetOrbital(k)->energy;
    r += nq * e;
  }
  return r;
}

double ZerothResidualConfig(CONFIG *cfg) {
  int i, n, nq, kappa, k;
  double r, e;

  r = 0.0;
  for (i = 0; i < cfg->n_shells; i++) {
    n = (cfg->shells[i]).n;
    nq = (cfg->shells[i]).nq;
    kappa = (cfg->shells[i]).kappa;
    k = OrbitalIndex(n, kappa, 0.0);
    ResidualPotential(&e, k, k);
    r += nq * e;
  }
  return r;
}
  
static double FKB(int ka, int kb, int k) {
  int ja, jb, ia, ib;
  double a, b;
  
  GetJLFromKappa(GetOrbital(ka)->kappa, &ja, &ia);
  GetJLFromKappa(GetOrbital(kb)->kappa, &jb, &ib);

  if (!Triangle(ia, k, ia) || !Triangle(ib, k, ib)) return 0.0;
  a = W3j(ja, k, ja, 1, 0, -1)*W3j(jb, k, jb, 1, 0, -1);
  if (fabs(a) < EPS30) return 0.0; 
  Slater(&b, ka, kb, ka, kb, k/2, 0);
  
  b *= a*(ja+1.0)*(jb+1.0);
  if (IsEven((ja+jb)/2)) b = -b;

  return b;
}

static double GKB(int ka, int kb, int k) {
  int ja, jb, ia, ib;
  double a, b;
  
  GetJLFromKappa(GetOrbital(ka)->kappa, &ja, &ia);
  GetJLFromKappa(GetOrbital(kb)->kappa, &jb, &ib);
  
  if (IsOdd((ia+k+ib)/2) || !Triangle(ia, k, ib)) return 0.0;
  a = W3j(ja, k, jb, 1, 0, -1);
  if (fabs(a) < EPS30) return 0.0;
  Slater(&b, ka, kb, kb, ka, k/2, 0);
  
  b *= a*a*(ja+1.0)*(jb+1.0);
  if (IsEven((ja+jb)/2)) b = -b;
  return b;
}  
  
static double ConfigEnergyVarianceParts0(SHELL *bra, int ia, int ib, 
					 int m2, int p) {
  int ja, jb, k, kp, k0, k1, kp0, kp1, ka, kb;
  double a, b, c, d, e;

  ja = GetJFromKappa(bra[ia].kappa);
  ka = OrbitalIndex(bra[ia].n, bra[ia].kappa, 0);
  if (p > 0) {
    jb = GetJFromKappa(bra[ib].kappa);
    kb = OrbitalIndex(bra[ib].n, bra[ib].kappa, 0);
  }
  e = 0.0;
  switch (p) {
  case 0:
    k0 = 4;
    k1 = 2*ja;
    a = 1.0/(ja*(ja+1.0));
    for (k = k0; k <= k1; k += 4) {
      for (kp = k0; kp <= k1; kp += 4) {
	b = -a + W6j(ja, ja, k, ja, ja, kp);
	if (k == kp) b += 1.0/(k+1.0);
	b *= a*FKB(ka, ka, k)*FKB(ka, ka, kp);
	e += b;
      }
    }
    break;
  case 1:
    k0 = 4;
    k1 = 2*ja;
    kp0 = 4;
    kp1 = 2*jb;
    kp1 = Min(kp1, k1);
    a = 1.0/(ja*(ja+1.0));
    for (k = k0; k <= k1; k += 4) {
      for (kp = kp0; kp <= kp1; kp += 4) {
	b = a - W6j(ja, ja, kp, ja, ja, k);
	if (k == kp) b -= 1.0/(k+1.0);
	b *= W6j(ja, ja, kp, jb, jb, m2);
	b /= 0.5*ja;
	b *= FKB(ka, ka, k)*FKB(ka, kb, kp);
	e += b;
      }
    }
    if (IsOdd((ja+jb+m2)/2)) e = -e;
    break;
  case 2:
    k0 = 4;
    k1 = 2*ja;
    kp0 = abs(ja-jb);
    kp1 = ja + jb;
    a = 1.0/(ja*(ja+1.0));
    for (k = k0; k <= k1; k += 4) {
      for (kp = kp0; kp <= kp1; kp += 2) {
	b = W6j(k, kp, m2, jb, ja, ja);
	b = -b*b;
	d=W6j(jb, jb, k, ja, ja, m2)*W6j(jb, jb, k, ja, ja, kp);
	if(IsOdd((m2+kp)/2)) b -= d;
	else b += d;
	if (m2 == kp) {
	  b -= (1.0/(m2+1.0)-1.0/(jb+1.0))*a;
	} else {
	  b += a/(jb+1.0);
	}
	b /= 0.5*ja;
	b *= FKB(ka, ka, k)*GKB(ka, kb, kp);
	e += b;
      }
    }
    if (IsOdd((ja+jb)/2+1)) e = -e;
    break;
  case 3:
    k0 = 4;
    k1 = 2*ja;
    k = 2*jb;
    k1 = Min(k, k1);
    for (k = k0; k <= k1; k += 4) {
      for (kp = k0; kp <= k1; kp += 4) {
	b = 0.0;
	if (k == kp) b += 1.0/((k+1.0)*(jb+1.0));
	b -= W9j(ja, ja, k, ja, m2, jb, kp, jb, jb);
	b -= W6j(ja, jb, m2, jb, ja, k)*W6j(ja, jb, m2, jb, ja, kp)/ja;
	b /= ja;
	b *= FKB(ka, kb, k)*FKB(ka, kb, kp);
	e += b;
      }
    }
    break;
  case 4:
    k0 = abs(ja-jb);
    k1 = ja+jb;
    for (k = k0; k <= k1; k += 2) {
      for (kp = k0; kp <= k1; kp += 2) {
	b = 0.0;
	if (k == kp) b += 1.0/((k+1.0)*(jb+1.0));
	b -= W9j(ja, jb, k, jb, m2, ja, kp, ja, jb);
	c = -1.0/(jb+1.0);	
	d = c;
	if (k == m2) {
	  c += 1.0/(m2+1.0);
	}
	if (kp == m2) {	  
	  d += 1.0/(m2+1.0);
	}
	b -= c*d/ja;
	b /= ja;
	b *= GKB(ka, kb, k)*GKB(ka, kb, kp);
	e += b;
      }
    }
    break;
  case 5:
    k0 = 4;
    k1 = 2*ja;
    k = 2*jb;
    k1 = Min(k1, k);
    kp0 = abs(ja-jb);
    kp1 = ja+jb;
    for (k = k0; k <= k1; k += 4) {
      for (kp = kp0; kp <= kp1; kp += 2) {
      	b=W6j(jb, jb, k, ja, ja, kp)/(jb+1.0);
      	if(IsOdd(kp/2)) b=-b;
	c = W6j(k, kp, m2, ja, jb, jb)*W6j(k, kp, m2, jb, ja, ja);
	if (IsOdd((ja+jb+kp+m2)/2)) b += c;
	else b -= c;
	c = -1.0/(jb+1.0);
	if (kp == m2) c += 1.0/(m2+1.0);
	c *= W6j(ja, jb, m2, jb, ja, k)/ja;
	if(IsOdd(m2/2)) b += c;
	else b -= c;
	b /= 0.5*ja;
	b *= FKB(ka, kb, k)*GKB(ka, kb, kp);
	e += b;
      }
    }
    break;
  case 6:
    k0 = 4;
    k1 = 2*ja;
    k = 2*jb;
    k1 = Min(k, k1);
    a = 1.0/((ja+1.0)*(jb+1.0));
    for (k = k0; k <= k1; k += 4) {
      b = a/(k+1.0);
      c = FKB(ka, kb, k);
      b *= c*c;
      e += b;
    }
    break;
  case 7:
    k0 = abs(ja-jb);
    k1 = ja+jb;
    a = 1.0/((ja+1.0)*(jb+1.0));
    for (k = k0; k <= k1; k += 2) {
      for (kp = k0; kp <= k1; kp += 2) {
	b = 0;
	if (k == kp) {
	  b += 1.0/(k+1.0);
	}
	b -= a;
	c = GKB(ka, kb, k);
	d = GKB(ka, kb, kp);
	e += a*b*c*d;
      }
    }
    break;
  case 8:    
    k0 = 4;
    k1 = 2*ja;
    k = 2*jb;
    k1 = Min(k, k1);
    kp0 = abs(ja-jb);
    kp1 = ja+jb;
    a = 1.0/((ja+1.0)*(jb+1.0));
    for (k = k0; k <= k1; k += 4) {
      for (kp = kp0; kp <= kp1; kp += 2) {
	b = W6j(jb, ja, kp, ja, jb, k);
	if (fabs(b) < EPS30) continue;
	b *= 2.0*a;
	if (IsOdd(kp/2)) b = -b;
	c = FKB(ka, kb, k);
	d = GKB(ka, kb, kp);
	e += b*c*d;
      }
    }
    break;
  }

  return e;
}

static double ConfigEnergyVarianceParts1(SHELL *bra, int i, 
					 int ia, int ib, int m2, int p) {  
  int js, ja, jb, k, kp, k0, k1, kp0, kp1, ka, kb, ks;
  double a, b, e;
  
  js = GetJFromKappa(bra[i].kappa);
  ks = OrbitalIndex(bra[i].n, bra[i].kappa, 0);
  ja = GetJFromKappa(bra[ia].kappa);
  ka = OrbitalIndex(bra[ia].n, bra[ia].kappa, 0);
  jb = GetJFromKappa(bra[ib].kappa);
  kb = OrbitalIndex(bra[ib].n, bra[ib].kappa, 0);
  e = 0.0;

  switch (p) {
  case 0:
    k0 = 4;
    k1 = 2*ja;
    k = 2*jb;
    k1 = Min(k, k1);
    k = 2*js;
    k1 = Min(k, k1);
    for (k = k0; k <= k1; k += 4) {
      b = W6j(ja, ja, k, jb, jb, m2);
      if (fabs(b) < EPS30) continue;
      b *= 2.0/((k+1.0)*(js+1.0));
      b *= FKB(ks, ka, k)*FKB(ks, kb, k);
      if (IsEven((ja+jb+m2)/2)) b = -b;
      e += b;
    }
    break;
  case 1:
    k0 = abs(js-ja);
    k1 = js+ja;
    kp0 = abs(js-jb);
    kp1 = js+jb;
    a = 1.0/((js+1.0)*(ja+1.0)*(jb+1.0));
    for (k = k0; k <= k1; k += 2) {
      for (kp = kp0; kp <= kp1; kp += 2) {
	b = W6j(k, kp, m2, jb, ja, js);
	b = -b*b + a;
	b /= (js+1.0);
	b *= 2.0*GKB(ks, ka, k)*GKB(ks, kb, kp);
	if (IsOdd((ja+jb)/2+1)) b = -b;
	e += b;
      }
    }
    break;
  case 2:
    k0 = 4;
    k1 = 2*js;    
    k = 2*ja;
    k1 = Min(k, k1);
    kp0 = abs(js-jb);
    kp1 = js+jb;
    for (k = k0; k <= k1; k += 4) {
      for (kp = kp0; kp <= kp1; kp += 2) {
	b = W6j(jb, jb, k, ja, ja, m2);
	b *= W6j(jb, jb, k, js, js, kp);
	if (fabs(b) < EPS30) continue;
	b /= (js+1.0);
	if (IsEven((ja+jb+kp+m2)/2)) b = -b;
	b *= 2.0*FKB(ks, ka, k)*GKB(ks, kb, kp);
	e += b;
      }
    }
    break;
  }

  return e;
}

double ConfigEnergyVariance(int ns, SHELL *bra, int ia, int ib, int m2) {
  int i, js, p;
  double e, a, b, c;
  
  e = 0.0;
  for (i = 0; i < ns; i++) {
    js = GetJFromKappa(bra[i].kappa);
    a = bra[i].nq;
    b = js+1.0 - bra[i].nq;
    if (i == ia) {
      a -= 1.0;      
    }
    if (i == ib) {
      b -= 1.0;
    }
    if (a == 0.0 || b == 0.0) continue;
    a = a*b;
    b = 0.0;
    if (i == ia) {
      for (p = 0; p < 6; p++) {
	c = ConfigEnergyVarianceParts0(bra, ia, ib, m2, p);
	b += c;
      }
      b /= js-1.0;
    } else if (i == ib) {
      for (p = 0; p < 6; p++) {
	c = ConfigEnergyVarianceParts0(bra, ib, ia, m2, p);
	b += c;
      }
      b /= js-1.0;
    } else {
      for (p = 6; p < 9; p++) {
	c = ConfigEnergyVarianceParts0(bra, i, ia, m2, p);
	b += c;
	c = ConfigEnergyVarianceParts0(bra, i, ib, m2, p);
	b += c;
      }
      c = ConfigEnergyVarianceParts1(bra, i, ia, ib, m2, 0);
      b += c;
      c = ConfigEnergyVarianceParts1(bra, i, ia, ib, m2, 1);
      b += c;
      c = ConfigEnergyVarianceParts1(bra, i, ia, ib, m2, 2);
      b += c;
      c = ConfigEnergyVarianceParts1(bra, i, ib, ia, m2, 2);
      b += c;
      b /= js;
    }

    e += a*b;
  }
    
  if (e < 0.0) e = 0.0;
  return e;
}

double ConfigEnergyShiftCI(int nrs0, int nrs1) {
  int n0, k0, q0, n1, k1, q1;
  int j0, j1, s0, s1, kmax, k;
  double a, g, pk, w, sd, q;

  UnpackNRShell(&nrs0, &n0, &k0, &q0);
  UnpackNRShell(&nrs1, &n1, &k1, &q1);
  q = (q0-1.0)/(2.0*k0+1.0) - q1/(2.0*k1+1.0);
  if (q == 0.0) return q;
  g = 0.0;
  for (j0 = k0-1; j0 <= k0+1; j0 += 2) {
    if (j0 < 0) continue;
    s0 = OrbitalIndex(n0, GetKappaFromJL(j0, k0), 0);    
    for (j1 = k1-1; j1 <= k1+1; j1 += 2) {
      if (j1 < 0) continue;
      s1 = OrbitalIndex(n1, GetKappaFromJL(j1, k1), 0);
      w = W6j(j0, k0, 1, k1, j1, 2);
      w = w*w*(j0+1.0)*(j1+1.0)*0.5;
      pk = ((j0+1.0)*(j1+1.0))/((k0+1.0)*(k1+1.0)*4.0) - w;
      Slater(&sd, s0, s1, s0, s1, 0, 0);
      g += sd*pk;
      kmax = 2*j0;
      k = 2*j1;
      kmax = Min(k, kmax);
      for (k = 4; k <= kmax; k += 4) {	
	pk = W3j(k0, k, k0, 0, 0, 0);
	if (fabs(pk) > 0) {
	  pk *= W3j(k1, k, k1, 0, 0, 0);
	  if (fabs(pk) > 0) {
	    pk *= 0.25;
	    a = W6j(j0, 2, j1, k1, 1, k0);
	    if (fabs(a) > 0) {
	      a = a*a;
	      a *= W6j(j0, k, j0, k0, 1, k0);
	      if (fabs(a) > 0) {
		a *= W6j(j1, k, j1, k1, 1, k1);
		if (fabs(a) > 0) {
		  a *= W6j(j0, k, j0, j1, 2, j1);
		  a *= 2.0*(j0+1.0)*(j1+1.0)*(k0+1.0)*(k1+1.0);
		}
	      }
	    }
	    a += W6j(k0, k, k0, k1, 2, k1);
	    pk *= a;
	  }
	  if (fabs(pk) > 0) {
	    Slater(&sd, s0, s1, s0, s1, k/2, 0);
	    g += pk*sd*(j0+1.0)*(j1+1.0);
	  }
	}
      }
      pk = W3j(k0, 2, k1, 0, 0, 0);
      if (fabs(pk) > 0) {
	pk = pk*pk;
	a = W6j(j0, 2, j1, k1, 1, k0);
	a *= a;
	pk *= a;
	pk *= (k0+1.0)*(k1+1.0)/3.0;
	pk *= (1.0 - 0.5*(j0+1.0)*(j1+1.0)*a);
	if (fabs(pk) > 0) {
	  Slater(&sd, s0, s1, s1, s0, 1, 0);
	  g += pk*sd*(j0+1.0)*(j1+1.0);
	}
      }
    }
  }

  return q*g;
}
	  
double ConfigEnergyShift(int ns, SHELL *bra, int ia, int ib, int m2) {
  double qa, qb, a, b, c, sd, e;
  int ja, jb, k, kmin, kmax;
  int k0, k1;

  qa = bra[ia].nq;
  qb = bra[ib].nq;
  ja = GetJFromKappa(bra[ia].kappa);
  jb = GetJFromKappa(bra[ib].kappa);
  if (qa == 1 && qb == 0) e = 0.0;
  else {
    e = (qa-1.0)/ja - qb/jb;
    if (e != 0.0) {
      kmin = 4;
      kmax = 2*ja;
      k = 2*jb;
      kmax = Min(kmax, k);
      a = 0.0;
      k0 = OrbitalIndex(bra[ia].n, bra[ia].kappa, 0);
      k1 = OrbitalIndex(bra[ib].n, bra[ib].kappa, 0);
      for (k = kmin; k <= kmax; k += 4) {
	b = W6j(k, ja, ja, m2, jb, jb);	
	if (fabs(b) > EPS30) {
	  a += b*FKB(k0, k1, k);
	}
      }
      
      if(IsOdd(m2/2)) a=-a;
      
      kmin = abs(ja-jb);
      kmax = ja + jb;
      c = 1.0/((ja+1.0)*(jb+1.0));
      for (k = kmin; k <= kmax; k += 2) {
	if (k == m2) {	  
	  b = 1.0/(m2+1.0)-c;
	} else {
	  b = -c;
	}
	a += b*GKB(k0, k1, k);
      }
      if (IsEven((ja+jb)/2)) a = -a;
      e *= a;
    }
  }

  return e;
}

/*
** when enabled, qr_norm stores the shift of orbital energies for MBPT
*/
void ShiftOrbitalEnergy(CONFIG *cfg) {
  int i, j;
  ORBITAL *orb;
  CONFIG c;
  double e0, e1;

  c.n_shells = cfg->n_shells+1;
  c.shells = malloc(sizeof(SHELL)*c.n_shells);
  memcpy(c.shells+1, cfg->shells, sizeof(SHELL)*cfg->n_shells);
  c.shells[0].nq = 1;
  c.shells[1].nq -= 1;
  for (i = 0; i < n_orbitals; i++) {
    orb = GetOrbital(i);
    orb->qr_norm = 0.0;
    /* this disables the shift */
    continue;
    for (j = 0; j < cfg->n_shells; j++) {
      if (cfg->shells[j].n == orb->n && cfg->shells[j].kappa == orb->kappa) {
	break;
      }
    }
    ResidualPotential(&e0, i, i);
    if (j < cfg->n_shells) {
      e1 = CoulombEnergyShell(cfg, j);
    } else {
      c.shells[0].n = orb->n;
      c.shells[0].kappa = orb->kappa;
      e1 = CoulombEnergyShell(&c, 0);
    }
    orb->qr_norm = 2*e1 + e0;
  }
  free(c.shells);
}
    
double CoulombEnergyShell(CONFIG *cfg, int i) {
  int n, kappa, kl, j2, nq, k, np, kappap, klp, j2p;
  int nqp, kp, kk, j, kkmin, kkmax;
  double y, t, q, b, a, r;

  n = (cfg->shells[i]).n;
  kappa = (cfg->shells[i]).kappa;
  kl = GetLFromKappa(kappa);
  j2 = GetJFromKappa(kappa);
  nq = (cfg->shells[i]).nq;
  k = OrbitalIndex(n, kappa, 0.0);
    
  if (nq > 1) {
    t = 0.0;
    for (kk = 2; kk <= j2; kk += 2) {
      Slater(&y, k, k, k, k, kk, 0);
      q = W3j(j2, 2*kk, j2, -1, 0, 1);
      t += y * q * q ;
    }
    Slater(&y, k, k, k, k, 0, 0);
    b = (nq-1.0) * (y - (1.0 + 1.0/j2)*t);

  } else {
    b = 0.0;
  }

  t = 0.0;
  for (j = 0; j < cfg->n_shells; j++) {
    if (j == i) continue;
    nqp = (cfg->shells[j]).nq;
    if (nqp == 0) continue;
    np = (cfg->shells[j]).n;
    kappap = (cfg->shells[j]).kappa;
    klp = GetLFromKappa(kappap);
    j2p = GetJFromKappa(kappap);
    kp = OrbitalIndex(np, kappap, 0.0);
    
    kkmin = abs(j2 - j2p);
    kkmax = (j2 + j2p);
    if (IsOdd((kkmin + kl + klp)/2)) kkmin += 2;
    a = 0.0;
    for (kk = kkmin; kk <= kkmax; kk += 4) {
      Slater(&y, k, kp, kp, k, kk/2, 0);
      q = W3j(j2, kk, j2p, -1, 0, 1);
      a += y * q * q;
    }
    Slater(&y, k, kp, k, kp, 0, 0);
    t += nqp * (y - a);
  }
  
  r = 0.5*(b + t);
  
  return r;
}

double AverageEnergyConfig(CONFIG *cfg) {
  return AverageEnergyConfigMode(cfg, 0);
}

/* calculate the average energy of a configuration */
double AverageEnergyConfigMode(CONFIG *cfg, int md) {
  int i, j, n, kappa, nq, np, kappap, nqp;
  int k, kp, kk, kl, klp, kkmin, kkmax, j2, j2p;
  double x, y, t, q, a, b, r, e, *v1, *v2;
 
  x = 0.0;
  for (i = 0; i < cfg->n_shells; i++) {
    n = (cfg->shells[i]).n;
    kappa = (cfg->shells[i]).kappa;
    kl = GetLFromKappa(kappa);
    j2 = GetJFromKappa(kappa);
    nq = (cfg->shells[i]).nq;
    k = OrbitalIndex(n, kappa, 0.0);
    if (md < 10) {
      double am = AMU * GetAtomicMass();
      if (nq > 1) {
	t = 0.0;
	for (kk = 1; kk <= j2; kk += 1) {
	  y = 0;
	  if (kk == 1 && qed.sms) {
	    v1 = Vinti(k, k);
	  } else {
	    v1 = NULL;
	  }
	  if (IsEven(kk)) {
	    Slater(&y, k, k, k, k, kk, 0);
	  }
	  if (v1 && qed.sms == 3) {
	    double a1 = ReducedCL(j2, 2*kk, j2);
	    if (fabs(a1) > 1e-10) {
	      y += 0.25*v1[2]*v1[2]/(am*a1*a1);
	    }	  
	    y += 0.25*v1[1]*v1[1]/am;
	  }
	  if (qed.br < 0 || n <= qed.br) {
	    if (qed.minbr <= 0 || n <= qed.minbr) {
	      int mbr = qed.mbr;
	      if (qed.nbr > 0 && n > qed.nbr) mbr = 0;
	      y += Breit(k, k, k, k, kk, kappa, kappa, kappa, kappa,
			 kl, kl, kl, kl, mbr);
	    }
	  }
	  if (y) {
	    q = W3j(j2, 2*kk, j2, -1, 0, 1);
	    t += y * q * q ;
	  }
	}
	Slater(&y, k, k, k, k, 0, 0);
	b = ((nq-1.0)/2.0) * (y - (1.0 + 1.0/j2)*t);
      } else {
	b = 0.0;
      }
      
      t = 0.0;
      for (j = 0; j < i; j++) {
	np = (cfg->shells[j]).n;
	kappap = (cfg->shells[j]).kappa;
	klp = GetLFromKappa(kappap);
	j2p = GetJFromKappa(kappap);
	nqp = (cfg->shells[j]).nq;
	kp = OrbitalIndex(np, kappap, 0.0);
	kkmin = abs(j2 - j2p);
	kkmax = (j2 + j2p);
	int maxn, minn;
	if (n < np) {
	  minn = n;
	  maxn = np;
	} else {
	  minn = np;
	  maxn = n;
	}
	//if (IsOdd((kkmin + kl + klp)/2)) kkmin += 2;
	a = 0.0;
	for (kk = kkmin; kk <= kkmax; kk += 2) {
	  y = 0;
	  int kk2 = kk/2;
	  if (kk == 2 && qed.sms) {
	    v1 = Vinti(k, kp);
	    v2 = Vinti(kp, k);
	  } else {
	    v1 = NULL;
	    v2 = NULL;
	  }
	  if (IsEven((kl+klp+kk)/2)) {
	    Slater(&y, k, kp, kp, k, kk2, 0);
	    if (v1 && v2) {	    
	      y -= (v1[0]+v1[2])*v2[0]/am;
	    }
	  }
	  if (v1 && v2) {
	    double a1 = ReducedCL(j2, kk, j2p);
	    double a2 = ReducedCL(j2p, kk, j2);
	    if (IsEven((kl+klp+kk)/2)) {
	      y -= v1[1]*v2[0]/am;
	      if (qed.sms == 3 && fabs(a1) > 1e-10) {
		y += 0.5*(0.5+a2)*v1[2]*v2[2]/(am*a1*a2);
	      }
	    } else {
	      if (qed.sms == 3 && fabs(a1) > 1e-10) {
		y += 0.25*v1[2]*v2[2]/(am*a1*a2);
	      }
	    }
	    if (qed.sms == 3) {
	      y += 0.25*v1[1]*v2[1]/am;
	    }
	  }
	  if (qed.br < 0 || maxn <= qed.br) {
	    if (qed.minbr <= 0 || minn <= qed.minbr) {
	      int mbr = qed.mbr;
	      if (qed.nbr > 0 && maxn > qed.nbr) mbr = 0;
	      y += Breit(k, kp, kp, k, kk2, kappa, kappap, kappap, kappa,
			 kl, klp, klp, kl, mbr);
	    }
	  }
	  if (y) {
	    q = W3j(j2, kk, j2p, -1, 0, 1);
	    a += y * q * q;
	  }
	}
	y = 0;
	Slater(&y, k, kp, k, kp, 0, 0);
	t += nqp * (y - a);
      }    
      ResidualPotential(&y, k, k);
      e = GetOrbital(k)->energy;
      a = QED1E(k, k);
      r = nq * (b + t + e + a + y);      
    } else {
      ORBITAL *orb = GetOrbitalSolved(k);
      a = SelfEnergy(orb, orb);
      r = nq * a;
    }
    x += r;
  }
  return x;
}

void DiExConfig(CONFIG *cfg, double *d0, double *d1) {
  int i, j, n, kappa, np, kappap;
  int k, kp, kk, kl, klp, kkmin, kkmax, j2, j2p;
  double y, t, q, a, b, nq, nqp;
  SHELL *s1, *s2;
 
  *d0 = 0.0;
  *d1 = 0.0;
  for (i = 0; i < cfg->n_shells; i++) {
    s1 = &cfg->shells[i];
    n = s1->n;
    kappa = s1->kappa;
    kl = GetLFromKappa(kappa);
    j2 = GetJFromKappa(kappa);
    nq = s1->nq;
    k = OrbitalIndex(n, kappa, 0.0);

    t = 0.0;
    for (kk = 2; kk <= j2; kk += 2) {
      Slater(&y, k, k, k, k, kk, 0);
      q = W3j(j2, 2*kk, j2, -1, 0, 1);
      t += y * q * q ;
    }
    Slater(&y, k, k, k, k, 0, 0);
    b = ((nq-1.0)/2.0);
    *d0 += nq*b*(y - (1.0+1.0/j2)*t);
    
#if FAC_DEBUG
      fprintf(debug_log, "\nAverage Radial: %lf\n", y);
#endif
      
    for (j = 0; j < i; j++) {
      s2 = &cfg->shells[j];
      np = s2->n;
      kappap = s2->kappa;
      klp = GetLFromKappa(kappap);
      j2p = GetJFromKappa(kappap);
      nqp = s2->nq;
      kp = OrbitalIndex(np, kappap, 0.0);

      kkmin = abs(j2 - j2p);
      kkmax = (j2 + j2p);
      if (IsOdd((kkmin + kl + klp)/2)) kkmin += 2;
      a = 0.0;
      for (kk = kkmin; kk <= kkmax; kk += 4) {
	Slater(&y, k, kp, kp, k, kk/2, 0);
	q = W3j(j2, kk, j2p, -1, 0, 1);
	a += y * q * q;
#if FAC_DEBUG
	fprintf(debug_log, "exchange rank: %d, q*q: %lf, Radial: %lf\n", 
		kk/2, q*q, y);
#endif

      }
      Slater(&y, k, kp, k, kp, 0, 0);

#if FAC_DEBUG
      fprintf(debug_log, "direct: %lf\n", y);
#endif
      *d0 += nq*nqp*y;
      *d1 += -nq*nqp*a;
    }
  }
}

/* calculate the average energy of an average configuration, 
** with seperate direct and exchange contributions */
void DiExAvgConfig(AVERAGE_CONFIG *cfg, double *d0, double *d1) {
  int i, j, n, kappa, np, kappap;
  int k, kp, kk, kl, klp, kkmin, kkmax, j2, j2p;
  double y, t, q, a, b, nq, nqp;
 
  *d0 = 0.0;
  *d1 = 0.0;
  for (i = 0; i < cfg->n_shells; i++) {
    n = cfg->n[i];
    kappa = cfg->kappa[i];
    kl = GetLFromKappa(kappa);
    j2 = GetJFromKappa(kappa);
    nq = cfg->nq[i];
    k = OrbitalIndex(n, kappa, 0.0);

    t = 0.0;
    for (kk = 2; kk <= j2; kk += 2) {
      Slater(&y, k, k, k, k, kk, 0);
      q = W3j(j2, 2*kk, j2, -1, 0, 1);
      t += y * q * q ;
    }
    Slater(&y, k, k, k, k, 0, 0);
    b = ((nq-1.0)/2.0);
    *d0 += nq*b*(y - (1.0+1.0/j2)*t);
    
#if FAC_DEBUG
      fprintf(debug_log, "\nAverage Radial: %lf\n", y);
#endif
      
    for (j = 0; j < i; j++) {
      np = cfg->n[j];
      kappap = cfg->kappa[j];
      klp = GetLFromKappa(kappap);
      j2p = GetJFromKappa(kappap);
      nqp = cfg->nq[j];
      kp = OrbitalIndex(np, kappap, 0.0);

      kkmin = abs(j2 - j2p);
      kkmax = (j2 + j2p);
      if (IsOdd((kkmin + kl + klp)/2)) kkmin += 2;
      a = 0.0;
      for (kk = kkmin; kk <= kkmax; kk += 4) {
	Slater(&y, k, kp, kp, k, kk/2, 0);
	q = W3j(j2, kk, j2p, -1, 0, 1);
	a += y * q * q;
#if FAC_DEBUG
	fprintf(debug_log, "exchange rank: %d, q*q: %lf, Radial: %lf\n", 
		kk/2, q*q, y);
#endif

      }
      Slater(&y, k, kp, k, kp, 0, 0);

#if FAC_DEBUG
      fprintf(debug_log, "direct: %lf\n", y);
#endif
      *d0 += nq*nqp*y;
      *d1 += -nq*nqp*a;
    }
  }
}

/* calculate the average energy of an average configuration */
double AverageEnergyAvgConfig(AVERAGE_CONFIG *cfg) {
  int i, j, n, kappa, np, kappap;
  int k, kp, kk, kl, klp, kkmin, kkmax, j2, j2p;
  double x, y, t, q, a, b, r, nq, nqp, r0, r1;
 
  r0 = 0.0;
  r1 = 0.0;
  x = 0.0;
  for (i = 0; i < cfg->n_shells; i++) {
    n = cfg->n[i];
    kappa = cfg->kappa[i];
    kl = GetLFromKappa(kappa);
    j2 = GetJFromKappa(kappa);
    nq = cfg->nq[i];
    k = OrbitalIndex(n, kappa, 0.0);
    
    t = 0.0;
    for (kk = 2; kk <= j2; kk += 2) {
      Slater(&y, k, k, k, k, kk, 0);
      q = W3j(j2, 2*kk, j2, -1, 0, 1);
      t += y * q * q ;
    }
    Slater(&y, k, k, k, k, 0, 0);
    b = ((nq-1.0)/2.0) * (y - (1.0 + 1.0/j2)*t);
    
#if FAC_DEBUG
      fprintf(debug_log, "\nAverage Radial: %lf\n", y);
#endif
      
    t = 0.0;
    for (j = 0; j < i; j++) {
      np = cfg->n[j];
      kappap = cfg->kappa[j];
      klp = GetLFromKappa(kappap);
      j2p = GetJFromKappa(kappap);
      nqp = cfg->nq[j];
      kp = OrbitalIndex(np, kappap, 0.0);

      kkmin = abs(j2 - j2p);
      kkmax = (j2 + j2p);
      if (IsOdd((kkmin + kl + klp)/2)) kkmin += 2;
      a = 0.0;
      for (kk = kkmin; kk <= kkmax; kk += 4) {
	Slater(&y, k, kp, kp, k, kk/2, 0);
	q = W3j(j2, kk, j2p, -1, 0, 1);
	a += y * q * q;
#if FAC_DEBUG
	fprintf(debug_log, "exchange rank: %d, q*q: %lf, Radial: %lf\n", 
		kk/2, q*q, y);
#endif

      }
      Slater(&y, k, kp, k, kp, 0, 0);

#if FAC_DEBUG
      fprintf(debug_log, "direct: %lf\n", y);
#endif

      t += nqp * (y - a);
    }

    ResidualPotential(&y, k, k);
    a = GetOrbital(k)->energy;
    r = nq * (b + t + a + y);
    r0 += nq*y;
    r1 += nq*(b+t);
    x += r;
  }

  /*printf("%12.5E %12.5E %15.8E\n", r0, r1, x);*/
  return x;
}

double ZerothHamilton(ORBITAL *orb0, ORBITAL *orb1) {
  double *p0, *p1, *q0, *q1;
  if (orb0->kappa != orb1->kappa) return 0.0;
  
  p0 = Large(orb0);
  p1 = Large(orb1);
  q0 = Small(orb0);
  q1 = Small(orb1);
  int ilast = Min(orb0->ilast, orb1->ilast);
  int i;
  double v, a, a2, b;
  Differential(p1, _xk, 0, ilast, potential->dr_drho);
  Differential(q1, _zk, 0, ilast, potential->dr_drho);
  a = 1/FINE_STRUCTURE_CONST;
  a2 = 1/FINE_STRUCTURE_CONST2;
  for (i = 0; i <= ilast; i++) {
    v = potential->Vc[i] + potential->U[i];
    b = orb1->kappa/potential->rad[i];
    _yk[i] = p0[i]*p1[i]*(v) + a*p0[i]*(-_zk[i]+q1[i]*b);
    _yk[i] += q0[i]*q1[i]*(v-2*a2) + a*q0[i]*(_xk[i]+p1[i]*b);
    _yk[i] *= potential->dr_drho[i];
    
  }
  _xk[0] = 0;
  NewtonCotes(_xk, _yk, 0, ilast, 1, 0);
  b = _xk[ilast];
  //b = Simpson(_yk, 0, ilast);
  return b;
}

double ConfigHamilton(CONFIG *cfg, ORBITAL *orb0, ORBITAL *orb1, double xdf) {
  double t, y, a, nqp, q, kappa;
  int k0, k1, j2, kl, j2p, klp, minn, maxn;
  int np, kp, kappap, kk, kkmin, kkmax, kk2, j;

  if (orb0->kappa != orb1->kappa) return 0;
  k0 = orb0->idx;
  k1 = orb1->idx;
  kappa = orb0->kappa;  
  GetJLFromKappa(kappa, &j2, &kl);
  if (orb1->n > orb0->n) {
    maxn = orb1->n;
    minn = orb0->n;
  } else {
    maxn = orb0->n;
    minn = orb1->n;
  }
  t = 0.0;
  double te = 0.0;
  double td = 0.0;
  for (j = 0; j < cfg->n_shells; j++) {
    np = (cfg->shells[j]).n;
    kappap = (cfg->shells[j]).kappa;
    kp = OrbitalIndex(np, kappap, 0);
    GetJLFromKappa(kappap, &j2p, &klp);
    nqp = (cfg->shells[j]).nq;
    if (kappap == kappa) {
      if (np == orb0->n) nqp -= xdf;
      if (np == orb1->n) nqp -= xdf;
    }
    if (nqp <= 0) continue;
    kkmin = abs(j2 - j2p);
    kkmax = (j2 + j2p);
    if (np > maxn) maxn = np;
    if (np < minn) minn = np;
    a = 0.0;
    for (kk = kkmin; kk <= kkmax; kk += 2) {
      y = 0;
      int kk2 = kk/2;
      if (IsEven((kl + klp + kk)/2)) {
	Slater(&y, k0, kp, kp, k1, kk2, 0);
      }
      /*
      if (qed.br < 0 || maxn <= qed.br) {
	if (qed.minbr <= 0 || minn <= qed.minbr) {
	  int mbr = qed.mbr;
	  if (qed.nbr > 0 && maxn > qed.nbr) mbr = 0;
	  y += Breit(k0, kp, kp, k0, kk, kappa, kp, kp, kappa,
		     kl, klp, klp, kl, mbr);
	}
      }
      */
      if (y) {
	q = W3j(j2, kk, j2p, -1, 0, 1);
	a += y*q*q;
      }
    }
    te -= nqp*a;
    y = 0;
    Slater(&y, k0, kp, k1, kp, 0, 0);
    td += nqp*y;
  }
  y = 0;
  ResidualPotential(&y, k0, k1);
  t = td + te + y;

  return t;
}

/* calculate the expectation value of the residual potential:
   -Z/r - v0(r), where v0(r) is central potential used to solve 
   dirac equations. the orbital index must be valid, i.e., their 
   radial equations must have been solved. */
int ResidualPotential(double *s, int k0, int k1) {
  int i;
  ORBITAL *orb1, *orb2;
  int index[2];
  LOCK *lock = NULL;
  double *p, z, *p1, *p2, *q1, *q2;

  orb1 = GetOrbitalSolved(k0);
  orb2 = GetOrbitalSolved(k1);
  if (!orb1 || !orb2) return -1;
  if (orb1->wfun == NULL || orb2->wfun == NULL) {
    *s = 0.0;
    return 0;
  }

  if (k0 > k1) {
    index[0] = k1;
    index[1] = k0;
  } else {
    index[0] = k0;
    index[1] = k1;
  }
  int myrank = MyRankMPI()+1;
  p = (double *) MultiSet(residual_array, index, NULL, &lock,
			  InitDoubleData, NULL);
  int locked = 0;
  if (lock && !(p && *p)) {
    SetLock(lock);
    locked = 1;
  }
  if (p && *p) {
    *s = *p;
    if (locked) ReleaseLock(lock);
#pragma omp atomic
    residual_array->iset -= myrank;
#pragma omp flush
    return 0;
  } 

  *s = 0.0;
  if (orb1->n < 0 || orb2->n < 0) {
    p1 = Large(orb1);
    p2 = Large(orb2);
    q1 = Small(orb1);
    q2 = Small(orb2);
    for (i = potential->ib; i <= potential->ib1; i++) {
      z = potential->Vc[i] + potential->U[i];
      _yk[i] = -(potential->Z[i]/potential->rad[i]) - z;
      _yk[i] *= potential->dr_drho[i];
      _yk[i] *= p1[i]*p2[i] + q1[i]*q2[i];
    }
    *s = Simpson(_yk, potential->ib, potential->ib1);
  } else {
    for (i = 0; i < potential->maxrp; i++) {
      z = potential->Vc[i] + potential->U[i];
      _yk[i] = -(potential->Z[i]/potential->rad[i]) - z;
    }
    Integrate(_yk, orb1, orb2, 1, s, -1);
    if (potential->nfrozen || potential->npseudo) {
      z = ZerothHamilton(orb1, orb2);
      *s += z;
      if (orb1->n == orb2->n) *s -= orb1->energy;
    }
  }
  *p = *s;
  if (locked) ReleaseLock(lock);
#pragma omp atomic
  residual_array->iset -= myrank;
#pragma omp flush
  return 0;
}

double MeanPotential(int k0, int k1) {
  int i;
  ORBITAL *orb1, *orb2;
  double z, *p1, *p2, *q1, *q2;

  orb1 = GetOrbitalSolved(k0);
  orb2 = GetOrbitalSolved(k1);
  if (!orb1 || !orb2) return -1;
  if (orb1->wfun == NULL || orb2->wfun == NULL) {
    return 0.0;
  }

  if (orb1->n < 0 || orb2->n < 0) {
    p1 = Large(orb1);
    p2 = Large(orb2);
    q1 = Small(orb1);
    q2 = Small(orb2);
    for (i = potential->ib; i <= potential->ib1; i++) {
      z = potential->U[i];
      z += potential->Vc[i];
      _yk[i] = z;
      _yk[i] *= potential->dr_drho[i];
      _yk[i] *= p1[i]*p2[i] + q1[i]*q2[i];
    }
    z = Simpson(_yk, potential->ib, potential->ib1);
  } else {
    for (i = 0; i < potential->maxrp; i++) {
      z = potential->U[i];
      z += potential->Vc[i];
      _yk[i] = z;
    }
    Integrate(_yk, orb1, orb2, 1, &z, -1);
  }

  return z;
}

double RadialMoments(int m, int k1, int k2) {
  int index[3];
  LOCK *lock = NULL;
  int npts, i0, i;
  ORBITAL *orb1, *orb2;
  double *q, r, z, *p1, *p2, *q1, *q2;
  int n1, n2;
  int kl1, kl2;
  int nh, klh;
  
  orb1 = GetOrbitalSolved(k1);
  orb2 = GetOrbitalSolved(k2);
  n1 = orb1->n;
  n2 = orb2->n;
  kl1 = orb1->kappa;
  kl2 = orb2->kappa;
  kl1 = GetLFromKappa(kl1);
  kl2 = GetLFromKappa(kl2);
  kl1 /= 2;
  kl2 /= 2;	

  GetHydrogenicNL(&nh, &klh, NULL, NULL);

  if (n1 > 0 && n2 > 0 && potential->ib <= 0) {
    if ((n1 > nh && n2 > nh) || 
	(kl1 > klh && kl2 > klh) ||
	orb1->wfun == NULL || 
	orb2->wfun == NULL) {
      if (n1 == n2 && kl1 == kl2) {
	z = GetResidualZ();
	r = HydrogenicExpectation(z, m, n1, kl1);
	if (r) {
	  return r;
	}
      } else if (m == 1) {
	z = GetResidualZ();
	if (n1 < n2) {
	  r = HydrogenicDipole(z, n1, kl1, n2, kl2);
	  return r;
	} else if (n1 < n2) {
	  r = HydrogenicDipole(z, n2, kl2, n1, kl1);
	  return r;
	}
      }
    }
  }

  if (potential->ib <= 0 && n1 == n2 && m > 0 && n1 > potential->nmax) {
    return 0.0;
  }
  if (orb1->wfun == NULL || orb2->wfun == NULL) {
    return 0.0;
  }

  if (m >= 0) {
    index[0] = 2*m;
  } else {
    index[0] = -2*m-1;
  }
    
  if (k1 < k2) {
    index[1] = k1;
    index[2] = k2;
  } else {
    index[1] = k2;
    index[2] = k1;
  }
  int myrank = MyRankMPI()+1;
  q = (double *) MultiSet(moments_array, index, NULL, &lock,
			  InitDoubleData, NULL);
  int locked = 0;
  if (lock && !(*q)) {
    SetLock(lock);
    locked = 1;
  }
  if (*q) {
    if (locked) ReleaseLock(lock);
#pragma omp atomic
    moments_array->iset -= myrank;
#pragma omp flush
    return *q;
  } 
  if (n1 < 0 || n2 < 0) {
    i0 = potential->ib;
    npts = potential->ib1;
    p1 = Large(orb1);
    q1 = Small(orb1);
    p2 = Large(orb2);
    q2 = Small(orb2);
    for (i = i0; i <= npts; i++) {
      r = p1[i]*p2[i] + q1[i]*q2[i];
      r *= potential->dr_drho[i];
      _yk[i] = r;
      if (m != 0) _yk[i] *= pow(potential->rad[i], m);
    }
    r = Simpson(_yk, i0, npts);
    *q = r;
  } else {    
    npts = potential->maxrp-1;
    if (n1 != 0) npts = Min(npts, orb1->ilast);
    if (n2 != 0) npts = Min(npts, orb2->ilast);
    if (m == 0) {
      for (i = 0; i <= npts; i++) {
	_yk[i] = 1.0;
      }
    } else {
      for (i = 0; i <= npts; i++) {
	_yk[i] = pow(potential->rad[i], m);
      }
    }
    r = 0.0;
    Integrate(_yk, orb1, orb2, 1, &r, m);
    *q = r;
  }
  if (locked) ReleaseLock(lock);
#pragma omp atomic
    moments_array->iset -= myrank;
#pragma omp flush
  return r;
}

double MultipoleRadialNR(int m, int k1, int k2, int gauge) {
  int i, p, t;
  ORBITAL *orb1, *orb2;
  double r;
  int kappa1, kappa2;

#ifdef PERFORM_STATISTICS
  clock_t start, stop; 
  start = clock();
#endif

  orb1 = GetOrbital(k1);
  orb2 = GetOrbital(k2);
  kappa1 = orb1->kappa;
  kappa2 = orb2->kappa;

  r = 0.0;
  if (m == 1) {
    /* use the relativistic version is just simpler */
    printf("should call MultipoleRadialFR instead\n");
  } else if (m > 1) {
    t = kappa1 + kappa2;
    p = m - t;
    if (p && t) {
      r = RadialMoments(m-1, k1, k2);
      r *= p*t;
      r /= sqrt(m*(m+1.0));
      r *= -0.5 * FINE_STRUCTURE_CONST;
      for (i = 2*m - 1; i > 0; i -= 2) r /= (double) i;
      r *= ReducedCL(GetJFromKappa(kappa1), 2*m, GetJFromKappa(kappa2));
    }
  } else if (m < 0) {
    m = -m;
    if (gauge == G_BABUSHKIN) {
      r = RadialMoments(m, k1, k2);
      r *= sqrt((m+1.0)/m);
      for (i = 2*m - 1; i > 1; i -= 2) r /= (double) i;
    } else {
      /* the velocity form is not implemented yet. 
	 the following is still the length form */
      r = RadialMoments(m, k1, k2);
      r *= sqrt((m+1.0)/m);
      for (i = 2*m - 1; i > 1; i -= 2) r /= (double) i;
    }
    r *= ReducedCL(GetJFromKappa(kappa1), 2*m, GetJFromKappa(kappa2));
  }

#ifdef PERFORM_STATISTICS 
  stop = clock();
  rad_timing.radial_1e += stop - start;
#endif
  
  return r;
}

int MultipoleRadialFRGrid(double **p0, int m, int k1, int k2, int gauge) {
  double q, ip, ipm, im, imm;
  int kappa1, kappa2;
  int am, t;
  int index[4], s;
  ORBITAL *orb1, *orb2, *orb;
  double x, a, r, rp, ef, **p1;
  int jy, n, i, j, npts;
  double rcl;

#ifdef PERFORM_STATISTICS 
  clock_t start, stop;
  start = clock();
#endif

  if (m == 0) return 0;
  
  if (m >= 0) {
    index[0] = 2*m;
    am = m;
  } else {
    index[0] = -2*m-1;
    am = -m;
  }
  index[1] = k1;
  index[2] = k2;
  LOCK *lock = NULL;
  int myrank = MyRankMPI()+1;
  p1 = (double **) MultiSet(multipole_array, index, NULL, &lock,
			    InitPointerData, FreeMultipole);
  int locked = 0;
  if (lock && !(*p1)) {
    SetLock(lock);
    locked = 1;
  }
  if (*p1) {
    *p0 = *p1;
    if (locked) ReleaseLock(lock);
#pragma omp atomic
    multipole_array->iset -= myrank;
#pragma omp flush
    return n_awgrid;
  }

  orb1 = GetOrbitalSolved(k1);
  orb2 = GetOrbitalSolved(k2);
  
  kappa1 = orb1->kappa;
  kappa2 = orb2->kappa;
  rcl = ReducedCL(GetJFromKappa(kappa1), abs(2*m), 
		  GetJFromKappa(kappa2));
  double *pt = (double *) malloc(sizeof(double)*n_awgrid);
  if (fabs(rcl) < EPS10) {
    for (i = 0; i < n_awgrid; i++) {
      pt[i] = 0;
    }
  } else {
    ef = 0;
    if (orb1->n == 0 || orb2->n == 0) {
      ef = Max(orb1->energy, orb2->energy);
      if (ef > 0.0) {
	ef *= FINE_STRUCTURE_CONST;
      } else {
	ef = 0.0;
      }
    }
    
    npts = potential->maxrp-1;
    if (orb1->n > 0) npts = Min(npts, orb1->ilast);
    if (orb2->n > 0) npts = Min(npts, orb2->ilast);
    r = 0.0;
    jy = 1;
    
    for (i = 0; i < n_awgrid; i++) {
      r = 0.0;
      a = awgrid[i];
      pt[i] = 0.0;
      if (ef > 0.0) a += ef;
      if (m > 0) {
	t = kappa1 + kappa2;
	if (t) {
	  for (j = 0; j <= npts; j++) {
	    x = a*potential->rad[j];
	    n = m;
	    _yk[j] = BESLJN(jy, n, x);
	  }
	  Integrate(_yk, orb1, orb2, 4, &r, 0);
	  r *= t;
	  r *= (2*m + 1.0)/sqrt(m*(m+1.0));
	  r /= pow(a, m);
	  pt[i] = r*rcl;
	}
      } else {
	if (gauge == G_COULOMB) {
	  t = kappa1 - kappa2;
	  q = sqrt(am/(am+1.0));
	  for (j = 0; j <= npts; j++) {
	    x = a*potential->rad[j];
	    n = am+1;
	    _yk[j] = BESLJN(jy, n, x);
	    n = am-1;
	    _zk[j] = BESLJN(jy, n, x);
	  }
	  r = 0.0;
	  rp = 0.0;
	  if (t) {
	    Integrate(_yk, orb1, orb2, 4, &ip, 0);
	    Integrate(_zk, orb1, orb2, 4, &ipm, 0);
	    r = t*ip*q - t*ipm/q;
	  }
	  if (k1 != k2) {
	    Integrate(_yk, orb1, orb2, 5, &im, 0);
	    Integrate(_zk, orb1, orb2, 5, &imm, 0);
	    rp = (am + 1.0)*im*q + am*imm/q;
	  }
	  r += rp;
	  if (am > 1) r /= pow(a, am-1);
	  pt[i] = r*rcl;
	} else if (gauge == G_BABUSHKIN) {
	  t = kappa1 - kappa2;
	  for (j = 0; j <= npts; j++) {
	    x = a*potential->rad[j];
	    n = am+1;
	    _yk[j] = BESLJN(jy, n, x);
	    n = am;
	    _zk[j] = BESLJN(jy, n, x);
	  }
	  if (t) {
	    Integrate(_yk, orb1, orb2, 4, &ip, 0);
	    r = t*ip;
	  }
	  if (k1 != k2) {
	    Integrate(_yk, orb1, orb2, 5, &im, 0);
	  } else {
	    im = 0.0;
	  }
	  Integrate(_zk, orb1, orb2, 1, &imm, 0);
	  rp = (am + 1.0) * (imm + im);
	  q = (2*am + 1.0)/sqrt(am*(am+1.0));
	  q /= pow(a, am);
	  r *= q;
	  rp *= q;
	  pt[i] = (r+rp)*rcl;
	}
      }
    }
  }
#ifdef PERFORM_STATISTICS 
  stop = clock();
  rad_timing.radial_1e += stop - start;
#endif

  *p0 = *p1 = pt;
  if (locked) ReleaseLock(lock);
#pragma omp atomic
    multipole_array->iset -= myrank;
#pragma omp flush
  return n_awgrid;
}

double MultipoleRadialFR(double aw, int m, int k1, int k2, int gauge) {
  int n;
  ORBITAL *orb1, *orb2;
  double *y, ef, r;

  orb1 = GetOrbitalSolved(k1);
  orb2 = GetOrbitalSolved(k2);
  if (orb1->wfun == NULL || orb2->wfun == NULL) {
    if (m == -1) {
      return MultipoleRadialNR(m, k1, k2, gauge);
    } else {
      return 0.0;
    }
  }
  
  n = MultipoleRadialFRGrid(&y, m, k1, k2, gauge);
  if (n == 0) return 0.0;
  ef = 0;
  if (orb1->n == 0 || orb2->n == 0) {
    ef = Max(orb1->energy, orb2->energy);
    if (ef > 0.0) {
      ef *= FINE_STRUCTURE_CONST;
    } else {
      ef = 0.0;
    }
  }
  if (n_awgrid > 1) {
    if (ef > 0) aw += ef;
  }

  r = InterpolateMultipole(aw, n, awgrid, y);
  if (gauge == G_COULOMB && m < 0) r /= aw;

  return r;
}

/* fully relativistic multipole operator, 
   see Grant, J. Phys. B. 1974. Vol. 7, 1458. */ 
double MultipoleRadialFR0(double aw, int m, int k1, int k2, int gauge) {
  double q, ip, ipm, im, imm;
  int kappa1, kappa2;
  int am, t;
  int index[4];
  ORBITAL *orb1, *orb2, *orb;
  double x, a, r, rp, **p1, ef;
  int jy, n, i, j, npts;
  double rcl;

#ifdef PERFORM_STATISTICS 
  clock_t start, stop;
  start = clock();
#endif

  if (m == 0) return 0.0;
  
  if (m >= 0) {
    index[0] = 2*m;
    am = m;
  } else {
    index[0] = -2*m-1;
    am = -m;
  }

  orb1 = GetOrbitalSolved(k1);
  orb2 = GetOrbitalSolved(k2);
  if (orb1->wfun == NULL || orb2->wfun == NULL) {
    if (m == -1) {
      return MultipoleRadialNR(m, k1, k2, gauge);
    } else {
      return 0.0;
    }
  }
  
  index[1] = k1;
  index[2] = k2;
  kappa1 = orb1->kappa;
  kappa2 = orb2->kappa;
  rcl = ReducedCL(GetJFromKappa(kappa1), abs(2*m), 
		  GetJFromKappa(kappa2));
  ef = 0;
  if (orb1->n == 0 || orb2->n == 0) {
    ef = Max(orb1->energy, orb2->energy);
    if (ef > 0.0) {
      ef *= FINE_STRUCTURE_CONST;
    } else {
      ef = 0.0;
    }
  }
  if (n_awgrid > 1) {
    if (ef > 0) aw += ef;
  }
  LOCK *lock = NULL;
  int myrank = MyRankMPI()+1;
  p1 = (double **) MultiSet(multipole_array, index, NULL, &lock,
			    InitPointerData, FreeMultipole);
  int locked = 0;
  if (lock && !(*p1)) {
    SetLock(lock);
    locked = 1;
  }
  if (*p1) {
    r = InterpolateMultipole(aw, n_awgrid, awgrid, *p1);
    if (gauge == G_COULOMB && m < 0) r /= aw;
    r *= rcl;
    if (locked) ReleaseLock(lock);
#pragma omp atomic
    multipole_array->iset -= myrank;
#pragma omp flush
    return r;
  }
  double *pt = (double *) malloc(sizeof(double)*n_awgrid);  
  npts = potential->maxrp-1;
  if (orb1->n > 0) npts = Min(npts, orb1->ilast);
  if (orb2->n > 0) npts = Min(npts, orb2->ilast);
  r = 0.0;
  jy = 1;

  for (i = 0; i < n_awgrid; i++) {
    r = 0.0;
    a = awgrid[i];
    pt[i] = 0.0;
    if (ef > 0.0) a += ef;
    if (m > 0) {
      t = kappa1 + kappa2;
      if (t) {
	for (j = 0; j <= npts; j++) {
	  x = a*potential->rad[j];
	  n = m;
	  _yk[j] = BESLJN(jy, n, x);
	}
	Integrate(_yk, orb1, orb2, 4, &r, 0);
	r *= t;
	r *= (2*m + 1.0)/sqrt(m*(m+1.0));
	r /= pow(a, m);
	pt[i] = r;
      }
    } else {
      if (gauge == G_COULOMB) {
	t = kappa1 - kappa2;
	q = sqrt(am/(am+1.0));
	for (j = 0; j <= npts; j++) {
	  x = a*potential->rad[j];
	  n = am+1;
	  _yk[j] = BESLJN(jy, n, x);
	  n = am-1;
	  _zk[j] = BESLJN(jy, n, x);
	}
	r = 0.0;
	rp = 0.0;
	if (t) {
	  Integrate(_yk, orb1, orb2, 4, &ip, 0);
	  Integrate(_zk, orb1, orb2, 4, &ipm, 0);
	  r = t*ip*q - t*ipm/q;
	}
	if (k1 != k2) {
	  Integrate(_yk, orb1, orb2, 5, &im, 0);
	  Integrate(_zk, orb1, orb2, 5, &imm, 0);
	  rp = (am + 1.0)*im*q + am*imm/q;
	}
	r += rp;
	if (am > 1) r /= pow(a, am-1);
	pt[i] = r;
      } else if (gauge == G_BABUSHKIN) {
	t = kappa1 - kappa2;
	for (j = 0; j < npts; j++) {
	  x = a*potential->rad[j];
	  n = am+1;
	  _yk[j] = BESLJN(jy, n, x);
	  n = am;
	  _zk[j] = BESLJN(jy, n, x);
	}
	if (t) {
	  Integrate(_yk, orb1, orb2, 4, &ip, 0);
	  r = t*ip;
	}
	if (k1 != k2) {
	  Integrate(_yk, orb1, orb2, 5, &im, 0);
	} else {
	  im = 0.0;
	}
	Integrate(_zk, orb1, orb2, 1, &imm, 0);
	rp = (am + 1.0) * (imm + im);
	q = (2*am + 1.0)/sqrt(am*(am+1.0));
	q /= pow(a, am);
	r *= q;
	rp *= q;
	pt[i] = r+rp;
      }
    }
  }

  r = InterpolateMultipole(aw, n_awgrid, awgrid, pt);
  if (gauge == G_COULOMB && m < 0) r /= aw;
  r *= rcl;

#ifdef PERFORM_STATISTICS 
  stop = clock();
  rad_timing.radial_1e += stop - start;
#endif
  *p1 = pt;
  if (locked) ReleaseLock(lock);
#pragma omp atomic
    multipole_array->iset -= myrank;
#pragma omp flush
  return r;
}

double *GeneralizedMoments(int k1, int k2, int m) {
  ORBITAL *orb1, *orb2;
  int n1, i, jy, nk;
  double x, r, r0;
  double *p1, *p2, *q1, *q2;
  int index[3], t;
  double **p, k, *kg;
  double amin, amax, kmin, kmax;
  
  index[0] = m;
  if (k1 > k2) {
    index[1] = k2;
    index[2] = k1;
    orb1 = GetOrbitalSolved(k2);
    orb2 = GetOrbitalSolved(k1);
  } else {
    index[1] = k1;
    index[2] = k2;
    orb1 = GetOrbitalSolved(k1);
    orb2 = GetOrbitalSolved(k2);
  }
  LOCK *lock = NULL;
  int myrank = MyRankMPI()+1;
  p = (double **) MultiSet(gos_array, index, NULL, &lock,
			   InitPointerData, FreeMultipole);
  int locked = 0;
  if (lock && !(*p)) {
    SetLock(lock);
    locked = 1;
  }
  if (*p) {
    if (locked) ReleaseLock(lock);
#pragma omp atomic
    gos_array->iset -= myrank;
#pragma omp flush
    return *p;
  }

  nk = NGOSK;
  double *pt = (double *) malloc(sizeof(double)*nk*2);
  kg = pt + nk;

  if (orb1->wfun == NULL || orb2->wfun == NULL || 
      (orb1->n <= 0 && orb2->n <= 0)) {
    for (t = 0; t < nk*2; t++) {
      pt[t] = 0.0;
    }
    *p = pt;
    if (locked) ReleaseLock(lock);
#pragma omp atomic
    gos_array->iset -= myrank;
#pragma omp flush
    return *p;
  }
  
  jy = 1;
  p1 = Large(orb1);
  p2 = Large(orb2);
  q1 = Small(orb1);
  q2 = Small(orb2);

  amin = sqrt(2.0*fabs(orb1->energy));
  amax = sqrt(2.0*fabs(orb2->energy));
  if (amin < amax) {
    kmin = amin;
    kmax = amax;
  } else {
    kmin = amax;
    kmax = amin;
  }
  kmin = log(0.05*kmin);
  kmax = log(50.0*kmax);
  r = (kmax - kmin)/(nk-1.0);
  kg[0] = kmin;
  for (i = 1; i < nk; i++) {
    kg[i] = kg[i-1] + r;
  }
  
  if (orb1->n > 0 && orb2->n > 0) {
    n1 = Min(orb1->ilast, orb2->ilast);
  
    for (i = 0; i <= n1; i++) {
      _phase[i] = (p1[i]*p2[i] + q1[i]*q2[i])*potential->dr_drho[i];
    }
    
    if (m == 0) {
      if (k1 == k2) r0 = 1.0;
      else if (orb1->n != orb2->n) r0 = 0.0;
      else {
	if (orb1->kappa + orb2->kappa != -1) r0 = 0.0;
	else {
	  r0 = Simpson(_phase, 0, n1);
	}
      }
    } else {
      r0 = 0.0;
    }
    
    for (t = 0; t < nk; t++) {
      k = exp(kg[t]);
      for (i = 0; i <= n1; i++) {
	x = k * potential->rad[i];
	_dphase[i] = BESLJN(jy, m, x);
	_dphase[i] *= _phase[i];
      }
      r = Simpson(_dphase, 0, n1);
      
      pt[t] = (r - r0)/k;
    }
  } else {
    if (orb1->n > 0) n1 = orb1->ilast;
    else n1 = orb2->ilast;
    for (t = 0; t < nk; t++) {
      k = exp(kg[t]);      
      for (i = 0; i <= n1; i++) {
	x = k * potential->rad[i];
	_yk[i] = BESLJN(jy, m, x);
      }
      Integrate(_yk, orb1, orb2, 1, &r, 0);
      pt[t] = r/k;
    }
  }
  *p = pt;
  if (locked) ReleaseLock(lock);
#pragma omp atomic
  gos_array->iset -= myrank;
#pragma omp flush
  return *p;
}

void PrintGeneralizedMoments(char *fn, int m, int n0, int k0, 
			     int n1, int k1, double e1) {
  FILE *f;
  int i0, i1, i;
  double *g, *x;

  i0 = OrbitalIndex(n0, k0, 0.0);
  e1 /= HARTREE_EV;
  i1 = OrbitalIndex(n1, k1, e1);
  f = fopen(fn, "w");
  if (f == NULL) {
    printf("cannot open file %s\n", fn);
    return;
  }
  g = GeneralizedMoments(i0, i1, m);
  x = g + NGOSK;
  for (i = 0; i < NGOSK; i++) {
    fprintf(f, "%15.8E %15.8E %15.8E\n", x[i], exp(x[i]), g[i]);
  }
  fclose(f);
}
  
double InterpolateMultipole(double aw, int n, double *x, double *y) {
  double r;
  int np, nd;

  if (n == 1) {
    r = y[0];
  } else {
    np = 3;
    nd = 1;
    UVIP3P(np, n, x, y, nd, &aw, &r);
  }

  return r;
}

int SlaterTotal(double *sd, double *se, int *j, int *ks, int k, int mode) {
  int t, kk, tt, maxn, minn;
  int tmin, tmax;
  double e, a, d, a1, a2, am, *v1, *v2;
  int err;
  int kl0, kl1, kl2, kl3;
  int k0, k1, k2, k3;
  int js[4];
  ORBITAL *orb0, *orb1, *orb2, *orb3;

#ifdef PERFORM_STATISTICS 
  clock_t start, stop;
  start = clock();
#endif

  k0 = ks[0];
  k1 = ks[1];
  k2 = ks[2];
  k3 = ks[3];
  kk = k/2;

  maxn = 0;
  minn = 1000000;
  orb0 = GetOrbitalSolved(k0);
  orb1 = GetOrbitalSolved(k1);
  orb2 = GetOrbitalSolved(k2);
  orb3 = GetOrbitalSolved(k3);

  if (orb0->n <= 0) {
    maxn = -1;
  } else {
    if (orb0->n > maxn) {
      maxn = orb0->n;
    }
    if (orb0->n < minn) {
      minn = orb0->n;
    }
  }
  if (orb1->n <= 0) {
    maxn = -1;
  } else {
    if (orb1->n > maxn) {
      maxn = orb1->n;
    }
    if (orb1->n < minn) {
      minn = orb1->n;
    }
  }
  if (orb2->n <= 0) {
    maxn = -1;
  } else {
    if (orb2->n > maxn) {
      maxn = orb2->n;
    }
    if (orb2->n < minn) {
      minn = orb2->n;
    }
  }
  if (orb3->n <= 0) {
    maxn = -1;
  } else {
    if (orb3->n > maxn) {
      maxn = orb3->n;
    }
    if (orb3->n < minn) {
      minn = orb3->n;
    }
  }

  if (orb0->wfun == NULL || orb1->wfun == NULL ||
      orb2->wfun == NULL || orb3->wfun == NULL) {
    if (sd) *sd = 0.0;
    if (se) *se = 0.0;
    return 0;
  }

  kl0 = GetLFromKappa(orb0->kappa);
  kl1 = GetLFromKappa(orb1->kappa);
  kl2 = GetLFromKappa(orb2->kappa);
  kl3 = GetLFromKappa(orb3->kappa);

  if (orb1->n < 0 || orb3->n < 0) {
    mode = 2;
  } else {
    if (kl1 > slater_cut.kl0 && kl3 > slater_cut.kl0) {
      if (se) {
	*se = 0.0;
	se = NULL;
      }
    }
    if (kl0 > slater_cut.kl0 && kl2 > slater_cut.kl0) {
      if (se) {
	*se = 0.0;
	se = NULL;
      }
    }  
    if (kl1 > slater_cut.kl1 && kl3 > slater_cut.kl1) {
      mode = 2;
    }
    if (kl0 > slater_cut.kl1 && kl2 > slater_cut.kl1) {
      mode = 2;
    }
  }
  if (qed.br == 0 && IsOdd((kl0+kl1+kl2+kl3)/2)) {
    if (sd) *sd = 0.0;
    if (se) *se = 0.0;
    return 0;
  }

  if (j) {
    memcpy(js, j, sizeof(int)*4);
  } else {
    js[0] = 0;
    js[1] = 0;
    js[2] = 0;
    js[3] = 0;
  }

  if (js[0] <= 0) js[0] = GetJFromKappa(orb0->kappa);
  if (js[1] <= 0) js[1] = GetJFromKappa(orb1->kappa);
  if (js[2] <= 0) js[2] = GetJFromKappa(orb2->kappa);
  if (js[3] <= 0) js[3] = GetJFromKappa(orb3->kappa);  

  am = AMU * potential->atom->mass;
  if (sd) {
    d = 0.0;
    if (Triangle(js[0], js[2], k) && Triangle(js[1], js[3], k)) {
      if (kk == 1 && qed.sms && maxn > 0) {
	v1 = Vinti(k0, k2);
	v2 = Vinti(k1, k3);
      } else {
	v1 = NULL;
	v2 = NULL;
      }
      if (IsEven((kl0+kl2)/2+kk) && IsEven((kl1+kl3)/2+kk)) {	
	err = Slater(&d, k0, k1, k2, k3, kk, mode);
	if (v1 && v2) {
	  d -= (v1[0] + v1[2]) * v2[0]/am;
	}
      }
      a1 = ReducedCL(js[0], k, js[2]);
      a2 = ReducedCL(js[1], k, js[3]); 
      if (v1 && v2) {
	if (IsEven((kl1+kl3)/2+kk)) {
	  d -= v1[1]*v2[0]/am;
	  if (qed.sms == 3 &&
	      fabs(a1)>1e-10 &&
	      fabs(a2)>1e-10 &&
	      kl0 == kl2 && js[0] == js[2]) {
	    d += ((0.25*(kl1==kl3 && js[1]==js[3]))+0.5*a2) *
	      v1[2]*v2[2]/(am*a1*a2);
	  }
	} else {
	  if (qed.sms == 3 &&
	      fabs(a1)>1e-10 &&
	      fabs(a2)>1e-10 &&
	      kl0 == kl2 && js[0] == js[2] &&
	      kl1 == kl3 && js[1] == js[3]) {
	    d += 0.25*v1[2]*v2[2]/(am*a1*a2);
	  }
	}
	if (qed.sms == 3) {
	  d += 0.25*v1[1]*v2[1]/am;
	}
      }
      if (qed.br < 0 || (maxn > 0 && maxn <= qed.br)) {
	if (qed.minbr <= 0 || minn <= qed.minbr) {
	  int mbr = qed.mbr;
	  if (qed.nbr > 0 && (maxn <= 0 || maxn > qed.nbr)) mbr = 0;
	  d += Breit(k0, k1, k2, k3, kk,
		     orb0->kappa, orb1->kappa, orb2->kappa, orb3->kappa,
		     kl0, kl1, kl2, kl3, mbr);
	}
      }
      if (d) {
	d *= a1*a2;      
	if (k0 == k1 && k2 == k3) d *= 0.5;
      }
    }
    *sd = d;
  }
  
  if (!se) goto EXIT;

  if (abs(mode) == 2) {
    *se = 0.0;
    goto EXIT;
  }
  *se = 0.0;
  if (k0 == k1 && (orb0->n > 0 || orb1->n > 0)) goto EXIT;
  if (k2 == k3 && (orb2->n > 0 || orb3->n > 0)) goto EXIT;
  tmin = abs(js[0] - js[3]);
  tt = abs(js[1] - js[2]);
  tmin = Max(tt, tmin);
  tmax = js[0] + js[3];
  tt = js[1] + js[2];
  tmax = Min(tt, tmax);
  tmax = Min(tmax, GetMaxRank());
  if (IsOdd(tmin)) tmin++;
  
  for (t = tmin; t <= tmax; t += 2) {
    a = W6j(js[0], js[2], k, js[1], js[3], t);
    if (fabs(a) > EPS30) {
      e = 0.0;
      if (t == 2 && qed.sms && maxn > 0) {
	v1 = Vinti(k0, k3);
	v2 = Vinti(k1, k2);
      } else {
	v1 = NULL;
	v2 = NULL;
      }
      if (IsEven((kl0+kl3+t)/2) && IsEven((kl1+kl2+t)/2)) {
	err = Slater(&e, k0, k1, k3, k2, t/2, mode);
	if (v1 && v2) {
	  e -= (v1[0]+v1[2])*v2[0]/am;
	}
      }
      a1 = ReducedCL(js[0], t, js[3]); 
      a2 = ReducedCL(js[1], t, js[2]);
      if (v1 && v2) {
	if (IsEven((kl1+kl2+t)/2)) {
	  e -= v1[2]*v2[0]/am;
	  if (qed.sms == 3 &&
	      fabs(a1)>1e-10 &&
	      fabs(a2)>1e-10 &&
	      kl0 == kl3 && js[0] == js[3]) {
	    e += ((0.25*(kl1==kl2&&js[1]==js[2]))+0.5*a2) *
	      v1[2]*v2[2]/(am*a1*a2);
	  }
	} else {
	  if (qed.sms == 3 &&
	      fabs(a1)>1e-10 &&
	      fabs(a2)>1e-10 &&
	      kl0 == kl3 && js[0] == js[3] &&
	      kl1 == kl2 && js[1] == js[2]) {
	    e += 0.25*v1[2]*v2[2]/(am*a1*a2);
	  }
	}
	if (qed.sms == 3) {
	  e += 0.25*v1[1]*v2[1]/am;
	}
      }
      if (qed.br < 0 || (maxn > 0 && maxn <= qed.br)) {
	if (qed.minbr <= 0 || minn <= qed.minbr) {
	  int mbr = qed.mbr;
	  if (qed.nbr > 0 && (maxn <= 0 || maxn > qed.nbr)) mbr = 0;
	  e += Breit(k0, k1, k3, k2, t/2,
		     orb0->kappa, orb1->kappa, orb3->kappa, orb2->kappa,
		     kl0, kl1, kl3, kl2, mbr);
	}
      }
      if (e) {
	e *= a * (k + 1.0) * a1 * a2;
	if (IsOdd(t/2 + kk)) e = -e;
	*se += e;
      }
    }
  }

 EXIT:
#ifdef PERFORM_STATISTICS 
    stop = clock();
    rad_timing.radial_slater += stop - start;
#endif
  return 0;
}

double SelfEnergyRatio(ORBITAL *orb, ORBITAL *horb) {
  int i, k, m, npts;
  double *p, *q, e, z;
  double *large, *small;
  double a, b;
  
  if (orb->wfun == NULL) return 1.0;
  m = qed.mse%10;
  if (m >= 2) {
    //modqed
    k = IdxVT(orb->kappa)-1;
    z = potential->atom->atomic_number;
    if (z < 10 || z > 120 || k < 0 || k >= NKSEP) {
      m = 1;
    } else {
      b = HydrogenicSelfEnergy(51, qed.pse, 1.0, potential, orb, orb);
      a = HydrogenicSelfEnergy(51, qed.pse, 1.0, potential, horb, horb);
      return b/a;
    }
  }
  npts = potential->maxrp;
  p = Large(horb);
  q = Small(horb);
  large = Large(orb);
  small = Small(orb);
  for (i = 0; i < npts; i++) {
    if (m == 0) {
      //uehling potential
      a = potential->ZVP[i]/potential->rad[i];
    } else if (m == 1) {
      //welton scaling
      a = potential->qdist[i];
    } else {
      a = 1.0;
    }
    a *= potential->dr_drho[i];
    _dwork11[i] = (p[i]*p[i] + q[i]*q[i])*a;
    _dwork12[i] = (large[i]*large[i] + small[i]*small[i])*a;
  }
  if (i < 5) return 1.0;
  a = Simpson(_dwork11, 0, i-1);
  b = Simpson(_dwork12, 0, i-1);
  if (1+a == 1) return 1.0;
  return b/a;
}

ORBITAL *SolveAltOrbital(ORBITAL *orb, POTENTIAL *p) {
  ORBITAL *norb;
  int ierr;

  if (!p->hlike && !p->hpvs) {
    return orb;
  }
  norb = (ORBITAL *) malloc(sizeof(ORBITAL));
  InitOrbitalData(norb, 1);
  norb->n = orb->n;
  norb->kappa = orb->kappa;
  p->flag = -1;
  ierr = RadialSolver(norb, p);
  if (ierr < 0) {
    MPrintf(-1, "error in SolveAltOrbital: %d %d %d\n",
	    ierr, norb->n, norb->kappa);
  }
  return norb;
}

double SelfEnergy(ORBITAL *orb1, ORBITAL *orb2) {
  double a, c;
  double an = 0.0;
  ORBITAL *orb0 = NULL;
  int msc = qed.mse%10;
  int ksc = qed.mse/10;
  
  if (qed.se == -1000000) return 0.0;
  if (orb1->n <= 0 || orb2->n <= 0) return 0.0;
  int nm = Min(orb1->n, orb2->n);
  if (orb1->energy > 0 && orb2->energy > 0) {
    return 0.0;
  } else {
    double ae, ose, eb0;
    int idx, nb0;
    int kv = IdxVT(orb1->kappa);
    if (kv <= 0) return 0;
    int kv1 = kv-1;
    if (potential->ib > 0 && potential->nb > 0) {
      eb0 = fabs(Min(orb1->energy, orb2->energy));
      ae = qed.ose1;
      eb0 *= ae;
      ose = qed.ose0;
      if (qed.ase) {
	ose /= pow(nm, qed.ase);
      }
      nb0 = potential->nb+1;
      int nk = 1+GetLFromKappa(orb1->kappa)/2;
      if (nb0 < nk) nb0 = nk;
      for (;;nb0++) {
	if (potential->mpseudo > 0 && nb0 > potential->mpseudo) break;
	idx = OrbitalIndex(nb0, orb1->kappa, 0);
	orb0 = GetOrbitalSolved(idx);
	if (orb0->energy > eb0) break;
	if (orb0->rfn < potential->rfn[kv1]) break;
      }
      potential->nfn[kv1] = nb0;
      orb0 = NULL;
      if (orb1->n < nb0 && orb2->n >= nb0) {
	if (ose >= 0) {
	  an = eb0/orb2->energy;
	  an = Min(1.0, an);
	  an = ose>0?pow(an, ose):1.0;
	  if (orb2->n > nb0) {
	    int idx = OrbitalIndex(nb0, orb2->kappa, 0);
	    orb2 = GetOrbitalSolved(idx);
	  }
	  orb0 = orb1;
	} else {
	  if (orb2->rfn < potential->rfn[kv1] &&
	      orb2->energy > eb0) return 0.0;
	}
      } else if (orb1->n >= nb0 && orb2->n < nb0) {
	if (ose >= 0) {
	  an = eb0/orb1->energy;
	  an = Min(1.0, an);
	  an = ose>0?pow(an, ose):1.0;
	  if (orb1->n > nb0) {
	    int idx = OrbitalIndex(nb0, orb1->kappa, 0);
	    orb1 = GetOrbitalSolved(idx);
	  }
	  orb0 = orb2;
	} else {
	  if (orb1->rfn < potential->rfn[kv1] &&
	      orb1->energy > eb0) return 0.0;
	}
      }
    }
  }
  if (potential->ib > 0 &&
      (orb1->n > potential->nb || orb2->n > potential->nb)) {
    if (ksc < 4) return 0;
  }
  if (orb1 != orb2) {
    if (ksc < 6) return 0.0;
    if (orb1->n <= 0 || orb2->n <= 0) return 0.0;
    if (qed.se >= 0 && (orb1->n > qed.se || orb2->n > qed.se)) return 0.0;
    if (potential->nb > nm &&
	(orb1->n > potential->nb || orb2->n > potential->nb)) {
      c = potential->nb-nm;
      c = pow(c, qed.ise);
      c = 1 + qed.cse0*c/(1.0 + qed.cse1*c);
    } else {
      c = 0.0;
    }
    if (orb0 != NULL && orb0->ose < 1e30) {
      a = orb0->ose;
      if (an > 0) a *= an;
      if (c > 0) a *= c;
      return a;
    }
    if (orb1->rorb == NULL) {
      orb1->rorb = SolveAltOrbital(orb1, rpotential);
    }
    if (orb2->rorb == NULL) {
      orb2->rorb = SolveAltOrbital(orb2, rpotential);
    }
    a = HydrogenicSelfEnergy(qed.mse, qed.pse, c, potential,
			     orb1->rorb, orb2->rorb);
    if (orb0 != NULL) orb0->ose = a;
    if (an > 0) a *= an;
    if (c > 0) a *= c;
    return a;
  }  
  if (orb1->se < 0.999e31) return orb1->se;
  if (orb1->n <= 0) {
    orb1->se = 0.0;
    return 0.0;
  }
  if (!(qed.se < 0 || orb1->n <= qed.se)) {
    orb1->se = 0.0;
    return 0.0;
  }
  if (qed.sse > 0 && orb1->n > qed.sse) {
    int k = GetLFromKappa(orb1->kappa)/2;
    if (k < orb1->n-1) {
      int idx = OrbitalIndex(k+1, orb1->kappa, 0);
      ORBITAL *orb = GetOrbitalSolved(idx);
      a = SelfEnergy(orb, orb);
      c = SelfEnergyRatio(orb1, orb);
      orb1->se = a*c;
      if (qed.pse) {
	MPrintf(-1, "SE: z=%g, n=%d, kappa=%2d, e0=%11.4E, md=%d, screen=%11.4E, final=%11.4E\n", potential->Z[potential->maxrp-1], orb1->n, orb1->kappa, orb1->energy, qed.mse, c, orb1->se);
      }
      return orb1->se;
    }
  }
  if (orb1->rorb == NULL) {
    orb1->rorb = SolveAltOrbital(orb1, rpotential);
  }
  if (msc == 9 || fabs(potential->N-1)<1e-5) {
    c = 1.0;
  } else {
    if (orb1->horb == NULL) {
      orb1->horb = SolveAltOrbital(orb1, hpotential);
    }
    if (orb1->horb->wfun == NULL) {
      MPrintf(-1, "hlike orbital for SE screening null wfun: %d %d\n",
	      orb1->horb->n, orb1->horb->kappa);
      c = 1.0;
    }
    c = SelfEnergyRatio(orb1->rorb, orb1->horb);
  }
  orb1->se = HydrogenicSelfEnergy(qed.nms*100+qed.mse, qed.pse, c, potential, orb1->rorb, NULL);
  return orb1->se;
}

double RadialNMS(ORBITAL *orb1, ORBITAL *orb2, int kv) {
  int i, k0, k1;
  double a, r, r2, mass;
  
  if (qed.se == -1000000) return 0.0;  
  if (qed.nms <= 0) return 0.0;
  if (orb1->n <= 0 || orb2->n <= 0) return 0.0;
  
  mass = AMU*potential->atom->mass;
  r = 0.0;
  
  if (qed.nms == 1) {
    k0 = orb1->idx;
    k1 = orb2->idx;
    for (i = 0; i < potential->maxrp; i++) {
      _yk[i] = potential->VT[kv][i];
    }
    a = 0.0;
    Integrate(_yk, orb1, orb2, 1, &a, -1);
    a = -a;
    if (k0 == k1) a += orb1->energy;
    a /= mass;
    r = a;
    for (i = 0; i < potential->maxrp; i++) {
      _yk[i] = orb1->energy - potential->VT[kv][i];
      _yk[i] *= orb2->energy - potential->VT[kv][i];
    }
    Integrate(_yk, orb1, orb2, 1, &a, -1);
    a *= FINE_STRUCTURE_CONST2/(2.0 * mass);
    r += a;
  } else {
    double *p1, *p2, *q1, *q2, az;
    int m1, j2, k, kt;
    m1 = Min(orb1->ilast, orb2->ilast);
    p1 = Large(orb1);
    p2 = Large(orb2);
    q1 = Small(orb1);
    q2 = Small(orb2);
    GetJLFromKappa(orb1->kappa, &j2, &k);
    k /= 2;
    kt = j2-k;
    DrLargeSmall(orb1, potential, _dwork1, _dwork3);
    DrLargeSmall(orb2, potential, _dwork2, _dwork4);
    for (i = 0; i <= m1; i++) {
      _dwork[i] = _dwork1[i]*_dwork2[i];
      _dwork[i] += _dwork3[i]*_dwork4[i];
      r2 = potential->rad[i];
      r2 *= r2;
      a = p1[i]*p2[i]*k*(k+1.0) + q1[i]*q2[i]*kt*(kt+1.0);
      _dwork[i] += a/r2;
      if (qed.nms > 2) {
	az = FINE_STRUCTURE_CONST*potential->VT[kv][i];
	a = 2*az*(q1[i]*_dwork2[i]+q2[i]*_dwork1[i]);
	_dwork[i] += a;
	a = az*(orb1->kappa-1)*(q1[i]*p2[i]+q2[i]*p1[i]);
	_dwork[i] += a/potential->rad[i];
      }
      _dwork[i] *= potential->dr_drho[i];
    }
    r = Simpson(_dwork, 0, m1);
    r /= 2*mass;
  }  
  if (qed.pse > 1) {
    MPrintf(-1, "NMS: %d %d %d %d %d %18.10E %18.10E\n", qed.nms, orb1->n, orb1->kappa, orb2->n, orb2->kappa, mass, r);
  }
  return r;
}

double QED1E(int k0, int k1) {
  int i, kv;
  ORBITAL *orb1, *orb2;
  int index[2];
  double *p, r, a, c;

  if (qed.se == -1000000) return 0.0;  
  if (qed.nms == 0 && qed.se == 0 && qed.vp%100 == 0) {
      return 0.0;
  }

  orb1 = GetOrbitalSolved(k0);
  orb2 = GetOrbitalSolved(k1);
  
  kv = orb1->kv;
  if (kv < 0 || kv > NKSEP) {
    MPrintf(-1, "invalid orbital kv in RadialNMS: %d %d %d\n",
	    orb1->n, orb1->kappa, orb1->kv);
    return 0.0;
  }
  
  if (orb1->n <= 0 || orb2->n <= 0) {
    return 0.0;
  }
  if (orb1->wfun == NULL || orb2->wfun == NULL) {
    return 0.0;
  }
  
  if (k0 > k1) {
    index[0] = k1;
    index[1] = k0;
  } else {
    index[0] = k0;
    index[1] = k1;
  }
  
  LOCK *lock = NULL;
  int myrank = MyRankMPI()+1;
  p = (double *) MultiSet(qed1e_array, index, NULL, &lock,
			  InitDoubleData, NULL);
  int locked = 0;
  if (lock && !(p && *p)) {
    SetLock(lock);
    locked = 1;
  }
  if (p && *p) {
    if (locked) ReleaseLock(lock);
#pragma omp atomic
    qed1e_array->iset -= myrank;
#pragma omp flush
    return *p;
  }
  r = 0.0;
  if (orb1->kappa != orb2->kappa) {
    printf("dk: %d %d %d %d\n", orb1->n, orb1->kappa, orb2->n, orb2->kappa);
  }
  a = RadialNMS(orb1, orb2, kv);
  r += a;
  if (!potential->pvp && potential->mvp) {
    for (i = 0; i < potential->maxrp; i++) {
      _yk[i] = -potential->ZVP[i]/potential->rad[i];
    }
    Integrate(_yk, orb1, orb2, 1, &a, 0);
    if (qed.pse > 2) {
      MPrintf(-1, "VP: %d %d %d %d %18.10E\n", orb1->n, orb1->kappa, orb2->n, orb2->kappa, a);
    }
    r += a;      
  }
  if (optimize_control.mce < 20) {
    a = SelfEnergy(orb1, orb2);
    r += a;
    if (k0 == k1) {
      orb1->qed = r;
    }
  }
  *p = r;
  if (locked) ReleaseLock(lock);
#pragma omp atomic
  qed1e_array->iset -= myrank;
#pragma omp flush
  return r;
}

double *Vinti(int k0, int k1) {
  int i;
  ORBITAL *orb1, *orb2;
  int index[2];
  double **p, *large0, *small0, *large1, *small1;
  int ka0, ka1, m1;
  double a, b;

  if (qed.se == -1000000) {    
    return NULL;
  }

  orb1 = GetOrbitalSolved(k0);
  orb2 = GetOrbitalSolved(k1);

  if (orb1->n <= 0 || orb2->n <= 0) {
    return NULL;
  }
  if (orb1->wfun == NULL || orb2->wfun == NULL) {
    return NULL;
  }
  
  index[0] = k0;
  index[1] = k1;
  LOCK *lock = NULL;
  int myrank = MyRankMPI()+1;
  p = (double **) MultiSet(vinti_array, index, NULL, &lock,
			   InitPointerData, FreeMultipole);
  int locked = 0;
  if (lock && !(*p)) {
    SetLock(lock);
    locked = 1;
  }
  if (*p) {
    if (locked) ReleaseLock(lock);
#pragma omp atomic
    vinti_array->iset -= myrank;
#pragma omp flush
    return *p;
  }

  double *r;
  r = (double *) malloc(sizeof(double)*3);
  m1 = Min(orb1->ilast, orb2->ilast);
  ka0 = orb1->kappa;
  ka1 = orb2->kappa;
  large0 = Large(orb1);
  large1 = Large(orb2);
  small0 = Small(orb1);
  small1 = Small(orb2);
  a = 0.5*(ka0*(ka0+1.0) - ka1*(ka1+1.0));
  b = 0.5*(-ka0*(-ka0+1.0) + ka1*(-ka1+1.0));
  DrLargeSmall(orb2, potential, _dwork1, _dwork2);  
  for (i = 0; i <= m1; i++) {
    _yk[i] = large0[i]*_dwork1[i] - a*large0[i]*large1[i]/potential->rad[i];
    _yk[i] += small0[i]*_dwork2[i] - b*small0[i]*small1[i]/potential->rad[i];
    _yk[i] *= potential->dr_drho[i];
  }
  r[0] = Simpson(_yk, 0, m1);
  r[1] = r[2] = 0;  
  if (qed.sms == 1) {
    *p = r;
    if (locked) ReleaseLock(lock);
#pragma omp atomic
    vinti_array->iset -= myrank;
#pragma omp flush
    return r;
  }

  int j0, j1, kl0, kl1, kt0, kt1, kv;
  double sd, sd1, sd2, az, ac;
  kv = orb1->kv;
  if (kv < 0 || kv > NKSEP) {
    MPrintf(-1, "invalid orbital kv in Vinti: %d %d %d\n",
	    orb1->n, orb1->kappa, orb1->kv);
    kv = 0;
  }
  GetJLFromKappa(orb1->kappa, &j0, &kl0);
  GetJLFromKappa(orb2->kappa, &j1, &kl1);
  kt0 = 2*j0-kl0;
  kt1 = 2*j1-kl1;
  r[1] = 0;
  if (kt0 == kl1 || kl0 == kt1) {
    ac = ReducedCL(j0, 2, j1);
    if (fabs(ac) > 1e-10) {
      for (i = 0; i < potential->maxrp; i++) {
	_yk[i] = 0;
      }
      sd = sqrt(6*(j0+1.0)*(j1+1.0));
      if (kt0 == kl1) {
	sd1 = sd*W6j(1, j0, kt0, j1, 1, 2);
	if (IsEven((kt0+1+j0)/2)) sd1 = -sd1;
	for (i = 0; i < potential->maxrp; i++) {
	  _yk[i] -= small0[i]*large1[i]*sd1;
	}
      }
      if (kl0 == kt1) {
	sd2 = sd*W6j(1, j0, kl0, j1, 1, 2);
	if (IsEven((kl0+1+j0)/2)) sd2 = -sd2;
	for (i = 0; i <= m1; i++) {
	  _yk[i] += small1[i]*large0[i]*sd2;
	}
      }    
      for (i = 0; i <= m1; i++) {
	az = FINE_STRUCTURE_CONST*potential->VT[kv][i];
	_yk[i] *= az*potential->dr_drho[i];
      }
      r[1] = -Simpson(_yk, 0, m1);
      r[1] /= ac;
    }
  }
  
  for (i = 0; i <= m1; i++) {
    _yk[i] = small0[i]*large1[i] - small1[i]*large0[i];
    az = FINE_STRUCTURE_CONST*potential->VT[kv][i];
    _yk[i] *= az*potential->dr_drho[i];
  }
  r[2] = -Simpson(_yk, 0, m1);
  *p = r;
  if(locked) ReleaseLock(lock);
#pragma omp atomic
  vinti_array->iset -= myrank;
#pragma omp flush
  return *p;
}

double BreitC(int n, int m, int k, int k0, int k1, int k2, int k3) {
  int ka0, ka1, ka2, ka3, kb, kp;
  double r, b, c;
  
  ka0 = GetOrbital(k0)->kappa;
  ka1 = GetOrbital(k1)->kappa;
  ka2 = GetOrbital(k2)->kappa;
  ka3 = GetOrbital(k3)->kappa;
  if (k == m) {
    r = -(ka0 + ka2)*(ka1 + ka3);
    if (r) r /= (m*(m+1.0));
  } else if (k == (m + 1)) {
    kb = ka2 - ka0;
    kp = ka3 - ka1;
    b = (m + 2.0)/(2.0*(2.0*m + 1.0));
    c = -(m - 1.0)/((2.0*m+1.0)*(2.0*m+2.0));
    switch (n) {
    case 0:
      r = (k + kb)*(b + c*kp);
      break;
    case 1:
      r = (k + kp)*(b + c*kb);
      break;
    case 2:
      r = (k - kb)*(b - c*kp);
      break;
    case 3:
      r = (k - kp)*(b - c*kb);
      break;
    case 4:
      r = -(k + kb)*(b - c*kp);
      break;
    case 5:
      r = -(k - kp)*(b + c*kb);
      break;
    case 6:
      r = -(k - kb)*(b + c*kp);
      break;
    case 7:
      r = -(k + kp)*(b - c*kb);
      break;
    default:
      r = 0;
    }
  } else {
    kb = ka2 - ka0;
    kp = ka3 - ka1;
    b = (m - 1.0)/(2.0*(2.0*m + 1.0));
    c = (m + 2.0)/(2.0*m*(2.0*m + 1.0));
    switch (n) {
    case 0:
      r = (b + c*kb)*(kp - k - 1.0);
      break;
    case 1:
      r = (b + c*kp)*(kb - k - 1.0);
      break;
    case 2:
      r = (b - c*kb)*(-kp - k - 1.0);
      break;
    case 3:
      r = (b - c*kp)*(-kb - k - 1.0);
      break;
    case 4:
      r = -(b + c*kb)*(-kp - k - 1.0);
      break;
    case 5:
      r = -(b - c*kp)*(kb - k - 1.0);
      break;
    case 6:
      r = -(b - c*kb)*(kp - k - 1.0);
      break;
    case 7:
      r = -(b + c*kp)*(-kb - k - 1.0);
      break;
    default:
      r = 0;
    }
  }

  return r;
}

double BreitS(int k0, int k1, int k2, int k3, int k) {
  ORBITAL *orb0, *orb1, *orb2, *orb3;
  int index[5], i;
  double *p0, *z, r;
  FLTARY *byk;
  LOCK *lock = NULL;
  int locked = 0;
  int myrank = MyRankMPI()+1;
  if (breit_array->maxsize != 0) {
    index[0] = k0;
    index[1] = k1;
    index[2] = k2;
    index[3] = k3;
    index[4] = k;
    p0 = (double *) MultiSet(breit_array, index, NULL, &lock,
			     InitDoubleData, NULL);
    if (lock && !(p0 && *p0)) {
      SetLock(lock);
      locked = 1;
    }
    if (p0 && *p0) {
      r = *p0;
      if (locked) ReleaseLock(lock);
#pragma omp atomic
      breit_array->iset -= myrank;
#pragma omp flush
      return r;
    }
  }
  byk = NULL;
  z = _zk;
  LOCK *xlock = NULL;
  if (xbreit_array[4]->maxsize != 0) {
    index[0] = k0;
    index[1] = k1;
    index[2] = k;
    byk = (FLTARY *) MultiSet(xbreit_array[4], index, NULL, &xlock,
			      InitFltAryData, FreeFltAryData);
    if (xlock) SetLock(xlock);
    if (byk->npts > 0) {
      for (i = 0; i < byk->npts; i++) {
	z[i] = byk->yk[i];
      }
    }
  }
  int npts;
  if (byk == NULL || byk->npts < 0) {
    orb0 = GetOrbitalSolved(k0);
    orb1 = GetOrbitalSolved(k1);
    for (i = 0; i < potential->maxrp; i++) {
      _dwork1[i] = pow(potential->rad[i], k);
    }    
    Integrate(_dwork1, orb0, orb1, -6, z, 0);    
    for (i = 0; i < potential->maxrp; i++) {
      if (z[i]) z[i] /= _dwork1[i]*potential->rad[i];
    }
    for (i = potential->maxrp-1; i >= 0; i--) {
      if (z[i]) break;
    }
    npts = i+1;
    if (byk != NULL) {
      int size = sizeof(float)*npts;
      byk->yk = malloc(size);
      AddMultiSize(xbreit_array[4], size);
      for (i = 0; i < npts; i++) {
	byk->yk[i] = z[i];
      }
      byk->npts = npts;
    }
  }
  if (xlock) ReleaseLock(xlock);
  orb2 = GetOrbitalSolved(k2);
  orb3 = GetOrbitalSolved(k3);
  Integrate(z, orb2, orb3, 6, &r, 0);
  if (breit_array->maxsize != 0) {
    if (!r) r = 1e-100;
    *p0 = r;
    if (locked) ReleaseLock(lock);
#pragma omp atomic
    breit_array->iset -= myrank;
  }
  if (xbreit_array[4]->maxsize != 0) {
#pragma omp atomic
    xbreit_array[4]->iset -= myrank;
  }
#pragma omp flush
  return r;
}

double BreitI(int n, int k0, int k1, int k2, int k3, int m) {
  double r;

  switch (n) {
  case 0:
    r = BreitS(k0, k2, k1, k3, m);
    break;
  case 1:
    r = BreitS(k1, k3, k0, k2, m);
    break;
  case 2:
    r = BreitS(k2, k0, k3, k1, m);
    break;
  case 3:
    r = BreitS(k3, k1, k2, k0, m);
    break;
  case 4:
    r = BreitS(k0, k2, k3, k1, m);
    break;
  case 5:
    r = BreitS(k3, k1, k0, k2, m);
    break;
  case 6:
    r = BreitS(k2, k0, k1, k3, m);
    break;
  case 7:
    r = BreitS(k1, k3, k2, k0, m);
    break;
  default:
    r = 0.0;
  }

  return r;
}

int BreitX(ORBITAL *orb0, ORBITAL *orb1, int k, int m, int w, int mbr,
	   double e, double *y) {
  int i;
  double kf = 1.0;
  double x, r, b;
  int k2 = 2*k;
  int jy, k1, fd, rs;
  int index[3];  
  FLTARY *byk;

  if (y == NULL) y = _xk;
  if (e < 0) e = fabs(orb0->energy-orb1->energy);
  byk = NULL;
  LOCK *lock = NULL;
  int locked = 0;
  int myrank = MyRankMPI()+1;
  if ((mbr == 2 || (m < 3 && w == 0) || (m == 3 && w == 1))
      && xbreit_array[m]->maxsize != 0) {
    index[0] = orb0->idx;
    index[1] = orb1->idx;
    index[2] = k;
    byk = (FLTARY *) MultiSet(xbreit_array[m], index, NULL, &lock,
			      InitFltAryData, FreeFltAryData);
    if (lock && byk->npts <= 0) {
      SetLock(lock);
      locked = 1;
    }
    if (byk->npts > 0) {
      for (i = 0; i < byk->npts; i++) {
	if (e > 0) {
	  _dwork1[i] = FINE_STRUCTURE_CONST*e*potential->rad[i];
	} else {
	  _dwork1[i] = 0;
	}
	_dwork2[i] = pow(potential->rad[i], k);
	y[i] = byk->yk[i];
      }
      if (locked) ReleaseLock(lock);
#pragma omp atomic
      xbreit_array[m]->iset -= myrank;
#pragma omp flush
      return byk->npts;
    }     
  }

  for (i = 1; i < k2; i += 2) {
    kf *= i;
  }
  double ef = FINE_STRUCTURE_CONST*e;
  double efk = pow(ef, k);

  for (i = 0; i < potential->maxrp; i++) {    
    if (e > 0) {
      x = FINE_STRUCTURE_CONST*e*potential->rad[i];
      _dwork1[i] = x;
    } else {
      x = 0;
      _dwork1[i] = 0;
    }
    _dwork2[i] = pow(potential->rad[i], k);
    switch (m) {
    case 0:
      if (x < qed.xbr) {
	_dwork[i] = (1 - 0.5*x*x/(k2+3.0))*_dwork2[i];
      } else {
	jy = 1;
	_dwork[i] = BESLJN(jy, k, x);
	_dwork[i] *= kf*(k2+1.0)/efk;	 
      }
      break;
    case 1:
      _dwork[i] = _dwork2[i]/potential->rad[i];
      break;
    case 2:
      if (x < qed.xbr) {
	_dwork[i] = 0.5/(1.0+k2);
      } else {
	jy = 3;
	k1 = k - 1;
	_dwork[i] = BESLJN(jy, k1, x);	
      }
      b = _dwork2[i]*potential->rad[i];
      _dwork[i] *= b;
      break;
    case 3:
      if (x < qed.xbr) {
	_dwork[i] = (1 - 0.5*x*x/(5+k2))*_dwork2[i]*potential->rad[i];
      } else {
	jy = 1;
	k1 = k + 1;
	_dwork[i] = BESLJN(jy, k1, x);
	_dwork[i] *= kf*(k2+1.0)*(k2+3.0)/(efk*ef);
      }
      break;
    default:
      _dwork[i] = 0;
      break;
    }
  }

  Integrate(_dwork, orb0, orb1, -6, y, 0);
  switch (m) {
  case 0:
    for (i = 0; i < potential->maxrp; i++) {
      y[i] /= _dwork2[i]*potential->rad[i];
    }
    break;
  case 1:
    for (i = 0; i < potential->maxrp; i++) {
      y[i] /= _dwork2[i];
    }
    break;
  case 2:
  case 3:
    for (i = 0; i < potential->maxrp; i++) {
      y[i] /= _dwork2[i]*potential->rad[i]*potential->rad[i];
    }
    break;
  }
  int npts;
  for (i = potential->maxrp-1; i >= 0; i--) {
    if (y[i]) break;
  }
  npts = i+1;
  if (byk) {
    int size = sizeof(float)*npts;
    byk->yk = malloc(size);
    AddMultiSize(xbreit_array[m], size);
    for (i = 0; i < npts; i++) {
      byk->yk[i] = y[i];
    }
    byk->npts = npts;
#pragma omp atomic
    xbreit_array[m]->iset -= myrank;
  }
  if (locked) ReleaseLock(lock);
#pragma omp flush
  return npts;
}

double BreitRW(int k0, int k1, int k2, int k3, int k, int w, int mbr) {
  double e, x, r, b, kf;
  ORBITAL *orb0, *orb1, *orb2, *orb3;
  int maxrp, i, jy, kd;
  
  orb0 = GetOrbital(k0);
  orb1 = GetOrbital(k1);
  orb2 = GetOrbital(k2);
  orb3 = GetOrbital(k3);
  if (mbr == 2) {
    e = 0;
  } else {
    if (w == 0) {
      e = fabs(orb0->energy - orb2->energy);
    } else {
      e = fabs(orb1->energy - orb3->energy);
    }
  }
  double ef = FINE_STRUCTURE_CONST*e;
  double efk = pow(ef, k);
  kf = 1;
  kd = 2*k;
  for (i = 1; i < kd; i += 2) {
    kf *= i;
  }
  int npts = BreitX(orb0, orb2, k, 0, w, mbr, e, _xk);
  for (i = 0; i < npts; i++) {
    x = _dwork1[i];
    if (x < qed.xbr) {
      b = 1 - 0.5*x*x/(1.0-kd);
    } else {
      jy = 2;
      b = -BESLJN(jy, k, x)*(_dwork2[i]*potential->rad[i]*efk*ef)/kf;
    }
    _xk[i] *= b;
  }
  for (; i < potential->maxrp; i++) _xk[i] = 0.0;
  Integrate(_xk, orb1, orb3, 6, &r, 0);
  return r;
}

double BreitRK(int k0, int k1, int k2, int k3, int k, int mbr) {
  double r1, r2, r3, r4;
  r1 = BreitRW(k0, k1, k2, k3, k, 0, mbr);
  if ((k0 == k1 && k2 == k3) || (k0 == k3 && k2 == k1)) {
    r2 = r1;
  } else {
    r2 = BreitRW(k0, k1, k2, k3, k, 1, mbr);
  }
  if (k0 == k1 && k2 == k3) {
    r3 = r1;
  } else {
    r3 = BreitRW(k1, k0, k3, k2, k, 0, mbr);
  }
  if ((k0 == k1 && k2 == k3) || (k0 == k3 && k2 == k1)) {
    r4 = r3;
  } else {
    r4 = BreitRW(k1, k0, k3, k2, k, 1, mbr);
  }
  return 0.5*(r1+r2+r3+r4);
}

double BreitSW(int k0, int k1, int k2, int k3, int k, int w, int mbr) {  
  double e, x, xk, r, b, kf, s1, s2;
  ORBITAL *orb0, *orb1, *orb2, *orb3;
  int i, jy, kd, kj;
  
  orb0 = GetOrbital(k0);
  orb1 = GetOrbital(k1);
  orb2 = GetOrbital(k2);
  orb3 = GetOrbital(k3);
  if (mbr == 2) {
    e = 0.0;
  } else {
    if (w == 0) {
      e = fabs(orb0->energy - orb2->energy);
    } else {
      e = fabs(orb1->energy - orb3->energy);
    }
  }
  double ef = FINE_STRUCTURE_CONST*e;
  double efk = pow(ef, k);
  kf = 1;
  kd = 2*k;
  for (i = 1; i < kd; i += 2) {
    kf *= i;
  }

  int npts1 = BreitX(orb0, orb2, k, 1, w, mbr, e, _dwork12);
  int npts2 = BreitX(orb0, orb2, k, 2, w, mbr, e, _dwork13);
  int npts = Min(npts1, npts2);
  for (i = 0; i < npts; i++) {
    x = _dwork1[i];
    xk = _dwork2[i]*efk;
    if (x < qed.xbr) {
      b = -0.5/(1+kd);
    } else {
      jy = 4;
      kj = k + 1;
      b = BESLJN(jy, kj, x);
    }
    _dwork12[i] = _dwork12[i]*b;
    _dwork13[i] = _dwork13[i]*(1-x*x*b);
    _yk[i] = _dwork12[i] + _dwork13[i];
    b = (kd+1);
    _yk[i] *= b*b;
  }
  for (; i < potential->maxrp; i++) _yk[i] = 0.0;
  Integrate(_yk, orb1, orb3, 6, &s1, 0);
  if (k >= 0) {
    npts = BreitX(orb1, orb3, k, 3, w, mbr, e, _dwork12);
    for (i = 0; i < npts; i++) {
      x = _dwork1[i];
      if (x < qed.xbr) {
	b = 1 - 0.5*x*x/(3.0-kd);
      } else {
	jy = 2;
	kj = k-1;
	xk = _dwork2[i]*efk;
	b = -BESLJN(jy, kj, x)*xk*(kd-1.0)/kf;
      }
      _dwork12[i] *= b;    
      b = x*x/((kd-1.0)*(kd+3));
      _dwork12[i] *= b;
    }
    for (; i < potential->maxrp; i++) _dwork12[i] = 0;
    Integrate(_dwork12, orb0, orb2, 6, &s2, 0);
  } else {
    s2 = 0;
  }  
  r = s1 - s2;
  return r;
}

double BreitSK(int k0, int k1, int k2, int k3, int k, int mbr) {
  double r1, r2;
  r1 = BreitSW(k0, k1, k2, k3, k, 0, mbr);
  if ((k0 == k1 && k2 == k3) || (k0 == k3 && k1 == k2)) {
    r2 = r1;
  } else {      
    r2 = BreitSW(k0, k1, k2, k3, k, 1, mbr);
  }
  return 0.5*(r1+r2);
}

double BreitWW(int k0, int k1, int k2, int k3, int k,
	       int kp0, int kp1, int kp2, int kp3,
	       int kl0, int kl1, int kl2, int kl3, int mbr) {
  int m, m0, m1, n, ka, kap, kd;
  double a, b, c, r, c1, c2, c3, c4;

  if (k <= 0) return 0;

  int index[5];
  double *p;
  LOCK *lock = NULL;
  int locked = 0;
  int myrank = MyRankMPI()+1;
  if (wbreit_array->maxsize != 0) {
    index[0] = k0;
    index[1] = k1;
    index[2] = k2;
    index[3] = k3;
    index[4] = k;
    p = (double *) MultiSet(wbreit_array, index, NULL, &lock,
			    InitDoubleData, NULL);
    if (lock && !(p && *p)) {
      SetLock(lock);
      locked = 1;
    }
    if (p && *p) {
      r = *p;
      if (locked) ReleaseLock(lock);
#pragma omp atomic
      wbreit_array->iset -= myrank;
#pragma omp flush
      return r;
    }
  }
  m0 = k - 1;
  if (m0 < 0) m0 = 0;
  m1 = k + 1;
  kd = 2*k;
  ka = kp2 - kp0;
  kap = kp3 - kp1;
  r = 0.0;  
  for (m = m0; m <= m1; m++) {
    if (IsEven((kl0+kl2)/2 + m) || IsEven((kl1+kl3)/2 + m)) continue;
    if (m < k) {
      a = (k+1.0)/(k*(kd-1.0)*(kd+1.0));
      c1 = a*(ka+k)*(kap+k);
      c2 = a*(ka-k)*(kap-k);
      c3 = a*(ka+k)*(kap-k);
      c4 = a*(ka-k)*(kap+k);
    } else if (m == k) {
      a = -(kp0 + kp2)*(kp1 + kp3)/(k*(k+1.0));
      c1 = c2 = c3 = c4 = a;
    } else {
      a = k/((k+1.0)*(kd+1.0)*(kd+3.0));
      c1 = a*(ka-k-1)*(kap-k-1);
      c2 = a*(ka+k+1)*(kap+k+1);
      c3 = a*(ka-k-1)*(kap+k+1);
      c4 = a*(ka+k+1)*(kap-k-1);
    }

    if (fabs(c1) > 1e-30) {
      a = BreitRK(k0, k1, k2, k3, m, mbr);
      r += a*c1;
      //printf("rk1: %d %d %d %d %d %d %g %g %g %d %d %d %d %d %d %d %d %d %d %d\n", k0, k1, k2, k3, k, m, a, c1, r, kl0, kl1, kl2, kl3, kp0, kp1, kp2, kp3, ka, kap, kd);
    }

    if (fabs(c2) > 1e-30) {
      a = BreitRK(k2, k3, k0, k1, m, mbr);
      r += a*c2;
      //printf("rk2: %d %d %d %d %d %d %g %g %g\n", k0, k1, k2, k3, k, m, a, c2, r);
    }

    if (fabs(c3) > 1e-30) {
      a = BreitRK(k0, k3, k2, k1, m, mbr);
      r += a*c3;
      //printf("rk3: %d %d %d %d %d %d %g %g %g\n", k0, k1, k2, k3, k, m, a, c3, r);
    }

    if (fabs(c4) > 1e-30) {
      a = BreitRK(k2, k1, k0, k3, m, mbr);
      r += a*c4;
      //printf("rk4: %d %d %d %d %d %d %g %g %g\n", k0, k1, k2, k3, k, m, a, c4, r);
    }
  }
  
  if (k >= 0 && IsEven((kl0+kl2)/2+k) && IsEven((kl1+kl3)/2+k)) {
    b = 1.0/((kd+1.0)*(kd+1.0));
    c = b*(ka+k)*(kap-k-1.0);
    if (fabs(c) > 1e-30) {
      a = BreitSK(k0, k1, k2, k3, k, mbr);      
      r += a*c;
      //printf("sk1: %d %d %d %d %d %g %g %g %g\n", k0, k1, k2, k3, k, b, a, c, r);
    }
    
    c = b*(kap+k)*(ka-k-1.0);
    if (fabs(c) > 1e-30) {
      a = BreitSK(k1, k0, k3, k2, k, mbr);
      r += a*c;
      //printf("sk2: %d %d %d %d %d %g %g %g %g\n", k0, k1, k2, k3, k, b, a, c, r);
    }
    
    c = b*(ka-k)*(kap+k+1.0);
    if (fabs(c) > 1e-30) {
      a = BreitSK(k2, k3, k0, k1, k, mbr);
      r += a*c;
      //printf("sk3: %d %d %d %d %d %g %g %g %g\n", k0, k1, k2, k3, k, b, a, c, r);
    }
    
    c = b*(kap-k)*(ka+k+1.0);
    if (fabs(c) > 1e-30) {
      a = BreitSK(k3, k2, k1, k0, k, mbr);
      r += a*c;
      //printf("sk4: %d %d %d %d %d %g %g %g %g\n", k0, k1, k2, k3, k, b, a, c, r);
    }
    
    c = b*(ka+k)*(kap+k+1.0);
    if (fabs(c) > 1e-30) {
      a = BreitSK(k0, k3, k2, k1, k, mbr);
      r += a*c;
      //printf("sk5: %d %d %d %d %d %g %g %g %g\n", k0, k1, k2, k3, k, b, a, c, r);
    }
    
    c = b*(kap-k)*(ka-k-1.0);
    if (fabs(c) > 1e-30) {
      a = BreitSK(k3, k0, k1, k2, k, mbr);
      r += a*c;
      //printf("sk6: %d %d %d %d %d %g %g %g %g\n", k0, k1, k2, k3, k, b, a, c, r);
    }
    
    c = b*(ka-k)*(kap-k-1.0);
    if (fabs(c) > 1e-30) {
      a = BreitSK(k2, k1, k0, k3, k, mbr);
      r += a*c;
      //printf("sk7: %d %d %d %d %d %g %g %g %g\n", k0, k1, k2, k3, k, b, a, c, r);
    }
    
    c = b*(kap+k)*(ka+k+1.0);
    if (fabs(c) > 1e-30) {
      a = BreitSK(k1, k2, k3, k0, k, mbr);
      r += a*c;
      //printf("sk8: %d %d %d %d %d %g %g %g %g\n", k0, k1, k2, k3, k, b, a, c, r);
    }
  }  
  if (wbreit_array->maxsize != 0) {
    if (!r) r = 1e-100;
    *p = r;
    if (locked) ReleaseLock(lock);
#pragma omp atomic
    wbreit_array->iset -= myrank;
  }
#pragma omp flush
  return r;
}

double BreitNW(int k0, int k1, int k2, int k3, int k,
	       int kl0, int kl1, int kl2, int kl3) {
  int m, m0, m1, n;
  double a, c, r;

  m0 = k - 1;
  if (m0 < 0) m0 = 0;
  m1 = k + 1;
  r = 0.0;
  for (m = m0; m <= m1; m++) {
    if (IsEven((kl0+kl2)/2 + m) || IsEven((kl1+kl3)/2 + m)) continue;
    for (n = 0; n < 8; n++) {
      c = BreitC(n, m, k, k0, k1, k2, k3);
      if (c) {
	a = BreitI(n, k0, k1, k2, k3, m);
	r += a*c;
      }
    }
  }

  return r;
}

double Breit(int k0, int k1, int k2, int k3, int k,
	     int kp0, int kp1, int kp2, int kp3,
	     int kl0, int kl1, int kl2, int kl3, int mbr) {
  double r;
  if (k <= 0) return 0;
  if (mbr == 0) {
    r = BreitNW(k0, k1, k2, k3, k, kl0, kl1, kl2, kl3);
  } else {
    r = BreitWW(k0, k1, k2, k3, k, kp0, kp1, kp2, kp3, kl0, kl1, kl2, kl3, mbr);
  }
  //printf("bri: %d %d %d %d %d %g %d %d %d %d\n", k0, k1, k2, k3, k, r, kp0, kp1, kp2, kp3);
  return r;
}

/* calculate the slater integral of rank k */
int Slater(double *s, int k0, int k1, int k2, int k3, int k, int mode) {
  int index[5];
  double *p;
  int ilast, i, npts, m;
  ORBITAL *orb0, *orb1, *orb2, *orb3;
  double norm;
#ifdef PERFORM_STATISTICS
  clock_t start, stop; 
  start = clock();
#endif

  index[0] = k0;
  index[1] = k1;
  index[2] = k2;
  index[3] = k3;
  index[4] = k;  

  LOCK *lock = NULL;
  int locked = 0;
  int myrank = MyRankMPI()+1;
  if (abs(mode) < 2) {
    SortSlaterKey(index);
    p = (double *) MultiSet(slater_array, index, NULL, &lock,
			    InitDoubleData, NULL);
    if (lock && !(p && *p)) {
      SetLock(lock);
      locked = 1;
    }
  } else {
    p = NULL;
  }
  if (p && *p) {
    *s = *p;
  } else {
    orb0 = GetOrbitalSolved(k0);
    orb1 = GetOrbitalSolved(k1);
    orb2 = GetOrbitalSolved(k2);
    orb3 = GetOrbitalSolved(k3);
    *s = 0.0; 
    npts = potential->maxrp;
    switch (mode) {
    case 0: /* fall through to case 1 */
    case 1: /* full relativistic with distorted free orbitals */
      GetYk(k, _yk, orb0, orb2, k0, k2, -1); 
      if (orb1->n > 0) ilast = orb1->ilast;
      else ilast = npts-1;
      if (orb3->n > 0) ilast = Min(ilast, orb3->ilast);
      for (i = 0; i <= ilast; i++) {
	_yk[i] = (_yk[i]/potential->rad[i]);
      }
      Integrate(_yk, orb1, orb3, 1, s, 0);
      break;
    
    case -1: /* quasi relativistic with distorted free orbitals */
      GetYk(k, _yk, orb0, orb2, k0, k2, -2);
      if (orb1->n > 0) ilast = orb1->ilast;
      else ilast = npts-1;
      if (orb3->n > 0) ilast = Min(ilast, orb3->ilast);
      for (i = 0; i <= ilast; i++) {
	_yk[i] /= potential->rad[i];
      }
      Integrate(_yk, orb1, orb3, 2, s, 0);

      norm  = orb0->qr_norm;
      norm *= orb1->qr_norm;
      norm *= orb2->qr_norm;
      norm *= orb3->qr_norm;

      *s *= norm;
      break;

    case 2: /* separable coulomb interaction, orb0, orb2 is inner part */
      m = k;
      *s = RadialMoments(m, k0, k2);
      if (*s != 0.0) {
	m = -m-1;
	*s *= RadialMoments(m, k1, k3);
      }
      break;

    case -2: /* separable coulomb interaction, orb1, orb3 is inner part  */
      m = k;
      *s = RadialMoments(m, k1, k3);
      if (*s != 0.0) {
	m = -m-1;
	*s *= RadialMoments(m, k0, k2);      
      }
      break;
      
    default:
      break;
    }      
    if (p) *p = *s;
  }
  if (locked) ReleaseLock(lock);
  if (p) {
#pragma omp atomic
    slater_array->iset -= myrank;
  }
#pragma omp flush
#ifdef PERFORM_STATISTICS 
    stop = clock();
    rad_timing.radial_2e += stop - start;
#endif
  return 0;
}

/* reorder the orbital index appears in the slater integral, so that it is
   in a form: a <= b <= d, a <= c, and if (a == b), c <= d. */ 
void SortSlaterKey(int *kd) {
  int i;

  if (kd[0] > kd[2]) {
    i = kd[0];
    kd[0] = kd[2];
    kd[2] = i;
  }

  if (kd[1] > kd[3]) {
    i = kd[1];
    kd[1] = kd[3];
    kd[3] = i;
  }
  
  if (kd[0] > kd[1]) {
    i = kd[0];
    kd[0] = kd[1];
    kd[1] = i;
    i = kd[2];
    kd[2] = kd[3];
    kd[3] = i;
  } else if (kd[0] == kd[1]) {
    if (kd[2] > kd[3]) {
      i = kd[2];
      kd[2] = kd[3];
      kd[3] = i;
    }
  }
}

void PrepSlater(int ib0, int iu0, int ib1, int iu1,
		int ib2, int iu2, int ib3, int iu3) {
#pragma omp parallel default(shared)
  {
  int k, kmax, kk, i, j, p, q, m, ilast;
  int j0, j1, j2, j3, k0, k1, k2, k3;
  int index[6];
  double *dp;
  ORBITAL *orb0, *orb1, *orb2, *orb3;
  int c = 0;

  int myrank = MyRankMPI()+1;
  double wt0 = WallTime();
  kmax = GetMaxRank();
  for (kk = 0; kk <= kmax; kk += 2) {
    k = kk/2;
    for (i = ib0; i <= iu0; i++) {
      orb0 = GetOrbital(i);
      GetJLFromKappa(orb0->kappa, &j0, &k0);
      for (p = ib2; p <= iu2; p++) {
	if (p < i) continue;
	orb2 = GetOrbital(p);
	GetJLFromKappa(orb2->kappa, &j2, &k2);
	if (k0 > slater_cut.kl0 || k2 > slater_cut.kl0) continue;
	int skip = SkipMPI();
	if (skip) continue;	     
	GetYk(k, _yk, orb0, orb2, i, p, -1);
	ilast = potential->maxrp-1;
	for (m = 0; m <= ilast; m++) {
	  _yk[m] /= potential->rad[m];
	}
	for (j = ib1; j <= iu1; j++) {
	  if (j < i) continue;
	  orb1 = GetOrbital(j);
	  GetJLFromKappa(orb1->kappa, &j1, &k1);
	  for (q = ib3; q <= iu3; q++) {
	    if (q < j) continue;
	    orb3 = GetOrbital(q);
	    GetJLFromKappa(orb3->kappa, &j3, &k3);
	    if (i == j && q < p) continue;
	    if (IsOdd((k0+k2)/2+k) || 
		IsOdd((k1+k3)/2+k) ||
		!Triangle(j0, j2, kk) ||
		!Triangle(j1, j3, kk)) continue;
	    index[0] = i;
	    index[1] = j;
	    index[2] = p;
	    index[3] = q;
	    index[4] = k;
	    index[5] = 0;
	    LOCK *lock = NULL;
	    dp = MultiSet(slater_array, index, NULL, &lock,
			  InitDoubleData, NULL);
	    c++;
	    //if (lock) SetLock(lock);
	    if (*dp == 0) {
	      Integrate(_yk, orb1, orb3, 1, dp, 0);
	    }
	    //if (lock) ReleaseLock(lock);
	    slater_array->iset -= myrank;
	  }
	}
      }
    }
  }
  double wt1 = WallTime();
  MPrintf(-1, "PrepSlater: %d %g\n", c, wt1-wt0);
  }
}
      
int GetYk0(int k, double *yk, ORBITAL *orb1, ORBITAL *orb2, int type) {
  int i, ilast, i0;
  double a, max;

  for (i = 0; i < potential->maxrp; i++) {
    _dwork1[i] = pow(potential->rad[i], k);
  }
  Integrate(_dwork1, orb1, orb2, type, _zk, 0);
  for (i = 0; i < potential->maxrp; i++) {
    _zk[i] /= _dwork1[i];
    yk[i] = _zk[i];
    _zk[i] = _dwork1[i];
  }
  if (k > 2) {
    max = 0.0;
    for (i = 0; i < potential->maxrp; i++) {
      a = fabs(yk[i]);
      if (max < a) max = a;
    }
    max *= 1E-3;
    ilast = Min(orb1->ilast, orb2->ilast);
    for (i = 0; i < ilast; i++) {
      a = Large(orb1)[i]*Large(orb2)[i]*potential->rad[i];
      if (fabs(a) > max) break;
      _dwork1[i] = 0.0;
    }
    i0 = i;
  } else i0 = 0;
  for (i = i0; i < potential->maxrp; i++) {
    _dwork1[i] = pow(potential->rad[i0]/potential->rad[i], k+1);
  }
  Integrate(_dwork1, orb1, orb2, type, _xk, 0);
  ilast = potential->maxrp - 1;    
  for (i = i0; i < potential->maxrp; i++) {
    _xk[i] = (_xk[ilast] - _xk[i])/_dwork1[i];
    yk[i] += _xk[i];
  }

  return 0;
}  

/*
** this is a better version of Yk than GetYk0.
** note that on exit, _zk consists r^k, which is used in GetYk
*/      
int GetYk1(int k, double *yk, ORBITAL *orb1, ORBITAL *orb2, int type) {
  int i, ilast;
  double r0, a;
  
  ilast = Min(orb1->ilast, orb2->ilast);
  r0 = sqrt(potential->rad[0]*potential->rad[ilast]);  
  for (i = 0; i < potential->maxrp; i++) {
    _dwork1[i] = pow(potential->rad[i]/r0, k);
  }
  Integrate(_dwork1, orb1, orb2, type, _zk, 0);
  a = pow(r0, k);
  for (i = 0; i < potential->maxrp; i++) {
    _zk[i] /= _dwork1[i];
    yk[i] = _zk[i];
    _zk[i] = _dwork1[i]*a;
  }  
  for (i = 0; i < potential->maxrp; i++) {
    _dwork1[i] = (r0/potential->rad[i])/_dwork1[i];
  }
  Integrate(_dwork1, orb1, orb2, type, _xk, -1);
  for (i = 0; i < potential->maxrp; i++) {
    yk[i] += _xk[i]/_dwork1[i];
  }
      
  return 0;
}
      
int GetYk(int k, double *yk, ORBITAL *orb1, ORBITAL *orb2, 
	  int k1, int k2, int type) {
  int i, i0, i1, n, npts, ic0, ic1;
  double a, b, a2, b2, max, max1;
  int index[3];
  FLTARY *syk;

  syk = NULL;
  LOCK *lock = NULL;
  int locked = 0;
  int myrank = MyRankMPI()+1;
  if (yk_array->maxsize != 0) {
    if (k1 <= k2) {
      index[0] = k1;
      index[1] = k2;
    } else {
      index[0] = k2;
      index[1] = k1;
    }
    index[2] = k;
    syk = (FLTARY *) MultiSet(yk_array, index, NULL, &lock,
			      InitFltAryData, FreeFltAryData);
    if (lock && syk->npts <= 0) {
      SetLock(lock);
      locked = 1;
    }
    if (syk->npts > 0) {
      npts = syk->npts-2;
      for (i = npts-1; i < potential->maxrp; i++) {
	_dwork1[i] = pow(potential->rad[i], k);
      }
      for (i = 0; i < npts; i++) {
	yk[i] = syk->yk[i];
      }
      ic0 = npts;
      ic1 = npts+1;
      i0 = npts-1;
      a = syk->yk[i0]*_dwork1[i0];
      for (i = npts; i < potential->maxrp; i++) {
	b = potential->rad[i] - potential->rad[i0];
	b = syk->yk[ic1]*b;
	if (b < -20) {
	  yk[i] = syk->yk[ic0];
	} else {
	  yk[i] = (a - syk->yk[ic0])*exp(b);
	  yk[i] += syk->yk[ic0];
	}
	yk[i] /= _dwork1[i];
      }    
    }
  }
  if (syk == NULL || syk->npts <= 0) {
    GetYk1(k, yk, orb1, orb2, type);
    max = 0;
    for (i = 0; i < potential->maxrp; i++) {
      _zk[i] *= yk[i];
      a = fabs(_zk[i]); 
      if (a > max) max = a;
    }
    max1 = max*EPS5;
    max = max*EPS4;
    int maxrp = potential->maxrp;
    a = _zk[potential->maxrp-1];
    for (i = potential->maxrp-2; i >= 0; i--) {
      if (fabs(_zk[i] - a) > max1) {
	break;
      }
    }
    i1 = i;
    for (i = i1; i >= 0; i--) {      
      b = fabs(a - _zk[i]);
      _zk[i] = log(b);
      if (b > max) {
	break;
      }
    }
    i0 = i;
    if (i0 == i1) {
      i0--;
      b = fabs(a - _zk[i0]);
      _zk[i0] = log(b);
    } 
    npts = i0+1;
    ic0 = npts;
    ic1 = npts+1;
    if (syk != NULL) {
      int size = sizeof(float)*(npts+2);
      syk->yk = malloc(size);
      AddMultiSize(yk_array, size);
      for (i = 0; i < npts ; i++) {
	syk->yk[i] = yk[i];
      }
      syk->yk[ic0] = a;
      n = i1 - i0 + 1;
      a = 0.0;
      b = 0.0;
      a2 = 0.0;
      b2 = 0.0;
      for (i = i0; i <= i1; i++) {      
	max = (potential->rad[i]-potential->rad[i0]);
	a += max;
	b += _zk[i];
	a2 += max*max;
	b2 += _zk[i]*max;
      }
      syk->yk[ic1] = (a*b - n*b2)/(a*a - n*a2);       
      if (syk->yk[ic1] >= 0) {
	i1 = i0 + (i1-i0)*0.3;
	if (i1 == i0) i1 = i0 + 1;
	for (i = i0; i <= i1; i++) {      
	  max = (potential->rad[i]-potential->rad[i0]);
	  a += max;
	  b += _zk[i];
	  a2 += max*max;
	  b2 += _zk[i]*max;
	}
	syk->yk[ic1] = (a*b - n*b2)/(a*a - n*a2);
      }
      if (syk->yk[ic1] >= 0) {
	syk->yk[ic1] = -10.0/(potential->rad[i1]-potential->rad[i0]);
      }
      syk->npts = npts+2;
    }
  }
  if (locked) ReleaseLock(lock);  
  if (yk_array->maxsize != 0) {
#pragma omp atomic
    yk_array->iset -= myrank;
  }
#pragma omp flush
  return 0;
}  

/* integrate a function given by f with two orbitals. */
/* type indicates the type of integral */
/* type = 1,    P1*P2 + Q1*Q2 */
/* type = 2,    P1*P2 */
/* type = 3,    Q1*Q2 */ 
/* type = 4:    P1*Q2 + Q1*P2 */
/* type = 5:    P1*Q2 - Q1*P2 */
/* type = 6:    P1*Q2 */
/* if type is positive, only the end point is returned, */
/* otherwise, the whole function is returned */
/* id indicate whether integrate inward (-1) or outward (0) */
int Integrate(double *f, ORBITAL *orb1, ORBITAL *orb2, 
	      int t, double *x, int id) {
  int i1, i2, ilast;
  int i, type;
  double *r, ext;

  if (t == 0) t = 1;
  if (t < 0) {
    r = x;
    type = -t;
  } else {
    r = _dwork;
    type = t;
  }
  for (i = 0; i < potential->maxrp; i++) {
    r[i] = 0.0;
  }

  ext = 0.0;
  ilast = Min(orb1->ilast, orb2->ilast);
  if (id >= 0) {
    IntegrateSubRegion(0, ilast, f, orb1, orb2, t, r, 0, NULL);
  } else {
    IntegrateSubRegion(0, ilast, f, orb1, orb2, t, r, -1, NULL);
  }
  i2 = ilast;
  if (orb1->ilast == ilast && orb1->n == 0) {
    i1 = ilast + 1;
    i2 = orb2->ilast;
    if (i2 > i1) {
      if (type == 6) {
	i = 7;
	if (t < 0) i = -i;
	IntegrateSubRegion(i1, i2, f, orb2, orb1, i, r, 1, NULL);
      } else {
	IntegrateSubRegion(i1, i2, f, orb1, orb2, t, r, 2, NULL);
      }
      i2--;
    }
    if (orb2->n == 0) {
      i1 = orb2->ilast + 1;
      i2 = potential->maxrp - 1;
      IntegrateSubRegion(i1, i2, f, orb1, orb2, t, r, 3, &ext);
      i2--;
    }
  } else if (orb2->ilast == ilast && orb2->n == 0) {
    i1 = ilast + 1;
    i2 = orb1->ilast;
    if (i2 > i1) {
      IntegrateSubRegion(i1, i2, f, orb1, orb2, t, r, 1, NULL);
      i2--;
    }
    if (orb1->n == 0) {
      i1 = orb1->ilast + 1;
      i2 = potential->maxrp - 1;
      IntegrateSubRegion(i1, i2, f, orb1, orb2, t, r, 3, &ext);
      i2--;
    }
  }
  if (t >= 0) {
    *x = r[i2] + ext;
  } else {
    if (id >= 0) {
      for (i = i2+1; i < potential->maxrp; i++) {
	r[i] = r[i2];
      }
    } else {
      if (i2 > ilast) {
	ext += r[i2];
	for (i = ilast + 1; i <= i2; i++) {
	  r[i] = ext - r[i];
	}
	for (i = i2+1; i < potential->maxrp; i++) {
	  r[i] = 0.0;
	}
	for (i = 0; i <= ilast; i++) {
	  r[i] = r[ilast+1] + r[i];
	}
      }
    }
  }
  
  return 0;
}

void AddEvenPoints(double *r, double *r1, int i0, int i1, int t) {
  int i;

  if (t < 0) {
    for (i = i0; i < i1; i++) {
      r[i] += r1[i];
    }
  } else {
    r[i1-1] += r1[i1-1];
  }
}

int IntegrateSubRegion(int i0, int i1, 
		       double *f, ORBITAL *orb1, ORBITAL *orb2,
		       int t, double *r, int m, double *ext) {
  int i, j, ip, i2, type;
  ORBITAL *tmp;
  double *large1, *large2, *small1, *small2;
  double *x, *y, *r1, *x1, *x2, *y1, *y2;
  double a, b, e1, e2, a2, r0;
  int kv1 = orb1->kv;
  int kv2 = orb2->kv;
  
  if (i1 <= i0) return 0;
  type = abs(t);

  x = _dwork3;
  y = _dwork4;
  x1 = _dwork5;
  x2 = _dwork6;
  y1 = _dwork7;
  y2 = _dwork8;
  r1 = _dwork9;
  i2 = i1;

  switch (m) {
  case -1: /* m = -1 same as m = 0 but integrate inward */
  case 0:  /* m = 0 */
    large1 = Large(orb1);
    large2 = Large(orb2);
    small1 = Small(orb1);
    small2 = Small(orb2);
    switch (type) {
    case 1: /* type = 1 */
      for (i = i0; i <= i1; i++) {
	x[i] = large1[i] * large2[i];
	x[i] += small1[i] * small2[i];
	x[i] *= f[i]*potential->dr_drho[i];
      }
      if (i1 == orb1->ilast && orb1->n == 0 && i < potential->maxrp) {
	if (i <= orb2->ilast) {
	  ip = i+1;
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  x[i] = large1[i]*a*large2[i];
	  x[i] += (small1[i]*b + small1[ip]*a)*small2[i];
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	} else if (orb2->n == 0) {
	  ip = i+1;
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  e1 = sin(large2[ip]);
	  e2 = cos(large2[ip]);
	  x[i] = large1[i]*a*large2[i]*e1;
	  x[i] += (small1[i]*b+small1[ip]*a)*(small2[i]*e2+small2[ip]*e1);
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	}
      } else if (i1 == orb2->ilast && orb2->n == 0 && i < potential->maxrp) {
	if (i <= orb1->ilast) {
	  ip = i+1;
	  a = sin(large2[ip]);
	  b = cos(large2[ip]);
	  x[i] = large2[i]*a*large1[i];
	  x[i] += (small2[i]*b + small2[ip]*a)*small1[i];
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	} else if (orb1->n == 0) {
	  ip = i+1;
	  a = sin(large2[ip]);
	  b = cos(large2[ip]);
	  e1 = sin(large1[ip]);
	  e2 = cos(large1[ip]);
	  x[i] = large2[i]*a*large1[i]*e1;
	  x[i] += (small2[i]*b+small2[ip]*a)*(small1[i]*e2+small1[ip]*e1);
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	}
      }
      break;
    case 2: /* type = 2 */
      for (i = i0; i <= i1; i++) {
	x[i] = large1[i] * large2[i];
	x[i] *= f[i]*potential->dr_drho[i];
      }
      if (i1 == orb1->ilast && orb1->n == 0 && i < potential->maxrp) {
	if (i <= orb2->ilast) {
	  ip = i+1;
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  x[i] = large1[i]*a*large2[i];
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	} else if (orb2->n == 0) {
	  ip = i+1;
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  e1 = sin(large2[ip]);
	  e2 = cos(large2[ip]);
	  x[i] = large1[i]*a*large2[i]*e1;
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	}
      } else if (i1 == orb2->ilast && orb2->n == 0 && i < potential->maxrp) {
	if (i <= orb1->ilast) {
	  ip = i+1;
	  a = sin(large2[ip]);
	  b = cos(large2[ip]);
	  x[i] = large2[i]*a*large1[i];
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	} else if (orb1->n == 0) {
	  ip = i+1;
	  a = sin(large2[ip]);
	  b = cos(large2[ip]);
	  e1 = sin(large1[ip]);
	  e2 = cos(large1[ip]);
	  x[i] = large2[i]*a*large1[i]*e1;
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	}
      }
      break;
    case 3: /*type = 3 */
      for (i = i0; i <= i1; i++) {
	x[i] = small1[i] * small2[i];
	x[i] *= f[i]*potential->dr_drho[i];
      }
      if (i1 == orb1->ilast && orb1->n == 0 && i < potential->maxrp) {
	if (i <= orb2->ilast) {
	  ip = i+1;
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  x[i] = (small1[i]*b + small1[ip]*a)*small2[i];
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	} else if (orb2->n == 0) {
	  ip = i+1;
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  e1 = sin(large2[ip]);
	  e2 = cos(large2[ip]);
	  x[i] = (small1[i]*b+small1[ip]*a)*(small2[i]*e2+small2[ip]*e1);
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	}
      } else if (i1 == orb2->ilast && orb2->n == 0 && i < potential->maxrp) {
	if (i <= orb1->ilast) {
	  ip = i+1;
	  a = sin(large2[ip]);
	  b = cos(large2[ip]);
	  x[i] = (small2[i]*b + small2[ip]*a)*small1[i];
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	} else if (orb1->n == 0) {
	  ip = i+1;
	  a = sin(large2[ip]);
	  b = cos(large2[ip]);
	  e1 = sin(large1[ip]);
	  e2 = cos(large1[ip]);
	  x[i] = (small2[i]*b+small2[ip]*a)*(small1[i]*e2+small1[ip]*e1);
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	}
      }
      break;
    case 4: /*type = 4 */
      for (i = i0; i <= i1; i++) {
	x[i] = large1[i] * small2[i];
	x[i] += small1[i] * large2[i];
	x[i] *= f[i]*potential->dr_drho[i];
      }
      if (i1 == orb1->ilast && orb1->n == 0 && i < potential->maxrp) {
	if (i <= orb2->ilast) {
	  ip = i+1;
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  x[i] = large1[i]*a*small2[i];
	  x[i] += (small1[i]*b + small1[ip]*a)*large2[i];
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	} else if (orb2->n == 0) {
	  ip = i+1;
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  e1 = sin(large2[ip]);
	  e2 = cos(large2[ip]);
	  x[i] = large1[i]*a*(small2[i]*e2+small2[ip]*e1);
	  x[i] += (small1[i]*b+small1[ip]*a)*large2[i]*e1;
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	}
      } else if (i1 == orb2->ilast && orb2->n == 0 && i < potential->maxrp) {
	if (i <= orb1->ilast) {
	  ip = i+1;
	  a = sin(large2[ip]);
	  b = cos(large2[ip]);
	  x[i] = (small2[i]*b + small2[ip]*a)*large1[i];
	  x[i] += large2[i]*a*small1[i];
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	} else if (orb1->n == 0) {
	  ip = i+1;
	  a = sin(large2[ip]);
	  b = cos(large2[ip]);
	  e1 = sin(large1[ip]);
	  e2 = cos(large1[ip]);
	  x[i] = (small2[i]*b+small2[ip]*a)*large1[i]*e1;
	  x[i] += large2[i]*a*(small1[i]*e2+small1[ip]*e1);
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	}
      }
      break;
    case 5: /* type = 5 */
      for (i = i0; i <= i1; i++) {
	x[i] = large1[i] * small2[i];
	x[i] -= small1[i] * large2[i];
	x[i] *= f[i]*potential->dr_drho[i];
      } 
      if (i1 == orb1->ilast && orb1->n == 0 && i < potential->maxrp) {
	if (i <= orb2->ilast) {
	  ip = i+1;
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  x[i] = large1[i]*a*small2[i];
	  x[i] -= (small1[i]*b + small1[ip]*a)*large2[i];
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	} else if (orb2->n == 0) {
	  ip = i+1;
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  e1 = sin(large2[ip]);
	  e2 = cos(large2[ip]);
	  x[i] = large1[i]*a*(small2[i]*e2+small2[ip]*e1);
	  x[i] -= (small1[i]*b+small1[ip]*a)*large2[i]*e1;
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	}
      } else if (i1 == orb2->ilast && orb2->n == 0 && i < potential->maxrp) {
	if (i <= orb1->ilast) {
	  ip = i+1;
	  a = sin(large2[ip]);
	  b = cos(large2[ip]);
	  x[i] = (small2[i]*b + small2[ip]*a)*large1[i];
	  x[i] -= large2[i]*a*small1[i];
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	} else if (orb1->n == 0) {
	  ip = i+1;
	  a = sin(large2[ip]);
	  b = cos(large2[ip]);
	  e1 = sin(large1[ip]);
	  e2 = cos(large1[ip]);
	  x[i] = (small2[i]*b+small2[ip]*a)*large1[i]*e1;
	  x[i] -= large2[i]*a*(small1[i]*e2+small1[ip]*e1);
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	}
      }
      break;
    case 6: /* type = 6 */
      for (i = i0; i <= i1; i++) {
	x[i] = large1[i] * small2[i];
	x[i] *= f[i]*potential->dr_drho[i];
      }
      if (i1 == orb1->ilast && orb1->n == 0 && i < potential->maxrp) {
	if (i <= orb2->ilast) {
	  ip = i+1;
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  x[i] = large1[i]*a*small2[i];
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	} else if (orb2->n == 0) {
	  ip = i+1;
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  e1 = sin(large2[ip]);
	  e2 = cos(large2[ip]);
	  x[i] = large1[i]*a*(small2[i]*e2+small2[ip]*e1);
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	}
      } else if (i1 == orb2->ilast && orb2->n == 0 && i < potential->maxrp) {
	if (i <= orb1->ilast) {
	  ip = i+1;
	  a = sin(large2[ip]);
	  b = cos(large2[ip]);
	  x[i] = (small2[i]*b + small2[ip]*a)*large1[i];
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	} else if (orb1->n == 0) {
	  ip = i+1;
	  a = sin(large2[ip]);
	  b = cos(large2[ip]);
	  e1 = sin(large1[ip]);
	  e2 = cos(large1[ip]);
	  x[i] = (small2[i]*b+small2[ip]*a)*large1[i]*e1;
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	}
      }
      break;
    default: /* error */
      return -1;
    }      
    NewtonCotes(r, x, i0, i2, t, m);
    break;

  case 1: /* m = 1 */
    if (type == 6) { /* type 6 needs special treatments */
      large1 = Large(orb1);
      large2 = Large(orb2);
      small1 = Small(orb1);
      small2 = Small(orb2);
      e1 = orb1->energy;
      e2 = orb2->energy;
      a2 = 0.5*FINE_STRUCTURE_CONST2;
      j = 0;
      for (i = i0; i <= i1; i+= 2) {
	ip = i+1;
	x[j] = large1[i] * small2[ip];
	x[j] *= f[i];
	y[j] = large1[i] * small2[i] * f[i];
	_phase[j] = large2[ip];
	_dphase[j] = (1+a2*e2-potential->VT[kv2][i])/(large2[i]*large2[i]);
	j++;
      }
      if (i < potential->maxrp) {
	ip = i+1;
	if (i > orb1->ilast && orb1->n == 0) { 
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  b = small1[i]*b + small1[ip]*a;
	  a = large1[i]*a;
	  x[j] = a * small2[ip];
	  x[j] *= f[i];
	  y[j] = a * small2[i] * f[i];
	  _phase[j] = large2[ip];
	  _dphase[j] = (1+a2*(e2-potential->VT[kv2][i]))/(large2[i]*large2[i]);
	  j++;
	  i2 = i;
	}
      }
      IntegrateSinCos(j, x, y, _phase, _dphase, i0, r, t, ext);
      if (IsOdd(i2)) r[i2] = r[i2-1];
      break;
    } else if (type == 7) {
      large1 = Large(orb1);
      large2 = Large(orb2);
      small1 = Small(orb1);
      small2 = Small(orb2);
      e1 = orb1->energy;
      e2 = orb2->energy;
      a2 = 0.5*FINE_STRUCTURE_CONST2;
      j = 0;
      for (i = i0; i <= i1; i+= 2) {
	ip = i+1;
	x[j] = small1[i] * large2[i];
	x[j] *= f[i];
	_phase[j] = large2[ip];
	_dphase[j] = (1+a2*(e2-potential->VT[kv2][i]))/(large2[i]*large2[i]);
	j++;	
      }
      if (i < potential->maxrp) {
	ip = i+1;
	if (i > orb1->ilast && orb1->n == 0) { 
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  b = small1[i]*b + small1[ip]*a;
	  x[j] = b * large2[i];
	  x[j] *= f[i];
	  _phase[j] = large2[ip];
	  _dphase[j] = (1+a2*(e2-potential->VT[kv2][i]))/(large2[i]*large2[i]);
	  j++;
	  i2 = i;
	}
      }
      IntegrateSinCos(j, x, NULL, _phase, _dphase, i0, r, t, ext);
      if (IsOdd(i2)) r[i2] = r[i2-1];
      break;
    }
  case 2: /* m = 1, 2 are essentially the same */
    /* type 6 is treated in m = 1 */
    if (m == 2) {
      tmp = orb1;
      orb1 = orb2;
      orb2 = tmp;
    }
    large1 = Large(orb1);
    large2 = Large(orb2);
    small1 = Small(orb1);
    small2 = Small(orb2);
    e1 = orb1->energy;
    e2 = orb2->energy;
    a2 = 0.5*FINE_STRUCTURE_CONST2;
    switch (type) {
    case 1: /* type = 1 */
      j = 0;
      for (i = i0; i <= i1; i+= 2) {
	ip = i+1;
	x[j] = large1[i] * large2[i];
	x[j] += small1[i] * small2[ip];
	x[j] *= f[i];
	y[j] = small1[i] * small2[i] * f[i];
	_phase[j] = large2[ip];
	_dphase[j] = (1+a2*(e2-potential->VT[kv2][i]))/(large2[i]*large2[i]);
	j++;
      }
      if (i < potential->maxrp) {
	ip = i+1;
	if (i > orb1->ilast && orb1->n == 0) {
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  b = small1[i]*b + small1[ip]*a;
	  a = large1[i]*a;
	  x[j] = a * large2[i];
	  x[j] += b * small2[ip];
	  x[j] *= f[i];
	  y[j] = b * small2[i] * f[i];
	  _phase[j] = large2[ip];
	  _dphase[j] = (1+a2*(e2-potential->VT[kv2][i]))/(large2[i]*large2[i]);
	  j++;
	  i2 = i;
	}
      }
      IntegrateSinCos(j, x, y, _phase, _dphase, i0, r, t, ext);
      if (IsOdd(i2)) r[i2] = r[i2-1];
      break;
    case 2: /* type = 2 */
      j = 0;
      for (i = i0; i <= i1; i+= 2) {
	ip = i+1;
	x[j] = large1[i] * large2[i];
	x[j] *= f[i]; 
	_phase[j] = large2[ip];
	_dphase[j] = (1+a2*(e2-potential->VT[kv2][i]))/(large2[i]*large2[i]);
	j++;
      }
      if (i < potential->maxrp) {
	ip = i+1;
	if (i > orb1->ilast && orb1->n == 0) {
	  a = sin(large1[ip]);
	  a = large1[i]*a;
	  x[j] = a * large2[i];
	  x[j] *= f[i];
	  _phase[j] = large2[ip];
	  _dphase[j] = (1+a2*(e2-potential->VT[kv2][i]))/(large2[i]*large2[i]);
	  j++;
	  i2 = i;
	}
      }
      IntegrateSinCos(j, x, NULL, _phase, _dphase, i0, r, t, ext);
      if (IsOdd(i2)) r[i2] = r[i2-1];
      break;
    case 3: /* type = 3 */
      j = 0;
      for (i = i0; i <= i1; i+= 2) {
	ip = i+1;
	x[j] = small1[i] * small2[ip];
	x[j] *= f[i];
	y[j] = small1[i] * small2[i] * f[i];
	_phase[j] = large2[ip];
	_dphase[j] = (1+a2*(e2-potential->VT[kv2][i]))/(large2[i]*large2[i]);
	j++;
      }
      if (i < potential->maxrp) {
	ip = i+1;
	if (i > orb1->ilast && orb1->n == 0) {
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  b = small1[i]*b + small1[ip]*a;
	  x[j] += b * small2[ip];
	  x[j] *= f[i];
	  y[j] = b * small2[i] * f[i];
	  _phase[j] = large2[ip];
	  _dphase[j] = (1+a2*(e2-potential->VT[kv2][i]))/(large2[i]*large2[i]);
	  j++;
	  i2 = i;
	}
      }
      IntegrateSinCos(j, x, y, _phase, _dphase, i0, r, t, ext);
      if (IsOdd(i2)) r[i2] = r[i2-1];
      break;  
    case 4: /* type = 4 */
      j = 0;
      for (i = i0; i <= i1; i+= 2) {
	ip = i+1;
	x[j] = large1[i] * small2[ip];
	x[j] += small1[i] * large2[i];
	x[j] *= f[i];
	y[j] = large1[i] * small2[i] * f[i];
	_phase[j] = large2[ip];
	_dphase[j] = (1+a2*(e2-potential->VT[kv2][i]))/(large2[i]*large2[i]);
	j++;
      }
      if (i < potential->maxrp) {
	ip = i+1;
	if (i > orb1->ilast && orb1->n == 0) {
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  b = small1[i]*b + small1[ip]*a;
	  a = large1[i]*a;
	  x[j] = a * small2[ip];
	  x[j] += b * large2[i];
	  x[j] *= f[i];
	  y[j] = a * small2[i] * f[i];
	  _phase[j] = large2[ip];
	  _dphase[j] = (1+a2*(e2-potential->VT[kv2][i]))/(large2[i]*large2[i]);
	  j++;
	  i2 = i;
	}
      }
      IntegrateSinCos(j, x, y, _phase, _dphase, i0, r, t, ext);
      if (IsOdd(i2)) r[i2] = r[i2-1];
      break;
    case 5: /* type = 5 */
      j = 0;
      if (m == 2) {
	r0 = r[i0];
	r[i0] = 0.0;
      }
      for (i = i0; i <= i1; i+= 2) {
	ip = i+1;
	x[j] = large1[i] * small2[ip];
	x[j] -= small1[i] * large2[i];
	x[j] *= f[i];
	y[j] = large1[i] * small2[i] * f[i];
	_phase[j] = large2[ip];
	_dphase[j] = (1+a2*(e2-potential->VT[kv2][i]))/(large2[i]*large2[i]);
	j++;
      }
      if (i < potential->maxrp) {
	ip = i+1;
	if (i > orb1->ilast && orb1->n == 0) { 
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  b = small1[i]*b + small1[ip]*a;
	  a = large1[i]*a;
	  x[j] = a * small2[ip];
	  x[j] -= b * large2[i];
	  x[j] *= f[i];
	  y[j] = a * small2[i] * f[i];
	  _phase[j] = large2[ip];
	  _dphase[j] = (1+a2*(e2-potential->VT[kv2][i]))/(large2[i]*large2[i]);
	  j++;
	  i2 = i;
	}
      }
      IntegrateSinCos(j, x, y, _phase, _dphase, i0, r, t, ext);
      if (m == 2) {
	if (t < 0) {
	  for (i = i0; i <= i2; i++) {
	    r[i] = r0 - r[i];
	  }
	} else {
	  if (IsOdd(i2)) r[i2-1] = r0 - r[i2-1];
	  else r[i2] = r0 - r[i2];
	}
      }
      if (IsOdd(i2)) r[i2] = r[i2-1];
      break;
    default:
      return -1;
    }
    break;

  case 3: /* m = 3 */
    large1 = Large(orb1);
    large2 = Large(orb2);
    small1 = Small(orb1);
    small2 = Small(orb2);
    e1 = orb1->energy;
    e2 = orb2->energy;
    a2 = 0.5*FINE_STRUCTURE_CONST2;
    r1[i0] = 0.0;
    switch (type) {
    case 1: /* type = 1 */
      j = 0;
      for (i = i0; i <= i1; i += 2) {
	ip = i+1;	
	x1[j] = small1[i] * small2[ip];
	x2[j] = small1[ip] * small2[i];	
	x[j] = x1[j] + x2[j];
	x[j] *= 0.5*f[i];
	y1[j] = small1[i] * small2[i];
	y2[j] = small1[ip] * small2[ip];
	y2[j] += large1[i] * large2[i];
	y[j] = y1[j] - y2[j];
	y[j] *= 0.5*f[i];
	_phase[j] = large1[ip] + large2[ip];
	a = (1+a2*(e1-potential->VT[kv1][i]))/(large1[i]*large1[i]);
	b = (1+a2*(e2-potential->VT[kv2][i]))/(large2[i]*large2[i]);
	_dphase[j] = a + b;
	_dphasep[j] = a - b;
	j++;
      }
      IntegrateSinCos(j, x, y, _phase, _dphase, i0, r, t, ext);
      j = 0;
      for (i = i0; i <= i1; i += 2) {
	ip = i+1;
	x[j] = -x1[j] + x2[j];
	x[j] *= 0.5*f[i];
	y[j] = y1[j] + y2[j];
	y[j] *= 0.5*f[i];
	_phase[j] = large1[ip] - large2[ip];
	j++;
      }
      IntegrateSinCos(j, NULL, y, _phase, _dphasep, i0, r1, t, NULL);
      AddEvenPoints(r, r1, i0, i1, t);
      break;	
    case 2: /* type = 2 */
      j = 0;
      for (i = i0; i <= i1; i += 2) {
	ip = i+1;
	y2[j] = -large1[i] * large2[i];
	y2[j] *= 0.5*f[i];
	y[j] = y2[j];
	_phase[j] = large1[ip] + large2[ip];
	a = (1+a2*(e1-potential->VT[kv1][i]))/(large1[i]*large1[i]);
	b = (1+a2*(e2-potential->VT[kv2][i]))/(large2[i]*large2[i]);
	_dphase[j] = a + b;
	_dphasep[j] = a - b;
	j++;
      } 
      IntegrateSinCos(j, NULL, y, _phase, _dphase, i0, r, t, ext);
      j = 0;
      for (i = i0; i <= i1; i += 2) {
	ip = i+1;
	y[j] = -y2[j];
	_phase[j] = large1[ip] - large2[ip];
	j++;
      }
      IntegrateSinCos(j, NULL, y, _phase, _dphasep, i0, r1, t, NULL);
      AddEvenPoints(r, r1, i0, i1, t);
      break;
    case 3: /* type = 3 */
      j = 0;
      for (i = i0; i <= i1; i += 2) {
	ip = i+1;
	x1[j] = small1[i] * small2[ip];
	x2[j] = small1[ip] * small2[i];
	x[j] = x1[j] + x2[j];
	x[j] *= 0.5*f[i];
	y1[j] = small1[i] * small2[i];
	y2[j] = small1[ip] * small2[ip];
	y[j] = y1[j] - y2[j];
	y[j] *= 0.5*f[i];
	_phase[j] = large1[ip] + large2[ip];
	a = (1+a2*(e1-potential->VT[kv1][i]))/(large1[i]*large1[i]);
	b = (1+a2*(e2-potential->VT[kv2][i]))/(large2[i]*large2[i]);
	_dphase[j] = a + b;
	_dphasep[j] = a - b;
	j++;
      }
      IntegrateSinCos(j, x, y, _phase, _dphase, i0, r, t, ext);
      j = 0;
      for (i = i0; i <= i1; i += 2) {
	ip = i+1;
	x[j] = -x1[j] + x2[j];
	x[j] *= 0.5*f[i];
	y[j] = y1[j] + y2[j];
	y[j] *= 0.5*f[i];
	_phase[j] = large1[ip] - large2[ip];
	j++;
      }
      IntegrateSinCos(j, x, y, _phase, _dphasep, i0, r1, t, NULL);
      AddEvenPoints(r, r1, i0, i1, t);
      break;
    case 4: /* type = 4 */
      j = 0;
      for (i = i0; i <= i1; i += 2) {
	ip = i+1;
	x1[j] = large1[i] * small2[i];
	x2[j] = small1[i] * large2[i];
	x[j] = x1[j] + x2[j];
	x[j] *= 0.5*f[i];
	y[j] = -large1[i] * small2[ip];
	y[j] -= small1[ip] * large2[i];
	y[j] *= 0.5*f[i];
	y2[j] = y[j];
	_phase[j] = large1[ip] + large2[ip];	
	a = (1+a2*(e1-potential->VT[kv1][i]))/(large1[i]*large1[i]);
	b = (1+a2*(e2-potential->VT[kv2][i]))/(large2[i]*large2[i]);
	_dphase[j] = a + b;
	_dphasep[j] = a - b;
	j++;
      }
      IntegrateSinCos(j, x, y, _phase, _dphase, i0, r, t, ext);
      j = 0;
      for (i = i0; i <= i1; i += 2) {
	ip = i+1;
	x[j] = x1[j] - x2[j];
	x[j] *= 0.5*f[i];
	y[j] = -y2[j];
	_phase[j] = large1[ip] - large2[ip];
	j++;
      }
      IntegrateSinCos(j, x, y, _phase, _dphasep, i0, r1, t, NULL);
      AddEvenPoints(r, r1, i0, i1, t);
      break;
    case 5: /* type = 5 */
      j = 0;
      for (i = i0; i <= i1; i += 2) {
	ip = i+1;
	x1[j] = large1[i] * small2[i];
	x2[j] = small1[i] * large2[i];
	x[j] = x1[j] - x2[j];
	x[j] *= 0.5*f[i];
	y[j] = -large1[i] * small2[ip];
	y[j] += small1[ip] * large2[i];
	y[j] *= 0.5*f[i];
	y2[j] = y[j];
	_phase[j] = large1[ip] + large2[ip];
	a = (1+a2*(e1-potential->VT[kv1][i]))/(large1[i]*large1[i]);
	b = (1+a2*(e2-potential->VT[kv2][i]))/(large2[i]*large2[i]);
	_dphase[j] = a + b;
	_dphasep[j] = a - b;
	j++;
      }
      IntegrateSinCos(j, x, y, _phase, _dphase, i0, r, t, ext);
      j = 0;
      for (i = i0; i <= i1; i += 2) {
	ip = i+1;
	x[j] = x1[j] + x2[j];
	x[j] *= 0.5*f[i];
	y[j] = -y2[j];
	_phase[j] = large1[ip] - large2[ip];
	j++;
      }
      IntegrateSinCos(j, x, y, _phase, _dphasep, i0, r1, t, NULL);
      AddEvenPoints(r, r1, i0, i1, t);
      break;
    case 6: /* type = 6 */
      j = 0;
      for (i = i0; i <= i1; i += 2) {
	ip = i+1;
	x1[j] = large1[i] * small2[i];
	x1[j] *= 0.5*f[i];
	x[j] = x1[j];
	y[j] = -large1[i] * small2[ip];
	y[j] *= 0.5*f[i];
	y2[j] = y[j];
	_phase[j] = large1[ip] + large2[ip];
	a = (1+a2*(e1-potential->VT[kv1][i]))/(large1[i]*large1[i]);
	b = (1+a2*(e2-potential->VT[kv2][i]))/(large2[i]*large2[i]);
	_dphase[j] = a + b;
	_dphasep[j] = a - b;
	j++;
      }
      IntegrateSinCos(j, x, y, _phase, _dphase, i0, r, t, ext);
      j = 0;
      for (i = i0; i <= i1; i += 2) {
	ip = i+1;
	x[j] = x1[j];
	y[j] = -y2[j];
	_phase[j] = large1[ip] - large2[ip];
	j++;
      }
      IntegrateSinCos(j, x, y, _phase, _dphasep, i0, r1, t, NULL);
      AddEvenPoints(r, r1, i0, i1, t);
      break;
    default:
      return -1;
    }
  default:
    return -1; 
  }

  return 0;
}

int IntegrateSinCos(int j, double *x, double *y, 
		    double *phase, double *dphase, 
		    int i0, double *r, int t, double *ext) {
  int i, k, m, n, q, s, i1;
  double si0, si1, cs0, cs1;
  double is0, is1, is2, is3;
  double ic0, ic1, ic2, ic3;
  double his0, his1, his2, his3;
  double hic0, hic1, hic2, hic3;
  double a0, a1, a2, a3;
  double d, p, h, dr;
  double *z, *u, *w, *u1, *w1;

  w = _dwork2;
  z = _dwork10;
  u = _dwork11;
  if (phase[j-1] < 0) {
    for (i = 0; i < j; i++) {
      phase[i] = -phase[i];
      dphase[i] = -dphase[i];
      if (x) x[i] = -x[i];
    }
  }
  for (i = 1, k = i0+2; i < j; i++, k += 2) {
    h = phase[i] - phase[i-1];
    z[k] = 0.0;
    if (x != NULL) z[k] += x[i]*sin(phase[i]);
    if (y != NULL) z[k] += y[i]*cos(phase[i]);
    z[k] *= potential->dr_drho[k];
    if (i < 4) {
      if (h > 0.8) break;
    } else {
      if (h > 0.4) break;
    }
  }
  if (i > 1) {
    z[i0] = 0.0;
    if (x != NULL) z[i0] += x[0]*sin(phase[0]);
    if (y != NULL) z[i0] += y[0]*cos(phase[0]);
    z[i0] *= potential->dr_drho[i0];
    if (i == j) {
      i1 = i;
    } else {
      i1 = i + 1;
    }
    u1 = u + i1;
    w1 = w + i1;
    for (m = 0, n = i0; m < i; m++, n += 2) {
      u[m] = potential->rad[n];
      w[m] = potential->rad[n+1];
      z[n+1] = 0.0;
    }
    if (i1 > i) {
      u[i] = potential->rad[n];
    }
    UVIP3P(3, i1, u, phase, i, w, w1);
    if (x) {
      UVIP3P(3, i1, u, x, i, w, u1);
      for (m = 0, n = i0+1; m < i; m++, n += 2) {
        z[n] += u1[m]*sin(w1[m]);
      }
    }
    if (y) {
      UVIP3P(3, i1, u, y, i, w, u1);
      for (m = 0, n = i0+1; m < i; m++, n += 2) {
        z[n] += u1[m]*cos(w1[m]);
      }
    }
    for (m = 0, n = i0+1; m < i; m++, n += 2) {
      z[n] *= potential->dr_drho[n];
    }
    NewtonCotes(r, z, i0, k-2, t, 0);
  }

  q = i-1;
  m = j-q;
  if (m < 2) {
    return 0;
  }

  for (n = 1; n <= 5; n++) {
    i1 = q-n;
    if (i1 < 0 || phase[i1] >= phase[i1+1]) {
      i1++;
      break;
    }
  }
  if (t < 0) {
    for (n = i1, s = i0; n < j; n++, s += 2) {
      u[n] = potential->rad[s];
      z[n] = potential->rad[s+1];
    }
    UVIP3P(3, j-i1, u+i1, phase+i1, m, z+q, w+q);
  }

  if (x != NULL) {
    for (n = i1; n < j; n++) {
      x[n] /= dphase[n];
    }
    UVIP3C(3, j-i1, phase+i1, x+i1, z+i1, z+j+i1, x+j+i1);
  }
  if (y != NULL) {
    for (n = i1; n < j; n++) {
      y[n] /= dphase[n];
    }
    UVIP3C(3, j-i1, phase+i1, y+i1, u+i1, u+j+i1, y+j+i1);
  }

  si0 = sin(phase[i-1]);
  cs0 = cos(phase[i-1]);
  for (; i < j; i++, k += 2) {
    if (t < 0) {
      dr = w[i-1] - phase[i-1];
      si1 = sin(w[i-1]);
      cs1 = cos(w[i-1]);
      his0 = -(cs1 - cs0);
      hic0 = si1 - si0;
      p = dr;
      his1 = -p * cs1 + hic0;
      hic1 = p * si1 - his0;
      p *= dr; 
      his2 = -p * cs1 + 2.0*hic1;
      hic2 = p * si1 - 2.0*his1;
      p *= dr;
      his3 = -p * cs1 + 3.0*hic2;
      hic3 = p * si1 - 3.0*his2;
      r[k-1] = r[k-2];
    }
    d = phase[i] - phase[i-1];
    si1 = sin(phase[i]);
    cs1 = cos(phase[i]);
    is0 = -(cs1 - cs0);
    ic0 = si1 - si0;
    p = d;
    is1 = -p * cs1 + ic0;
    ic1 = p * si1 - is0;
    p *= d;
    is2 = -p * cs1 + 2.0*ic1;
    ic2 = p * si1 - 2.0*is1; 
    p *= d;
    is3 = -p * cs1 + 3.0*ic2;
    ic3 = p * si1 - 3.0*is2;
    r[k] = r[k-2];
    if (x != NULL) {
      a0 = x[i-1];
      a1 = z[i-1]; 
      a2 = z[j+i-1];
      a3 = x[j+i-1];
      h = a0*is0 + a1*is1 + a2*is2 + a3*is3;
      r[k] += h;
      if (t < 0) {
        h = a0*his0 + a1*his1 + a2*his2 + a3*his3;
        r[k-1] += h;
      }
    }
    if (y != NULL) {
      a0 = y[i-1];
      a1 = u[i-1];
      a2 = u[j+i-1];
      a3 = y[j+i-1];
      h = a0*ic0 + a1*ic1 + a2*ic2 + a3*ic3;
      r[k] += h;
      if (t < 0) {
        h = a0*hic0 + a1*hic1 + a2*hic2 + a3*hic3;
        r[k-1] += h;
      }
    }
    si0 = si1;
    cs0 = cs1;
  }

  return 0;
}

void LimitArrayRadial(int m, double n) {
  int i;
  
  if (m < 100) {
    n *= 1e6;
  }
  switch (m) {
  case -100:
    if (n < 0) n = ARYCTH;
    yk_array->cth = n;
    slater_array->cth = n;
    breit_array->cth = n;
    wbreit_array->cth = n;
    gos_array->cth = n;
    moments_array->cth = n;
    multipole_array->cth = n;
    residual_array->cth = n;
    vinti_array->cth = n;
    for (i = 0; i <= 4; i++) {
      xbreit_array[i]->cth = n;
    }
    break;
  case -1:
    LimitMultiSize(NULL, n);
    break;
  case 0:
    LimitMultiSize(yk_array, n);
    break;
  case 100:
    yk_array->cth = n;
    break;
  case 1:
    LimitMultiSize(slater_array, n);
    break;
  case 101:
    slater_array->cth = n;
    break;
  case 2:
    LimitMultiSize(breit_array, n);
    break;
  case 102:
    breit_array->cth = n;
    break;
  case 3:
    LimitMultiSize(wbreit_array, n);
    break;
  case 103:
    wbreit_array->cth = n;
    break;
  case 4:
    LimitMultiSize(gos_array, n);
    break;
  case 104:
    gos_array->cth = n;
    break;
  case 5:
    LimitMultiSize(moments_array, n);
    break;
  case 105:
    moments_array->cth = n;
    break;
  case 6:
    LimitMultiSize(multipole_array, n);
    break;
  case 106:
    multipole_array->cth = n;
    break;
  case 7:
    LimitMultiSize(residual_array, n);
    break;
  case 107:
    residual_array->cth = n;
    break;
  case 8:
    LimitMultiSize(vinti_array, n);
    break;
  case 108:
    vinti_array->cth = n;
    break;
  case 9:
    LimitMultiSize(qed1e_array, n);
    break;
  case 109:
    qed1e_array->cth = n;
    break;
  case 20:
  case 21:
  case 22:
  case 23:
  case 24:
    LimitMultiSize(xbreit_array[m-20], n);    
    break;
  case 120:
  case 121:
  case 122:
  case 123:
  case 124:
    xbreit_array[m-120]->cth = n;    
    break;
  default:
    printf("nothing is done\n");
    break;
  }
}

int InitRadial(void) {
  int ndim, i;
  int blocks[5] = {MULTI_BLOCK6,MULTI_BLOCK6,MULTI_BLOCK6,
		   MULTI_BLOCK6,MULTI_BLOCK6};
  potential = malloc(sizeof(POTENTIAL));
  hpotential = malloc(sizeof(POTENTIAL));
  rpotential = malloc(sizeof(POTENTIAL));
  potential->nfrozen = 0;
  potential->npseudo = 0;
  potential->mpseudo = 0;
  potential->dpseudo = 1;
  potential->mode = POTMODE;
  potential->hxs = POTHXS;
  potential->ihx = POTIHX;
  potential->hx0 = POTHX0;
  potential->hx1 = POTHX1;
  potential->hlike = 0;
  potential->nse = qed.se;
  if (qed.mse >= 1000) {
    potential->hpvs = 1;
    qed.mse = qed.mse%1000;
  } else {
    potential->hpvs = 0;
  }
  potential->pse = qed.mse >= 140;
  qed.mse = qed.mse%100;
  potential->mse = qed.mse;
  potential->mvp = qed.vp%100;
  potential->pvp = qed.vp > 100;  
  potential->flag = 0;
  potential->rb = 0;
  potential->atom = GetAtomicNucleus();
  SetBoundaryMaster(0, 1.0, -1.0);
  n_orbitals = 0;
  n_continua = 0;

  double cth = ARYCTH;
  orbitals = malloc(sizeof(ARRAY));
  if (!orbitals) return -1;
  if (ArrayInit(orbitals, sizeof(ORBITAL), ORBITALS_BLOCK) < 0) return -1;
  ndim = 5;
  slater_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(slater_array, sizeof(double), ndim, blocks, "slater_array");
  slater_array->cth = cth;
  
  ndim = 5;
  for (i = 0; i < ndim; i++) blocks[i] = MULTI_BLOCK5;
  breit_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(breit_array, sizeof(double *), ndim, blocks, "breit_array");
  breit_array->cth = cth;
  
  wbreit_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(wbreit_array, sizeof(double), ndim, blocks, "wbreit_array");
  wbreit_array->cth = cth;

  ndim = 3;
  for (i = 0; i < ndim; i++) blocks[i] = MULTI_BLOCK3;
  for (i = 0; i < 5; i++) {
    xbreit_array[i] = (MULTI *) malloc(sizeof(MULTI));
    char id[MULTI_IDLEN];
    sprintf(id, "xbreit_array%d", i);
    MultiInit(xbreit_array[i], sizeof(FLTARY), ndim, blocks, id);
    xbreit_array[i]->cth = cth;
  }
  
  ndim = 2;
  for (i = 0; i < ndim; i++) blocks[i] = MULTI_BLOCK2;
  residual_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(residual_array, sizeof(double), ndim, blocks, "residual_array");

  vinti_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(vinti_array, sizeof(double *), ndim, blocks, "vinti_array");

  qed1e_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(qed1e_array, sizeof(double), ndim, blocks, "qed1e_array");
  
  ndim = 3;
  for (i = 0; i < ndim; i++) blocks[i] = MULTI_BLOCK3;
  multipole_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(multipole_array, sizeof(double *), ndim, blocks, "multipole_array");
  multipole_array->cth = cth;

  ndim = 3;
  for (i = 0; i < ndim; i++) blocks[i] = MULTI_BLOCK3;
  moments_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(moments_array, sizeof(double), ndim, blocks, "moments_array"); 

  gos_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(gos_array, sizeof(double *), ndim, blocks, "gos_array");
  gos_array->cth = cth;

  yk_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(yk_array, sizeof(FLTARY), ndim, blocks, "yk_array");
  yk_array->cth = cth;

  n_awgrid = 1;
  awgrid[0]= EPS3;
  
  SetRadialGrid(DMAXRP, -1.0, -1.0, -1.0, -1.0);
  SetSlaterCut(-1, -1);

  SetOrbMap(0, 0, 0, 0);

  return 0;
}

void SetRadialCleanFlags(void) {  
  int i;
  SetMultiCleanFlag(slater_array);
  SetMultiCleanFlag(yk_array);
  SetMultiCleanFlag(breit_array);
  for (i = 0; i < 5; i++) {
    SetMultiCleanFlag(xbreit_array[i]);
  }
  SetMultiCleanFlag(wbreit_array);
  ReportMultiStats();
}

int ReinitRadial(int m) {
  if (m < 0) return 0;
#pragma omp barrier
#pragma omp master
  SetSlaterCut(-1, -1);
  ClearOrbitalTable(m);
  FreeSimpleArray(slater_array);
  FreeBreitArray();
  FreeSimpleArray(residual_array);
  FreeSimpleArray(qed1e_array);
  FreeMultipoleArray();
  FreeMomentsArray();
  FreeVintiArray();
  FreeYkArray();
  if (m < 2) {
    FreeGOSArray();
    if (m == 0) {
      if (optimize_control.n_screen > 0) {
	free(optimize_control.screened_n);
	optimize_control.n_screen = 0;
      }
      potential->flag = 0;
      n_awgrid = 1;
      awgrid[0] = EPS3;
      SetRadialGrid(DMAXRP, -1.0, -1.0, -1.0, -1.0);
    }
  }
#pragma omp barrier
  return 0;
}

void ElectronDensity(char *ofn, int n, int *ilev, int t) {
  FILE *f;
  ANGULAR_ZMIX *ang;
  LEVEL *lev;
  ORBITAL *orb0, *orb1;
  int i, j, k, m, nz, iz;
  double a, *p0, *q0, *p1, *q1;

  f = fopen(ofn, "w");
  if (f == NULL) {
    printf("cannot open file %s\n", ofn);
    return;
  }

  for (i = 0; i < n; i++) {
    nz = AngularZMix(&ang, ilev[i], ilev[i], 0, 0, NULL, NULL);
    if (nz <= 0) {
      continue;
    }
    for (k = 0; k < potential->maxrp; k++) {
      _dwork15[k] = 0;
      for (iz = 0; iz < nz; iz++) {
	if (ang[iz].k != 0) continue;
	orb0 = GetOrbitalSolved(ang[iz].k0);
	orb1 = GetOrbitalSolved(ang[iz].k1);
	if (k > orb0->ilast || k > orb1->ilast) continue;
	j = GetJFromKappa(orb0->kappa);
	p0 = Large(orb0);
	q0 = Small(orb0);
	p1 = Large(orb1);
	q1 = Small(orb1);
	switch (t) {
	case 1:
	  a = p0[k]*p1[k] + q0[k]*q1[k];	  
	  break;
	case 2:
	  a = p0[k]*p1[k];
	  break;
	case 3:
	  a = q0[k]*q1[k];
	  break;
	case 4:
	  a = p0[k]*q1[k] + q0[k]*p1[k];
	  break;
	case 5:
	  a = p0[k]*q1[k] - q0[k]*p1[k];	  
	  break;
	case 6:
	  a = p0[k]*q1[k];
	  break;
	default:
	  a = 0.0;
	  break;
	}
	a *= ang[iz].coeff*sqrt(j+1.0);
	_dwork15[k] += a;
      }
    }
    if (nz > 0) free(ang);
    lev = GetLevel(ilev[i]);
    DecodePJ(lev->pj, NULL, &j);
    a = 1.0/sqrt(j+1.0);
    for (m = potential->maxrp-1; m >= 0; m--) {
      if (_dwork15[m]) break;
    }
    fprintf(f, "# %4d %4d %4d %d\n", i, n, m, t);
    for (k = 0; k <= m; k++) {
      fprintf(f, "%6d %4d %15.8E %15.8E\n",
	      i, k, potential->rad[k], _dwork15[k]*a);
    }
  }
  fclose(f);
}

void ExpectationValue(char *ifn, char *ofn, int n, int *ilev,
		      double a, int t0) {
  int ni, i, np, k, iz, nz, j, t;
  double *ri, *vi, r, v;
  FILE *f0, *f1;
  char buf[1024], *p;
  ORBITAL *orb0, *orb1;
  ANGULAR_ZMIX *ang;
  LEVEL *lev;

  t = abs(t0);
  f0 = fopen(ifn, "r");
  if (f0 == NULL) {
    printf("cannot open file %s\n", ifn);
    return;
  }
  ni = 0;
  while (1) {
    if (fgets(buf, 1024, f0) == NULL) break;
    p = buf;
    while(*p == ' ') p++;
    if (*p == '#') continue;
    k = sscanf(p, "%lg %lg", &r, &v);
    if (k != 2) continue;
    ni++;
  }  
  if (ni < 2) {
    printf("too few points in ExpectationValue input: %d\n", ni);
    fclose(f0);
    return;
  }
  ri = malloc(sizeof(double)*ni);
  vi = malloc(sizeof(double)*ni);
  fseek(f0, 0, SEEK_SET);
  i = 0;
  while (1) {
    if (fgets(buf, 1024, f0) == NULL) break;
    p = buf;
    while(*p == ' ') p++;
    if (*p == '#') continue;
    k = sscanf(buf, "%lg %lg", &r, &v);
    if (k != 2) continue;
    ri[i] = potential->ar*pow(r, potential->qr) + potential->br*log(r);
    vi[i] = v;
    i++;
  }
  fclose(f0);
  np = 3;
  for (i = 0; i < potential->maxrp; i++) {
    _dwork14[i] = potential->ar*pow(potential->rad[i], potential->qr) +
      potential->br*log(potential->rad[i]);
  }
  UVIP3P(np, ni, ri, vi, potential->maxrp, _dwork14, _dwork15);
  if (t0 < 0) {
    for (i = 0; i < potential->maxrp; i++) {
      _dwork15[i] = exp(_dwork15[i]);
    }
  }
  if (a) {
    for (i = 0; i < potential->maxrp; i++) {
      _dwork15[i] *= pow(potential->rad[i], a);
    }
  }
  f1 = fopen(ofn, "w");
  if (f1 == NULL) {
    printf("cannot open file: %s\n", ofn);
    free(ri);
    free(vi);
    return;
  }
  fprintf(f1, "# %4d %g %d\n", n, a, t);
  for (i = 0; i < n; i++) {
    nz = AngularZMix(&ang, ilev[i], ilev[i], 0, 0, NULL, NULL);
    if (nz <= 0) {
      continue;
    }
    v = 0;
    for (iz = 0; iz < nz; iz++) {
      if (ang[iz].k != 0) continue;
      orb0 = GetOrbitalSolved(ang[iz].k0);
      orb1 = GetOrbitalSolved(ang[iz].k1);
      Integrate(_dwork15, orb0, orb1, t, &r, 0);
      j = GetJFromKappa(orb0->kappa);
      v += r * ang[iz].coeff * sqrt(j+1.0);
    }
    lev = GetLevel(ilev[i]);
    DecodePJ(lev->pj, NULL, &j);
    v /= sqrt(j+1.0);
    fprintf(f1, "%6d %20.12E\n", ilev[i], v);
    free(ang);
  }  
  fclose(f1);  
}

int TestIntegrate(void) {
  int n, m;
  double r, r0;  

  m = 20;
  for (n = 0; n < potential->maxrp; n++) {
    _xk[n] = pow(potential->rad[n], m)*potential->dr_drho[n];
  }
  _yk[0] = 0.0;
  NewtonCotes(_yk, _xk, 0, potential->maxrp-1, -1, 0);
  for (n = 0; n < potential->maxrp; n++) {
    if (n > 10) break;
    if (n > 0 && n < potential->maxrp-1) {
      if (n%2 == 1) r = 0.5*(_yk[n-1]+_yk[n+1]);
      else r = _yk[n];
    }
    r0 = (pow(potential->rad[n],m+1)-pow(potential->rad[0],m+1))/(m+1.0);
    printf("%4d %12.5E %12.5E %12.5E %12.5E %12.5E\n",
	   n, potential->rad[n], _xk[n], _yk[n], r, r0);
  }

  return 0;
}

void RemoveOrbitalLock(void) {
  if (orbitals->lock) {
    DestroyLock(orbitals->lock);
    free(orbitals->lock);
    orbitals->lock = NULL;
  }
}

int TestIntegrate0(void) {
  ORBITAL *orb1, *orb2, *orb3, *orb4;
  int k1, k2, k3, k4, k, i, i0=1500;
  double r, a, s[6];
 
  k1 = OrbitalIndex(2, -2, 0);
  k2 = OrbitalIndex(3, -1, 0);
  k3 = OrbitalIndex(0, -5, 3.30265398e3/HARTREE_EV);
  k4 = OrbitalIndex(0, -4, 2.577e3/HARTREE_EV);
  printf("# %d %d %d %d\n", k1, k2, k3, k4);
  orb1 = GetOrbital(k1);
  orb2 = GetOrbital(k2);
  orb3 = GetOrbital(k3);
  orb4 = GetOrbital(k4);

  for (k = 1; k < 7; k++) {    
    for (i = 0; i < potential->maxrp; i++) {
      _yk[i] = pow(potential->rad[i],k);
      _yk[i+i0] = 1.0/_yk[i];
    }
    Integrate(_yk+i0, orb3, orb4, -1, _zk, 0);
    Integrate(_yk+i0, orb3, orb4, -1, _zk+i0, -1);
    for (i = 0; i < potential->maxrp; i++) {
      printf("%d %4d %15.8E %15.8E %15.8E %15.8E %15.8E\n", 
	     k, i, potential->rad[i], _yk[i], _zk[i], _yk[i+i0], _zk[i+i0]);
    }
    Integrate(_yk+i0, orb3, orb4, 1, &r, 0);
    Integrate(_yk+i0, orb3, orb4, 1, &a, -1);
    printf("# %15.8E %15.8E\n", r, a);
    printf("\n\n");
  }
  
  return 0;
}

void OptimizeModSE(int n, int ka, double dr, int ni) {
  double e0, eps, e, r0, r1, r, z, rs, rms;
  ORBITAL *orb;
  int k, i;

  z = potential->atom->atomic_number;
  if (z < 60) eps = 1e-4;
  else eps = 1e-5;
  r = potential->atom->rmse;
  rms = potential->atom->rms*1e5*RBOHR;
  if (r <= 0) {
    r = rms;
  }
  dr = r*dr;
  k = OrbitalIndex(n, ka, 0);
  orb = GetOrbitalSolved(k);
  e0 = HydrogenicSelfEnergy(31, 0, 1.0, potential, orb, orb);
  if (potential->N > 1) {
    if (orb->rorb == NULL) {
      orb->rorb = SolveAltOrbital(orb, rpotential);
    }
    if (orb->horb == NULL) {
      orb->horb = SolveAltOrbital(orb, hpotential);
    }
    rs = SelfEnergyRatio(orb->rorb, orb->horb);
    e0 *= rs;
  } else {
    rs = 1.0;
  }
  INIQED(z, 9, 1, r);
  e = HydrogenicSelfEnergy(51, 0, 1.0, potential, orb, orb);
  if (ni > 0 && fabs(e0-e) > eps) {
    i = 0;
    if (e < e0) {
      r0 = r;
      r1 = r0;
      while (i < ni && e < e0) {
	r1 = r1 + dr;
	INIQED(z, 9, 1, r1);
	e = HydrogenicSelfEnergy(51, 0, 1.0, potential, orb, orb);
	if (fabs(e0-e) < eps) {
	  r = r1;
	  break;
	}
	i++;
      }
    } else {
      r1 = r;
      r0 = r1;
      while (i < ni && e > e0) {
	r0 = r0 - dr;
	INIQED(z, 9, 1, r0);
	e = HydrogenicSelfEnergy(51, 0, 1.0, potential, orb, orb);
	if (fabs(e0-e) < eps) {
	  r = r0;
	  break;
	}
	i++;
      }
    }
    i = 0;
    while(i < ni) {
      r = 0.5*(r0 + r1);
      INIQED(z, 9, 1, r);
      e = HydrogenicSelfEnergy(51, 0, 1.0, potential, orb, orb);
      if (fabs(e-e0) < eps) break;
      if (r1-r0 < 1e-5) break;
      if (e < e0) {
	r0 = r;
      } else if (e > e0) {
	r1 = r;
      }
      i++;
    }
  }
  
  printf("modqed rms: %3.0f %18.10E %18.10E %18.10E %18.10E %18.10E %18.10E\n",
	 z, rms, r, rs, e, e0, e-e0);
}

int AddNewConfigToList(int k, int ni, int *kc,
		       CONFIG *c0, int nb, int **kcb,
		       int nc, SHELL_RESTRICTION *sr) {
  CONFIG *cfg = ConfigFromIList(ni, kc);
  int r;
  double sth;
  if (nc > 0) {
    r = ApplyRestriction(1, cfg, nc, sr);
    if (r <= 0) return -1;
  }
  if (ConfigExists(cfg)) return -1;
  sth = c0->sth;
  if (sth > 0) {
    int i0, i1, i2, i3;
    int n0, k0, n1, k1, n2, k2, n3, k3;
    int i, j;
    double s = 0;
    int sms0 = qed.sms;
    int br0 = qed.br;
    qed.sms = 0;
    qed.br = 0;
    for (i = 0; i < nb; i++) {
      i0 = -1;
      i1 = -1;
      i2 = -1;
      i3 = -1;
      for (j = 0; j < ni; j++) {
	if (kc[j] == kcb[i][j]+2) {
	  if (i1 >= 0 || i3 >= 0) break;	  
	  i1 = j;
	  i3 = j;
	} else if (kc[j] == kcb[i][j]-2) {
	  if (i0 >= 0 || i2 >= 0) break;
	  i0 = j;
	  i2 = j;
	} else if (kc[j] == kcb[i][j]+1) {
	  if (i1 < 0) i1 = j;
	  else if (i3 < 0) i3 = j;
	  else break;
	} else if (kc[j] == kcb[i][j]-1) {
	  if (i0 < 0) i0 = j;
	  else if (i2 < 0) i2 = j;
	  else break;
	} else if (abs(kc[j]-kcb[i][j]) > 2) {
	  break;
	}
      }
      if (j < ni) continue;
      
      IntToShell(i0, &n0, &k0);
      IntToShell(i1, &n1, &k1);
      int i2s0, i2s1, i2s;
      int i3s0, i3s1, i3s;
      if (i2 >= 0) {
	IntToShell(i2, &n2, &k2);
	i2s0 = 0;
	i2s1 = 1;
      } else {
	i2s0 = 0;
	i2s1 = cfg->n_shells;
	//n2 = cfg->shells[0].n;
	//k2 = cfg->shells[0].kappa;
      }
      if (i3 >= 0) {
	i3s0 = 0;
	i3s1 = 1;
	IntToShell(i3, &n3, &k3);
      } else {
	i3s0 = 0;
	i3s1 = cfg->n_shells;
	//n3 = cfg->shells[0].n;
	//k3 = cfg->shells[0].kappa;
      }
      int ko0, ko1, ko2, ko3;
      ORBITAL *o0, *o1, *o2, *o3;
      int j0, kl0, j1, kl1, j2, kl2, j3, kl3;
      ko0 = OrbitalIndex(n0, k0, 0);
      ko1 = OrbitalIndex(n1, k1, 0);
      GetJLFromKappa(k0, &j0, &kl0);
      GetJLFromKappa(k1, &j1, &kl1);
      o0 = GetOrbital(ko0);
      o1 = GetOrbital(ko1);
      for (i2s = i2s0; i2s < i2s1; i2s++) {
	if (i2 < 0) {
	  n2 = cfg->shells[i2s].n;
	  k2 = cfg->shells[i2s].kappa;
	  if (n2 == n0 && k2 == k0 && cfg->shells[i2s].nq < 2) continue;
	}
	ko2 = OrbitalIndex(n2, k2, 0);
	GetJLFromKappa(k2, &j2, &kl2);
	o2 = GetOrbital(ko2);
	for (i3s = i3s0; i3s < i3s1; i3s++) {
	  if (i3 < 0) {
	    if (i2 < 0 && i2s != i3s) continue;
	    n3 = cfg->shells[i3s].n;
	    k3 = cfg->shells[i3s].kappa;
	    if (n3 == n1 && k3 == k1 && cfg->shells[i3s].nq < 2) continue;
	  }
	  ko3 = OrbitalIndex(n3, k3, 0);
	  GetJLFromKappa(k3, &j3, &kl3);
	  if (IsOdd((kl0+kl1+kl2+kl3)/2)) continue;
	  o3 = GetOrbital(ko3);
	  double de = fabs(o1->energy-o0->energy + o3->energy-o2->energy);
	  if (de < EPS10) {
	    s = 1e10;
	    break;
	  } else {
	    double s2 = 0, s1 = 0, sd = 0, se = 0;
	    int kk, ks[4];
	    ks[0] = ko0;
	    ks[1] = ko2;
	    ks[2] = ko1;
	    ks[3] = ko3;
	    for (kk = 0; kk < 5; kk += 2) {
	      SlaterTotal(&sd, &se, NULL, ks, kk, 0);
	      sd = fabs(sd+se);
	      if (s2 < sd) s2 = sd;
	    }
	    if (i2 < 0 && i3 < 0 && k0 == k1) {
	      ResidualPotential(&s1, ko0, ko1);
	      s1 = fabs(s1);
	      if (s2 < s1) s2 = s1;
	    }
	    s2 /= de;
	    if (s < s2) s = s2;
	  }
	}
      }
    }
    qed.sms = sms0;
    qed.br = br0;
    if (s < sth) {
      return -10;
    }
  }
  if (Couple(cfg) < 0) return -1;
  r = AddConfigToList(k, cfg);
  free(cfg);
  return r;
}

int ConfigSD(int m0, int ng, int *kg, char *s, char *gn1, char *gn2,
	     int n0, int n1, int n0d, int n1d, int k0, int k1,
	     int ngb, int *kgb, double sth) {
  int ni, nr, *kc, nb, **kcb, i, j, k, ir, ns, ks, ks2, ka, is, js;
  int t, ird, nd, kd, kd2, jd, id, ig1, ig2;
  CONFIG_GROUP *g;
  CONFIG *c, *cr;
  SHELL_RESTRICTION *sr;
  int m, mar, *kcr, nc, *kcrn, nnr, nn, km, kt, km0;
  double sth0;

  if (s) {
    nc = GetRestriction(s, &sr, 0);
  } else {
    nc = 0;
    sr = NULL;
  }
  
  if (m0 == 0) {
    ig1 = GroupIndex(gn1);
    if (ig1 < 0) {
      printf("invalid config group name: %s\n", gn1);
      if (nc > 0) {
	for (i = 0; i < nc; i++) {
	  free(sr[i].shells);
	}
	free(sr);
      }
      return -1;
    }
    nr = GetConfigFromString(&cr, s);
    if (nr <= 0) {
      if (nc > 0) {
	for (i = 0; i < nc; i++) {
	  free(sr[i].shells);
	}
	free(sr);
      }
      g = GetGroup(ig1);
      if (g != NULL && g->n_cfgs == 0) RemoveGroup(ig1);
      return 0;
    }
    ni = 0;
    for (i = 0; i < nr; i++) {
      cr[i].sth = sth;
      if (ni < cr[i].shells[0].n) ni = cr[i].shells[0].n;
    }
    for (i = 0; i < ng; i++) {
      g = GetGroup(kg[i]);
      if (ni < g->nmax) ni = g->nmax;
    }
    ni = ni*ni;
    kc = malloc(sizeof(int)*ni);
    nb = 0;
    if (ngb > 0) {
      if (kgb) {
	for (i = 0; i < ngb; i++) {
	  g = GetGroup(kgb[i]);
	  nb += g->n_cfgs;
	}
	kcb = NULL;
	if (nb > 0) {
	  kcb = malloc(sizeof(int *)*nb);
	  for (i = 0; i < nb; i++) {
	    kcb[i] = malloc(sizeof(int)*ni);
	  }
	}
	k = 0;
	for (i = 0; i < ngb; i++) {
	  g = GetGroup(kgb[i]);
	  for (j = 0; j < g->n_cfgs; j++, k++) {
	    c = GetConfigFromGroup(kgb[i], j);
	    ConfigToIList(c, ni, kcb[k]);
	  }
	}
      } 
    }
    for (i = 0; i < nr; i++) {
      ConfigToIList(&cr[i], ni, kc);
      AddNewConfigToList(ig1, ni, kc, &cr[i], nb, kcb, nc, sr);
    }
    
    free(kc);
    if (nb > 0) {
      for (i = 0; i < nb; i++) {
	free(kcb[i]);
      }
      free(kcb);
    }
    if (nc > 0) {
      for (i = 0; i < nc; i++) {
	free(sr[i].shells);
      }
      free(sr);
    }
    if (sth > 0) ReinitRadial(2);
    g = GetGroup(ig1);
    if (g != NULL && g->n_cfgs == 0) RemoveGroup(ig1);
    return 0;
  }
  
  if (n0d <= 0) n0d = n0;
  if (n1d <= 0) n1d = n1;
  if (k0 < 0) k0 = 0;
  if (k1 < 0) k1 = Max(n1,n1d)-1;
  if (gn2 == NULL || strlen(gn2) == 0) gn2 = gn1;
  m = m0%10;
  mar = m0/10;
  if (m < 1 || m > 3 || mar > 2) {
    printf("invalid mode: %d %d %d\n", m0, m, mar);
    return -1;
  }
  if (mar > 0) sth = 0;
  if (mar == 2) {
    nr = 1;
    cr = NULL;
  } else {
    nr = GetConfigFromString(&cr, s);
  }
  if (mar == 1) {
    n0 = 1;
    n1 = 1;
    k0 = 0;
    k1 = 0;
    n0d = 1;
    n1d = 1;    
  }
  if (nr <= 0) {
    printf("invalid reference shell spec in ConfigSD: %s\n", s);
    return -1;
  }
  ni = Max(n1, n1d);
  for (i = 0; i < ng; i++) {
    g = GetGroup(kg[i]);
    if (ni < g->nmax) ni = g->nmax;
  }
  if (kgb && kgb != kg) {
    for (i = 0; i < ngb; i++) {
      g = GetGroup(kgb[i]);
      if (ni < g->nmax) ni = g->nmax;
    }
  }
  ni = ni*ni;
  kc = malloc(sizeof(int)*ni);
  kcr = NULL;  
  if (cr) {
    for (i = 0; i < ni; i++) kc[i] = 0;
    for (i = 0; i < nr; i++) {
      for (j = 0; j < cr[i].n_shells; j++) {
	k = ShellToInt(cr[i].shells[j].n, cr[i].shells[j].kappa);
	if (k >= 0) kc[k] = 1;
      }
      free(cr[i].shells);
    }  
    free(cr);
    nr = 0;
    for (i = 0; i < ni; i++) {
      if (kc[i]) nr++;
    }
    kcr = malloc(sizeof(int)*nr);
    j = 0;
    for (i = 0; i < ni; i++) {
      if (kc[i]) {
	kcr[j] = i;
	j++;
      }
    }
  }
  ig1 = GroupIndex(gn1);
  if (ig1 < 0) {
    printf("invalid config group name: %s\n", gn1);
    return -1;
  }
  ig2 = GroupIndex(gn2);
  if (ig2 < 0) {
    printf("invalid config group name: %s\n", gn2);
    return -1;
  }
  nb = 0;
  if (ngb > 0) {
    if (kgb && sth >= 0) {
      for (i = 0; i < ngb; i++) {
	g = GetGroup(kgb[i]);
	nb += g->n_cfgs;
      }
      if (nb > 0) {
	kcb = malloc(sizeof(int *)*nb);
	for (i = 0; i < nb; i++) {
	  kcb[i] = malloc(sizeof(int)*ni);
	}
      }
      k = 0;
      for (i = 0; i < ngb; i++) {
	g = GetGroup(kgb[i]);
	for (j = 0; j < g->n_cfgs; j++, k++) {
	  c = GetConfigFromGroup(kgb[i], j);
	  ConfigToIList(c, ni, kcb[k]);
	}
      }
    } else {
      nb = 1;
      kcb = malloc(sizeof(int *)*nb);
      kcb[0] = malloc(sizeof(int)*ni);
    }
  }
  nnr = nr;
  kcrn = malloc(sizeof(int)*(nnr+1));
  if (kcr) {
    km0 = -1;
    nnr = 0;
    for (k = 0; k < nr; k++) {
      ir = kcr[k];
      IntToShell(ir, &nn, &km);
      km = GetLFromKappa(km);
      if (km != km0) {
	kcrn[nnr] = k;
	km0 = km;
	nnr++;
      }
    }
    kcrn[nnr] = nr;
  } else {
    kcrn[0] = 0;
    kcrn[1] = 1;
  }
  if (ng > 0 && kg) {
    for (i = 0; i < ng; i++) {
      g = GetGroup(kg[i]);
      for (j = 0; j < g->n_cfgs; j++) {
	c = GetConfigFromGroup(kg[i], j);
	c->sth = sth;
      }
    }
    if (sth < -EPS10) {
      if (ngb <= 0 || !kgb) {
	printf("sth < 0 without kgb: %g %d %d\n", sth, ngb, ng);
	return -1;
      }
      GetInteractConfigs(ngb, kgb, ng, kg, -sth);
    }
  }
  if (m != 2) {
    for (i = 0; i < ng; i++) {
      g = GetGroup(kg[i]);
      for (j = 0; j < g->n_cfgs; j++) {
	c = GetConfigFromGroup(kg[i], j);
	if (c->sth < -EPS10) continue;
	ConfigToIList(c, ni, kc);
	if ((sth < 0 || kgb == NULL) && nb == 1) {
	  memcpy(kcb[0], kc, sizeof(int)*ni);
	}
	for (km = 0; km < nnr; km++) {	
	  for (ns = n0; ns <= n1; ns++) {
	    for (ks = k0; ks <= k1; ks++) {
	      if (ks >= ns) break;
	      int pr = 1000000;
	      for (k = kcrn[km]; k < kcrn[km+1]; k++) {
		if (kcr) {
		  ir = kcr[k];
		  if (kc[ir] <= 0) continue;
		} else {
		  ir = -1;
		}
		ks2 = 2*ks;
		for (js = ks2-1; js <= ks2+1; js += 2) {
		  if (js < 0) continue;
		  if (mar != 1) {
		    ka = GetKappaFromJL(js, ks2);
		    is = ShellToInt(ns, ka);
		    if (kc[is] == js+1) continue;
		  } else {
		    is = -1;
		  }
		  if (ir >= 0) kc[ir]--;
		  if (is >= 0) kc[is]++;
		  if (pr == 1000000) {
		    pr = AddNewConfigToList(ig1, ni, kc, c, nb, kcb, nc, sr);
		  } else if (pr > -10) {
		    double sth0 = c->sth;
		    c->sth = 0.0;
		    AddNewConfigToList(ig1, ni, kc, c, nb, kcb, nc, sr);
		    c->sth = sth0;
		  }
		  if (ir >= 0) kc[ir]++;
		  if (is >= 0) kc[is]--;
		}
	      }
	    }
	  }
	}
      }
    }
  }
  if (m != 1) {
    for (i = 0; i < ng; i++) {
      g = GetGroup(kg[i]);
      for (j = 0; j < g->n_cfgs; j++) {
	c = GetConfigFromGroup(kg[i], j);
	if (c->sth < -EPS10) continue;
	ConfigToIList(c, ni, kc);
	if ((sth < 0 || kgb == NULL) && nb == 1) {
	  memcpy(kcb[0], kc, sizeof(int)*ni);
	}
	for (km = 0; km < nnr; km++) {	
	  for (ns = n0; ns <= n1; ns++) {
	    for (kt = 0; kt < nnr; kt++) {
	      for (nd = ns; nd <= n1d; nd++) {
		if (nd < n0d) continue;
		for (ks = k0; ks <= k1; ks++) {
		  if (ks >= ns) break;
		  for (kd = k0; kd <= k1; kd++) {
		    if (kd >= nd) break;
		    int pr = 1000000;
		    ks2 = 2*ks;
		    for (js = ks2-1; js <= ks2+1; js += 2) {
		      if (js < 0) continue;
		      if (mar != 1) {
			ka = GetKappaFromJL(js, ks2);
			is = ShellToInt(ns, ka);
			if (kc[is] == js+1) continue;
		      } else {
			is = -1;
		      }	
		      kd2 = 2*kd;
		      for (jd = kd2-1; jd <= kd2+1; jd += 2) {
			if (jd < 0) continue;	  
			for (k = kcrn[km]; k < kcrn[km+1]; k++) {
			  if (kcr) {
			    ir = kcr[k];
			    if (kc[ir] <= 0) continue;
			  } else {
			    ir = -1;
			  }
			  for (t = kcrn[kt]; t < kcrn[kt+1]; t++) {
			    if (kcr) {
			      ird = kcr[t];
			      if (kc[ird] <= 0) continue;
			    } else {
			      ird = -1;
			    }
			    if (mar != 1) {
			      ka = GetKappaFromJL(jd, kd2);
			      id = ShellToInt(nd, ka);
			      if (kc[id] == jd+1) continue;
			    } else {
			      id = -1;
			    }
			    if (ir >= 0) kc[ir]--;			      
			    if (is >= 0) kc[is]++;
			    if (ird >= 0) kc[ird]--;
			    if (id >= 0) kc[id]++;
			    if (kc[ir] >= 0 && kc[ird] >= 0 &&
				kc[is] <= js+1 && kc[id] <= jd+1) {
			      if (pr == 1000000) {
				pr = AddNewConfigToList(ig2, ni, kc, c,
							nb, kcb, nc, sr);
			      } else if (pr > -10) {
				double sth0 = c->sth;
				c->sth = 0.0;
				AddNewConfigToList(ig2, ni, kc, c,
						   nb, kcb, nc, sr);
				c->sth = sth0;
			      }
			    }
			    if (ird >= 0) kc[ird]++;
			    if (id >= 0) kc[id]--;
			    if (ir >= 0) kc[ir]++;
			    if (is >= 0) kc[is]--;
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  free(kc);
  if (kcr) free(kcr);
  free(kcrn);
  if (nc > 0) {
    for (i = 0; i < nc; i++) {
      free(sr[i].shells);
    }
    free(sr);
  }
  if (nb > 0) {
    for (i = 0; i < nb; i++) {
      free(kcb[i]);
    }
    free(kcb);
  }
  if (sth > 0) ReinitRadial(2);
  g = GetGroup(ig2);
  if (g != NULL && g->n_cfgs == 0) RemoveGroup(ig2);
  if (ig1 != ig2) {
    g = GetGroup(ig1);
    if (g != NULL && g->n_cfgs == 0) RemoveGroup(ig1);
  }
  return 0;
}

