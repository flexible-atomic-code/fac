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

#include "rmatrix.h"
#include "cf77.h"

static char *rcsid="$Id: rmatrix.c,v 1.8 2004/12/23 21:31:03 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

static RBASIS rbasis;
static int ntg, *tg, nts, *ts;
static int ncg, *cg, ncs, *cs;
static DCFG dcfg, dcfg0;
static int _nrefine=4;
static int _refiter=1;
static int _mrefine=16;
static double _mineref=1e-8;
static double _rrefine=0.25;
static int fmode;
static int _rmx_dk = 3;
static int _rmx_acs = 1;
static char _rmx_bfn[256] = "";
static char _rmx_efn[256] = "";
static int _rmx_pj = -1;
static int _rmx_mpj = -1;
static int _stark_amp = 0;
static int _stark_nts = 0;
static int _stark_pw = 0;
static int *_stark_lower = NULL;
static int *_stark_upper = NULL;
static IDXARY _stark_idx;
static int _stark_eadj = 0;
static int _gailitis_expni = 10;
static double _gailitis_exprf = 2.0;
static double _gailitis_exprt = 0.25;
static int _rmx_isave = 1;
static double _rmx_fmaxe = 1.5;
static int _minsp = 3;

#pragma omp threadprivate(dcfg)

int SkipSym(int pj) {
  int i0, i1;
  if (_rmx_pj < 0 && _rmx_mpj < 0) return 0;
  i0 = 0;
  if (_rmx_pj >= 0) i0 = _rmx_pj;
  i1 = i0;
  if (_rmx_mpj >= 0) i1 = _rmx_mpj;
  else if (_rmx_mpj == -2) i1 = MAX_SYMMETRIES;
  if (pj < i0) return 1;
  if (pj > i1) return 1;
  return 0;
}

void SetStarkIDs(char *ids) {
  char buf[8192];
  
  if (_stark_nts > 0) {
    free(_stark_lower);
    free(_stark_upper);
    _stark_nts = 0;
  }
  FreeIdxAry(&_stark_idx, 0);

  if (ids == NULL || strlen(ids) == 0) return;
  strncpy(buf, ids, 8192);
  int n = StrSplit(buf, ',');
  if (n <= 0) return;

  int *ilo = malloc(sizeof(int)*n);
  int *iup = malloc(sizeof(int)*n);
  char *p, *r;
  int m, k, u, i;
  p = buf;
  k = 0;
  m = 0;
  for (i = 0; i < n; i++) {
    while (*p == ' ') p++;
    r = p;
    while (*r && *r != '-') r++;
    if (*r == '-') {
      r++;
      ilo[k] = atoi(p);
      iup[k] = atoi(r);
      if (ilo[k] > m) m = ilo[k];
      if (iup[k] > m) m = iup[k];
      k++;
    }
    while (*p) p++;
    p++;
  }
  m++;
  int *id = malloc(sizeof(int)*(m+1));
  for (i = 0; i < m; i++) id[i] = 0;
  for (i = 0; i < k; i++) {
    id[ilo[i]] = 1;
    id[iup[i]] = 1;
  }
  u = 0;
  for (i = 0; i < m; i++) {
    if (id[i]) {
      id[u++] = i;
    }
  }
  
  InitIdxAry(&_stark_idx, u, id);
  _stark_lower = ilo;
  _stark_upper = iup;
  _stark_nts = k;
}

void RMatrixRefine(int n, int m, double r) {
  if (n >= 0) _nrefine = n;
  if (m > 0) _mrefine = m;
  if (r > 0) _rrefine = r;
}

void RMatrixFMode(int m) {
  fmode = m;
}

int InitRMatrix(void) {
  ntg = 0;
  ncg = 0;
  nts = 0;
  ncs = 0;
  rbasis.nkappa = 0;
  dcfg.nmultipoles = 2;
  dcfg.ngailitis = 0;
  dcfg.rgailitis = 0.0;
  dcfg.degenerate = EPS4;
  dcfg.accuracy = EPS4;
  dcfg.pdirection = 1;

  fmode = 0;
  InitIdxAry(&_stark_idx, 0, NULL);
  _stark_nts = 0;
  return 0;
}

void ClearRMatrixBasis(RBASIS *rbs) {
  int i, k;

  for (i = 0; i < rbs->nkappa; i++) {
    for (k = 0; k < NBTERMS; k++) {
      free(rbs->cbuttle[k][i]);
    }
    free(rbs->ebuttle[i]);
    free(rbs->basis[i]);
    free(rbs->bnode[i]);
    free(rbs->w0[i]);
    free(rbs->w1[i]);
    free(rbs->ek[i]);
  }
  if (rbs->nkappa > 0) {
    for (k = 0; k < NBTERMS; k++) {
      free(rbs->cbuttle[k]);
    }
    free(rbs->kappa);
    free(rbs->ebuttle);
    free(rbs->basis);
    free(rbs->bnode);
    free(rbs->w0);
    free(rbs->w1);
    free(rbs->ek);
  }
  rbs->nkappa = 0;
}

void RMatrixExpansion(int n, double d, double a, double r) {
  if (n < 0) n = 0;
  dcfg.ngailitis = n;
  if (d < 1E-4) d = 1E-4;
  dcfg.degenerate = d;
  dcfg.accuracy = a;
  dcfg.rgailitis = r;
}

void PropogateDirection(int m) {
  dcfg.pdirection = m;
}

void RMatrixNMultipoles(int n) {
  dcfg.nmultipoles = n;
}

void ReadRMatrixBasis(char *fn, RBASIS *rbs, int fmt) {
  FILE *f;
  int i, j, k, n, kappa, m, nr;

  f = fopen(fn, "r");
  if (fmt == 0) {
    nr = fread(&(rbs->ib0), sizeof(int), 1, f);
    nr = fread(&(rbs->rb0), sizeof(double), 1, f);
    nr = fread(&(rbs->ib1), sizeof(int), 1, f);
    nr = fread(&(rbs->rb1), sizeof(double), 1, f);
    nr = fread(&(rbs->bqp), sizeof(double), 1, f);
    nr = fread(&(rbs->emin), sizeof(double), 1, f);
    nr = fread(&(rbs->kmax), sizeof(int), 1, f);
    rbs->ki = rbs->kmax/1000000;
    rbs->kmax = rbs->kmax%1000000;
    rbs->kmin = rbs->kmax/1000;
    rbs->kmax = rbs->kmax%1000;
    rbs->kapmin = 0;
    if (rbs->kmin > 0) {
      rbs->kapmin = 2*(rbs->kmin-1) + (rbs->kmin>1?2:1);
    }
    nr = fread(&(rbs->nbk), sizeof(int), 1, f);
    rbs->nbi = rbs->nbk/1000;
    rbs->nbk = rbs->nbk%1000;
    nr = fread(&(rbs->nkappa), sizeof(int), 1, f);
    nr = fread(&(rbs->nbuttle), sizeof(int), 1, f);    
    rbs->kappa = malloc(sizeof(int)*rbs->nkappa);
    nr = fread(rbs->kappa, sizeof(int), rbs->nkappa, f);
    rbs->basis = malloc(sizeof(int *)*rbs->nkappa);
    rbs->bnode = malloc(sizeof(int *)*rbs->nkappa);
    rbs->ebuttle = malloc(sizeof(double *)*rbs->nkappa);
    for (k = 0; k < NBTERMS; k++) {
      rbs->cbuttle[k] = malloc(sizeof(double *)*rbs->nkappa);
    }
    rbs->w0 = malloc(sizeof(double *)*rbs->nkappa);
    rbs->w1 = malloc(sizeof(double *)*rbs->nkappa);
    rbs->ek = malloc(sizeof(double *)*rbs->nkappa);
    for (i = 0; i < rbs->nkappa; i++) {
      rbs->basis[i] = malloc(sizeof(int)*rbs->nbk);
      rbs->bnode[i] = malloc(sizeof(int)*rbs->nbk);
      rbs->ebuttle[i] = malloc(sizeof(double)*rbs->nbuttle);
      for (k = 0; k < NBTERMS; k++) {
	rbs->cbuttle[k][i] = malloc(sizeof(double)*rbs->nbuttle);
      }
      rbs->w0[i] = malloc(sizeof(double)*rbs->nbk);    
      rbs->w1[i] = malloc(sizeof(double)*rbs->nbk);
      rbs->ek[i] = malloc(sizeof(double)*rbs->nbk);
    }
    for (i = 0; i < rbs->nkappa; i++) {
      nr = fread(rbs->basis[i], sizeof(int), rbs->nbk, f);
      nr = fread(rbs->bnode[i], sizeof(int), rbs->nbk, f);
      nr = fread(rbs->ek[i], sizeof(double), rbs->nbk, f);
      nr = fread(rbs->w0[i], sizeof(double), rbs->nbk, f);
      nr = fread(rbs->w1[i], sizeof(double), rbs->nbk, f);
      nr = fread(rbs->ebuttle[i], sizeof(double), rbs->nbuttle, f);
      for (k = 0; k < NBTERMS; k++) {
	nr = fread(rbs->cbuttle[k][i], sizeof(double), rbs->nbuttle, f);
      }
    }
  } else {
    fscanf(f, "%d %lf %d %lf %lf\n",
	   &(rbs->ib0), &(rbs->rb0), &(rbs->ib1), &(rbs->rb1),
	   &(rbs->bqp));
    fscanf(f, "%d %d %d %d\n",
	   &(rbs->kmax), &(rbs->nbk), &(rbs->nkappa), &(rbs->nbuttle));
    rbs->ki = rbs->kmax/1000000;
    rbs->kmax = rbs->kmax%1000000;
    rbs->kmin = rbs->kmax/1000;
    rbs->kmax = rbs->kmax%1000;
    rbs->nbi = rbs->nbk/1000;
    rbs->nbk = rbs->nbk%1000;    
    rbs->kappa = malloc(sizeof(int)*rbs->nkappa);
    rbs->basis = malloc(sizeof(int *)*rbs->nkappa);
    rbs->bnode = malloc(sizeof(int *)*rbs->nkappa);
    rbs->ebuttle = malloc(sizeof(double *)*rbs->nkappa);
    for (k = 0; k < NBTERMS; k++) {
      rbs->cbuttle[k] = malloc(sizeof(double *)*rbs->nkappa);
    }
    rbs->w0 = malloc(sizeof(double *)*rbs->nkappa);
    rbs->w1 = malloc(sizeof(double *)*rbs->nkappa);
    rbs->ek = malloc(sizeof(double *)*rbs->nkappa);
    for (i = 0; i < rbs->nkappa; i++) {
      rbs->basis[i] = malloc(sizeof(int)*rbs->nbk);
      rbs->bnode[i] = malloc(sizeof(int)*rbs->nbk);
      rbs->ebuttle[i] = malloc(sizeof(double)*rbs->nbuttle);
      for (k = 0; k < NBTERMS; k++) {
	rbs->cbuttle[k][i] = malloc(sizeof(double)*rbs->nbuttle);
      }
      rbs->w0[i] = malloc(sizeof(double)*rbs->nbk);    
      rbs->w1[i] = malloc(sizeof(double)*rbs->nbk);
      rbs->ek[i] = malloc(sizeof(double)*rbs->nbk);
    }  
    for (i = 0; i < rbs->nkappa; i++) {
      for (n = 0; n < rbs->nbk; n++) {
	fscanf(f, "%d %d %d %lf %lf %lf\n", 
	       &(rbs->basis[i][n]), &kappa,
	       &(rbs->bnode[i][n]), &(rbs->ek[i][n]),
	       &(rbs->w0[i][n]), &(rbs->w1[i][n]));
      }
      rbs->kappa[i] = kappa;
      for (n = 0; n < rbs->nbuttle; n++) {
	fscanf(f, "%d %lf", &kappa, &(rbs->ebuttle[i][n]));
	for (k = 0; k < NBTERMS; k++) {
	  fscanf(f, "%lf", &(rbs->cbuttle[k][i][n]));
	}
      }
      fscanf(f, "\n");
    }
  }
  fclose(f);
}

void WriteRMatrixBasis(char *fn, int fmt) {
  FILE *f;
  int i, n, k, ka, nr;
  ORBITAL *orb;

  f = fopen(fn, "w");
  if (f == NULL) return;

  if (fmt == 0) {
    nr = fwrite(&(rbasis.ib0), sizeof(int), 1, f);
    nr = fwrite(&(rbasis.rb0), sizeof(double), 1, f);
    nr = fwrite(&(rbasis.ib1), sizeof(int), 1, f);
    nr = fwrite(&(rbasis.rb1), sizeof(double), 1, f);
    nr = fwrite(&(rbasis.bqp), sizeof(double), 1, f);
    nr = fwrite(&(rbasis.emin), sizeof(double), 1, f);
    k = rbasis.kmin*1000000 + rbasis.ki*1000 + rbasis.kmax;
    nr = fwrite(&k, sizeof(int), 1, f);
    k = rbasis.nbi*1000 + rbasis.nbk;
    nr = fwrite(&k, sizeof(int), 1, f);
    nr = fwrite(&(rbasis.nkappa), sizeof(int), 1, f);
    nr = fwrite(&(rbasis.nbuttle), sizeof(int), 1, f);
    nr = fwrite(rbasis.kappa, sizeof(int), rbasis.nkappa, f);
    for (i = 0; i < rbasis.nkappa; i++) {
      nr = fwrite(rbasis.basis[i], sizeof(int), rbasis.nbk, f);
      nr = fwrite(rbasis.bnode[i], sizeof(int), rbasis.nbk, f);
      nr = fwrite(rbasis.ek[i], sizeof(double), rbasis.nbk, f);
      nr = fwrite(rbasis.w0[i], sizeof(double), rbasis.nbk, f);
      nr = fwrite(rbasis.w1[i], sizeof(double), rbasis.nbk, f);
      nr = fwrite(rbasis.ebuttle[i], sizeof(double), rbasis.nbuttle, f);
      for (k = 0; k < NBTERMS; k++) {
	nr = fwrite(rbasis.cbuttle[k][i], sizeof(double), rbasis.nbuttle, f);
      }
    }
  } else {
    fprintf(f, "%4d %15.8E %4d %15.8E %15.8E\n", 
	    rbasis.ib0, rbasis.rb0, rbasis.ib1, rbasis.rb1, rbasis.bqp);
    fprintf(f, "%2d %3d %3d %3d\n",
	    rbasis.kmin*1000000 + rbasis.ki*1000 + rbasis.kmax,
	    rbasis.nbi*1000+rbasis.nbk,
	    rbasis.nkappa, rbasis.nbuttle);
    for (i = 0; i < rbasis.nkappa; i++) {
      ka = rbasis.kappa[i];
      for (n = 0; n < rbasis.nbk; n++) {
	fprintf(f, "%4d %3d %3d %15.8E %15.8E %15.8E\n",
		rbasis.basis[i][n], ka, rbasis.bnode[i][n], rbasis.ek[i][n],
		rbasis.w0[i][n], rbasis.w1[i][n]);
      }
      for (n = 0; n < rbasis.nbuttle; n++) {
	fprintf(f, "%4d %15.8E",
		ka, rbasis.ebuttle[i][n]);
	for (k = 0; k < NBTERMS; k++) {
	  fprintf(f, " %15.8E", rbasis.cbuttle[k][i][n]);
	}
	fprintf(f, "\n");
      }
    }
  }
  fclose(f);
}

void RMatrixBoundary(double r0, double r1, double b) {
  int i, n, nmax;
  ORBITAL *orb;
  POTENTIAL *pot;

  if (r0 == 0) {
    n = GetNumOrbitals();
    nmax = 0;
    for (i = 0; i < n; i++) {
      orb = GetOrbital(i);
      if (nmax < orb->n) nmax = orb->n;
    }
    SetBoundary(nmax, r1, b, 0, 0, 0, NULL, NULL, 0, 0, 0, 0, 0, 0);
    pot = RadialPotential();
    pot->ib1 = pot->ib;
    rbasis.ib0 = 0;    
    rbasis.ib1 = pot->ib;
  } else {
    pot = RadialPotential();
    if (pot->ib1 > pot->ib) {
      ReinitRadial(2);
    }
    n = 0;
    for (i = 0; i < nts; i++) {
      if (ts[i] > n) n = ts[i];
    }
    for (i = 0; i < ncs; i++) {
      if (cs[i] > n) n = cs[i];
    }
    ClearRMatrixLevels(n+1);
    pot->bqp = b;
    if (r1 < 0) {
      if (pot->ib1 > pot->ib) {
	r1 = (1-r1)*pot->rad[pot->ib1] - pot->rad[pot->ib];
      } else {
	r1 = (1-r1)*pot->rad[pot->ib];
      }
    }
    pot->ib = pot->ib1;
    SetBoundary(-100, r1, r0, 0, 0, 0, NULL, NULL, 0, 0, 0, 0, 0, 0);
    rbasis.ib0 = pot->ib;
    rbasis.ib1 = pot->ib1;
  }
  printf("rmatrix boundary: %d %g %d %g\n",
	 rbasis.ib0, pot->rad[rbasis.ib0],
	 rbasis.ib1, pot->rad[rbasis.ib1]);
}

void ExtrapolateButtle(RBASIS *rbs, int t, int m, double *e, 
		       double *r0, double *r1, double *r2) {
  int nb, i, n, k;
  double cp0[2], cp1[2], cp2[3];
  double xb1[NBFIT], xb2[NBFIT], yb0[NBFIT], yb1[NBFIT], eb[NBFIT];
  double a0, a1, b, ep;
  
  nb = rbs->nbk;
  if (nb <= NBFIT) {
    for (k = 0; k < m; k++) {
      r0[k] = 0.0;
      r1[k] = 0.0;
      r2[k] = 0.0;
    }
    return;
  }
  for (i = 0; i < NBFIT; i++) {
    n = nb - NBFIT + i;
    xb1[i] = n;
    xb2[i] = 1.0/xb1[i];
    yb0[i] = fabs(rbs->w0[t][n]);
    yb1[i] = fabs(rbs->w1[t][n]);
    eb[i] = rbs->ek[t][n];
  }
  if (rbs->ib0 > 0) {
    PolyFit(2, cp0, NBFIT, xb2, yb0);
  }
  PolyFit(2, cp1, NBFIT, xb2, yb1);
  PolyFit(3, cp2, NBFIT, xb1, eb);  
  for (k = 0; k < m; k++) {
    r0[k] = 0.0;
    r1[k] = 0.0;
    r2[k] = 0.0;
    for (i = 0; i < 25000; i++) {
      n = i + nb;
      b = 1.0/n;
      a1 = cp1[0] + b*cp1[1];
      ep = cp2[0] + cp2[1]*n + cp2[2]*n*n;
      if (rbs->ib0 > 0) {
	a0 = cp0[0] + b*cp0[1];
	r0[k] += 0.5*a0*a0/(ep - e[k]);
	b = 0.5*a0*a1;
	if (IsOdd(i)) b = -b;
	b /= ep - e[k];
	r2[k] += b;
      }
      r1[k] += 0.5*a1*a1/(ep - e[k]);
    }
    n++;
    if (rbs->ib0 > 0) {
      a0 = 0.5*cp0[0]*cp0[0]/cp2[2];
      r0[k] += a0/n;
      a0 = 0.5*cp0[0]*cp1[0]/cp2[2];
      r2[k] += 0.5*(a0/n - a0/(n+1.0));
      if (rbs->w0[t][nb-1]*rbs->w1[t][nb-1] > 0) r2[k] = -r2[k];
    }
    a0 = 0.5*cp1[0]*cp1[0]/cp2[2];
    r1[k] += a0/n;
  }
}
    
int RMatrixBasis(char *fn, int kmax, int nb) {
  int k, i, j, k2, ib0, ib1, kappa, t, in;
  int nkb0, nkb1, n, n0, nmax, kb, ki, kmin, ni;
  double e0, e1, ep, a0, a1, b, rb0, rb1, bb, c0, c1;
  double r01, r10, r0, r1, r2, p0, p1, q0, q1, x0, x1;
  ORBITAL *orb, orbf;
  POTENTIAL *pot;

  double wt0 = WallTime();
  ClearRMatrixBasis(&rbasis);
  pot = RadialPotential();  
  nmax = pot->nb;
  ib0 = rbasis.ib0;
  ib1 = rbasis.ib1;
  rb0 = pot->rad[ib0];
  rb1 = pot->rad[ib1];
  bb = pot->bqp;

  rbasis.rb0 = rb0;
  rbasis.rb1 = rb1;
  rbasis.bqp = bb;
  rbasis.kmin = kmax/1000000;
  kmax = kmax%1000000;
  rbasis.ki = kmax/1000;
  rbasis.kmax = kmax%1000;
  kmin = rbasis.kmin;
  kmax = rbasis.kmax;
  rbasis.kapmin = 0;
  if (rbasis.kmin > 0) {
    rbasis.kapmin = 2*(rbasis.kmin-1) + (rbasis.kmin>1?2:1);
  }
  ni = nb/1000;
  nb = nb%1000;
  rbasis.nbk = nb;
  rbasis.nbi = ni;
  rbasis.nkappa = 2*(kmax-kmin) + (kmin==0?1:2);
  rbasis.nbuttle = nb;
  rbasis.basis = malloc(sizeof(int *)*rbasis.nkappa);
  rbasis.bnode = malloc(sizeof(int *)*rbasis.nkappa);
  rbasis.kappa = malloc(sizeof(int)*rbasis.nkappa);
  rbasis.ebuttle = malloc(sizeof(double *)*rbasis.nkappa);
  for (kb = 0; kb < NBTERMS; kb++) {
    rbasis.cbuttle[kb] = malloc(sizeof(double *)*rbasis.nkappa);
  }
  rbasis.w0 = malloc(sizeof(double *)*rbasis.nkappa);
  rbasis.w1 = malloc(sizeof(double *)*rbasis.nkappa);
  rbasis.ek = malloc(sizeof(double *)*rbasis.nkappa);
  for (i = 0; i < rbasis.nkappa; i++) {
    rbasis.basis[i] = malloc(sizeof(int)*rbasis.nbk);
    rbasis.bnode[i] = malloc(sizeof(int)*rbasis.nbk);
    rbasis.ebuttle[i] = malloc(sizeof(double)*rbasis.nbuttle);
    for (kb = 0; kb < NBTERMS; kb++) {
      rbasis.cbuttle[kb][i] = malloc(sizeof(double)*rbasis.nbuttle);
      for (j = 0; j < rbasis.nbuttle; j++) {
	rbasis.cbuttle[kb][i][j] = 0.0;
      }
    }
    rbasis.w0[i] = malloc(sizeof(double)*rbasis.nbk);    
    rbasis.w1[i] = malloc(sizeof(double)*rbasis.nbk);
    rbasis.ek[i] = malloc(sizeof(double)*rbasis.nbk);
  }
  double wt1 = WallTime();
  nkb0 = GetNumOrbitals();
  for (k = kmin; k <= kmax; k++) {
    n0 = k+1;
    if (rbasis.ib0 == 0 && n0 <= nmax) n0 = nmax + 1;
    k2 = 2*k;
    t = k2 - 1;
    for (j = k2-1; j <= k2+1; j += 2, t++) {     
      if (j < 0) continue;
      kappa = GetKappaFromJL(j, k2);
      rbasis.kappa[t] = kappa;
      for (in = 0; in < nb; in++) {
	if (rbasis.ib0 == 0) {
	  n = in + n0;
	} else {
	  n = -(in + n0);
	}
	int ix = OrbitalExistsNoLock(n, kappa, 0.0);
	if (ix < 0) {
	  orb = GetNewOrbitalNoLock(n, kappa, 0.0);
	  ix = orb->idx;
	}
	rbasis.basis[t][in] = ix;
      }
    }
  }
  ResetWidMPI();
#pragma omp parallel default(shared) private(k, k2, t, j, in, n, n0, kappa)
  {
  for (k = kmin; k <= kmax; k++) {
    k2 = 2*k;
    t = k2 - 1;
    for (j = k2-1; j <= k2+1; j += 2, t++) {  
      int skip = SkipMPI();
      if (skip) continue;   
      if (j < 0) continue;
      double wtt0 = WallTime();
      double emax = -1e30;
      for (in = 0; in < nb; in++) {
	int ix = rbasis.basis[t][in];
	orb = GetOrbitalSolvedNoLock(ix);
	if (orb->energy > emax) emax = orb->energy;
      }
      double wtt1 = WallTime();
      MPrintf(-1, "Basis: %d %d %d %12.5E %10.3E\n",
	      k2, j, nb, emax*HARTREE_EV, wtt1-wtt0);
      fflush(stdout);
    }
  }
  }
  double wt2 = WallTime();
  e1 = 0.0;
  e0 = 0.0;
  double emin = 1E31;
  for (k = kmin; k <= kmax; k++) {
    k2 = 2*k;
    t = k2 - 1;
    for (j = k2-1; j <= k2+1; j += 2, t++) {     
      if (j < 0) continue;
      for (in = 0; in < nb; in++) {
	orb = GetOrbital(rbasis.basis[t][in]);
	if (in == nb-1) {
	  if (orb->energy < emin) emin = orb->energy;
	}
	rbasis.bnode[t][in] = orb->n;
	rbasis.ek[t][in] = orb->energy;
	rbasis.w1[t][in] = WLarge(orb)[ib1];
	rbasis.w0[t][in] = WLarge(orb)[ib0];
	e0 = e1;
	e1 = orb->energy;
	if (in > 0) {
	  rbasis.ebuttle[t][in] = 0.5*(e1 + e0);
	  if (in == 1) {
	    if (e0 > 0) {
	      rbasis.ebuttle[t][0] = -(e1-e0);
	    } else {
	      rbasis.ebuttle[t][0] = e0 - (e1-e0);
	    }
	  }
	}
      }
    }
  }
  ResetWidMPI();
#pragma omp parallel default(shared) private(k, k2, t, j, kappa, i, orbf, r0, r1, r2, in, b, p1, q1, x1, a1, p0, q0, x0, a0, r01, r10, c0, c1)
  {
  for (k = kmin; k <= kmax; k++) {
    k2 = 2*k;
    t = k2 - 1;
    for (j = k2-1; j <= k2+1; j += 2, t++) {     
      if (j < 0) continue;
      int skip = SkipMPI();
      if (skip) continue;
      kappa = GetKappaFromJL(j, k2);
      ExtrapolateButtle(&rbasis, t, rbasis.nbuttle, rbasis.ebuttle[t], 
			rbasis.cbuttle[3][t], rbasis.cbuttle[2][t],
			rbasis.cbuttle[4][t]);
      for (i = 0; i < rbasis.nbuttle; i++) {
	InitOrbitalData(&orbf, 1);
	orbf.n = 1000000;
	orbf.kappa = kappa;
	orbf.energy = rbasis.ebuttle[t][i];
	RadialSolver(&orbf, pot);
	r0 = 0.0;
	r1 = 0.0;
	r2 = 0.0;
	for (in = 0; in < nb; in++) {
	  b = rbasis.w1[t][in];
	  b = b*b/(rbasis.ek[t][in] - orbf.energy);
	  r1 += b;
	  if (ib0 > 0) {
	    b = rbasis.w0[t][in];
	    b = b*b/(rbasis.ek[t][in] - orbf.energy);
	    r0 += b;
	    b = rbasis.w0[t][in]*rbasis.w1[t][in];
	    r2 += b/(rbasis.ek[t][in] - orbf.energy);
	  }
	}
	p1 = WLarge(&orbf)[ib1];
	q1 = WSmall(&orbf)[ib1];
	x1 = 2.0*rb1*q1/FINE_STRUCTURE_CONST - (bb+kappa)*p1;
	a1 = x1/rb1;
	r1 /= 2.0*rb1;
	if (ib0 > 0) {
	  p0 = WLarge(&orbf)[ib0];
	  q0 = WSmall(&orbf)[ib0];
	  x0 = 2.0*rb0*q0/FINE_STRUCTURE_CONST - (bb+kappa)*p0;
	  a0 = x0/rb0;
	  r0 /= 2.0*rb0;
	  r01 = r2/(2.0*rb1);
	  r10 = r2/(2.0*rb0);
	  c0 = p0 - r01*x1 + r0*x0;
	  c1 = p1 - r1*x1 + r10*x0;
	  b = rbasis.cbuttle[4][t][i];
	  rbasis.cbuttle[0][t][i] = (c1 + a0*b)/a1;
	  rbasis.cbuttle[1][t][i] = -(c0 - a1*b)/a0;
	} else {
	  c1 = p1 - r1*x1;
	  rbasis.cbuttle[0][t][i] = c1/a1;
	}
	free(orbf.wfun);
      }
    }
  }
  }
  rbasis.emin = emin;
  WriteRMatrixBasis(fn, fmode);
  nkb1 = GetNumOrbitals();
  double wt3 = WallTime();
  MPrintf(-1, "rmx bas: %d %d %12.5E %g %g %g\n",
	  nkb0, nkb1, emin*HARTREE_EV, wt1-wt0, wt2-wt1, wt3-wt2);
  if (nts > 0 && rbasis.ib0 == 0) {
    PrepSlater(0, nkb0-1, nkb0, nkb1-1, 0, nkb0-1, nkb0, nkb1-1);
    PrepSlater(0, nkb0-1, 0, nkb0-1, nkb0, nkb1-1, nkb0, nkb1-1);
  }

  return 0;
}

int IndexFromKappa(int ka, int k0) {
  int j, k;

  GetJLFromKappa(ka, &j, &k);
  if (j > k) return k-k0;
  else return (k-1)-k0;
}

int KappaFromIndex(int k, int k0) {
  int j;
  
  k += k0;
  if (IsOdd(k)) {
    k++;
    j = k-1;
  } else {
    j = k+1;
  }
  return GetKappaFromJL(j, k);
}

void RMatrixTargets(int nt, int *kt, int nc, int *kc) {
  int nlevs, i, m;
  SYMMETRY *sym;
  STATE *s;
  LEVEL *lev;
  
  if (nts > 0) free(ts);
  if (ncs > 0) free(cs);
  nlevs = SetFrozenTargets(nt, kt, nc, kc, &nts, &ts, &ncs, &cs);
  ntg = nt;
  tg = malloc(sizeof(int)*ntg);
  memcpy(tg, kt, sizeof(int)*ntg);
  ncg = nc;
  cg = malloc(sizeof(int)*ncg);
  memcpy(cg, kc, sizeof(int)*ncg);
}

void ClearRMatrixSurface(RMATRIX *rmx) {
  int i, j;

  if (rmx->nts > 0) {
    free(rmx->et);
    free(rmx->ts);
    free(rmx->pts);
    free(rmx->jts);
  } 
  if (rmx->ncs > 0) {
    free(rmx->ec);
    free(rmx->cs);
    free(rmx->pcs);
    free(rmx->jcs);
  }
  if (rmx->ndim > 0) free(rmx->ek);
  if (rmx->nchan0 > 0) {
    for (i = 0; i < rmx->nchan0; i++) {
      free(rmx->w0[rmx->chans[i]]);
      free(rmx->w1[rmx->chans[i]]);
    }
    free(rmx->chans);
    for (i = 0; i < rmx->nlam; i++) {
      free(rmx->aij[i]);
    }
    for (i = 0; i < 3; i++) {
      for (j = 0; j < dcfg.nr; j++) {
	free(rmx->rmatrix[i][j]);
      }
      free(rmx->rmatrix[i]);
    }
  }
  if (rmx->nlam > 0) free(rmx->aij);
  rmx->nlam = 0;
  rmx->nchan0 = 0;
  rmx->nchan = 0;
  rmx->ndim = 0;
  rmx->ncs = 0;
  rmx->nts = 0;
}  
  
int ReadRMatrixSurface(FILE *f, RMATRIX *rmx, int m, int fmt) {
  int nkappa, isym, nsym, p, j;
  int n, nchan, mchan, nchan0, ndim, i, k, t, ilam;
  int k1, k2, k3, k4, ierr;
  double a, b, z;

  if (fmt == 0) {
    if (m == 0) {
      ierr = fread(&nsym, sizeof(int), 1, f);
      ierr = fread(&mchan, sizeof(int), 1, f);
      ierr = fread(&nts, sizeof(int), 1, f);
      ierr = fread(&ncs, sizeof(int), 1, f);
      ierr = fread(&nkappa, sizeof(int), 1, f);
      ierr = fread(&z, sizeof(double), 1, f);
      ierr = fread(&(rmx->nlam), sizeof(int), 1, f);
      rmx->nsym = nsym;
      rmx->mchan = mchan;
      rmx->nts = nts;
      rmx->ncs = ncs;
      rmx->nkappa = nkappa;
      rmx->z = z;
      nchan = nts*nkappa;
      rmx->et = malloc(sizeof(double)*nts);
      rmx->ec = malloc(sizeof(double)*ncs);
      rmx->w0 = malloc(sizeof(double *)*nchan);
      rmx->w1 = malloc(sizeof(double *)*nchan);
      if (rmx->nlam > 0) {
	rmx->aij = malloc(sizeof(double *)*rmx->nlam);
      }
      for (i = 0; i < nchan; i++) {
	rmx->w0[i] = NULL;
	rmx->w1[i] = NULL;
      }
      rmx->ts = malloc(sizeof(int)*nts);
      rmx->pts = malloc(sizeof(int)*nts);
      rmx->jts = malloc(sizeof(int)*nts);
      rmx->cs = malloc(sizeof(int)*ncs);
      rmx->pcs = malloc(sizeof(int)*ncs);
      rmx->jcs = malloc(sizeof(int)*ncs);
      rmx->et0 = 1E30;
      for (i = 0; i < nts; i++) {
	ierr = fread(&(rmx->ts[i]), sizeof(int), 1, f);
	ierr = fread(&(rmx->pts[i]), sizeof(int), 1, f);
	ierr = fread(&(rmx->jts[i]), sizeof(int), 1, f);
	ierr = fread(&(rmx->et[i]), sizeof(double), 1, f);
	if (rmx->et[i] < rmx->et0) rmx->et0 = rmx->et[i];
      }
      for (i = 0; i < ncs; i++) {
	ierr = fread(&(rmx->cs[i]), sizeof(int), 1, f);
	ierr = fread(&(rmx->pcs[i]), sizeof(int), 1, f);
	ierr = fread(&(rmx->jcs[i]), sizeof(int), 1, f);
	ierr = fread(&(rmx->ec[i]), sizeof(double), 1, f);
      }
      rmx->ndim = 0;
      rmx->nchan0 = 0;
      rmx->ek = NULL;
      rmx->chans = NULL;
      return 0;
    }
    
    if (rmx->ndim > 0) free(rmx->ek);
    if (rmx->nchan0 > 0) {
      for (i = 0; i < rmx->nchan0; i++) {
	free(rmx->w0[rmx->chans[i]]);
	free(rmx->w1[rmx->chans[i]]);
	rmx->w0[rmx->chans[i]] = NULL;
	rmx->w1[rmx->chans[i]] = NULL;
      }
      free(rmx->chans);
      free(rmx->ilev);
      free(rmx->kappa);
      for (i = 0; i < 3; i++) {	
	for (j = 0; j < dcfg.nr; j++) {
	  free(rmx->rmatrix[i][j]);
	}
	free(rmx->rmatrix[i]);
      }
      for (i = 0; i < rmx->nlam; i++) {
	free(rmx->aij[i]);
      }
    }
    
    ierr = fread(&isym, sizeof(int), 1, f);
    ierr = fread(&p, sizeof(int), 1, f);
    ierr = fread(&j, sizeof(int), 1, f);
    ierr = fread(&ndim, sizeof(int), 1, f);
    ierr = fread(&nchan0, sizeof(int), 1, f);
    rmx->isym = isym;
    rmx->p = p;
    rmx->j = j;
    rmx->ndim = ndim;
    rmx->ek = malloc(sizeof(double)*ndim);
    rmx->nchan0 = nchan0;
    rmx->chans = malloc(sizeof(int)*nchan0);
    rmx->kappa = malloc(sizeof(int)*nchan0);
    rmx->ilev = malloc(sizeof(int)*nchan0);
    for (i = 0; i < 3; i++) {
      rmx->rmatrix[i] = malloc(sizeof(double *)*dcfg.nr);
      for (k = 0; k < dcfg.nr; k++) {
	rmx->rmatrix[i][k] = malloc(sizeof(double)*nchan0*nchan0);
      }
    }
    for (i = 0; i < rmx->nlam; i++) {
      rmx->aij[i] = malloc(sizeof(double)*nchan0*nchan0);
    }
    ierr = fread(rmx->ek, sizeof(double), ndim, f);
    for (i = 0; i < nchan0; i++) {
      ierr = fread(&(rmx->chans[i]), sizeof(int), 1, f);
      ierr = fread(&(rmx->ilev[i]), sizeof(int), 1, f);
      ierr = fread(&(rmx->kappa[i]), sizeof(int), 1, f);
    }
    for (i = 0; i < nchan0; i++) {
      rmx->w0[rmx->chans[i]] = malloc(sizeof(double)*ndim);
      rmx->w1[rmx->chans[i]] = malloc(sizeof(double)*ndim);
    }
    for (ilam = 0; ilam < rmx->nlam; ilam++) {
      for (i = 0; i < nchan0; i++) {
	for (t = 0; t <= i; t++) {
	  ierr = fread(&a, sizeof(double), 1, f);	
	  k = t*nchan0 + i;
	  rmx->aij[ilam][k] = a;
	  k = i*nchan0 + t;
	  rmx->aij[ilam][k] = a;
	}
      }
    }
    for (i = 0; i < nchan0; i++) {
      ierr = fread(rmx->w0[rmx->chans[i]], sizeof(double), ndim, f);
      ierr = fread(rmx->w1[rmx->chans[i]], sizeof(double), ndim, f);
    }
  } else {
    if (m == 0) {
      ierr = fscanf(f, "%d %d %d %d %d %lf %d\n", 
		    &nsym, &mchan, &nts, &ncs, &nkappa, &z, &(rmx->nlam));
      if (ierr == EOF) return -1;
      rmx->nsym = nsym;
      rmx->mchan = mchan;
      rmx->nts = nts;
      rmx->ncs = ncs;
      rmx->nkappa = nkappa;
      rmx->z = z;
      nchan = nts*nkappa;
      rmx->et = malloc(sizeof(double)*nts);
      rmx->ec = malloc(sizeof(double)*ncs);
      rmx->w0 = malloc(sizeof(double *)*nchan);
      rmx->w1 = malloc(sizeof(double *)*nchan);
      if (rmx->nlam > 0) {
	rmx->aij = malloc(sizeof(double *)*rmx->nlam);
      }
      for (i = 0; i < nchan; i++) {
	rmx->w0[i] = NULL;
	rmx->w1[i] = NULL;
      }
      rmx->ts = malloc(sizeof(int)*nts);
      rmx->pts = malloc(sizeof(int)*nts);
      rmx->jts = malloc(sizeof(int)*nts);
      rmx->cs = malloc(sizeof(int)*ncs);
      rmx->pcs = malloc(sizeof(int)*ncs);
      rmx->jcs = malloc(sizeof(int)*ncs);
      rmx->et0 = 1E30;
      for (i = 0; i < nts; i++) {
	ierr = fscanf(f, "T %d %d %d %lf\n", 
		      &(rmx->ts[i]), &(rmx->pts[i]),
		      &(rmx->jts[i]), &(rmx->et[i]));
	if (rmx->et[i] < rmx->et0) rmx->et0 = rmx->et[i];
      }
      for (i = 0; i < ncs; i++) {
	ierr = fscanf(f, "C %d %d %d %lf\n", &k, &k1, &k2, &a);
	rmx->cs[i] = k;
	rmx->pcs[i] = k1;
	rmx->jcs[i] = k2;
	rmx->ec[i] = a;
      }    
      rmx->ndim = 0;
      rmx->nchan0 = 0;
      rmx->ek = NULL;
      rmx->chans = NULL;
      return 0;
    }

    if (rmx->ndim > 0) free(rmx->ek);
    if (rmx->nchan0 > 0) {
      for (i = 0; i < rmx->nchan0; i++) {
	free(rmx->w0[rmx->chans[i]]);
	free(rmx->w1[rmx->chans[i]]);
	rmx->w0[rmx->chans[i]] = NULL;
	rmx->w1[rmx->chans[i]] = NULL;
      }
      free(rmx->chans);
      free(rmx->ilev);
      free(rmx->kappa);
      for (i = 0; i < 3; i++) {
	for (k = 0; k < dcfg.nr; k++) {
	  free(rmx->rmatrix[i][k]);
	}
	free(rmx->rmatrix[i]);
      }
      for (i = 0; i < rmx->nlam; i++) {
	free(rmx->aij[i]);
      }
    }
  
    ierr = fscanf(f, "%d %d %d %d %d\n", &isym, &p, &j, &ndim, &nchan0);
    if (ierr == EOF) {
      return -1;
    }
    rmx->isym = isym;
    rmx->p = p;
    rmx->j = j;
    rmx->ndim = ndim;
    rmx->ek = malloc(sizeof(double)*ndim);
    rmx->nchan0 = nchan0;
    rmx->chans = malloc(sizeof(int)*nchan0);
    rmx->kappa = malloc(sizeof(int)*nchan0);
    rmx->ilev = malloc(sizeof(int)*nchan0);
    for (i = 0; i < 3; i++) {
      rmx->rmatrix[i] = malloc(sizeof(double *)*dcfg.nr);
      for (k = 0; k < dcfg.nr; k++) {
	rmx->rmatrix[i][k] = malloc(sizeof(double)*nchan0*nchan0);
      }
    }
    for (i = 0; i < rmx->nlam; i++) {
      rmx->aij[i] = malloc(sizeof(double)*nchan0*nchan0);
    }
    for (i = 0; i < ndim; i++) {
      fscanf(f, "%d %lf\n", &k, &(rmx->ek[i]));
    }  
    for (i = 0; i < nchan0; i++) {
      fscanf(f, "%d %d %d %d %d\n", 
	     &k, &k1, &k2, &k3, &k4);
      rmx->chans[i] = k;
      rmx->ilev[i] = k1;
      rmx->kappa[i] = k4;
      rmx->w0[rmx->chans[i]] = malloc(sizeof(double)*ndim);
      rmx->w1[rmx->chans[i]] = malloc(sizeof(double)*ndim);
    }
    for (ilam = 0; ilam < rmx->nlam; ilam++) {
      for (i = 0; i < nchan0; i++) {
	for (t = 0; t <= i; t++) {
	  fscanf(f, "%d %d %d %d %d %d %d %lf\n",
		 &k, &k1, &k2, &k3, &k1, &k2, &k3, &a);
	  k = t*nchan0 + i;
	  rmx->aij[ilam][k] = a;
	  k = i*nchan0 + t;
	  rmx->aij[ilam][k] = a;
	}
      }
    }
    for (i = 0; i < nchan0; i++) {
      for (n = 0; n < ndim; n++) {
	fscanf(f, "%d %d %lf %lf\n", &k, &t, &a, &b);
	rmx->w0[k][t] = a;
	rmx->w1[k][t] = b;
      }
    }
  }
  return 0;
}  

int WriteRMatrixSurface(FILE *f, double **wik0, double **wik1, int m,
			int fmt, RMATRIX *rmx, HAMILTON *h) {
  int nchan, t, ic, nchan0, i, k, ilev, ka;
  int p, j, jc, i1, k1, ilev1, ka1, ilam, nr;
  double z, a;
  LEVEL *lev;
  
  if (m == 0) {
    if (rmx == NULL) {
      z = GetAtomicNumber();
      z -= GetNumElectrons(ts[0]);
      if (fmt == 0) {
	nr = fwrite(&m, sizeof(int), 1, f);
	nr = fwrite(&m, sizeof(int), 1, f);
	nr = fwrite(&nts, sizeof(int), 1, f);
	nr = fwrite(&ncs, sizeof(int), 1, f);
	nr = fwrite(&(rbasis.nkappa), sizeof(int), 1, f);
	nr = fwrite(&z, sizeof(double), 1, f);
	nr = fwrite(&(dcfg.nmultipoles), sizeof(int), 1, f);
	for (t = 0; t < nts; t++) {
	  ilev = ts[t];
	  lev = GetLevel(ilev);
	  DecodePJ(lev->pj, &k1, &k);
	  nr = fwrite(&ilev, sizeof(int), 1, f);
	  nr = fwrite(&k1, sizeof(int), 1, f);
	  nr = fwrite(&k, sizeof(int), 1, f);
	  nr = fwrite(&(lev->energy), sizeof(double), 1, f);
	}
	for (t = 0; t < ncs; t++) {
	  ilev = cs[t];
	  lev = GetLevel(ilev);
	  DecodePJ(lev->pj, &k1, &k);
	  nr = fwrite(&ilev, sizeof(int), 1, f);
	  nr = fwrite(&k1, sizeof(int), 1, f);
	  nr = fwrite(&k, sizeof(int), 1, f);
	  nr = fwrite(&(lev->energy), sizeof(double), 1, f);	  
	}
      } else {
	fprintf(f, "%3d %3d %3d %3d %3d %10.3E %2d\n",
		m, m, nts, ncs, rbasis.nkappa, z, dcfg.nmultipoles);
	for (t = 0; t < nts; t++) {
	  ilev = ts[t];
	  lev = GetLevel(ilev);
	  DecodePJ(lev->pj, &k1, &k);
	  fprintf(f, "T %3d %3d %3d %15.8E\n", ilev, k1, k, lev->energy);
	}
	for (t = 0; t < ncs; t++) {
	  ilev = cs[t];
	  lev = GetLevel(ilev);
	  DecodePJ(lev->pj, &k1, &k);
	  fprintf(f, "C %3d %d %3d %15.8E\n", ilev, k1, k, lev->energy);
	}
      }
    } else {
      z = rmx->z;
      nts = rmx->nts;
      ncs = rmx->ncs;
      if (fmt == 0) {
	nr = fwrite(&(rmx->nsym), sizeof(int), 1, f);
	nr = fwrite(&(rmx->mchan), sizeof(int), 1, f);
	nr = fwrite(&nts, sizeof(int), 1, f);
	nr = fwrite(&ncs, sizeof(int), 1, f);	
	nr = fwrite(&(rmx->nkappa), sizeof(int), 1, f);
	nr = fwrite(&z, sizeof(double), 1, f);
	nr = fwrite(&(rmx->nlam), sizeof(int), 1, f);
	for (t = 0; t < nts; t++) {
	  nr = fwrite(&(rmx->ts[t]), sizeof(int), 1, f);
	  nr = fwrite(&(rmx->pts[t]), sizeof(int), 1, f);
	  nr = fwrite(&(rmx->jts[t]), sizeof(int), 1, f);
	  nr = fwrite(&(rmx->et[t]), sizeof(double), 1, f);
	}
	for (t = 0; t < ncs; t++) {
	  nr = fwrite(&(rmx->cs[t]), sizeof(int), 1, f);
	  nr = fwrite(&(rmx->pcs[t]), sizeof(int), 1, f);
	  nr = fwrite(&(rmx->jcs[t]), sizeof(int), 1, f);
	  nr = fwrite(&(rmx->ec[t]), sizeof(double), 1, f);
	}
      } else {
	fprintf(f, "%3d %3d %3d %3d %3d %10.3E %2d\n",
		rmx->nsym, rmx->mchan, nts, ncs, rmx->nkappa, z, rmx->nlam);
	for (t = 0; t < nts; t++) {
	  ilev = rmx->ts[t];
	  k = rmx->jts[t];
	  k1 = rmx->pts[t];
	  fprintf(f, "T %3d %3d %3d %15.8E\n", ilev, k1, k, rmx->et[t]);
	}
	for (t = 0; t < ncs; t++) {
	  ilev = rmx->cs[t];
	  k = rmx->jcs[t];
	  k1 = rmx->pcs[t];
	  fprintf(f, "C %3d %3d %3d %15.8E\n", ilev, k1, k, rmx->ec[t]);
	}
      }	
    }
    return 0;
  }

  if (rmx == NULL) {
    DecodePJ(h->pj, &p, &j);
    nchan = rbasis.nkappa * nts;
    nchan0 = 0;
    for (ic = 0; ic < nchan; ic++) {
      if (wik1[ic]) nchan0++;
    }
    if (nchan0 == 0) return 0;
    if (fmt == 0) {
      nr = fwrite(&(h->pj), sizeof(int), 1, f);
      nr = fwrite(&p, sizeof(int), 1, f);
      nr = fwrite(&j, sizeof(int), 1, f);
      nr = fwrite(&(h->dim), sizeof(int), 1, f);
      nr = fwrite(&(nchan0), sizeof(int), 1, f);
      nr = fwrite(h->mixing, sizeof(double), h->dim, f);
      for (ic = 0; ic < nchan; ic++) {
	if (wik1[ic]) {
	  i = ic/rbasis.nkappa;
	  k = ic%rbasis.nkappa;
	  ilev = ts[i];
	  lev = GetLevel(ilev);
	  ka = rbasis.kappa[k];
	  nr = fwrite(&ic, sizeof(int), 1, f);
	  nr = fwrite(&i, sizeof(int), 1, f);
	  nr = fwrite(&ka, sizeof(int), 1, f);
	}
      }
      for (ilam = 0; ilam < dcfg.nmultipoles; ilam++) {
	for (ic = 0; ic < nchan; ic++) {
	  if (wik1[ic] == NULL) continue;
	  i = ic/rbasis.nkappa;
	  ilev = ts[i];
	  k = ic%rbasis.nkappa;
	  ka = rbasis.kappa[k];
	  for (jc = 0; jc <= ic; jc++) {
	    if (wik1[jc] == NULL) continue;
	    i1 = jc/rbasis.nkappa;
	    ilev1 = ts[i1];
	    k1 = jc%rbasis.nkappa;
	    ka1 = rbasis.kappa[k1];
	    a = MultipoleCoeff(h->pj, ilev1, ka1, ilev, ka, ilam+1);
	    nr = fwrite(&a, sizeof(double), 1, f);
	  }
	}
      }
      for (ic = 0; ic < nchan; ic++) {
	if (wik1[ic]) {
	  nr = fwrite(wik0[ic], sizeof(double), h->dim, f);
	  nr = fwrite(wik1[ic], sizeof(double), h->dim, f);
	}
	free(wik0[ic]);
	free(wik1[ic]);
	wik0[ic] = NULL;
	wik1[ic] = NULL;
      }
    } else {
      fprintf(f, "%3d %3d %3d %6d %3d\n", h->pj, p, j, h->dim, nchan0);
      int *im = malloc(sizeof(int)*h->dim);
      ArgSort(h->dim, h->mixing, im);
      for (t = 0; t < h->dim; t++) {
	fprintf(f, "%6d %17.10E\n", t, h->mixing[im[t]]);
      }
      for (ic = 0; ic < nchan; ic++) {
	if (wik1[ic]) {
	  i = ic/rbasis.nkappa;
	  k = ic%rbasis.nkappa;
	  ilev = ts[i];
	  lev = GetLevel(ilev);
	  ka = rbasis.kappa[k];
	  fprintf(f, "%6d %3d %3d %3d %3d\n",
		  ic, i, ilev, k, ka);
	}
      }
      for (ilam = 0; ilam < dcfg.nmultipoles; ilam++) {
	for (ic = 0; ic < nchan; ic++) {
	  if (wik1[ic] == NULL) continue;
	  i = ic/rbasis.nkappa;
	  ilev = ts[i];
	  k = ic%rbasis.nkappa;
	  ka = rbasis.kappa[k];
	  for (jc = 0; jc <= ic; jc++) {
	    if (wik1[jc] == NULL) continue;
	    i1 = jc/rbasis.nkappa;
	    ilev1 = ts[i1];
	    k1 = jc%rbasis.nkappa;
	    ka1 = rbasis.kappa[k1];
	    a = MultipoleCoeff(h->pj, ilev1, ka1, ilev, ka, ilam+1);
	    fprintf(f, "%2d %6d %3d %3d %6d %3d %3d %15.8E\n",
		    ilam+1, ic, ilev, ka, jc, ilev1, ka1, a);
	  }
	}
      }
      for (ic = 0; ic < nchan; ic++) {
	if (wik1[ic]) {
	  for (t = 0; t < h->dim; t++) {
	    fprintf(f, "%6d %6d %15.8E %15.8E\n", 
		    ic, t, wik0[ic][im[t]], wik1[ic][im[t]]);
	  }
	  free(wik0[ic]);
	  free(wik1[ic]);
	  wik0[ic] = NULL;
	  wik1[ic] = NULL;
	}
      }
      free(im);
    }
  } else {
    nchan = rmx->nkappa * rmx->nts;
    nchan0 = rmx->nchan0;
    if (fmt == 0) {
      nr = fwrite(&(rmx->isym), sizeof(int), 1, f);
      nr = fwrite(&(rmx->p), sizeof(int), 1, f);
      nr = fwrite(&(rmx->j), sizeof(int), 1, f);
      nr = fwrite(&(rmx->ndim), sizeof(int), 1, f);
      nr = fwrite(&(nchan0), sizeof(int), 1, f);
      nr = fwrite(rmx->ek, sizeof(double), rmx->ndim, f);
      for (ic = 0; ic < nchan0; ic++) {
	nr = fwrite(&(rmx->chans[ic]), sizeof(int), 1, f);
	nr = fwrite(&(rmx->ilev[ic]), sizeof(int), 1, f);
	nr = fwrite(&(rmx->kappa[ic]), sizeof(int), 1, f);
      }
      for (ilam = 0; ilam < rmx->nlam; ilam++) {
	for (i = 0; i < nchan0; i++) {
	  for (t = 0; t <= i; t++) {
	    k = t*nchan0 + i;
	    nr = fwrite(&(rmx->aij[ilam][k]), sizeof(double), 1, f);
	  }
	}
      }
      for (i = 0; i < nchan0; i++) {
	nr = fwrite(rmx->w0[rmx->chans[i]], sizeof(double), rmx->ndim, f);
	nr = fwrite(rmx->w1[rmx->chans[i]], sizeof(double), rmx->ndim, f);
      }
    } else {
      fprintf(f, "%3d %3d %3d %6d %3d\n", 
	      rmx->isym, rmx->p, rmx->j, rmx->ndim, nchan0);
      int *im = malloc(sizeof(int)*rmx->ndim);
      ArgSort(rmx->ndim, rmx->ek, im);
      for (t = 0; t < rmx->ndim; t++) {
	fprintf(f, "%6d %17.10E\n", t, rmx->ek[im[t]]);
      }
      for (ic = 0; ic < nchan0; ic++) {
	i = rmx->ilev[ic];
	ka = rmx->kappa[ic];
	k = IndexFromKappa(ka, rbasis.kapmin);
	fprintf(f, "%6d %3d %3d %3d %3d\n",
		rmx->chans[ic], i, rmx->ts[i], k, ka);
      }
      for (ilam = 0; ilam < rmx->nlam; ilam++) {
	for (i = 0; i < nchan0; i++) {
	  for (t = 0; t <= i; t++) {
	    k = t*nchan0 + i;
	    a = rmx->aij[ilam][k];
	    ic = rmx->chans[i];
	    jc = rmx->chans[t];
	    ilev = rmx->ts[rmx->ilev[i]];
	    ilev1 = rmx->ts[rmx->ilev[t]];
	    ka = rmx->kappa[i];
	    ka1 = rmx->kappa[t];
	    fprintf(f, "%2d %6d %3d %3d %6d %3d %3d %15.8E\n",
		    ilam+1, ic, ilev, ka, jc, ilev1, ka1, a);
	  }
	}
      }
      for (ic = 0; ic < nchan0; ic++) {
	i = rmx->chans[ic];
	for (t = 0; t < rmx->ndim; t++) {
	  fprintf(f,  "%6d %6d %15.8E %15.8E\n", 
		  i, t, rmx->w0[i][im[t]], rmx->w1[i][im[t]]);
	}	
      }
      free(im);
    }
  }

  return nchan0;
}

int RMatrixSurface(char *fn) {
  int i, m, k, t, kb, ilev, ika, nsym, nchm, nchan0, nb, nbi;
  int j0, p0, j1, kk1, j, p, jmin, jmax, nchan, q, ic, ib;
  int nb0, nb1, na, n, ki, k0, k1, nk, ik, ibk, kk, nbk, nbs1;
  double *mix, **wik0, **wik1;
  SYMMETRY *sym;
  STATE *s;
  LEVEL *lev;
  ORBITAL *orb;
  HAMILTON *h, **hs;
  FILE *f;

  f = fopen(fn, "w");
  if (f == NULL) return -1;
  WriteRMatrixSurface(f, NULL, NULL, 0, fmode, NULL, NULL);
  nchan = nts*rbasis.nkappa;
  wik0 = malloc(sizeof(double *)*nchan);
  wik1 = malloc(sizeof(double *)*nchan);
  for (t = 0; t < nchan; t++) {
    wik0[t] = NULL;
    wik1[t] = NULL;
  }
  nsym = 0;
  nchm = 0;

  n = 0;
  for (i = 0; i < nts; i++) {
    if (ts[i] > n) n = ts[i];
  }
  for (i = 0; i < ncs; i++) {
    if (cs[i] > n) n = cs[i];
  }
  n++;
    
  nbi = rbasis.nbi;
  if (nbi > 0) {
    nb = rbasis.nbk+rbasis.kmax - nbi + 2;
  } else {
    nb = 1;
    nbi = rbasis.nbk+rbasis.kmax+1;
    nb0 = 0;
    nb1 = 0;
  }
  ki = rbasis.ki;
  if (ki > rbasis.kmin) {
    nk = rbasis.kmax-ki+2;
  } else {
    nk = 1;
    ki = rbasis.kmax+1;
    k0 = -1;
    k1 = -1;
  }
  nbk = nb*nk;
  hs = malloc(sizeof(HAMILTON *)*nbk);
  for (ib = 0; ib < nb; ib++) {
    if (rbasis.nbi > 0) {
      if (ib == 0) {
	nb0 = 1;
	nb1 = nbi-1;    
      } else {
	nb0 = nbi+ib-1;
	nb1 = nb0;
      }
    }
    for (ik = 0; ik < nk; ik++) {
      ibk = ib*nk +ik;
      hs[ibk] = malloc(sizeof(HAMILTON)*MAX_SYMMETRIES);
      if (rbasis.ki > rbasis.kmin) {
	if (ik == 0) {
	  k0 = 0;
	  k1 = ki-1;
	} else {
	  k0 = ki+ik-1;
	  k1 = k0;
	}
      }
      int nlevs = GetNumLevels();
      /*
      printf("rmatrix surface ibk: %d %d %d %d %d %d %d %d %d\n",
	     ib, nb, ik, nk, ibk, nb0, nb1, k0, k1);
      */
      if (_rmx_acs || (ib == 0 && ik == 0)) {
	for (i = 0; i < ncs; i++) {
	  lev = GetLevel(cs[i]);
	  DecodePJ(lev->pj, &p0, &j0);
	  AddStateToSymmetry(-(cs[i]+1), -1, j0, p0, j0);
	}
      }
      for (i = 0; i < nts; i++) {
	lev = GetLevel(ts[i]);
	m = lev->pb;
	sym = GetSymmetry(lev->pj);
	s = ArrayGet(&(sym->states), m);
	DecodePJ(lev->pj, &p0, &j0);
	for (k = 0; k < rbasis.nkappa; k++) {
	  for (t = 0; t < rbasis.nbk; t++) {
	    kb = rbasis.basis[k][t];
	    orb = GetOrbital(kb);
	    if (rbasis.nbi > 0) {
	      na = abs(orb->n);	  
	      if (na < nb0) continue;
	      if (na > nb1) continue;
	    }
	    if (rbasis.ki > rbasis.kmin) {
	      kk = GetLFromKappa(orb->kappa)/2;
	      if (kk < k0-_rmx_dk) continue;
	      if (kk > k1+_rmx_dk) continue;
	    }
	    GetJLFromKappa(orb->kappa, &j1, &kk1);
	    p = p0 + kk1/2;
	    jmin = abs(j0 - j1);
	    jmax = j0 + j1;
	    for (j = jmin; j <= jmax; j += 2) {
	      AddStateToSymmetry(-(ts[i]+1), kb, j, p, j);
	    }
	  }
	}
      }

      for (i = 0; i < MAX_SYMMETRIES; i++) {
	double wt0 = WallTime();
	h = GetHamilton(i);
	if (SkipSym(i)) {
	  AllocHamMem(h, -1, -1);
	  AllocHamMem(h, 0, 0);
	  hs[ibk][i].dim = 0;
	  continue;
	}	
	if (rbasis.ib0 == 0) {
	  if (ib == 0 && ik == 0) {
	    nbs1 = ncg;
	  } else {
	    if (_rmx_acs) {
	      nbs1 = -ncg;
	    } else {
	      nbs1 = 0;
	    }
	  }
	  k = ConstructHamiltonFrozen(i, ntg, tg, 0, nbs1, cg,
				      nb0, nb1, k0, k1);
	} else {
	  if (nb1 == 0) {
	    nbs1 = -10000000;
	  } else {
	    nbs1 = -nb1;
	  }
	  k = ConstructHamiltonFrozen(i, ntg, tg, -10000000, 0, NULL,
				      -nb0, nbs1, k0, k1);
	}
	if (k < 0) {
	  AllocHamMem(h, -1, -1);
	  AllocHamMem(h, 0, 0);
	  hs[ibk][i].dim = 0;
	  continue;
	}
	double wt1= WallTime();
	MPrintf(-1, "construct hamilton: %d %d %d %g\n",
		i, h->dim, h->n_basis, wt1-wt0);
	fflush(stdout);
      }
      ResetWidMPI();
#pragma omp parallel default(shared) private(i, h, sym, j, k, s, ilev, kb, orb, ika, ic, nbs1)
      {
	for (i = 0; i < MAX_SYMMETRIES; i++) {
	  int skip = SkipMPI();
	  if (skip) continue;
	  double wt0 = WallTime();
	  h = GetHamilton(i);
	  sym = GetSymmetry(i);
	  if (h->dim <= 0) continue;
	  DiagnolizeHamilton(h);
	  double wt1 = WallTime();
	  MPrintf(-1, "diagonalize hamilton: %d %d %d %g\n",
		  i, h->dim, h->n_basis, wt1-wt0);
	  fflush(stdout);
	  memcpy(&(hs[ibk][i]), h, sizeof(HAMILTON));
	  hs[ibk][i].basis = malloc(sizeof(int)*(size_t)(h->n_basis*2));
	  hs[ibk][i].mixing = malloc(sizeof(double)*(size_t)h->msize);
	  memcpy(hs[ibk][i].mixing, h->mixing, sizeof(double)*(size_t)h->msize);
	  nbs1 = h->n_basis;
	  for (j = 0; j < nbs1; j++) {
	    hs[ibk][i].basis[j] = -1;
	    hs[ibk][i].basis[j+nbs1] = -1;
	    k = h->basis[j];
	    s = ArrayGet(&(sym->states), k);
	    k = -(s->kgroup+1);
	    ilev = IBisect(k, nts, ts);
	    if (ilev >= 0) {
	      kb = s->kcfg;
	      orb = GetOrbital(kb);
	      ika = IndexFromKappa(orb->kappa, rbasis.kapmin);
	      ic = ilev*rbasis.nkappa + ika;	    
	      hs[ibk][i].basis[j] = ic;
	      hs[ibk][i].basis[j+nbs1] = kb;
	    }
	  }
	  if (strlen(_rmx_efn) > 0) {
	    AddToLevels(h, 0, NULL);
	  }
	  AllocHamMem(h, -1, -1);
	  AllocHamMem(h, 0, 0);
	}
      }
      if (strlen(_rmx_efn) > 0) {
	SortLevels(nlevs, -1, 0);
	SaveLevels(_rmx_efn, nlevs, -1);
      }
      if (strlen(_rmx_bfn) > 0) {
	char bfn[256];
	sprintf(bfn, "%s.b%dk%d", _rmx_bfn, ib, ik);
	GetBasisTableLR(bfn, 0, 0, nlevs, -1);
      }
      ClearRMatrixLevels(n+1);
    }
  }
  for (i = 0; i < MAX_SYMMETRIES; i++) {
    double wt0 = WallTime();
    int dim = 0;
    for (ib = 0; ib < nb; ib++) {
      for (ik = 0; ik < nk; ik++) {
	ibk = ib*nk + ik;
	dim += hs[ibk][i].dim;
      }
    }
    if (dim == 0) continue;
    HAMILTON *h0 = GetHamilton(i);
    int dim0 = 0;
    h0->pj = i;
    h0->dim = dim;
    h0->mixing = malloc(sizeof(double)*(size_t)dim);
    for (ib = 0; ib < nb; ib++) {
      for (ik = 0; ik < nk; ik++) {
	ibk = ib*nk + ik;
	h = &(hs[ibk][i]);
	if (h->dim <= 0) continue;
	mix = h->mixing+h->dim;
	nbs1 = h->n_basis;
	memcpy(h0->mixing+dim0, h->mixing, sizeof(double)*(size_t)h->dim);
	for (t = 0; t < h->dim; t++) {
	  int tm = t+dim0;
	  for (q = 0; q < nbs1; q++) {
	    ic = h->basis[q];
	    if (ic >= 0) {
	      kb = h->basis[q+nbs1];
	      orb = GetOrbital(kb);
	      if (wik1[ic] == NULL) {
		wik1[ic] = malloc(sizeof(double)*dim);
		for (p = 0; p < dim; p++) {
		  wik1[ic][p] = 0.0;
		}
	      }
	      wik1[ic][tm] += mix[q]*WLarge(orb)[rbasis.ib1];
	      if (wik0[ic] == NULL) {
		wik0[ic] = malloc(sizeof(double)*dim);
		for (p = 0; p < dim; p++) {
		  wik0[ic][p] = 0.0;
		}
	      }
	      wik0[ic][tm] += mix[q]*WLarge(orb)[rbasis.ib0];
	    }	  
	  }
	  mix += h->n_basis;
	}
	dim0 += h->dim;
      }
    }
    nchan0 = WriteRMatrixSurface(f, wik0, wik1, 1, fmode, NULL, h0);
    if (nchan0 > nchm) nchm = nchan0;
    if (nchan0 > 0) nsym++;
    double wt1 = WallTime();
    printf("rmx sym: %d %d %d %d %g\n", i, h0->dim, nchan0, nchan, wt1-wt0);
    h0->dim = 0;
    free(h0->mixing);
    h0->mixing = NULL;
    fflush(stdout);
  }
  
  fseek(f, 0, SEEK_SET);
  if (fmode == 0) {
    fwrite(&nsym, sizeof(int), 1, f);
    fwrite(&nchm, sizeof(int), 1, f);
  } else {
    fprintf(f, "%3d %3d", nsym, nchm);
  }
  for (ib = 0; ib < nb; ib++) {
    for (ik = 0; ik < nk; ik++) {
      ibk = ib*nk + ik;
      for (i = 0; i < MAX_SYMMETRIES; i++) {
	if (hs[ibk][i].dim > 0) {
	  free(hs[ibk][i].basis);
	  free(hs[ibk][i].mixing);
	}
      }
      free(hs[ibk]);
    }
  }
  free(hs);
  free(wik0);
  free(wik1);
  fclose(f);
  return 0;
}

int RMatrix(double e, RMATRIX *rmx, RBASIS *rbs, int m) {
  int i, j, p, q, nb, k, s;
  double *si0, *si1, *sj0, *sj1;
  double de, a00, a11, a01, a10;
  double *x, *y0, *y1, *y2, b;

  for (i = 0; i < rmx->nchan0; i++) {
    si0 = rmx->w0[rmx->chans[i]];
    si1 = rmx->w1[rmx->chans[i]];
    for (j = 0; j <= i; j++) {
      sj0 = rmx->w0[rmx->chans[j]];
      sj1 = rmx->w1[rmx->chans[j]];
      p = i*rmx->nchan0 + j;
      a00 = 0.0;
      a11 = 0.0;
      a01 = 0.0;
      a10 = 0.0;
      for (k = 0; k < rmx->ndim; k++) {
	de = rmx->ek[k] - rmx->et0;
	de -= e;
	a11 += si1[k]*sj1[k]/de;
	if (rbs->ib0 > 0) {
	  a00 += si0[k]*sj0[k]/de;
	  a01 += si0[k]*sj1[k]/de;
	  a10 += si1[k]*sj0[k]/de;
	}
      }
      a11 *= 0.5;
      if (rbs->ib0 > 0) {
	a00 *= 0.5;
	a01 *= 0.5;
	a10 *= 0.5;
      }
      rmx->rmatrix[0][dcfg.mr][p] = a11;
      rmx->rmatrix[1][dcfg.mr][p] = a00;
      rmx->rmatrix[2][dcfg.mr][p] = a01;
      if (i != j) {	
	p = j*rmx->nchan0 + i;
	rmx->rmatrix[0][dcfg.mr][p] = a11;
	rmx->rmatrix[1][dcfg.mr][p] = a00;
	rmx->rmatrix[2][dcfg.mr][p] = a10;
      } else if (m >= 0) {
	q = rmx->chans[i]/rbs->nkappa;
	k = rmx->chans[i]%rbs->nkappa;
	nb = rbs->nbuttle;
	x = rbs->ebuttle[k];
	if (m == 0) {
	  y0 = rbs->cbuttle[0][k];
	  y1 = rbs->cbuttle[1][k];
	} else {
	  y0 = rbs->cbuttle[2][k];
	  y1 = rbs->cbuttle[3][k];
	}
	de = e - (rmx->et[q] - rmx->et0);
	if (de <= x[nb-1]) {
	  y2 = rbs->cbuttle[4][k];
	  UVIP3P(3, nb, x, y0, 1, &de, &b);
	  rmx->rmatrix[0][dcfg.mr][p] += b;	
	  if (rbs->ib0 > 0) {
	    UVIP3P(3, nb, x, y1, 1, &de, &b);
	    rmx->rmatrix[1][dcfg.mr][p] += b;
	    UVIP3P(3, nb, x, y2, 1, &de, &b);
	    rmx->rmatrix[2][dcfg.mr][p] += b;
	  }
	} else {
	  ExtrapolateButtle(rbs, k, 1, &de, &a00, &a01, &a11);
	  rmx->rmatrix[0][dcfg.mr][p] += a01;
	  if (rbs->ib0 > 0) {
	    rmx->rmatrix[1][dcfg.mr][p] += a00;
	    rmx->rmatrix[2][dcfg.mr][p] += a11;
	  }
	}
      }
    }
  }

  return 0;
}    

int RMatrixPropogate(double *r0, double *r1, RMATRIX *rmx1) {
  int i, j, k, m, p, q, r;
  double a, b;
  int *iwork = dcfg.iwork;

  for (i = 0; i < rmx1->nchan0; i++) {
    for (j = 0; j <= i; j++) {
      p = i*rmx1->nchan0 + j;
      r0[p] += rmx1->rmatrix[1][dcfg.mr][p];
      if (i != j) {
	q = j*rmx1->nchan0 + i;
	r0[q] += rmx1->rmatrix[1][dcfg.mr][q];
	r1[q] = 0.0;
	r1[p] = 0.0;
      } else {
	r1[p] = 1.0;
      }
    }
  }
  m = rmx1->nchan0;
  DGESV(m, m, r0, m, iwork, r1, m, &k);
  memcpy(r0, r1, sizeof(double)*m*m);
  for (i = 0; i < rmx1->nchan0; i++) {
    for (j = 0; j <= i; j++) {
      a = 0.0;
      for (k = 0; k < rmx1->nchan0; k++) {
	p = k*rmx1->nchan0 + j;
	for (m = 0; m < rmx1->nchan0; m++) {
	  q = m*rmx1->nchan0 + i;
	  r = m*rmx1->nchan0 + k;
	  b = r0[r];
	  b *= rmx1->rmatrix[2][dcfg.mr][q]*rmx1->rmatrix[2][dcfg.mr][p];
	  a += b;
	}
      }
      p = i*rmx1->nchan0 + j;
      r1[p] = rmx1->rmatrix[0][dcfg.mr][p] - a;
      if (i != j) {
	q = j*rmx1->nchan0 + i;
	r1[q] = r1[p];
      }
    }
  }
  
  m = rmx1->nchan0;
  memcpy(r0, r1, sizeof(double)*m*m);  

  return 0;
}

int RMatrixKMatrix(RMATRIX *rmx0, RBASIS *rbs, double *r0) {
  double *fs, *fc, *gs, *gc, *a, *b, x;
  int nch, i, j, k, nop, ierr, m, ij, im, mj;
  int *iwork;
  
  iwork = dcfg.iwork;
  a = rmx0->rmatrix[1][dcfg.mr];
  b = rmx0->rmatrix[2][dcfg.mr];
  nop = dcfg.nop;
  if (dcfg.diag) {
    fs = dcfg.fs0;
    fc = dcfg.fc0;
    gs = dcfg.gs0;
    gc = dcfg.gc0;
  } else {
    fs = dcfg.fs;
    fc = dcfg.fc;
    gs = dcfg.gs;
    gc = dcfg.gc;
  }
  nch = rmx0->nchan0;
  if (dcfg.diag) {
    for (i = 0; i < nch; i++) {
      for (j = 0; j < nch; j++) {
	k = rmx0->kappa[j];
	x = (rbs->bqp + k)/rbs->rb1;
	k = j*nch + i;
	a[k] = (r0[k]*x)*fc[j] - 2.0*r0[k]*gc[j]/FINE_STRUCTURE_CONST;
	if (i == j) a[k] += fc[j];
	if (j < nop) {
	  b[k] = (r0[k]*x)*fs[j] - 2.0*r0[k]*gs[j]/FINE_STRUCTURE_CONST;
	  if (i == j) b[k] += fs[j];
	  b[k] = -b[k];
	}
      }
    }
  } else {
    for (i = 0; i < nch; i++) {
      for (j = 0; j < nch; j++) {
	ij = j*nch + i;
	a[ij] = fc[ij];
	if (j < nop) {
	  b[ij] = -fs[ij];
	} else {
	  b[ij] = 0.0;
	}
	for (m = 0; m < nch; m++) {
	  k = rmx0->kappa[m];
	  x = (rbs->bqp + k)/rbs->rb1;
	  im = m*nch + i;
	  mj = j*nch + m;
	  a[ij] += (r0[im]*x)*fc[mj] - 2.0*r0[im]*gc[mj]/FINE_STRUCTURE_CONST;
	  if (j < nop) {
	    b[ij] -= (r0[im]*x)*fs[mj] - 2.0*r0[im]*gs[mj]/FINE_STRUCTURE_CONST;
	  }
	}
      }
    }
  }
  DGESV(nch, nop, a, nch, iwork, b, nch, &ierr);
  return ierr;
}

int SMatrix(RMATRIX *rmx0) {
  double *a, *b, *c, x;
  int i, j, k, p, q;
  int nop, nch, ierr;
  int *iwork = dcfg.iwork;
  
  c = rmx0->rmatrix[2][dcfg.mr];
  a = rmx0->rmatrix[0][dcfg.mr];
  b = rmx0->rmatrix[1][dcfg.mr];
  nch = rmx0->nchan0;
  nop = dcfg.nop;

  for (i = 0; i < nop; i++) {
    for (j = 0; j <= i; j++) {
      x = 0.0;
      for (k = 0; k < nop; k++) {
	p = k*nch + i;
	q = j*nch + k;
	x += c[p]*c[q];
      }
      p = j*nch + i;
      b[p] = -2.0*c[p];
      a[p] = x;
      if (i != j) {
	q = i*nch + j;
	a[q] = x;
	b[q] = -2.0*c[q];
      } else {
	a[p] += 1.0;
      }
    }
  }
  DGESV(nop, nop, a, nch, iwork, b, nch, &ierr);
  for (i = 0; i < nop; i++) {
    for (j = 0; j <= i; j++) {
      x = 0.0;
      for (k = 0; k < nop; k++) {
	p = k*nch + i;
	q = j*nch + k;
	x += b[p]*c[q];
      }
      p = j*nch + i;
      a[p] = x;
      if (j != i) {
	q = i*nch + j;
	a[q] = x;
      }
    }
  }

  return 0;
}

void DCPQ(int *neq, double *t, double *y, double *yd) {
  double r, z, ka, e, a, b;
  
  r = *t;
  z = y[2];
  ka = y[3];
  e = y[4];
  a = ka/r;
  yd[0] = -a*y[0];
  yd[1] = a*y[1];
  a = (e + z/r)*FINE_STRUCTURE_CONST;
  b = 2.0/FINE_STRUCTURE_CONST;
  yd[0] += (b+a)*y[1];
  yd[1] -= a*y[0];
}
FCALLSCSUB4(DCPQ, DCPQ, dcpq, PINT, PDOUBLE, DOUBLEV, DOUBLEV)

void IntegrateDiracCoulomb(double z, int ka, double r, double rt,
			   double e, double *p, double *q) {
  int lrw, liw, itask, iopt, istate, mf, neq, *iwork, itol;
  double y[5], rtol, atol, *rwork, rs;

  y[2] = z;
  y[3] = ka;
  y[4] = e;
  y[0] = *p;
  y[1] = *q;
  neq = 2;
  rwork = dcfg.rwork;
  iwork = dcfg.iwork;
  itol = 1;
  rtol = EPS10;
  atol = EPS10;
  itask = 1;
  istate = 1;
  iopt = 0;
  mf = 10;
  lrw = dcfg.lrw;
  liw = dcfg.liw;
  rs = rt;
  while (rs != r) {
    LSODE(C_FUNCTION(DCPQ, dcpq), neq, y, &rs, r,
	  itol, rtol, &atol, itask, &istate, iopt, rwork,
	  lrw, iwork, liw, NULL, mf);
    if (istate == -1) istate = 2;
    else if (istate < 0) {
      MPrintf(-1, "LSODE Error in IntegrateDiracCoulomb %d\n", istate);
      Abort(1);
    }
  }
  *p = y[0];
  *q = y[1];
}

void ExtDPQ(int *neq, double *t, double *y, double *ydot) {
  double a, b, c, r;
  int i, j, k, nch, ilam;
  double *p, *q, *dp, *dq;
  RMATRIX *rmx;

  rmx = dcfg.rmx;
  r = *t;
  nch = rmx->nchan0;
  p = y;
  q = y + nch;
  dp = ydot;
  dq = ydot + nch;
  for (i = 0; i < nch; i++) {
    a = 0.0;
    b = 0.0;
    for (j = 0; j < nch; j++) {
      k = j*nch + i;
      c = 1.0/(r*r);
      for (ilam = 0; ilam < dcfg.nlam; ilam++) {
	a += rmx->aij[ilam][k]*q[j]*c;
	b += rmx->aij[ilam][k]*p[j]*c;
	c /= r;
      }
    }
    a *= FINE_STRUCTURE_CONST;
    b *= FINE_STRUCTURE_CONST;    
    
    dp[i] = -a;
    dq[i] = b;
    a = rmx->kappa[i]/r;
    dp[i] -= a*p[i];
    dq[i] += a*q[i];
    a = (dcfg.e[i]+rmx->z/r)*FINE_STRUCTURE_CONST;
    b = 2.0/FINE_STRUCTURE_CONST;
    dp[i] += (b + a)*q[i];
    dq[i] -= a*p[i];
  }
}
FCALLSCSUB4(ExtDPQ, EXTDPQ, extdpq, PINT, PDOUBLE, DOUBLEV, DOUBLEV)
    
int IntegrateExternal(RMATRIX *rmx, double r1, double r0) {
  int i, j, k, lrw, liw, itask, iopt, istate, mf, neq, *iwork;
  int itol;
  double *y, rtol, atol, *rwork, rs;

  dcfg.rmx = rmx;
  y = dcfg.p;
  for (i = 0; i < dcfg.nop; i++) {
    for (j = 0; j < rmx->nchan0; j++) {
      k = i*rmx->nchan0 + j;
      y[j] = dcfg.fs[k];
      y[j+rmx->nchan0] = dcfg.gs[k];
    }
    neq = 2*rmx->nchan0;
    rwork = dcfg.rwork;
    iwork = dcfg.iwork;
    itol = 1;
    rtol = EPS10;
    atol = EPS10;
    itask = 1;
    istate = 1;
    iopt = 0;
    mf = 10;
    lrw = dcfg.lrw;
    liw = dcfg.liw;
    rs = r0;
    while (rs != r1) {
      LSODE(C_FUNCTION(EXTDPQ, extdpq), neq, y, &rs, r1,
	    itol, rtol, &atol, itask, &istate, iopt, rwork,
	    lrw, iwork, liw, NULL, mf);
      if (istate == -1) istate = 2;
      else if (istate < 0) {
	printf("LSODE Error in IntegrateExternal %d\n", istate);
	exit(1);
      }
    }
    for (j = 0; j < rmx->nchan0; j++) {
      k = i*rmx->nchan0 + j;
      dcfg.fs[k] = y[j];
      dcfg.gs[k] = y[j+rmx->nchan0];
    }
  }
  for (i = 0; i < rmx->nchan0; i++) {
    for (j = 0; j < rmx->nchan0; j++) {
      k = i*rmx->nchan0 + j;
      y[j] = dcfg.fc[k];
      y[j+rmx->nchan0] = dcfg.gc[k];
    }
    neq = 2*rmx->nchan0;
    rwork = dcfg.rwork;
    iwork = dcfg.iwork;
    itol = 1;
    rtol = EPS10;
    atol = EPS10;
    itask = 1;
    istate = 1;
    iopt = 0;
    mf = 10;
    lrw = dcfg.lrw;
    liw = dcfg.liw;
    rs = r0;
    while (rs != r1) {
      LSODE(C_FUNCTION(EXTDPQ, extdpq), neq, y, &rs, r1,
	    itol, rtol, &atol, itask, &istate, iopt, rwork,
	    lrw, iwork, liw, NULL, mf);
      if (istate == -1) istate = 2;
      else if (istate < 0) {
	printf("LSODE Error in IntegrateExternal %d\n", istate);
	exit(1);
      }
    }
    for (j = 0; j < rmx->nchan0; j++) {
      k = i*rmx->nchan0 + j;
      dcfg.fc[k] = y[j];
      dcfg.gc[k] = y[j+rmx->nchan0];
    }
  }
  dcfg.diag = 0;
  return 0;
}

void TransformQ(RMATRIX *rmx, double b, double r, int m) {
  int i, j, ij;

  if (m == 0) {
    for (i = 0; i < rmx->nchan0; i++) {
      for (j = 0; j < rmx->nchan0; j++) {
	ij = i*rmx->nchan0 + j;
	if (i < dcfg.nop) {
	  dcfg.gs[ij] = 2.0*dcfg.gs[ij]/FINE_STRUCTURE_CONST;
	  dcfg.gs[ij] -= ((b + rmx->kappa[j])/r)*dcfg.fs[ij];
	}
	dcfg.gc[ij] = 2.0*dcfg.gc[ij]/FINE_STRUCTURE_CONST;
	dcfg.gc[ij] -= ((b + rmx->kappa[j])/r)*dcfg.fc[ij];
      }
    }
  } else {
    for (i = 0; i < rmx->nchan0; i++) {
      for (j = 0; j < rmx->nchan0; j++) {
	ij = i*rmx->nchan0 + j;
	if (i < dcfg.nop) {
	  dcfg.gs[ij] += ((b + rmx->kappa[j])/r)*dcfg.fs[ij];
	  dcfg.gs[ij] /= 2.0/FINE_STRUCTURE_CONST;
	}
	dcfg.gc[ij] += ((b + rmx->kappa[j])/r)*dcfg.fc[ij];
	dcfg.gc[ij] /= 2.0/FINE_STRUCTURE_CONST;
      }
    }
  }
}

int PropogateExternal(RMATRIX *rmx, RBASIS *rbs) {
  int i, j, ji, m, jm, mi, nch, nch2, nop, ierr;
  double *x, *a, *y, b;

  nch = rmx->nchan0;
  nch2 = nch*nch;
  nop = dcfg.nop;
  TransformQ(rmx, rbs->bqp, rbs->rb1, 0);
  a = dcfg.rwork;
  x = dcfg.rwork + nch2;
  y = dcfg.p;
  memcpy(a, rmx->rmatrix[2][dcfg.mr], sizeof(double)*nch2);
  for (i = 0; i < nch; i++) {
    for (j = 0; j < nch; j++) {
      ji = j*nch + i;
      if (i == j) {
	x[ji] = 1.0;
      } else {
	x[ji] = 0.0;
      }
    } 
  }
  DGESV(nch, nch, a, nch, dcfg.iwork, x, nch, &ierr);
  for (i = 0; i < nop; i++) {
    for (j = 0; j < nch; j++) {
      ji = i*nch + j;
      b = 0.0;
      for (m = 0; m < nch; m++) {
	jm = m*nch + j;
	mi = i*nch + m;
	b += rmx->rmatrix[0][dcfg.mr][jm]*dcfg.gs[mi];
      }
      b -= dcfg.fs[ji];
      a[ji] = b;
    }
  }
  for (i = 0; i < nop; i++) {    
    for (j = 0; j < nch; j++) {
      ji = i*nch + j;
      b = 0.0;
      for (m = 0; m < nch; m++) {
	jm = m*nch + j;
	mi = i*nch + m;
	b += x[jm]*a[mi];
      }
      y[j] = b;
    }
    for (j = 0; j < nch; j++) {
      ji = i*nch + j;
      b = 0.0;      
      for (m = 0; m < nch; m++) {
	jm = j*nch + m;
	mi = i*nch + m;
	b += rmx->rmatrix[2][dcfg.mr][jm]*dcfg.gs[mi];
	b -= rmx->rmatrix[1][dcfg.mr][jm]*y[m];
      }
      dcfg.fs[ji] = b;
    }
    for (j = 0; j < nch; j++) {
      dcfg.gs[i*nch+j] = y[j];
    }
  }

  for (i = 0; i < nch; i++) {
    for (j = 0; j < nch; j++) {
      ji = i*nch + j;
      b = 0.0;
      for (m = 0; m < nch; m++) {
	jm = m*nch + j;
	mi = i*nch + m;
	b += rmx->rmatrix[0][dcfg.mr][jm]*dcfg.gc[mi];
      }
      b -= dcfg.fc[ji];
      a[ji] = b;
    }
  }
  for (i = 0; i < nch; i++) {    
    for (j = 0; j < nch; j++) {
      ji = i*nch + j;
      b = 0.0;
      for (m = 0; m < nch; m++) {
	jm = m*nch + j;
	mi = i*nch + m;
	b += x[jm]*a[mi];
      }
      y[j] = b;
    }
    for (j = 0; j < nch; j++) {
      ji = i*nch + j;
      b = 0.0;      
      for (m = 0; m < nch; m++) {
	jm = j*nch + m;
	mi = i*nch + m;
	b += rmx->rmatrix[2][dcfg.mr][jm]*dcfg.gc[mi];
	b -= rmx->rmatrix[1][dcfg.mr][jm]*y[m];
      }
      dcfg.fc[ji] = b;
    }
    for (j = 0; j < nch; j++) {
      dcfg.gc[i*nch+j] = y[j];
    }
  }

  TransformQ(rmx, rbs->bqp, rbs->rb0, -1);
  dcfg.diag = 0;

  return 0;
}

void PrepDiracCoulomb(RMATRIX *rmx, RBASIS *rbs, double r) {
  ResetWidMPI();
  double wt0 = WallTime();  
#pragma omp parallel default(shared)
  {
    double e, a, rt, t1, c1, t2, c2;
    int i, j, k, ij, ka, ierr, iter;
    int w = 0;
    for (j = 0; j < rbs->nkappa; j++) {
      for (i = 0; i < rmx->nts; i++) {
	for (k = 0; k < dcfg0.nke; k++) {
	  if (SkipWMPI(w++)) continue;
	  e = dcfg0.ek[k] - (rmx->et[i]-rmx->et0);
	  ij = (i*rbs->nkappa + j)*dcfg0.nke + k;
	  ka = rbs->kappa[j];
	  if (e > 0) {
	    a = 0.5*ka*(ka+1.0);
	    rt = a/(rmx->z + sqrt(rmx->z + 4*a*e));
	    if (_gailitis_exprt > 0) rt *= _gailitis_exprt;
	    if (rt < r) rt = r;
	  } else {
	    rt = r;
	  }
	  ierr = 0;
	  DCOUL(rmx->z, e, ka, rt, &t1, &c1, &t2, &c2, &ierr);
	  if (e > 0) {
	    iter = 0;
	    while (ierr) {
	      iter++;
	      if (iter > _gailitis_expni) {
		MPrintf(-1, "DCOUL error: %d %d %d %d %g %g %g %g %g\n",
			iter, ka, ierr, k, rmx->z,
			e, dcfg0.ek[k], rmx->et[i]-rmx->et0, r);
		break;
	      }
	      rt *= _gailitis_exprf;
	      ierr = 0;
	      DCOUL(rmx->z, e, ka, rt, &t1, &c1, &t2, &c2, &ierr);
	    }
	    if (iter > _gailitis_expni) {
	      ierr = 0;
	      DCOUL(rmx->z, e-(a/r-rmx->z)/r, ka, r, &t1, &c1, &t2, &c2, &ierr);
	      ierr = -9999;
	    } else if (rt > r) {
	      IntegrateDiracCoulomb(rmx->z, ka, r, rt, e, &t1, &c1);
	      IntegrateDiracCoulomb(rmx->z, ka, r, rt, e, &t2, &c2);
	    }
	  }
	  if (e > 0 && ierr != -9999) {
	    dcfg0.afs[ij] = t1;
	    dcfg0.ags[ij] = c1;
	    dcfg0.afc[ij] = t2;
	    dcfg0.agc[ij] = c2;
	  } else {
	    dcfg0.afs[ij] = 0;
	    dcfg0.ags[ij] = 0;
	    dcfg0.afc[ij] = t1;
	    dcfg0.agc[ij] = c1;
	  }
	}      
      }
      MPrintf(0, "PrepDiracCoulomb: %3d %3d %3d %12.5E %11.4E %11.4E\n",
	      j, rbs->kappa[j], rmx->nts, r, WallTime()-wt0, TotalSize());
    }
  }
}

int GailitisExp(RMATRIX *rmx, RBASIS *rbs, double r) {
  int nlam, i, j, n, m, mi, mj, mp, ij, ierr, ilam, nop;
  double a, d, pa, pb, t1, t2, t3, t4;
  double c1, c2, c3, c4, x1, x2, x3, x4;
  double *e, *p2, *p, deps;
  
  e = dcfg.e;
  p2 = dcfg.p2;
  p = dcfg.p;
  
  if (dcfg.nmultipoles > 0) {
    nlam = dcfg.nmultipoles;
    if (nlam > rmx->nlam) nlam = rmx->nlam;
  } else if (dcfg.nmultipoles == 0) {
    nlam = rmx->nlam;
  } else {
    nlam = dcfg.nmultipoles;
  }
  dcfg.nlam = nlam;

  nop = 0;
  for (i = 0; i < rmx->nchan0; i++) {
    e[i] = dcfg.energy - (rmx->et[rmx->ilev[i]] - rmx->et0);
    if (nlam > 0 && dcfg.ngailitis > 1) { 
      p2[i] = 2.0*e[i]*(1.0+0.5*FINE_STRUCTURE_CONST2*e[i]);
      p[i] = sqrt(fabs(p2[i]));
    }
    j = IndexFromKappa(rmx->kappa[i], rbs->kapmin);
    ij = (rmx->ilev[i]*rmx->nkappa + j)*dcfg0.nke + dcfg.ike;
    dcfg.fs0[i] = dcfg0.afs[ij];
    dcfg.gs0[i] = dcfg0.ags[ij];
    dcfg.fc0[i] = dcfg0.afc[ij];
    dcfg.gc0[i] = dcfg0.agc[ij];
    if (e[i] > 0) nop++;
    for (j = 0; j < rmx->nchan0; j++) {
      ij = j*rmx->nchan0 + i;
      if (i == j) {
	dcfg.fs[ij] = dcfg.fs0[i];
	dcfg.gs[ij] = dcfg.gs0[i];
	dcfg.fc[ij] = dcfg.fc0[i];
	dcfg.gc[ij] = dcfg.gc0[i];
      } else {
	dcfg.fs[ij] = 0;
	dcfg.gs[ij] = 0;
	dcfg.fc[ij] = 0;
	dcfg.gc[ij] = 0;
      }
    }
  }
  
  rmx->nop = nop;
  dcfg.nop = nop;

  if (nlam <= 0 || dcfg.ngailitis <= 1) {
    dcfg.diag = 1;
    return nop;
  }
  dcfg.diag = 0;

  a = 1.0;
  for (m = 0; m < dcfg.ngailitis; m++) {
    dcfg.rm[m] = a;
    a /= r;
  }
  if (rmx->z > 1) {
    deps = dcfg.degenerate*rmx->z*rmx->z;
  } else {
    deps = dcfg.degenerate;
  }

  for (n = 0; n < rmx->nchan0; n++) {
    for (i = 0; i < rmx->nchan0; i++) {
      dcfg.a[i] = 0.0;
      dcfg.b[i] = 0.0;
      dcfg.c[i] = 0.0;
      dcfg.d[i] = 0.0;
    }
    dcfg.a[n] = 1.0;
    d = e[n]*FINE_STRUCTURE_CONST/p[n];
    dcfg.d[n] = d;
    a = rmx->z*FINE_STRUCTURE_CONST;
    pa = a*d;
    pb = rmx->z*p[n]/e[n];
    if (e[n] > 0) {
      for (m = 1; m < dcfg.ngailitis; m++) {
	for (i = 0; i < rmx->nchan0; i++) {
	  mi = (m-1)*rmx->nchan0 + i;
	  if (fabs(e[i] - e[n]) > deps) {
	    t1 = -pa*dcfg.a[mi];
	    t1 += (m-1.0-rmx->kappa[i]-rmx->kappa[n])*dcfg.b[mi];
	    t1 += a*dcfg.d[mi];
	    t2 = -pb*dcfg.b[mi];
	    t2 -= (m-1.0-rmx->kappa[i]+rmx->kappa[n])*dcfg.a[mi];
	    t2 += a*dcfg.c[mi];
	  } else {
	    t1 = 0.0;
	    t2 = 0.0;
	  }
	  t3 = -pb*dcfg.d[mi];
	  t3 += (m-1.0+rmx->kappa[i]+rmx->kappa[n])*dcfg.c[mi];
	  t3 += a*dcfg.a[mi];
	  t4 = -pa*dcfg.c[mi];
	  t4 -= (m-1.0+rmx->kappa[i]-rmx->kappa[n])*dcfg.d[mi];
	  t4 += a*dcfg.b[mi];
	  c1 = 0.0;
	  c2 = 0.0;
	  c3 = 0.0;
	  c4 = 0.0;
	  for (ilam = 0; ilam < nlam; ilam++) {
	    mj = m - ilam - 2;
	    if (mj < 0) break;
	    mj = mj*rmx->nchan0;
	    for (j = 0; j < rmx->nchan0; j++) {
	      mp = mj + j;
	      ij = i*rmx->nchan0+j;
	      if (fabs(e[i] - e[n]) > deps) {
		c1 += rmx->aij[ilam][ij]*dcfg.d[mp];
		c2 += rmx->aij[ilam][ij]*dcfg.c[mp];
	      }
	      c3 += rmx->aij[ilam][ij]*dcfg.a[mp];
	      c4 += rmx->aij[ilam][ij]*dcfg.b[mp];
	    }
	  }
	  if (fabs(e[i] - e[n]) > deps) {
	    c1 *= FINE_STRUCTURE_CONST;
	    c2 *= FINE_STRUCTURE_CONST;
	    t1 -= c1;
	    t2 -= c2;
	  } else {
	    for (ilam = 0; ilam < nlam; ilam++) {
	      mj = m - ilam - 1;
	      if (mj < 0) break;
	      mj = mj*rmx->nchan0;
	      for (j = 0; j < rmx->nchan0; j++) {
		mp = mj + j;
		ij = i*rmx->nchan0 + j;
		t1 += rmx->aij[ilam][ij]*dcfg.a[mp];
		t1 += rmx->aij[ilam][ij]*d*dcfg.d[mp];
		t2 += rmx->aij[ilam][ij]*dcfg.b[mp];
		t2 += rmx->aij[ilam][ij]*d*dcfg.c[mp];
	      }
	    }
	    t1 *= p[n]*FINE_STRUCTURE_CONST;
	    t2 *= p[n]*FINE_STRUCTURE_CONST;
	  }
	  c3 *= FINE_STRUCTURE_CONST;
	  c4 *= FINE_STRUCTURE_CONST;
	  t3 -= c3;
	  t4 -= c4;
	  mi = m*rmx->nchan0 + i;
	  if (fabs(e[i] - e[n]) > deps) {
	    x1 = e[i]*FINE_STRUCTURE_CONST;
	    x2 = 2.0/FINE_STRUCTURE_CONST;
	    x3 = p2[n] - p2[i];
	    dcfg.a[mi] = (p[n]*t1 + (x1+x2)*t3)/x3;
	    dcfg.b[mi] = (p[n]*t2 + (x1+x2)*t4)/x3;
	    dcfg.c[mi] = (p[n]*t4 + x1*t2)/x3;
	    dcfg.d[mi] = (p[n]*t3 + x1*t2)/x3;
	  } else {
	    x1 = 2.0*m*e[n]*FINE_STRUCTURE_CONST;
	    x2 = 2.0*m*p[n];
	    x3 = 2.0*rmx->z/p[n];
	    dcfg.a[mi] = -(t2 + (m+rmx->kappa[i]-rmx->kappa[n])*t3)/x1;
	    dcfg.b[mi] = (t1+x3*t3-(m+rmx->kappa[i]+rmx->kappa[n])*t4)/x1;
	    dcfg.c[mi] = (t1+x3*t3+(m-rmx->kappa[i]-rmx->kappa[n])*t4)/x2;
	    dcfg.d[mi] = (-t2+(m-rmx->kappa[i]+rmx->kappa[n])*t3)/x2;
	  }
	}
      }
    } else {
      for (m = 1; m < dcfg.ngailitis; m++) {
	for (i = 0; i < rmx->nchan0; i++) {
	  mi = (m-1)*rmx->nchan0 + i;
	  if (fabs(e[i] - e[n]) > deps) {
	    t1 = pa*dcfg.a[mi];
	    t1 += (m-1.0-rmx->kappa[i]-rmx->kappa[n])*dcfg.b[mi];
	    t1 += a*dcfg.d[mi];
	    t2 = pb*dcfg.b[mi];
	    t2 -= (m-1.0-rmx->kappa[i]+rmx->kappa[n])*dcfg.a[mi];
	    t2 += a*dcfg.c[mi];
	  } else {
	    t1 = 0.0;
	    t2 = 0.0;
	  }
	  t3 = pb*dcfg.d[mi];
	  t3 += (m-1.0+rmx->kappa[i]+rmx->kappa[n])*dcfg.c[mi];
	  t3 += a*dcfg.a[mi];
	  t4 = pa*dcfg.c[mi];
	  t4 -= (m-1.0+rmx->kappa[i]-rmx->kappa[n])*dcfg.d[mi];
	  t4 += a*dcfg.b[mi];
	  c1 = 0.0;
	  c2 = 0.0;
	  c3 = 0.0;
	  c4 = 0.0;
	  for (ilam = 0; ilam < nlam; ilam++) {
	    mj = m - ilam - 2;
	    if (mj < 0) break;
	    mj = mj*rmx->nchan0;
	    for (j = 0; j < rmx->nchan0; j++) {
	      mp = mj + j;
	      ij = i*rmx->nchan0+j;
	      if (fabs(e[i] - e[n]) > deps) {
		c1 += rmx->aij[ilam][ij]*dcfg.d[mp];
		c2 += rmx->aij[ilam][ij]*dcfg.c[mp];
	      }
	      c3 += rmx->aij[ilam][ij]*dcfg.a[mp];
	      c4 += rmx->aij[ilam][ij]*dcfg.b[mp];
	    }
	  }
	  if (fabs(e[i] - e[n]) > deps) {
	    c1 *= FINE_STRUCTURE_CONST;
	    c2 *= FINE_STRUCTURE_CONST;
	    t1 -= c1;
	    t2 -= c2;
	  } else {
	    for (ilam = 0; ilam < nlam; ilam++) {
	      mj = m - ilam - 1;
	      if (mj < 0) break;
	      mj = mj*rmx->nchan0;
	      for (j = 0; j < rmx->nchan0; j++) {
		mp = mj + j;
		ij = i*rmx->nchan0 + j;
		t1 += rmx->aij[ilam][ij]*dcfg.a[mp];
		t1 += rmx->aij[ilam][ij]*d*dcfg.d[mp];
		t2 -= rmx->aij[ilam][ij]*dcfg.b[mp];
		t2 += rmx->aij[ilam][ij]*d*dcfg.c[mp];
	      }
	    }
	    t1 *= p[n]*FINE_STRUCTURE_CONST;
	    t2 *= p[n]*FINE_STRUCTURE_CONST;
	  }
	  c3 *= FINE_STRUCTURE_CONST;
	  c4 *= FINE_STRUCTURE_CONST;
	  t3 -= c3;
	  t4 -= c4;
	  mi = m*rmx->nchan0 + i;
	  if (fabs(e[i] - e[n]) > deps) {
	    x1 = e[i]*FINE_STRUCTURE_CONST;
	    x2 = 2.0/FINE_STRUCTURE_CONST;
	    x3 = p2[i] - p2[n];
	    dcfg.a[mi] = (p[n]*t1 - (x1+x2)*t3)/x3;
	    dcfg.b[mi] = -(p[n]*t2 + (x1+x2)*t4)/x3;
	    dcfg.c[mi] = (p[n]*t4 - x1*t2)/x3;
	    dcfg.d[mi] = -(p[n]*t3 + x1*t2)/x3;
	  } else {
	    x1 = 2.0*m*e[n]*FINE_STRUCTURE_CONST;
	    x2 = 2.0*m*p[n];
	    x3 = 2.0*rmx->z/p[n];
	    dcfg.a[mi] = -(t2 + (m+rmx->kappa[i]-rmx->kappa[n])*t3)/x1;
	    dcfg.b[mi] = (t1-x3*t3-(m+rmx->kappa[i]+rmx->kappa[n])*t4)/x1;
	    dcfg.c[mi] = (t1-x3*t3+(m-rmx->kappa[i]-rmx->kappa[n])*t4)/x2;
	    dcfg.d[mi] = (t2-(m-rmx->kappa[i]+rmx->kappa[n])*t3)/x2;
	  }
	}
      }
    }
    for (i = 0; i < rmx->nchan0; i++) {
      c1 = 0.0;
      c2 = 0.0;
      c3 = 0.0;
      c4 = 0.0;
      for (m = 0; m < dcfg.ngailitis; m++) {
	mi = m*rmx->nchan0 + i;
	x1 = fabs(t1);
	x2 = fabs(t2);
	x3 = fabs(t3);
	x4 = fabs(t4);
	t1 = dcfg.a[mi]*dcfg.rm[m];
	t2 = dcfg.b[mi]*dcfg.rm[m];
	t3 = dcfg.c[mi]*dcfg.rm[m];
	t4 = dcfg.d[mi]*dcfg.rm[m];
	/*
	if (m > 4) {
	  if ((x1 > 0 && fabs(t1) > x1) ||
	      (x2 > 0 && fabs(t2) > x2) ||
	      (x3 > 0 && fabs(t3) > x3) ||
	      (x4 > 0 && fabs(t4) > x4)) break;
	}
	*/
	c1 += t1;
	c2 += t2;
	c3 += t3;
	c4 += t4;	
	/*	
	printf("%2d %2d %2d %2d %10.3E %10.3E %10.3E %10.3E %10.3E\n",
	       rmx->isym, n, i, m, t1, t2, t3, t4, dcfg.rm[m]);
	*/
      }
      mp = n*rmx->nchan0 + i;
      if (e[n] > 0) {
	dcfg.fs[mp] = dcfg.fs0[n]*c1 + (dcfg.gs0[n]/d)*c2;
	dcfg.gs[mp] = (dcfg.gs0[n]/d)*c4 - dcfg.fs0[n]*c3;
	dcfg.fc[mp] = dcfg.fc0[n]*c1 + (dcfg.gc0[n]/d)*c2;
	dcfg.gc[mp] = (dcfg.gc0[n]/d)*c4 - dcfg.fc0[n]*c3;
      } else {
	dcfg.fs[mp] = 0.0;
	dcfg.gs[mp] = 0.0;
	dcfg.fc[mp] = dcfg.fc0[n]*c1 - (dcfg.gc0[n]/d)*c2;
	dcfg.gc[mp] = -(dcfg.gc0[n]/d)*c4 - dcfg.fc0[n]*c3;
      }
    }
  }

  return nop;
}

static void InitDCFG(int nw, int nts, int nk,
		     int nke, double *ek, int n0, double *eo) {
  if (nw > 0) {
    dcfg.nts = 0;
    dcfg.nka = 0;
    dcfg.nke = 0;
    dcfg.ntk = 0;
    dcfg.ike = 0;
    dcfg.ek = NULL;
    memcpy(&dcfg0, &dcfg, sizeof(dcfg));  
#pragma omp parallel default(shared)
    {
      int nw2, ns, n;
      memcpy(&dcfg, &dcfg0, sizeof(dcfg));
      nw2 = nw*nw;
      ns = nw*dcfg.ngailitis;
      dcfg.lrw = 20 + 16*nw*2;
      if (dcfg.lrw < nw2*2) dcfg.lrw = nw2*2;
      dcfg.liw = 20+nw;
      dcfg.nts = 0;
      dcfg.nka = 0;
      dcfg.nke = 0;
      dcfg.ntk = 0;
      n = nw2*4+nw*7+4*ns+dcfg.ngailitis+dcfg.lrw;
      dcfg.dwork = malloc(sizeof(double)*n);
      dcfg.fc = dcfg.dwork;
      dcfg.gc = dcfg.fc + nw2;
      dcfg.fs = dcfg.gc + nw2;
      dcfg.gs = dcfg.fs + nw2;
      dcfg.fc0 = dcfg.gs + nw2;
      dcfg.gc0 = dcfg.fc0 + nw;
      dcfg.fs0 = dcfg.gc0 + nw;
      dcfg.gs0 = dcfg.fs0 + nw;
      dcfg.e = dcfg.gs0 + nw;
      dcfg.p = dcfg.e + nw;
      dcfg.p2 = dcfg.p + nw;
      dcfg.a = dcfg.p2 + nw;
      dcfg.b = dcfg.a + ns;
      dcfg.c = dcfg.b + ns;
      dcfg.d = dcfg.c + ns;
      dcfg.rm = dcfg.d + ns;
      dcfg.rwork = dcfg.rm + dcfg.ngailitis;
      dcfg.iwork = malloc(sizeof(int)*dcfg.liw);
      dcfg.mr = MPIRank(&dcfg.nr);
    }
  }
  if (dcfg0.ntk > 0) {
    free(dcfg0.afc);
    dcfg0.ntk = 0;
    dcfg0.ek = NULL;
  }
  dcfg0.nts = nts;
  dcfg0.nka = nk;
  dcfg0.nke = nke;
  dcfg0.ntk = nts*nk*nke;
  dcfg0.n0 = n0;
  dcfg0.eo = eo;
  if (eo) {
    dcfg0.eo0 = eo[0];
    dcfg0.eo1 = eo[n0-1];
  }
  if (dcfg0.ntk > 0) {
    dcfg0.afc = malloc(sizeof(double)*dcfg0.ntk*4);
    dcfg0.agc = dcfg0.afc + dcfg0.ntk;
    dcfg0.afs = dcfg0.agc + dcfg0.ntk;
    dcfg0.ags = dcfg0.afs + dcfg0.ntk;
    dcfg0.ek = ek;
  }
}

static void ClearDCFG(void) {
#pragma omp parallel default(shared)
  {
    free(dcfg.dwork);
    free(dcfg.iwork);
  }
  if (dcfg0.ntk > 0) free(dcfg.afc);
}

void SortGroupEnergy(RMXCE *rs, RBASIS *rbs, RMATRIX *rmx) {
  double *y, e;
  int *k;
  int i, t, q, p, j, m, n;

  double wt0 = WallTime();
  y = malloc(sizeof(double)*rs->nke);
  k = malloc(sizeof(int)*rs->nke);
  for (i = 0; i < rs->nke; i++) {
    y[i] = rs->e[i];
  }
  ArgSort(rs->nke, y, k);
  for (i = 0; i < rs->nke; i++) {
    y[i] = rs->e[k[i]];
  }
  for (i = 0; i < rs->nke; i++) {
    rs->e[i] = y[i];
  }
  int ns = rmx[0].nts*(rmx[0].nts+1)/2;
  for (t = 0; t < ns; t++) {
    for (i = 0; i < rs->nke; i++) {
      y[i] = rs->s[t][k[i]];
    }
    for (i = 0; i < rs->nke; i++) {
      rs->s[t][i] = y[i];
    }
  }
  
  if (_stark_amp && _stark_idx.n > 0) {
    int nkw = 2*(rbs->kmax-rbs->kmin+1);
    for (t = 0; t < _stark_idx.n; t++) {
      j = IdxGet(&rs->its, _stark_idx.d[t]);
      if (j < 0) continue;
      n = 2*(rmx->jts[j]+1);
      n *= n;      
      for (q = 0; q < n; q++) {	
	for (p = 0; p < nkw; p++) {
	  for (i = 0; i < rs->nke; i++) {
	    y[i] = rs->ap[t][q][p][k[i]];
	  }
	  for (i = 0; i < rs->nke; i++) {
	    rs->ap[t][q][p][i] = y[i];
	  }
	}
      }
    }
  }

  for (m = 0; m < rmx->nsym; m++) {
    if (SkipSym(rmx->isym)) continue;
    for (t = 0; t < _stark_idx.n; t++) {
      j = IdxGet(&rs->its, _stark_idx.d[t]);
      if (j < 0) continue;
      for (q = 0; q < rs->smx[m].nk[t]; q++) {
	for (i = 0; i < rs->nke; i++) {
	  y[i] = rs->smx[m].rp[t][q][k[i]];
	}
	for (i = 0; i < rs->nke; i++) {
	  rs->smx[m].rp[t][q][i] = y[i];	  
	}
	for (i = 0; i < rs->nke; i++) {
	  y[i] = rs->smx[m].ip[t][q][k[i]];
	}
	for (i = 0; i < rs->nke; i++) {
	  rs->smx[m].ip[t][q][i] = y[i];	  
	}
      }
    }
  }
  double de = (rs->e[rs->nke-1]-rs->e[0])/(1.5*rs->nke);
  if (_minsp > 0) {
    while (1) {
      for (i = 0; i < rs->nke-_minsp; i++) {
	if (rs->e[i+_minsp]-rs->e[i] < de) {
	  de /= 1.5;
	  break;
	}
      }
      if (i >= rs->nke-_minsp) break;
    }
  }
  rs->nsp = (int)(1+(rs->e[rs->nke-1]-rs->e[0])/de);
  if (rs->nsp < rs->nke) rs->nsp = rs->nke;
  rs->de = (rs->e[rs->nke-1]-rs->e[0])/(rs->nsp-1);
  free(k);
  free(y);
  rs->isp = malloc(sizeof(int)*rs->nsp);
  t = 0;
  for (i = 0; i < rs->nsp; i++) {
    e = rs->e[0]+i*rs->de;
    for (q = t; q < rs->nke; q++) {
      if (e < rs->e[q]) {
	rs->isp[i] = q-1;
	break;
      }
    }
    t = q;
  }
  rs->nes = 0;
  rs->es = NULL;
  MPrintf(-1, "SortGroupEnergy: %d %d %d %12.5E %12.5E %12.5E %12.5E %11.4E %11.4E\n",
	  rs->nke, rs->nsp, dcfg0.n0,
	  rs->e[0]*HARTREE_EV, rs->e[rs->nke-1]*HARTREE_EV,
	  rs->de*HARTREE_EV, dcfg0.eo1*HARTREE_EV,
	  WallTime()-wt0, TotalSize());
}

int *IdxEgy(int n, double *e0, int n0, double *eo) {
  int i, j, p;
  double et, em1;

  if (eo == NULL) return NULL;
  int *io = malloc(sizeof(int)*n);
  for (i = 0; i < n; i++) {
    io[i] = 0;
  }
  for (i = 0; i < n0; i++) {
    et = 1e31;
    for (j = 0; j < n; j++) {
      em1 = fabs(e0[j]-eo[i]);
      if (em1 < et) {
	p = j;
	et = em1;
      }
    }
    io[p] = 1;
  }
  return io;
}

void SaveRMatrixCE(RMXCE *rs, RBASIS *rbs, RMATRIX *rmx,
		   int iter, char *ofn, double wt0) {
  int *isp, nsp, i, t, q;
  FILE  *f1[4];
  char fn[1024];

  if (_stark_nts > 0) {
    if (dcfg0.n0 > 0) {
      rs->nes = dcfg0.n0;
      rs->es = malloc(sizeof(double)*rs->nes);
      memcpy(rs->es, dcfg0.eo, sizeof(double)*rs->nes);
    } else {
      rs->nes = 2*rs->nke*_stark_nts;
      rs->es = malloc(sizeof(double)*rs->nes);
      t = 0;
      for (q = 0; q < _stark_nts; q++) {
	int ilo = IdxGet(&rs->its, _stark_lower[q]);
	int iup = IdxGet(&rs->its, _stark_upper[q]);
	for (i = 0; i < rs->nke; i++) {
	  rs->es[t++] = fabs(rs->e[i]-(rmx->et[ilo]-rmx->et0));
	  rs->es[t++] = fabs(rs->e[i]-(rmx->et[iup]-rmx->et0));
	}
      }
      for (i = 0; i < rs->nes; i++) {
	if (rs->es[i] <= 0) {
	  rs->es[i] = rs->de*1e-3;
	}
      }
      rs->nes = SortUniqueDouble(rs->nes, rs->es, rs->de*0.5);
      for (i = rs->nes-1; i >= 0; i--) {
	if (rs->es[i] <= dcfg0.eo1) break;
      }
      rs->nes = i;
    }
  }
  
  if (iter >= 0) {
    sprintf(fn, "%s%02d", ofn, iter);
  } else {
    sprintf(fn, "%s", ofn);
  }
  MPrintf(-1, "SaveRMatrixCE Beg: %s %d %d %d %11.4E %11.4E\n",
	  fn, _stark_nts, _stark_amp, _stark_pw, WallTime()-wt0, TotalSize());
  
  f1[0] = fopen(fn, "w");
  for (i = 1; i < 4; i++) {
    f1[i] = NULL;
  }
  char buf[1024];
  if (_stark_nts) {
    sprintf(buf, "%s.ws", fn);
    f1[1] = fopen(buf, "w");
    if (_stark_amp > 1) {
      sprintf(buf, "%s.amp", fn);
      f1[2] = fopen(buf, "w");
    }
    if (_stark_pw) {
      sprintf(buf, "%s.wsp", fn);
      f1[3] = fopen(buf, "w");
    }
  }
  int npw = rbs->kmax-rbs->kmin+1;
  int nkw = 2*npw;
  int ns = rmx[0].nts*(rmx[0].nts+1)/2;
  int nsw = _stark_nts*3+3;
  int ns0 = 0;
  int ns1 = _stark_nts;
  int ns2 = _stark_nts*2;
  int ns3 = _stark_nts*3;
  int its0, its1, st0, k, ika0, ika1, j, p;
  isp = rs->isp;
  nsp = rs->nsp;
  int nes = rs->nes;
  double **sw = NULL;
  double ***swp = NULL;
  if (_stark_nts > 0) {
    sw = malloc(sizeof(double *)*nsw);
    for (i = 0; i < nsw; i++) {
      sw[i] = malloc(sizeof(double)*nes);
      for (j = 0; j < nes; j++) sw[i][j] = 0;
    }
    if (_stark_pw) {
      swp = malloc(sizeof(double **)*_stark_nts*2);
      for (i = 0; i < _stark_nts*2; i++) {
	swp[i] = malloc(sizeof(double *)*nes);
	for (k = 0; k < nes; k++) {
	  swp[i][k] = malloc(sizeof(double)*npw);
	  for (t = 0; t < npw; t++) {
	    swp[i][k][t] = 0.0;
	  }
	}
      }
    }
  }
  SMATRIX *smx = rs->smx;
  double ****ap = rs->ap;
  int nke = rs->nke;
  double *e = rs->e;  
  double **s = rs->s;
  double et, ek;

  int *ik0 = NULL;
  int *ik1 = NULL;
  if (dcfg0.n0 > 0 && dcfg0.eo) {
    ik0 = IdxEgy(nke, rs->e, dcfg0.n0, dcfg0.eo);
    ik1 = IdxEgy(rs->nes, rs->es, dcfg0.n0, dcfg0.eo);
  }
  MPrintf(-1, "save collision strength: %d %11.4E %11.4E\n",
	  rmx[0].nts, WallTime()-wt0, TotalSize());
  for (its0 = 0; its0 < rmx[0].nts; its0++) {
    for (its1 = 0; its1 < rmx[0].nts; its1++) {
      if (rmx[0].et[its1] < rmx[0].et[its0]) continue;
      st0 = its1*(its1+1)/2;
      if (its0 != its1) {
	for (i = 0; i < _stark_nts; i++) {	
	  if (rmx[0].ts[its0] == _stark_lower[i] ||
	      rmx[0].ts[its0] == _stark_upper[i]) {
	    double w0 = 1.0/(1.0+rmx[0].jts[its0]);
	    ResetWidMPI();
#pragma omp parallel default(shared) private(t, et)
	    {
	      int w = 0;
	      for (t = 0; t < rs->nes; t++) {
		if (SkipWMPI(w++)) continue;
		et = rs->es[t] + rmx[0].et[its0] - rmx[0].et0;
		sw[ns2+i][t] += InterpLinear(rs->de, isp, nke, e,
					     s[st0+its0], et)*w0;
	      }
	    }
	  }
	  if (rmx[0].ts[its1] == _stark_lower[i] ||
	      rmx[0].ts[its1] == _stark_upper[i]) {
	    double w0 = 1.0/(1.0+rmx[0].jts[its1]);
	    ResetWidMPI();
#pragma omp parallel default(shared) private(t, et)
	    {
	      int w = 0;
	      for (t = 0; t < nes; t++) {
		if (SkipWMPI(w++)) continue;
		et = rs->es[t] + rmx[0].et[its1] - rmx[0].et0;
		sw[ns2+i][t] += InterpLinear(rs->de, isp, nke, e,
					     s[st0+its0], et)*w0;
	      }
	    }
	  }
	}
      }
      fprintf(f1[0], "# %3d %3d %3d %3d %12.5E %6d %3d\n", 
	      rmx[0].ts[its0], rmx[0].jts[its0],
	      rmx[0].ts[its1], rmx[0].jts[its1],	      
	      (rmx[0].et[its1]-rmx[0].et[its0])*HARTREE_EV, 
	      nke, npw);
      for (k = 0; k < nke; k++) {
	if (ik0 && !ik0[k]) continue;
	et = (e[k]-rmx[0].et[its0]+rmx[0].et0)*HARTREE_EV;
	if (et < 0.0) continue;
	fprintf(f1[0], "%3d %3d %14.8E %15.8E %15.8E %12.5E\n",
		rmx[0].ts[its0], rmx[0].ts[its1],
		e[k]*HARTREE_EV, et,
		(e[k]-rmx[0].et[its1]+rmx[0].et0)*HARTREE_EV,
		s[st0+its0][k]);
      }
      fprintf(f1[0], "\n");
    }
  }

  if (_stark_idx.n > 0) {
    MPrintf(-1, "compute stark width: %d %d %11.4E %11.4E\n",
	    _stark_idx.n, _stark_nts, WallTime()-wt0, TotalSize());
    ResetWidMPI();
#pragma omp parallel default(shared) private(k, ika0, ika1, et, ek)
    {      
      int w = 0;
      int iup, ilo, j0, j1, kl0, kl1, ikup, iklo, is0, is1;
      int xup, xlo, ix;
      double s0r, s0i, s1r, s1i, ssr, ssi, af;
      for (k = 0; k < nes; k++) {
	if (SkipWMPI(w++)) continue;
	if (ik1 && !ik1[k]) continue;
	ek = rs->es[k];
	for (ix = 0; ix < _stark_nts; ix++) {
	  xup = IdxGet(&_stark_idx, _stark_upper[ix]);
	  iup = IdxGet(&rs->its, _stark_upper[ix]);
	  if (iup < 0) continue;
	  xlo = IdxGet(&_stark_idx, _stark_lower[ix]);
	  ilo = IdxGet(&rs->its, _stark_lower[ix]);
	  if (ilo < 0) continue;
	  for (is0 = 0; is0 < rmx[0].nsym; is0++) {
	    if (smx[is0].isym < 0 || smx[is0].nj[xup] == 0) continue;	    
	    for (is1 = 0; is1 < rmx[0].nsym; is1++) {
	      if (smx[is1].isym < 0 || smx[is1].nj[xlo] == 0) continue;
	      int dj = (smx[is0].jmin[xup] - smx[is1].jmin[xlo])/2;
	      //if (!Triangle(smx[is0].jj, smx[is1].jj, 2)) continue;		
	      for (ika0 = 0; ika0 < smx[is0].nj[xup]; ika0++) {
		j0 = smx[is0].jmin[xup] + 2*ika0;
		if (j0 < smx[is1].jmin[xlo] ||
		    j0 > smx[is1].jmax[xlo]) {
		  continue;
		}
		//if (!Triangle(smx[is0].jj, rmx[0].jts[iup], j0)) continue;
		//if (!Triangle(smx[is1].jj, rmx[0].jts[ilo], j0)) continue;
		kl0 = GetChanL(j0, smx[is0].jmin[xup], smx[is0].sp[xup]);
		if (kl0 < rbs[0].kmin || kl0 > rbs[0].kmax) continue;
		int kl0p = GetChanL(j0, smx[is1].jmin[xlo],
				    smx[is1].sp[xlo]);
		if (kl0 != kl0p) continue;
		kl0 -= rbs[0].kmin;
		for (ika1 = 0; ika1 < smx[is0].nj[xup]; ika1++) {
		  j1 = smx[is0].jmin[xup] + 2*ika1;
		  if (j1 < smx[is1].jmin[xlo] ||
		      j1 > smx[is1].jmax[xlo]) {
		    continue;
		  }
		  //if (!Triangle(smx[is0].jj, rmx[0].jts[iup], j1)) continue;
		  //if (!Triangle(smx[is1].jj, rmx[0].jts[ilo], j1)) continue;
		  kl1 = GetChanL(j1, smx[is0].jmin[xup], smx[is0].sp[xup]);
		  if (kl1 < rbs[0].kmin || kl1 > rbs[0].kmax) continue;
		  int kl1p = GetChanL(j1, smx[is1].jmin[xlo],
				      smx[is1].sp[xlo]);
		  if (kl1 != kl1p) continue;
		  kl1 -= rbs[0].kmin;
		  ikup = ika1*smx[is0].nj[xup] + ika0;
		  iklo = (ika1+dj)*smx[is1].nj[xlo] + ika0+dj;
		  double eup, elo;
		  eup = ek + rmx[0].et[iup]-rmx[0].et0;
		  elo = ek + rmx[0].et[ilo]-rmx[0].et0;
		  s0r = InterpLinear(rs->de, isp, nke, e,
				     smx[is0].rp[xup][ikup], eup);
		  s0i = InterpLinear(rs->de, isp, nke, e,
				     smx[is0].ip[xup][ikup], eup);
		  s1r = InterpLinear(rs->de, isp, nke, e,
				     smx[is1].rp[xlo][iklo], elo);
		  s1i = InterpLinear(rs->de, isp, nke, e,
				     smx[is1].ip[xlo][iklo], elo);
		  if (ika1 == ika0) {
		    s0r += 1.0;
		    s1r += 1.0;
		  }
		  ssr = -(s0r*s1r + s0i*s1i);
		  ssi = -(s0i*s1r - s0r*s1i);
		  if (ika1 == ika0) ssr += 1.0;
		  if (!ssr && !ssi) continue;
		  af = (smx[is0].jj+1.0)*(smx[is1].jj+1.0);
		  af *= W6j(smx[is1].jj, smx[is0].jj, 2,
			    rmx[0].jts[iup], rmx[0].jts[ilo], j0);
		  af *= W6j(smx[is1].jj, smx[is0].jj, 2,
			    rmx[0].jts[iup], rmx[0].jts[ilo], j1);
		  int ph = (j0+j1)/2 + rmx[0].jts[iup] + smx[is1].jj;
		  if (IsOdd(ph)) {		    
		    af = -af;
		  }
		  ssr *= af;
		  ssi *= af;
		  sw[ns0+ix][k] += ssr;
		  sw[ns1+ix][k] += ssi;
		  /*
		  if (k == 0) {
		    printf("%2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E\n",
			   is0, is1, smx[is0].isym, smx[is1].isym,
			   iup, ilo, kl0, j0, kl1, j1, ikup, iklo, dj, smx[is0].sp[xup], smx[is1].sp[xlo],
			   s0r, s1r, s0i, s1i, ssr, ssi, af, sw[ns0+ix][k], sw[ns1+ix][k]);
		  }
		  */
		  if (_stark_pw && kl0 < npw) {
		    swp[ns0+ix][k][kl0] += ssr;
		    swp[ns1+ix][k][kl0] += ssi;
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  
  if (_stark_amp > 1) {
    MPrintf(-1, "save stark amplitude: %d %d %11.4E %11.4E\n",
	    _stark_idx.n, _stark_nts, WallTime()-wt0, TotalSize());
    for (i = 0; i < _stark_idx.n; i++) {
      j = IdxGet(&rs->its, _stark_idx.d[i]);
      if (j < 0) continue;
      int mi0, mi1, ms0, ms1, mid;
      mid = 0;
      for (ms0 = -1; ms0 <= 1; ms0 += 2) {
	for (ms1 = -1; ms1 <= 1; ms1 += 2) {
	  for (mi0 = -rmx[0].jts[j]; mi0 <= rmx[0].jts[j]; mi0 += 2) {
	    for (mi1 = -rmx[0].jts[j]; mi1 <= rmx[0].jts[j]; mi1 += 2) {
	      int dm = ((mi0+ms0)-(mi1+ms1))/2;
	      for (k = 0; k < nke; k++) {
		if (ik0 && !ik0[k]) continue;
		et = (e[k] - (rmx[0].et[j]-rmx[0].et0))*HARTREE_EV;
		for (p = 0; p < npw; p++) {
		  q = p + rbs[0].kmin;
		  t = k + nke*p;
		  double rr = ap[i][mid][p][k];
		  double ri = ap[i][mid][p+npw][k];
		  fprintf(f1[2],
			  "%3d %3d %2d %2d %3d %3d %5d %3d %3d %12.5E %12.5E %12.5E %12.5E\n",
			  i, _stark_idx.d[i], ms0, ms1, mi0, mi1, k, q, dm,
			  e[k]*HARTREE_EV, et, rr, ri);
		}
	      }
	      mid++;
	    }
	  }
	}
      }
      fprintf(f1[2], "\n");
    }
  }

  if (_stark_nts > 0) {
    MPrintf(-1, "save stark width: %d %d %11.4E %11.4E\n",
	    _stark_idx.n, _stark_nts, WallTime()-wt0, TotalSize());
    for (i = 0; i < _stark_nts; i++) {
      its0 = IdxGet(&rs->its, _stark_lower[i]);
      if (its0 < 0) continue;
      its1 = IdxGet(&rs->its, _stark_upper[i]);
      if (its1 < 0) continue;
      st0 = its1*(its1+1)/2;
      double *edt = NULL;
      if (_stark_amp) {
	int st0p = its0*(its0+1)/2;
	int js0 = IdxGet(&_stark_idx, _stark_lower[i]);
	int js1 = IdxGet(&_stark_idx, _stark_upper[i]);
	double w0 = rmx[0].jts[its0]+1;
	double w1 = rmx[0].jts[its1]+1;    
	edt = sw[ns3];
	ResetWidMPI();
#pragma omp parallel default(shared) private(k, q, et, ek, p, t)
	{
	  int w = 0;
	  double et0, et1;
	  for (k = 0; k < nes; k++) {
	    if (SkipWMPI(w++)) continue;
	    if (ik1 && !ik1[k]) continue;
	    ek = rs->es[k];
	    et0 = ek + (rmx[0].et[its0] - rmx[0].et0);
	    sw[ns3+1][k] = InterpLinear(rs->de, isp, nke, e,
				       s[st0p+its0], et0)/w0;
	    //printf("edt0: %d %g %g %g\n", k, et, e[k], edt[k]);
	    et1 = ek + (rmx[0].et[its1] - rmx[0].et0);
	    sw[ns3+2][k] = InterpLinear(rs->de, isp, nke, e,
					s[st0+its1], et1)/w1;
	    //printf("edt1: %d %g %g %g\n", k, et, e[k], edt[k]);
	    edt[k] = sw[ns3+1][k] + sw[ns3+2][k];
	    //edt[k] = 0;
	    int nm0 = rmx[0].jts[its0]+1;
	    int nm0s = nm0*nm0;
	    int nm1 = rmx[0].jts[its1]+1;
	    int nm1s = nm1*nm1;
	    int ms0, ms1, mi0, mi1;
	    for (ms0 = -1; ms0 <= 1; ms0 += 2) {
	      int ims0 = ((ms0+1)/2)*(2*nm0s);
	      int ims1 = ((ms0+1)/2)*(2*nm1s);
	      for (ms1 = -1; ms1 <= 1; ms1 += 2) {
		int jms0 = ((ms1+1)/2)*(nm0s) + ims0;
		int jms1 = ((ms1+1)/2)*(nm1s) + ims1;
		for (q = -2; q <= 2; q += 2) {
		  for (mi0 = -rmx[0].jts[its1];
		       mi0 <= rmx[0].jts[its1]; mi0 += 2) {
		    int mi0p = mi0 + q;
		    if (mi0p < -rmx[0].jts[its0] || mi0p > rmx[0].jts[its0]) {
		      continue;
		    }
		    for (mi1 = -rmx[0].jts[its1];
			 mi1 <= rmx[0].jts[its1]; mi1 += 2) {
		      int mi1p = mi1 + q;
		      if (mi1p < -rmx[0].jts[its0] || mi1p > rmx[0].jts[its0]) {
			continue;
		      }
		      /*
		      double f1 = ClebschGordan(rmx[0].jts[its0], mi0p, 2, q,
						rmx[0].jts[its1], mi0);
		      double f2 = ClebschGordan(rmx[0].jts[its0], mi1p, 2, q,
						rmx[0].jts[its1], mi1);
		      double ff = f1*f2/w1;
		      */
		      double f1 = W3j(rmx[0].jts[its0], 2, rmx[0].jts[its1],
				      -mi0p, q, mi0);
		      double f2 = W3j(rmx[0].jts[its0], 2, rmx[0].jts[its1],
				      -mi1p, q, mi1);
		      double ff = f1*f2;
		      if (IsOdd(rmx[0].jts[its1]+(mi0+mi1)/2)) ff = -ff;
		      int ik1 =  jms1 + (mi0+rmx[0].jts[its1])/2*nm1 +
			(mi1+rmx[0].jts[its1])/2;
		      int ik0 = jms0 + (mi0p+rmx[0].jts[its0])/2*nm0 +
			(mi1p+rmx[0].jts[its0])/2; 
		      for (p = 0; p < npw; p++) {
			double *ra1 = ap[js1][ik1][p];
			double *ia1 = ap[js1][ik1][p+npw];
			double *ra0 = ap[js0][ik0][p];
			double *ia0 = ap[js0][ik0][p+npw];
			double rr0 = InterpLinear(rs->de, isp, nke, e,
						  ra0, et0);
			double ri0 = InterpLinear(rs->de, isp, nke, e,
						  ia0, et0);
			double rr1 = InterpLinear(rs->de, isp, nke, e,
						  ra1, et1);
			double ri1 = InterpLinear(rs->de, isp, nke, e,
						  ia1, et1);
			double f = (rr1*rr0+ri1*ri0)*ff;
			edt[k] -= f;
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }    
      fprintf(f1[1], "# %3d %3d %3d %3d %12.5E %6d %6d %3d\n", 
	      rmx[0].ts[its0], rmx[0].jts[its0],
	      rmx[0].ts[its1], rmx[0].jts[its1],	      
	      (rmx[0].et[its1]-rmx[0].et[its0])*HARTREE_EV, 
	      nsp, nes, npw);
      for (k = 0; k < nes; k++) {
	if (ik1 && !ik1[k]) continue;
	ek = rs->es[k];
	fprintf(f1[1],
		"%3d %3d %14.8E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E\n",
		rmx[0].ts[its0], rmx[0].ts[its1],
		ek*HARTREE_EV, sw[ns0+i][k], sw[ns1+i][k],
		sw[ns2+i][k], sw[ns3][k], sw[ns3+1][k], sw[ns3+2][k]);
      }
      fprintf(f1[1], "\n");
      if (_stark_pw) {
	fprintf(f1[3], "# %3d %3d %3d %3d %12.5E %6d %6d %3d\n", 
		rmx[0].ts[its0], rmx[0].jts[its0],
		rmx[0].ts[its1], rmx[0].jts[its1],	      
		(rmx[0].et[its1]-rmx[0].et[its0])*HARTREE_EV, 
		nsp, nes, npw);
	for (k = 0; k < nes; k++) {
	  if (ik1 && !ik1[k]) continue;
	  ek = rs->es[k];
	  for (t = 0; t < npw; t++) {
	    fprintf(f1[3],
		    "%3d %3d %14.8E %3d %12.5E %12.5E\n",
		    rmx[0].ts[its0], rmx[0].ts[its1],
		    ek*HARTREE_EV, t+rbs[0].kmin,
		    swp[ns0+i][k][t], swp[ns1+i][k][t]);
	  }
	}
	fprintf(f1[3], "\n");
      }
    }
  }
  free(isp);
  if (rs->nes > 0) {
    free(rs->es);
  }
  if (_stark_nts > 0) {
    for (i = 0; i < nsw; i++) {
      free(sw[i]);
    }
    free(sw);
    if (_stark_pw) {
      for (i = 0; i < _stark_nts*2; i++) {
	for (k = 0; k < nes; k++) {
	  free(swp[i][k]);
	}
	free(swp[i]);
      }
      free(swp);
    }
  }
  for (i = 0; i < 4; i++) {
    if (f1[i]) fclose(f1[i]);
  }
  if (ik0) free(ik0);
  if (ik1) free(ik1);
  MPrintf(-1, "SaveRMatrixCE End: %s %11.4E %11.4E\n",
	  fn, WallTime()-wt0, TotalSize());
}

int RMatrixCE(char *fn, int np, char *bfn[], char *rfn[], 
	      double emin, double emax,
	      int nst, double *sde, int m, int mb) {
  RBASIS *rbs;
  RMATRIX *rmx;
  FILE **f, *f1[2];
  int ns, i, j, p, q, i0, i1, t, t0, n, nke;
  double **s, *e, *e0, et;
  RMXCE rs;
  char buf[1024];
  
#if USE_MPI == 2
  if (!MPIReady()) InitializeMPI(0, 0);
#endif
  double wt0 = WallTime();
  dcfg.mr = MPIRank(&dcfg.nr);
  rbs = malloc(sizeof(RBASIS)*np);
  f = malloc(sizeof(FILE *)*np);
  rmx = malloc(sizeof(RMATRIX)*np);
  f1[0] = NULL;
  f1[1] = NULL;
  if (m & 1) {
    sprintf(buf, "%s.pw", fn);
    f1[0] = fopen(buf, "w");
  }
  if (m & 2) {
    sprintf(buf, "%s.smx", fn);
    f1[1] = fopen(buf, "w");
  }
  
  double ebmin = 1e31;
  for (i = 0; i < np; i++) {
    ReadRMatrixBasis(bfn[i], &(rbs[i]), fmode);
    if (rbs[i].emin < ebmin) ebmin = rbs[i].emin;
    f[i] = fopen(rfn[i], "r");
    if (f[i] == NULL) return -1;
    ReadRMatrixSurface(f[i], &(rmx[i]), 0, fmode);
  }
  printf("maximum energy vs basis set max: %12.5E %12.5E\n",
	 emax, ebmin*HARTREE_EV);  
  for (i = 0; i < nst; i++) {
    sde[i] /= HARTREE_EV;
  }
  emin /= HARTREE_EV;
  emax /= HARTREE_EV;

  n = 1;
  double em0 = emin;
  double em1;
  int nsp = 1 + (nst-1)/2;
  int *isp = malloc(sizeof(int)*nsp);
  for (i = 0; i < nst; i += 2) {
    if (i < nst-1) {
      em1 = sde[i+1];
    } else {
      em1 = emax;
    }    
    t = (int)(0.1+(em1-em0)/sde[i]);
    if (t <= 0) t = 1;
    sde[i] = (em1-em0)/t;
    isp[i/2] = n-1;
    n += t;
    em0 += t*sde[i];
    if (i < nst-1) {
      sde[i+1] = em0;
    }
  }  
  int na = n;
  if (_stark_eadj && _stark_idx.n > 0) {
    na += n*_stark_idx.n;
  }
  e0 = malloc(sizeof(double)*na);
  e0[0] = emin;
  i = 0;
  for (t = 0; t < nst; t += 2) {
    for (i++; i < n; i++) {    
      e0[i] = e0[i-1] + sde[t];
      if (t < nst-1) {
	if (i == isp[1+t/2]) break;
      }
    }
  }
  int n0 = n;
  double *eo = NULL;
  if (na > n) {
    n0 = n;
    eo = malloc(sizeof(double)*n);
    for (i = 0; i < n; i++) {
      eo[i] = e0[i];
    }
    t = n;
    for (j = 0; j < n; j++) {
      for (i = 0; i < _stark_idx.n; i++) {
	p = _stark_idx.d[i];
	e0[t++] = e0[j] + rmx[0].et[p] - rmx[0].et0;
      }
    }
    qsort(e0, na, sizeof(double), CompareDouble);
    n = na;
    i0 = 0;
    for (i = 1; i < na; i++) {
      if (e0[i]-e0[i0] < 1e-8) {
	e0[i] = 1e31;
	n--;
      } else {
	i0 = i;
      }
    }
    qsort(e0, na, sizeof(double), CompareDouble);
  }
  double de = sde[0];
  for (t = 0; t < rmx[0].nts; t++) {
    et = rmx[0].et[t] - rmx[0].et0;
    i0 = (et - emin)/de;
    i = i0-10;
    if (i < 0) i = 0;
    for (; i <= i0+10 && i < n; i++) {
      if (fabs(et-e0[i]) < 1E-10) {
	if (et < e0[i]) e0[i] += 1E-10;
	else e0[i] -= 1E-10;
      }
    }
  }

  ns = rmx[0].nts*(rmx[0].nts+1)/2;
  e = e0;
  InitDCFG(rmx[0].mchan, 0, 0, 0, NULL, 0, NULL);
  int *ts = malloc(sizeof(int)*rmx[0].nts);
  memcpy(ts, rmx[0].ts, sizeof(int)*rmx[0].nts);
  InitIdxAry(&rs.its, rmx[0].nts, ts);
  rs.nke = 0;
  nke = n;
  if (eo) {
    InitDCFG(0, rmx[0].nts, rbs[0].nkappa, nke, e, n0, eo);
  } else {
    InitDCFG(0, rmx[0].nts, rbs[0].nkappa, nke, e, nke, e);
  }
  RMatrixCEW(np, rbs, rmx, f, f1, fn, &rs, nke, e, m, mb, 0);
  SaveRMatrixCE(&rs, &rbs[0], &rmx[0], -1, fn, wt0);
  SMATRIX *smx = rs.smx;
  double ****ap = rs.ap;
  int npw = rbs[0].kmax-rbs[0].kmin+1;
  int nkw = 2*npw;
  if (_stark_idx.n > 0) {
    for (i = 0; i < rmx[0].nsym; i++) {
      if (smx[i].isym < 0) continue;
      for (j = 0; j < _stark_idx.n; j++) {
	if (smx[i].nj[j] <= 0) continue;
	for (t = 0; t < smx[i].nk[j]; t++) {
	  free(smx[i].rp[j][t]);
	  free(smx[i].ip[j][t]);
	}
	free(smx[i].rp[j]);
	free(smx[i].ip[j]);
      }
      free(smx[i].rp);
      free(smx[i].ip);
      free(smx[i].jmin);
      free(smx[i].jmax);
      free(smx[i].nj);
      free(smx[i].nk);
      free(smx[i].sp);
    }
    free(smx);
    if (_stark_amp) {
      for (i = 0; i < _stark_idx.n; i++) {
	j = IdxGet(&rs.its, _stark_idx.d[i]);
	if (j >= 0) {
	  p = 2*(rmx[0].jts[j]+1);
	  p *= p;
	  for (q = 0; q < p; q++) {
	    for (t = 0; t < nkw; t++) {	      
	      free(ap[i][q][t]);
	    }
	    free(ap[i][q]);
	  }
	  free(ap[i]);
	}
      }
      free(ap);
    }
  }
  for (i = 0; i < np; i++) {
    fclose(f[i]);
    ClearRMatrixBasis(&(rbs[i]));
    ClearRMatrixSurface(&(rmx[i]));
  }
  for (i = 0; i < ns; i++) {
    free(rs.s[i]);
  }
  free(rs.s);
  free(e0);
  if (eo) free(eo);
  free(rs.e);
  free(f);
  free(rbs);
  free(rmx);
  free(isp);
  FreeIdxAry(&rs.its, 0);
  ClearDCFG();
  for (i = 0; i < 2; i++) {
    if (f1[i]) fclose(f1[i]);
  }
  return 0;
}

int RefineRMatrixEGrid(double **er, int idep,
		       RMXCE *rs, RBASIS *rbs, RMATRIX *rmx) {
  double rde;
  int i, j, k, nke, nkr, nde;
  
  SortGroupEnergy(rs, rbs, rmx);
  if (idep > _mrefine) return 0;
  nde = _nrefine;
  if (nde > 0 && idep >= _refiter) {
    nde = 2;
  }
  nke = rs->nke;
  nkr = nke*nde;
  double *erp = malloc(sizeof(double)*nkr);
  double *e = rs->e;
  double **s = rs->s;
  int t = 0;
  for (i = 0; i < nke-1; i++) {
    rde = (e[i+1]-e[i])/nde;
    if (rde < _mineref) continue;
    if (idep > 0) {
      int ir = 0;
      for (k = 0; k < rmx->nts; k++) {
	j = (k*(k+1))/2 + k;
	double ds = fabs(s[j][i+1]-s[j][i]);
	double ms = 0.5*fabs(s[j][i+1]+s[j][i])*_rrefine;
	double a;
	if (ds > ms) {
	  ir = 1;
	  break;
	}
	if (i > 0) {
	  ds = fabs(s[j][i]-s[j][i-1]);	  
	  a = (e[i+1]-e[i])/(e[i]-e[i-1]);
	  ds *= Min(a,10.0);
	  if (ds > ms) {
	    ir = 1;
	    break;
	  }
	}
	if (i < nke-2) {
	  ds = fabs(s[j][i+2]-s[j][i+1]);
	  a = (e[i+1]-e[i])/(e[i+2]-e[i+1]);
	  ds *= Min(a, 10.0);
	  if (ds > ms) {
	    ir = 1;
	    break;
	  }
	}
      }    
      if (!ir) continue;
    }
    for (k = 1; k < nde; k++,t++) {
      erp[t] = e[i]+k*rde;
    }
  }
  erp = realloc(erp, sizeof(double)*t);
  *er = erp;
  MPrintf(0, "RefineRMatrixEGrid: %d %d %d %d\n",
	  idep, nke, nde, t);
  return t;
}

double InterpLinear(double de, int *isp, int nke, double *ea,
		    double *ra, double e) {
  int i, i0;
  double f;
  if (nke == 1) return ra[0];
  if (e <= ea[0]) return ra[0];
  if (e >= ea[nke-1]) return ra[nke-1];
  i = (int)((e-ea[0])/de);
  for (i0 = isp[i]; i0 < nke-1; i0++) {
    if (e <= ea[i0+1]) break;
  }
  f = (e-ea[i0])/(ea[i0+1]-ea[i0]);
  if (f < 0 || f > 1) {
    MPrintf(-1,"invalid interplinear: %d %d %d %d %g %g %g %g %g\n", i, isp[i], i0, nke, de, e, ea[i0], ea[i0+1], f);
  }
  return (1-f)*ra[i0] + f*ra[i0+1];
}

int GetChanL(int j, int j0, int k0) {
  int i = (j - j0)/2;
  int k = k0+i;
  int k2 = 2*k0;
  if (IsOdd(i)) {
    if (j0 > k2) {
      k++;
    } else {
      k--;
    }
  }
  return k;
}

int RMatrixCEW(int np, RBASIS *rbs, RMATRIX *rmx,
	       FILE **f, FILE **f1, char *fn, RMXCE *rs,
	       int nke, double *e, int m, int mb, int idep) {
  int i, j, k, t, p, q, h, i0, ns, *iwork, *ir, nkr, *ipr;
  double et, ***sp, *r0, *r1, x, y, *er;
  int pp, jj, its0, its1, st0, ka0, ka1, mka0, mka1, npw, ika0, ika1;
  SMATRIX *smx;
  
  double wt00 = WallTime();
  for (i = 0; i < np; i++) {
    ClearRMatrixSurface(&(rmx[i]));
    fseek(f[i], 0, SEEK_SET); 
    ReadRMatrixSurface(f[i], &(rmx[i]), 0, fmode);
  }
  ns = rmx[0].nts*(rmx[0].nts+1)/2;
  smx = NULL;
  if (_stark_idx.n > 0) {
    smx = malloc(sizeof(SMATRIX)*rmx[0].nsym);
  }
  
  npw = rbs[0].kmax-rbs[0].kmin+1;
  if (m & 1) {
    sp = malloc(sizeof(double **)*ns);
    for (i = 0; i < ns; i++) {
      sp[i] = malloc(sizeof(double *)*nke);
      for (j = 0; j < nke; j++) {
	sp[i][j] = malloc(sizeof(double)*npw);
	for (k = 0; k < npw; k++) {
	  sp[i][j][k] = 0.0;
	}
      }
    }
  }

  double ****ap = NULL;
  int nkw = 2*npw;
  if (_stark_amp && _stark_idx.n > 0) {
    ap = malloc(sizeof(double **)*_stark_idx.n);
    for (i = 0; i < _stark_idx.n; i++) {
      j = IdxGet(&rs->its, _stark_idx.d[i]);
      if (j < 0) {
	ap[i] = NULL;
      } else {
	p = 2*(rmx[0].jts[j]+1);
	p *= p;
	ap[i] = malloc(sizeof(double *)*p);
	for (q = 0; q < p; q++) {
	  ap[i][q] = malloc(sizeof(double *)*nkw);
	  for (t = 0; t < nkw; t++) {
	    ap[i][q][t] = malloc(sizeof(double)*nke);
	    for (k = 0; k < nke; k++) {
	      ap[i][q][t][k] = 0.0;
	    }
	  }
	}
      }
    }
  }
	 
  if (rs->nke == 0) {
    rs->e = malloc(sizeof(double)*nke);
    for (k = 0; k < nke; k++) {
      rs->e[k] = e[k];
    }
    rs->s = malloc(sizeof(double *)*ns);
    for (i = 0; i < ns; i++) {
      rs->s[i] = malloc(sizeof(double)*nke);
    }
    rs->ap = ap;
    rs->smx = smx;
  } else {
    t = rs->nke+nke;
    rs->e = realloc(rs->e, sizeof(double)*t);
    for (k = rs->nke; k < t; k++) {
      p = k - rs->nke;
      rs->e[k] = e[p];
    }
    for (i = 0; i < ns; i++) {
      rs->s[i] = realloc(rs->s[i], sizeof(double)*t);
    }
  }
  
  double **s = malloc(sizeof(double *)*ns);
  for (i = 0; i < ns; i++) {
    s[i] = rs->s[i]+rs->nke;
    for (k = 0; k < nke; k++) {
      s[i][k] = 0;
    }
  }
  
  double rg = dcfg.rgailitis;
  if (rg < rbs[np-1].rb1) rg = rbs[np-1].rb1;
  PrepDiracCoulomb(&(rmx[0]), &(rbs[0]), rg);
  for (i = 0; i < rmx[0].nsym;  i++) {
    double wt0 = WallTime();
    for (j = 0; j < np; j++) {
      ReadRMatrixSurface(f[j], &(rmx[j]), 1, fmode);
      if (j > 0 && rmx[j].nchan0 != rmx[0].nchan0) {
	printf("inconsistent rmatrix nchan0: %d %d %d %d %d\n",
	       i, rmx[0].isym, j, rmx[j].nchan0, rmx[0].nchan0);
	Abort(1);
      }
    }
    smx[i].isym = -1;
    smx[i].nj = NULL;
    if (SkipSym(rmx[0].isym)) continue;
    DecodePJ(rmx[0].isym, &pp, &jj);
    if (_stark_idx.n > 0) {
      smx[i].pp = pp;
      smx[i].jj = jj;
      smx[i].isym = rmx[0].isym;
      smx[i].jmin = malloc(sizeof(int)*_stark_idx.n);
      smx[i].jmax = malloc(sizeof(int)*_stark_idx.n);
      smx[i].nj = malloc(sizeof(int)*_stark_idx.n);
      smx[i].nk = malloc(sizeof(int)*_stark_idx.n);
      smx[i].sp = malloc(sizeof(int)*_stark_idx.n);
      smx[i].rp = malloc(sizeof(double **)*_stark_idx.n);
      smx[i].ip = malloc(sizeof(double **)*_stark_idx.n);
      for (j = 0; j < _stark_idx.n; j++) {
	its0 = IdxGet(&rs->its, _stark_idx.d[j]);
	if (its0 < 0) {
	  smx[i].sp[j] = 0;	  
	  smx[i].nj[j] = 0;
	  smx[i].nk[j] = 0;
	  continue;
	}
	  
	smx[i].jmin[j] = abs(jj-rmx[0].jts[its0]);
	smx[i].jmax[j] = abs(jj+rmx[0].jts[its0]);
	smx[i].sp[j] = (smx[i].jmin[j]-1)/2;
	if (rmx[0].pts[its0] == smx[i].pp) {
	  if (IsOdd(smx[i].sp[j])) smx[i].sp[j]++;
	} else {
	  if (IsEven(smx[i].sp[j])) smx[i].sp[j]++;
	}
	smx[i].nj[j] = 1+(smx[i].jmax[j]-smx[i].jmin[j])/2;
	smx[i].nk[j] = smx[i].nj[j]*smx[i].nj[j];
	smx[i].rp[j] = malloc(sizeof(double *)*smx[i].nk[j]);
	smx[i].ip[j] = malloc(sizeof(double *)*smx[i].nk[j]);
	for (t = 0; t < smx[i].nk[j]; t++) {
	  smx[i].rp[j][t] = malloc(sizeof(double)*nke);
	  smx[i].ip[j][t] = malloc(sizeof(double)*nke);
	  for (k = 0; k < nke; k++) {
	    smx[i].rp[j][t][k] = 0.0;
	    smx[i].ip[j][t][k] = 0.0;
	  }
	}
      }
    }
    //if (i != 0) continue;
    double wjj = 0.5*(jj+1.0);
    ResetWidMPI();
#pragma omp parallel default(shared) private(k, j, r0, r1, p, its1, st0, ka1, ika1, mka1, q, h, x, its0, ka0, ika0, mka0, et)
    {
      int w = 0;
      for (k = 0; k < nke; k++) {
	if (SkipWMPI(w++)) continue;
	dcfg.energy = e[k];
	for (j = 0; j < np; j++) {
	  RMatrix(e[k], &(rmx[j]), &(rbs[j]), mb);
	}
	r0 = rmx[0].rmatrix[0][dcfg.mr];
	r1 = rmx[0].rmatrix[1][dcfg.mr];
	dcfg.ike = k;
	GailitisExp(&(rmx[0]), &(rbs[0]), rg);
	if (rg > rbs[np-1].rb1) {
	  IntegrateExternal(&(rmx[0]), rbs[np-1].rb1, rg);
	}
	if (dcfg.pdirection >= 0) {
	  for (j = 1; j < np; j++) {
	    RMatrixPropogate(r0, r1, &(rmx[j]));
	  }
	  RMatrixKMatrix(&(rmx[0]), &(rbs[np-1]), r0);
	} else {
	  for (j = np-1; j > 0; j--) {
	    PropogateExternal(&(rmx[j]), &(rbs[j]));
	  }
	  RMatrixKMatrix(&(rmx[0]), &(rbs[0]), r0);
	}
	SMatrix(&(rmx[0]));
	for (p = 0; p < dcfg.nop; p++) {
	  its1 = rmx[0].ilev[p];	      
	  st0 = its1*(its1+1)/2;
	  ka1 = rmx[0].kappa[p];
	  GetJLFromKappa(ka1, &ika1, &mka1);
	  et = e[k] - (rmx[0].et[its1]-rmx[0].et0);
	  if (_stark_idx.n > 0) {
	    for (q = 0; q < dcfg.nop; q++) {
	      its0 = rmx[0].ilev[q];
	      if (its0 != its1) continue;
	      ka0 = rmx[0].kappa[q];
	      GetJLFromKappa(ka0, &ika0, &mka0);
	      h = q*rmx[0].nchan0 + p;
	      int ix = IdxGet(&_stark_idx, rmx[0].ts[its0]);
	      if (ix >= 0) {
		if (ika1 >= smx[i].jmin[ix] &&
		    ika1 <= smx[i].jmax[ix] &&
		    ika0 >= smx[i].jmin[ix] &&
		    ika0 <= smx[i].jmax[ix]) {
		  int jka1 = (ika1 - smx[i].jmin[ix])/2;
		  int jka0 = (ika0 - smx[i].jmin[ix])/2;
		  int ika = jka1*smx[i].nj[ix] + jka0;
		  smx[i].rp[ix][ika][k] = r0[h];
		  smx[i].ip[ix][ika][k] = -r1[h];
		  if (_stark_amp) {
		    int mi0, mi1, ms0, ms1, ii0, ii1, mid;
		    /*
		    double sig0, sig1, sr, si;
		    double sigs = CoulombPhaseShift(rmx[0].z, et, -1);
		    sig0 = CoulombPhaseShift(rmx[0].z, et, ka0) - sigs;
		    sig1 = CoulombPhaseShift(rmx[0].z, et, ka1) - sigs;
		    sig1 += sig0 + (PI/2)*((mka0-mka1)/2);
		    sr = cos(sig1);
		    si = sin(sig1);
		    */
		    mid = 0;
		    for (ms0 = -1; ms0 <= 1; ms0 += 2) {
		      for (ms1 = -1; ms1 <= 1; ms1 += 2) {
			for (mi0 = -rmx[0].jts[its0];
			     mi0 <= rmx[0].jts[its0]; mi0 += 2) {
			  for (mi1 = -rmx[0].jts[its0];
			       mi1 <= rmx[0].jts[its0]; mi1 += 2) {
			    int mt = mi0+ms0;
			    int dm = mt - (mi1+ms1);
			    if (dm > mka1 || dm < -mka1) {
			      mid++;
			      continue;
			    }
			    int mf = dm+ms1;
			    double f1 = ClebschGordan(mka0, 0, 1, ms0,
						      ika0, ms0);
			    double f2 = ClebschGordan(mka1, dm, 1, ms1,
						      ika1, mf);
			    double f3 = ClebschGordan(rmx[0].jts[its0], mi0,
						      ika0, ms0, jj, mt);
			    double f4 = ClebschGordan(rmx[0].jts[its0], mi1,
						      ika1, mf, jj, mt);
			    double ff = f1*f2*f3*f4*sqrt(mka0+1.0);
			    int ikw = (mka1/2-rbs[0].kmin);
			    /*
			    ap[ix][mid][ikw][k] += ff*(r0[h]*si-r1[h]*sr);
			    ap[ix][mid][ikw+npw][k] += ff*(-r0[h]*sr-r1[h]*si);
			    */
			    ap[ix][mid][ikw][k] += -ff*r0[h];
			    ap[ix][mid][ikw+npw][k] += ff*r1[h];
			    mid++;
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	  mka1 /= 2;
	  mka1 -= rbs[0].kmin;
	  for (q = 0; q < dcfg.nop; q++) {
	    its0 = rmx[0].ilev[q];
	    if (its0 > its1) continue;
	    h = q*rmx[0].nchan0 + p;
	    x = r0[h]*r0[h] + r1[h]*r1[h];
	    x *= wjj;
	    ka0 = rmx[0].kappa[q];
	    mka0 = GetLFromKappa(ka0);
	    mka0 /= 2;
	    mka0 -= rbs[0].kmin;
	    if (rbs[0].kmin > 0 && mka1 == rbs[0].kmin) continue;
	    s[st0+its0][k] += x;
	    if (m & 1) {
	      sp[st0+its0][k][mka1] += x;
	    }
	    if (m & 2) {
	      et = (e[k]-rmx[0].et[its0]+rmx[0].et0)*HARTREE_EV;
	      fprintf(f1[1], 
		      "%3d %4d %4d %4d %4d %12.5E %12.5E %12.5E %12.5E %12.5E\n",
		      rmx[0].isym, its0, ka0, its1, ka1, et, 
		      r0[h], r1[h], 
		      x, rmx[0].rmatrix[2][dcfg.mr][h]);
	    }
	  }
	}
      }
    }
    double wt1 = WallTime();
    MPrintf(-1, "sym: %3d %1d %3d %3d %3d %11.4E %11.4E %10.3E %10.3E %10.3E\n",
	    rmx[0].isym, pp, jj, idep, nke,
	    (nke>1?(e[1]-e[0]):0.0)*HARTREE_EV, e[0]*HARTREE_EV,
	    wt1-wt0, wt1-wt00, TotalSize());
    fflush(stdout);
  }
  if (rs->nke > 0) {
    free(s);
    t = rs->nke + nke;
    if (_stark_idx.n > 0) {
      if (_stark_amp) {
	for (i = 0; i < _stark_idx.n; i++) {
	  j = IdxGet(&rs->its, _stark_idx.d[i]);
	  if (j < 0) continue;
	  p = 2*(rmx[0].jts[j]+1);
	  p *= p;
	  for (q = 0; q < p; q++) {
	    for (h = 0; h < nkw; h++) {
	      rs->ap[i][q][h] = realloc(rs->ap[i][q][h], sizeof(double)*t);
	      for (k = rs->nke; k < t; k++) {
		rs->ap[i][q][h][k] = ap[i][q][h][k-rs->nke];
	      }
	    }
	  }
	}
      }
      for (i = 0; i < rmx[0].nsym; i++) {
	for (j = 0; j < _stark_idx.n; j++) {
	  its0 = IdxGet(&rs->its, _stark_idx.d[j]);
	  if (its0 < 0) continue;
	  for (q = 0; q < smx[i].nk[j]; q++) {
	    rs->smx[i].rp[j][q] = realloc(rs->smx[i].rp[j][q],
					  sizeof(double)*t);
	    rs->smx[i].ip[j][q] = realloc(rs->smx[i].ip[j][q],
					  sizeof(double)*t);
	    for (k = rs->nke; k < t; k++) {
	      p = k - rs->nke;
	      rs->smx[i].rp[j][q][k] = smx[i].rp[j][q][p];
	      rs->smx[i].ip[j][q][k] = smx[i].ip[j][q][p];
	    }
	  }
	}
      }
    }
  }

  if (m & 1) {
    for (its0 = 0; its0 < rmx[0].nts; its0++) {
      for (its1 = 0; its1 < rmx[0].nts; its1++) {
	if (rmx[0].et[its1] < rmx[0].et[its0]) continue;
	st0 = its1*(its1+1)/2;
	fprintf(f1[0], "# %3d %3d %3d %3d %12.5E %6d %3d\n", 
		rmx[0].ts[its0], rmx[0].jts[its0],
		rmx[0].ts[its1], rmx[0].jts[its1],	      
		(rmx[0].et[its1]-rmx[0].et[its0])*HARTREE_EV, 
		nke, npw);
	for (k = 0; k < nke; k++) {
	  et = (e[k]-rmx[0].et[its0]+rmx[0].et0)*HARTREE_EV;
	  if (et < 0.0) continue;
	  for (j = 0; j < npw; j++) {
	    fprintf(f1[0], "%3d %3d %14.8E %15.8E %15.8E %3d %11.5E\n",
		    rmx[0].ts[its0], rmx[0].ts[its1],
		    e[k]*HARTREE_EV, et,
		    (e[k]-rmx[0].et[its1]+rmx[0].et0)*HARTREE_EV,
		    j+rbs[0].kmin, sp[st0+its0][k][j]);
	  }
	}
	fprintf(f1[0], "\n");
      }
    }

    for (i = 0; i < ns; i++) {
      for (j = 0; j < nke; j++) {
	free(sp[i][j]);
      }
      free(sp[i]);
    }
    free(sp);
  }
  int onke = rs->nke;
  rs->nke += nke;
  int n = RefineRMatrixEGrid(&er, idep, rs, rbs, rmx);
  if (n > 0) {
    if (_rmx_isave) {
      SaveRMatrixCE(rs, rbs, rmx, idep, fn, wt00);
    }
    InitDCFG(0, rmx[0].nts, rbs[0].nkappa, n, er, 0, NULL);
    RMatrixCEW(np, rbs, rmx, f, f1, fn, rs, n, er, m, mb, idep+1);
    free(er);
  }

  if (onke == 0) return 0;
  if (_stark_idx.n > 0) {
    for (i = 0; i < rmx[0].nsym; i++) {
      if (smx[i].isym < 0) continue;
      for (j = 0; j < _stark_idx.n; j++) {
	if (smx[i].nj[j] <= 0) continue;
	for (t = 0; t < smx[i].nk[j]; t++) {
	  free(smx[i].rp[j][t]);
	  free(smx[i].ip[j][t]);
	}
	free(smx[i].rp[j]);
	free(smx[i].ip[j]);
      }
      free(smx[i].rp);
      free(smx[i].ip);
      free(smx[i].jmin);
      free(smx[i].jmax);
      free(smx[i].nj);
      free(smx[i].nk);
      free(smx[i].sp);
    }
    free(smx);
    if (_stark_amp) {
      for (i = 0; i < _stark_idx.n; i++) {
	j = IdxGet(&rs->its, _stark_idx.d[i]);
	if (j >= 0) {
	  p = 2*(rmx[0].jts[j]+1);
	  p *= p;
	  for (q = 0; q < p; q++) {
	    for (t = 0; t < nkw; t++) {	      
	      free(ap[i][q][t]);
	    }
	    free(ap[i][q]);
	  }
	  free(ap[i]);
	}
      }
      free(ap);
    }
  }

  return 0;
}

int RMatrixConvert(char *ifn, char *ofn, int m) {
  int i;
  FILE *f0, *f1;
  RMATRIX rmx;

  if (m == 0) {
    ReadRMatrixBasis(ifn, &rbasis, 0);
    WriteRMatrixBasis(ofn, 1);
    ClearRMatrixBasis(&rbasis);
    return 0;
  } else if (m == 1) {
    ReadRMatrixBasis(ifn, &rbasis, 1);
    WriteRMatrixBasis(ofn, 0);
    ClearRMatrixBasis(&rbasis);
    return 0;
  } else if (m == 2) {
    f0 = fopen(ifn, "r");
    if (!f0) {
      printf("cannot open file %s\n", ifn);
      return -1;
    }
    f1 = fopen(ofn, "w");
    if (!f1) {
      printf("cannot open file %s\n", ofn);
      fclose(f0);
      return -1;
    }
    ReadRMatrixSurface(f0, &rmx, 0, 0);
    WriteRMatrixSurface(f1, NULL, NULL, 0, 1, &rmx, NULL);
    for (i = 0; i < rmx.nsym; i++) {
      ReadRMatrixSurface(f0, &rmx, 1, 0);
      if (SkipSym(rmx.isym)) continue;
      WriteRMatrixSurface(f1, NULL, NULL, 1, 1, &rmx, NULL);
    }
    ClearRMatrixSurface(&rmx);
    fclose(f0);
    fclose(f1);
    return 0;
  } else if (m == 3) {
    f0 = fopen(ifn, "r");
    if (!f0) {
      printf("cannot open file %s\n", ifn);
      return -1;
    }
    f1 = fopen(ofn, "w");
    if (!f1) {
      printf("cannot open file %s\n", ofn);
      fclose(f0);
      return -1;
    }
    ReadRMatrixSurface(f0, &rmx, 0, 1);
    WriteRMatrixSurface(f1, NULL, NULL, 0, 0, &rmx, NULL);
    for (i = 0; i < rmx.nsym; i++) {
      ReadRMatrixSurface(f0, &rmx, 1, 1);
      if (SkipSym(rmx.isym)) continue;
      WriteRMatrixSurface(f1, NULL, NULL, 1, 0, &rmx, NULL);
    }
    ClearRMatrixSurface(&rmx);
    fclose(f0);
    fclose(f1);
    return 0;
  }
  return 0;
}

void TestRMatrix(double e, int m, char *fn1, char *fn2, char *fn3) {
  FILE *f, *f1;
  RMATRIX rmx;
  RBASIS rbs;
  int i, j, p, k;
  double ei, a1, a2;

  e /= HARTREE_EV;
  f = fopen(fn2, "r");
  f1 = fopen(fn3, "w");
  ReadRMatrixBasis(fn1, &rbs, fmode);
  ReadRMatrixSurface(f, &rmx, 0, fmode);
  InitDCFG(rmx.mchan, rmx.nts, rbs.nkappa, 1, &e, 0, NULL);
  dcfg.ike = 0;
  for (k = 0; k < rmx.nsym; k++) {
    ReadRMatrixSurface(f, &rmx, 1, fmode);
    RMatrix(e, &rmx, &rbs, 1);
    if (rmx.isym == 0) {
      if (m >= 0) {
	if (dcfg.rgailitis > rbs.rb1) {      
	  GailitisExp(&rmx, &rbs, dcfg.rgailitis);
	  IntegrateExternal(&rmx, rbs.rb1, dcfg.rgailitis);
	} else {
	  GailitisExp(&rmx, &rbs, rbs.rb1);
	}
      } else {
	GailitisExp(&rmx, &rbs, rbs.rb1);
	PropogateExternal(&rmx, &rbs);
      }
      for (i = 0; i < rmx.nchan0; i++) {
	ei = (e-rmx.et[rmx.ilev[i]]+rmx.et0);
	fprintf(f1, "%2d %2d %2d %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E\n",
		rmx.isym, i, i, ei, a1, dcfg.fs0[i], 
		dcfg.gs0[i], dcfg.fc0[i], dcfg.gc0[i]);      	  
	for (j = 0; j < rmx.nchan0; j++) {
	  p = i*rmx.nchan0 + j; 
	  fprintf(f1, "%2d %2d %2d %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E\n",
		  rmx.isym, i, j, ei, a1, dcfg.fs[p], 
		  dcfg.gs[p], dcfg.fc[p], dcfg.gc[p]);
	}
      }
      fprintf(f1, "\n\n");
      break;
    }
  }
  ClearDCFG();
  fclose(f);
  fclose(f1);
}

void SetOptionRMatrix(char *s, char *sp, int ip, double dp) {
  if (0 == strcmp("rmatrix:nrefine", s)) {
    _nrefine = ip;
    return;
  }
  if (0 == strcmp("rmatrix:isave", s)) {
    _rmx_isave = ip;
    return;
  }
  if (0 == strcmp("rmatrix:mrefine", s)) {
    _mrefine = ip;
    return;
  }
  if (0 == strcmp("rmatrix:refiter", s)) {
    _refiter = ip;
    return;
  }  
  if (0 == strcmp("rmatrix:mineref", s)) {
    _mineref = dp/HARTREE_EV;
    return;
  }
  if (0 == strcmp("rmatrix:rrefine", s)) {
    _rrefine = dp;
    return;
  }
  if (0 == strcmp("rmatrix:fmaxe", s)) {
    _rmx_fmaxe = dp;
    return;
  }
  if (0 == strcmp("rmatrix:fmode", s)) {
    fmode = ip;
    return;
  }
  if (0 == strcmp("rmatrix:dk", s)) {
    _rmx_dk = ip;
    return;
  }
  if (0 == strcmp("rmatrix:acs", s)) {
    _rmx_acs = ip;
    return;
  }
  if (0 == strcmp("rmatrix:bfn", s)) {
    strncpy(_rmx_bfn, sp, 256);
    return;
  }
  if (0 == strcmp("rmatrix:efn", s)) {
    strncpy(_rmx_efn, sp, 256);
    return;
  }
  if (0 == strcmp("rmatrix:pj", s)) {
    _rmx_pj = ip;
    return;
  }
  if (0 == strcmp("rmatrix:mpj", s)) {
    _rmx_mpj = ip;
    return;
  }
  if (0 == strcmp("rmatrix:stark_ids", s)) {
    SetStarkIDs(sp);
    return;
  }
  if (0 == strcmp("rmatrix:stark_amp", s)) {
    _stark_amp = ip;
    return;
  }
  if (0 == strcmp("rmatrix:stark_pw", s)) {
    _stark_pw = ip;
    return;
  }
  if (0 == strcmp("rmatrix:stark_eadj", s)) {
    _stark_eadj = ip;
    return;
  }
  if (0 == strcmp("rmatrix:gailitis_expni", s)) {
    _gailitis_expni = ip;
    return;
  }
  if (0 == strcmp("rmatrix:gailitis_exprf", s)) {
    _gailitis_exprf = dp;
    return;
  }
  if (0 == strcmp("rmatrix:gailitis_exprt", s)) {
    _gailitis_exprt = dp;
    return;
  }
  if (0 == strcmp("rmatrix:minsp", s)) {
    _minsp = ip;
    return;
  }
}
