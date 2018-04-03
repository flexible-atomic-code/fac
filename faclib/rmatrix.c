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
static DCFG dcfg;
static int nbatch;
static int fmode;

void RMatrixNBatch(int n) {
  if (n > 0) {
    nbatch = n; 
  } else {
    nbatch = 100;
  }
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

  RMatrixNBatch(0);
  fmode = 0;
  
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
    nr = fread(&(rbs->kmax), sizeof(int), 1, f);
    nr = fread(&(rbs->nbk), sizeof(int), 1, f);
    nr = fread(&(rbs->nkappa), sizeof(int), 1, f);
    nr = fread(&(rbs->nbuttle), sizeof(int), 1, f);
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
	       &(rbs->basis[i][n]), &kappa, &(rbs->bnode[i][n]), &(rbs->ek[i][n]),
	       &(rbs->w0[i][n]), &(rbs->w1[i][n]));
      }
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
    nr = fwrite(&(rbasis.kmax), sizeof(int), 1, f);
    nr = fwrite(&(rbasis.nbk), sizeof(int), 1, f);
    nr = fwrite(&(rbasis.nkappa), sizeof(int), 1, f);
    nr = fwrite(&(rbasis.nbuttle), sizeof(int), 1, f);
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
	    rbasis.kmax, rbasis.nbk, rbasis.nkappa, rbasis.nbuttle);
    for (i = 0; i < rbasis.nkappa; i++) {
      ka = KappaFromIndex(i);
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
    SetBoundary(nmax, r1, b);
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
	r1 = 2.0*pot->rad[pot->ib1] - pot->rad[pot->ib];
      } else {
	r1 = 2.0*pot->rad[pot->ib];
      }
    }
    pot->ib = pot->ib1;
    SetBoundary(-100, r1, r0);
    rbasis.ib0 = pot->ib;
    rbasis.ib1 = pot->ib1;
  }
}

void ExtrapolateButtle(RBASIS *rbs, int t, int m, double *e, 
		       double *r0, double *r1, double *r2) {
  int nb, i, n, k;
  double cp0[2], cp1[2], cp2[3];
  double xb1[NBFIT], xb2[NBFIT], yb0[NBFIT], yb1[NBFIT], eb[NBFIT];
  double a0, a1, b, ep;
  
  nb = rbs->nbk;
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
  int nkb0, nkb1, n, n0, nmax, kb;
  double e0, e1, ep, a0, a1, b, rb0, rb1, bb, c0, c1;
  double r01, r10, r0, r1, r2, p0, p1, q0, q1, x0, x1;
  ORBITAL *orb, orbf;
  POTENTIAL *pot;

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
  rbasis.kmax = kmax;
  rbasis.nbk = nb;
  rbasis.nkappa = (2*kmax + 1);
  rbasis.nbuttle = nb;
  rbasis.basis = malloc(sizeof(int *)*rbasis.nkappa);
  rbasis.bnode = malloc(sizeof(int *)*rbasis.nkappa);
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

  nkb0 = GetNumOrbitals();
  e1 = 0.0;
  e0 = 0.0;
  for (k = 0; k <= kmax; k++) {
    n0 = k+1;
    if (rbasis.ib0 == 0 && n0 <= nmax) n0 = nmax + 1;
    k2 = 2*k;
    t = k2 - 1;
    for (j = k2-1; j <= k2+1; j += 2, t++) {     
      if (j < 0) continue;
      kappa = GetKappaFromJL(j, k2);
      printf("Basis: %3d\n", kappa);
      fflush(stdout);
      for (in = 0; in < nb; in++) {
	if (rbasis.ib0 == 0) {
	  n = in + n0;
	} else {
	  n = -(in + n0);
	}
	rbasis.basis[t][in] = OrbitalIndex(n, kappa, 0.0);
	rbasis.bnode[t][in] = n;
	orb = GetOrbital(rbasis.basis[t][in]);
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

      ExtrapolateButtle(&rbasis, t, rbasis.nbuttle, rbasis.ebuttle[t], 
			rbasis.cbuttle[3][t], rbasis.cbuttle[2][t],
			rbasis.cbuttle[4][t]);
      for (i = 0; i < rbasis.nbuttle; i++) {
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
  WriteRMatrixBasis(fn, fmode);
  if (nts > 0 && rbasis.ib0 == 0) {
    nkb1 = GetNumOrbitals();
    PrepSlater(0, nkb0-1, nkb0, nkb1-1, 0, nkb0-1, nkb0, nkb1-1);
    PrepSlater(0, nkb0-1, 0, nkb0-1, nkb0, nkb1-1, nkb0, nkb1-1);
  }
  return 0;
}

int IndexFromKappa(int ka) {
  int j, k;

  GetJLFromKappa(ka, &j, &k);
  if (j > k) return k;
  else return k-1;
}

int KappaFromIndex(int k) {
  int j;

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
  
  nlevs = GetNumLevels();
  if (nts > 0) free(ts);
  if (ncs > 0) free(cs);
  ts = malloc(sizeof(int)*nlevs);
  cs = malloc(sizeof(int)*nlevs);
  nts = 0;
  ncs = 0;
  for (i = 0; i < nlevs; i++) {
    lev = GetLevel(i);
    m = lev->pb;
    sym = GetSymmetry(lev->pj);
    s = ArrayGet(&(sym->states), m);
    if (InGroups(s->kgroup, nt, kt)) {
      ts[nts] = i;
      nts++;
    } else if (InGroups(s->kgroup, nc, kc)) {
      cs[ncs] = i;
      ncs++;
    }
  }
  if (ncs == 0) {
    free(cs);
  }
  
  AngularFrozen(nts, ts, ncs, cs);
  ntg = nt;
  tg = malloc(sizeof(int)*ntg);
  memcpy(tg, kt, sizeof(int)*ntg);
  ncg = nc;
  cg = malloc(sizeof(int)*ncg);
  memcpy(cg, kc, sizeof(int)*ncg);
}

void ClearRMatrixSurface(RMATRIX *rmx) {
  int i;

  if (rmx->nts > 0) {
    free(rmx->et);
    free(rmx->ts);
    free(rmx->jts);
  } 
  if (rmx->ncs > 0) {
    free(rmx->ec);
    free(rmx->cs);
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
      free(rmx->rmatrix[i]);
    }
  }
  if (rmx->nlam > 0) free(rmx->aij);
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
      rmx->jts = malloc(sizeof(int)*nts);
      rmx->cs = malloc(sizeof(int)*ncs);
      rmx->jcs = malloc(sizeof(int)*ncs);
      rmx->et0 = 1E30;
      for (i = 0; i < nts; i++) {
	ierr = fread(&(rmx->ts[i]), sizeof(int), 1, f);
	ierr = fread(&(rmx->jts[i]), sizeof(int), 1, f);
	ierr = fread(&(rmx->et[i]), sizeof(double), 1, f);
	if (rmx->et[i] < rmx->et0) rmx->et0 = rmx->et[i];
      }
      for (i = 0; i < ncs; i++) {
	ierr = fread(&(rmx->cs[i]), sizeof(int), 1, f);
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
      rmx->rmatrix[i] = malloc(sizeof(double)*nchan0*nchan0);
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
      rmx->jts = malloc(sizeof(int)*nts);
      rmx->cs = malloc(sizeof(int)*ncs);
      rmx->jcs = malloc(sizeof(int)*ncs);
      rmx->et0 = 1E30;
      for (i = 0; i < nts; i++) {
	ierr = fscanf(f, "T %d %d %lf\n", 
		      &(rmx->ts[i]), &(rmx->jts[i]), &(rmx->et[i]));
	if (rmx->et[i] < rmx->et0) rmx->et0 = rmx->et[i];
      }
      for (i = 0; i < ncs; i++) {
	ierr = fscanf(f, "C %d %d %lf\n", &k, &k1, &a);
	rmx->cs[i] = k;
	rmx->jcs[i] = k1;
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
      rmx->rmatrix[i] = malloc(sizeof(double)*nchan0*nchan0);
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
	  DecodePJ(lev->pj, NULL, &k);
	  nr = fwrite(&ilev, sizeof(int), 1, f);
	  nr = fwrite(&k, sizeof(int), 1, f);
	  nr = fwrite(&(lev->energy), sizeof(double), 1, f);
	}
	for (t = 0; t < ncs; t++) {
	  ilev = cs[t];
	  lev = GetLevel(ilev);
	  DecodePJ(lev->pj, NULL, &k);
	  nr = fwrite(&ilev, sizeof(int), 1, f);
	  nr = fwrite(&k, sizeof(int), 1, f);
	  nr = fwrite(&(lev->energy), sizeof(double), 1, f);	  
	}
      } else {
	fprintf(f, "%3d %3d %3d %3d %3d %10.3E %2d\n",
		m, m, nts, ncs, rbasis.nkappa, z, dcfg.nmultipoles);
	for (t = 0; t < nts; t++) {
	  ilev = ts[t];
	  lev = GetLevel(ilev);
	  DecodePJ(lev->pj, NULL, &k);
	  fprintf(f, "T %3d %3d %15.8E\n", ilev, k, lev->energy);
	}
	for (t = 0; t < ncs; t++) {
	  ilev = cs[t];
	  lev = GetLevel(ilev);
	  DecodePJ(lev->pj, NULL, &k);
	  fprintf(f, "C %3d %3d %15.8E\n", ilev, k, lev->energy);
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
	  nr = fwrite(&(rmx->jts[t]), sizeof(int), 1, f);
	  nr = fwrite(&(rmx->et[t]), sizeof(double), 1, f);
	}
	for (t = 0; t < ncs; t++) {
	  nr = fwrite(&(rmx->cs[t]), sizeof(int), 1, f);
	  nr = fwrite(&(rmx->jcs[t]), sizeof(int), 1, f);
	  nr = fwrite(&(rmx->ec[t]), sizeof(double), 1, f);
	}
      } else {
	fprintf(f, "%3d %3d %3d %3d %3d %10.3E %2d\n",
		rmx->nsym, rmx->mchan, nts, ncs, rmx->nkappa, z, rmx->nlam);
	for (t = 0; t < nts; t++) {
	  ilev = rmx->ts[t];
	  k = rmx->jts[t];
	  fprintf(f, "T %3d %3d %15.8E\n", ilev, k, rmx->et[t]);
	}
	for (t = 0; t < ncs; t++) {
	  ilev = rmx->cs[t];
	  k = rmx->jcs[t];
	  fprintf(f, "C %3d %3d %15.8E\n", ilev, k, rmx->ec[t]);
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
	  ka = KappaFromIndex(k);
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
	  ka = KappaFromIndex(k);
	  for (jc = 0; jc <= ic; jc++) {
	    if (wik1[jc] == NULL) continue;
	    i1 = jc/rbasis.nkappa;
	    ilev1 = ts[i1];
	    k1 = jc%rbasis.nkappa;
	    ka1 = KappaFromIndex(k1);
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
      for (t = 0; t < h->dim; t++) {
	fprintf(f, "%6d %17.10E\n", t, h->mixing[t]);
      }
      for (ic = 0; ic < nchan; ic++) {
	if (wik1[ic]) {
	  i = ic/rbasis.nkappa;
	  k = ic%rbasis.nkappa;
	  ilev = ts[i];
	  lev = GetLevel(ilev);
	  ka = KappaFromIndex(k);
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
	  ka = KappaFromIndex(k);
	  for (jc = 0; jc <= ic; jc++) {
	    if (wik1[jc] == NULL) continue;
	    i1 = jc/rbasis.nkappa;
	    ilev1 = ts[i1];
	    k1 = jc%rbasis.nkappa;
	    ka1 = KappaFromIndex(k1);
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
		    ic, t, wik0[ic][t], wik1[ic][t]);
	  }
	  free(wik0[ic]);
	  free(wik1[ic]);
	  wik0[ic] = NULL;
	  wik1[ic] = NULL;
	}
      }
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
      for (t = 0; t < rmx->ndim; t++) {
	fprintf(f, "%6d %17.10E\n", t, rmx->ek[t]);
      }
      for (ic = 0; ic < nchan0; ic++) {
	i = rmx->ilev[ic];
	ka = rmx->kappa[ic];
	k = IndexFromKappa(ka);
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
		  i, t, rmx->w0[i][t], rmx->w1[i][t]);
	}	
      }
    }
  }

  return nchan0;
}

int RMatrixSurface(char *fn) {
  int i, m, k, t, kb, ilev, ika, nsym, nchm, nchan0;
  int j0, p0, j1, k1, j, p, jmin, jmax, nchan, q, ic;
  double *ek, *mix, **wik0, **wik1;
  SYMMETRY *sym;
  STATE *s;
  LEVEL *lev;
  ORBITAL *orb;
  HAMILTON *h;
  FILE *f;

  f = fopen(fn, "w");
  if (f == NULL) return -1;

  for (i = 0; i < ncs; i++) {
    lev = GetLevel(cs[i]);
    DecodePJ(lev->pj, &p0, &j0);
    AddStateToSymmetry(-(cs[i]+1), -1, j0, p0, j0);
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
	GetJLFromKappa(orb->kappa, &j1, &k1);
	p = p0 + k1/2;
	jmin = abs(j0 - j1);
	jmax = j0 + j1;
	for (j = jmin; j <= jmax; j += 2) {
	  AddStateToSymmetry(-(ts[i]+1), kb, j, p, j);
	}
      }
    }
  }

  WriteRMatrixSurface(f, NULL, NULL, 0, fmode, NULL, NULL);

  nchan = nts*rbasis.nkappa;
  printf("%d %d %d %d\n", nchan, nts, rbasis.nkappa, ncs);
  wik0 = malloc(sizeof(double *)*nchan);
  wik1 = malloc(sizeof(double *)*nchan);
  for (t = 0; t < nchan; t++) {
    wik0[t] = NULL;
    wik1[t] = NULL;
  }
  nsym = 0;
  nchm = 0;
  for (i = 0; i < MAX_SYMMETRIES; i++) {
    h = GetHamilton(i);
    if (rbasis.ib0 == 0) {
      k = ConstructHamiltonFrozen(i, ntg, tg, 0, ncg, cg);
    } else {
      k = ConstructHamiltonFrozen(i, ntg, tg, -1, 0, NULL);
    }
    if (k < 0) continue;
    DiagnolizeHamilton(h);
    ek = h->mixing;
    mix = h->mixing+h->dim;
    sym = GetSymmetry(i);
    printf("sym: %d %d\n", i, h->dim);
    fflush(stdout);
    for (t = 0; t < h->dim; t++) {
      for (q = 0; q < h->dim; q++) {
	k = h->basis[q];
	s = ArrayGet(&(sym->states), k);
	k = -(s->kgroup+1);
	ilev = IBisect(k, nts, ts);
	if (ilev >= 0) {
	  kb = s->kcfg;
	  orb = GetOrbital(kb);
	  ika = IndexFromKappa(orb->kappa);
	  ic = ilev*rbasis.nkappa + ika;
	  if (wik1[ic] == NULL) {
	    wik1[ic] = malloc(sizeof(double)*h->dim);
	    for (p = 0; p < h->dim; p++) {
	      wik1[ic][p] = 0.0;
	    }
	  }
	  wik1[ic][t] += mix[q]*WLarge(orb)[rbasis.ib1];
	  if (wik0[ic] == NULL) {
	    wik0[ic] = malloc(sizeof(double)*h->dim);
	    for (p = 0; p < h->dim; p++) {
	      wik0[ic][p] = 0.0;
	    }
	  }
	  wik0[ic][t] += mix[q]*WLarge(orb)[rbasis.ib0];
	}
      }
      mix += h->dim;
    }
    nchan0 = WriteRMatrixSurface(f, wik0, wik1, 1, fmode, NULL, h);
    if (nchan0 > nchm) nchm = nchan0;
    nsym++;
  }
  
  fseek(f, 0, SEEK_SET);
  if (fmode == 0) {
    fwrite(&nsym, sizeof(int), 1, f);
    fwrite(&nchm, sizeof(int), 1, f);
  } else {
    fprintf(f, "%3d %3d", nsym, nchm);
  }

  free(wik0);
  free(wik1);
  fclose(f);
  return 0;
}

int RMatrix(double e, RMATRIX *rmx, RBASIS *rbs, int m) {
  int i, j, p, q, nb, k;
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
      rmx->rmatrix[0][p] = a11;
      rmx->rmatrix[1][p] = a00;
      rmx->rmatrix[2][p] = a01;
      if (i != j) {	
	p = j*rmx->nchan0 + i;
	rmx->rmatrix[0][p] = a11;
	rmx->rmatrix[1][p] = a00;
	rmx->rmatrix[2][p] = a10;
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
	  rmx->rmatrix[0][p] += b;	
	  if (rbs->ib0 > 0) {
	    UVIP3P(3, nb, x, y1, 1, &de, &b);
	    rmx->rmatrix[1][p] += b;
	    UVIP3P(3, nb, x, y2, 1, &de, &b);
	    rmx->rmatrix[2][p] += b;
	  }
	} else {
	  ExtrapolateButtle(rbs, k, 1, &de, &a00, &a01, &a11);
	  rmx->rmatrix[0][p] += a01;
	  if (rbs->ib0 > 0) {
	    rmx->rmatrix[1][p] += a00;
	    rmx->rmatrix[2][p] += a11;
	  }
	}
      }
    }
  }

  rmx->energy = e;
  return 0;
}    

int RMatrixPropogate(double *r0, double *r1, RMATRIX *rmx1) {
  int i, j, k, m, p, q, r;
  double a, b;
  int *iwork = dcfg.iwork;

  for (i = 0; i < rmx1->nchan0; i++) {
    for (j = 0; j <= i; j++) {
      p = i*rmx1->nchan0 + j;
      r0[p] += rmx1->rmatrix[1][p];
      if (i != j) {
	q = j*rmx1->nchan0 + i;
	r0[q] += rmx1->rmatrix[1][q];
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
	  b *= rmx1->rmatrix[2][q]*rmx1->rmatrix[2][p];
	  a += b;
	}
      }
      p = i*rmx1->nchan0 + j;
      r1[p] = rmx1->rmatrix[0][p] - a;
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
  a = rmx0->rmatrix[1];
  b = rmx0->rmatrix[2];
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
  
  c = rmx0->rmatrix[2];
  a = rmx0->rmatrix[0];
  b = rmx0->rmatrix[1];
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
  memcpy(a, rmx->rmatrix[2], sizeof(double)*nch2);
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
	b += rmx->rmatrix[0][jm]*dcfg.gs[mi];
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
	b += rmx->rmatrix[2][jm]*dcfg.gs[mi];
	b -= rmx->rmatrix[1][jm]*y[m];
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
	b += rmx->rmatrix[0][jm]*dcfg.gc[mi];
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
	b += rmx->rmatrix[2][jm]*dcfg.gc[mi];
	b -= rmx->rmatrix[1][jm]*y[m];
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

int GailitisExp(RMATRIX *rmx, double r) {
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
    e[i] = rmx->energy - (rmx->et[rmx->ilev[i]] - rmx->et0);    
    if (nlam > 0 && dcfg.ngailitis > 1) { 
      p2[i] = 2.0*e[i]*(1.0+0.5*FINE_STRUCTURE_CONST2*e[i]);
      p[i] = sqrt(fabs(p2[i]));
    }
    ierr = 0;
    DCOUL(rmx->z, e[i], rmx->kappa[i], r, &t1, &c1, &t2, &c2, &ierr);
    if (e[i] > 0) {
      dcfg.fs0[i] = t1;
      dcfg.gs0[i] = c1;
      dcfg.fc0[i] = t2;
      dcfg.gc0[i] = c2;
      nop++;
    } else {
      dcfg.fs0[i] = 0.0;
      dcfg.gs0[i] = 0.0;
      dcfg.fc0[i] = t1;
      dcfg.gc0[i] = c1;
    }
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

static void InitDCFG(int nw) {
  int nw2, ns, n;

  nw2 = nw*nw;
  ns = nw*dcfg.ngailitis;
  dcfg.lrw = 20 + 16*nw*2;
  if (dcfg.lrw < nw2*2) dcfg.lrw = nw2*2;
  dcfg.liw = 20+nw;
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
}

static void ClearDCFG(void) {
  free(dcfg.dwork);
  free(dcfg.iwork);
}

int RMatrixCE(char *fn, int np, char *bfn[], char *rfn[], 
	      double emin, double emax, double de, int m, int mb) {
  RBASIS *rbs;
  RMATRIX *rmx;
  FILE **f, *f1;
  int i, j, k, t, p, q, h, n, i0, ns, *iwork;
  double *e0, *e, et, **s, ***sp, *r0, *r1, x;
  int pp, jj, its0, its1, st0, ka0, ka1, mka0, mka1, npw;
  int nke, npe;

  emin /= HARTREE_EV;
  emax /= HARTREE_EV;
  de /= HARTREE_EV;

  n = (emax - emin)/de + 1;
  e0 = malloc(sizeof(double)*n);
  e0[0] = emin;
  for (i = 1; i < n; i++) {
    e0[i] = e0[i-1] + de;
  }

  f1 = fopen(fn, "w");
  rbs = malloc(sizeof(RBASIS)*np);
  f = malloc(sizeof(FILE *)*np);
  rmx = malloc(sizeof(RMATRIX)*np);
  
  for (i = 0; i < np; i++) {
    ReadRMatrixBasis(bfn[i], &(rbs[i]), fmode);
    f[i] = fopen(rfn[i], "r");
    if (f[i] == NULL) return -1;
    ReadRMatrixSurface(f[i], &(rmx[i]), 0, fmode);
  }
  
  for (t = 0; t < rmx[0].nts; t++) {
    et = rmx[0].et[t] - rmx[0].et0;
    i0 = (et - emin)/de;
    i = i0-10;
    if (i < 0) i = 0;
    for (; i <= i0+10 && i < n; i++) {
      if (fabs(et-e0[i]) < 1E-5) {
	if (et < e0[i]) e0[i] += 1E-5;
	else e0[i] -= 1E-5;
      }
    }
  }

  InitDCFG(rmx[0].mchan);

  e = e0;
  npe = 0;
  while (npe < n) {
    if (npe + nbatch > n) {
      nke = n - npe;
    } else {
      nke = nbatch;
    }
    for (i = 0; i < np; i++) {
      ClearRMatrixSurface(&(rmx[i]));      
      fseek(f[i], 0, SEEK_SET);
      ReadRMatrixSurface(f[i], &(rmx[i]), 0, fmode);
    }
    ns = rmx[0].nts*(rmx[0].nts+1)/2;
    s = malloc(sizeof(double *)*ns);
    for (i = 0; i < ns; i++) {
      s[i] = malloc(sizeof(double)*nke);
      for (k = 0; k < nke; k++) {
	s[i][k] = 0.0;
      }
    }
    npw = 0;
    if (m & 1) {
      npw = rbs[0].kmax+1;
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
    fprintf(f1, "### %6d %6d %4d %12.5E %12.5E\n", 
	    n, npe, nke, e[0]*HARTREE_EV, e[nke-1]*HARTREE_EV);
    fflush(f1);
    for (i = 0; i < rmx[0].nsym;  i++) {
      for (j = 0; j < np; j++) {
	ReadRMatrixSurface(f[j], &(rmx[j]), 1, fmode);
      }
      DecodePJ(rmx[0].isym, &pp, &jj);
      printf("sym: %d %d %d\n", rmx[0].isym, pp, jj);
      fflush(stdout);
      for (k = 0; k < nke; k++) {
	for (j = 0; j < np; j++) {
	  RMatrix(e[k], &(rmx[j]), &(rbs[j]), mb);
	}
	r0 = rmx[0].rmatrix[0];
	r1 = rmx[0].rmatrix[1];
	if (dcfg.rgailitis > rbs[np-1].rb1) {
	  GailitisExp(&(rmx[0]), dcfg.rgailitis);
	  IntegrateExternal(&(rmx[0]), rbs[np-1].rb1, dcfg.rgailitis);
	} else {
	  GailitisExp(&(rmx[0]), rbs[np-1].rb1);
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
	  mka1 = GetLFromKappa(ka1);
	  mka1 /= 2;
	  for (q = 0; q <= p; q++) {
	    h = q*rmx[0].nchan0 + p;
	    x = r0[h]*r0[h] + r1[h]*r1[h];
	    x *= 0.5*(jj+1.0);
	    its0 = rmx[0].ilev[q];
	    ka0 = rmx[0].kappa[q];
	    mka0 = GetLFromKappa(ka0);
	    mka0 /= 2;
	    s[st0+its0][k] += x;
	    if (m & 1) {
	      sp[st0+its0][k][mka0] += x;
	    }
	    if (m & 2) {
	      et = (e[k]-rmx[0].et[its0]+rmx[0].et0)*HARTREE_EV;
	      fprintf(f1, 
		      "%3d %4d %4d %4d %4d %12.5E %12.5E %12.5E %12.5E %12.5E\n",
		      rmx[0].isym, its0, ka0, its1, ka1, et, 
		      r0[h], r1[h], 
		      x, rmx[0].rmatrix[2][h]);
	    }
	  }
	}
      }
    }

    if (m & 2) {
      fprintf(f1, "\n\n");
    }
    for (its0 = 0; its0 < rmx[0].nts; its0++) {
      for (its1 = its0; its1 < rmx[0].nts; its1++) {
	st0 = its1*(its1+1)/2;
	fprintf(f1, "# %3d %3d %3d %3d %12.5E %6d %3d\n", 
		rmx[0].ts[its0], rmx[0].jts[its0],
		rmx[0].ts[its1], rmx[0].jts[its1], 
		(rmx[0].et[its1]-rmx[0].et[its0])*HARTREE_EV, 
		nke, npw);
	for (k = 0; k < nke; k++) {
	  et = (e[k]-rmx[0].et[its0]+rmx[0].et0)*HARTREE_EV;
	  if (et < 0.0) continue;
	  fprintf(f1, "%15.8E %15.8E\n", et, s[st0+its0][k]);
	}
	fprintf(f1, "\n\n");
	if (m & 1) {
	  for (k = 0; k < nke; k++) {
	    et = (e[k]-rmx[0].et[its0]+rmx[0].et0)*HARTREE_EV;
	    if (et < 0.0) continue;
	    for (j = 0; j < npw; j++) {
	      fprintf(f1, "%15.8E %3d %15.8E\n", et, j, sp[st0+its0][k][j]);
	    }
	  }
	  fprintf(f1, "\n\n");
	}
      }
    }
    
    for (i = 0; i < ns; i++) {
      free(s[i]);
      if (m & 1) {
	for (j = 0; j < nke; j++) {
	  free(sp[i][j]);
	}
	free(sp[i]);
      }
    }
    free(s);
    if (m & 1) {
      free(sp);
    }
    e += nke;
    npe += nke;
    fflush(f1);
  }
      
  for (i = 0; i < np; i++) {
    fclose(f[i]);
    ClearRMatrixBasis(&(rbs[i]));
    ClearRMatrixSurface(&(rmx[i]));
  }
  free(f);
  free(rbs);
  free(rmx);
  fclose(f1);
  free(e0);
  ClearDCFG();
  
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
      printf("sym: %3d\n", i);
      ReadRMatrixSurface(f0, &rmx, 1, 0);
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
      printf("sym: %3d\n", i);
      ReadRMatrixSurface(f0, &rmx, 1, 1);
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
  InitDCFG(rmx.mchan);
  for (k = 0; k < rmx.nsym; k++) {
    ReadRMatrixSurface(f, &rmx, 1, fmode);
    RMatrix(e, &rmx, &rbs, 1);
    /*
    if (rmx.isym == 20) {
      GailitisExp(&rmx, dcfg.rgailitis);
      a2 = (dcfg.rgailitis-rbs.rb1)/3000;
      a1 = dcfg.rgailitis;
      while (a1-a2 - rbs.rb1 > -0.1*a2) {
	IntegrateExternal(&rmx, a1-a2, a1);
	a1 -= a2;
	for (i = 0; i < 1; i++) {
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
      }
      break;
    }  
    */
    if (rmx.isym == 0) {
      if (m >= 0) {
	if (dcfg.rgailitis > rbs.rb1) {      
	  GailitisExp(&rmx, dcfg.rgailitis);
	  IntegrateExternal(&rmx, rbs.rb1, dcfg.rgailitis);
	} else {
	  GailitisExp(&rmx, rbs.rb1);
	}
      } else {
	GailitisExp(&rmx, rbs.rb1);
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
    /*
    for (i = 0; i < rmx.nchan0; i++) {
      for (j = 0; j < rmx.nchan0; j++) {
	p = i*rmx.nchan0+j;
	if (rmx.nlam >= 2) {
	  a1 = rmx.aij[0][p];
	  a2 = rmx.aij[1][p];
	} else {
	  a1 = 0.0;
	  a2 = 0.0;
	}
	fprintf(f1, "%3d %2d %2d %15.8E %15.8E %15.8E %15.8E %15.8E\n",
		rmx.isym, i, j, a1, a2,
		rmx.rmatrix[0][p], rmx.rmatrix[1][p], rmx.rmatrix[2][p]);
      }
    }
    fprintf(f1, "\n\n");
    */
  }
  ClearDCFG();
  fclose(f);
  fclose(f1);
}
