#include "rmatrix.h"

static RBASIS rbasis;
static int ntg, *tg, nts, *ts;
static int ncg, *cg, ncs, *cs;

int InitRMatrix(void) {
  ntg = 0;
  ncg = 0;
  nts = 0;
  ncs = 0;
  rbasis.nkappa = 0;
  
  return 0;
}

void ClearRMatrixBasis(void) {
  int i, k;

  for (i = 0; i < rbasis.nkappa; i++) {
    for (k = 0; k < NBTERMS; k++) {
      free(rbasis.cbuttle[k][i]);
    }
    free(rbasis.ebuttle[i]);
    free(rbasis.basis[i]);
    free(rbasis.w0[i]);
    free(rbasis.w1[i]);
    free(rbasis.ek[i]);
  }
  if (rbasis.nkappa > 0) {
    for (k = 0; k < NBTERMS; k++) {
      free(rbasis.cbuttle[k]);
    }
    free(rbasis.ebuttle);
    free(rbasis.basis);
    free(rbasis.w0);
    free(rbasis.w1);
    free(rbasis.ek);
  }
  rbasis.nkappa = 0;
}
  
void WriteRMatrixBasis(char *fn) {
  FILE *f;
  int i, n, k;
  ORBITAL *orb;

  f = fopen(fn, "w");
  if (f == NULL) return;

  fprintf(f, "%4d %15.8E %4d %15.8E %15.8E\n", 
	  rbasis.ib0, rbasis.rb0, rbasis.ib1, rbasis.rb1, rbasis.bqp);
  fprintf(f, "%2d %3d %3d %3d\n",
	  rbasis.kmax, rbasis.nbk, rbasis.nkappa, rbasis.nbuttle);
  for (i = 0; i < rbasis.nkappa; i++) {
    for (n = 0; n < rbasis.nbk; n++) {
      orb = GetOrbital(rbasis.basis[i][n]);
      fprintf(f, "%4d %3d %3d %15.8E %15.8E %15.8E\n",
	      rbasis.basis[i][n], orb->kappa, orb->n, rbasis.ek[i][n],
	      rbasis.w0[i][n], rbasis.w1[i][n]);
    }
    for (n = 0; n < rbasis.nbuttle; n++) {
      fprintf(f, "%4d %15.8E",
	      orb->kappa, rbasis.ebuttle[i][n]);
      for (k = 0; k < NBTERMS; k++) {
	fprintf(f, " %15.8E", rbasis.cbuttle[k][i][n]);
      }
      fprintf(f, "\n");
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
    pot->bqp = b;
    pot->ib = pot->ib1;
    SetBoundary(-100, r1, r0);
    rbasis.ib0 = pot->ib;
    rbasis.ib1 = pot->ib1;
  }
}

int RMatrixBasis(char *fn, int kmax, int nb) {
  int k, i, j, k2, ib0, ib1, kappa, t, in;
  int nkb0, nkb1, n, n0, nmax, kb;
  double e0, e1, ep, a0, a1, b, rb0, rb1, bb, c0, c1;
  double r01, r10, r0, r1, r2, p0, p1, q0, q1, x0, x1;
  double xb1[NBFIT], xb2[NBFIT];
  double yb0[NBFIT], yb1[NBFIT], eb[NBFIT];
  double cp0[2], cp1[2], cp2[3];
  ORBITAL *orb, orbf;
  POTENTIAL *pot;

  ClearRMatrixBasis();
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
  rbasis.nbuttle = (nb-1);
  rbasis.basis = malloc(sizeof(int *)*rbasis.nkappa);
  rbasis.ebuttle = malloc(sizeof(double *)*rbasis.nkappa);
  for (kb = 0; kb < NBTERMS; kb++) {
    rbasis.cbuttle[kb] = malloc(sizeof(double *)*rbasis.nkappa);
  }
  rbasis.w0 = malloc(sizeof(double *)*rbasis.nkappa);
  rbasis.w1 = malloc(sizeof(double *)*rbasis.nkappa);
  rbasis.ek = malloc(sizeof(double *)*rbasis.nkappa);
  for (i = 0; i < rbasis.nkappa; i++) {
    rbasis.basis[i] = malloc(sizeof(int)*rbasis.nbk);
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
      for (in = 0; in < nb; in++) {
	if (rbasis.ib0 == 0) {
	  n = in + n0;
	} else {
	  n = -(in + n0);
	}
	rbasis.basis[t][in] = OrbitalIndex(n, kappa, 0.0);	
	orb = GetOrbital(rbasis.basis[t][in]);
	rbasis.ek[t][in] = orb->energy;
	rbasis.w1[t][in] = WLarge(orb)[ib1];
	rbasis.w0[t][in] = WLarge(orb)[ib0];
	e0 = e1;
	e1 = orb->energy;
	if (in > 0) {
	  /*
	  rbasis.ebuttle[t][in*2-2] = e0 + (e1 - e0)*0.25;
	  rbasis.ebuttle[t][in*2-1] = e1 - (e1 - e0)*0.25;
	  */
	  rbasis.ebuttle[t][in-1] = 0.5*(e1 + e0);
	}
      }
      for (in = 0; in < NBFIT; in++) {
	n = nb - NBFIT + in;
	xb1[in] = n+n0;
	xb2[in] = 1.0/xb1[in];
	yb0[in] = fabs(rbasis.w0[t][n]);
	yb1[in] = fabs(rbasis.w1[t][n]);
	eb[in] = rbasis.ek[t][n];
      }
      if (rbasis.ib0 > 0) {
	PolyFit(2, cp0, NBFIT, xb2, yb0);
      }
      PolyFit(2, cp1, NBFIT, xb2, yb1);
      PolyFit(3, cp2, NBFIT, xb1, eb);
      for (i = 0; i < rbasis.nbuttle; i++) {
	r0 = 0.0;
	r1 = 0.0;
	r2 = 0.0;
	for (in = 0; in < 25000; in++) {
	  n = in + nb + n0;
	  b = 1.0/n;
	  if (rbasis.ib0 > 0) {
	    a0 = cp0[0] + cp0[1]*b;
	  }
	  a1 = cp1[0] + cp1[1]*b;
	  ep = cp2[0] + cp2[1]*n + cp2[2]*n*n;
	  if (rbasis.ib0 > 0) {
	    r0 += 0.5*a0*a0/(ep - rbasis.ebuttle[t][i]);
	    b = 0.5*a0*a1;
	    if (IsOdd(in)) b = -b;
	    b /= ep - rbasis.ebuttle[t][i];
	    r2 += b;
	  }
	  r1 += 0.5*a1*a1/(ep - rbasis.ebuttle[t][i]);
	}
	n++;
	if (rbasis.ib0 > 0) {
	  a0 = 0.5*cp0[0]*cp0[0]/cp2[2];
	  r0 += a0/n;
	  a0 = 0.5*cp0[0]*cp1[0]/cp2[2];
	  r2 += 0.5*(a0/n - a0/(n+1.0));
	  if (rbasis.w0[t][nb-1]*rbasis.w1[t][nb-1] > 0) r2 = -r2;
	}
	a0 = 0.5*cp1[0]*cp1[0]/cp2[2];
	r1 += a0/n;
	rbasis.cbuttle[2][t][i] = r1;
	rbasis.cbuttle[3][t][i] = r0;
	rbasis.cbuttle[4][t][i] = r2;
	orbf.n = 1000000;
	orbf.kappa = kappa;
	orbf.energy = rbasis.ebuttle[t][i];
	RadialFreeInner(&orbf, pot);
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
  if (rbasis.ib0 == 0) {
    nkb1 = GetNumOrbitals();
    PrepSlater(0, nkb0-1, nkb0, nkb1-1, 0, nkb0-1, nkb0, nkb1-1);
    PrepSlater(0, nkb0-1, 0, nkb0-1, nkb0, nkb1-1, nkb0, nkb1-1);
  }
  WriteRMatrixBasis(fn);
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

void WriteRMatrixSurface(FILE *f, double **wik0, double **wik1) {
  HAMILTON *h;
  int nchan, t, ic, nchan0, i, k, ilev, ka;
  int p, j;
  LEVEL *lev;
  
  h = GetHamilton();
  nchan = rbasis.nkappa * nts;
  DecodePJ(h->pj, &p, &j);
  fprintf(f, "%3d %3d %3d %3d %2d %2d %6d\n",
	  nts, ncs, rbasis.nkappa, h->pj, p, j, h->dim);
  for (t = 0; t < h->dim; t++) {
    fprintf(f, "%6d %17.10E\n", t, h->mixing[t]);
  }
  nchan0 = 0;
  for (ic = 0; ic < nchan; ic++) {
    if (wik1[ic]) nchan0++;
  }
  fprintf(f, "%6d\n", nchan0);
  for (ic = 0; ic < nchan; ic++) {
    if (wik1[ic]) {
      i = ic/rbasis.nkappa;
      k = ic%rbasis.nkappa;
      ilev = ts[i];
      lev = GetLevel(ilev);
      ka = KappaFromIndex(k);
      fprintf(f, "%6d %3d %15.8E %3d %2d\n",
	      ic, ilev, lev->energy, k, ka);
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

int RMatrixSurface(char *fn) {
  int i, m, k, t, kb, ilev, ika;
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
  if (rbasis.ib0 == 0) {
    for (i = 0; i < ncs; i++) {
      lev = GetLevel(cs[i]);
      DecodePJ(lev->pj, &p0, &j0);
      AddStateToSymmetry(-(cs[i]+1), -1, j0, p0, j0);
    }
  }

  nchan = nts*rbasis.nkappa;
  printf("%d %d %d %d\n", nchan, nts, rbasis.nkappa, ncs);
  wik0 = malloc(sizeof(double *)*nchan);
  wik1 = malloc(sizeof(double *)*nchan);
  for (t = 0; t < nchan; t++) {
    wik0[t] = NULL;
    wik1[t] = NULL;
  }
  for (i = 0; i < MAX_SYMMETRIES; i++) {
    if (rbasis.ib0 == 0) {
      k = ConstructHamiltonFrozen(i, ntg, tg, 0, ncg, cg);
    } else {
      k = ConstructHamiltonFrozen(i, ntg, tg, -1, 0, NULL);
    }
    if (k < 0) continue;
    DiagnolizeHamilton();
    h = GetHamilton();
    ek = h->mixing;
    mix = h->mixing+h->dim;
    sym = GetSymmetry(i);
    printf("sym: %d %d\n", i, h->dim);
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
    WriteRMatrixSurface(f, wik0, wik1);
  }
  
  free(wik0);
  free(wik1);
  fclose(f);
  return 0;
}
