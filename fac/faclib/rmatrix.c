#include "rmatrix.h"

static RBASIS rbasis;

int RMatrixBasis(int kmax, int nb) {
  int k, i, j, k2, ib, kappa, t, in;
  int nkb0, nkb1, n, n0, nmax;
  double e0, e1, a, b, c, rb, bb;
  ORBITAL *orb, orbf;
  
  ib = GetBoundary(&rb, &bb, &nmax);
  rbasis.kmax = kmax;
  rbasis.nbk = nb;
  rbasis.nkappa = (2*kmax + 1);
  rbasis.nbuttle = 2*(nb-1);
  rbasis.basis = malloc(sizeof(int *)*rbasis.nkappa);
  rbasis.ebuttle = malloc(sizeof(double *)*rbasis.nkappa);
  rbasis.cbuttle = malloc(sizeof(double *)*rbasis.nkappa);
  for (i = 0; i < rbasis.nkappa; i++) {
    rbasis.basis[i] = malloc(sizeof(int)*rbasis.nbk);
    rbasis.ebuttle[i] = malloc(sizeof(double)*rbasis.nbuttle);
    rbasis.cbuttle[i] = malloc(sizeof(double)*rbasis.nbuttle);
  }

  nkb0 = GetNumOrbitals();
  e1 = 0.0;
  e0 = 0.0;
  for (k = 0; k <= kmax; k++) {
    n0 = k+1;
    if (n0 <= nmax) n0 = nmax + 1;
    k2 = 2*k;
    t = k2 - 1;
    for (j = k2-1; j <= k2+1; j += 2, t++) {     
      if (j < 0) continue;
      kappa = GetKappaFromJL(j, k2);
      for (in = 0; in < nb; in++) {
	n = in + n0;
	rbasis.basis[t][in] = OrbitalIndex(n, kappa, 0.0);	
	orb = GetOrbital(rbasis.basis[t][in]);	
	e0 = e1;
	e1 = orb->energy;
	printf("%2d %2d %3d %3d %15.8E %15.8E %15.8E\n",
	       kappa, orb->n, orb->ilast, ib, orb->energy,
	       rb, WLarge(orb)[orb->ilast]);
	if (in > 0) {
	  rbasis.ebuttle[t][in*2-2] = e0 + (e1 - e0)*0.25;
	  rbasis.ebuttle[t][in*2-1] = e1 - (e1 - e0)*0.25;
	}
      }
      for (i = 0; i < rbasis.nbuttle; i++) {
	rbasis.cbuttle[t][i] = 0.0;
	orbf.n = -1;
	orbf.kappa = kappa;
	orbf.energy = rbasis.ebuttle[t][i];
	SolveDirac(&orbf);
	a = 0.0;
	for (in = 0; in < nb; in++) {
	  orb = GetOrbital(rbasis.basis[t][in]);
	  b = WLarge(orb)[orb->ilast];
	  b = b*b/(orb->energy - orbf.energy);
	  a += b;
	}
	a /= 2.0*rb;
	c = WSmall(&orbf)[ib]/WLarge(&orbf)[ib];
	c = c*2.0*rb/FINE_STRUCTURE_CONST - bb - kappa;
	rbasis.cbuttle[t][i] = 1/c - a;
	printf("%2d %15.8E %11.4E %11.4E %11.4E\n",
	       kappa, orbf.energy, rbasis.cbuttle[t][i], 1/c, a);
	free(orbf.wfun);
      }
    }
  }
  
  nkb1 = GetNumOrbitals();
  PrepSlater(0, nkb0-1, nkb0, nkb1-1, 0, nkb0-1, nkb0, nkb1-1);
  PrepSlater(0, nkb0-1, 0, nkb0-1, nkb0, nkb1-1, nkb0, nkb1-1);

  return 0;
}

int IndexKappa(int ka) {
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
    
int RMatrix(char *fn, int nt, int *kt, int nc, int *kc) {
  int nlevs, nts, ncs, i, m, k, t, kb, ilev, ika, *ts, *cs;
  int j0, p0, j1, k1, j, p, jmin, jmax, nchan, q, ic;
  double *ek, *mix, **wik;
  SYMMETRY *sym;
  STATE *s;
  LEVEL *lev;
  ORBITAL *orb;
  HAMILTON *h;
  FILE *f;

  nlevs = GetNumLevels();
  printf("%d %d %d\n", nlevs, nt, nc);
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
	    AddStateToSymmetry(-(i+1), kb, j, p, j);
	  }
	}
      }
    } else if (InGroups(s->kgroup, nc, kc)) {
      cs[ncs] = i;
      ncs++;
      DecodePJ(lev->pj, &p0, &j0);
      AddStateToSymmetry(-(i+1), -1, j0, p0, j0);
    }
  }
  if (ncs == 0) {
    free(cs);
  }

  f = fopen(fn, "w");

  AngularFrozen(nts, ts, ncs, cs);
  nchan = nts*rbasis.nkappa;
  printf("%d %d %d %d\n", nchan, nts, rbasis.nkappa, ncs);
  wik = malloc(sizeof(double *)*nchan);
  for (t = 0; t < nchan; t++) {
    wik[t] = NULL;
  }
  for (i = 0; i < MAX_SYMMETRIES; i++) {
    k = ConstructHamiltonFrozen(i, nt, kt, 0, nc, kc);
    if (k < 0) continue;
    DiagnolizeHamilton();
    h = GetHamilton();
    ek = h->mixing;
    mix = h->mixing+h->dim;
    sym = GetSymmetry(i);
    printf("%d %d\n", i, h->dim);
    fprintf(f, "ISYM: %d\n", i);
    for (t = 0; t < h->dim; t++) {
      for (q = 0; q < h->dim; q++) {
	k = h->basis[q];
	s = ArrayGet(&(sym->states), k);
	k = -(s->kgroup+1);
	ilev = IBisect(k, nts, ts);
	if (ilev >= 0) {
	  kb = s->kcfg;
	  orb = GetOrbital(kb);
	  ika = IndexKappa(orb->kappa);
	  ic = ilev*rbasis.nkappa + ika;
	  if (wik[ic] == NULL) {
	    wik[ic] = malloc(sizeof(double)*h->dim);
	    for (p = 0; p < h->dim; p++) {
	      wik[ic][p] = 0.0;
	    }
	  }
	  wik[ic][t] += mix[q]*WLarge(orb)[orb->ilast];
	}
      }
      mix += h->dim;
    }
    if (i == 29) {
      AddToLevels(0, NULL);
    }
    for (t = 0; t < h->dim; t++) {
      fprintf(f, "%15.8E ", ek[t]);
    }
    fprintf(f, "\n");
    for (ic = 0; ic < nchan; ic++) {
      if (wik[ic]) {
	fprintf(f, "%2d %2d %2d\n", ic,
		ic/rbasis.nkappa, KappaFromIndex(ic%rbasis.nkappa));
	for (t = 0; t < h->dim; t++) {
	  fprintf(f, "%15.8E ", wik[ic][t]);
	}
	fprintf(f, "\n");
	free(wik[ic]);
	wik[ic] = NULL;
      }
    }
  }
  
  SaveLevels("tb.en", nlevs, -1);
  free(wik);
  ClearAngularFrozen();
  fclose(f);

  return 0;
}
