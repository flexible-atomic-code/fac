#include "structure.h"
#include <time.h>

#if (FAC_DEBUG >= DEBUG_STRUCTURE)
#define debug_integral(s, ne, r) \
        {int ks; \
         for (ks = 0; ks < (2*(ne)); ks++) { \
         fprintf(debug_log, "%d %d %d %d %d\n", \
	         (s)[ks].index, (s)[ks].n, (s)[ks].kappa,\
                 (s)[ks].nq_bra, (s)[ks].nq_ket); \
        } \
         fprintf(debug_log, "%d electron: %lf\n\n", (ne), (r)); \
        }
#endif /* (FAC_DEBUG >= DEBUG_STRUCTURE) */

static ARRAY *hamiltons;
static int n_hamiltons = 0;

static ARRAY *levels;
static int n_levels = 0;

static MULTI *angz_array;
static MULTI *angzxz_array;

static int rydberg_ignored = 0;
static double angz_cut = EPS4;
static double mix_cut = EPS3;

static STRUCT_TIMING timing = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

int GetStructTiming(STRUCT_TIMING *t) {
  memcpy(t, &timing, sizeof(timing));
  return 0;
}

int SetAngZOptions(int n, double mix_cut, double cut) {
  rydberg_ignored = n;
  mix_cut = mix_cut;
  angz_cut = cut;
}

HAMILTON *GetNewHamilton() {
  int i;
  HAMILTON *h;

  h = (HAMILTON *) ArrayAppend(hamiltons, NULL);
  if (!h) {
    printf("Not sufficient memory for hamiltons array\n");
    abort();
  }
 
  h->basis = NULL;
  h->mixing = NULL;
  n_hamiltons++;
  return h;
}

HAMILTON *GetHamilton(int ih) {
  return (HAMILTON *) ArrayGet(hamiltons, ih);
}

int ConstructHamilton(int isym, int k, int *kg) {
  int i, j, t;
  CONFIG *cfg_i, *cfg_j;
  HAMILTON *h;
  ARRAY *st;
  STATE *s;
  SYMMETRY *sym;
  double r;
  clock_t start, stop;

#if (FAC_DEBUG >= DEBUG_STRUCTURE) 
  char name[LEVEL_NAME_LEN];
#endif

#ifdef PERFORM_STATISTICS
  start = clock();
#endif

  j = 0;
  sym = GetSymmetry(isym);
  if (sym == NULL) return -1;
  st = &(sym->states);
  for (t = 0; t < sym->n_states; t++) {
    s = (STATE *) ArrayGet(st, t);
    if (InGroups(s->kgroup, k, kg)) j++;
  }
  
  if (j == 0) return -1;

  h = GetNewHamilton();
  if (!h) return -1;
  h->pj = isym;

  h->dim = j;
  t = j*(j+1)/2;
  /*  if (h->basis != NULL) free(h->basis); */
  h->basis = (int *) malloc(sizeof(int)*j);
  if (!(h->basis)) goto ERROR;
  h->hamilton = (double *) malloc(sizeof(double)*t);
  if (!(h->hamilton)) goto ERROR;

  j = 0;
  for (t = 0; t < sym->n_states; t++) {
    s = (STATE *) ArrayGet(st, t);
    if (InGroups(s->kgroup, k, kg)) {
      h->basis[j] = t;
      j++;
    }
  }
#if (FAC_DEBUG >= DEBUG_STRUCTURE)
  fprintf(debug_log, "%d %d %X \n", h->dim, n_hamiltons-1, h->pj);
#endif /* (FAC_DEBUG >= DEBUG_STRUCTURE) */

  for (j = 0; j < h->dim; j++) {
    t = j*(j+1)/2;
    for (i = 0; i <= j; i++) {
      r = HamiltonElement(isym, h->basis[i], h->basis[j]);
      h->hamilton[i+t] = r;
    }
  }

#ifdef PERFORM_STATISTICS
  stop = clock();
  timing.set_ham += stop-start;
#endif

  return n_hamiltons-1;

 ERROR:
  if (h->basis) free(h->basis);
  if (h->hamilton) free(h->hamilton);
#ifdef PERFORM_STATISTICS
  stop = clock();
  timing.set_ham += stop-start;
#endif
  return -1;
}

int ValidBasis(STATE *s, int k, int *kg, int n) {
  int t, i, m, kb;
  LEVEL *lev;
  HAMILTON *h;
  STATE *sp;
  SYMMETRY *sym;
  
  t = s->kgroup;
  if (t >= 0) return 0;
  
  kb = s->kcfg;
  kb = GetOrbital(kb)->n;
  if (kb != n) return 0;

  t = -t-1;
  if (kg) {
    lev = GetLevel(t);
    i = lev->ham_index;
    h = GetHamilton(i);
    m = lev->major_component;
    sym = GetSymmetry(h->pj);
    sp = (STATE *) ArrayGet(&(sym->states), h->basis[m]);
    t = sp->kgroup;
    return InGroups(t, k, kg);
  } else {
    if (t == k) return 1;
    else return 0;
  }
}


int ConstructHamiltonFrozen(int isym, int k, int *kg, int n) {
  int i, j, t;
  HAMILTON *h;
  ARRAY *st;
  STATE *s;
  SYMMETRY *sym;
  double r;
  clock_t start, stop;

#ifdef PERFORM_STATISTICS
  start = clock();
#endif

  j = 0;
  sym = GetSymmetry(isym);
  st = &(sym->states);
  for (t = 0; t < sym->n_states; t++) { 
    s = (STATE *) ArrayGet(st, t);
    if (i = ValidBasis(s, k, kg, n)) {
      j++;
    } 
  }
  if (j == 0) return -1;
  
  h = GetNewHamilton();
  if (!h) return -1;
  h->pj = isym;
  
  h->dim = j;
  t = j*(j+1)/2;
  h->basis = (int *) malloc(sizeof(int)*j); 
  if (!(h->basis)) goto ERROR;
  h->hamilton = (double *) malloc(sizeof(double)*t); 
  if (!(h->hamilton)) goto ERROR;
      
  j = 0;
  for (t = 0; t < sym->n_states; t++) { 
    s = (STATE *) ArrayGet(st, t);
    if (ValidBasis(s, k, kg, n)) {
      h->basis[j] = t;
      j++;
    }
  }

  /*  printf("%d %d\n", h->dim, h->pj); */
  for (j = 0; j < h->dim; j++) {
    t = j*(j+1)/2;
    for (i = 0; i <= j; i++) {
      r = HamiltonElementFrozen(isym, h->basis[i], h->basis[j]);
      h->hamilton[i+t] = r;
    }
  }

#ifdef PERFORM_STATISTICS
  stop = clock();
  timing.set_ham += stop-start;
#endif

  return n_hamiltons-1;

 ERROR:
  if (h->basis) free(h->basis);
  if (h->hamilton) free(h->hamilton);

#ifdef PERFORM_STATISTICS
  stop = clock();
  timing.set_ham += stop-start;
#endif
  return -1;
}

double HamiltonElementFrozen(int isym, int isi, int isj) {
  STATE *si, *sj;
  double r, r0, sd, se, a;
  int i, ti, tj, ji1, ji2, jj1, jj2, ki2, kj2, j, nz;
  int ks[4];
  ORBITAL *orbi, *orbj;
  ANGULAR_ZMIX *ang;
  SYMMETRY *sym;
  HAMILTON *h1, *h2;
  LEVEL *lev1, *lev2;

  sym = GetSymmetry(isym);
  si = (STATE *) ArrayGet(&(sym->states), isi);
  sj = (STATE *) ArrayGet(&(sym->states), isj);
  r = 0.0;
  ti = si->kgroup;
  tj = sj->kgroup;
  ti = -ti-1;
  tj = -tj-1;
  orbi = GetOrbital(si->kcfg);
  orbj = GetOrbital(sj->kcfg);
  lev1 = GetLevel(ti);
  lev2 = GetLevel(tj);
  h1 = GetHamilton(lev1->ham_index);
  h2 = GetHamilton(lev2->ham_index);
  ji1 = h1->pj;
  jj1 = h2->pj;
  DecodePJ(ji1, NULL, &ji1);
  DecodePJ(jj1, NULL, &jj1);
  GetJLFromKappa(orbi->kappa, &ji2, &ki2);
  GetJLFromKappa(orbj->kappa, &jj2, &kj2);
  j = si->kstate;
  if (si->kgroup == sj->kgroup) { 
    if (ji2 == jj2 && ki2 == kj2) {
      ResidualPotential(&r0, si->kcfg, sj->kcfg);
      r += r0;
    } 
  
    if (si->kcfg == sj->kcfg) {
      r += lev1->energy;
      r += orbi->energy;
    }
  }
 
  ks[1] = si->kcfg;
  ks[3] = sj->kcfg;

  nz = AngularZMix(&ang, ti, tj, -1, -1);
  a = 0.0;
  for (i = 0; i < nz; i++) {
    if (fabs(ang[i].coeff) < EPS10) continue;
    r0 = W6j(ji1, ji2, j, jj2, jj1, ang[i].k);
    if (fabs(r0) < EPS10) continue;
    ks[0] = ang[i].k0;
    ks[2] = ang[i].k1;
    SlaterTotal(&sd, &se, NULL, ks, ang[i].k, 0);
    r0 *= ang[i].coeff*(sd+se);
    a += r0;
  }
  
  /*  a /= sqrt(j+1.0);*/
  if (IsOdd((ji2 + jj1 + j)/2)) a = -a;
  r += a;

  if (nz > 0) {
    free(ang);
  } 

  /*
  if (ti == tj && si->kcfg == sj->kcfg) {
    printf("%10.3E %10.3E\n", a, r);
  } 
  */
  return r;
} 
    

double HamiltonElement(int isym, int isi, int isj) { 
  CONFIG *ci, *cj;
  int ki, kj;
  SYMMETRY *sym;
  STATE *si, *sj;
  SHELL_STATE *sbra, *sket;
  SHELL *bra;
  int n_shells, i, j;
  INTERACT_SHELL s[4];
  double x, r;
  int phase;

  sym = GetSymmetry(isym);
  si = (STATE *) ArrayGet(&(sym->states), isi);
  sj = (STATE *) ArrayGet(&(sym->states), isj);
  ci = GetConfig(si);
  if (ci->n_shells == 0) return 0.0;
  cj = GetConfig(sj);
  if (cj->n_shells == 0) return 0.0;
  
  ki = si->kstate;
  kj = sj->kstate;
  for (i = 0; i < 4; i++) s[i].index = -1;
  n_shells = GetInteract(&phase, s, &bra, &sbra, &sket, 
			 ci, ki, cj, kj, si, sj);
  if (n_shells <= 0) return 0.0;
  
  x = 0.0;
  if (s[0].index >= 0 && s[3].index >= 0) {
    r = Hamilton2E(n_shells, sbra, sket, s);
    x += r;

#if (FAC_DEBUG >= DEBUG_STRUCTURE)
    debug_integral(s, 2, r);
#endif

  } else if( s[0].index >= 0) {
    r = Hamilton1E(n_shells, sbra, sket, s);
    x += r;

#if (FAC_DEBUG >= DEBUG_STRUCTURE)
    debug_integral(s, 1, r);
#endif

    for (i = 0; i < n_shells; i++) {
      s[2].index = n_shells - i - 1;
      s[3].index = s[2].index;
      s[2].n = bra[i].n;
      s[3].n = s[2].n;
      s[2].kappa = bra[i].kappa;
      s[3].kappa = s[2].kappa;
      s[2].j = GetJ(bra+i);
      s[3].j = s[2].j;
      s[2].kl = GetL(bra+i);
      s[3].kl = s[2].kl;
      s[2].nq_bra = GetNq(bra+i);
      if (s[2].index == s[0].index) {
	s[2].nq_ket = s[2].nq_bra - 1;
      } else if (s[2].index == s[1].index) {
	s[2].nq_ket = s[2].nq_bra + 1;
      } else {
	s[2].nq_ket = s[2].nq_bra;
      }
      if (s[2].nq_bra < 0 || s[2].nq_ket < 0 ||
	  s[2].nq_bra > s[2].j+1 || s[2].nq_ket > s[2].j+1.0) {
	continue;
      }
      s[3].nq_bra = s[2].nq_bra;
      s[3].nq_ket = s[2].nq_ket;
      r = Hamilton2E(n_shells, sbra, sket, s);
      x += r;
#if (FAC_DEBUG >= DEBUG_STRUCTURE)
      debug_integral(s, 2, r);
#endif
    }
  } else {
    for (i = 0; i < n_shells; i++) {
      s[0].index = n_shells - i - 1;
      s[1].index = s[0].index;
      s[0].n = bra[i].n;
      s[1].n = s[0].n;
      s[0].kappa = bra[i].kappa;
      s[1].kappa = s[0].kappa;
      s[0].j = GetJ(bra+i);
      s[1].j = s[0].j;
      s[0].kl = GetL(bra+i);
      s[1].kl = s[0].kl;
      s[0].nq_bra = GetNq(bra+i);
      s[0].nq_ket = s[0].nq_bra;
      s[1].nq_bra = s[0].nq_bra;
      s[1].nq_ket = s[1].nq_bra;

      r = Hamilton1E(n_shells, sbra, sket, s);
      x += r;
#if (FAC_DEBUG >= DEBUG_STRUCTURE)
      debug_integral(s, 1, r);
#endif
      for (j = 0; j <= i; j++) {
	s[2].index = n_shells - j - 1;
	s[3].index = s[2].index;
	s[2].n = bra[j].n;
	s[3].n = s[2].n;
	s[2].kappa = bra[j].kappa;
	s[3].kappa = s[2].kappa;
	s[2].j = GetJ(bra+j);
	s[3].j = s[2].j;
	s[2].kl = GetL(bra+j);
	s[3].kl = s[2].kl;
	s[2].nq_bra = GetNq(bra+j);
	s[2].nq_ket = s[2].nq_bra;
	s[3].nq_bra = s[2].nq_bra;
	s[3].nq_ket = s[3].nq_bra;
	r = Hamilton2E(n_shells, sbra, sket, s);
	x += r;

#if (FAC_DEBUG >= DEBUG_STRUCTURE)
	debug_integral(s, 2, r);
#endif

      }
    }
  }
  /* the prefactor in the Wigner-Eckart theorem should be included in 
     the matrix element. for a scalar operator, this is [J]^{-1/2} */
  x /= sqrt(sbra[0].totalJ + 1.0);
  if (IsOdd(phase)) x = -x;

  free(sbra);
  free(sket);

  return x;
}

double Hamilton1E(int n_shells, SHELL_STATE *sbra, SHELL_STATE *sket,
		  INTERACT_SHELL *s) {
  int nk0, k;
  int *k0;
  double *x, z0, r0;
  int k1, k2;

  if (s[0].j != s[1].j ||
      s[0].kl != s[1].kl) return 0.0;
  nk0 = 1;
  k = 0;
  k0 = &k;
  x = &z0;
  nk0 = AngularZ(&x, &k0, nk0, n_shells, sbra, sket, s, s+1);
  k1 = OrbitalIndex(s[0].n, s[0].kappa, 0.0);
  k2 = OrbitalIndex(s[1].n, s[1].kappa, 0.0);
  if (fabs(z0) < EPS10) return 0.0;
  ResidualPotential(&r0, k1, k2);

#if (FAC_DEBUG >= DEBUG_STRUCTURE)
  fprintf(debug_log, "Radial1E: %lf\n", r0);
#endif

  if (k1 == k2) r0 += (GetOrbital(k1))->energy;

#if (FAC_DEBUG >= DEBUG_STRUCTURE)
  fprintf(debug_log, "Angular: %lf, Radial: %lf\n", z0, r0);
#endif

  return r0 * z0 * sqrt(s[0].j + 1.0);
}

double Hamilton2E(int n_shells, SHELL_STATE *sbra, SHELL_STATE *sket, 
		  INTERACT_SHELL *s) {
  int nk0, nk, *kk, k, *kk0, i;
  double *ang;
  double se, sd, x;
  double z0, *y;
  int ks[4], js[4];

  js[0] = 0;
  js[1] = 0;
  js[2] = 0;
  js[3] = 0;
  
  ks[0] = OrbitalIndex(s[0].n, s[0].kappa, 0.0);
  ks[1] = OrbitalIndex(s[2].n, s[2].kappa, 0.0);
  ks[2] = OrbitalIndex(s[1].n, s[1].kappa, 0.0);
  ks[3] = OrbitalIndex(s[3].n, s[3].kappa, 0.0);
 
  z0 = 0.0;
  nk0 = 0;
  if (ks[1] == ks[2]) {
    nk0 = 1;
    k = 0;
    kk0 = &k;
    y = &z0;
    nk0 = AngularZ(&y, &kk0, nk0, n_shells, sbra, sket, s, s+3);
    if (nk0 > 0) {
      z0 /= sqrt(s[0].j + 1.0);
      if (IsOdd((s[0].j - s[1].j)/2)) z0 = -z0;
    }
  }

  x = 0.0;
    
  nk = AngularZxZ0(&ang, &kk, 0, n_shells, sbra, sket, s);
  for (i = 0; i < nk; i++) {

#if (FAC_DEBUG >= DEBUG_STRUCTURE)
    fprintf(debug_log, "rank: %d, Angular: %lf \n", kk[i], ang[i]);
#endif
    if (fabs(ang[i]) > EPS10) {
      SlaterTotal(&sd, &se, js, ks, kk[i], 0);
      x += ang[i] * (sd+se);
    } else if (nk0 > 0) {
      SlaterTotal(&sd, NULL, js, ks, kk[i], 0);
    }
    if (nk0 > 0) x -= z0 * sd;

  }
    
  if (nk > 0) {
    free(ang);
    free(kk);
  }

  return x;
}
    
int DiagnolizeHamilton(int ih) {
  double *ap;
  double *w;
  double *z;
  char jobz[] = "V";
  char uplo[] = "U";
  int n;
  int ldz;
  double *work;
  int lwork;
  int *iwork;
  int liwork;
  int info;
  HAMILTON *h;
  clock_t start, stop;

#ifdef PERFOM_SATISTICS
  start = clock();
#endif

  h = GetHamilton(ih);
  n = h->dim;
  ap = h->hamilton;
  ldz = n;
  lwork = 1 + 6*n + n*n;
  liwork = 3 + 5*n;
  work = malloc(sizeof(double)*lwork);
  iwork = malloc(sizeof(int)*liwork);

  h->mixing = malloc(sizeof(double)*n*(n+1));
  if (!work ||
      !iwork ||
      !(h->mixing)) goto ERROR;
  w = h->mixing;
  z = h->mixing + n;

  dspevd_(jobz, uplo, &n, ap, w, z, &ldz, work, &lwork,
	  iwork, &liwork, &info);
  if (info) goto ERROR;
  free(work);
  free(iwork);
  return 0;
  
 ERROR:
  if(work) free(work);
  if(iwork) free(iwork);
  if(h->mixing) free(h->mixing);

#ifdef PERFORM_STATISTICS
  stop = clock();
  timing.diag_ham += stop-start;
#endif

  return -1;
}

int AddToLevels(int ih) {
  int i, d, j, k;
  HAMILTON *h;
  LEVEL lev;
  double *mix;
  
  h = GetHamilton(ih);
  if (h->basis == NULL ||
      h->mixing == NULL) return -1;
  d = h->dim;
  mix = h->mixing + d;
  j = n_levels;
  for (i = 0; i < d; i++) {
    lev.energy = h->mixing[i];
    lev.ham_index = ih;
    lev.mix_index = i;
    k = GetPrincipleBasis(mix, d);
    lev.major_component = k;
    if (ArrayAppend(levels, &lev) == NULL) {
      printf("Not enough memory for levels array\n");
      abort();
    }
    j++;
    mix += d;
  }
  
  n_levels = j;
  if (i < d-1) return -2;
  return 0;
}

int CorrectEnergy(int n, int *k, double *e) {
  int i;
  LEVEL *lev;
  double e0;

  lev = GetLevel(0);
  e0 = lev->energy;
  for (i = 0; i < n; i++) {
    lev = GetLevel(k[i]);
    if (lev == NULL) {
      printf("Level %d does not exist when correcting energy\n", k[i]);
      continue;
    }
    lev->energy = e[i] + e0;
  }
  return 0;
}

LEVEL *GetLevel(int k) {
  return (LEVEL *) ArrayGet(levels, k);
}

int LevelTotalJ(int k) {
  int i;
  i = GetLevel(k)->ham_index;
  i = GetHamilton(i)->pj;
  DecodePJ(i, NULL, &i);
  return i;
}

int GetNumLevels() {
  return n_levels;
}

int GetPrincipleBasis(double *mix, int d) {
  int i, k;
  double c = 0.0;
  double fm;
  for (i = 0; i < d; i++) {
    fm = fabs(mix[i]);
    if (fm > c) {
      c = fm;
      k = i;
    }
  }

  return k;
}

int SortLevels(int start, int n) {
  int i, j, i0, j0;
  LEVEL tmp, *lev1, *lev2, *levp;
  if (n < 0) n = n_levels-start;
  while (1 < n) {
    i = start;
    j = start + n - 1;
    lev1 = GetLevel(i);
    lev2 = GetLevel(j);
    levp = lev2;
    
    while (i < j) {
      while (i < j) {
	if (lev1->energy > levp->energy) break;
	i++;
	lev1 = GetLevel(i);
      }
      while (i < j) {
	if (levp->energy > lev2->energy) break;
	j--;
	lev2 = GetLevel(j);
      }
      if (i < j) {
	memcpy(&tmp, lev1, sizeof(LEVEL));
	memcpy(lev1, lev2, sizeof(LEVEL));
	memcpy(lev2, &tmp, sizeof(LEVEL));
	i++;
	lev1 = GetLevel(i);
      }
    }
    if (lev1 != levp) {
      memcpy(&tmp, lev1, sizeof(LEVEL));
      memcpy(lev1, levp, sizeof(LEVEL));
      memcpy(levp, &tmp, sizeof(LEVEL));
    }

    i0 = i - start;
    j0 = n - i0 - 1;
    if (j0 < i0) {
      if (1 < j0) {
	SortLevels(i+1, j0);
      }
      n = i0;
    } else {
      if (1 < i0) {
	SortLevels(start, i0);
      }
      start = i+1;
      n = j0;
    }
  }
}

int SaveLevelsToAscii(char *fn, int m, int n) {
  FILE *f;
  char name[LEVEL_NAME_LEN];
  char sname[LEVEL_NAME_LEN];
  int i, j, k;
  int hi, si;
  double mix;
  STATE *s, *sp;
  ORBITAL *orb;
  CONFIG *c;
  SYMMETRY *sym;
  HAMILTON *h;
  LEVEL *lev;
  char *asym;
  double e0;

  f = fopen(fn, "w");
  if (f == NULL) return -1;
  
  e0 = GetLevel(0)->energy;
  asym = GetAtomSymbol();
  fprintf(f, "%s,  Ground State Energy: %10.4E\n\n", 
	  asym, e0*HARTREE_EV);
 
  for (i = 0; i < n_levels; i++) {
    lev = GetLevel(i);
    hi = lev->ham_index;
    si = lev->major_component;
    h = GetHamilton(hi);
    sym = GetSymmetry(h->pj);
    s = (STATE *) ArrayGet(&(sym->states), h->basis[si]);
    if (n > 0) {
      if (s->kgroup >= 0) {
	c = GetConfig(s);
	if (n != c->shells[0].n) continue;
      } else {
	orb = GetOrbital(s->kcfg);
	if (n != orb->n) continue;
      }
    }
    ConstructLevelName(name, sname, s);
    fprintf(f, "%-5d %11.4E  %-20s  %-30s\n", i, 
	    (lev->energy-e0)*HARTREE_EV, sname, name);
    if (m == 0) continue;

    si = lev->mix_index;
    k = h->dim * (si + 1);
    for (j = 0; j < h->dim; j++) {
      mix = h->mixing[k];
      sp = (STATE *) ArrayGet(&(sym->states), h->basis[j]);
      fprintf(f, "  (%2d %2d %2d) %10.4E \n", 
	      sp->kgroup,
	      sp->kcfg,
	      sp->kstate, mix);
      k++;
    }
  }
  fclose(f);
  
  return 0;
}

int ConstructLevelName(char *name, char *sname, STATE *basis) {
  int n, nq, kl, j;
  int i, len;
  char symbol[20];
  char jsym;
  char ashell[16];
  CONFIG *c;
  SHELL_STATE *s;
  ORBITAL *orb;
  LEVEL *lev;
  HAMILTON *h;
  SYMMETRY *sym;
  int hi, si;
  int n0, kl0, nq0;

  symbol[0] = '\0';
  if (basis->kgroup < 0) {
    i = basis->kgroup;
    i = -(i + 1);
    if (name) {
      orb = GetOrbital(basis->kcfg);
      GetJLFromKappa(orb->kappa, &j, &kl);
      if (j < kl) jsym = '-';
      else jsym = '+';
      
      kl /= 2;
      SpecSymbol(symbol, kl);
      sprintf(name, "%5d + %d%s%c1(%d)%d ", 
	      i, orb->n, symbol, jsym, j, basis->kstate);
    }
    if (sname) {
      lev = GetLevel(i);
      hi = lev->ham_index;
      si = lev->major_component;
      h = GetHamilton(hi);
      sym = GetSymmetry(h->pj);
      basis = (STATE *) ArrayGet(&(sym->states), h->basis[si]);
      ConstructLevelName(NULL, sname, basis);
    }
    return 0;
  }

  c = GetConfig(basis);
  s = c->csfs + basis->kstate;

  len = 0;
  name[0] = '\0';
  sname[0] = '\0';
  n0 = 0;
  kl0= -1;
  nq0 = 0;
  for (i = c->n_shells-1; i >= 0; i--) {
    UnpackShell(c->shells+i, &n, &kl, &j, &nq);
    if (name) {
      if (j < kl) jsym = '-';
      else jsym = '+';
      kl = kl/2;
      SpecSymbol(symbol, kl);
      sprintf(ashell, "%1d%s%c%1d(%1d)%1d ", 
	      n, symbol, jsym, nq, s[i].shellJ, s[i].totalJ); 
      len += strlen(ashell);
      if (len >= LEVEL_NAME_LEN) return -1;
      strcat(name, ashell);
    }
    if (sname) {
      if (n == n0 && kl == kl0) {
	nq0 += nq;
      } else {
	if (n0 > 0) {
	  SpecSymbol(symbol, kl0);
	  sprintf(ashell, "%1d%s%1d ", n0, symbol, nq0);
	  strcat(sname, ashell);
	}
	n0 = n;
	kl0 = kl;
	nq0 = nq;
      }
    }
  }
  
  if (n0 > 0) {
    SpecSymbol(symbol, kl0);
    sprintf(ashell, "%1d%s%1d ", n0, symbol, nq0);
    strcat(sname, ashell);
  }
  return len;
}
    

int GetBasisTable(char *fn) {
  FILE *f;
  int i, p, j, k, nsym;
  char name[LEVEL_NAME_LEN];
  char sname[LEVEL_NAME_LEN];
  ARRAY *st;
  STATE *s;
  SYMMETRY *sym;

  f = fopen(fn, "w");
  if (!f) return -1;

  nsym = MAX_SYMMETRIES;
  for (i = 0; i < nsym; i++) {
    sym = GetSymmetry(i);
    DecodePJ(i, &p, &j);
    st = &(sym->states);
    fprintf(f, "%d: J = %d, Parity = %d\n-------------------\n", i, j, p);
    for (k = 0; k < sym->n_states; k++) {
      s = (STATE *) ArrayGet(st, k);
      ConstructLevelName(name, sname, s);
      fprintf(f, "%4d (%2d %2d %2d) %-20s %-50s \n",
	      k, s->kgroup, s->kcfg, s->kstate, sname, name);
    }
    fprintf(f, "\n");
  }
  
  fclose(f);
}

int AngularZMixStates(ANGULAR_ZMIX **ang, STATE *s1, STATE *s2) {
  CONFIG *c1, *c2;
  int nz, n, p, q;
  int n_shells, *k, nkk;
  INTERACT_SHELL s[4];
  double *r;
  int phase, im;
  int orb0, orb1;
  SHELL *bra;
  SHELL_STATE *sbra, *sket;
  STATE *stmp;
  int j1, j2;
  ANGZ_DATUM *angz_datum;
  int index[6];
  clock_t start, stop;

#ifdef PERFORM_STATISTICS
  start = clock();
#endif

  index[0] = s2->kgroup;
  index[1] = s1->kgroup;
  index[2] = s2->kcfg;
  index[3] = s1->kcfg;
  index[4] = s2->kstate;
  index[5] = s1->kstate;
  angz_datum = (ANGZ_DATUM *) MultiSet(angz_array, index, NULL);
  nz = angz_datum->nz;
  if (nz < 0) {
#ifdef PERFORM_STATISTICS
    stop = clock();
    timing.angz_states += stop-start;
#endif
    return -1;
  }
  
  if (nz > 0) {
    (*ang) = (ANGULAR_ZMIX *) (angz_datum->angz);
#ifdef PERFORM_STATISTICS
    stop = clock();
    timing.angz_states += stop-start;
#endif
    return nz;
  }

  c1 = GetConfig(s1);
  c2 = GetConfig(s2);
  n_shells = GetInteract(&phase, s, &bra, &sbra, &sket, 
			 c1, s1->kstate, c2, s2->kstate, s1, s2);
			
  if (n_shells <= 0) return -1;
  if (s[2].index >= 0 && s[0].index >= 0) {
    free(sbra);
    free(sket);
    return -1;
  }

  nz = ANGZ_BLOCK;
  n = 0;
  (*ang) = malloc(sizeof(ANGULAR_ZMIX)*nz);
  if (!(*ang)) return -1;
  
  if (s[0].index >= 0) {
    nkk = AngularZ(&r, &k, 0, n_shells, sbra, sket, s, s+1);
    if (nkk > 0) {
      orb0 = OrbitalIndex(s[0].n, s[0].kappa, 0.0);
      orb1 = OrbitalIndex(s[1].n, s[1].kappa, 0.0);
      for (p = 0; p < nkk; p++) {
	if (fabs(r[p]) > angz_cut) {
	  if (IsOdd(phase)) r[p] = -r[p];
	  im = AddToAngularZMix(&n, &nz, ang, k[p], orb0, orb1, r[p]);
	  if (im < 0) {
	    free(r);
	    free(k);
	    free(*ang);
	    return -1;
	  }
	}
      }
      free(r);
      free(k);
    }
  } else {
    for (q = 0; q < n_shells; q++) {	    
      s[0].index = n_shells - q - 1;      
      s[1].index = s[0].index;
      s[0].n = bra[q].n;
      s[1].n = s[0].n;
      s[0].kappa = bra[q].kappa;
      s[1].kappa = s[0].kappa;
      s[0].j = GetJ(bra+q);
      s[1].j = s[0].j;
      s[0].kl = GetL(bra+q);
      s[1].kl = s[0].kl;
      s[0].nq_bra = GetNq(bra+q);
      s[0].nq_ket = s[0].nq_bra;
      s[1].nq_bra = s[0].nq_bra;
      s[1].nq_ket = s[1].nq_bra;
      nkk = AngularZ(&r, &k, 0, n_shells, sbra, sket, s, s+1);
      if (nkk > 0) {
	orb0 = OrbitalIndex(s[0].n, s[0].kappa, 0.0);
	orb1 = orb0;
	for (p = 0; p < nkk; p++) {
	  if (fabs(r[p]) > angz_cut) {
	    if (IsOdd(phase)) r[p] = -r[p];
	    im = AddToAngularZMix(&n, &nz, ang, k[p], orb0, orb1, r[p]);
	    if (im < 0) {
	      free(r);
	      free(k);
	      free(*ang);
	      return -1;
	    }
	  }
	}	    
	free(r);
	free(k);
      }
    }
  }
      
  free(sbra);
  free(sket);  
  
  if (n <= 0) {
    angz_datum->nz = -1;
    free(*ang);
    return -1;
  } else {
    if (n < nz) (*ang) = realloc(*ang, sizeof(ANGULAR_ZMIX)*n);
    angz_datum->angz = (void *) (*ang);
    angz_datum->nz = n;
  }

#ifdef PERFORM_STATISTICS
  stop =clock();
  timing.angz_states += stop - start;
#endif
  return n;
}

int AngZSwapBraKet(int nz, ANGULAR_ZMIX *ang, int p) {
  int i;
  int k0, k1;
  for (i = 0; i < nz; i++) {
    k0 = ang[i].k0;
    k1 = ang[i].k1;
    k0 = GetOrbital(k0)->kappa;
    k1 = GetOrbital(k1)->kappa;
    k0 = GetJFromKappa(k0);
    k1 = GetJFromKappa(k1);
    if (IsOdd((k1-k0)/2+p)) ang[i].coeff = - ang[i].coeff;
  }
  return 0;
}
    
int AngularZFreeBoundStates(ANGULAR_ZFB **ang, STATE *slow, STATE *sup) {  
  int nz, j1, j2, ia, kb;
  int n_shells, *k, k0, nkk;
  INTERACT_SHELL s[4];
  SHELL *bra;
  SHELL_STATE *sbra, *sket;
  CONFIG *clow, *cup, clowf;
  int phase;
  int orb0, orb1;
  double *r, r0;
  SHELL sh[200];
  SHELL_STATE cs[200];  
  ANGULAR_ZFB afb;
  ANGZ_DATUM *angz_datum;
  int index[6];
  clock_t start, stop;
  
#ifdef PERFORM_STATISTICS
  start = clock();
#endif
  
  index[1] = slow->kgroup;
  index[0] = sup->kgroup;
  index[3] = slow->kcfg;
  index[2] = sup->kcfg;
  index[5] = slow->kstate;
  index[4] = sup->kstate;
  angz_datum = (ANGZ_DATUM *) MultiSet(angz_array, index, NULL);
  nz = angz_datum->nz;
  if (nz < 0) {
#ifdef PERFORM_STATISTICS
    stop = clock();
    timing.angzfb_states += stop-start;
#endif
    return -1;
  }
  if (nz == 1) {
    (*ang) = (ANGULAR_ZFB *) angz_datum->angz;
#ifdef PERFORM_STATISTICS
    stop = clock();
    timing.angzfb_states += stop-start;
#endif
    return nz;
  }

  k = &k0;
  r = &r0;
  k0 = 0;

  clowf.shells = sh;
  clowf.csfs = cs;
    
  clow = GetConfig(slow);
  j1 = clow->csfs[slow->kstate].totalJ;

  clowf.n_shells = clow->n_shells+1;
  memcpy(clowf.shells+1, clow->shells, sizeof(SHELL)*clow->n_shells);
  clowf.shells[0].n = 1000;
  clowf.shells[0].nq = 1;
  clowf.shells[0].kappa = -1;
  clowf.n_csfs = 1;
  memcpy(clowf.csfs+1, clow->csfs+slow->kstate, 
	 sizeof(SHELL_STATE)*clow->n_shells);	
  cup = GetConfig(sup);
  j2 = cup->csfs[sup->kstate].totalJ;

  clowf.csfs[0].shellJ = 1;
  clowf.csfs[0].totalJ = j2;
  clowf.csfs[0].nu = 1;
  clowf.csfs[0].Nr = 0; 

  n_shells = GetInteract(&phase, s, &bra, &sbra, &sket,
			 &clowf, 0, cup, sup->kstate, 
			 slow, sup);
  if (n_shells <= 0) {
    nz = -1;
    goto END;
  }
  if (s[0].index < 0 || (s[0].index >= 0 && s[2].index >= 0)) {
    free(sbra);
    free(sket);
    nz = -1;
    goto END;
  }

  s[0].kappa = s[1].kappa;
  s[0].j = s[1].j;
  s[0].kl = s[1].kl;
  sbra[0].shellJ = s[0].j;
  
  nz = -1;
  nkk = AngularZ(&r, &k, 1, n_shells, sbra, sket, s, s+1);
  if (nkk == 1) {
    if (IsOdd(phase+(j1+s[0].j-j2)/2)) *r = -(*r);
    (*r) *= sqrt(s[0].j + 1.0);
    if (fabs(*r) > angz_cut) {
      kb = OrbitalIndex(s[1].n, s[1].kappa, 0.0);
      afb.kb = kb;
      afb.coeff = *r;
      nz = 1;
    }
  }

  free(sbra);
  free(sket);

 END:
  if (nz != 1) {
    angz_datum->nz = -1;
    nz = -1;
  } else {
    angz_datum->angz = malloc(sizeof(ANGULAR_ZFB));
    angz_datum->nz = 1;
    memcpy(angz_datum->angz, &afb, sizeof(ANGULAR_ZFB));
    (*ang) = (ANGULAR_ZFB *) angz_datum->angz;
    nz = 1;
  }
 
#ifdef PERFORM_STATISTICS
  stop = clock();
  timing.angzfb_states += stop-start;
#endif
  return nz;

}

int AngularZxZMixStates(ANGULAR_ZxZMIX **ang, STATE *slow, STATE *sup) {
  int nz, n, i, j, im, p;
  int n_shells;
  INTERACT_SHELL s[4];
  SHELL *bra;
  SHELL_STATE *sbra, *sket;
  CONFIG *clow, *cup;
  int phase; 

  clow = GetConfig(slow);
  cup = GetConfig(sup);
  n_shells = GetInteract(&phase, s, &bra, &sbra, &sket,
			 clow, slow->kstate, cup, sup->kstate,
			 slow, sup);

  if (n_shells <= 0) return -1;

  n = 0;
  nz = ANGZ_BLOCK;
  (*ang) = malloc(sizeof(ANGULAR_ZxZMIX)*nz);
  
  if (s[0].index >= 0 && s[2].index >= 0) {
    AddToAngularZxZ(&n, &nz, ang, n_shells, phase, sbra, sket, s, 0);
  } else if (s[0].index >= 0) {
    for (i = 0; i < n_shells; i++) {      
      s[2].index = n_shells - i - 1;
      if (s[2].index == s[0].index) continue;
      s[3].index = s[2].index;
      s[2].n = bra[i].n;
      s[3].n = s[2].n;
      s[2].kappa = bra[i].kappa;
      s[3].kappa = s[2].kappa;
      s[2].j = GetJ(bra+i);
      s[3].j = s[2].j;
      s[2].kl = GetL(bra+i);
      s[3].kl = s[2].kl;
      s[2].nq_bra = GetNq(bra+i);
      if (s[2].index == s[0].index) {
	s[2].nq_ket = s[2].nq_bra - 1;
      } else if (s[2].index == s[1].index) {
	s[2].nq_ket = s[2].nq_bra + 1;
      } else {
	s[2].nq_ket = s[2].nq_bra;
      }
      if (s[2].nq_bra < 0 || s[2].nq_ket < 0 ||
	  s[2].nq_bra > s[2].j+1 || s[2].nq_ket > s[2].j+1.0) {
	continue;
      }
      s[3].nq_bra = s[2].nq_bra;
      s[3].nq_ket = s[2].nq_ket;
      AddToAngularZxZ(&n, &nz, ang, n_shells, phase, sbra, sket, s, 0);
    }
  } else {   
    for (i = 0; i < n_shells; i++) {
      s[0].index = n_shells - i - 1;
      s[1].index = s[0].index;
      s[0].n = bra[i].n;
      s[1].n = s[0].n;
      s[0].kappa = bra[i].kappa;
      s[1].kappa = s[0].kappa;
      s[0].j = GetJ(bra+i);
      s[1].j = s[0].j;
      s[0].kl = GetL(bra+i);
      s[1].kl = s[0].kl;
      s[0].nq_bra = GetNq(bra+i);
      s[0].nq_ket = s[0].nq_bra;
      s[1].nq_bra = s[0].nq_bra;
      s[1].nq_ket = s[1].nq_bra;  
      for (j = 0; j <= i; j++) {
	s[2].index = n_shells - j - 1;
	if (s[2].index == s[0].index) continue;
	s[3].index = s[2].index;
	s[2].n = bra[j].n;
	s[3].n = s[2].n;
	s[2].kappa = bra[j].kappa;
	s[3].kappa = s[2].kappa;
	s[2].j = GetJ(bra+j);
	s[3].j = s[2].j;
	s[2].kl = GetL(bra+j);
	s[3].kl = s[2].kl;
	s[2].nq_bra = GetNq(bra+j);
	s[2].nq_ket = s[2].nq_bra;
	s[3].nq_bra = s[2].nq_bra;
	s[3].nq_ket = s[3].nq_bra;
	AddToAngularZxZ(&n, &nz, ang, n_shells, phase, sbra, sket, s, 0);
      }
    }
  }
  free(sbra);
  free(sket);
  if (n < nz) (*ang) = realloc(*ang, n*sizeof(ANGULAR_ZxZMIX));

  return n;
}

int AngularZxZFreeBoundStates(ANGULAR_ZxZMIX **ang, 
			      STATE *slow, STATE *sup) {  
  int nz, n, j1, j2, i, im, p;
  int n_shells;
  INTERACT_SHELL s[4], sp[4];
  SHELL *bra;
  SHELL_STATE *sbra, *sket;
  CONFIG *clow, *cup, clowf;
  int phase, jmin, jmax, jf;
  SHELL sh[200];
  SHELL_STATE cs[200];
  ANGZ_DATUM *angz_datum;
  int index[6];
  clock_t start, stop; 
  
#ifdef PERFORM_STATISTICS
  start = clock();
#endif

  index[1] = slow->kgroup;
  index[0] = sup->kgroup;
  index[3] = slow->kcfg;
  index[2] = sup->kcfg;
  index[5] = slow->kstate;
  index[4] = sup->kstate;
  angz_datum = (ANGZ_DATUM *) MultiSet(angzxz_array, index, NULL);
  nz = angz_datum->nz;
  if (nz < 0) { 
#ifdef PERFORM_STATISTICS
    stop = clock();
    timing.angzxzfb_states += stop-start;
#endif
    return -1;
  }
  
  if (nz > 0) {
    (*ang) = (ANGULAR_ZxZMIX *) angz_datum->angz;
#ifdef PERFORM_STATISTICS
    stop = clock();
    timing.angzxzfb_states += stop-start;
#endif
    return nz;
  }

  clowf.shells = sh;
  clowf.csfs = cs;
  clow = GetConfig(slow);
  j1 = clow->csfs[slow->kstate].totalJ;
  
  clowf.n_shells = clow->n_shells+1;
  memcpy(clowf.shells+1, clow->shells, sizeof(SHELL)*clow->n_shells);
  clowf.shells[0].n = 1000;
  clowf.shells[0].nq = 1;
  clowf.shells[0].kappa = -1;
  clowf.n_csfs = 1;
  memcpy(clowf.csfs+1, clow->csfs+slow->kstate, 
	 sizeof(SHELL_STATE)*clow->n_shells);	
  cup = GetConfig(sup);
  j2 = cup->csfs[sup->kstate].totalJ;

  clowf.csfs[0].shellJ = 1;
  clowf.csfs[0].totalJ = j2;
  clowf.csfs[0].nu = 1;
  clowf.csfs[0].Nr = 0;   
  n_shells = GetInteract(&phase, sp, &bra, &sbra, &sket,
			 &clowf, 0, cup, sup->kstate,
			 slow, sup);

  
  if (n_shells <= 0) {
#ifdef PERFORM_STATISTICS
  stop = clock();
  timing.angzfb_states += stop-start;
#endif
    return -1;
  } 

  n = 0;
  nz = ANGZ_BLOCK;
  (*ang) = malloc(sizeof(ANGULAR_ZxZMIX)*nz);

  jmin = abs(j2-j1);
  jmax = j1+j2;

  for (jf = jmin; jf <= jmax; jf += 2) {    
    memcpy(s, sp, sizeof(INTERACT_SHELL)*4);
    s[0].j = jf;
    s[0].kl = jf+1;
    s[0].kappa = GetKappaFromJL(s[0].j, s[0].kl);
    sbra[0].shellJ = s[0].j;

    if (s[2].index >= 0) {
      AddToAngularZxZ(&n, &nz, ang, n_shells, phase, sbra, sket, s, 1);
    } else {      
      for (i = 0; i < n_shells; i++) {      
	s[2].index = n_shells - i - 1;
	if (s[2].index == s[0].index) continue;
	s[3].index = s[2].index;
	s[2].n = bra[i].n;
	s[3].n = s[2].n;
	s[2].kappa = bra[i].kappa;
	s[3].kappa = s[2].kappa;
	s[2].j = GetJ(bra+i);
	s[3].j = s[2].j;
	s[2].kl = GetL(bra+i);
	s[3].kl = s[2].kl;
	s[2].nq_bra = GetNq(bra+i);
	if (s[2].index == s[0].index) {
	  s[2].nq_ket = s[2].nq_bra - 1;
	} else if (s[2].index == s[1].index) {
	  s[2].nq_ket = s[2].nq_bra + 1;
	} else {
	  s[2].nq_ket = s[2].nq_bra;
	}
	s[3].nq_bra = s[2].nq_bra;
	s[3].nq_ket = s[2].nq_ket;
	AddToAngularZxZ(&n, &nz, ang, n_shells, phase, sbra, sket, s, 1);
      }
    }
  }
  
  
  free(sbra);
  free(sket);

  if (n <= 0) {
    angz_datum->nz = -1;
    free(*ang);
    return -1;
  } else {
    if (n < nz) (*ang) = realloc(*ang, sizeof(ANGULAR_ZxZMIX)*n);
    angz_datum->angz = (void *) *ang;
    angz_datum->nz = n;
  }

#ifdef PERFORM_STATISTICS
  stop = clock();
  timing.angzfb_states += stop-start;
#endif
  return n;
}

int AddToAngularZxZ(int *n, int *nz, ANGULAR_ZxZMIX **ang, 
		    int n_shells, int phase, SHELL_STATE *sbra, 
		    SHELL_STATE *sket, INTERACT_SHELL *s, int m) {
  int nkk, *k, p, im;
  double *r;
  int orb0, orb1;
  int kk0, kk1;
  
  nkk = AngularZxZ0(&r, &k, 0, n_shells, sbra, sket, s);
  if (nkk > 0) {
    if (m == 0) {
      orb0 = OrbitalIndex(s[0].n, s[0].kappa, 0.0);
    } else {
      orb0 = s[0].j;
    }
    orb1 = OrbitalIndex(s[1].n, s[1].kappa, 0.0);
    kk0 = OrbitalIndex(s[2].n, s[2].kappa, 0.0);
    kk1 = OrbitalIndex(s[3].n, s[3].kappa, 0.0);
    for (p = 0; p < nkk; p++) {
      if (fabs(r[p]) > angz_cut) {
	if (IsOdd(phase)) r[p] = -r[p];
	im = AddToAngularZxZMix(n, nz, ang, k[p], 
				orb0, orb1, kk0, kk1, r[p]);
	if (im < 0) {
	  free(r);
	  free(k);
	  free(*ang);
	  return -1;
	}
      }
    }
    free(r);
    free(k);
  }
  
  return 0;
}

int AddToAngularZxZMix(int *n, int *nz, ANGULAR_ZxZMIX **ang, 
		       int k, int k0, int k1, int k2, int k3, double r) {
  int im;
  clock_t start, stop;

#ifdef PERFORM_STATISTICS
  start = clock();
#endif

  for (im = 0; im < *n; im++) {
    if ((*ang)[im].k == k &&
	(*ang)[im].k0 == k0 &&
	(*ang)[im].k1 == k1 &&
	(*ang)[im].k2 == k2 &&
	(*ang)[im].k3 == k3) break;
  }
  if (im < *n) {
    (*ang)[im].coeff += r;
  } else {
    (*n)++;
    if (*n > *nz) {
      *nz += ANGZ_BLOCK;
      *ang = realloc((*ang), (*nz)*sizeof(ANGULAR_ZxZMIX));
      if (!(*ang)) {
	printf("Can't enlarge AngularZxZMix array\n");
	return -1;
      }
    }
    (*ang)[im].k = k;
    (*ang)[im].k0 = k0;
    (*ang)[im].k1 = k1;
    (*ang)[im].k2 = k2;
    (*ang)[im].k3 = k3;
    (*ang)[im].coeff = r;
  }
#ifdef PERFORM_STATISTICS
  stop = clock();
  timing.add_angzxz += stop -start;
#endif

  return 0;
}

int AngularZFreeBound(ANGULAR_ZFB **ang, int lower, int upper) {
  int i, j, m; 
  int nz, n;
  double r0;
  STATE *slow, *sup;
  SYMMETRY *sym1, *sym2;
  int hlow, hup;
  HAMILTON *h1, *h2;
  LEVEL *lev1, *lev2;
  double mix1, mix2;
  int kg, jf, i1, i2, kb, ia, j1, j2;
  ANGULAR_ZFB *ang_sub;
  clock_t start, stop;

#ifdef PERFORM_STATISTICS
  start = clock();
#endif

  lev1 = GetLevel(lower);
  lev2 = GetLevel(upper);
  hlow = lev1->ham_index;
  hup = lev2->ham_index;
  h1 = GetHamilton(hlow);
  h2 = GetHamilton(hup);
  sym1 = GetSymmetry(h1->pj);
  sym2 = GetSymmetry(h2->pj);
  j1 = h1->pj;
  j2 = h2->pj;
  DecodePJ(j1, NULL, &j1);
  DecodePJ(j2, NULL, &j2);
  
  sup = (STATE *) ArrayGet(&(sym2->states), h2->basis[lev2->major_component]);
  if (sup->kgroup < 0) {
    n = 1;
    nz = 1;
    (*ang) = malloc(sizeof(ANGULAR_ZFB));
    (*ang)->kb = sup->kcfg;
    (*ang)->coeff = 0.0;
    i2 = lev2->mix_index;
    i2 = h2->dim * (i2 + 1);
    for (j = 0; j < h2->dim; j++) {
      mix2 = h2->mixing[i2];
      if (fabs(mix2) < mix_cut) {
	i2++;
	continue;
      }
      sup = (STATE *) ArrayGet(&(sym2->states), h2->basis[j]);
      kg = sup->kgroup;
      kg = -kg-1;
      if (kg == lower) {
	(*ang)->kb = sup->kcfg;
	kb = GetOrbital(sup->kcfg)->kappa;
	jf = GetJFromKappa(kb);
	mix2 *= sqrt(j2 + 1.0);
	if (IsEven((j2+jf-j1)/2)) mix2 = -mix2;
	(*ang)->coeff = mix2;
	break;
      }
      i2++;
    }
    
    if (fabs((*ang)->coeff) < angz_cut) {
      nz = 0;
      n = 0;
      free(*ang);
    }
  } else {
    n = 0;
    nz = ANGZ_BLOCK;
    (*ang) = malloc(sizeof(ANGULAR_ZFB)*nz);
    i1 = lev1->mix_index;
    i1 = h1->dim * (i1 + 1);
    
    for (i = 0; i < h1->dim; i++) {
      mix1 = h1->mixing[i1];
      if (fabs(mix1) < mix_cut) {
	i1++;
	continue;
      }
      i2 = lev2->mix_index;
      i2 = h2->dim * (i2 + 1);
      slow = (STATE *) ArrayGet(&(sym1->states), h1->basis[i]);
      for (j = 0; j < h2->dim; j++) {
	mix2 = h2->mixing[i2];
	if (fabs(mix2) < mix_cut) {
	  i2++;
	  continue;
	}
	sup = (STATE *) ArrayGet(&(sym2->states), h2->basis[j]);
	r0 = mix1*mix2;
	m = AngularZFreeBoundStates(&ang_sub, slow, sup);
	if (m == 1) {
	  kb = ang_sub->kb;
	  r0 *= ang_sub->coeff;
	  if (fabs(r0) < angz_cut) {
	    i2++;
	    continue;
	  }
	  for (ia = 0; ia < n; ia++) {
	    if ((*ang)[ia].kb == kb) break;
	  }
	  if (ia == n) {
	    n++;
	    (*ang)[ia].kb = kb;
	    (*ang)[ia].coeff = r0;
	  } else {
	    (*ang)[ia].coeff += r0;
	  }
	}
	i2++;
      }
      i1++;
    }
    if (n == 0) {
      free(*ang);
    }
  }

#ifdef PERFORM_STATISTICS
    stop = clock();
    timing.angz_fb += stop-start;
#endif
  return n;
}

int GetBaseJ(STATE *s) {
  int ih;
  if (s->kgroup >= 0) return -1;
  ih = s->kgroup;
  ih = -ih-1;
  ih = GetLevel(ih)->ham_index;
  ih = GetHamilton(ih)->pj;
  DecodePJ(ih, NULL, &ih);
  return ih;
}

int AngularZMix(ANGULAR_ZMIX **ang, int lower, int upper, 
		int mink, int maxk) {
  int i, j, j1, j2, jb1, jb2, kg1, kg2, jlow, jup, kb1, kb2;
  int nz, n;
  int i1, i2, im;
  double r0;
  STATE *slow, *sup;
  SYMMETRY *sym1, *sym2;
  HAMILTON *h1, *h2;
  LEVEL *lev1, *lev2;
  int hlow, hup;
  int ik, kmin, kmax, m;
  int nz_sub, nfb;
  ANGULAR_ZMIX *ang_sub;
  ANGULAR_ZFB *afb;
  double mix1, mix2, sqrt_j12, a;
  int ignore_ryd = 0;
  clock_t start, stop;

#ifdef PERFORM_STATISTICS
  start = clock();
#endif

  lev1 = GetLevel(lower);
  lev2 = GetLevel(upper);
  hlow = lev1->ham_index;
  hup = lev2->ham_index;
  h1 = GetHamilton(hlow);
  h2 = GetHamilton(hup);
  sym1 = GetSymmetry(h1->pj);
  sym2 = GetSymmetry(h2->pj);
  j1 = h1->pj;
  j2 = h2->pj;
  DecodePJ(j1, NULL, &j1);
  DecodePJ(j2, NULL, &j2);
  sqrt_j12 = sqrt((j1+1.0)*(j2+1.0));
  kmin = abs(j1-j2);
  kmax = j1+j2;
  if (IsOdd(kmin)) kmin++;
  if (mink >= 0) kmin = Max(kmin, mink);
  if (maxk >= 0) kmax = Min(kmax, maxk);
  if (kmax < kmin) {
#ifdef PERFORM_STATISTICS
    stop = clock();
    timing.angz_mix += stop-start;
#endif
    return 0;
  }

  slow = (STATE *) ArrayGet(&(sym1->states), h1->basis[lev1->major_component]);
  sup = (STATE *) ArrayGet(&(sym2->states), h2->basis[lev2->major_component]);
  kg1 = slow->kgroup;
  kg2 = sup->kgroup;

  if (kg1 < 0) {
    kb1 = slow->kcfg;
    n = GetOrbital(kb1)->n;
    if (rydberg_ignored > 0 && rydberg_ignored < n) ignore_ryd =1;
  } else if (kg2 < 0) {
    kb2 = sup->kcfg;
    n = GetOrbital(kb2)->n;
    if (rydberg_ignored > 0 && rydberg_ignored < n) ignore_ryd =1;
  }

  if (kg1 < 0 && kg2 < 0) {
    n = 0;
    nz = ANGZ_BLOCK;
    (*ang) = malloc(sizeof(ANGULAR_ZMIX)*nz);
    
    i1 = lev1->mix_index;
    i1 = h1->dim * (i1 + 1);
    
    for (i = 0; i < h1->dim; i++) {
      mix1 = h1->mixing[i1];
      if (fabs(mix1) < mix_cut) {
	i1++;
	continue;
      }
      i2 = lev2->mix_index;
      i2 = h2->dim * (i2 + 1);
      slow = (STATE *) ArrayGet(&(sym1->states), h1->basis[i]);
      jlow = GetBaseJ(slow);
      kg1 = slow->kgroup;
      kg1 = -kg1-1;
      kb1 = slow->kcfg;
      jb1 = GetOrbital(kb1)->kappa;
      jb1 = GetJFromKappa(jb1);
      for (j = 0; j < h2->dim; j++) {
	mix2 = h2->mixing[i2];
	if (fabs(mix2) < mix_cut) {
	  i2++;
	  continue;
	}
	a = mix1*mix2;
	sup = (STATE *) ArrayGet(&(sym2->states), h2->basis[j]);
	jup = GetBaseJ(sup);
	kg2 = sup->kgroup;
	kg2 = -kg2-1;
	kb2 = sup->kcfg;
	jb2 = GetOrbital(kb2)->kappa;
	jb2 = GetJFromKappa(jb2);

	if (ignore_ryd && kb1 != kb2) {
	  i2++;
	  continue;
	}

	if (kg1 == kg2) {
	  for (ik = kmin; ik <= kmax; ik += 2) {
	    r0 = W6j(jlow, jb1, j1, ik, j2, jb2);
	    if (fabs(r0) < angz_cut) continue;
	    r0 *= a*sqrt_j12;
	    if (IsEven((j1+jb2-jlow-ik)/2+j2)) 
	      r0 = -r0;
	    im = AddToAngularZMix(&n, &nz, ang, ik, kb1, kb2, r0);
	  }
	}
	if (kb1 == kb2){	  
	  nz_sub = AngularZMix(&ang_sub, kg1, kg2, kmin, kmax);
	  if (nz_sub <= 0) {
	    i2++;
	    continue;
	  }
	  for (m = 0; m < nz_sub; m++) {
	    r0 = W6j(jlow, jb1, j1, j2, ang_sub[m].k, jup);
	    if (fabs(r0) < angz_cut) continue;
	    r0 *= a*ang_sub[m].coeff;
	    r0 *= sqrt_j12;
	    if (IsOdd((jlow+jb1+j2+ang_sub[m].k)/2)) r0 = -r0;
	    im = AddToAngularZMix(&n, &nz, ang, ang_sub[m].k, 
				  ang_sub[m].k0, ang_sub[m].k1, r0);
	  }
	  if (nz_sub > 0) free(ang_sub);
	}
	i2++;
      }
      i1++;
    }
  } else if (kg1 < 0 && !ignore_ryd) {        
    nz = ANGZ_BLOCK;
    n = 0;
    (*ang) = malloc(sizeof(ANGULAR_ZMIX)*nz);
    i1 = lev1->mix_index;
    i1 = h1->dim * (i1 + 1); 
    for (i = 0; i < h1->dim; i++) {
      mix1 = h1->mixing[i1];
      if (fabs(mix1) < mix_cut) {
	i1++;
	continue;
      }
      slow = (STATE *) ArrayGet(&(sym1->states), h1->basis[i]);
      kb1 = slow->kcfg;
      jb1 = GetOrbital(kb1)->kappa;
      jb1 = GetJFromKappa(jb1);
      jlow = GetBaseJ(slow);
      kg1 = slow->kgroup;
      kg1 = -kg1-1;
      nfb = AngularZFreeBound(&afb, kg1, upper);
      for (m = 0; m < nfb; m++) {
	jb2 = GetOrbital(afb[m].kb)->kappa;
	jb2 = GetJFromKappa(jb2);
	for (ik = kmin; ik <= kmax; ik += 2) {
	  r0 = W6j(jlow, jb1, j1, ik, j2, jb2);
	  r0 *= mix1*afb[m].coeff*sqrt(j1+1.0);
	  if (fabs(r0) < angz_cut) continue;
	  if (IsOdd((j1+j2-ik)/2)) r0 = -r0;
	  im = AddToAngularZMix(&n, &nz, ang, ik, kb1, afb[m].kb, r0);
	}
      }
      if (nfb > 0) free(afb);
      i1++;
    }
  } else if (kg2 < 0 && !ignore_ryd) {    
    nz = ANGZ_BLOCK;
    n = 0;
    (*ang) = malloc(sizeof(ANGULAR_ZMIX)*nz);
    i2 = lev2->mix_index;
    i2 = h2->dim * (i2 + 1);
    for (j = 0; j < h2->dim; j++) {
      mix2 = h2->mixing[i2];
      if (fabs(mix2) < mix_cut) {
	i2++;
	continue;
      }
      sup = (STATE *) ArrayGet(&(sym2->states), h2->basis[j]);
      kb2 = sup->kcfg;
      jb2 = GetOrbital(kb2)->kappa;
      jb2 = GetJFromKappa(jb2);
      jup = GetBaseJ(sup);
      kg2 = sup->kgroup;
      kg2 = -kg2-1;
      nfb = AngularZFreeBound(&afb, kg2, lower);
      for (m = 0; m < nfb; m++) {
	jb1 = GetOrbital(afb[m].kb)->kappa;
	jb1 = GetJFromKappa(jb1);
	for (ik = kmin; ik <= kmax; ik += 2) {
	  r0 = W6j(jup, jb2, j2, ik, j1, jb1);
	  r0 *= mix2*afb[m].coeff*sqrt(j2+1.0);
	  if (fabs(r0) < angz_cut) continue;
	  if (IsOdd((2*j1-ik+jb1-jb2)/2)) r0 = -r0;
	  im = AddToAngularZMix(&n, &nz, ang, ik, afb[m].kb, kb2, r0);
	}
      }
      if (nfb > 0) free(afb);
      i2++;
    }
  } else {
    nz = ANGZ_BLOCK;
    n = 0;
    (*ang) = malloc(sizeof(ANGULAR_ZMIX)*nz);
    if (!(*ang)) return -1;
  
    i1 = lev1->mix_index;
    i1 = h1->dim * (i1 + 1);

    for (i = 0; i < h1->dim; i++) {
      mix1 = h1->mixing[i1];
      if (fabs(mix1) < mix_cut) {
	i1++;
	continue;
      }
      i2 = lev2->mix_index;
      i2 = h2->dim * (i2 + 1);
      slow = (STATE *) ArrayGet(&(sym1->states), h1->basis[i]);
      for (j = 0; j < h2->dim; j++) {
	mix2 = h2->mixing[i2];
	if (fabs(mix2) < mix_cut) {
	  i2++;
	  continue;
	}
	a = mix1*mix2;
	sup = (STATE *) ArrayGet(&(sym2->states), h2->basis[j]);
	nz_sub = AngularZMixStates(&ang_sub, slow, sup);
	for (m = 0; m < nz_sub; m++) {
	  r0 = ang_sub[m].coeff*a;
	  if (fabs(r0) < angz_cut) continue;
	  if (ang_sub[m].k > kmax || ang_sub[m].k < kmin) continue;
	  im = AddToAngularZMix(&n, &nz, ang, ang_sub[m].k, 
				ang_sub[m].k0, ang_sub[m].k1, r0);
	}
	/*if (nz_sub > 0) free(ang_sub);*/
	i2++;
      }
      i1++;
    }
  }

#ifdef PERFORM_STATISTICS
    stop = clock();
    timing.angz_mix += stop - start;
#endif

  if (n == 0) free(*ang);
  return n;
}

int AngularZxZFreeBound(ANGULAR_ZxZMIX **ang, int lower, int upper) {
  int i, j, j1, j2;
  int nz, n;
  int kg, jf, i1, i2;
  double mix1, mix2;
  int hlow, hup;
  STATE *slow, *sup;
  SYMMETRY *sym1, *sym2;
  HAMILTON *h1, *h2;
  LEVEL *lev1, *lev2;
  int kb, jb, jup, jmin, jmax, orb0, orb1;
  double r, r0, sqrt_j2;
  int nz_sub;
  ANGULAR_ZxZMIX *ang_sub;
  ANGULAR_ZMIX *ang_z;
  int im, m;
  clock_t start, stop;
  	
#ifdef PERFORM_STATISTICS
  start = clock();
#endif

  lev1 = GetLevel(lower);
  lev2 = GetLevel(upper);
  hlow = lev1->ham_index;
  hup = lev2->ham_index;
  h1 = GetHamilton(hlow);
  h2 = GetHamilton(hup);
  sym1 = GetSymmetry(h1->pj);
  sym2 = GetSymmetry(h2->pj);
  j1 = h1->pj;
  j2 = h2->pj;
  DecodePJ(j1, NULL, &j1);
  DecodePJ(j2, NULL, &j2);
  sqrt_j2 = sqrt(j2+1.0);

  nz = ANGZ_BLOCK;
  n = 0;
  (*ang) = malloc(sizeof(ANGULAR_ZxZMIX)*nz);
  if (!(*ang)) return -1;
    
  sup = (STATE *) ArrayGet(&(sym2->states), h2->basis[lev2->major_component]);
  if (sup->kgroup < 0) {  
    i2 = lev2->mix_index;
    i2 = h2->dim * (i2 + 1);
    for (j = 0; j < h2->dim; j++) {
      mix2 = h2->mixing[i2];
      if (fabs(mix2) < mix_cut) {
	i2++;
	continue;
      }
      sup = (STATE *) ArrayGet(&(sym2->states), h2->basis[j]);
      kb = sup->kcfg;
      jb = GetOrbital(kb)->kappa;
      jb = GetJFromKappa(jb);
      jup = GetBaseJ(sup);
      kg = sup->kgroup;
      kg = -kg-1;
      nz_sub = AngularZMix(&ang_z, lower, kg, -1, -1);
      /*      printf("%d %10.3E %d %d\n", j, mix2, nz_sub, kg); */
      if (nz_sub <= 0) {
	i2++;
	continue;
      }
      for (i = 0; i < nz_sub; i++) {
	r = ang_z[i].coeff*mix2*sqrt_j2;
	/*	printf("%d %10.3E %10.3E\n", i, ang_sub[i].coeff, r); */
	if (fabs(r) < angz_cut) continue;	
	orb0 = ang_z[i].k0;
	orb1 = ang_z[i].k1;
	jmin = abs(j1-j2);
	jmin = Max(jmin, abs(jb-ang_z[i].k));
	jmax = j1+j2;
	jmax = Min(jmax, jb+ang_z[i].k);
	/*	printf("%d %d %d %d\n", i, jmin, jmax, ang_sub[i].k); */
	for (jf = jmin; jf <= jmax; jf += 2) {
	  r0 = W6j(j1, jf, j2, jb, jup, ang_z[i].k);
	  if (fabs(r0) < angz_cut) continue;
	  if (IsOdd((jup+jf+j2)/2)) r0 = -r0;
	  r0 *= r;
	  im = AddToAngularZxZMix(&n, &nz, ang, ang_z[i].k, 
				  jf, kb, orb0, orb1, r0);
	}
      }
      i2++;
      free(ang_z);
    }
  } else { 

    i1 = lev1->mix_index;
    i1 = h1->dim * (i1 + 1);

    for (i = 0; i < h1->dim; i++) {
      mix1 = h1->mixing[i1];
      if (fabs(mix1) < mix_cut) {
	i1++;
	continue;
      }
      i2 = lev2->mix_index;
      i2 = h2->dim * (i2 + 1);
      slow = (STATE *) ArrayGet(&(sym1->states), h1->basis[i]);
      for (j = 0; j < h2->dim; j++) {
	mix2 = h2->mixing[i2];
	if (fabs(mix2) < mix_cut) {
	  i2++;
	  continue;
	}
	r = mix1*mix2;
	sup = (STATE *) ArrayGet(&(sym2->states), h2->basis[j]); 

	nz_sub = AngularZxZFreeBoundStates(&ang_sub, slow, sup); 

	for (m = 0; m < nz_sub; m++) {
	  r0 = ang_sub[m].coeff*r;
	  im = AddToAngularZxZMix(&n, &nz, ang, 
				  ang_sub[m].k, ang_sub[m].k0,
				  ang_sub[m].k1, ang_sub[m].k2,
				  ang_sub[m].k3, r0); 	  
	}
	/*if (nz_sub > 0) free(ang_sub);*/

	i2++;
      }
      i1++;
    }
  }

  if (n == 0) free(*ang);

#ifdef PERFORM_STATISTICS
    stop = clock();
    timing.angzxz_fb += stop - start;
#endif
	
  return n;
}
  
int CompareAngularZMix(const void *c1, const void *c2) {
  ANGULAR_ZMIX *a1, *a2;
  double coeff1, coeff2;
  int r;

  a1 = (ANGULAR_ZMIX *) c1;
  a2 = (ANGULAR_ZMIX *) c2;
  coeff1 = fabs(a1->coeff);
  coeff2 = fabs(a2->coeff);
  
  if (coeff1 > coeff2) r = -1;
  else if (coeff1 < coeff2) r = 1;
  else r = 0;

  return r;
}

int AddToAngularZMix(int *n, int *nz, ANGULAR_ZMIX **ang,
		     int k, int k0, int k1, double coeff) {
  int im;
  clock_t start, stop;

#ifdef PERFORM_STATISTICS
  start = clock();
#endif

  for (im = 0; im < *n; im++) {
    if ((*ang)[im].k == k &&
	(*ang)[im].k0 == k0 &&
	(*ang)[im].k1 == k1) break;
  }
  if (im < *n) {
    (*ang)[im].coeff += coeff;
  } else {
    (*n)++;
    if (*n > *nz) {
      *nz += ANGZ_BLOCK;
      *ang = realloc((*ang), (*nz)*sizeof(ANGULAR_ZMIX));
      if (!(*ang)) {
	printf("Can't enlarge AngularZMix array\n");
	return -1;
      }
    }
    (*ang)[im].k = k;
    (*ang)[im].k0 = k0;
    (*ang)[im].k1 = k1;
    (*ang)[im].coeff = coeff;
  }
  
#ifdef PERFORM_STATISTICS
  stop = clock();
  timing.add_angz += stop -start;
#endif

  return 0;
}

void _FreeAngZDatum(void *p) {
  ANGZ_DATUM *ap;
  ap = (ANGZ_DATUM *) p;
  if (ap->nz > 0) free(ap->angz);
  ap->nz = 0;
}

int _FreeAngZ(int g, MULTI *ma) {  
  ARRAY *a, *b;  
  int ndim;
  
  a = ma->array;
  if (a == NULL) return 0;
  ndim = ma->ndim;
  if (g < 0) {
    MultiFreeData(a, ndim, _FreeAngZDatum);
  } else {
    b = (ARRAY *) ArrayGet(a, g);
    if (b) {
      MultiFreeData(b, ndim - 1, _FreeAngZDatum);
    }
  }
  return 0;
}

int FreeAngZ(int g, int which_array) {
  ARRAY *a, *b;  
  int ndim;

  if (which_array == 0) {
    _FreeAngZ(g, angz_array);
  } else if (which_array > 0) {
    _FreeAngZ(g, angzxz_array);
  } else {
    _FreeAngZ(g, angz_array);
    _FreeAngZ(g, angzxz_array);
  }
  return 0;
}

int InitStructure() {
  int i;
  int ndim = 6;
  int blocks[6] = {10, 5, 10, 5, 20, 50};
  n_hamiltons = 0;
  hamiltons = malloc(sizeof(ARRAY));
  if (!hamiltons) return -1;
  ArrayInit(hamiltons, sizeof(HAMILTON), HAMILTONS_BLOCK);

  n_levels = 0;
  levels = malloc(sizeof(ARRAY));
  if (!levels) return -1;
  ArrayInit(levels, sizeof(LEVEL), LEVELS_BLOCK);

  angz_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(angz_array, sizeof(ANGZ_DATUM), ndim, blocks);

  angzxz_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(angzxz_array, sizeof(ANGZ_DATUM), ndim, blocks);

  return 0;
}


