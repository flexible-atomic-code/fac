#include <time.h>

#include "structure.h"
#include "cf77.h"

static char *rcsid="$Id: structure.c,v 1.47 2003/05/14 02:31:27 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

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

static HAMILTON _ham = {0, 0, 0, 0, 0, 0, 0, 0, 0, 
			NULL, NULL, NULL, NULL, NULL};

static ARRAY *levels;
static int n_levels = 0;

static MULTI *angz_array;
static MULTI *angzxz_array;

static int ncorrections = 0;
static ARRAY *ecorrections;

static int rydberg_ignored = 0;
static double angz_cut = ANGZCUT;
static double mix_cut = MIXCUT;


#ifdef PERFORM_STATISTICS 
static STRUCT_TIMING timing = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
int GetStructTiming(STRUCT_TIMING *t) {
  memcpy(t, &timing, sizeof(timing));
  return 0;
}
#endif

static void InitAngzDatum(void *p, int n) {
  ANGZ_DATUM *d;
  int i;
  
  d = (ANGZ_DATUM *) p;
  for (i = 0; i < n; i++) {
    d[i].ns = 0;
  }
}

int SetAngZCut(double cut) {
  angz_cut = cut;
  return 0;
}

int SetMixCut(double cut) {
  mix_cut = cut;
  return 0;
}

int SetAngZOptions(int n, double mix, double cut) {
  rydberg_ignored = n;
  mix_cut = mix;
  angz_cut = cut;
  return 0;
}

int ConstructHamilton(int isym, int k, int *kg, int kp, int *kgp) {
  int i, j, t, jp;
  HAMILTON *h;
  ARRAY *st;
  STATE *s;
  SYMMETRY *sym;
  double r;
#if (FAC_DEBUG >= DEBUG_STRUCTURE) 
  char name[LEVEL_NAME_LEN];
#endif
#ifdef PERFORM_STATISTICS
  clock_t start, stop;
  start = clock();
#endif

  j = 0;
  jp = 0;
  sym = GetSymmetry(isym);
  if (sym == NULL) return -1;
  st = &(sym->states);
  for (t = 0; t < sym->n_states; t++) {
    s = (STATE *) ArrayGet(st, t);
    if (InGroups(s->kgroup, k, kg)) j++;
    else if (kp > 0 && InGroups(s->kgroup, kp, kgp)) jp++;
  }
  
  if (j == 0) return -1;

  h = &_ham;
  h->pj = isym;

  h->dim = j;
  h->n_basis = jp+j;
  t = j*(j+1)/2;
  h->hsize = t + (h->dim*jp) + jp;

  if (h->basis == NULL) {
    h->n_basis0 = h->n_basis;
    h->basis = (int *) malloc(sizeof(int)*(h->n_basis));
  } else if (h->n_basis > h->n_basis0) {
    h->n_basis0 = h->n_basis;
    h->basis = (int *) realloc(h->basis, sizeof(int)*h->n_basis);
  }
  if (!(h->basis)) goto ERROR;

  if (h->hamilton == NULL) {
    h->hsize0 = h->hsize;
    h->hamilton = (double *) malloc(sizeof(double)*h->hsize);
  } else if (h->hsize > h->hsize0) {
    h->hsize0 = h->hsize;
    h->hamilton = (double *) realloc(h->hamilton, sizeof(double)*h->hsize);
  }
  if (!(h->hamilton)) goto ERROR;

  j = 0;  
  for (t = 0; t < sym->n_states; t++) {
    s = (STATE *) ArrayGet(st, t);
    if (InGroups(s->kgroup, k, kg)) {
      h->basis[j] = t;
      j++;
    }
  }
  if (jp > 0) {  
    for (t = 0; t < sym->n_states; t++) {
      s = (STATE *) ArrayGet(st, t);
      if (kp > 0 && InGroups(s->kgroup, kp, kgp)) {
	h->basis[j] = t;
	j++;
      }
    }
  }

  for (j = 0; j < h->dim; j++) {
    t = j*(j+1)/2;
    for (i = 0; i <= j; i++) {
      r = HamiltonElement(isym, h->basis[i], h->basis[j]);
      h->hamilton[i+t] = r;
    }
  }

  if (jp > 0) {
    t = ((h->dim+1)*(h->dim))/2;
    for (i = 0; i < h->dim; i++) {
      for (j = h->dim; j < h->n_basis; j++) {
	r = HamiltonElement(isym, h->basis[i], h->basis[j]);
	h->hamilton[t++] = r;
      }
    }
    for (j = h->dim; j < h->n_basis; j++) {
      r = HamiltonElement(isym, h->basis[j], h->basis[j]);
      h->hamilton[t++] = r;
    }
  }
#ifdef PERFORM_STATISTICS
  stop = clock();
  timing.set_ham += stop-start;
#endif

  return 0;

 ERROR:
#ifdef PERFORM_STATISTICS
  stop = clock();
  timing.set_ham += stop-start;
#endif
  printf("ConstructHamilton Error\n");
  return -1;
}

int ValidBasis(STATE *s, int k, int *kg, int n) {
  int t, m, kb;
  LEVEL *lev;
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
    m = lev->pb;
    sym = GetSymmetry(lev->pj);
    sp = (STATE *) ArrayGet(&(sym->states), m);
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
  double r, delta;
#ifdef PERFORM_STATISTICS
  clock_t start, stop;
  start = clock();
#endif

  j = 0;
  sym = GetSymmetry(isym);
  st = &(sym->states);
  for (t = 0; t < sym->n_states; t++) { 
    s = (STATE *) ArrayGet(st, t);
    if ((i = ValidBasis(s, k, kg, n))) {
      j++;
    } 
  }
  
  if (j == 0) return -1;

  h = &_ham;
  h->pj = isym;

  h->dim = j;
  h->n_basis = j;
  t = j*(j+1)/2;
  h->hsize = t;
  h->msize = h->dim*h->n_basis + h->dim;

  if (h->basis == NULL) {
    h->n_basis0 = h->n_basis;
    h->basis = (int *) malloc(sizeof(int)*(h->n_basis));
  } else if (h->n_basis > h->n_basis0) {
    h->n_basis0 = h->n_basis;
    h->basis = (int *) realloc(h->basis, sizeof(int)*h->n_basis);
  }
  if (!(h->basis)) goto ERROR;

  if (h->hamilton == NULL) {
    h->hsize0 = h->hsize;
    h->hamilton = (double *) malloc(sizeof(double)*h->hsize);
  } else if (h->hsize > h->hsize0) {
    h->hsize0 = h->hsize;
    h->hamilton = (double *) realloc(h->hamilton, sizeof(double)*h->hsize);
  }
  if (!(h->hamilton)) goto ERROR;
      
  j = 0;
  for (t = 0; t < sym->n_states; t++) { 
    s = (STATE *) ArrayGet(st, t);
    if (ValidBasis(s, k, kg, n)) {
      h->basis[j] = t;
      j++;
    }
  }

  for (j = 0; j < h->dim; j++) {
    t = j*(j+1)/2;
    for (i = 0; i <= j; i++) {
      r = HamiltonElementFrozen(isym, h->basis[i], h->basis[j]);
      h->hamilton[i+t] = r;
    }
    for (i = 0; i < j; i++) {
      delta = fabs(h->hamilton[i+t]/h->hamilton[j+t]);
      if (delta < EPS10) h->hamilton[i+t] = 0.0;
    }
  }

#ifdef PERFORM_STATISTICS
  stop = clock();
  timing.set_ham += stop-start;
#endif

  return 0;

 ERROR:
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
  ji1 = lev1->pj;
  jj1 = lev2->pj;
  DecodePJ(ji1, NULL, &ji1);
  DecodePJ(jj1, NULL, &jj1);
  GetJLFromKappa(orbi->kappa, &ji2, &ki2);
  GetJLFromKappa(orbj->kappa, &jj2, &kj2);
  j = si->kstate;
  if (si->kgroup == sj->kgroup) { 
    if (ji2 == jj2 && ki2 == kj2) {
      ResidualPotential(&r0, si->kcfg, sj->kcfg);
      r += r0;
      r += QED1E(si->kcfg, sj->kcfg);
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
  
  if (IsOdd((ji2 + jj1 + j)/2)) a = -a;
  r += a;

  if (nz > 0) {
    free(ang);
  } 

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
  INTERACT_DATUM *idatum;
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

  idatum = NULL;
  n_shells = GetInteract(&idatum, &sbra, &sket, 
			 si->kgroup, sj->kgroup,
			 si->kcfg, sj->kcfg,
			 ki, kj, 0);
  if (n_shells <= 0) return 0.0;
  memcpy(s, idatum->s, sizeof(INTERACT_SHELL)*4);
  bra = idatum->bra;
  phase = idatum->phase;

  x = 0.0;
  if (s[0].index >= 0 && s[3].index >= 0) {
    r = Hamilton2E(n_shells, sbra, sket, s);
    x += r;
  } else if( s[0].index >= 0) {
    r = Hamilton1E(n_shells, sbra, sket, s);
    x += r;
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

  if (isi == isj) {
    x += ci->delta;
  }

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
  if (k1 == k2) r0 += (GetOrbital(k1))->energy;
  r0 += QED1E(k1, k2);

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
    
int DiagnolizeHamilton(void) {
  double *ap;
  double *w;
  double *z, *x, *y, *b, *d, *ep;
  double *mixing = NULL;
  char jobz[] = "V";
  char uplo[] = "U";
  char trans[] = "N";
  int n, m, np;
  int ldz;
  int lwork;
  int liwork;
  int info;
  HAMILTON *h;
  int i, j, t, t0, k, one;
  double d_one, d_zero, a;
#ifdef PERFORM_STATISTICS
  clock_t start, stop;
  start = clock();
#endif

  h = &_ham;
  n = h->dim;
  m = h->n_basis;
  ldz = n;
  t0 = n*(n+1);

  lwork = 1 + 5*n + t0;
  liwork = 3 + 5*n;  
  if (h->work == NULL) {
    h->dim0 = h->dim;
    h->work = (double *) malloc(sizeof(double)*(lwork+2*t0));
    h->iwork = (int *) malloc(sizeof(int)*liwork);
  } else if (h->dim > h->dim0) {
    h->dim0 = h->dim;
    h->work = (double *) realloc(h->work, sizeof(double)*(lwork+2*t0));
    h->iwork = (int *) realloc(h->iwork, sizeof(int)*liwork);
  }

  h->msize = h->dim * h->n_basis + h->dim;  
  if (h->mixing == NULL) {
    h->msize0 = h->msize;
    h->mixing = (double *) malloc(sizeof(double)*h->msize);
  } else if (h->msize > h->msize0) {
    h->msize0 = h->msize;
    h->mixing = (double *) realloc(h->mixing, sizeof(double)*h->msize);
  }
  if (!(h->mixing)) {
    printf("Allocating Mixing Error\n");
    goto ERROR;
  }  
  ap = h->hamilton;
  if (m > n) {
    mixing = h->work + lwork;
  } else {
    mixing = h->mixing;
  }
  
  w = mixing;
  z = mixing + n;
  DSPEVD(jobz, uplo, n, ap, w, z, ldz, h->work, lwork,
	 h->iwork, liwork, &info);
  if (info) {
    printf("dspevd Error: %d\n", info);
    goto ERROR;
  }
  if (m > n) {
    np = m-n;
    b = h->hamilton + t0/2;
    ep = b + n*np;
    y = h->mixing+n;
    one = 1;
    d_one = 1.0;
    d_zero = 0.0;
    for (i = 0; i < n; i++) {
      x = y+n;
      DGEMV(trans, np, n, d_one, b, np, z, one, d_zero, x, one);
      y += m;
      z += n;
    }
    y = h->mixing + 2*n;
    for (j = 0; j < n; j++) {
      t = j*(j+1)/2;
      x = h->mixing + 2*n;
      for (i = 0; i <= j; i++) {
	a = 0.0;
	for (k = 0; k < np; k++) {
	  a += x[k]*y[k]/(w[j] - ep[k]);
	}
	if (i == j) a += w[j];
	h->hamilton[i+t] = a;
	x += m;
      }
      y += m;
    }
    w = h->mixing;
    d = h->work+lwork+t0;
    DSPEVD(jobz, uplo, n, ap, w, d, ldz, h->work, lwork,
	    h->iwork, liwork, &info);
    y = h->mixing+n;
    z = mixing+n;
    for (i = 0; i < n; i++) {
      x = y+n;
      DGEMV(trans, n, n, d_one, z, n, d, one, d_zero, y, one);
      DGEMV(trans, np, n, d_one, b, np, y, one, d_zero, x, one);
      for (j = 0; j < np; j++) {
	x[j] *= 1.0/(w[i]-ep[j]);
      } 
      a = DDOT(np, x, one, x, one);
      a = 1.0/sqrt(1.0+a);
      DSCAL(m, a, y, one);
      y += m;
      d += n;
    }
  }

  return 0;
  
 ERROR: 
#ifdef PERFORM_STATISTICS
  stop = clock();
  timing.diag_ham += stop-start;
#endif

  return -1;
}

int AddToLevels(int ng, int *kg) {
  int i, d, j, k, t, m;
  HAMILTON *h;
  LEVEL lev;
  SYMMETRY *sym;
  STATE *s;
  CONFIG *c;
  int g0, p0;
  double *mix, a;
  
  h = &_ham;
  if (h->basis == NULL ||
      h->mixing == NULL) return -1;
  d = h->dim;
  mix = h->mixing + d;
  j = n_levels;
  sym = GetSymmetry(h->pj);
  for (i = 0; i < d; i++) {
    k = GetPrincipleBasis(mix, d);
    s = (STATE *) ArrayGet(&(sym->states), h->basis[k]);
    if (ng > 0) {
      if (!InGroups(s->kgroup, ng, kg)) {
	mix += h->n_basis;
	continue;
      }
    }
    lev.energy = h->mixing[i];
    lev.pj = h->pj;
    lev.pb = h->basis[k];
    lev.basis = (int *) malloc(sizeof(int)*h->n_basis);
    lev.mixing = (double *) malloc(sizeof(double)*h->n_basis);
    a = fabs(mix_cut * mix[k]);
    for (t = 0, m = 0; t < h->n_basis; t++) {
      if (fabs(mix[t]) < a) continue;
      lev.basis[m] = h->basis[t];
      lev.mixing[m] = mix[t];
      m++;
    }
    lev.n_basis = m;
    if (m < t) {
      lev.basis = (int *) realloc(lev.basis, sizeof(int)*m);
      lev.mixing = (double *) realloc(lev.mixing, sizeof(double)*m);
    }
    SortMixing(0, m, lev.basis, lev.mixing, sym);

    if (s->kgroup >= 0) {
      lev.ngp = 0;
      lev.igp = (int *) malloc(sizeof(int)*(m+1));
      g0 = -1;
      p0 = -1;
      for (t = 0; t < m; t++) {
	s = (STATE *) ArrayGet(&(sym->states), lev.basis[t]);
	c = GetConfig(s);
	if (s->kgroup != g0 || c->ipart != p0) {
	  lev.igp[lev.ngp] = t;
	  lev.ngp++;
	  g0 = s->kgroup;
	  p0 = c->ipart;
	}
      }
      lev.igp[lev.ngp] = m;
      if (lev.ngp < m) {
	lev.igp = (int *) realloc(lev.igp, sizeof(int)*(lev.ngp+1));
      }
    } else {
      lev.ngp = 0;
      lev.igp = NULL;
    }

    if (ArrayAppend(levels, &lev, NULL) == NULL) {
      printf("Not enough memory for levels array\n");
      exit(1);
    }
    j++;
    mix += h->n_basis;
  }
  
  n_levels = j;
  if (i < d-1) return -2;

  return 0;
}

static int CompareBasis(int b1, int b2, 
			double m1, double m2, SYMMETRY *sym) {
  STATE *s1, *s2;
  
  s1 = (STATE *) ArrayGet(&(sym->states), b1);
  s2 = (STATE *) ArrayGet(&(sym->states), b2);
  if (s1->kgroup >= 0 && s2->kgroup >= 0) {
    if (s1->kgroup > s2->kgroup) return 1;
    else if (s1->kgroup < s2->kgroup) return -1;
    else {
      if (s1->kcfg > s2->kcfg) return 1;
      else if (s1->kcfg < s2->kcfg) return -1;
      else {
	if (fabs(m1) > fabs(m2)) return 1;
	else if (fabs(m1) < fabs(m2)) return -1;
	else return 0;
      }
    }
  } else {
    if (fabs(m1) > fabs(m2)) return 1;
    else if (fabs(m1) < fabs(m2)) return -1;
    else return 0;
  }
}  

int SortMixing(int start, int n, 
	       int *basis, double *mix, SYMMETRY *sym) {
  int i, j, i0, j0, t;
  int *b1, *b2, *bp;
  double *m1, *m2, *mp, tmp;

  while (1 < n) {
    i = start;
    j = start + n - 1;
    m1 = mix + i;
    m2 = mix + j;
    b1 = basis + i;
    b2 = basis + j;
    mp = m2;    
    bp = b2;
    
    while (i < j) {
      while (i < j) {
	if (CompareBasis(*b1, *bp, *m1, *mp, sym) < 0) break;
	i++;
	m1 = mix + i;
	b1 = basis + i;
      }
      while (i < j) {
	if (CompareBasis(*bp, *b2, *mp, *m2, sym) < 0) break;
	j--;
	m2 = mix + j;
	b2 = basis + j;
      }
      if (i < j) {
	tmp = *m1;
	*m1 = *m2;
	*m2 = tmp;
	t = *b1;
	*b1 = *b2;
	*b2 = t;
	i++;
	m1 = mix + i;
	b1 = basis + i;
      }
    }
    if (CompareBasis(*b1, *bp, *m1, *mp, sym)) {
      tmp = *m1;
      *m1 = *mp;
      *mp = tmp;
      t = *b1;
      *b1 = *bp;
      *bp = t;
    }

    i0 = i - start;
    j0 = n - i0 - 1;
    if (j0 < i0) {
      if (1 < j0) {
	SortMixing(i+1, j0, basis, mix, sym);
      }
      n = i0;
    } else {
      if (1 < i0) {
	SortMixing(start, i0, basis, mix, sym);
      }
      start = i+1;
      n = j0;
    }
  }
  return 0;
}
  
int AddECorrection(int iref, int ilev, double e, int nmin) {
  ECORRECTION c;

  c.iref = iref;
  c.ilev = ilev;
  c.e = e;
  c.nmin = nmin;
  ArrayAppend(ecorrections, &c, NULL);
  ncorrections += 1;

  return 0;
}

LEVEL *GetLevel(int k) {
  return (LEVEL *) ArrayGet(levels, k);
}

int LevelTotalJ(int k) {
  int i;
  i = GetLevel(k)->pj;
  DecodePJ(i, NULL, &i);
  return i;
}

int GetNumLevels(void) {
  return n_levels;
}

int GetPrincipleBasis(double *mix, int d) {
  int i, k;
  double c = 0.0;
  double fm;

  k = 0;
  for (i = 0; i < d; i++) {
    fm = fabs(mix[i]);
    if (fm > c) {
      c = fm;
      k = i;
    }
  }

  return k;
}

int CompareLevels(LEVEL *lev1, LEVEL *lev2) {
  STATE *s1, *s2;
  SYMMETRY *sym1, *sym2;
  ORBITAL *orb;
  int i1, i2;
  int p1, p2, j1, j2;

  i1 = lev1->pb;
  i2 = lev2->pb;
  sym1 = GetSymmetry(lev1->pj);
  sym2 = GetSymmetry(lev2->pj);
  s1 = (STATE *) ArrayGet(&(sym1->states), i1);
  s2 = (STATE *) ArrayGet(&(sym2->states), i2);
  if (s1->kgroup < 0 && s2->kgroup < 0) {
    orb = GetOrbital(s1->kcfg);
    GetJLFromKappa(orb->kappa, &p1, &j1);
    orb = GetOrbital(s2->kcfg);
    GetJLFromKappa(orb->kappa, &p2, &j2);
    i1 = p1 - p2;
    if (i1) return i1;
    i1 = j1 - j2;
    if (i1) return i1;
    i1 = (s2->kgroup - s1->kgroup);
    if (i1) return i1;
    DecodePJ(lev1->pj, &p1, &j1);
    DecodePJ(lev2->pj, &p2, &j2);
    i1 = p1 - p2;
    if (i1) return i1;
    return (j1 - j2);
  } else {
    if (lev1->energy > lev2->energy) return 1;
    else if (lev1->energy < lev2->energy) return -1;
    else return 0;
  }
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
	if (CompareLevels(lev1, levp) > 0) break;
	i++;
	lev1 = GetLevel(i);
      }
      while (i < j) {
	if (CompareLevels(levp, lev2) > 0) break;
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
  return 0;
}

int GetNumElectrons(int k) {
  LEVEL *lev;
  SYMMETRY *sym;
  STATE *s;
  int nele;
  
  lev = GetLevel(k);
  sym = GetSymmetry(lev->pj);
  s = (STATE *) ArrayGet(&(sym->states), 0);
  nele = ConstructLevelName(NULL, NULL, NULL, NULL, s);
  return nele;
}

int SaveLevels(char *fn, int m, int n) {
  STATE *s;
  SYMMETRY *sym;
  CONFIG *cfg, *cfg1;
  SHELL_STATE *csf, *csf1;
  LEVEL *lev;
  EN_RECORD r;
  EN_HEADER en_hdr;
  F_HEADER fhdr;
  ECORRECTION *ec;
  double e0;
  char name[LEVEL_NAME_LEN];
  char sname[LEVEL_NAME_LEN];
  char nc[LEVEL_NAME_LEN];
  FILE *f;
  int i, k, p, j0;
  int nele, nele0, vnl;
  int si, ms, mst;
#ifdef PERFORM_STATISTICS
  STRUCT_TIMING structt;
  ANGULAR_TIMING angt;
  RECOUPLE_TIMING recouplet;
  RAD_TIMING radt;
#endif
 
  f = NULL;
  nele0 = -1;
  if (n < 0) n = n_levels - m;
  fhdr.type = DB_EN;
  strcpy(fhdr.symbol, GetAtomicSymbol());
  fhdr.atom = GetAtomicNumber();
  f = OpenFile(fn, &fhdr);
  for (k = 0; k < n; k++) {
    i = m + k;
    lev = GetLevel(i);
    si = lev->pb;
    sym = GetSymmetry(lev->pj);
    s = (STATE *) ArrayGet(&(sym->states), si);
    if (ncorrections > 0) {
      for (p = 0; p < ecorrections->dim; p++) {
	ec = (ECORRECTION *) ArrayGet(ecorrections, p);
	if (ec->ilev == i) {
	  if (ec->ilev == ec->iref) {
	    e0 = lev->energy;
	  } else {
	    e0 = GetLevel(ec->iref)->energy;
	  }
	  ec->e = e0 + ec->e - lev->energy;
	  lev->energy += ec->e;
	  ec->s = s;
	  ec->ilev = -(ec->ilev+1);
	  ncorrections -= 1;
	  break;
	}
      }
    }

    if (s->kgroup > 0) {
      cfg = GetConfig(s);
      if ((cfg->shells[0]).n > (cfg->shells[1]).n &&
	  (cfg->shells[0]).nq == 1) {
	for (p = 0; p < ecorrections->dim; p++) {
	  ec = (ECORRECTION *) ArrayGet(ecorrections, p);
	  if (-(i+1) == ec->ilev) break;
	}
	if (p == ecorrections->dim) {
	  csf = cfg->csfs + s->kstate;
	  for (p = 0; p < ecorrections->dim; p++) {
	    ec = (ECORRECTION *) ArrayGet(ecorrections, p);
	    if (ec->ilev >= 0) continue;
	    if ((cfg->shells[0]).n < ec->nmin) continue;
	    cfg1 = GetConfig(ec->s);
	    ms = cfg1->n_shells*sizeof(SHELL);
	    mst = cfg1->n_shells*sizeof(SHELL_STATE);
	    csf1 = cfg1->csfs + ec->s->kstate;
	    if (cfg->n_electrons == cfg1->n_electrons+1 &&
		cfg->n_shells == cfg1->n_shells+1 &&
		memcmp(cfg->shells+1, cfg1->shells, ms) == 0 &&
		memcmp(csf+1, csf1, mst) == 0) {
	      lev->energy += ec->e;
	      break;
	    }
	  }
	}
      }  
    }
 
    DecodePJ(lev->pj, &p, &j0);
    r.ilev = i;
    r.p = p;
    r.j = j0;
    r.energy = lev->energy;

    nele = ConstructLevelName(name, sname, nc, &vnl, s);
    strncpy(r.name, name, LNAME);
    strncpy(r.sname, sname, LSNAME);
    strncpy(r.ncomplex, nc, LNCOMPLEX);
    r.name[LNAME-1] = '\0';
    r.sname[LSNAME-1] = '\0';
    r.ncomplex[LNCOMPLEX-1] = '\0';
    if (r.p == 0) {
      r.p = vnl;
    } else {
      r.p = -vnl;
    }
    if (nele != nele0) {
      if (nele0 >= 0) DeinitFile(f, &fhdr);
      nele0 = nele;
      en_hdr.nele = nele;
      InitFile(f, &fhdr, &en_hdr);
    }
    WriteENRecord(f, &r);
  }

  DeinitFile(f, &fhdr);
  CloseFile(f, &fhdr);

#ifdef PERFORM_STATISTICS
  GetStructTiming(&structt);
  fprintf(perform_log, "AngZMix: %6.1E, AngZFB: %6.1E, AngZxZFB: %6.1E, SetH: %6.1E DiagH: %6.1E\n",
	  ((double) (structt.angz_mix))/CLOCKS_PER_SEC,
	  ((double) (structt.angz_fb))/CLOCKS_PER_SEC,
	  ((double) (structt.angzxz_fb))/CLOCKS_PER_SEC,
	  ((double) (structt.set_ham))/CLOCKS_PER_SEC,
	  ((double) (structt.diag_ham))/CLOCKS_PER_SEC);
  fprintf(perform_log, "AngZS: %6.1E, AngZFBS: %6.1E, AngZxZFBS: %6.1E, AddZ: %6.1E, AddZxZ: %6.1E\n",
	  ((double) (structt.angz_states))/CLOCKS_PER_SEC,
	  ((double) (structt.angzfb_states))/CLOCKS_PER_SEC,
	  ((double) (structt.angzxzfb_states))/CLOCKS_PER_SEC,
	  ((double) (structt.add_angz))/CLOCKS_PER_SEC,
	  ((double) (structt.add_angzxz))/CLOCKS_PER_SEC);

  GetAngularTiming(&angt);
  fprintf(perform_log, "W3J: %6.1E, W6J: %6.1E, W9J: %6.1E\n", 
	  ((double)angt.w3j)/CLOCKS_PER_SEC, 
	  ((double)angt.w6j)/CLOCKS_PER_SEC, 
	  ((double)angt.w9j)/CLOCKS_PER_SEC);
  GetRecoupleTiming(&recouplet);
  fprintf(perform_log, "AngZ: %6.1E, AngZxZ: %6.1E, Interact: %6.1E\n",
	  ((double)recouplet.angz)/CLOCKS_PER_SEC,
	  ((double)recouplet.angzxz)/CLOCKS_PER_SEC,
	  ((double)recouplet.interact)/CLOCKS_PER_SEC);
  GetRadTiming(&radt);
  fprintf(perform_log, "Dirac: %d, %6.1E, 1E: %6.1E, Slater: %6.1E, 2E: %6.1E\n", 
	  GetNumContinua(),
	  ((double)radt.dirac)/CLOCKS_PER_SEC, 
	  ((double)radt.radial_1e)/CLOCKS_PER_SEC,
	  ((double)radt.radial_slater)/CLOCKS_PER_SEC,
	  ((double)radt.radial_2e)/CLOCKS_PER_SEC);
  fprintf(perform_log, "\n");
  fflush(perform_log);
#endif /* PERFORM_STATISTICS */
  
  return 0;
}

int ConstructLevelName(char *name, char *sname, char *nc, 
		       int *vnl, STATE *basis) {
  int n, nq, kl, j;
  int nele, i, len;
  char symbol[20];
  char jsym;
  char ashell[16];
  CONFIG *c;
  SHELL_STATE *s;
  ORBITAL *orb;
  LEVEL *lev;
  SYMMETRY *sym;
  int si;
  int n0, kl0, nq0;

  symbol[0] = '\0';
  if (basis->kgroup < 0) {
    i = basis->kgroup;
    i = -(i + 1);    
    orb = GetOrbital(basis->kcfg);
    GetJLFromKappa(orb->kappa, &j, &kl);
    if (vnl) {
      *vnl = (kl/2) + 100*(orb->n);
    }
    if (name) {
      if (j < kl) jsym = '-';
      else jsym = '+';
      
      kl /= 2;
      SpecSymbol(symbol, kl);
      sprintf(name, "%5d + %d%s%c1(%d)%d ", 
	      i, orb->n, symbol, jsym, j, basis->kstate);
    }
    lev = GetLevel(i);
    si = lev->pb;
    sym = GetSymmetry(lev->pj);
    basis = (STATE *) ArrayGet(&(sym->states), si);
    if (sname || nc) {
      nele = ConstructLevelName(NULL, sname, nc, NULL, basis);
      if (nc) {
	if (nele == 0) {
	  nc[0] = '\0';
	}
	sprintf(ashell, "%1d*1", orb->n);
	strcat(nc, ashell);
      }
    } else {
      nele = ConstructLevelName(NULL, NULL, NULL, NULL, basis);
    }
    return nele+1;
  }

  c = GetConfig(basis);
  nele = c->n_electrons;
  if (!name && !sname && !nc) return nele;

  s = c->csfs + basis->kstate;

  len = 0;
  if (name) name[0] = '\0';
  if (sname) sname[0] = '\0';
  if (nc) nc[0] = '\0';
  n0 = 0;
  kl0= -1;
  nq0 = 0;
  for (i = c->n_shells-1; i >= 0; i--) {
    UnpackShell(c->shells+i, &n, &kl, &j, &nq);
    if (j < kl) jsym = '-';
    else jsym = '+';
    kl = kl/2;
    if (name) {
      if (((nq < j+1) && nq > 0) || (i == 0 && name[0] == '\0')) {
	SpecSymbol(symbol, kl);
	sprintf(ashell, "%1d%s%c%1d(%1d)%1d ", 
		n, symbol, jsym, nq, s[i].shellJ, s[i].totalJ); 
	len += strlen(ashell);
	if (len >= LEVEL_NAME_LEN) return -1;
	strcat(name, ashell);
      }
    }
    if (sname) {
      if (n == n0 && kl == kl0) {
	nq0 += nq;
      } else {
	if (nq0 > 0 && nq0 < 2*(2*kl0+1)) {
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
  
  if (sname && n0 > 0) {
    if ((nq0 > 0 && nq0 < 2*(2*kl0+1)) || sname[0] == '\0') {
      SpecSymbol(symbol, kl0);
      sprintf(ashell, "%1d%s%1d ", n0, symbol, nq0);
      strcat(sname, ashell);
    }
  }

  if (nc) {
    n0 = 0;
    nq0 = 0;
    for (i = c->n_shells-1; i >= 0; i--) {
      UnpackShell(c->shells+i, &n, &kl, &j, &nq);
      if (n == n0) {
	nq0 += nq;
      } else {
	if (nq0 > 0) {
	  sprintf(ashell, "%1d*%1d ", n0, nq0);
	  strcat(nc, ashell);
	}
	n0 = n;
	nq0 = nq;
      }
    }
    if (n0 > 0) {
      sprintf(ashell, "%1d*%1d ", n0, nq0);
      strcat(nc, ashell);
    }
  }

  if (vnl) {
    UnpackShell(c->shells, &n, &kl, &j, &nq);
    *vnl = (kl/2) + 100*n;
  }
      
  return nele;
}
    

int GetBasisTable(char *fn) {
  FILE *f;
  int i, p, j, k, si, nsym;
  char name[LEVEL_NAME_LEN];
  char sname[LEVEL_NAME_LEN];
  ARRAY *st;
  STATE *s;
  LEVEL *lev;
  SYMMETRY *sym;

  f = fopen(fn, "w");
  if (!f) return -1;

  nsym = MAX_SYMMETRIES;
  fprintf(f, "============Basis Table===================\n");
  for (i = 0; i < nsym; i++) {
    sym = GetSymmetry(i);
    DecodePJ(i, &p, &j);
    st = &(sym->states);
    if (sym->n_states <= 0) continue;
    for (k = 0; k < sym->n_states; k++) {
      s = (STATE *) ArrayGet(st, k);
      ConstructLevelName(name, sname, NULL, NULL, s);
      fprintf(f, "%6d   %2d %2d   %3d %3d %3d %3d   %-20s %-20s\n",
	      i, p, j, k, s->kgroup, s->kcfg, s->kstate, sname, name);
    }
    fprintf(f, "\n");
  }

  fprintf(f, "============Mixing Coefficients===================\n");
  for (i = 0; i < n_levels; i++) {
    lev = GetLevel(i);
    sym = GetSymmetry(lev->pj);
    DecodePJ(lev->pj, &p, &j);
    for (k = 0; k < lev->n_basis; k++) {
      si = lev->basis[k];
      s = (STATE *) ArrayGet(&(sym->states), si);
      fprintf(f, "%6d   %2d %2d   %3d %3d %3d %3d   %15.8E\n", 
	      i, p, j, si, s->kgroup, s->kcfg, s->kstate, lev->mixing[k]);
    }
    fprintf(f, "\n");
  }
      
  fclose(f);
  return 0;
}

int AngularZMixStates(ANGZ_DATUM **ad, 
		      int kg1, int kg2, 
		      int kp1, int kp2) {
  int ns, n, p, q, nz, iz;
  int ns1, ns2, *pnz, *pic;
  int nc1, nc2;
  int kc1, kc2;
  int ks1, ks2;
  int n_shells, *k, nkk;
  double *r;
  int phase, im;
  int orb0, orb1;
  CONFIG *c1, *c2;
  PARTITION *part1, *part2;
  INTERACT_DATUM *idatum;
  INTERACT_SHELL s[4];
  SHELL *bra;
  SHELL_STATE *sbra, *sket;
  ANGULAR_ZMIX **a, *ang;
  int index[4];
#ifdef PERFORM_STATISTICS
  clock_t start, stop;
  start = clock();
#endif

  index[0] = kg1;
  index[1] = kg2;
  index[2] = kp1;
  index[3] = kp2;

  *ad = (ANGZ_DATUM *) MultiSet(angz_array, index, NULL, InitAngzDatum);
  ns = (*ad)->ns;
  if (ns < 0) {
#ifdef PERFORM_STATISTICS
    stop = clock();
    timing.angz_states_load += stop-start;
    timing.n_angz_states_load++;
#endif
    return -1;
  }
  
  if (ns > 0) {
#ifdef PERFORM_STATISTICS
    stop = clock();
    timing.angz_states_load += stop-start;
    timing.n_angz_states_load++;
#endif
    return ns;
  }

  part1 = GetPartition(kg1, kp1);
  part2 = GetPartition(kg2, kp2);
  nc1 = part1->icfg2 - part1->icfg1 + 1;
  nc2 = part2->icfg2 - part2->icfg1 + 1;
  (*ad)->ns = part1->n_csfs * part2->n_csfs;
  if ((*ad)->ns == 0) {
    (*ad)->ns = -1;
    return (*ad)->ns;
  }
  ns = (*ad)->ns;
  (*ad)->angz = malloc(sizeof(ANGULAR_ZMIX *)*ns);
  (*ad)->nz = (int *) malloc(sizeof(int)*ns);
  (*ad)->ic = (int *) malloc(sizeof(int)*nc1*nc2);
  
  iz = 0;
  a = (ANGULAR_ZMIX **) (*ad)->angz;
  pnz = (*ad)->nz;
  pic = (*ad)->ic;

  nz = ANGZ_BLOCK;
  for (kc1 = part1->icfg1; kc1 <= part1->icfg2; kc1++) {
    c1 = GetConfigFromGroup(kg1, kc1);
    ns1 = c1->n_csfs * c1->n_shells; 
    for (kc2 = part2->icfg1; kc2 <= part2->icfg2; kc2++) {
      *pic = iz;
      pic++;
      c2 = GetConfigFromGroup(kg2, kc2);  
      if (abs(c1->n_shells - c2->n_shells) > 1) {
	for (ks1 = 0; ks1 < c1->n_csfs; ks1++) {
	  for (ks2 = 0; ks2 < c2->n_csfs; ks2++) {
	    a[iz] = NULL;
	    pnz[iz] = 0;
	    iz++;
	  }
	}
	continue;
      }
      ns2 = c2->n_csfs * c2->n_shells;
      idatum = NULL; 
      for (ks1 = 0; ks1 < ns1; ks1 += c1->n_shells) {  
	for (ks2 = 0; ks2 < ns2; ks2 += c2->n_shells) {
	  n_shells = GetInteract(&idatum, &sbra, &sket, 
				 kg1, kg2, kc1, kc2, 
				 ks1, ks2, 0);
	  n = 0;
	  ang = NULL;
	  if (n_shells <= 0) goto OUT;
	  memcpy(s, idatum->s, sizeof(INTERACT_SHELL)*4);
	  phase = idatum->phase;
	  bra = idatum->bra;
	  if (s[2].index >= 0 && s[0].index >= 0) {
	    free(sbra);
	    free(sket);
	    goto OUT;
	  }

	  ang = malloc(sizeof(ANGULAR_ZMIX)*nz);
	  if (s[0].index >= 0) {
	    nkk = AngularZ(&r, &k, 0, n_shells, sbra, sket, s, s+1);
	    if (nkk > 0) {
	      orb0 = OrbitalIndex(s[0].n, s[0].kappa, 0.0);
	      orb1 = OrbitalIndex(s[1].n, s[1].kappa, 0.0);
	      for (p = 0; p < nkk; p++) {
		if (IsOdd(phase)) r[p] = -r[p];
		im = AddToAngularZMix(&n, &nz, &ang, k[p], orb0, orb1, r[p]);
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
		  if (fabs(r[p]) < EPS10) continue;
		  if (IsOdd(phase)) r[p] = -r[p];
		  im = AddToAngularZMix(&n, &nz, &ang, k[p], orb0, orb1, r[p]);
		}	    
		free(r);
		free(k);
	      }
	    }
	  }
	  PackAngularZMix(&n, &ang, nz);
	  free(sbra);
	  free(sket);
	OUT:
	  a[iz] = ang;
	  pnz[iz] = n;
	  iz++;
	}
      }
    }
  }
#ifdef PERFORM_STATISTICS
  stop = clock();
  timing.angz_states += stop - start;
  timing.n_angz_states++;
#endif
  return (*ad)->ns;
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
    
int AngularZFreeBoundStates(ANGZ_DATUM **ad, 
			    int kg1, int kg2, 
			    int kp1, int kp2) {
  int n_shells;
  int phase;
  int j1, j2, kb;
  int *k, k0, nkk, kmax;
  int jf, jp, tf;
  double *r, r0;
  int ns1, ns2, *pnz, iz, *pic;
  int kc1, kc2, nc1, nc2;
  int ns, ks1, ks2, n;
  PARTITION *part1, *part2;
  INTERACT_DATUM *idatum;
  INTERACT_SHELL s[4];
  SHELL_STATE *sbra, *sket;
  CONFIG *c1, *c2;
  ANGULAR_ZFB *ang, **a;
  
  int index[4];
#ifdef PERFORM_STATISTICS
  clock_t start, stop;
  start = clock();
#endif
  
  index[0] = kg1;
  index[1] = kg2;
  index[2] = kp1;
  index[3] = kp2;
  *ad = (ANGZ_DATUM *) MultiSet(angz_array, index, NULL, InitAngzDatum);
  ns = (*ad)->ns;
  if (ns < 0) {
#ifdef PERFORM_STATISTICS
    stop = clock();
    timing.angzfb_states += stop-start;
#endif
    return -1;
  }
  if (ns > 0) {
#ifdef PERFORM_STATISTICS
    stop = clock();
    timing.angzfb_states += stop-start;
#endif
    return ns;
  }

  part1 = GetPartition(kg1, kp1);
  part2 = GetPartition(kg2, kp2);
  nc1 = part1->icfg2 - part1->icfg1 + 1;
  nc2 = part2->icfg2 - part2->icfg1 + 1;
  (*ad)->ns = part1->n_csfs * part2->n_csfs;
  if ((*ad)->ns == 0) {
    (*ad)->ns = -1;
    return (*ad)->ns;
  }
  ns = (*ad)->ns;
  (*ad)->angz = malloc(sizeof(ANGULAR_ZMIX *)*ns);
  (*ad)->nz = (int *) malloc(sizeof(int)*ns);
  (*ad)->ic = (int *) malloc(sizeof(int)*nc1*nc2);
  
  kmax = GetMaxRank();

  iz = 0;
  a = (ANGULAR_ZFB **) (*ad)->angz;
  pnz = (*ad)->nz;
  pic = (*ad)->ic;

  for (kc1 = part1->icfg1; kc1 <= part1->icfg2; kc1++) {
    c1 = GetConfigFromGroup(kg1, kc1);
    ns1 = c1->n_csfs * c1->n_shells; 
    for (kc2 = part2->icfg1; kc2 <= part2->icfg2; kc2++) {
      *pic = iz;
      pic++;
      c2 = GetConfigFromGroup(kg2, kc2);    
      if (abs(c1->n_shells+1 - c2->n_shells) > 1) {
	for (ks1 = 0; ks1 < c1->n_csfs; ks1++) {
	  for (ks2 = 0; ks2 < c2->n_csfs; ks2++) {
	    a[iz] = NULL;
	    pnz[iz] = 0;
	    iz++;
	  }
	}
	continue;
      }      
      ns2 = c2->n_csfs * c2->n_shells;    
      idatum = NULL;
      for (ks1 = 0; ks1 < ns1; ks1 += c1->n_shells) {
	j1 = (c1->csfs[ks1]).totalJ;
	for (ks2 = 0; ks2 < ns2; ks2 += c2->n_shells) {
	  j2 = (c2->csfs[ks2]).totalJ;
	  n_shells = GetInteract(&idatum, &sbra, &sket,
				 kg1, kg2, kc1, kc2, 
				 ks1, ks2, 1);
	  n = 0;
	  ang = NULL;
	  if (n_shells <= 0) goto END;
	  memcpy(s, idatum->s, sizeof(INTERACT_SHELL)*4);
	  phase = idatum->phase;
	  if (s[0].index < 0 || (s[0].index >= 0 && s[2].index >= 0)) {
	    free(sbra);
	    free(sket);
	    goto END;
	  }
      
	  tf = 0;
	  for (k0 = 0; k0 <= kmax; k0 += 2) {
	    for (jf = abs(k0-s[1].j); jf <= k0+s[1].j; jf += 2) {
	      for (jp = abs(jf-j1); jp <= jf+j1; jp += 2) {
		if (Triangle(jp, j2, k0) && Triangle(j1, j2, s[1].j)) {
		  tf = 1;
		  goto TF;
		}
	      }
	    }
	  }	  
	TF:	  
	  if (tf == 1) {
	    s[0].j = jf;
	    s[0].kl = jf + 1;
	    s[0].kappa = GetKappaFromJL(s[0].j, s[0].kl);
	    sbra[0].shellJ = jf;
	    sbra[0].totalJ = jp; 
	    k = &k0;
	    r = &r0;
	    nkk = AngularZ(&r, &k, 1, n_shells, sbra, sket, s, s+1);
	    if (fabs(*r) < EPS10) goto END;
	    if (IsOdd(phase+(jp+j2-k0)/2)) *r = -(*r);
	    *r /= sqrt(jp+1.0)*W6j(j1, jf, jp, k0, j2, s[1].j);
	    kb = OrbitalIndex(s[1].n, s[1].kappa, 0.0);
	    ang = (ANGULAR_ZFB *) malloc(sizeof(ANGULAR_ZFB));
	    ang->kb = kb;
	    ang->coeff = *r;
	    n = 1;
	  }
	  free(sbra);
	  free(sket);
	END:
	  a[iz] = ang;
	  pnz[iz] = n;
	  iz++;
	}
      }
    }
  }
  
#ifdef PERFORM_STATISTICS
  stop = clock();
  timing.angzfb_states += stop-start;
#endif
  return (*ad)->ns;
}

int AngularZxZMixStates(ANGZ_DATUM **ad, 
			int kg1, int kg2, 
			int kp1, int kp2) {
  return 0;
}

int AngularZxZFreeBoundStates(ANGZ_DATUM **ad,  
			      int kg1, int kg2, 
			      int kp1, int kp2) {
  int n_shells;
  int phase;
  int j1, j2, i, n, nz;
  int jmin, jmax, jf;
  int ns1, ns2, *pnz, iz, *pic;
  int kc1, kc2, nc1, nc2;
  int ns, ks1, ks2;
  PARTITION *part1, *part2;
  INTERACT_SHELL s[4];
  SHELL *bra;
  SHELL_STATE *sbra, *sket;
  INTERACT_DATUM *idatum;
  CONFIG *c1, *c2;
  ANGULAR_ZxZMIX **a, *ang;
  int index[4];

#ifdef PERFORM_STATISTICS
  clock_t start, stop; 
  start = clock();
#endif

  index[0] = kg1;
  index[1] = kg2;
  index[2] = kp1;
  index[3] = kp2;
  *ad = (ANGZ_DATUM *) MultiSet(angzxz_array, index, NULL, InitAngzDatum);

  ns = (*ad)->ns;
  if (ns < 0) { 
#ifdef PERFORM_STATISTICS
    stop = clock();
    timing.angzxzfb_states += stop-start;
#endif
    return -1;
  }
  
  if (ns > 0) {
#ifdef PERFORM_STATISTICS
    stop = clock();
    timing.angzfb_states += stop-start;
#endif
    return ns;
  }

  part1 = GetPartition(kg1, kp1);
  part2 = GetPartition(kg2, kp2);
  nc1 = part1->icfg2 - part1->icfg1 + 1;
  nc2 = part2->icfg2 - part2->icfg1 + 1;
  (*ad)->ns = part1->n_csfs * part2->n_csfs;
  if ((*ad)->ns == 0) {
    (*ad)->ns = -1;
    return (*ad)->ns;
  }
  ns = (*ad)->ns;
  (*ad)->angz = malloc(sizeof(ANGULAR_ZxZMIX *)*ns);
  (*ad)->nz = (int *) malloc(sizeof(int)*ns);
  (*ad)->ic = (int *) malloc(sizeof(int)*nc1*nc2);
  
  iz = 0;
  a = (ANGULAR_ZxZMIX **) (*ad)->angz;
  pnz = (*ad)->nz;
  pic = (*ad)->ic;

  nz = ANGZ_BLOCK;

  for (kc1 = part1->icfg1; kc1 <= part1->icfg2; kc1++) {
    c1 = GetConfigFromGroup(kg1, kc1);
    ns1 = c1->n_csfs * c1->n_shells; 
    for (kc2 = part2->icfg1; kc2 <= part2->icfg2; kc2++) {
      *pic = iz;
      pic++;
      c2 = GetConfigFromGroup(kg2, kc2);    
      if (abs(c1->n_shells+1 - c2->n_shells) > 2) {
	for (ks1 = 0; ks1 < c1->n_csfs; ks1++) {
	  for (ks2 = 0; ks2 < c2->n_csfs; ks2++) {
	    a[iz] = NULL;
	    pnz[iz] = 0;
	    iz++;
	  }
	}
	continue;
      }
      ns2 = c2->n_csfs * c2->n_shells;    
      idatum = NULL;
      for (ks1 = 0; ks1 < ns1; ks1 += c1->n_shells) {
	j1 = (c1->csfs[ks1]).totalJ;
	for (ks2 = 0; ks2 < ns2; ks2 += c2->n_shells) {
	  j2 = (c2->csfs[ks2]).totalJ;
	  n_shells = GetInteract(&idatum, &sbra, &sket,
				 kg1, kg2, kc1, kc2, 
				 ks1, ks2, 1);
	  n = 0;
	  ang = NULL;
	  if (n_shells <= 0) goto END;

	  ang = malloc(sizeof(ANGULAR_ZxZMIX)*nz);
	  phase = idatum->phase;
	  bra = idatum->bra;
	  sbra[0].totalJ = j2;
	  jmin = abs(j2 - j1);
	  jmax = j1 + j2;
	  for (jf = jmin; jf <= jmax; jf += 2) {
	    memcpy(s, idatum->s, sizeof(INTERACT_SHELL)*4);
	    s[0].j = jf;
	    s[0].kl = jf+1;
	    s[0].kappa = GetKappaFromJL(s[0].j, s[0].kl);
	    sbra[0].shellJ = s[0].j;
	    
	    if (s[2].index >= 0) {
	      AddToAngularZxZ(&n, &nz, &ang, n_shells, phase, 
			      sbra, sket, s, 1);
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
		AddToAngularZxZ(&n, &nz, &ang, n_shells, phase, 
				sbra, sket, s, 1);
	      }
	    }
	  }

	  PackAngularZxZMix(&n, &ang, nz);
	  free(sbra);
	  free(sket);

	END:	  
	  a[iz] = ang;
	  pnz[iz] = n;
	  iz++;
	}
      }
    }
  }

#ifdef PERFORM_STATISTICS
  stop = clock();
  timing.angzfb_states += stop-start;
#endif

  return (*ad)->ns;
}


int AngularZFreeBound(ANGULAR_ZFB **ang, int lower, int upper) {
  int i, j, m; 
  int nz, n;
  double r0;
  int kg1, kg2, kc1, kc2, ks1, ks2;
  int igp,jgp, imin, imax, jmin, jmax;
  int ns, isz, kp1, kp2;
  PARTITION *part1, *part2;
  CONFIG *c1, *c2;
  STATE *slow, *sup;
  SYMMETRY *sym1, *sym2;
  LEVEL *lev1, *lev2;
  ANGZ_DATUM *ad;
  ANGULAR_ZFB *ang_sub;
  double mix1, mix2, sqrt_j2;
  int kg, jf, kb, ia, j1, j2;
#ifdef PERFORM_STATISTICS
  clock_t start, stop;
  start = clock();
#endif

  lev1 = GetLevel(lower);
  lev2 = GetLevel(upper);
  sym1 = GetSymmetry(lev1->pj);
  sym2 = GetSymmetry(lev2->pj);
  j1 = lev1->pj;
  j2 = lev2->pj;
  DecodePJ(j1, NULL, &j1);
  DecodePJ(j2, NULL, &j2);
  
  sup = (STATE *) ArrayGet(&(sym2->states), lev2->pb);
  if (sup->kgroup < 0) {
    sqrt_j2 = sqrt(j2 + 1.0);
    n = 0;
    nz = ANGZ_BLOCK;
    (*ang) = malloc(sizeof(ANGULAR_ZFB)*nz);
    for (j = 0; j < lev2->n_basis; j++) {
      mix2 = lev2->mixing[j];
      if (fabs(mix2) < angz_cut) {
	break;
      }
      sup = (STATE *) ArrayGet(&(sym2->states), lev2->basis[j]);
      kg = sup->kgroup;
      kg = -kg-1;
      if (kg == lower) {
	kb = sup->kcfg;
	jf = GetOrbital(kb)->kappa;
	jf = GetJFromKappa(jf);
	r0 = mix2*sqrt_j2;
	if (IsEven((j2+jf-j1)/2)) r0 = -r0;
	ia = AddToAngularZFB(&n, &nz, ang, kb, r0);
      }
    }    
  } else {
    n = 0;
    nz = ANGZ_BLOCK;
    (*ang) = malloc(sizeof(ANGULAR_ZFB)*nz);
    ns = -1;
    for (igp = 0; igp < lev1->ngp; igp++) {
      imin = lev1->igp[igp];
      imax = lev1->igp[igp+1];
      slow = (STATE *) ArrayGet(&(sym1->states), lev1->basis[imin]);
      kg1 = slow->kgroup;
      kc1 = slow->kcfg;
      c1 = GetConfigFromGroup(kg1, kc1);
      kp1 = c1->ipart;
      part1 = GetPartition(kg1, kp1);
      for (jgp = 0; jgp < lev2->ngp; jgp++) {
	jmin = lev2->igp[jgp];
	jmax = lev2->igp[jgp+1];
	sup = (STATE *) ArrayGet(&(sym2->states), lev2->basis[jmin]);
	kg2 = sup->kgroup;
	kc2 = sup->kcfg;
	c2 = GetConfigFromGroup(kg2, kc2);
	kp2 = c2->ipart;
	part2 = GetPartition(kg2, kp2);
	ns = AngularZFreeBoundStates(&ad, kg1, kg2, kp1, kp2);
	if (ns <= 0) continue;
	for (i = imin; i < imax; i++) {
	  mix1 = lev1->mixing[i];
	  if (fabs(mix1) < angz_cut) continue;
	  slow = (STATE *) ArrayGet(&(sym1->states), lev1->basis[i]);
	  kc1 = slow->kcfg;
	  c1 = GetConfigFromGroup(kg1, kc1);
	  ks1 = slow->kstate/c1->n_shells;
	  for (j = jmin; j < jmax; j++) {
	    mix2 = lev2->mixing[j];	
	    r0 = mix1*mix2;
	    if (fabs(r0) < angz_cut) continue;
	    sup = (STATE *) ArrayGet(&(sym2->states), lev2->basis[j]);
	    kc2 = sup->kcfg;
	    c2 = GetConfigFromGroup(kg2, kc2);
	    ks2 = sup->kstate/c2->n_shells;
	    m = (kc1 - part1->icfg1)*(part2->icfg2 - part2->icfg1 + 1);
	    m += kc2 - part2->icfg1;
	    isz = (ad->ic)[m];
	    isz += ks1*c2->n_csfs + ks2;
	    m = (ad->nz)[isz];
	    if (m == 1) {
	      ang_sub = (ad->angz)[isz];
	      kb = ang_sub->kb;
	      r0 *= ang_sub->coeff;
	      ia = AddToAngularZFB(&n, &nz, ang, kb, r0);
	    }
	  }
	}
      }
    }
  }

  PackAngularZFB(&n, ang, nz);

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
  ih = GetLevel(ih)->pj;
  DecodePJ(ih, NULL, &ih);
  return ih;
}

int AngularZMix(ANGULAR_ZMIX **ang, 
		int lower, int upper, 
		int mink, int maxk) {
  int i, j, j1, j2, jb1, jb2;
  int kg1, kg2, kc1, kc2, kp1, kp2;
  int ks1, ks2, isz;
  int jlow, jup, kb1, kb2;
  int nz, n, ns, im;
  int igp, jgp, imin, imax, jmin, jmax;
  double r0;
  int ik, kmin, kmax, m, nmax;
  int nz_sub, nfb;
  PARTITION *part1, *part2;
  STATE *slow, *sup;
  SYMMETRY *sym1, *sym2;
  LEVEL *lev1, *lev2;
  CONFIG *c1, *c2;
  ANGZ_DATUM *ad;
  ANGULAR_ZMIX *ang_sub;
  ANGULAR_ZFB *afb;
  double mix1, mix2, sqrt_j12, a;
  int ignore_ryd = 0;
#ifdef PERFORM_STATISTICS
  clock_t start, stop;
  start = clock();
#endif

  lev1 = GetLevel(lower);
  lev2 = GetLevel(upper);
  sym1 = GetSymmetry(lev1->pj);
  sym2 = GetSymmetry(lev2->pj);
  j1 = lev1->pj;
  j2 = lev2->pj;
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

  slow = (STATE *) ArrayGet(&(sym1->states), lev1->pb);
  sup = (STATE *) ArrayGet(&(sym2->states), lev2->pb);
  kg1 = slow->kgroup;
  kg2 = sup->kgroup;
  kc1 = slow->kcfg;
  kc2 = sup->kcfg;

  if (kg1 < 0) {
    kb1 = slow->kcfg;
    n = GetOrbital(kb1)->n;
    if (rydberg_ignored > 0 && rydberg_ignored < n) ignore_ryd = 1;
    nmax = GetNMax();
    if (n >= nmax && kg2 < 0) {
      kb2 = sup->kcfg;
      m = GetOrbital(kb2)->n;
      if (m == n) {
	ignore_ryd = 1;
      }
    }
  } else if (kg2 < 0) {
    kb2 = sup->kcfg;
    n = GetOrbital(kb2)->n;
    if (rydberg_ignored > 0 && rydberg_ignored < n) ignore_ryd = 1;
  }

  if (kg1 < 0 && kg2 < 0) {
    n = 0;
    nz = ANGZ_BLOCK;
    (*ang) = malloc(sizeof(ANGULAR_ZMIX)*nz);    
    
    for (i = 0; i < lev1->n_basis; i++) {
      mix1 = lev1->mixing[i];
      if (fabs(mix1) < angz_cut) continue;
      slow = (STATE *) ArrayGet(&(sym1->states), lev1->basis[i]);
      jlow = GetBaseJ(slow);
      kg1 = slow->kgroup;
      kg1 = -kg1-1;
      kb1 = slow->kcfg;
      jb1 = GetOrbital(kb1)->kappa;
      jb1 = GetJFromKappa(jb1);
      for (j = 0; j < lev2->n_basis; j++) {
	mix2 = lev2->mixing[j];
	a = mix1*mix2;
	if (fabs(a) < angz_cut) continue;
	sup = (STATE *) ArrayGet(&(sym2->states), lev2->basis[j]);
	jup = GetBaseJ(sup);
	kg2 = sup->kgroup;
	kg2 = -kg2-1;
	kb2 = sup->kcfg;
	jb2 = GetOrbital(kb2)->kappa;
	jb2 = GetJFromKappa(jb2);

	if (ignore_ryd && kb1 != kb2) {
	  continue;
	}

	if (kg1 == kg2) {
	  for (ik = kmin; ik <= kmax; ik += 2) {
	    r0 = W6j(jlow, jb1, j1, ik, j2, jb2);
	    if (fabs(r0) < EPS10) continue;
	    r0 *= a*sqrt_j12;
	    if (IsEven((j1+jb2-jlow-ik)/2+j2)) 
	      r0 = -r0;
	    im = AddToAngularZMix(&n, &nz, ang, ik, kb1, kb2, r0);
	  }
	}
	if (kb1 == kb2){	  
	  nz_sub = AngularZMix(&ang_sub, kg1, kg2, kmin, kmax);
	  if (nz_sub <= 0) {
	    continue;
	  }
	  for (m = 0; m < nz_sub; m++) {
	    r0 = W6j(jlow, jb1, j1, j2, ang_sub[m].k, jup);
	    if (fabs(r0) < EPS10) continue;
	    r0 *= a*ang_sub[m].coeff;
	    r0 *= sqrt_j12;
	    if (IsOdd((jlow+jb1+j2+ang_sub[m].k)/2)) r0 = -r0;
	    im = AddToAngularZMix(&n, &nz, ang, ang_sub[m].k, 
				  ang_sub[m].k0, ang_sub[m].k1, r0);
	  }
	  if (nz_sub > 0) free(ang_sub);
	}
      }
    }
  } else if (kg1 < 0 && !ignore_ryd) {        
    nz = ANGZ_BLOCK;
    n = 0;
    (*ang) = malloc(sizeof(ANGULAR_ZMIX)*nz);
    for (i = 0; i < lev1->n_basis; i++) {
      mix1 = lev1->mixing[i];
      if (fabs(mix1) < angz_cut) continue;
      slow = (STATE *) ArrayGet(&(sym1->states), lev1->basis[i]);
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
	  if (fabs(r0) < EPS10) continue;
	  r0 *= mix1*afb[m].coeff*sqrt(j1+1.0);
	  if (IsOdd((j1+j2-ik)/2)) r0 = -r0;
	  im = AddToAngularZMix(&n, &nz, ang, ik, kb1, afb[m].kb, r0);
	}
      }
      if (nfb > 0) free(afb);
    }
  } else if (kg2 < 0 && !ignore_ryd) {    
    nz = ANGZ_BLOCK;
    n = 0;
    (*ang) = malloc(sizeof(ANGULAR_ZMIX)*nz);
    for (j = 0; j < lev2->n_basis; j++) {
      mix2 = lev2->mixing[j];
      if (fabs(mix2) < angz_cut) continue;
      sup = (STATE *) ArrayGet(&(sym2->states), lev2->basis[j]);
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
	  if (fabs(r0) < EPS10) continue;
	  r0 *= mix2*afb[m].coeff*sqrt(j2+1.0);
	  if (IsOdd((2*j1-ik+jb1-jb2)/2)) r0 = -r0;
	  im = AddToAngularZMix(&n, &nz, ang, ik, afb[m].kb, kb2, r0);
	}
      }
      if (nfb > 0) free(afb);
    }
  } else {
    nz = ANGZ_BLOCK;
    n = 0;
    (*ang) = malloc(sizeof(ANGULAR_ZMIX)*nz);
    for (igp = 0; igp < lev1->ngp; igp++) {
      imin = lev1->igp[igp];
      imax = lev1->igp[igp+1];
      slow = (STATE *) ArrayGet(&(sym1->states), lev1->basis[imin]);
      kg1 = slow->kgroup;
      kc1 = slow->kcfg;
      c1 = GetConfigFromGroup(kg1, kc1);
      kp1 = c1->ipart;
      part1 = GetPartition(kg1, kp1);
      for (jgp = 0; jgp < lev2->ngp; jgp++) {
	jmin = lev2->igp[jgp];
	jmax = lev2->igp[jgp+1];
	sup = (STATE *) ArrayGet(&(sym2->states), lev2->basis[jmin]);
	kg2 = sup->kgroup;
	kc2 = sup->kcfg;
	c2 = GetConfigFromGroup(kg2, kc2);
	kp2 = c2->ipart;
	part2 = GetPartition(kg2, kp2);
	ns = AngularZMixStates(&ad, kg1, kg2, kp1, kp2);
	if (ns <= 0) continue;
	for (i = imin; i < imax; i++) {
	  mix1 = lev1->mixing[i];
	  if (fabs(mix1) < angz_cut) continue;
	  slow = (STATE *) ArrayGet(&(sym1->states), lev1->basis[i]);
	  kc1 = slow->kcfg;
	  c1 = GetConfigFromGroup(kg1, kc1);
	  ks1 = slow->kstate/c1->n_shells;
	  for (j = jmin; j < jmax; j++) {
	    mix2 = lev2->mixing[j];
	    a = mix1*mix2;
	    if (fabs(a) < angz_cut) continue;
	    sup = (STATE *) ArrayGet(&(sym2->states), lev2->basis[j]);
	    kc2 = sup->kcfg;
	    c2 = GetConfigFromGroup(kg2, kc2);
	    ks2 = sup->kstate/c2->n_shells;
	    m = (kc1 - part1->icfg1)*(part2->icfg2 - part2->icfg1 + 1);
	    m += kc2 - part2->icfg1;
	    isz = (ad->ic)[m];
	    isz += ks1*c2->n_csfs + ks2;
	    nz_sub = (ad->nz)[isz];
	    if (nz_sub > 0) {
	      ang_sub = (ad->angz)[isz];
	      for (m = 0; m < nz_sub; m++) {
		r0 = ang_sub[m].coeff*a;
		if (ang_sub[m].k > kmax || ang_sub[m].k < kmin) continue;
		im = AddToAngularZMix(&n, &nz, ang, ang_sub[m].k, 
				      ang_sub[m].k0, ang_sub[m].k1, r0);
	      }
	    }
	  }
	}
      }
    }
  }

  PackAngularZMix(&n, ang, nz);

#ifdef PERFORM_STATISTICS
    stop = clock();
    timing.angz_mix += stop - start;
#endif

  return n;
}

int AngularZxZFreeBound(ANGULAR_ZxZMIX **ang, int lower, int upper) {
  int i, j, j1, j2;
  int nz, n;
  int kg, jf;
  int kg1, kg2, kc1, kc2, ks1, ks2;
  int igp, jgp, imin, imax, jmin, jmax;
  int ns, isz, kp1, kp2;
  double mix1, mix2;
  PARTITION *part1, *part2;
  CONFIG *c1, *c2;
  STATE *slow, *sup;
  SYMMETRY *sym1, *sym2;
  LEVEL *lev1, *lev2;
  ANGZ_DATUM *ad;
  ANGULAR_ZxZMIX *ang_sub;
  ANGULAR_ZMIX *ang_z;
  int kb, jb, jup, orb0, orb1;
  double r, r0, sqrt_j2;
  int nz_sub;
  int im, m;
#ifdef PERFORM_STATISTICS
  clock_t start, stop;
  start = clock();
#endif

  lev1 = GetLevel(lower);
  lev2 = GetLevel(upper);
  sym1 = GetSymmetry(lev1->pj);
  sym2 = GetSymmetry(lev2->pj);
  j1 = lev1->pj;
  j2 = lev2->pj;
  DecodePJ(j1, NULL, &j1);
  DecodePJ(j2, NULL, &j2);
  sqrt_j2 = sqrt(j2+1.0);

  nz = ANGZ_BLOCK;
  n = 0;
  (*ang) = malloc(sizeof(ANGULAR_ZxZMIX)*nz);
  if (!(*ang)) return -1;
    
  sup = (STATE *) ArrayGet(&(sym2->states), lev2->pb);
  if (sup->kgroup < 0) {  
    for (j = 0; j < lev2->n_basis; j++) {
      mix2 = lev2->mixing[j];
      if (fabs(mix2) < angz_cut) continue;
      sup = (STATE *) ArrayGet(&(sym2->states), lev2->basis[j]);
      kb = sup->kcfg;
      jb = GetOrbital(kb)->kappa;
      jb = GetJFromKappa(jb);
      jup = GetBaseJ(sup);
      kg = sup->kgroup;
      kg = -kg-1;
      nz_sub = AngularZMix(&ang_z, lower, kg, -1, -1);
      if (nz_sub <= 0) {
	continue;
      }
      for (i = 0; i < nz_sub; i++) {
	r = ang_z[i].coeff*mix2*sqrt_j2;
	orb0 = ang_z[i].k0;
	orb1 = ang_z[i].k1;
	jmin = abs(j1-j2);
	jmin = Max(jmin, abs(jb-ang_z[i].k));
	jmax = j1+j2;
	jmax = Min(jmax, jb+ang_z[i].k);
	for (jf = jmin; jf <= jmax; jf += 2) {
	  r0 = W6j(j1, jf, j2, jb, jup, ang_z[i].k);
	  if (fabs(r0) < EPS10) continue;
	  if (IsOdd((jup+jf+j2)/2)) r0 = -r0;
	  r0 *= r;
	  im = AddToAngularZxZMix(&n, &nz, ang, ang_z[i].k, 
				  jf, kb, orb0, orb1, r0);
	}
      }
      free(ang_z);
    }
  } else {
    for (igp = 0; igp < lev1->ngp; igp++) {
      imin = lev1->igp[igp];
      imax = lev1->igp[igp+1];
      slow = (STATE *) ArrayGet(&(sym1->states), lev1->basis[imin]);
      kg1 = slow->kgroup;
      kc1 = slow->kcfg;
      c1 = GetConfigFromGroup(kg1, kc1);
      kp1 = c1->ipart;
      part1 = GetPartition(kg1, kp1);
      for (jgp = 0; jgp < lev2->ngp; jgp++) {
	jmin = lev2->igp[jgp];
	jmax = lev2->igp[jgp+1];
	sup = (STATE *) ArrayGet(&(sym2->states), lev2->basis[jmin]);
	kg2 = sup->kgroup;
	kc2 = sup->kcfg;
	c2 = GetConfigFromGroup(kg2, kc2);
	kp2 = c2->ipart;
	part2 = GetPartition(kg2, kp2);
	ns = AngularZxZFreeBoundStates(&ad, kg1, kg2, kp1, kp2);
	if (ns <= 0) continue;
	for (i = imin; i < imax; i++) {
	  mix1 = lev1->mixing[i];
	  if (fabs(mix1) < angz_cut) continue;
	  slow = (STATE *) ArrayGet(&(sym1->states), lev1->basis[i]);
	  kc1 = slow->kcfg;
	  c1 = GetConfigFromGroup(kg1, kc1);
	  ks1 = slow->kstate/c1->n_shells;
	  for (j = jmin; j < jmax; j++) {
	    mix2 = lev2->mixing[j];
	    r = mix1*mix2;
	    if (fabs(r) < angz_cut) continue;
	    sup = (STATE *) ArrayGet(&(sym2->states), lev2->basis[j]); 
	    kc2 = sup->kcfg;
	    c2 = GetConfigFromGroup(kg2, kc2);
	    ks2 = sup->kstate/c2->n_shells;
	    m = (kc1 - part1->icfg1)*(part2->icfg2 - part2->icfg1 + 1);
	    m += kc2 - part2->icfg1;
	    isz = (ad->ic)[m];
	    isz += ks1*c2->n_csfs + ks2;
	    nz_sub = (ad->nz)[isz];
	    if (nz_sub > 0) {
	      ang_sub = (ad->angz)[isz];
	      for (m = 0; m < nz_sub; m++) {
		r0 = ang_sub[m].coeff*r;
		im = AddToAngularZxZMix(&n, &nz, ang, 
					ang_sub[m].k, ang_sub[m].k0,
					ang_sub[m].k1, ang_sub[m].k2,
					ang_sub[m].k3, r0); 	  
	      }
	    }
	  }
	}
      }
    }
  }

  PackAngularZxZMix(&n, ang, nz);

#ifdef PERFORM_STATISTICS
    stop = clock();
    timing.angzxz_fb += stop - start;
#endif
	
  return n;
}
  
int CompareAngularZxZMix(const void *c1, const void *c2) {
  ANGULAR_ZxZMIX *a1, *a2;

  a1 = (ANGULAR_ZxZMIX *) c1;
  a2 = (ANGULAR_ZxZMIX *) c2;
  
  if (a1->k > a2->k) return 1;
  else if (a1->k < a2->k) return -1;
  else {
    if (a1->k0 > a2->k0) return 1;
    else if (a1->k0 < a2->k0) return -1;
    else {
      if (a1->k1 > a2->k1) return 1;
      else if (a1->k1 < a2->k1) return -1;
      else {
	if (a1->k2 > a2->k2) return 1;
	else if (a1->k2 < a2->k2) return -1;
	else {
	  if (a1->k3 > a2->k3) return 1;
	  else if (a1->k3 < a2->k3) return -1;
	  else return 0;
	}
      }
    }
  }
}
  
int CompareAngularZMix(const void *c1, const void *c2) {
  ANGULAR_ZMIX *a1, *a2;

  a1 = (ANGULAR_ZMIX *) c1;
  a2 = (ANGULAR_ZMIX *) c2;
  
  if (a1->k > a2->k) return 1;
  else if (a1->k < a2->k) return -1;
  else {
    if (a1->k0 > a2->k0) return 1;
    else if (a1->k0 < a2->k0) return -1;
    else {
      if (a1->k1 > a2->k1) return 1;
      else if (a1->k1 < a2->k1) return -1;
      else return 0;
    }
  }
}
  
int CompareAngularZFB(const void *c1, const void *c2) {
  ANGULAR_ZFB *a1, *a2;

  a1 = (ANGULAR_ZFB *) c1;
  a2 = (ANGULAR_ZFB *) c2;
  
  if (a1->kb > a2->kb) return 1;
  else if (a1->kb < a2->kb) return -1;
  else return 0;
}

int PackAngularZxZMix(int *n, ANGULAR_ZxZMIX **ang, int nz) {
  int j, m;
  ANGULAR_ZxZMIX *p1, *p2;
  
  m = *n;
  if (*n <= 1) goto OUT;
  if (*n > 2) {
    qsort((void *)(*ang), *n, sizeof(ANGULAR_ZxZMIX), CompareAngularZxZMix);
  }
  m = 1;
  p1 = (*ang);
  j = 1;
  p2 = p1 + 1;
  while (j < *n) {
    if (CompareAngularZxZMix(p1, p2) == 0) {
      p1->coeff += p2->coeff;
    } else {
      p1++;
      m++;
      memcpy(p1, p2, sizeof(ANGULAR_ZxZMIX));
    }
    j++;
    p2++;
  }
  
 OUT:
  if (*n <= 0) {
    if (nz > 0) free(*ang);
  } else {
    if (m < nz) {
      (*ang) = realloc((*ang), m*sizeof(ANGULAR_ZxZMIX));
      *n = m;
    }
  }

  return 0;
}

int PackAngularZMix(int *n, ANGULAR_ZMIX **ang, int nz) {
  int j, m;
  ANGULAR_ZMIX *p1, *p2;

  m = *n;
  if (*n <= 1) goto OUT;
  if (*n > 2) {
    qsort((void *)(*ang), *n, sizeof(ANGULAR_ZMIX), CompareAngularZMix);
  }
  m = 1;
  p1 = (*ang);
  j = 1;
  p2 = p1 + 1;
  while (j < *n) {
    if (CompareAngularZMix(p1, p2) == 0) {
      p1->coeff += p2->coeff;
    } else {
      p1++;
      m++;
      memcpy(p1, p2, sizeof(ANGULAR_ZMIX));
    }
    j++;
    p2++;
  }
  
 OUT:
  if (*n <= 0) {
    if (nz > 0) free(*ang);
  } else {
    if (m < nz) {
      (*ang) = realloc((*ang), m*sizeof(ANGULAR_ZMIX));
      *n = m;
    }
  }

  return 0;
}

int PackAngularZFB(int *n, ANGULAR_ZFB **ang, int nz) {
  int j, m;
  ANGULAR_ZFB *p1, *p2;
  
  m = *n;
  if (*n <= 1) goto OUT;
  if (*n > 2) {
    qsort((void *)(*ang), *n, sizeof(ANGULAR_ZFB), CompareAngularZFB);
  }
  m = 1;
  p1 = (*ang);
  j = 1;
  p2 = p1 + 1;
  while (j < *n) {
    if (CompareAngularZFB(p1, p2) == 0) {
      p1->coeff += p2->coeff;
    } else {
      p1++;
      m++;
      memcpy(p1, p2, sizeof(ANGULAR_ZFB));
    }
    j++;
    p2++;
  }

 OUT:
  if (*n <= 0) {
    if (nz > 0) free(*ang);
  } else {
    if (m < nz) {
      (*ang) = realloc((*ang), m*sizeof(ANGULAR_ZFB));
      *n = m;
    }
  }
  return 0;
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
      if (fabs(r[p]) < EPS10) continue;
      if (IsOdd(phase)) r[p] = -r[p];
      im = AddToAngularZxZMix(n, nz, ang, k[p], 
			      orb0, orb1, kk0, kk1, r[p]);
    }
    free(r);
    free(k);
  }
  
  return 0;
}

int AddToAngularZxZMix(int *n, int *nz, ANGULAR_ZxZMIX **ang, 
		       int k, int k0, int k1, int k2, int k3, double r) {
  int im;
#ifdef PERFORM_STATISTICS
  clock_t start, stop;
  start = clock();
#endif

  im = *n;
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
#ifdef PERFORM_STATISTICS
  stop = clock();
  timing.add_angzxz += stop -start;
#endif

  return 0;
}

int AddToAngularZMix(int *n, int *nz, ANGULAR_ZMIX **ang,
		     int k, int k0, int k1, double coeff) {
  int im;
#ifdef PERFORM_STATISTICS
  clock_t start, stop;
  start = clock();
#endif
  
  im = *n;
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
  
#ifdef PERFORM_STATISTICS
  stop = clock();
  timing.add_angz += stop -start;
#endif

  return 0;
}

int AddToAngularZFB(int *n, int *nz, ANGULAR_ZFB **ang,
		    int kb, double coeff) {
  int im;

  im = *n;
  (*n)++;
  if (*n > *nz) {
    *nz += ANGZ_BLOCK;
    *ang = realloc((*ang), (*nz)*sizeof(ANGULAR_ZFB));
    if (!(*ang)) {
      printf("Cannot enlarge AngularZFB array\n");
      return -1;
    }
  }
  (*ang)[im].kb = kb;
  (*ang)[im].coeff = coeff;

  return 0;
}

void FreeAngZDatum(void *p) {
  ANGZ_DATUM *ap;
  int i;

  ap = (ANGZ_DATUM *) p;
  for (i = 0; i < ap->ns; i++) {
    if (ap->nz[i] > 0) free(ap->angz[i]);
  }
  if (ap->ns > 0) {
    free(ap->angz);
    free(ap->nz);
    free(ap->ic);
  }
}

int FreeAngZArray(int g, MULTI *ma) {  
  ARRAY *a, *b;  
  int ndim;
  
  a = ma->array;
  if (a == NULL) return 0;
  ndim = ma->ndim;
  if (g < 0) {
    MultiFreeData(a, ndim, FreeAngZDatum);
  } else {
    b = (ARRAY *) ArrayGet(a, g);
    if (b) {
      MultiFreeData(b, ndim - 1, FreeAngZDatum);
    }
  }
  return 0;
}

int FreeAngZ(int g, int which_array) {

  if (which_array == 0) {
    FreeAngZArray(g, angz_array);
  } else if (which_array > 0) {
    FreeAngZArray(g, angzxz_array);
  } else {
    FreeAngZArray(g, angz_array);
    FreeAngZArray(g, angzxz_array);
  }
  return 0;
}

void FreeLevelData(void *p) {
  LEVEL *lev;
  lev = (LEVEL *) p;
  if (lev->n_basis > 0) {
    free(lev->basis);
    free(lev->mixing);
    free(lev->igp);
    lev->n_basis = 0;
  }
}
   
int ClearLevelTable(void) {
  CONFIG_GROUP *g;
  CONFIG *cfg;
  int ng, i, k;

  n_levels = 0;
  ArrayFree(levels, FreeLevelData);

  ng = GetNumGroups();
  for (k = 0; k < ng; k++) {
    g = GetGroup(k);
    for (i = 0; i < g->n_cfgs; i++) {
      cfg = (CONFIG *) ArrayGet(&(g->cfg_list), i);
      cfg->energy = 0.0;
      cfg->delta = 0.0;
    }
  }
  
  ncorrections = 0;
  ArrayFree(ecorrections, NULL);

  return 0;
}

int InitStructure(void) {
  int i, ndim = 4;
  int blocks[4];

  n_levels = 0;
  levels = malloc(sizeof(ARRAY));
  if (!levels) return -1;
  ArrayInit(levels, sizeof(LEVEL), LEVELS_BLOCK);

  for (i = 0; i < ndim; i++) blocks[i] = MULTI_BLOCK4;
  angz_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(angz_array, sizeof(ANGZ_DATUM), ndim, blocks);

  angzxz_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(angzxz_array, sizeof(ANGZ_DATUM), ndim, blocks);
  
  ecorrections = malloc(sizeof(ARRAY));
  ArrayInit(ecorrections, sizeof(ECORRECTION), 512);
  ncorrections = 0;

  return 0;
}

int ReinitStructure(int m) {

  if (m < 0) {
    return 0;
  } else if (m == 0) {
    FreeAngZ(-1, -1);
    ClearLevelTable();
  } else if (m == 1) {
    FreeAngZ(-1, -1);
  } else {
    ClearLevelTable();
  }

  return 0;
}

