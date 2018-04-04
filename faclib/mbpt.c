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

#include "mbpt.h"
#include "cf77.h"
#include "mpiutil.h"

static char *rcsid="$Id$";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

static int mbpt_extra = 0;
static int mbpt_reinit_ncps = 0;
static double mbpt_reinit_mem = 0;
static int mbpt_nlev = 0;
static int *mbpt_ilev = NULL;
static double mbpt_mcut = EPS4;
static int mbpt_n3 = 0;
static int mbpt_3rd = 0;
static int mbpt_nsplit = 0;
static int mbpt_ne = 0;
static int *mbpt_se = NULL;
static int *mbpt_de = NULL;
static CONFIG **mbpt_cs = NULL;
static CONFIG mbpt_cfg;
static int *mbpt_bas0, *mbpt_bas0s, *mbpt_bas0d, *mbpt_bas1;
static IDXARY mbpt_ibas0, mbpt_ibas1;
static int **mbpt_rij = NULL;
static struct {
  int nj;
  int *jp;
  IDXARY ibs;
} mbptjp = {0, NULL};

#pragma omp threadprivate(mbpt_cs, mbpt_cfg, mbpt_bas0, mbpt_bas0s, mbpt_bas0d, mbpt_bas1, mbpt_ibas0, mbpt_ibas1, mbptjp)
  
static TR_OPT mbpt_tr;
  
void InitMBPT(void) {
  mbpt_tr.mktr = 0;
  mbpt_tr.naw = 0;
  mbpt_tr.awgrid = NULL;
  mbpt_tr.nlow = 0;
  mbpt_tr.nup = 0;
  mbpt_tr.low = NULL;
  mbpt_tr.up = NULL;
  mbpt_ne = 10;
  int n = mbpt_ne*mbpt_ne;
  int i;
  mbpt_se = (int *) malloc(sizeof(int)*n);
  mbpt_de = (int *) malloc(sizeof(int)*n);
  for (i = 0; i < n; i++) {
    mbpt_se[i] = -1;
    mbpt_de[i] = -1;
  }
}

int IdxSD(int n, int k) {
  int kl = GetLFromKappa(k);
  return (n-1)*(n-1) + (kl-1) + (k<0);
}

void TransitionMBPT(int mk, int n) {  
  if (mk < 0) mk = 0;
  if (n < 0) {
    if (mk > 0) n = 3;
    else n = 0;
  }
  if (n > 0 && n != mbpt_tr.naw) {
    if (mbpt_tr.naw > 0) free(mbpt_tr.awgrid);
    mbpt_tr.awgrid = malloc(sizeof(double)*n);
  }
  mbpt_tr.mktr = mk;
  mbpt_tr.naw = n;
}

int GetAWGridMBPT(double **awgrid) {
  if (awgrid) *awgrid = mbpt_tr.awgrid;
  return mbpt_tr.naw;
}

void SetAWGridMBPT(double emin, double emax) {
  int i;
  double *a;

  SetAWGrid(mbpt_tr.naw, emin, emax);
  GetAWGrid(&a);
  for (i = 0; i < mbpt_tr.naw; i++) {
    mbpt_tr.awgrid[i] = a[i];
  }
}

void TRTableMBPT(char *fn, int nlow, int *low, int nup, int *up) {
  if (nlow <= 0 || nup <= 0) {
    printf("No lower or upper levels %d %d\n", nlow, nup);
    return;
  }
  sprintf(mbpt_tr.tfn, "%s", fn);  
  mbpt_tr.nlow = nlow;
  mbpt_tr.nup = nup;
  mbpt_tr.low = malloc(sizeof(int)*nlow);
  mbpt_tr.up = malloc(sizeof(int)*nup);
  memcpy(mbpt_tr.low, low, sizeof(int)*nlow);
  memcpy(mbpt_tr.up, up, sizeof(int)*nup);
}

void SetOptMBPT(int i3rd, int n3, double c) {
  mbpt_3rd = i3rd;
  mbpt_n3 = n3;
  mbpt_mcut = c;
}

void SetExtraMBPT(int m) {
  int am = abs(m);
  mbpt_extra = am%10;
  if (am < 10) {
    mbpt_reinit_ncps = 0;
    mbpt_reinit_mem = 0;
  } else if (m < 0) {
    mbpt_reinit_mem = 0;
    mbpt_reinit_ncps = am/10;
  } else {
    mbpt_reinit_ncps = 0;
    mbpt_reinit_mem = 1e9*(am/10);
  }
}

void SetSymMBPT(int nlev, int *ilev) {
  if (mbpt_nlev > 0) free(mbpt_ilev);  
  if (nlev > 0 && ilev[0] < 0) nlev = 0;
  mbpt_nlev = nlev;
  if (nlev > 0) {
    mbpt_ilev = malloc(sizeof(int)*nlev);
    memcpy(mbpt_ilev, ilev, sizeof(int)*nlev);
    qsort(mbpt_ilev, nlev, sizeof(int), CompareInt);
  }
}

void SetExcMBPT(int nd, int ns, char *s) {
  int i, j, k, n, nc;
  char *p;
  char s0[512];
  CONFIG *cfg;
  
  strncpy(s0, s, 511);
  s = s0;
  while (*s == ' ') s++;
  if (*s == '\0') {
    n = mbpt_ne*mbpt_ne;    
    for (k = 0; k < n; k++) {
      mbpt_se[i] = nd;
      mbpt_de[k] = ns;
    }
    return;
  }
  n = StrSplit(s, ' ');
  p = s;
  for (i = 0; i < n; i++) {
    while (*p == ' ') p++;
    nc = GetConfigFromString(&cfg, p);
    for (j = 0; j < nc; j++) {
      if (cfg[j].n_shells != 1) {
	printf("incorrect mbpt excitation limit spec: %d %s\n", i, p);
	continue;
      }
      k = IdxSD((cfg[j].shells)[0].n, (cfg[j].shells)[0].kappa);
      if (k >= 0) {
	mbpt_se[k] = ns;
	mbpt_de[k] = nd;
      }
    }
    if (nc > 0) free(cfg);
    while (*p) p++;
    p++;
  }
}

int AddMBPTConfig(int kgp, CORR_CONFIG *ccp) {
  CONFIG *c1, cp;

  c1 = ccp->c;

  if (ccp->kp == 0) {
    cp.n_shells = c1->n_shells;
    cp.shells = malloc(sizeof(SHELL)*cp.n_shells);
    memcpy(cp.shells, c1->shells, sizeof(SHELL)*c1->n_shells);
    Couple(&cp);
    AddConfigToList(kgp, &cp);
    return 0;
  }
  if (ccp->kq == 0) {
    if (c1->shells[0].nq != 0) {
      cp.n_shells = c1->n_shells + 1;
      cp.shells = malloc(sizeof(SHELL)*cp.n_shells);
      memcpy(cp.shells+1, c1->shells, sizeof(SHELL)*c1->n_shells);
    } else {
      cp.n_shells = 1;
      cp.shells = malloc(sizeof(SHELL)*cp.n_shells);
    }
    cp.shells[0].nq = 1;
    cp.shells[0].n = ccp->np;
    cp.shells[0].kappa = ccp->kp;
    Couple(&cp);
    AddConfigToList(kgp, &cp);
  } else {
    if (ccp->np == ccp->nq && ccp->kp == ccp->kq) {
      if (c1->shells[0].nq != 0) {
	cp.n_shells = c1->n_shells + 1;
	cp.shells = malloc(sizeof(SHELL)*cp.n_shells);
	memcpy(cp.shells+1, c1->shells, 
	       sizeof(SHELL)*c1->n_shells);
      } else {
	cp.n_shells = 1;
	cp.shells = malloc(sizeof(SHELL)*cp.n_shells);
      }
      cp.shells[0].nq = 2;
      cp.shells[0].n = ccp->np;
      cp.shells[0].kappa = ccp->kp;
    } else {
      if (c1->shells[0].nq != 0) {
	cp.n_shells = c1->n_shells + 2;
	cp.shells = malloc(sizeof(SHELL)*cp.n_shells);
	memcpy(cp.shells+2, c1->shells, 
	       sizeof(SHELL)*c1->n_shells);
      } else {
	cp.n_shells = 2;
	cp.shells = malloc(sizeof(SHELL)*cp.n_shells);
      }
      cp.shells[1].nq = 1;
      cp.shells[1].n = ccp->np;
      cp.shells[1].kappa = ccp->kp;
      cp.shells[0].nq = 1;
      cp.shells[0].n = ccp->nq;
      cp.shells[0].kappa = ccp->kq;
    }
    Couple(&cp);
    AddConfigToList(kgp, &cp);
  }
  
  return 0;
}    

void ConfigChangeNE(int *m, CONFIG *c, int k, int *ik, int n, int *em) {
  int na, i, j, p;
  int nn, jj, ka;
  SHELL *s;

  na = abs(n);
  for (j = 0; j < na; j++) {
    if (n < 0) {
      ik[em[j]]--;
    } else {
      ik[em[j]]++;
    }
  }
  for (j = 0; j < na; j++) {
    if (ik[em[j]] < 0) {
      goto OUT;
    }
    IntToShell(em[j], &nn, &ka);
    jj = GetJFromKappa(ka);
    if (ik[em[j]] > jj+1) goto OUT;
  }
  p = 0;
  for (i = 0; i < k; i++) {
    if (ik[i] >= 0) p++;
  }
  s = malloc(sizeof(SHELL)*p);
  j = 0;
  for (i = k-1; i >= 0; i--) {
    if (ik[i] >= 0) {
      IntToShell(i, &nn, &ka);
      s[j].n = nn;
      s[j].kappa = ka;
      s[j].nq = ik[i];
      j++;
    }
  }
  nn = 0;
  for (i = 0; i < *m; i++) {
    if (p == c[i].n_shells) {
      ka = memcmp(c[i].shells, s, sizeof(SHELL)*p);
      if (ka == 0) {
	nn = 1;
	break;
      }
    }
  }
  if (nn == 1) {
    free(s);
  } else {
    c[*m].shells = s;
    c[*m].n_shells = p;
    (*m)++;
  }

 OUT:  
  for (j = 0; j < na; j++) {
    if (n < 0) ik[em[j]]++;
    else ik[em[j]]--;
  }
}
	
void RemoveEmpty(CONFIG *c) {
  int i, n;
  SHELL *s;
  
  n = 0;
  for (i = 0; i < c->n_shells; i++) {
    if (c->shells[i].nq > 0) {
      n++;
    }
  }
  if (n == 0) {
    n = 1;
    s = malloc(sizeof(SHELL));
    PackShell(s, 1, 0, 1, 0);
    goto END;
  }
  s = malloc(sizeof(SHELL)*n);
  n = 0;
  for (i = 0; i < c->n_shells; i++) {
    if (c->shells[i].nq > 0) {
      memcpy(s+n, c->shells+i, sizeof(SHELL));
      n++;
    }
  }
 END:
  free(c->shells);
  c->shells = s;
  c->n_shells = n;
}

void BaseConfig(int n, int *kg, int n3, int *n3g, int n4, int *n4g,
		int n1, int *n1g, int n2, int *n2g, int *nbc, CONFIG **bc, 
		int *nbc1, CONFIG **bc1, int *nb, int **bk, FILE *f) {
  int nmax, nsm, ncs, i, j, km, k, t, m, p, jp, nq;
  int mcs, m1e, m2e, mpe, mcs1, *ik, em[2], t1, t2;
  CONFIG *c, *c1, **c0;
  CONFIG_GROUP *g;
  char sc[2000];

  nmax = 0;
  nsm = 0;
  ncs = 0;
  for (i = 0; i < n; i++) {
    g = GetGroup(kg[i]);
    ncs += g->n_cfgs;
    for (j = 0; j < g->n_cfgs; j++) {
      c = GetConfigFromGroup(kg[i], j);
      if (nmax < c->shells[0].n) nmax = c->shells[0].n;
      if (nsm < c->n_shells) nsm = c->n_shells;
    }
  }

  k = nmax*nmax;
  ik = malloc(sizeof(int)*k);
  for (i = 0; i < k; i++) {
    ik[i] = -1;
  }
  c0 = malloc(sizeof(CONFIG *)*ncs);
  m = 0;
  for (i = 0; i < n; i++) {
    g = GetGroup(kg[i]);
    nq = g->n_electrons;
    for (j = 0; j < g->n_cfgs; j++) {
      c = GetConfigFromGroup(kg[i], j);
      for (km = 0; km < c->n_shells; km++) {
	t = ShellToInt(c->shells[km].n, c->shells[km].kappa);
	ik[t] = 0;
      }
      c0[m++] = c;
    }
  }

  *nb = 0;
  for (j = 0; j < k; j++) {
    if (ik[j] == 0) (*nb)++;
  }
  (*bk) = malloc(sizeof(int)*(*nb));
  i = 0;
  for (j = 0; j < k; j++) {
    if (ik[j] == 0) {
      IntToShell(j, &p, &m);
      t = OrbitalIndex(p, m, 0.0);
      (*bk)[i++] = t;
    }
  }
  qsort(*bk, *nb, sizeof(int), CompareInt);

  mcs = ncs*nsm*(1+nsm+nsm*nsm);
  (*bc) = malloc(sizeof(CONFIG)*mcs);
  c = *bc;

  m = 0;
  for (i = 0; i < ncs; i++) {    
    for (j = 0; j < k; j++) {
      if (ik[j] >= 0) ik[j] = 0;
    }
    for (j = 0; j < c0[i]->n_shells; j++) {
      t = ShellToInt(c0[i]->shells[j].n, c0[i]->shells[j].kappa);
      ik[t] = c0[i]->shells[j].nq;
    }
    for (j = 0; j < c0[i]->n_shells; j++) {
      t = IBisect(c0[i]->shells[j].n, n3, n3g);
      if (t < 0) continue;
      em[0] = ShellToInt(c0[i]->shells[j].n, c0[i]->shells[j].kappa);      
      ConfigChangeNE(&m, c, k, ik, -1, em);      
    }
  }
  m1e = m;
  c = c + m;
  m = 0;
  for (i = 0; i < ncs; i++) {
    for (j = 0; j < k; j++) {
      if (ik[j] >= 0) ik[j] = 0;
    }
    for (j = 0; j < c0[i]->n_shells; j++) {
      t = ShellToInt(c0[i]->shells[j].n, c0[i]->shells[j].kappa);
      ik[t] = c0[i]->shells[j].nq;
    }
    for (j = 0; j < c0[i]->n_shells; j++) {
      em[0] = ShellToInt(c0[i]->shells[j].n, c0[i]->shells[j].kappa);
      for (jp = 0; jp <= j; jp++) {
	t1 = IBisect(c0[i]->shells[j].n, n3, n3g);
	t2 = IBisect(c0[i]->shells[jp].n, n4, n4g);
	if (t1 < 0 || t2 < 0) {
	  t2 = IBisect(c0[i]->shells[j].n, n3, n3g);
	  t1 = IBisect(c0[i]->shells[jp].n, n4, n4g);
	  if (t1 < 0 || t2 < 0) continue;
	}
	em[1] = ShellToInt(c0[i]->shells[jp].n, c0[i]->shells[jp].kappa);
	ConfigChangeNE(&m, c, k, ik, -2, em);
      }
    }
  }
  m2e = m;
  m = m1e+m2e; 
  if (n2 > 1 || (n2 == 1 &&  n2g[0] > 0)) {
    c1 = *bc + m1e;
    c = *bc;
    for (i = 0; i < m2e; i++) {
      for (j = 0; j < k; j++) {
	if (ik[j] >= 0) ik[j] = 0;
      }
      for (j = 0; j < c1[i].n_shells; j++) {
	t = ShellToInt(c1[i].shells[j].n, c1[i].shells[j].kappa);
	ik[t] = c1[i].shells[j].nq;
      }
      for (j = 0; j < c1[i].n_shells; j++) {
	t = IBisect(c1[i].shells[j].n, n2, n2g);
	if (t < 0) continue;
	em[0] = ShellToInt(c1[i].shells[j].n, c1[i].shells[j].kappa);
	ConfigChangeNE(&m, c, k, ik, 1, em);
      }
    }
  }
  mpe = m - m1e - m2e;
  
  mcs1 = (m1e+mpe)*nsm + m2e*nsm*nsm;
  *bc1 = malloc(sizeof(CONFIG)*mcs1);
  c1 = *bc1;
  m = 0;
  c = *bc;
  t = 0;
  if (n2 > 0) {
    t = IBisect(0, n2, n2g);
  }
  if (t >= 0) {
    for (i = 0; i < m1e; i++) {
      for (j = 0; j < k; j++) {
	if (ik[j] >= 0) ik[j] = 0;
      }
      for (j = 0; j < c[i].n_shells; j++) {
	t = ShellToInt(c[i].shells[j].n, c[i].shells[j].kappa);
	ik[t] = c[i].shells[j].nq;
      }
      for (j = 0; j < c[i].n_shells; j++) {
	t = IBisect(c[i].shells[j].n, n1, n1g);
	if (t < 0) continue;
	em[0] = ShellToInt(c[i].shells[j].n, c[i].shells[j].kappa);
	ConfigChangeNE(&m, c1, k, ik, 1, em);
      }
    }
  }
  c += m1e;
  /*
  for (i = 0; i < m2e; i++) {
    for (j = 0; j < k; j++) {
      if (ik[j] >= 0) ik[j] = 0;
    }
    for (j = 0; j < c[i].n_shells; j++) {
      t = ShellToInt(c[i].shells[j].n, c[i].shells[j].kappa);
      ik[t] = c[i].shells[j].nq;
    }
    for (j = 0; j < c[i].n_shells; j++) {
      em[0] = ShellToInt(c[i].shells[j].n, c[i].shells[j].kappa);
      for (jp = 0; jp <= j; jp++) {
	t1 = IBisect(c[i].shells[j].n, n1, n1g);
	t2 = IBisect(c[i].shells[jp].n, n2, n2g);
	if (t1 < 0 || t2 < 0) {
	  t2 = IBisect(c[i].shells[j].n, n1, n2g);
	  t1 = IBisect(c[i].shells[jp].n, n1, n2g);
	  if (t1 < 0 || t2 < 0) continue;
	}
	em[1] = ShellToInt(c[i].shells[jp].n, c[i].shells[jp].kappa);
	ConfigChangeNE(&m, c1, k, ik, 2, em);
      }
    }
  }
  */
  c += m2e;
  for (i = 0; i < mpe; i++) {
    for (j = 0; j < k; j++) {
      if (ik[j] >= 0) ik[j] = 0;
    }
    for (j = 0; j < c[i].n_shells; j++) {
      t = ShellToInt(c[i].shells[j].n, c[i].shells[j].kappa);
      ik[t] = c[i].shells[j].nq;
    }
    for (j = 0; j < c[i].n_shells; j++) {
      t = IBisect(c[i].shells[j].n, n1, n1g);
      if (t < 0) continue;
      em[0] = ShellToInt(c[i].shells[j].n, c[i].shells[j].kappa);
      ConfigChangeNE(&m, c1, k, ik, 1, em);
    }
  }

  *nbc = m1e + m2e + mpe;  
  c = *bc;
  for (i = 0; i < *nbc; i++) {
    RemoveEmpty(c+i);
    if (i < m1e || i >= m1e+m2e) c[i].n_electrons = nq-1;
    else c[i].n_electrons = nq-2;
  }

  c1 = *bc1;
  *nbc1 = m;
  *bc1 = malloc(sizeof(CONFIG)*m);
  c = *bc1;
  k = 0;
  for (i = 0; i < *nbc1; i++) {
    RemoveEmpty(c1+i);
    p = 0;
    for (j = 0; j < ncs; j++) {      
      if (c1[i].n_shells == c0[j]->n_shells) {
	t = memcmp(c1[i].shells, c0[j]->shells, sizeof(SHELL)*c0[j]->n_shells);
	if (t == 0) {
	  p = 1;
	  break;
	}
      }
    }
    if (p == 0) {
      c[k].n_shells = c1[i].n_shells;
      c[k].shells = c1[i].shells;
      c[k].n_electrons = nq;
      k++;
    } else {
      free(c1[i].shells);
    }
  }
  *nbc1 = k;
  free(ik);  
  free(c1);

  if (f) {
    for (i = 0; i < ncs; i++) {
      ConstructConfigName(sc, 2000, c0[i]);
      fprintf(f, "0E: %4d   %s\n", i, sc);
    }
    fprintf(f, "\n");
    c = *bc;
    for (i = 0; i < m1e; i++) {
      ConstructConfigName(sc, 2000, c+i);
      fprintf(f, "1E: %4d   %s\n", i, sc);
    }
    fprintf(f, "\n");
    c += m1e;
    for (i = 0; i < m2e; i++) {
      ConstructConfigName(sc, 2000, c+i);
      fprintf(f, "2E: %4d   %s\n", i, sc);
    }
    fprintf(f, "\n");
    c += m2e;
    for (i = 0; i < mpe; i++) {
      ConstructConfigName(sc, 2000, c+i);
      fprintf(f, "2P: %4d   %s\n", i, sc);
    }
    fprintf(f, "\n");    
    c = *bc1;
    for (i = 0; i < *nbc1; i++) {
      ConstructConfigName(sc, 2000, c+i);
      fprintf(f, "EX: %4d   %s\n", i, sc);
    }
    fprintf(f, "\n");
    fflush(f);
  }
  free(c0);
}

int ConstructNGrid(int n, int **g) {
  int i, *p, *q, n0, n1, nt, m, di, dm;

  if (n <= 0) return n;
  p = *g;
  if (p[0] >= 0) {
    if (n > 1) {
      qsort(p, n, sizeof(int), CompareInt);
    }
    return n;
  }
  
  if (n != 5) return -1;
  n1 = -p[0];
  n0 = p[1];
  nt = p[2];
  m = p[3];
  dm = p[4];

  n = (n1-n0+1) + nt;
  q = malloc(sizeof(int)*n);
  for (i = n0; i < n1; i++) {
    q[i] = i;
  }
  di = 2;
  for (; i < n; i++) {
    q[i] = q[i-1] + di;
    if (di < dm) {
      if (m == 0) {
	di *= 2;
      } else {
	di += m;
      }
      if (di >= dm) di = dm;
    }
  }

  free(*g);
  *g = q;
  return n;
}
  
int StrongInteractConfig(int ncc, double *e1, CORR_CONFIG *ccp, 
			 int k1, double a, double de) {
  int k0, m, q;
  double b, d;

  GetSymmetrySet(&k0, &m);
  if (ccp->np > 0 && ccp->kp != 0) {
    a += GetOrbital(OrbitalIndex(ccp->np, ccp->kp, 0))->energy;
    k1 += GetLFromKappa(ccp->kp)/2;
  }
  if (ccp->nq > 0 && ccp->kq != 0) {
    a += GetOrbital(OrbitalIndex(ccp->nq, ccp->kq, 0))->energy;
    k1 += GetLFromKappa(ccp->kq)/2;
  }
  k1 = IsOdd(k1);
  if (k0 >= 0 && k0 != k1) return 0;

  if (de <= 0) return 1;  
  b = 1e30;
  for (q = 0; q < ncc; q++) {
    d = fabs(e1[q] - a);
    if (b > d) b = d;
  }

  if (b < de) {
    return 1;
  }

  return 0;
}
  
int StructureMBPT0(char *fn, double de, double ccut, int n, int *s0, int kmax,
		   int n1, int *nm, int n2, int *nmp, 
		   int n3, int *n3g, int n4, int *n4g, char *gn) {
  CONFIG_GROUP *g1, *g0;
  CONFIG *c1, *bc, *bc1;
  SYMMETRY *sym;
  STATE *st;
  LEVEL *lev;
  HAMILTON *ha;
  int nbc, nbc1, k;
  int nele0, nele1, nk, m, t, r;
  int i, p, np, nq, inp, inq, k0, k1, q;
  int jp, jq, kap, kaq, kp, kq, kp2, kq2;
  int ic, kgp, nb, *bk;
  int ncc, ncc0, ncc1, ncc2;
  int *bs1, nbs1, *bs0, nbs0, *s;
  double a1, a2, d, a, b, *e1, *ham;
  CORR_CONFIG cc, *ccp, *ccp1;
  ARRAY ccfg;
  int *icg, ncg, icg0, icg1, rg;
  typedef struct _MBPT_BASE_ {
    int isym;
    int nbasis, nb2, bmax;
    int *basis;
    double *ene;
  } MBPT_BASE;
  MBPT_BASE mb, *mbp;
  ARRAY base;
  char cname[2000];
  FILE *f;
  double t0, t1, t2;
  int sr, nr;
  
  sr = MyRankMPI();
  nr = NProcMPI();
  
  t0 = clock();
  t0 /= CLOCKS_PER_SEC;

  f = NULL;
  if (sr == 0) {
    f = fopen(fn, "w");
    if (f == NULL) {
      printf("cannot open file %s\n", fn);
      return -1;
    }
  }

  n1 = ConstructNGrid(n1, &nm);
  n2 = ConstructNGrid(n2, &nmp);

  BaseConfig(n, s0, n3, n3g, n4, n4g, n1, nm, n2, nmp,
	     &nbc, &bc, &nbc1, &bc1, &nb, &bk, f);

  for (inp = 0; inp < n1; inp++) {
    np = nm[inp];
    for (kp = 0; kp <= kmax; kp++) {
      if (kp >= np) break;
      kp2 = 2*kp;
      for (jp = kp2-1; jp <= kp2+1; jp += 2) {
	if (jp < 0) continue;
	kap = GetKappaFromJL(jp, kp2);
	i = OrbitalIndex(np, kap, 0.0);
	if (i < 0) return -1;
	for (inq = 0; inq < n2; inq++) {
	  nq = nmp[inq];
	  for (kq = 0; kq <= kmax; kq++) {
	    if (kq >= nq) break;
	    kq2 = 2*kq;
	    for (jq = kq2-1; jq <= kq2+1; jq += 2) {
	      if (jq < 0) continue;
	      kaq = GetKappaFromJL(jq, kq2);
	      i = OrbitalIndex(nq, kaq, 0.0);
	      if (i < 0) return -1;
	    }
	  }
	}
      }
    }
  }

  ArrayInit(&ccfg, sizeof(CORR_CONFIG), 10000);
  ArrayInit(&base, sizeof(MBPT_BASE), MAX_SYMMETRIES);
  ncc = ZerothEnergyConfigSym(n, s0, &e1);
  if (ncc == 0) goto ADDCFG;
  de /= HARTREE_EV;
  k = nbc;
  icg = malloc(sizeof(int)*(n1*(1+n2)*k+2));

  icg[0] = 0;
  t = 1;
  g0 = GetGroup(s0[0]);
  nele0 = g0->n_electrons;
  for (p = 0; p < nbc1; p++) {
    c1 = bc1 + p;
    cc.c = c1;
    cc.ig = p;
    cc.ncs = 0;
    cc.np = c1->shells[0].n;
    cc.inp = IBisect(cc.np, n1, nm);
    if (cc.inp < 0) continue;
    cc.inq = 0;
    cc.nq = 0;
    cc.kp = 0;
    cc.kq = 0;
    a = ZerothEnergyConfig(c1);
    r = ConfigParity(c1);
    if (StrongInteractConfig(ncc, e1, &cc, r, a, de)) {
      ccp = ArrayAppend(&ccfg, &cc, NULL);
    }
  } 
  icg[t++] = ccfg.dim;
  for (p = 0; p < nbc; p++) {
    c1 = bc + p;
    a = ZerothEnergyConfig(c1);
    r = ConfigParity(c1);
    nele1 = c1->n_electrons;
    if (nele1 == nele0-1) {
      for (inp = 0; inp < n1; inp++) {
	np = nm[inp];
	for (kp = 0; kp <= kmax; kp++) {
	  if (kp >= np) break;
	  kp2 = 2*kp;
	  for (jp = kp2 - 1; jp <= kp2 + 1; jp += 2) {
	    if (jp < 0) continue;
	    cc.c = c1;
	    cc.ncs = 0;
	    cc.ig = p;
	    cc.inp = inp;
	    cc.np = np;
	    cc.inq = 0;
	    cc.nq = 0;
	    cc.kp = GetKappaFromJL(jp, kp2);
	    k0 = OrbitalIndex(np, cc.kp, 0);
	    k1 = IBisect(k0, nb, bk);
	    if (k1 >= 0) continue;
	    cc.kq = 0;
	    if (StrongInteractConfig(ncc, e1, &cc, r, a, de)) {
	      ccp = ArrayAppend(&ccfg, &cc, NULL);
	    }
	  }
	}
	icg[t++] = ccfg.dim;
      }
    } else if (nele1 == nele0-2) {
      for (inp = 0; inp < n1; inp++) {
	np = nm[inp];
	for (inq = 0; inq < n2; inq++) {
	  nq = nmp[inq];
	  for (kp = 0; kp <= kmax; kp++) {
	    if (kp >= np) break;
	    kp2 = 2*kp;
	    for (jp = kp2-1; jp <= kp2+1; jp += 2) {
	      if (jp < 0) continue;
	      kap = GetKappaFromJL(jp, kp2);
	      k0 = OrbitalIndex(np, kap, 0);
	      k0 = IBisect(k0, nb, bk);
	      if (k0 >= 0) continue;
	      for (kq = 0; kq <= kmax; kq++) {
		if (kq >= nq) break;
		kq2 = 2*kq;
		for (jq = kq2-1; jq <= kq2+1; jq += 2) {
		  if (jq < 0) continue;
		  kaq = GetKappaFromJL(jq, kq2);
		  if (np == nq && kaq > kap) continue;
		  k1 = OrbitalIndex(nq, kaq, 0);
		  k1 = IBisect(k1, nb, bk);
		  if (k1 >= 0) continue;
		  cc.c = c1;
		  cc.ncs = 0;
		  cc.ig = p;
		  cc.inp = inp;
		  cc.np = np;
		  cc.inq = inq;
		  cc.nq = nq;
		  cc.kp = kap;
		  cc.kq = kaq;
		  if (StrongInteractConfig(ncc, e1, &cc, r, a, de)) {
		    ccp = ArrayAppend(&ccfg, &cc, NULL);
		  }
		}
	      }
	    }
	  }
	  icg[t++] = ccfg.dim;
	}
      }
    }
  }
  free(e1);
  if (ccut <= 0) goto ADDCFG;
 
  ncg = t-1;
  r = ncg/nr;
  icg0 = sr*r;
  if (sr < nr-1) {
    icg1 = (sr+1)*r;
  } else {
    icg1 = ncg;
  }

  for (i = 0; i < MAX_SYMMETRIES; i++) {
    ha = GetHamilton(i);
    nk = ConstructHamiltonDiagonal(i, n, s0, 0);    
    if (nk < 0) continue;
    mb.nbasis = ha->dim;
    mb.isym = i;
    mb.basis = malloc(sizeof(int)*mb.nbasis);
    mb.ene = malloc(sizeof(double)*mb.nbasis);
    memcpy(mb.basis, ha->basis, sizeof(int)*mb.nbasis);
    memcpy(mb.ene, ha->hamilton, sizeof(double)*mb.nbasis);
    mb.bmax = mb.basis[mb.nbasis-1];
    ArrayAppend(&base, &mb, NULL);
  }
  ncc1 = 0;  
  for (rg = icg0; rg < icg1; rg++) {
    ncc0 = icg[rg];
    ncc = 0;
    ccp = ArrayGet(&ccfg, icg[rg]);
    for (p = icg[rg]; p < icg[rg+1]; p++) {
      ccp = ArrayGet(&ccfg, p);
      kgp = GroupIndex(gn);
      AddMBPTConfig(kgp, ccp);
      ncc++;
      if (ncc == ccfg.block || p == icg[rg+1]-1) {
	t1 = t2;
	for (i = 0; i < base.dim; i++) {
	  mbp = ArrayGet(&base, i);
	  nk = ConstructHamiltonDiagonal(mbp->isym, 1, &kgp, 0);
	  if (nk < 0) continue;
	  nbs1 = ha->dim;
	  bs1 = ha->basis;
	  e1 = ha->hamilton;	
	  sym = GetSymmetry(mbp->isym);
	  nbs0 = mbp->nbasis;
	  bs0 = mbp->basis;
	  ham = malloc(sizeof(double)*nbs0*nbs1);	  
	  m = 0;
	  for (q = 0; q < nbs1; q++) {
	    k1 = bs1[q];
	    st = ArrayGet(&(sym->states), k1);
	    r = st->kcfg + ncc0;
	    ccp1 = ArrayGet(&ccfg, r);
	    for (t = 0; t < nbs0; t++) {
	      k0 = bs0[t];
	      ham[m] = HamiltonElement(mbp->isym, k0, k1);
	      m++;
	    }	    
	  }
	  for (q = 0; q < nbs1; q++) {
	    st = ArrayGet(&(sym->states), bs1[q]);
	    r = st->kcfg + ncc0;
	    ccp1 = ArrayGet(&ccfg, r);
	    if (ccp1->ncs) continue;
	    for (r = 0; r < nbs0; r++) {
	      m = q*nbs0 + r;
	      a1 = ham[m];
	      b = a1/(mbp->ene[r] - e1[q]);
	      b = b*b;
	      if (b >= ccut) {
		ccp1->ncs = 1;
		ncc1++;		
		break;
	      }
	    }
	  }
	  free(ham);
	}
	RemoveGroup(kgp);
	ReinitRecouple(0);
	ncc = 0;
	ncc0 += ccfg.block;
      }
    }
  }
#if USE_MPI == 1
  if (nr > 1) {
    ncc1 = 0;
    for (p = 0; p < ccfg.dim; p++) {
      ccp = ArrayGet(&ccfg, p);
      if (ccp->ncs) ncc1++;
    }
    MPI_Allreduce(&ncc1, &ncc0, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    int *ics0, *ics1;
    if (ncc0 > 0) {
      ics1 = malloc(sizeof(int)*ncc0);
      ics0 = malloc(sizeof(int)*ncc0*nr);
    }
    for (i = 0; i < ncc0; i++) ics1[i] = -1;
    i = 0;
    for (p = 0; p < ccfg.dim; p++) {
      ccp = ArrayGet(&ccfg, p);
      if (ccp->ncs) ics1[i++] = p;
    }
    MPI_Allgather(ics1, ncc0, MPI_INT, ics0, ncc0, MPI_INT, MPI_COMM_WORLD);
    for (p = 0; p < ncc0*nr; p++) {
      if (ics0[p] >= 0) {
	ccp = ArrayGet(&ccfg, ics0[p]);
	ccp->ncs = 1;
      }
    }
    if (ncc0) {
      free(ics0);
      free(ics1);
    }
  }
#endif

 ADDCFG:  
  if (sr == 0) {
    kgp = GroupIndex(gn);
    ncc = 0;
    for (p = 0; p < ccfg.dim; p++) {
      ccp = ArrayGet(&ccfg, p);
      if (ccut <= 0 || ccp->ncs) {
	AddMBPTConfig(kgp, ccp);
	ncc++;
      }
    }   
    g1 = GetGroup(kgp);    
    for (p = 0; p < g1->n_cfgs; p++) {
      c1 = GetConfigFromGroup(kgp, p);
      ConstructConfigName(cname, 2000, c1);
      fprintf(f, "CC: %4d   %s\n", p, cname);
    }
    fclose(f);
  }

 DONE:
  for (i = 0; i < nbc; i++) {
    if (bc[i].n_shells > 0) free(bc[i].shells);
  }
  for (i = 0; i < nbc1; i++) {
    if (bc1[i].n_shells > 0) free(bc1[i].shells);
  }
  if (nbc > 0) {
    free(bc);
  }
  if (nbc1 > 0) {
    free(bc1);
  }
  if (nb > 0) {
    free(bk);
  }
  free(icg);
  ArrayFree(&ccfg, NULL);
  for (i = 0; i < base.dim; i++) {
    mbp = ArrayGet(&base, i);
    free(mbp->basis);
    free(mbp->ene);
  }
  ArrayFree(&base, NULL);

  return 0;
}

int PrepRadialBasisMBPT(int nk, int *nkm, int n, int *ng, int **bas) {
  int nb, k, j, k2, i, m, ka;
  ORBITAL *orb;

  nb = 2*nk*n;
  *bas = malloc(sizeof(int)*nb);
  m = 0;
  for (k = 0; k < nk; k++) {
    k2 = 2*k;
    for (j = k2-1; j <= k2+1; j += 2) {
      if (j < 0) continue;
      ka = GetKappaFromJL(j, k2);
      for (i = 0; i < n; i++) {
	if (ng[i] <= k) continue;
	if (nkm && nkm[k] > 0) {
	  if (ng[i] > nkm[k]) continue;
	}
	int ix = OrbitalExistsNoLock(ng[i], ka, 0);
	if (ix < 0) {
	  orb = GetNewOrbitalNoLock(ng[i], ka, 0);
	  ix = orb->idx;
	}
	(*bas)[m] = ix;
	m++;
      }
    }
  }
  nb = m;
  *bas = realloc(*bas, sizeof(int)*nb);
  return nb;
}

void SolveRadialBasisMBPT(int nmax) {
  int n;
  
  n = GetNumOrbitals();
  ResetWidMPI();
#pragma omp parallel default(shared)
  {
    int i, ib, nb;
    double wt0 = WallTime();
    ORBITAL *orb;
    nb = 0;
    for (i = 0; i < n; i++) {
#if USE_MPI == 1
      if (mbpt_nsplit) {
	orb = GetOrbital(i);
	if (orb->wfun != NULL) {
	  continue;
	}
	if (orb->n <= nmax) {
	  orb = GetOrbitalSolved(i);
	  nb++;
	  continue;
	}	
	ib = IdxGet(&mbptjp.ibs, i);
	if (ib < 0) continue;
      }	
#elif USE_MPI == 2
      if (SkipMPI()) continue;
#endif
      orb = GetOrbitalSolved(i);
      nb++;
    }
    double wt1 = WallTime();
    MPrintf(-1, "RadialBasis Time=%11.4E nb=%d\n", wt1-wt0, nb);
  }
}

/* pad the shells in c1 and c2, so that the shell structures are the same */
int PadStates(CONFIG *c1, CONFIG *c2, SHELL **bra, SHELL **ket,
	      SHELL_STATE **sbra, SHELL_STATE **sket) {
  int i, j, m, k, ns;
  SHELL *bra1, *ket1;
  SHELL_STATE *sbra1, *sket1, *s1, *s2;

  s1 = c1->csfs;
  s2 = c2->csfs;

  ns = c1->n_shells + c2->n_shells + 2;
  *bra = malloc(sizeof(SHELL)*ns);
  *ket = malloc(sizeof(SHELL)*ns);
  *sbra = malloc(sizeof(SHELL_STATE)*ns);
  *sket = malloc(sizeof(SHELL_STATE)*ns);
  /* first 2 orbitals are reserved for virtual orbitals */
  for (i = 0; i < 2; i++) {
    (*bra)[i].nq = 0;
    (*ket)[i].nq = 0;
    (*sbra)[i].shellJ = 0;
    (*sbra)[i].totalJ = s1[0].totalJ;
    (*sket)[i].shellJ = 0;
    (*sket)[i].totalJ = s2[0].totalJ;
    (*sbra)[i].nu = 0;
    (*sbra)[i].Nr = 0;
    (*sket)[i].nu = 0;
    (*sket)[i].Nr = 0;
  }
  bra1 = *bra + 2;
  ket1 = *ket + 2;
  sbra1 = *sbra + 2;
  sket1 = *sket + 2;
  m = 0;
  i = 0;
  j = 0;
  while (i < c1->n_shells || j < c2->n_shells) {
    if (i == c1->n_shells) {
      k = -1;
    } else if (j == c2->n_shells) {
      k = 1;
    } else {
      k = CompareShell(c1->shells+i, c2->shells+j);
    }
    if (k == 0) {
      memcpy(bra1+m, c1->shells+i, sizeof(SHELL));
      memcpy(ket1+m, c2->shells+j, sizeof(SHELL));
      memcpy(sbra1+m, s1+i, sizeof(SHELL_STATE));
      memcpy(sket1+m, s2+j, sizeof(SHELL_STATE));
      i++;
      j++;
      m++;
    } else if (k < 0) {
      memcpy(ket1+m, c2->shells+j, sizeof(SHELL));
      memcpy(sket1+m, s2+j, sizeof(SHELL_STATE));
      memcpy(bra1+m, c2->shells+j, sizeof(SHELL));
      memcpy(sbra1+m, s2+j, sizeof(SHELL_STATE));
      bra1[m].nq = 0;
      sbra1[m].shellJ = 0;
      sbra1[m].nu = 0;
      sbra1[m].Nr = 0;
      if (i == c1->n_shells) {
	sbra1[m].totalJ = 0;
      } else {
	sbra1[m].totalJ = s1[i].totalJ;
      }
      j++;
      m++;
    } else if (k > 0) {
      memcpy(ket1+m, c1->shells+i, sizeof(SHELL));
      memcpy(sket1+m, s1+i, sizeof(SHELL_STATE));
      memcpy(bra1+m, c1->shells+i, sizeof(SHELL));
      memcpy(sbra1+m, s1+i, sizeof(SHELL_STATE));
      ket1[m].nq = 0;
      sket1[m].shellJ = 0;
      sket1[m].nu = 0;
      sket1[m].Nr = 0;
      if (j == c2->n_shells) {
	sket1[m].totalJ = 0;
      } else {
	sket1[m].totalJ = s2[j].totalJ;
      }
      i++;
      m++;
    }
  }
  ns = m;

  return ns;
}

int CheckInteraction(int ns, SHELL *bra, SHELL *ket, 
		     int np, int *op, int nm, int *om) {
  int i, k, nqk, ph, j;

  for (i = 0; i < ns; i++) {
    nqk = ket[i].nq;
    k = IsPresent(i, np, op);
    nqk += k;
    k = IsPresent(i, nm, om);
    nqk -= k;
    if (nqk != bra[i].nq) return -1;
  }
  ph = 0;
  for (i = 0; i < np; i++) {
    for (j = op[i]+1; j < ns; j++) {
      ph += ket[j].nq;
    }
  }
  for (i = 0; i < nm; i++) {
    for (j = om[i]+1; j < ns; j++) {
      ph += ket[j].nq;
    }
  }

  return ph;
}

int CheckConfig(int ns, SHELL *ket, int np, int *op, int nm, int *om,
		int nc, CONFIG **cs) {
  int i, k, m, nq;
  CONFIG *c;

  k = 0;
  if (nc > 0) {
    c = cs[nc];
  }

  for (i = 0; i < ns; i++) {
    nq = ket[i].nq;
    m = IsPresent(i, nm, om);
    nq -= m;
    if (nq < 0) return 0;
    m = IsPresent(i, np, op);
    nq += m;
    if (nq > 0) {
      if (nq > GetJFromKappa(ket[i].kappa)+1.0) return 0;
      if (nc > 0) {
	memcpy(c->shells+k, ket+i, sizeof(SHELL));
	c->shells[k].nq = nq;
	k++;
      }
    }
  }
   
  if (k > 0) {
    qsort(c->shells, k, sizeof(SHELL), CompareShellInvert);
    if (c->shells[0].nq > 1) {
      if (c->shells[0].n > cs[nc]->nnrs) return -1;
    } else if (c->n_shells > 1) {
      if (c->shells[1].n > cs[nc]->nnrs) return -1;
    }
    for (i = 0; i < nc; i++) {
      if (k != cs[i]->n_shells) continue;
      m = memcmp(c->shells, cs[i]->shells, sizeof(SHELL)*k);
      if (m == 0) return i;
    }
  }
  return -1;
}

double SumInterp1D(int n, double *z, double *x, double *t, double *y) {
  int i, k0, k1, nk;
  double r, a, b, c, d, e, f, g, h;
  double p1, p2, p3, q1, q2, q3;

  r = 0.0;
  for (i = 0; i < n; i++) {
    r += z[i];
  }
  if (1+r == 1) return 0.0;
  if (n == 1) return r;

  for (k0 = 1; k0 < n; k0++) {
    if (x[k0]-x[k0-1] > 1) break;
  }
  if (k0 == n) {
    k1 = n-2;
    if (k1 < 0) k1 = 0;
    for (i = k1; i < n; i++) {
      t[i] = log(x[i]);
      y[i] = log(fabs(z[i]));
    }
    nk = 2;
    goto END;
  }
  k0--;
  
  for (k1 = n-1; k1 >= k0; k1--) {
    if (z[n-1] > 0) {
      if (z[k1] <= 0 || 1+z[k1] == 1) {
	k1++;
	break;
      }
    } else {
      if (z[k1] >= 0 || 1+z[k1] == 1) {
	k1++;
	break;
      }
    }
  }
  if (k1 < k0) k1 = k0;

  for (i = k0; i < n; i++) {
    t[i] = log(x[i]);
  }
  nk = n - k0;
  for (i = k0; i < k1; i++) {    
    for (a = x[i]+1; a < x[i+1]; a += 1.0) {      
      d = log(a);
      UVIP3P(3, nk, t+k0, z+k0, 1, &d, &b);
      r += b;
    }
  }
  for (i = k1; i < n; i++) {
    y[i] = log(fabs(z[i]));
  }
  nk = n - k1;
  for (i = k1+1; i < n; i++) {
    for (a = x[i-1]+1; a < x[i]; a += 1.0) {     
      d = log(a);
      UVIP3P(3, nk, t+k1, y+k1, 1, &d, &b);
      b = exp(b);
      if (z[n-1] < 0) b = -b;
      r += b;
    }
  }
 END:
  h = 0.0;
  a = 0.0;
  b = 0.0;
  c = 0.0;
  if (mbpt_extra > 0) {
    if (nk > 2) {
      k0 = n-3;
      k1 = n-2;
      i = n-1;
      p1 = t[k0] - t[k1];
      p2 = t[k0]*t[k0] - t[k1]*t[k1];
      p3 = y[k0] - y[k1];
      q1 = t[k1] - t[i];
      q2 = t[k1]*t[k1] - t[i]*t[i];
      q3 = y[k1] - y[i];
      c = (p3*q1 - q3*p1)/(q1*p2 - p1*q2);
      b = (p3 - c*p2)/p1;
      a = y[k0] - b*t[k0] - c*t[k0]*t[k0];      
      if (c <= 0 && 2.0*log(x[i])*c+b < -1.2) {
	g = 1.0;
	d = x[i]+1.0;
	while (g > EPS3) {
	  e = log(d);
	  f = exp(a + b*e + c*e*e);
	  if (z[i] < 0) f = -f;
	  h += f;
	  g = f*d/(-(2.0*e*c+b)-1.0);
	  g = fabs(g/r);
	  d += 1.0;
	}
	g *= fabs(r);
	if (f < 0) g = -g;
	h += g-f;	
      }
    }
    if (h == 0 && nk > 1) {
      k0 = n-2;
      k1 = n-1;
      a = y[k1] - y[k0];
      b = t[k1] - t[k0];
      a = -a/b;      
      if (a > 1.2) {
	a -= 1.0;
	b = pow(x[k1]/(1.0+x[k1]), a);
	h = b*z[k1]*x[k1]/a;	
      }
      b = 0.0;
      c = 0.0;
    }
    if (mbpt_extra == 2) {
      for (i = 0; i < n; i++) {
	printf("  %12.5E %12.5E\n", x[i], z[i]);
      }
      printf("# %12.5E %12.5E %12.5E %12.5E %12.5E\n", a, b, c, h, r);
    }
    r += h;
  }
  return r;
}

double SumInterpH(int n, int *ng, int n2, int *ng2, 
		  double *h, double *h1, double *w) {
  double r, *p, *x, *y, *z;
  int m, i0, i1;
  
  i0 = n+n2;
  x = w + i0;
  y = x + i0;
  z = y + i0;
  for (i0 = 0; i0 < n; i0++) {
    x[i0] = ng[i0];
  }
  r = SumInterp1D(n, h1, x, w, y); 
  p = h;
  for (i0 = 0; i0 < n; i0++) {
    for (i1 = 0; i1 < n2; i1++) {
      x[i1] = ng2[i1] + ng[i0];
    }
    z[i0] = SumInterp1D(n2, p, x, w, y);
    p += n2;
  }
  for (i0 = 0; i0 < n; i0++) {
    x[i0] = ng[i0];
  }
  r += SumInterp1D(n, z, x, w, y);
  
  return r;
}

#define MKK 20

void FixTotalJ(int ns, SHELL_STATE *st, SHELL *s, CONFIG *c, int m) {
  SHELL_STATE *st0;
  int i, j;

  st0 = c->csfs + m*c->n_shells;
  i = 0;
  j = 0;
  while (i < ns && j < c->n_shells) {
    st[i].totalJ = st0[j].totalJ;
    if (s[i].nq > 0) {
      st[i].shellJ = st0[j].shellJ;
      st[i].nu = st0[j].nu;
      st[i].Nr = st0[j].Nr;
      j++;
    }
    i++;
  }
  for (; i < ns; i++) {
    st[i].shellJ = 0;
    st[i].nu = 0;
    st[i].Nr = 0;
    st[i].totalJ = 0;
  }
  if (j < c->n_shells) {
    printf("Error in FixTotalJ, aborting\n");
    exit(1);
  }
}

void H3rd0(MBPT_EFF *meff, int ia, int ib, double de, double h, int i0, int md) {
  int i, m, j;
  double *e0, *r0, *h0, **hab, **hba, a, c;

  if (mbpt_3rd == 0) return;

  if (md == 1) {
    hab = meff->hab1;
    hba = meff->hba1;
  } else {
    hab = meff->hab;
    hba = meff->hba;
  }
  
  if (meff == NULL) return;
  if (meff->nbasis <= 0) return;

  e0 = meff->e0;
  h0 = meff->h0;
  for (i = 0; i < ia; i++) {
    if (i != ib) {
      m = ia*(ia+1)/2 + i;
      if (meff->hab1[m] == NULL) continue;
      if (i < ib) {
	j = ib*(ib+1)/2 + i;
      } else {
	j = i*(i+1)/2 + ib;
      }
      a = h0[j];
      c = h*a/(de*(de+e0[i]-e0[ib]));
#if CPMEFF == 0
#pragma omp atomic
#endif
      hba[m][i0] -= c;
    }
  }
  for (i = ia; i < meff->nbasis; i++) {
    if (i != ib) {
      m = i*(i+1)/2 + ia;
      if (meff->hab1[m] == NULL) continue;
      if (i < ib) {
	j = ib*(ib+1)/2 + i;
      } else {
	j = i*(i+1)/2 + ib;
      }
      a = h0[j];
      c = h*a/(de*(de+e0[i]-e0[ib]));
#if CPMEFF == 0
#pragma omp atomic
#endif
      hab[m][i0] -= c;
    }
  }
}
    
void H22Term(MBPT_EFF **meff, CONFIG *c0, CONFIG *c1,
	     int ns, SHELL *bra, SHELL *ket, 
	     SHELL_STATE *sbra, SHELL_STATE *sket,
	     int mst, int *bst, int *kst,
	     INTERACT_SHELL *s, int ph, int *ks1, int *ks2,
	     FORMULA *fm, double **a, int i0) {
  int m, kk1, kk2, kmin1, kmin2, kmax1, kmax2;
  int mkk1, mkk2, mkk, k, i1, ng, i1g;
  int q0, q1, m0, m1, ms0, ms1, s0, s1;
  double c, y, sd1, sd2, se1, se2;
  double a1[MKK], a2[MKK], *h1, *h2, d1, d2;
  ORBITAL *orb;

  int md;
  if (fm->j1 < 0 || fm->j2 < 0) {
    md = 0;
    fm->j1 = s[1].j;
    fm->j2 = s[3].j;
  } else if (s[1].j != fm->j1 || s[3].j != fm->j2) {
    md = 1;
    fm->j1 = s[1].j;
    fm->j2 = s[3].j;
  } else {
    md = 2;
  }
  if (md == 0) {
    TriadsZ(2, 2, fm);	      
    RecoupleTensor(8, s, fm);
  }  
  kmin1 = abs(s[0].j-s[1].j);
  kk1 = abs(s[2].j-s[3].j);
  kmin1 = Max(kmin1, kk1);
  kmax1 = s[0].j + s[1].j;
  kk1 = s[2].j + s[3].j;
  kmax1 = Min(kmax1, kk1);
  kmin2 = abs(s[4].j-s[5].j);
  kk2 = abs(s[6].j-s[7].j);
  kmin2 = Max(kmin2, kk2);
  kmax2 = s[4].j + s[5].j;
  kk2 = s[6].j + s[7].j;
  kmax2 = Min(kmax2, kk2);
  
  kk1 = 2*(MKK-1);
  kmax1 = Min(kmax1, kk1);
  kmax2 = Min(kmax2, kk1);
  if (kmax1 < kmin1) return;
  if (kmax2 < kmin2) return;
  if (md <= 1) {
    FixJsZ(s, fm);
    for (kk1 = kmin1; kk1 <= kmax1; kk1 += 2) {
      mkk1 = kk1/2;
      mkk = mkk1*MKK;
      for (kk2 = kmin2; kk2 <= kmax2; kk2 += 2) {
	mkk2 = kk2/2;
	fm->js[9] = kk1;
	fm->js[10] = kk1;
	fm->js[11] = kk2;
	fm->js[12] = kk2;
	fm->js[13] = 0;
	fm->js[14] = 0;
	fm->js[15] = 0;	
	for (k = 0; k < mst; k++) {
	  q0 = bst[k];
	  q1 = kst[k];
	  FixTotalJ(ns, sbra, bra, c0, q0);
	  FixTotalJ(ns, sket, ket, c1, q1);
	  ms0 = c0->symstate[q0];
	  UnpackSymState(ms0, &s0, &m0);
	  ms1 = c1->symstate[q1];
	  UnpackSymState(ms1, &s1, &m1);	  
	  EvaluateTensor(ns, sbra, sket, s, 1, fm);
	  if (IsOdd(ph)) fm->coeff = -fm->coeff;
	  a[mkk+mkk2][k] = fm->coeff;
	}
      }
    }
  }
  /* 
  ** these conditions must be placed after the above block,
  ** so that the recouping coeff. are calculated the first
  ** time even these parity conditions are not satisfied.
  */
  if (IsOdd((s[0].kl+s[1].kl+s[2].kl+s[3].kl)/2)) return;
  if (IsOdd((s[4].kl+s[5].kl+s[6].kl+s[7].kl)/2)) return;
  for (kk1 = kmin1; kk1 <= kmax1; kk1 += 2) {
    m = 0;
    mkk1 = kk1/2;
    mkk = mkk1*MKK;
    for (kk2 = kmin2; kk2 <= kmax2; kk2 += 2) {
      mkk2 = kk2/2;
      for (k = 0; k < mst; k++) {
	if (fabs(a[mkk+mkk2][k]) > EPS30) {
	  m = 1;
	  break;
	}
      }
      if (m == 1) break;
    }
    if (m) {
      SlaterTotal(&sd1, &se1, NULL, ks1, kk1, 0);
      a1[mkk1] = sd1 + se1;
    } else {
      a1[mkk1] = 0.0;
    }
  }
  for (kk2 = kmin2; kk2 <= kmax2; kk2 += 2) {
    m = 0;
    mkk2 = kk2/2;
    for (kk1 = kmin1; kk1 <= kmax1; kk1 += 2) {
      mkk1 = kk1/2;
      mkk = mkk1*MKK+mkk2;
      for (k = 0; k < mst; k++) {
	if (fabs(a[mkk][k]) > EPS30) {
	  m = 1;
	  break;
	}
      }      
      if (m == 1) break;
    }
    if (m) {
      SlaterTotal(&sd2, &se2, NULL, ks2, kk2, 0);
      a2[mkk2] = sd2 + se2;
    } else {
      a2[mkk2] = 0.0;
    }
  }

  orb = GetOrbital(ks2[2]);
  d1 = orb->energy + orb->qed;
  orb = GetOrbital(ks2[3]);
  d1 += orb->energy + orb->qed;
  orb = GetOrbital(ks2[0]);
  d1 -= orb->energy + orb->qed;
  orb = GetOrbital(ks2[1]);
  d1 -= orb->energy + orb->qed;
  orb = GetOrbital(ks1[0]);
  d2 = orb->energy + orb->qed;
  orb = GetOrbital(ks1[1]);
  d2 += orb->energy + orb->qed;
  orb = GetOrbital(ks1[2]);
  d2 -= orb->energy + orb->qed;
  orb = GetOrbital(ks1[3]);
  d2 -= orb->energy + orb->qed;
  /*
  d1 = GetOrbital(ks2[2])->energy + GetOrbital(ks2[3])->energy;
  d1 -= GetOrbital(ks2[0])->energy + GetOrbital(ks2[1])->energy;
  d2 = GetOrbital(ks1[0])->energy + GetOrbital(ks1[1])->energy;
  d2 -= GetOrbital(ks1[2])->energy + GetOrbital(ks1[3])->energy;
  */
  for (k = 0; k < mst; k++) {
    q0 = bst[k];
    q1 = kst[k];
    ms0 = c0->symstate[q0];
    UnpackSymState(ms0, &s0, &m0);
    ms1 = c1->symstate[q1];
    UnpackSymState(ms1, &s1, &m1);
   
    c = 0.0;
    for (kk1 = kmin1; kk1 <= kmax1; kk1 += 2) {
      mkk1 = kk1/2;
      mkk = mkk1*MKK;
      for (kk2 = kmin2; kk2 <= kmax2; kk2 += 2) {
	mkk2 = kk2/2;
	y = a[mkk+mkk2][k];
	if (fabs(y) < EPS30) continue;
	if (fabs(a1[mkk1]) < EPS30) continue;
	if (fabs(a2[mkk2]) < EPS30) continue;	  
	y /= sqrt((kk1+1.0)*(kk2+1.0));
	if (IsOdd((kk1+kk2)/2)) y = -y;      
	c += y*a1[mkk1]*a2[mkk2];
      }
    }
    if (fabs(c) < EPS30) continue;
    m = m1*(m1+1)/2 + m0;
    if (i0 >= 0) {
      h1 = meff[s0]->hab[m];
      h2 = meff[s0]->hba[m];
      ng = meff[s0]->n * meff[s0]->n2;
      i1 = i0;
      H3rd0(meff[s0], m0, m1, d1, c, i1, 2);
      if (m0 != m1) {
	H3rd0(meff[s0], m1, m0, d2, c, i1, 2);
      }
    } else {
      h1 = meff[s0]->hab1[m];
      h2 = meff[s0]->hba1[m];
      ng = meff[s0]->n;
      i1 = -(i0+1);
      H3rd0(meff[s0], m0, m1, d1, c, i1, 1);
      if (m0 != m1) {
	H3rd0(meff[s0], m1, m0, d2, c, i1, 1);
      }
    }
    sd1 = c/d1;
    sd2 = c/d2;
#if CPMEFF == 0
#pragma omp atomic
#endif
    h1[i1] += sd1;
#if CPMEFF == 0
#pragma omp atomic
#endif
    h2[i1] += sd2;
    c /= d1*d2;
    i1g = i1+ng;
#if CPMEFF == 0
#pragma omp atomic
#endif
    h1[i1g] += c;
#if CPMEFF == 0
#pragma omp atomic
#endif
    h2[i1g] += c;
  }
}

void H12Term(MBPT_EFF **meff, CONFIG *c0, CONFIG *c1,
	     int ns, SHELL *bra, SHELL *ket, 
	     SHELL_STATE *sbra, SHELL_STATE *sket,
	     int mst, int *bst, int *kst,
	     INTERACT_SHELL *s, int ph, int *ks, int k0, int k1,
	     FORMULA *fm, double **a, int i0) {  
  int kk, kk2, kmin, kmax, ng, i0g;
  int q0, q1, k, m0, m1, s0, s1, m, ms0, ms1;  
  double c, r1, y, sd, se, yk[MKK], *h1, *h2, d1, d2;
  ORBITAL *orb;

  int md;
  if (fm->j1 < 0) {
    md = 0;
    fm->j1 = s[1].j;
  } else if (s[1].j != fm->j1) {
    md = 1;
    fm->j1 = s[1].j;
  } else {
    md = 2;
  }
  if (md == 0) {
    TriadsZ(1, 2, fm);
    RecoupleTensor(6, s, fm);
  }
  if (s[0].j != s[1].j) return;
  kmin = abs(s[2].j-s[3].j);
  kk = abs(s[4].j-s[5].j);
  kmin = Max(kmin, kk);
  kmax = s[2].j + s[3].j;
  kk = s[4].j + s[5].j;
  kmax = Min(kmax, kk);
  kk = 2*(MKK-1);
  kmax = Min(kmax, kk);
  if (kmax < kmin) return;
  if (md <= 1) {
    FixJsZ(s, fm);
    for (kk = kmin; kk <= kmax; kk += 2) {
      kk2 = kk/2;
      fm->js[7] = 0;
      fm->js[8] = kk;
      fm->js[9] = kk;
      fm->js[10] = 0;
      fm->js[11] = 0;
      for (k = 0; k < mst; k++) {
	q0 = bst[k];
	q1 = kst[k];
	FixTotalJ(ns, sbra, bra, c0, q0);
	FixTotalJ(ns, sket, ket, c1, q1);
	EvaluateTensor(ns, sbra, sket, s, 1, fm);
	if (IsOdd(ph)) fm->coeff = -fm->coeff;
	a[kk2][k] = fm->coeff;
      }
    }
  }
  /* 
  ** these conditions must be placed after the above block,
  ** so that the recouping coeff. are calculated the first
  ** time even these parity conditions are not satisfied.
  */
  if (IsOdd((s[0].kl+s[1].kl)/2)) return;
  if (IsOdd((s[2].kl+s[3].kl+s[4].kl+s[5].kl)/2)) return;

  orb = GetOrbital(ks[2]);
  d1 = orb->energy + orb->qed;
  orb = GetOrbital(ks[3]);
  d1 += orb->energy + orb->qed;
  orb = GetOrbital(ks[0]);
  d1 -= orb->energy + orb->qed;
  orb = GetOrbital(ks[1]);
  d1 -= orb->energy + orb->qed;
  orb = GetOrbital(k0);
  d2 = orb->energy + orb->qed;
  orb = GetOrbital(k1);
  d2 -= orb->energy + orb->qed;
  /*
  d1 = GetOrbital(ks[2])->energy + GetOrbital(ks[3])->energy;
  d1 -= GetOrbital(ks[0])->energy + GetOrbital(ks[1])->energy;
  d2 = GetOrbital(k0)->energy - GetOrbital(k1)->energy;
  */
  for (kk = kmin; kk <= kmax; kk += 2) {
    kk2 = kk/2;
    yk[kk2] = 0.0;
    for (k = 0; k < mst; k++) {
      if (a[kk2][k]) break;
    }
    if (k == mst) continue;
    SlaterTotal(&sd, &se, NULL, ks, kk, 0);  
    yk[kk2] = sd + se;
    if (IsOdd(kk2)) yk[kk2] = -yk[kk2];
  }
  for (k = 0; k < mst; k++) {
    q0 = bst[k];
    q1 = kst[k];
    ms0 = c0->symstate[q0];
    UnpackSymState(ms0, &s0, &m0);
    ms1 = c1->symstate[q1];
    UnpackSymState(ms1, &s1, &m1);
    c = 0.0;
    for (kk = kmin; kk <= kmax; kk += 2) {
      kk2 = kk/2;
      if (a[kk2][k] == 0) continue;
      if (yk[kk2] == 0) continue;
      y = a[kk2][k];
      y /= sqrt(kk+1.0);
      c += y*yk[kk2];
    }
    if (fabs(c) < EPS30) continue;
    ResidualPotential(&r1, k0, k1);
    r1 += QED1E(k0, k1);
    r1 *= sqrt(s[0].j+1.0);
    c *= r1;
    /* minus sign is from the definition of Z^k */
    c = -c;    
    ng = meff[s0]->n;
    if (m0 <= m1) {
      m = m1*(m1+1)/2 + m0;
      h1 = meff[s0]->hab1[m];
      h2 = meff[s0]->hba1[m];
      H3rd0(meff[s0], m0, m1, d1, c, i0, 1);
      if (m0 != m1) {
	H3rd0(meff[s0], m1, m0, d2, c, i0, 1);
      }
    } else {
      m = m0*(m0+1)/2 + m1;
      h1 = meff[s0]->hba1[m];
      h2 = meff[s0]->hab1[m];
      H3rd0(meff[s0], m0, m1, d1, c, i0, 1);
      H3rd0(meff[s0], m1, m0, d2, c, i0, 1);
    }
    sd = c/d1;
    se = c/d2;
#if CPMEFF == 0
#pragma omp atomic
#endif
    h1[i0] += sd;
#if CPMEFF == 0
#pragma omp atomic
#endif
    h2[i0] += se;
    c /= d1*d2;
    i0g = i0 + ng;
#if CPMEFF == 0
#pragma omp atomic
#endif
    h1[i0g] += c;
#if CPMEFF == 0
#pragma omp atomic
#endif
    h2[i0g] += c;
  }
}

void TR12Term(MBPT_TR *mtr, CONFIG *c0, CONFIG *c1,
	      int ns, SHELL *bra, SHELL *ket, 
	      SHELL_STATE *sbra, SHELL_STATE *sket,
	      int mst, int *bst, int *kst,
	      INTERACT_SHELL *s, int ph, int *ks, int k0, int k1,
	      FORMULA *fm, double **a, int i0, int ng) {  
  int kk, kk2, kmin, kmax, gauge, n, p1, j0, j1;
  int q0, q1, k, m0, m1, s0, s1, m, ms0, ms1, is0, is1, p, p0;  
  double c, *r1, sd, se, yk[MKK], d2, rt;
  ORBITAL *orb;  

  int md;
  if (fm->j1 < 0) {
    md = 0;
    fm->j1 = s[1].j;
  } else if (fm->j1 != s[1].j) {
    md = 1;
    fm->j1 = s[1].j;
  } else {
    md = 2;
  }
  if (md == 0) {
    TriadsZ(1, 2, fm);
    RecoupleTensor(6, s, fm);
  }
  if (abs(s[0].j - s[1].j) > 2*mbpt_tr.mktr) return;
  kmin = abs(s[2].j-s[3].j);
  kk = abs(s[4].j-s[5].j);
  kmin = Max(kmin, kk);
  kmax = s[2].j + s[3].j;
  kk = s[4].j + s[5].j;
  kmax = Min(kmax, kk);
  kk = 2*(MKK-1);
  kmax = Min(kmax, kk);
  if (kmax < kmin) return;
  if (md <= 1) {
    FixJsZ(s, fm);
    for (p = 1; p <= mbpt_tr.mktr; p++) {
      for (kk = kmin; kk <= kmax; kk += 2) {
	kk2 = kk/2;
	fm->js[7] = 2*p;
	fm->js[8] = kk;
	fm->js[9] = kk;
	fm->js[10] = 0;
	fm->js[11] = 2*p;
	for (k = 0; k < mst; k++) {
	  q0 = bst[k];
	  q1 = kst[k];
	  FixTotalJ(ns, sbra, bra, c0, q0);
	  FixTotalJ(ns, sket, ket, c1, q1);
	  EvaluateTensor(ns, sbra, sket, s, 1, fm);
	  if (IsOdd(ph)) fm->coeff = -fm->coeff;
	  if (IsOdd(kk2)) fm->coeff = -fm->coeff;
	  q0 = k*mbpt_tr.mktr + (p-1);
	  a[kk2][q0] = fm->coeff;
	  a[kk2][q0] /= -sqrt((2.0*p+1.0)*(kk+1.0));
	}
      }
    }
  }
  /* 
  ** these conditions must be placed after the above block,
  ** so that the recouping coeff. are calculated the first
  ** time even these parity conditions are not satisfied.
  */
  if (IsOdd((s[2].kl+s[3].kl+s[4].kl+s[5].kl)/2)) return;
  
  orb = GetOrbital(ks[2]);
  d2 = orb->energy + orb->qed;
  orb = GetOrbital(ks[3]);
  d2 += orb->energy + orb->qed;
  orb = GetOrbital(ks[0]);
  d2 -= orb->energy + orb->qed;
  orb = GetOrbital(ks[1]);
  d2 -= orb->energy + orb->qed;
  gauge = GetTransitionGauge();

  for (kk = kmin; kk <= kmax; kk += 2) {
    kk2 = kk/2;
    SlaterTotal(&sd, &se, NULL, ks, kk, 0);  
    yk[kk2] = sd + se;
  }
  for (p = 1; p <= mbpt_tr.mktr; p++) {
    if (abs(s[0].j - s[1].j) > 2*p) continue;
    p0 = (s[0].kl + s[1].kl)/2 + p;
    if (IsEven(p0)) {
      n = MultipoleRadialFRGrid(&r1, -p, k0, k1, gauge);
    } else {
      n = MultipoleRadialFRGrid(&r1, p, k0, k1, gauge);
    }
    for (k = 0; k < mst; k++) {
      q0 = bst[k];
      q1 = kst[k];
      ms0 = c0->symstate[q0];
      UnpackSymState(ms0, &s0, &m0);
      ms1 = c1->symstate[q1];
      UnpackSymState(ms1, &s1, &m1);
      is0 = 2*mbpt_tr.mktr*s0 + 2*(p-1);
      if (IsOdd(p0)) is0++;
      is1 = IBisect(s1, mtr[is0].nsym1, mtr[is0].isym1);
      if (is1 < 0) continue;
      q0 = k*mbpt_tr.mktr + (p-1);
      c = 0.0;
      for (kk = kmin; kk <= kmax; kk += 2) {
	kk2 = kk/2;
	c += a[kk2][q0]*yk[kk2];
      }
      for (p1 = 0; p1 < n; p1++) {
	q1 = ((m0*mtr[is0].sym1[is1]->n_states + m1)*n + p1)*ng + i0;
	rt = c*r1[p1]/d2;
#if CPMTR == 0
#pragma omp atomic
#endif
	mtr[is0].tma[is1][q1] += rt;
      }
    }
  }
	
  kk = ks[0];
  ks[0] = ks[3];
  ks[3] = kk;
  kk = ks[1];
  ks[1] = ks[2];
  ks[2] = kk;
  for (kk = kmin; kk <= kmax; kk += 2) {
    kk2 = kk/2;
    SlaterTotal(&sd, &se, NULL, ks, kk, 0);  
    yk[kk2] = sd + se;
  }
  for (p = 1; p <= mbpt_tr.mktr; p++) {
    if (abs(s[0].j - s[1].j) > 2*p) continue;
    p0 = (s[0].kl + s[1].kl)/2 + p;
    if (IsEven(p0)) {
      n = MultipoleRadialFRGrid(&r1, -p, k1, k0, gauge);
    } else {
      n = MultipoleRadialFRGrid(&r1, p, k1, k0, gauge);
    }
    for (k = 0; k < mst; k++) {
      q0 = bst[k];
      q1 = kst[k];
      ms0 = c0->symstate[q0];
      UnpackSymState(ms0, &s0, &m0);
      DecodePJ(s0, NULL, &j0);
      ms1 = c1->symstate[q1];
      UnpackSymState(ms1, &s1, &m1);
      DecodePJ(s1, NULL, &j1);
      is0 = 2*mbpt_tr.mktr*s0 + 2*(p-1);
      if (IsOdd(p0)) is0++;
      is1 = IBisect(s1, mtr[is0].nsym1, mtr[is0].isym1);
      if (is1 < 0) continue;
      q0 = k*mbpt_tr.mktr + (p-1);
      c = 0.0;
      for (kk = kmin; kk <= kmax; kk += 2) {
	kk2 = kk/2;
	c += a[kk2][q0]*yk[kk2];
      }
      if (IsOdd(abs(j0-j1+s[2].j-s[3].j+s[4].j-s[5].j+s[0].j-s[1].j)/2)) {
	c = -c;
      }
      for (p1 = 0; p1 < n; p1++) {
	q1 = ((m0*mtr[is0].sym1[is1]->n_states + m1)*n + p1)*ng + i0;
	rt = c*r1[p1]/d2;
#if CPMTR == 0
#pragma omp atomic
#endif
	mtr[is0].rma[is1][q1] += rt;
      }
    }
  }
}

void H11Term(MBPT_EFF **meff, CONFIG *c0, CONFIG *c1, 
	     int ns, SHELL *bra, SHELL *ket,
	     SHELL_STATE *sbra, SHELL_STATE *sket, 
	     int mst, int *bst, int *kst,
	     INTERACT_SHELL *s, int ph, int k0, int k1, int k2, int k3,
	     FORMULA *fm, double *a, int i0) {
  double y, d1, d2, r1, r2, *h1, *h2, cd1, cd2;
  int q0, q1, k, m0, m1, m, ms0, ms1, s0, s1, ng, i0g;
  ORBITAL *orb;

  /* setup recouple tensor */
  int md; 
  if (fm->j1 < 0) {
    md = 0;
    fm->j1 = s[1].j;
  } else if (fm->j1 != s[1].j) {
    md = 1;
    fm->j1 = s[1].j;
  } else {
    md = 2;
  }
  if (md == 0) {
    TriadsZ(1, 1, fm);
    RecoupleTensor(4, s, fm);
  }
  if (s[0].j != s[1].j) return;
  if (s[2].j != s[3].j) return;
  if (md <= 1) {
    fm->j1 = s[1].j;
    FixJsZ(s, fm);
    fm->js[5] = 0;
    fm->js[6] = 0;
    fm->js[7] = 0;
    for (k = 0; k < mst; k++) {
      q0 = bst[k];
      q1 = kst[k];
      FixTotalJ(ns, sbra, bra, c0, q0);
      FixTotalJ(ns, sket, ket, c1, q1);	
      EvaluateTensor(ns, sbra, sket, s, 1, fm);
      a[k] = fm->coeff;
      if (IsOdd(ph)) a[k] = -a[k];
    }
  }
  
  /* 
  ** these conditions must be placed after the above block,
  ** so that the recouping coeff. are calculated the first
  ** time even these parity conditions are not satisfied.
  */
  if (IsOdd((s[0].kl+s[1].kl)/2)) return;
  if (IsOdd((s[2].kl+s[3].kl)/2)) return;

  y = 0;
  for (k = 0; k < mst; k++) {
    if (a[k]) break;
  }
  if (k == mst) return;
  
  ResidualPotential(&r1, k0, k1);
  r1 += QED1E(k0, k1);
  r1 *= sqrt(s[0].j+1.0);
  if (k3 != k0 || k2 != k1) {
    ResidualPotential(&r2, k2, k3);
    r2 += QED1E(k2, k3);
    r2 *= sqrt(s[3].j+1.0);
  } else {
    r2 = r1;
  }

  orb = GetOrbital(k3);
  d1 = orb->energy + orb->qed;
  orb = GetOrbital(k2);
  d1 -= orb->energy + orb->qed;
  orb = GetOrbital(k0);
  d2 = orb->energy + orb->qed;
  orb = GetOrbital(k1);
  d2 -= orb->energy + orb->qed;
  for (k = 0; k < mst; k++) {
    q0 = bst[k];
    q1 = kst[k];
    ms0 = c0->symstate[q0];
    UnpackSymState(ms0, &s0, &m0);
    ms1 = c1->symstate[q1];
    UnpackSymState(ms1, &s1, &m1);
    if (a[k] == 0) continue;      
    m = m1*(m1+1)/2 + m0;
    h1 = meff[s0]->hab1[m];
    h2 = meff[s0]->hba1[m];
    ng = meff[s0]->n;
    y = r1*r2*a[k];
    cd1 = y/d1;
    cd2 = y/d2;
#if CPMEFF == 0
#pragma omp atomic
#endif
    h1[i0] += cd1;
#if CPMEFF == 0
#pragma omp atomic
#endif
    h2[i0] += cd2;
    H3rd0(meff[s0], m0, m1, d1, y, i0, 1);
    if (m0 != m1) {
      H3rd0(meff[s0], m1, m0, d2, y, i0, 1);
    }
    y /= d1*d2;
    i0g = i0 + ng;
#if CPMEFF == 0
#pragma omp atomic
#endif
    h1[i0g] += y;
#if CPMEFF == 0
#pragma omp atomic
#endif
    h2[i0g] += y;
  }
}

void TR11Term(MBPT_TR *mtr, CONFIG *c0, CONFIG *c1, 
	      int ns, SHELL *bra, SHELL *ket,
	      SHELL_STATE *sbra, SHELL_STATE *sket, 
	      int mst, int *bst, int *kst,
	      INTERACT_SHELL *s, int ph, int k0, int k1, int k2, int k3,
	      FORMULA *fm, double *a, int i0, int ng) {
  double d2, *r1, r2, c, d;
  int q0, q1, k, m0, m1, m, ms0, ms1, s0, s1, p, p0, n, is0, is1;
  int gauge, p1, j0, j1;
  ORBITAL *orb;
  
  /* setup recouple tensor */
  int md;  
  if (fm->j1 < 0) {
    md = 0;
    fm->j1 = s[1].j;
  } else if (fm->j1 != s[1].j) {
    md = 1;
    fm->j1 = s[1].j;
  } else {
    md = 2;
  }
  if (md == 0) {
    TriadsZ(1, 1, fm);
    RecoupleTensor(4, s, fm);
  }
  if (abs(s[0].j - s[1].j) > 2*mbpt_tr.mktr) return;
  if (s[2].j != s[3].j) return;
  if (md <= 1) {
    FixJsZ(s, fm);
    for (p = 1; p <= mbpt_tr.mktr; p++) {
      fm->js[5] = 2*p;
      fm->js[6] = 0;
      fm->js[7] = 2*p;
      for (k = 0; k < mst; k++) {
	q0 = bst[k];
	q1 = kst[k];
	FixTotalJ(ns, sbra, bra, c0, q0);
	FixTotalJ(ns, sket, ket, c1, q1);	
	EvaluateTensor(ns, sbra, sket, s, 1, fm);
	q0 = k*mbpt_tr.mktr + (p-1);
	a[q0] = fm->coeff;
	if (IsOdd(ph)) a[q0] = -a[q0];
	a[q0] /= sqrt(2.0*p+1.0);
      }
    }
  }
  
  /* 
  ** these conditions must be placed after the above block,
  ** so that the recouping coeff. are calculated the first
  ** time even these parity conditions are not satisfied.
  */
  if (IsOdd((s[2].kl+s[3].kl)/2)) return;

  ResidualPotential(&r2, k2, k3);
  r2 += QED1E(k2, k3);
  r2 *= sqrt(s[3].j+1.0);

  orb = GetOrbital(k3);
  d2 = orb->energy + orb->qed;
  orb = GetOrbital(k2);
  d2 -= orb->energy + orb->qed;
  /*
  d2 = GetOrbital(k3)->energy - GetOrbital(k2)->energy;
  */
  gauge = GetTransitionGauge();
  for (p = 1; p <= mbpt_tr.mktr; p++) {
    if (abs(s[0].j - s[1].j) > 2*p) continue;
    p0 = (s[0].kl+s[1].kl)/2 + p;
    if (IsEven(p0)) {
      n = MultipoleRadialFRGrid(&r1, -p, k0, k1, gauge);
    } else {
      n = MultipoleRadialFRGrid(&r1, p, k0, k1, gauge);
    }
    for (k = 0; k < mst; k++) {
      q0 = bst[k];
      q1 = kst[k];
      ms0 = c0->symstate[q0];
      UnpackSymState(ms0, &s0, &m0);
      ms1 = c1->symstate[q1];
      UnpackSymState(ms1, &s1, &m1);
      is0 = 2*mbpt_tr.mktr*s0 + 2*(p-1);
      if (IsOdd(p0)) is0++;
      is1 = IBisect(s1, mtr[is0].nsym1, mtr[is0].isym1);
      if (is1 < 0) continue;
      q0 = k*mbpt_tr.mktr + (p-1);
      for (p1 = 0; p1 < n; p1++) {
	q1 = ((m0*mtr[is0].sym1[is1]->n_states + m1)*n + p1)*ng + i0;
	d = a[q0]*r1[p1]*r2/d2;
#if CPMTR == 0
#pragma omp atomic
#endif
	mtr[is0].tma[is1][q1] += d;
      }
    }
    if (IsEven(p0)) {
      n = MultipoleRadialFRGrid(&r1, -p, k1, k0, gauge);
    } else {
      n = MultipoleRadialFRGrid(&r1, p, k1, k0, gauge);
    }
    for (k = 0; k < mst; k++) {
      q0 = bst[k];
      q1 = kst[k];
      ms0 = c0->symstate[q0];
      UnpackSymState(ms0, &s0, &m0);
      DecodePJ(s0, NULL, &j0);
      ms1 = c1->symstate[q1];
      UnpackSymState(ms1, &s1, &m1);
      DecodePJ(s1, NULL, &j1);
      is0 = 2*mbpt_tr.mktr*s0 + 2*(p-1);
      if (IsOdd(p0)) is0++;
      is1 = IBisect(s1, mtr[is0].nsym1, mtr[is0].isym1);
      if (is1 < 0) continue;
      q0 = k*mbpt_tr.mktr + (p-1);
      if (IsOdd((abs(j0-j1+s[0].j-s[1].j)/2))) c = -r2;
      else c = r2;
      for (p1 = 0; p1 < n; p1++) {
	q1 = ((m0*mtr[is0].sym1[is1]->n_states + m1)*n + p1)*ng + i0;	
	d = a[q0]*r1[p1]*c/d2;
#if CPMTR == 0
#pragma omp atomic
#endif
	mtr[is0].rma[is1][q1] += d;
      }
    }
  }
}

void DeltaH22M2Loop(MBPT_EFF **meff, CONFIG *c0, CONFIG *c1, int ns, 
		    SHELL *bra, SHELL *ket,
		    SHELL_STATE *sbra, SHELL_STATE *sket,
		    int mst, int *bst, int *kst,
		    IDXARY *ib0, IDXARY *ib1, IDXARY *ing, 
		    IDXARY *ing2, int ph,
		    int ia, int ib, int ic, int id, int ik, int im,
		    FORMULA *fm, double **a) {
  double c, d1, d2;
  int ip, iq;
  int ks1[4], ks2[4];
  int i, i1, i2, m1, m2;
  INTERACT_SHELL s[8];
  ORBITAL *o[8];
  
  ip = im;
  iq = ik;
  ks1[0] = ib0->d[ia];
  ks1[1] = ib0->d[ib];
  ks1[2] = ib1->d[ik];
  ks1[3] = ib1->d[im];
  ks2[0] = ib1->d[iq];
  ks2[1] = ib1->d[ip];
  ks2[2] = ib0->d[ic];
  ks2[3] = ib0->d[id];
  s[0].index = ia+2;
  s[1].index = 0;
  s[2].index = ib+2;
  s[3].index = 1;
  s[4].index = 0;
  s[5].index = ic+2;
  s[6].index = 1;
  s[7].index = id+2;
  if (ik == im) {
    s[1].index = 1;
    s[4].index = 1;
  }
  for (i = 0; i < 8; i++) {
    if (i == 1 || i == 4) {
      o[i] = GetOrbital(ib1->d[ik]);
    } else if (i == 3 || i == 6) {
      o[i] = GetOrbital(ib1->d[im]);
    } else {
      o[i] = GetOrbital(ib0->d[s[i].index-2]);
    }
    if (IsEven(i)) s[i].n = 1;
    else s[i].n = -1;
    s[i].kappa = o[i]->kappa;
    GetJLFromKappa(s[i].kappa, &(s[i].j), &(s[i].kl));
    s[i].nq_bra = bra[s[i].index].nq;
    s[i].nq_ket = ket[s[i].index].nq;
    s[i].index = ns-s[i].index-1;
  }
  
  m1 = o[1]->n;
  m2 = o[3]->n;
  int vk = ik;
  if (m1 > m2) {
    i = m1;
    m1 = m2;
    m2 = i;
    vk = im;
  }
  i1 = IdxGet(ing, m1);
  if (i1 < 0) return;
  i2 = IdxGet(ing2, m2-m1);
  if (i2 < 0) return;
  i = i1*ing2->n + i2;
  H22Term(meff, c0, c1, ns, bra, ket, sbra, sket, mst, bst, kst,
	  s, ph, ks1, ks2, fm, a, i);
}
    
void DeltaH22M2(MBPT_EFF **meff, int ns,
		SHELL *bra, SHELL *ket, SHELL_STATE *sbra, SHELL_STATE *sket,
		int mst, int *bst, int *kst,
		CONFIG *c0, CONFIG *c1, 
		IDXARY *ib0, int *ib0s, int *ib0d, IDXARY *ib1, IDXARY *ing, 
		IDXARY *ing2, int nc, CONFIG **cs) {
  int ia, ib, ic, id, ik, im;
  int m1, m2, i, k, j1, j2;
  int op[4], om[4], ph;
  double *a[MKK*MKK];
  ORBITAL *o;
  FORMULA fm;

  fm.ns = -1;
  if (ib1->n <= 0) return;
  k = MKK*MKK;
  for (i = 0; i < k; i++) {
    a[i] = malloc(sizeof(double)*mst);
  }
  fm.js[0] = 0;
  int nmb, nmk;
  for (ia = 0; ia < ib0->n; ia++) {    
    for (ib = 0; ib <= ia; ib++) {
      if (ia == ib) {
	nmb = ib0d[ia];
      } else {
	nmb = Min(ib0s[ia], ib0s[ib]);
      }
      if (nmb <= 0) continue;
      for (ic = 0; ic < ib0->n; ic++) {
	for (id = 0; id <= ic; id++) {
	  if (ic == id) {
	    nmk = ib0d[ic];
	  } else {
	    nmk = Min(ib0s[ic], ib0s[id]);
	  }
	  if (nmk <= 0) continue;
	  om[0] = id;
	  om[1] = ic;
	  k = CheckConfig(ns-2, ket+2, 0, op, 2, om, 0, NULL);
	  if (k >= 0) continue;
	  om[0] = ia;
	  om[1] = ib;
	  k = CheckConfig(ns-2, bra+2, 0, op, 2, om, 0, NULL);
	  if (k >= 0) continue;
	  op[0] = ib;
	  op[1] = ia;
	  om[0] = ic;
	  om[1] = id;
	  ph = CheckInteraction(ns-2, bra+2, ket+2, 2, op, 2, om);
	  if (ph < 0) continue;	  
	  if (ing2->m0 == 0) {
	    fm.j1 = -1;
	    fm.j2 = -1;
	    for (m1 = 0; m1 < mbptjp.nj; m1++) {
	      if (SkipMPI()) continue;
	      for (im = mbptjp.jp[m1]; im < mbptjp.jp[m1+1]; im++) {
		ik = im;
		o = GetOrbital(ib1->d[im]);
		if (o->n > nmb || o->n > nmk) continue;
		ket[1].n = o->n;
		ket[1].kappa = o->kappa;
		ket[1].nq = 0;
		if (ket[1].n <= cs[nc]->n_csfs) {
		  om[0] = id+1;
		  om[1] = ic+1;
		  op[0] = 0;
		  op[1] = 0;
		  k = CheckConfig(ns-1, ket+1, 2, op, 2, om, nc, cs);
		  if (k >= 0) continue;
		}
		DeltaH22M2Loop(meff, c0, c1, ns, bra, ket, sbra, sket, 
			       mst, bst, kst,
			       ib0, ib1, ing, ing2, ph,
			       ia, ib, ic, id, ik, im, &fm, a);
	      }
	    }
	  }	  
	  fm.j1 = -1;
	  fm.j2 = -1;
	  for (m1 = 0; m1 < mbptjp.nj; m1++) {
	    for (m2 = 0; m2 <= m1; m2++) {
	      if (SkipMPI()) continue;
	      for (im = mbptjp.jp[m1]; im < mbptjp.jp[m1+1]; im++) {
		o = GetOrbital(ib1->d[im]);
		if (o->n > nmb || o->n > nmk) continue;
		ket[0].n = o->n;		
		ket[0].kappa = o->kappa;
		ket[0].nq = 0;
		for (ik = mbptjp.jp[m2]; ik < mbptjp.jp[m2+1]; ik++) {
		  if (im <= ik) continue;
		  o = GetOrbital(ib1->d[ik]);
		  ket[1].n = o->n;
		  ket[1].kappa = o->kappa;
		  ket[1].nq = 0;
		  if (Max(ket[0].n, ket[1].n) <= cs[nc]->n_csfs) {
		    om[0] = id+2;
		    om[1] = ic+2;
		    op[0] = 0;
		    op[1] = 1;
		    k = CheckConfig(ns, ket, 2, op, 2, om, nc, cs);
		    if (k >= 0) continue;
		  }
		  DeltaH22M2Loop(meff, c0, c1, ns, bra, ket, sbra, sket, 
				 mst, bst, kst,
				 ib0, ib1, ing, ing2, ph,
				 ia, ib, ic, id, ik, im, &fm, a);
		}
	      }	    
	    }
	  }
	}
      }
    }
  }
  k = MKK*MKK;
  for (i = 0; i < k; i++) {
    free(a[i]);
  }
}

void DeltaH22M1(MBPT_EFF **meff, int ns,
		SHELL *bra, SHELL *ket, SHELL_STATE *sbra, SHELL_STATE *sket,
		int mst, int *bst, int *kst, CONFIG *c0, CONFIG *c1, 
		IDXARY *ib0, int *ib0s, int *ib0d, IDXARY *ib1, IDXARY *ing,
		int nc, CONFIG **cs) {
  int ia, ib, ic, id, ik, im, ip, iq;
  int op[4], om[4], ph, k, ks1[4], ks2[4];
  int i, m1, i1, ij;
  double *a[MKK*MKK];
  INTERACT_SHELL s[8];
  ORBITAL *o[8];
  FORMULA fm;

  fm.ns = -1;
  fm.js[0] = 0;
  if (ib1->n <= 0) return;
  k = MKK*MKK;
  for (i = 0; i < k; i++) {
    a[i] = malloc(sizeof(double)*mst);
  }
  int nmb, nmk;
  for (ia = 0; ia < ib0->n; ia++) {
    for (ib = 0; ib <= ia; ib++) {
      for (ic = 0; ic < ib0->n; ic++) {
	for (id = 0; id <= ic; id++) {
	  for (ik = 0; ik < ib0->n; ik++) {
	    if (ia == ib) {
	      if (ik == ia) {
		nmb = ib0s[ia];
	      } else {
		nmb = ib0d[ia];
	      }
	    } else {
	      if (ik == ia) {
		nmb = ib0s[ib];
	      } else if (ik == ib) {
		nmb = ib0s[ia];
	      } else {
		nmb = Min(ib0s[ia], ib0s[ib]);
	      } 
	    }
	    if (nmb <= 0) continue;	      
	    for (iq = 0; iq < ib0->n; iq++) {
	      if (ic == id) {
		if (iq == ic) {
		  nmk = ib0s[ic];
		} else {
		  nmk = ib0d[ic];
		}
	      } else {
		if (iq == ic) {
		  nmk = ib0s[id];
		} else if (iq == id) {
		  nmk = ib0s[ic];
		} else {
		  nmk = Min(ib0s[ic], ib0s[id]);
		} 
	      }
	      if (nmk <= 0) continue;
	      op[0] = iq;
	      om[0] = id;
	      om[1] = ic;
	      k = CheckConfig(ns-1, ket+1, 1, op, 2, om, 0, NULL);
	      if (k >= 0) continue;
	      op[0] = ik;
	      om[0] = ia;
	      om[1] = ib;
	      k = CheckConfig(ns-1, bra+1, 1, op, 2, om, 0, NULL);
	      if (k >= 0) continue;
	      op[0] = ia;
	      op[1] = ib;
	      op[2] = iq;
	      om[0] = ic;
	      om[1] = id;
	      om[2] = ik;	      
	      ph = CheckInteraction(ns-1, bra+1, ket+1, 3, op, 3, om);
	      if (ph < 0) continue;
	      fm.j1 = -1;
	      for (ij = 0; ij < mbptjp.nj; ij++) {
		if (SkipMPI()) continue;
		for (ip = mbptjp.jp[ij]; ip < mbptjp.jp[ij+1]; ip++) {
		  o[1] = GetOrbital(ib1->d[ip]);
		  if (o[1]->n > nmb || o[1]->n > nmk) continue;
		  ket[0].n = o[1]->n;
		  ket[0].kappa = o[1]->kappa;
		  ket[0].nq = 0;
		  if (ket[0].n <= cs[nc]->n_csfs) {
		    om[0] = id+1;
		    om[1] = ic+1;
		    op[0] = iq+1;
		    op[1] = 0;
		    k = CheckConfig(ns, ket, 2, op, 2, om, nc, cs);
		    if (k >= 0) continue;
		  }
		  im = ip;
		  ks1[0] = ib0->d[ia];
		  ks1[1] = ib0->d[ib];
		  ks1[2] = ib1->d[im];
		  ks1[3] = ib0->d[ik];
		  ks2[0] = ib0->d[iq];
		  ks2[1] = ib1->d[ip];
		  ks2[2] = ib0->d[ic];
		  ks2[3] = ib0->d[id];
		  s[0].index = ia+1;
		  s[1].index = 0;
		  s[2].index = ib+1;
		  s[3].index = ik+1;
		  s[4].index = iq+1;
		  s[5].index = ic+1;
		  s[6].index = 0;
		  s[7].index = id+1;
		  i1 = 0;
		  for (i = 0; i < 8; i++) {
		    if (i == 1 || i == 6) {
		      o[i] = GetOrbital(ib1->d[ip]);
		    if (i == 1) {
		      m1 = o[1]->n;
		      i1 = IdxGet(ing, m1);
		      if (i1 < 0) break;
		    }
		    } else {
		      o[i] = GetOrbital(ib0->d[s[i].index-1]);
		    }		  
		    if (IsEven(i)) s[i].n = 1;
		    else s[i].n = -1;
		    s[i].kappa = o[i]->kappa;
		    GetJLFromKappa(s[i].kappa, &(s[i].j), &(s[i].kl));
		    s[i].nq_bra = bra[s[i].index].nq;
		    s[i].nq_ket = ket[s[i].index].nq;
		    s[i].index = ns-s[i].index-1;
		  }
		  if (i1 < 0) continue;
		  H22Term(meff, c0, c1, ns, bra, ket, sbra, sket, 
			  mst, bst, kst, s, ph, 
			  ks1, ks2, &fm, a, -(i1+1));
		}
	      }
	    }
	  }
	}
      }
    }
  }
  k = MKK*MKK;
  for (i = 0; i < k; i++) {
    free(a[i]);
  }
}

void DeltaH22M0(MBPT_EFF **meff, int ns,
		SHELL *bra, SHELL *ket, SHELL_STATE *sbra, SHELL_STATE *sket,
		int mst, int *bst, int *kst, CONFIG *c0, CONFIG *c1, 
		IDXARY *ib0, int *ib0s, int *ib0d,
		IDXARY *ing, int nc, CONFIG **cs) {
  int ia, ib, ic, id, ik, im, ip, iq, i, i1;
  int op[4], om[4], ph, k, ks1[4], ks2[4];
  double *a[MKK*MKK];
  INTERACT_SHELL s[8];
  ORBITAL *o[8];  
  FORMULA fm;

  fm.ns = -1;
  i1 = bra[0].n;
  i1 = IdxGet(ing, i1);
  if (i1 < 0) return;
  k = MKK*MKK;
  for (i = 0; i < k; i++) {
    a[i] = malloc(sizeof(double)*mst);
  }
  fm.js[0] = 0;
  int nmb, nmk;
  for (ia = 0; ia < ib0->n; ia++) {
    for (ib = 0; ib <= ia; ib++) {
      for (ic = 0; ic < ib0->n; ic++) {
	for (id = 0; id <= ic; id++) {
	  for (ik = 0; ik < ib0->n; ik++) {
	    for (im = 0; im <= ik; im++) {
	      if (ia == ib) {
		if (ik == ia || im == ia) {
		  nmb = ib0s[ia];
		} else {
		  nmb = ib0d[ia];
		}
	      } else {
		if (ik == ia || im == ia) {
		  nmb = ib0s[ib];
		} else if (ik == ib || im == ib) {
		  nmb = ib0s[ia];
		} else {
		  nmb = Min(ib0s[ia], ib0s[ib]);
		}
	      }
	      if (nmb <= 0) continue;
	      for (ip = 0; ip < ib0->n; ip++) {
		for (iq = 0; iq <= ip; iq++) {
		  if (ic == id) {
		    if (ip == ic || iq == ic) {
		      nmk = ib0s[ic];
		    } else {
		      nmk = ib0d[ic];
		    }
		  } else {
		    if (ip == ic || iq == ic) {
		      nmk = ib0s[id];
		    } else if (ip == id || iq == id) {
		      nmk = ib0s[ic];
		    } else {
		      nmk = Min(ib0s[ic], ib0s[id]);
		    }
		  }
		  if (nmk <= 0) continue;
		  op[0] = ip;
		  op[1] = iq;
		  om[0] = ic;
		  om[1] = id;
		  k = CheckConfig(ns, ket, 2, op, 2, om, nc, cs);
		  if (k >= 0) continue;
		  op[0] = ik;
		  op[1] = im;
		  om[0] = ia;
		  om[1] = ib;
		  k = CheckConfig(ns, bra, 2, op, 2, om, 0, NULL);
		  if (k >= 0) continue;
		  op[0] = ia;
		  op[1] = ib;
		  op[2] = ip;
		  op[3] = iq;
		  om[0] = ic;
		  om[1] = id;
		  om[2] = ik;
		  om[3] = im;
		  ph = CheckInteraction(ns, bra, ket, 4, op, 4, om);
		  if (ph < 0) continue;
		  ks1[0] = ib0->d[ia];
		  ks1[1] = ib0->d[ib];
		  ks1[2] = ib0->d[ik];
		  ks1[3] = ib0->d[im];
		  ks2[0] = ib0->d[ip];
		  ks2[1] = ib0->d[iq];
		  ks2[2] = ib0->d[ic];
		  ks2[3] = ib0->d[id];
		  s[0].index = ia;
		  s[1].index = ik;
		  s[2].index = ib;
		  s[3].index = im;
		  s[4].index = ip;
		  s[5].index = ic;
		  s[6].index = iq;
		  s[7].index = id;
		  if (ib == ik) {
		    if (im != ik) {
		      i = ks1[2];
		      ks1[2] = ks1[3];
		      ks1[3] = i;
		      i = s[1].index;
		      s[1].index = s[3].index;
		      s[3].index = i;
		    } else {
		      i = ks1[0];
		      ks1[0] = ks1[1];
		      ks1[1] = i;
		      i = s[0].index;
		      s[0].index = s[2].index;
		      s[2].index = i;
		    }
		  }
		  if (iq == ic) {
		    if (ip != iq) {
		      i = ks2[0];
		      ks2[0] = ks2[1];
		      ks2[1] = i;
		      i = s[4].index;
		      s[4].index = s[6].index;
		      s[6].index = i;
		    } else {
		      i = ks2[2];
		      ks2[2] = ks2[3];
		      ks2[3] = i;
		      i = s[5].index;
		      s[5].index = s[7].index;
		      s[7].index = i;
		    }
		  }
		  for (i = 0; i < 8; i++) {
		    o[i] = GetOrbital(ib0->d[s[i].index]);
		    if (IsEven(i)) s[i].n = 1;
		    else s[i].n = -1;
		    s[i].kappa = o[i]->kappa;
		    GetJLFromKappa(s[i].kappa, &(s[i].j), &(s[i].kl));
		    s[i].nq_bra = bra[s[i].index].nq;
		    s[i].nq_ket = ket[s[i].index].nq;
		    s[i].index = ns-s[i].index-1;
		  }
		  fm.j1 = -1;
		  if (SkipMPI()) continue;
		  H22Term(meff, c0, c1, ns, bra, ket, sbra, sket, 
			  mst, bst, kst, s, ph,
			  ks1, ks2, &fm, a, -(i1+1));
		}
	      }
	    }
	  }
	}
      }
    }
  }
  k = MKK*MKK;
  for (i = 0; i < k; i++) {
    free(a[i]);
  }
}

void DeltaH12M1(void *mptr, int ns,
		SHELL *bra, SHELL *ket, SHELL_STATE *sbra, SHELL_STATE *sket,
		int mst, int *bst, int *kst, CONFIG *c0, CONFIG *c1, 
		IDXARY *ib0, int *ib0s, int *ib0d,
		IDXARY *ib1, IDXARY *ing,
		int nc, CONFIG **cs, int mode) {
  int ia, ib, ic, id, ik, im, ij;
  int m1, i1, k, i, op[3], om[3], ks[4], ph;
  double *a[MKK];
  INTERACT_SHELL s[6];
  ORBITAL *o[6];
  FORMULA fm;

  fm.ns = -1;
  fm.js[0] = 0;
  if (ib1->n <= 0) return;
  k = MKK;
  for (i = 0; i < k; i++) {
    if (mode == 0) {
      a[i] = malloc(sizeof(double)*mst);
    } else {
      a[i] = malloc(sizeof(double)*mst*mbpt_tr.mktr);
    }
  }
  int nmb, nmk;
  for (ia = 0; ia < ib0->n; ia++) {
    if (bra[ia+1].nq == 0) continue;
    nmb = ib0s[ia];
    if (nmb <= 0) continue;
    for (ib = 0; ib < ib0->n; ib++) {
      for (ic = 0; ic <= ib; ic++) {
	for (id = 0; id < ib0->n; id++) {
	  if (ib == ic) {
	    if (id == ib) {
	      nmk = ib0s[ib];
	    } else {
	      nmk = ib0d[ib];
	    }
	  } else {
	    if (id == ib) {
	      nmk = ib0s[ic];
	    } else if (id == ic) {
	      nmk = ib0s[ib];
	    } else {
	      nmk = Min(ib0s[ib], ib0s[ic]);
	    }
	  }
	  if (nmk <= 0) continue;
	  op[0] = id;
	  om[0] = ic;
	  om[1] = ib;
	  k = CheckConfig(ns-1, ket+1, 1, op, 2, om, 0, NULL);
	  if (k >= 0) continue;
	  op[0] = ia;
	  op[1] = id;
	  om[0] = ic;
	  om[1] = ib;
	  ph = CheckInteraction(ns-1, bra+1, ket+1, 2, op, 2, om);
	  if (ph < 0) continue;
	  fm.j1 = -1;
	  for (ij = 0; ij < mbptjp.nj; ij++) {	  
	    if (SkipMPI()) continue;
	    for (ik = mbptjp.jp[ij]; ik < mbptjp.jp[ij+1]; ik++) {
	      o[1] = GetOrbital(ib1->d[ik]);
	      if (o[1]->n > nmb || o[1]->n > nmk) continue;
	      ket[0].n = o[1]->n;
	      ket[0].kappa = o[1]->kappa;
	      ket[0].nq = 0;
	      if (ket[0].n <= cs[nc]->n_csfs) {
		om[0] = ic+1;
		om[1] = ib+1;
		op[0] = id+1;
		op[1] = 0;
		k = CheckConfig(ns, ket, 2, op, 2, om, nc, cs);
		if (k >= 0) continue;
	      }
	      im = ik;
	      ks[0] = ib0->d[id];
	      ks[1] = ib1->d[im];
	      ks[2] = ib0->d[ib];
	      ks[3] = ib0->d[ic];
	      s[0].index = ia+1;
	      s[1].index = 0;
	      s[2].index = id+1;
	      s[3].index = ib+1;
	      s[4].index = 0;
	      s[5].index = ic+1;
	      i1 = 0;
	      for (i = 0; i < 6; i++) {
		if (i == 1 || i == 4) {		
		  o[i] = GetOrbital(ib1->d[ik]);
		  if (i == 1) {
		    m1 = o[1]->n;
		    i1 = IdxGet(ing, m1);
		    if (i1 < 0) break;
		  }
		} else {
		  o[i] = GetOrbital(ib0->d[s[i].index-1]);
		}
		if (IsEven(i)) s[i].n = 1;
		else s[i].n = -1;
		s[i].kappa = o[i]->kappa;
		GetJLFromKappa(s[i].kappa, &(s[i].j), &(s[i].kl));
		s[i].nq_bra = bra[s[i].index].nq;
		s[i].nq_ket = ket[s[i].index].nq;
		s[i].index = ns-s[i].index-1;
	      }
	      if (i1 < 0) continue;
	      if (mode == 0) {
		H12Term((MBPT_EFF **) mptr, c0, c1, ns, bra, ket, sbra, sket, 
			mst, bst, kst, s, ph, 
			ks, ib0->d[ia], ib1->d[ik], &fm, a, i1);
	      } else {
		TR12Term((MBPT_TR *) mptr, c0, c1, ns, bra, ket, sbra, sket, 
			 mst, bst, kst, s, ph, 
			 ks, ib0->d[ia], ib1->d[ik], &fm, a, i1, ing->n);
	      }
	    }
	  }
	}
      }
    }
  }
  k = MKK;
  for (i = 0; i < k; i++) {
    free(a[i]);
  }
}

void DeltaH12M0(void *mptr, int ns,
		SHELL *bra, SHELL *ket, SHELL_STATE *sbra, SHELL_STATE *sket,
		int mst, int *bst, int *kst, CONFIG *c0, CONFIG *c1, 
		IDXARY *ib0, int *ib0s, int *ib0d,
		IDXARY *ing, int nc, CONFIG **cs, int mode) {
  int ia, ib, ic, id, ik, im;
  int i1, k, i, op[3], om[3], ks[4], ph;
  double *a[MKK];
  INTERACT_SHELL s[6];
  ORBITAL *o[6];
  FORMULA fm;

  fm.ns = -1;
  i1 = bra[0].n;
  i1 = IdxGet(ing, i1);
  if (i1 < 0) return;
  k = MKK;
  for (i = 0; i < k; i++) {
    if (mode == 0) {
      a[i] = malloc(sizeof(double)*mst);
    } else {
      a[i] = malloc(sizeof(double)*mst*mbpt_tr.mktr);
    }
  }
  fm.js[0] = 0;
  int nmb, nmk;
  for (ia = 0; ia < ib0->n; ia++) {
    if (bra[ia].nq == 0) continue;
    for (ib = 0; ib < ib0->n; ib++) {
      if (ket[ib].nq == 0) continue;      
      for (ic = 0; ic <= ib; ic++) {
	if (ket[ic].nq == 0) continue;
	for (ik = 0; ik < ib0->n; ik++) {
	  if (ik == ia) continue;
	  nmb = ib0s[ia];
	  if (nmb <= 0) continue;
	  for (id = 0; id < ib0->n; id++) {
	    for (im = 0; im <= id; im++) {
	      if (ib == ic) {
		if (id == ib || im == ib) {
		  nmk = ib0s[ib];
		} else {
		  nmk = ib0d[ib];
		}
	      } else {
		if (id == ib || im == ib) {
		  nmk = ib0s[ic];
		} else if (id == ic || im == ic) {
		  nmk = ib0s[ib];
		} else {
		  nmk = Min(ib0s[ib], ib0s[ic]);
		}
	      }
	      if (nmk <= 0) continue;
	      if (SkipMPI()) continue;
	      op[0] = ia;
	      op[1] = im;
	      op[2] = id;
	      om[0] = ic;
	      om[1] = ib;
	      om[2] = ik;
	      ph = CheckInteraction(ns, bra, ket, 3, op, 3, om);
	      if (ph < 0) continue;
	      op[0] = im;
	      op[1] = id;
	      om[0] = ic;
	      om[1] = ib;
	      k = CheckConfig(ns, ket, 2, op, 2, om, nc, cs);
	      if (k >= 0) continue;
	      op[0] = ik;
	      om[0] = ia;
	      k = CheckConfig(ns, bra, 1, op, 1, om, 0, NULL);
	      if (k >= 0) continue;
	      ks[0] = ib0->d[id];
	      ks[1] = ib0->d[im];
	      ks[2] = ib0->d[ib];
	      ks[3] = ib0->d[ic];
	      s[0].index = ia;
	      s[1].index = ik;
	      s[2].index = id;
	      s[3].index = ib;
	      s[4].index = im;
	      s[5].index = ic;
	      /* make sure that in Z(a,b)\dot Z(c,d), b!=c, by rearrange orbs */
	      if (im == ib) {
		if (id != im) {
		  i = ks[0];
		  ks[0] = ks[1];
		  ks[1] = i;
		  i = s[2].index;
		  s[2].index = s[4].index;
		  s[4].index = i;
		} else {
		  i = ks[2];
		  ks[2] = ks[3];
		  ks[3] = i;
		  i = s[3].index;
		  s[3].index = s[5].index;
		  s[5].index = i;
		}
	      }
	      for (i = 0; i < 6; i++) {
		o[i] = GetOrbital(ib0->d[s[i].index]);
		if (IsEven(i)) s[i].n = 1;
		else s[i].n = -1;
		s[i].kappa = o[i]->kappa;
		GetJLFromKappa(s[i].kappa, &(s[i].j), &(s[i].kl));
		s[i].nq_bra = bra[s[i].index].nq;
		s[i].nq_ket = ket[s[i].index].nq;
		s[i].index = ns-s[i].index-1;
	      }
	      fm.j1 = -1;
	      if (mode == 0) {
		H12Term((MBPT_EFF **) mptr, c0, c1, ns, bra, ket, sbra, sket,
			mst, bst, kst, s, ph,
			ks, ib0->d[ia], ib0->d[ik], &fm, a, i1);
	      } else {
		TR12Term((MBPT_TR *) mptr, c0, c1, ns, bra, ket, sbra, sket,
			 mst, bst, kst, s, ph,
			 ks, ib0->d[ia], ib0->d[ik], &fm, a, i1, ing->n);
	      }
	    }
	  }
	}
      }
    }
  }
  k = MKK;
  for (i = 0; i < k; i++) {
    free(a[i]);
  }
}

void DeltaH11M1(void *mptr, int ns, 
		SHELL *bra, SHELL *ket, SHELL_STATE *sbra, SHELL_STATE *sket,
		int mst, int *bst, int *kst, CONFIG *c0, CONFIG *c1, 
		IDXARY *ib0, int *ib0s, IDXARY *ib1, IDXARY *ing,
		int nc, CONFIG **cs, int mode) {
  int ia, ib, ik, im, k, k0, k1, k2, k3, i, i1, m1, ij;
  int op[2], om[2], ph;
  INTERACT_SHELL s[4];
  ORBITAL *o0, *o1, *o2, *o3;
  double *a;
  FORMULA fm;

  fm.ns = -1;
  fm.js[0] = 0;
  if (ib1->n <= 0) return;
  if (mode == 0) {
    a = malloc(sizeof(double)*mst);
  } else {
    a = malloc(sizeof(double)*mst*mbpt_tr.mktr);
  }
  int nmb, nmk;
  for (ia = 0; ia < ib0->n; ia++) {
    if (bra[ia+1].nq == 0) continue;
    nmb = ib0s[ia];
    if (nmb <= 0) continue;
    for (ib = 0; ib < ib0->n; ib++) { 
      if (ket[ib+1].nq == 0) continue;
      nmk = ib0s[ib];
      if (nmk <= 0) continue;
      op[0] = ia;
      om[0] = ib;
      /* check if bra and ket can interact with 1 virtual orb */
      ph = CheckInteraction(ns-1, bra+1, ket+1, 1, op, 1, om);
      if (ph < 0) continue;
      fm.j1 = -1;
      for (ij = 0; ij < mbptjp.nj; ij++) {
	if (SkipMPI()) continue;
	for (ik = mbptjp.jp[ij]; ik < mbptjp.jp[ij+1]; ik++) {
	  o1 = GetOrbital(ib1->d[ik]);
	  if (o1->n > nmb || o1->n > nmk) continue;
	  ket[0].n = o1->n;
	  ket[0].kappa = o1->kappa;
	  ket[0].nq = 0;
	  /* if n < maxn(real), check if ket(ib+1->ik) is in model space */
	  if (ket[0].n <= cs[nc]->n_csfs) {
	    om[0] = ib+1;
	    op[0] = 0;
	    k = CheckConfig(ns, ket, 1, op, 1, om, nc, cs);
	    if (k >= 0) continue;
	  }
	  im = ik;
	  k0 = ib0->d[ia];
	  k1 = ib1->d[ik];
	  k2 = ib1->d[im];
	  k3 = ib0->d[ib];
	  s[0].index = ia+1;
	  s[0].n = 1;
	  s[1].index = 0;
	  s[1].n = -1;
	  s[2].index = 0;
	  s[2].n = 1;
	  s[3].index = ib+1;
	  s[3].n = -1;	  
	  s[0].nq_bra = bra[s[0].index].nq;
	  s[0].nq_ket = ket[s[0].index].nq;
	  s[1].nq_bra = s[1].nq_ket = 0;
	  s[2].nq_bra = s[1].nq_ket = 0;
	  s[3].nq_bra = bra[s[3].index].nq;
	  s[3].nq_ket = ket[s[3].index].nq;
	  o0 = GetOrbital(k0);
	  o1 = GetOrbital(k1);
	  m1 = o1->n;	
	  i1 = IdxGet(ing, m1);
	  if (i1 < 0) continue;
	  o2 = GetOrbital(k2);
	  o3 = GetOrbital(k3);
	  s[0].kappa = o0->kappa;
	  s[1].kappa = o1->kappa;
	  s[2].kappa = o2->kappa;
	  s[3].kappa = o3->kappa;
	  for (i = 0; i < 4; i++) {
	    GetJLFromKappa(s[i].kappa, &(s[i].j), &(s[i].kl));
	    s[i].index = ns-s[i].index-1;
	  }
	  if (mode == 0) {
	    H11Term((MBPT_EFF **) mptr, c0, c1, ns, bra, ket, 
		    sbra, sket, mst, bst, kst, 
		    s, ph, k0, k1, k2, k3, &fm, a, i1);
	  } else {
	    TR11Term((MBPT_TR *) mptr, c0, c1, ns, bra, ket,
		     sbra, sket, mst, bst, kst, 
		     s, ph, k0, k1, k2, k3, &fm, a, i1, ing->n);
	  }
	}
      }
    }
  }
  free(a);
}
  
void DeltaH11M0(void *mptr, int ns, 
		SHELL *bra, SHELL *ket, SHELL_STATE *sbra, SHELL_STATE *sket,
		int mst, int *bst, int *kst,
		CONFIG *c0, CONFIG *c1,
		IDXARY *ib0, int *ib0s, IDXARY *ing,
		int nc, CONFIG **cs, int mode) {
  int ia, ib, ik, im, k, k0, k1, k2, k3, i, i1;
  int op[2], om[2], ph;
  INTERACT_SHELL s[4];
  ORBITAL *o0, *o1, *o2, *o3;
  double *a;
  FORMULA fm;

  fm.ns = -1;
  i1 = bra[0].n;
  i1 = IdxGet(ing, i1);
  if (i1 < 0) return;
  if (mode == 0) {
    a = malloc(sizeof(double)*mst);
  } else {
    a = malloc(sizeof(double)*mst*mbpt_tr.mktr);
  }
  fm.js[0] = 0;
  int nmb, nmk;
  for (ia = 0; ia < ib0->n; ia++) {
    if (bra[ia].nq == 0) continue;
    if (ib0s[ia] <= 0) continue;
    for (ib = 0; ib < ib0->n; ib++) { 
      if (ket[ib].nq == 0) continue;
      if (ib0s[ib] <= 0) continue;
      for (ik = 0; ik < ib0->n; ik++) {
	if (ik == ia) {
	  nmb = 1000000;
	} else {
	  nmb = ib0s[ia];
	}
	if (nmb <= 0) continue;
	for (im = 0; im < ib0->n; im++) {
	  if (ia == ik) {
	    continue;
	  }
	  if (ib == im) {
	    continue;
	  }
	  nmk = ib0s[ib];
	  if (nmk <= 0) continue;
	  if (SkipMPI()) continue;
	  /* op contains the index for the creation operators */
	  op[0] = ia;
	  op[1] = im;
	  /* om contains the index for the annihilation operators */
	  om[0] = ib;
	  om[1] = ik;
	  /* check if bra and ket states can interact with op om operators */
	  ph = CheckInteraction(ns, bra, ket, 2, op, 2, om);
	  /* ph returns the phase factor (-1)^ph for recoupling */
	  if (ph < 0) continue;
	  op[0] = im;
	  om[0] = ib;
	  /* check if ket(ib->im) is within the model space */
	  k = CheckConfig(ns, ket, 1, op, 1, om, nc, cs);	  
	  if (k >= 0) continue;
	  op[0] = ik;
	  om[0] = ia;	  
	  /* check if bra(ia->ik) is within the model space */
	  k = CheckConfig(ns, bra, 1, op, 1, om, 0, NULL);
	  if (k >= 0) continue;
	  k0 = ib0->d[ia];
	  k1 = ib0->d[ik];
	  k2 = ib0->d[im];
	  k3 = ib0->d[ib];
	  s[0].index = ia;
	  s[0].n = 1;
	  s[1].index = ik;
	  s[1].n = -1;
	  s[2].index = im;
	  s[2].n = 1;
	  s[3].index = ib;
	  s[3].n = -1;
	  for (i = 0; i < 4; i++) {
	    s[i].nq_bra = bra[s[i].index].nq;
	    s[i].nq_ket = ket[s[i].index].nq;
	  }
	  o0 = GetOrbital(k0);
	  o1 = GetOrbital(k1);
	  o2 = GetOrbital(k2);
	  o3 = GetOrbital(k3);
	  s[0].kappa = o0->kappa;
	  s[1].kappa = o1->kappa;
	  s[2].kappa = o2->kappa;
	  s[3].kappa = o3->kappa;
	  for (i = 0; i < 4; i++) {
	    GetJLFromKappa(s[i].kappa, &(s[i].j), &(s[i].kl));
	    s[i].index = ns-s[i].index-1;
	  }
	  fm.j1 = -1;
	  if (mode == 0) {
	    H11Term((MBPT_EFF **) mptr, c0, c1, ns, bra, ket, 
		    sbra, sket, mst, bst, kst,
		    s, ph, k0, k1, k2, k3, &fm, a, i1);
	  } else {
	    TR11Term((MBPT_TR *) mptr, c0, c1, ns, bra, ket,
		     sbra, sket, mst, bst, kst,
		     s, ph, k0, k1, k2, k3, &fm, a, i1, ing->n);
	  }
	}
      }
    }
  }
  free(a);
}
  
void FreeEffMBPT(MBPT_EFF **meff) {
  int i, j, k, m;

  for (i = 0; i < MAX_SYMMETRIES; i++) {
    if (meff[i] == NULL) continue;
    if (meff[i]->nbasis > 0) {
      k = meff[i]->nbasis;
      k = k*(k+1)/2;
      for (m = 0; m < k; m++) {
	if (meff[i]->hab1[m]) {
	  free(meff[i]->hab1[m]);
	  free(meff[i]->hba1[m]);
	  free(meff[i]->hab[m]);
	  free(meff[i]->hba[m]);
	}
      }
      free(meff[i]->hab1);
      free(meff[i]->hba1);
      free(meff[i]->hab);
      free(meff[i]->hba);
      free(meff[i]->h0);
      free(meff[i]->e0);
      free(meff[i]->heff);
      free(meff[i]->basis);
    }
    free(meff[i]);
  }
}

void InitTransitionMBPT(MBPT_TR **mtr0, int n) {
  MBPT_TR *mtr;
  int i, j, k, m, isym, k0, k1, m0, m1, ms0, mst;

  mtr = malloc(sizeof(MBPT_TR)*(2*mbpt_tr.mktr)*MAX_SYMMETRIES);
  *mtr0 = mtr;
  j = 0;
  for (isym = 0; isym < MAX_SYMMETRIES; isym++) {
    DecodePJ(isym, &k1, &k0);
    for (i = 1; i <= mbpt_tr.mktr; i++) {
      for (k = -1; k <= 1; k += 2) {
	mtr[j].m = k*i;
	mtr[j].sym0 = GetSymmetry(isym);
	mtr[j].nsym1 = 0;
	if (mtr[j].sym0->n_states > 0) {
	  mtr[j].nsym1 = 2*i + 1;
	  mtr[j].isym1 = malloc(sizeof(int)*mtr[j].nsym1);
	  mtr[j].sym1 = malloc(sizeof(SYMMETRY *)*mtr[j].nsym1);
	  mtr[j].tma = malloc(sizeof(double *)*mtr[j].nsym1);
	  mtr[j].rma = malloc(sizeof(double *)*mtr[j].nsym1);
	  m = 0;
	  for (m0 = k0 - 2*i; m0 <= k0 + 2*i; m0 += 2) {
	    if (m0 < 0) continue;
	    m1 = k1;
	    if (k < 0 && IsOdd(i)) m1 = 1-m1;
	    if (k > 0 && IsEven(i)) m1 = 1-m1;
	    mtr[j].isym1[m] = IsEven(m1)? 2*m0 : (2*m0+1);
	    if (mtr[j].isym1[m] >= MAX_SYMMETRIES) continue;
	    mtr[j].sym1[m] = GetSymmetry(mtr[j].isym1[m]);
	    mst = mtr[j].sym0->n_states * mtr[j].sym1[m]->n_states;
	    mst *=  n * mbpt_tr.naw;
	    if (mst > 0) {
	      mtr[j].tma[m] = malloc(sizeof(double)*mst);
	      mtr[j].rma[m] = malloc(sizeof(double)*mst);
	      for (ms0 = 0; ms0 < mst; ms0++) {
		mtr[j].tma[m][ms0] = 0.0;
		mtr[j].rma[m][ms0] = 0.0;
	      }
	      m++;
	    }
	  }
	  mtr[j].nsym1 = m;
	  if (mtr[j].nsym1 == 0) {
	    free(mtr[j].isym1);
	    free(mtr[j].sym1);
	    free(mtr[j].tma);
	    free(mtr[j].rma);
	  }
	}
	j++;
      }
    }
  }
}

void FreeTransitionMBPT(MBPT_TR *mtr) {
  int j, k, m, n;

  if (mbpt_tr.mktr <= 0) return;
  n = 2*mbpt_tr.mktr*MAX_SYMMETRIES;
  for (j = 0; j < n; j++) {
    if (mtr[j].nsym1 == 0) continue;
    for (m = 0; m < mtr[j].nsym1; m++) {
      k = mtr[j].sym0->n_states * mtr[j].sym1[m]->n_states;
      if (k > 0) {
	free(mtr[j].tma[m]);
	free(mtr[j].rma[m]);
      }
    }
    free(mtr[j].tma);
    free(mtr[j].rma);
    free(mtr[j].isym1);
    free(mtr[j].sym1);    
  }
  free(mtr);
}

int GetJpList(int n, int *bas, int *jp) {
  ORBITAL *orb;
  int i, j1, j2, nj;
  
  jp[0] = 0;
  nj = 1;
  orb = GetOrbital(bas[0]);
  j1 = GetJFromKappa(orb->kappa);
  for (i = 1; i < n; i++) {
    orb = GetOrbital(bas[i]);
    j2 = GetJFromKappa(orb->kappa);
    if (j2 != j1) {
      jp[nj] = i;
      nj++;
      j1 = j2;
    }
  }
  jp[nj] = n;
  return nj;
}

/*
** fn, the energy file
** fn1, the effective hamilton file
** nkg, number of groups or the model space
** kg, the group indexes of the model space
** nk, nkm, specifies the maximum n-value for each 
**     orbital angular momentum of the virtual orbital.
** n, the number of virtual orbitals of the first excitation.
** ng, the n-grid of the first excitation.
** n2, the number of virtual orbitals of the 2nd excitation.
** ng2, the n-grid (n2[i,j]=ng[i]+ng2[j]) of the 2nd excitation.
** nkg0, the number of groups in kg to be included for MBPT correction.
*/
int StructureMBPT1(char *fn, char *fn1, int nkg, int *kg, int nk, int *nkm, 
		   int n, int *ng, int n2, int *ng2, int nkg0) {
  int *bas, *bas0, *bas1, nlevels, nb0, nb, n3, q, nhab, nhab1;
  int i, j, k, i0, i1, n0, n1, isym, ierr, nc, m, mks, *ks;
  int pp, jj, nmax, na, *ga, k0, k1, m0, m1, nmax1, mst, ncps;
  int p0, p1, j0, j1, j2, q0, q1, ms0, ms1, *bst, *kst, *bst0, *kst0;
  char tfn[1024];
  SYMMETRY *sym;
  STATE *st;
  HAMILTON *h;
  SHELL *bra, *ket, *bra1, *ket1, *bra2, *ket2;
  SHELL_STATE *sbra, *sket, *sbra1, *sket1, *sbra2, *sket2;
  CONFIG **cs, *c0, *c1, *ct0, *ct1;
  CONFIG_GROUP *g;
  ORBITAL *orb0, *orb1;
  MBPT_EFF *meff[MAX_SYMMETRIES];
  MBPT_TR *mtr;
  double a, b, c, *mix, *hab, *hba, emin, emax;
  double *h0, *heff, *hab1, *hba1, *dw, tt0, tt1, tbg, dt, dtt;
  FILE *f;
  IDXARY ing, ing2;
  ORBITAL *orb;

  ing.n = ing.m = 0;
  ing2.n = ing2.m = 0;
  ierr = 0;
  n3 = mbpt_n3;
  if (nkg0 <= 0 || nkg0 > nkg) nkg0 = nkg;

  /* construct configurations in kg, determine the maximum n-value*/
  nc = 0;
  for (i = 0; i < nkg; i++) {
    g = GetGroup(kg[i]);
    nc += g->n_cfgs;
  }
  int nmaxm, nmaxm1;
  CONFIG **csm;
  ResetWidMPI();
#pragma omp parallel default(shared) private(k, i1, cs, i, j, nmax, nmax1, g)
  {
    cs = malloc(sizeof(CONFIG *)*(nc+2));
    k = 0;
    i1 = 0;
    nmax = 0;
    nmax1 = 0;
    for (i = 0; i < nkg; i++) {
      g = GetGroup(kg[i]);    
      for (j = 0; j < g->n_cfgs; j++) {
	cs[k] = GetConfigFromGroup(kg[i], j);
	if (cs[k]->n_shells > i1) i1 = cs[k]->n_shells;
	if (cs[k]->shells[0].n > nmax) nmax = cs[k]->shells[0].n;
	if (cs[k]->shells[0].nq > 1) {
	  if (cs[k]->shells[0].n > nmax1) nmax1 = cs[k]->shells[0].n;
	} else if (cs[k]->n_shells > 1) {
	  if (cs[k]->shells[1].n > nmax1) nmax1 = cs[k]->shells[1].n;
	}
	k++;
      }
    }  
    /* use cfg.n_csfs to store the maximum n-value of the configuration in kg */
    mbpt_cfg.n_csfs = nmax;
    mbpt_cfg.nnrs = nmax1;
    cs[k] = &mbpt_cfg;
    cs[k]->n_shells = i1;
    cs[k]->shells = malloc(sizeof(SHELL)*(i1+2));
    mbpt_cs = cs;
#pragma omp master
    {
      nmaxm = nmax;
      nmaxm1 = nmax1;
      csm = cs;
    }
#pragma omp flush
  }
  nmax = nmaxm;
  nmax1 = nmaxm1;
  cs = csm;
  tt0 = WallTime();
  tbg = WallTime();

  MPrintf(-1, "Construct Radial Basis, %d/%d.\n", MyRankMPI(), NProcMPI());
  fflush(stdout);
  n = ConstructNGrid(n, &ng);
  n2 = ConstructNGrid(n2, &ng2);
  InitIdxAry(&ing, n, ng);
  InitIdxAry(&ing2, n2, ng2);
  nhab1 = n*2;
  nhab = n*n2*2;
  na = n*n2+n+nmax;
  ga = malloc(sizeof(int)*na);
  k = 0;
  for (i = 1; i <= nmax; i++) {
    ga[k++] = i;
  }
  nb = PrepRadialBasisMBPT(nmax, NULL, k, ga, &bas0);
  nb0 = nb;
  k = 0;
  for (i = 0; i < n; i++) {
    ga[k++] = ng[i];
    for (j = 0; j < n2; j++) {
      if (ng2[j] == 0) continue;
      ga[k++] = ng[i]+ng2[j];      
    }
  }
  na = SortUnique(k, ga);
  nb = PrepRadialBasisMBPT(nk, nkm, na, ga, &bas);
  if (mbpt_nsplit) {
    mbpt_rij = malloc(sizeof(int *)*nb);
    m = 0;
    for (i = 0; i < nb; i++) {
      mbpt_rij[i] = malloc(sizeof(int)*nb);
      orb0 = GetOrbital(i);
      i0 = IdxGet(&ing, orb0->n);      
      for (j = 0; j < nb; j++) {
	orb1 = GetOrbital(j);
	j0 = IdxGet(&ing2, orb1->n-orb0->n);
	if (i0 < 0 || j0 < 0) {
	  mbpt_rij[i][j] = -1;
	  continue;
	}
	mbpt_rij[i][j] = m;
	m++;
	if (m == NProcMPI()) m = 0;
      }
    }
  }
  ResetWidMPI();
#pragma omp parallel default(shared) private(i,j,i0,j0)
  {
    mbpt_bas0 = malloc(sizeof(int)*nb0);  
    mbpt_bas0s = malloc(sizeof(int)*nb0);
    mbpt_bas0d = malloc(sizeof(int)*nb0);  
    mbpt_bas1 = malloc(sizeof(int)*nb);
    mbptjp.jp = malloc(sizeof(int)*(nb+1));
    if (mbpt_nsplit) {
      InitIdxAry(&mbptjp.ibs, nb, bas);
      for (i = 0; i < mbptjp.ibs.m; i++) {
	mbptjp.ibs.i[i] = -1-mbptjp.ibs.i[i];
      }
      for (i = 0; i < nb; i++) {
	i0 = bas[i] - mbptjp.ibs.m0;	
	for (j = 0; j < nb; j++) {
	  j0 = bas[j] - mbptjp.ibs.m0;
	  if (mbpt_rij[i][j] == MyRankMPI()) {
	    if (mbptjp.ibs.i[i0] < 0) {
	      mbptjp.ibs.i[i0] = -1-mbptjp.ibs.i[i0];
	    }
	    if (mbptjp.ibs.i[j0] < 0) {
	      mbptjp.ibs.i[j0] = -1-mbptjp.ibs.i[j0];
	    }
	  }
	}
      }
    }
  }
  free(bas0);
  free(ga);
  SolveRadialBasisMBPT(nmax);
  
  ShiftOrbitalEnergy(cs[0]);

  tt1 = WallTime();
  dt = tt1-tt0;
  tt0 = tt1;
  MPrintf(-1, "Time = %12.5E\n", dt);
  fflush(stdout);
  if (nb < 0) return -1;

  if (n3 >= 0 && MyRankMPI() == 0) {
    f = fopen(fn1, "w");
    if (f == NULL) {
      MPrintf(-1, "cannot open file %s\n", fn1);
      FreeIdxAry(&ing, 2);
      FreeIdxAry(&ing2, 2);
      ResetWidMPI();
#pragma omp parallel
      {
	free(mbptjp.jp);
	free(mbpt_bas0);
	free(mbpt_bas0s);
	free(mbpt_bas0d);
	free(mbpt_bas1);
	if (mbpt_nsplit) {
	  FreeIdxAry(&mbptjp.ibs, 2);
	}
      }
      free(bas);
      if (mbpt_nsplit) {
	for (i = 0; i < nb; i++) {
	  free(mbpt_rij[i]);
	}
	free(mbpt_rij);
      }
      return -1;
    }
    fwrite(&n, sizeof(int), 1, f);
    fwrite(ng, sizeof(int), n, f);
    fwrite(&n2, sizeof(int), 1, f);
    fwrite(ng2, sizeof(int), n2, f);
    fwrite(&n3, sizeof(int), 1, f);
  }

  dw = malloc(sizeof(double)*(n+n2)*4);
  MPrintf(-1, "CI Structure.\n");
  fflush(stdout);
  nlevels = GetNumLevels();
  emax = -1E31;
  emin = 1E31;
  for (isym = 0; isym < MAX_SYMMETRIES; isym++) {
    meff[isym] = NULL;
    /* the construction is such that h->dim = h->n_basis
    ** one also makes sure that all the states in one symmetry forms a
    ** single hamiltonian matrix */
    h = GetHamilton(i);
    k = ConstructHamilton(isym, nkg0, nkg, kg, 0, NULL, 110);
    if (k == -1) continue;
    meff[isym] = (MBPT_EFF *) malloc(sizeof(MBPT_EFF));
    if (k < 0) {
      meff[isym]->nbasis = 0;
      continue;
    }
    double mem0=TotalSize();
    MPrintf(-1, "sym: %3d %d %g\n", isym, h->dim, mem0);
    fflush(stdout);    
    sym = GetSymmetry(isym);
    meff[isym]->h0 = malloc(sizeof(double)*h->hsize);
    meff[isym]->e0 = malloc(sizeof(double)*h->n_basis);
    memcpy(meff[isym]->h0, h->hamilton, sizeof(double)*h->hsize);
    meff[isym]->basis = malloc(sizeof(int)*h->n_basis);
    meff[isym]->nbasis = h->n_basis;
    meff[isym]->hsize = h->hsize;
    memcpy(meff[isym]->basis, h->basis, sizeof(int)*h->n_basis);
    meff[isym]->heff = malloc(sizeof(double)*h->n_basis*h->n_basis);
    meff[isym]->hab1 = malloc(sizeof(double *)*h->hsize);
    meff[isym]->hba1 = malloc(sizeof(double *)*h->hsize);
    meff[isym]->hab = malloc(sizeof(double *)*h->hsize);
    meff[isym]->hba = malloc(sizeof(double *)*h->hsize);  
    meff[isym]->n = n;
    meff[isym]->n2 = n2;
    if (DiagnolizeHamilton(h) < 0) {
      MPrintf(-1, "Diagnolizing Hamiltonian Error\n");
      fflush(stdout);
      Abort(1);
    }    
    ks = malloc(sizeof(int)*h->dim);
    mks = 0;
    m = 0;
    mix = h->mixing + h->dim;
    for (k = 0; k < h->dim; k++) {
      /* determine the states that need MBPT correction */
      if (nkg0 != nkg) {
	q = GetPrincipleBasis(mix, h->dim, NULL);
	st = (STATE *) ArrayGet(&(sym->states), h->basis[q]);
	if (InGroups(st->kgroup, nkg0, kg)) {	  
	  if (mbpt_nlev <= 0 || IBisect(m, mbpt_nlev, mbpt_ilev) >= 0) {
	    ks[mks++] = k;
	  }
	  m++;
	  if (h->mixing[k] < emin) emin = h->mixing[k];
	  if (h->mixing[k] > emax) emax = h->mixing[k];
	}
      } else {
	if (mbpt_nlev <= 0 || IBisect(m, mbpt_nlev, mbpt_ilev) >= 0) {
	  ks[mks++] = k;
	}
	m++;
	if (h->mixing[k] < emin) emin = h->mixing[k];
	if (h->mixing[k] > emax) emax = h->mixing[k];
      }
      mix += h->dim;
    }
    int ki=0, ke=0;
    for (j = 0; j < h->dim; j++) {
      for (i = 0; i <= j; i++) {
	c = 0;
	for (m = 0; m < mks; m++) {
	  k = ks[m];
	  mix = h->mixing + h->dim*(k+1);
	  a = fabs(mix[i]*mix[j]);
	  if (a > c) c = a;
	}
	k = j*(j+1)/2 + i;
	a = mbpt_mcut;
	if (i == j) a *= 0.2;
	if (c >= a) {
	  meff[isym]->hab1[k] = malloc(sizeof(double)*nhab1);
	  meff[isym]->hba1[k] = malloc(sizeof(double)*nhab1);
	  meff[isym]->hab[k] = malloc(sizeof(double)*nhab);
	  meff[isym]->hba[k] = malloc(sizeof(double)*nhab);
	  for (i0 = 0; i0 < nhab1; i0++) {
	    meff[isym]->hab1[k][i0] = 0.0;
	    meff[isym]->hba1[k][i0] = 0.0;
	  } 
	  for (i0 = 0; i0 < nhab; i0++) {
	    meff[isym]->hab[k][i0] = 0.0;
	    meff[isym]->hba[k][i0] = 0.0;
	  }
	  ki++;
	} else {
	  meff[isym]->hab1[k] = NULL;
	  ke++;
	}
      }
    }
    free(ks);
    tt1 = WallTime();
    dt = tt1-tt0;
    tt0 = tt1;
    double mem1=TotalSize();
    MPrintf(-1, "Time = %12.5E %g %g %d %d\n",
	    dt, mem1, mem1-mem0, ki, ke);
  }
  if (mbpt_tr.mktr > 0) {
    double mem0 = TotalSize();
    emax = (emax-emin);
    emin = 0.1;
    emax *= FINE_STRUCTURE_CONST;
    emin *= FINE_STRUCTURE_CONST;
    SetAWGridMBPT(emin, emax);
    InitTransitionMBPT(&mtr, n);
    double mem1 = TotalSize();
    MPrintf(-1, "TR Mem = %g %g %g\n", mem0, mem1, mem1-mem0);
  }

  if (n3 >= 0) {
    MPrintf(-1, "Construct Effective Hamiltonian %d\n", nc);
    fflush(stdout);
    for (k0 = 0; k0 < nc; k0++) {
      c0 = cs[k0];
      a = ZerothEnergyConfig(c0);
      b = ZerothResidualConfig(c0);
      for (m0 = 0; m0 < c0->n_csfs; m0++) {
	ms0 = c0->symstate[m0];
	UnpackSymState(ms0, &i0, &q0);
	if (meff[i0] && meff[i0]->nbasis > 0) {
	  meff[i0]->e0[q0] = a + 0.5*b;
	}
      }
    }
    double ttskip = 0, ttlock=0;
    long long tnlock = 0;
    ncps = 0;
    ResetWidMPI();
#pragma omp parallel default(shared) private(isym,n0,bra,ket,sbra,sket,bra1,ket1,bra2,ket2,sbra1,sket1,sbra2,sket2,cs,dt,dtt,k0,k1,c0,p0,c1,p1,m,bst0,kst0,m0,m1,ms0,ms1,q,q0,q1,k,mst,i0,i1,ct0,ct1,bst,kst,n1,bas0,bas1)
    {
      MBPT_EFF *imeff[MAX_SYMMETRIES];
      int cpmeff = 0;
#if CPMEFF == 1      
#if USE_MPI == 2
      if (NProcMPI() > 1) {
	cpmeff = 1;
      }
#endif
#endif
      if (cpmeff) {
	for (isym = 0; isym < MAX_SYMMETRIES; isym++) {
	  if (meff[isym] == NULL) {
	    imeff[isym] = NULL;
	    continue;
	  }
	  imeff[isym] = malloc(sizeof(MBPT_EFF));
	  memcpy(imeff[isym], meff[isym], sizeof(MBPT_EFF));
	  q = meff[isym]->hsize;
	  imeff[isym]->hab1 = malloc(sizeof(double *)*q);
	  imeff[isym]->hba1 = malloc(sizeof(double *)*q);
	  imeff[isym]->hab = malloc(sizeof(double *)*q);
	  imeff[isym]->hba = malloc(sizeof(double *)*q);
	  for (k = 0; k < q; k++) {
	    if (meff[isym]->hab1[k] == NULL) {
	      imeff[isym]->hab1[k] = NULL;
	    } else {
	      imeff[isym]->hab1[k] = malloc(sizeof(double)*nhab1);
	      imeff[isym]->hba1[k] = malloc(sizeof(double)*nhab1);
	      imeff[isym]->hab[k] = malloc(sizeof(double)*nhab);
	      imeff[isym]->hba[k] = malloc(sizeof(double)*nhab);
	      for (i0 = 0; i0 < nhab1; i0++) {
		imeff[isym]->hab1[k][i0] = 0.0;
		imeff[isym]->hba1[k][i0] = 0.0;
	      }	  
	      for (i0 = 0; i0 < nhab; i0++) {
		imeff[isym]->hab[k][i0] = 0.0;
		imeff[isym]->hba[k][i0] = 0.0;
	      }
	    }
	  }
	}
      } else {
	for (isym = 0; isym < MAX_SYMMETRIES; isym++) {
	  imeff[isym] = meff[isym];
	}
      }
      double ptt0, ptt1, tskip, tlock;
      long long nlock;
      cs = mbpt_cs;
      bas0 = mbpt_bas0;
      bas1 = mbpt_bas1;
      mbpt_ibas0.n = mbpt_ibas0.m = 0;
      mbpt_ibas1.n = mbpt_ibas1.m = 0;
      for (k0 = 0; k0 < nc; k0++) {
	c0 = cs[k0];      
	p0 = ConfigParity(c0);
	for (k1 = k0; k1 < nc; k1++) {
	  c1 = cs[k1];
	  p1 = ConfigParity(c1);
	  if (p0 != p1) continue;
	  /* pair of bra and ket states */
	  m = 0;      
	  bst0 = malloc(sizeof(int)*c0->n_csfs*c1->n_csfs);
	  kst0 = malloc(sizeof(int)*c0->n_csfs*c1->n_csfs);
	  for (m0 = 0; m0 < c0->n_csfs; m0++) {
	    ms0 = c0->symstate[m0];
	    UnpackSymState(ms0, &i0, &q0);	
	    if (c0 == c1) q = m0;
	    else q = 0;
	    for (m1 = q; m1 < c1->n_csfs; m1++) {
	      ms1 = c1->symstate[m1];
	      UnpackSymState(ms1, &i1, &q1);
	      if (i0 != i1) continue;
	      if (q0 <= q1) {
		k = q1*(q1+1)/2 + q0;
	      } else {
		k = q0*(q0+1)/2 + q1;
	      }
	      if (meff[i0] && meff[i0]->nbasis > 0 && meff[i0]->hab1[k]) {
		bst0[m] = m0;
		kst0[m] = m1;
		m++;
	      }
	    }
	  }
	  if (m == 0) {
	    continue;
	  }
	  ptt0 = tt0;
	  /* mst pairs */
	  mst = m;
	  /* if q0 <= q1 for the 1st pair, so are for the rest pairs */
	  ms0 = c0->symstate[bst0[0]];
	  ms1 = c1->symstate[kst0[0]];
	  UnpackSymState(ms0, &i0, &q0);
	  UnpackSymState(ms1, &i1, &q1);
	  if (q0 <= q1) {
	    ct0 = c0;
	    ct1 = c1;
	    bst = bst0;
	    kst = kst0;
	  } else {
	    ct0 = c1;
	    ct1 = c0;
	    bst = kst0;
	    kst = bst0;
	  }
	  /* make sure ct0 and ct1 have the same set of shells */
	  n0 = PadStates(ct0, ct1, &bra, &ket, &sbra, &sket);
	  /* pointers 1 starts from 2nd virtual orb. */
	  /* pointers 2 starts from the real orb. */
	  bra1 = bra + 1;
	  ket1 = ket + 1;
	  bra2 = bra + 2;
	  ket2 = ket + 2;
	  sbra1 = sbra + 1;
	  sket1 = sket + 1;
	  sbra2 = sbra + 2;
	  sket2 = sket + 2;
	  /* determine all real orbs */
	  for (k = 0; k < n0; k++) {	  
	    bas0[k] = OrbitalIndex(bra2[k].n, bra2[k].kappa, 0.0);
	    mbpt_bas0s[k] = 1000000;
	    mbpt_bas0d[k] = 1000000;
	    if(bra2[k].n <= mbpt_ne) {
	      int idx = IdxSD(bra2[k].n, bra2[k].kappa);
	      mbpt_bas0s[k] = mbpt_se[idx];
	      mbpt_bas0d[k] = mbpt_de[idx];
	      if (mbpt_bas0s[k] < 0) mbpt_bas0s[k] = 1000000;
	      if (mbpt_bas0d[k] < 0) mbpt_bas0d[k] = 1000000;
	    }
	  }	    
	  FreeIdxAry(&mbpt_ibas0, 2);
	  InitIdxAry(&mbpt_ibas0, n0, bas0);
	  /* determine all virtual orbs */
	  n1 = 0;
	  for (m = 0; m < nb; m++) {
	    k = IdxGet(&mbpt_ibas0, bas[m]);
	    if (k >= 0) continue;
	    bas1[n1] = bas[m];
	    n1++;
	  }
	  FreeIdxAry(&mbpt_ibas1, 2);
	  InitIdxAry(&mbpt_ibas1, n1, bas1);
	  mbptjp.nj = GetJpList(n1, bas1, mbptjp.jp);
	  if (n3 != 2) {
	    /* 1-b 2-b term no virtual orb */
	    DeltaH12M0(imeff, n0, bra2, ket2, sbra2, sket2, mst, bst, kst,
	    	       ct0, ct1, &mbpt_ibas0, mbpt_bas0s, mbpt_bas0d,
	    	       &ing, nc, cs, 0);	      
	    /* 1-b 2-b term 1 virtual orb */
	    DeltaH12M1(imeff, n0+1, bra1, ket1, sbra1, sket1, mst, bst, kst,
	    	       ct0, ct1, &mbpt_ibas0, mbpt_bas0s, mbpt_bas0d,
	    	       &mbpt_ibas1, &ing, nc, cs, 0);	  
	    /* 2-b 1-b term no virtual orb */
	    DeltaH12M0(imeff, n0, ket2, bra2, sket2, sbra2, mst, kst, bst,
	    	       ct1, ct0, &mbpt_ibas0, mbpt_bas0s, mbpt_bas0d,
	    	       &ing, nc, cs, 0);
	    /* 2-b 1-b term 1 virtual orb */
	    DeltaH12M1(imeff, n0+1, ket1, bra1, sket1, sbra1, mst, kst, bst,
	    	       ct1, ct0, &mbpt_ibas0, mbpt_bas0s, mbpt_bas0d,
	    	       &mbpt_ibas1, &ing, nc, cs, 0);	  
	    /* 1-b 1-b term no virtual orb */
	    DeltaH11M0(imeff, n0, bra2, ket2, sbra2, sket2, mst, bst, kst,
		       ct0, ct1, &mbpt_ibas0, mbpt_bas0s, &ing, nc, cs, 0);  
	    /* 1-b 1-b term 1 virtual orb */	    
	    DeltaH11M1(imeff, n0+1, bra1, ket1, sbra1, sket1, mst, bst, kst, 
		       ct0, ct1, &mbpt_ibas0, mbpt_bas0s,
		       &mbpt_ibas1, &ing, nc, cs, 0);
	    /* 2-b 2-b term no virtual */
	    DeltaH22M0(imeff, n0, bra2, ket2, sbra2, sket2, mst, bst, kst,
	    	       ct0, ct1, &mbpt_ibas0,
		       mbpt_bas0s, mbpt_bas0d, &ing, nc, cs);
	    /* 2-b 2-b term 1 virtual */
	    DeltaH22M1(imeff, n0+1, bra1, ket1, sbra1, sket1, mst, bst, kst,
	    	       ct0, ct1, &mbpt_ibas0, mbpt_bas0s, mbpt_bas0d,
	    	       &mbpt_ibas1, &ing, nc, cs);
	  }
	  
	  if (n3 != 1) {
	    /* 2-b 2-b term 2 virtual */
	    DeltaH22M2(imeff, n0+2, bra, ket, sbra, sket, mst, bst, kst,
		       ct0, ct1, &mbpt_ibas0, mbpt_bas0s, mbpt_bas0d,
		       &mbpt_ibas1, &ing, &ing2, nc, cs);
	  }
	  free(bra);
	  free(ket);
	  free(sbra);
	  free(sket);
	  free(bst0);
	  free(kst0);
	  ptt1 = WallTime();
	  dt = ptt1-ptt0;
	  dtt = ptt1-tbg;
	  tt0 = ptt1;
	  tskip = TimeSkip();
	  tlock = TimeLock();
	  nlock = NumLock();
	  double tmem = TotalSize();
	  MPrintf(0, "%3d %3d %3d %3d %3d %3d ... %12.5E %12.5E %12.5E %12.5E %12.5E %ld %ld\n", 
		  k0, k1, nc, mst, n0, n1, dt, dtt, tmem,
		  tskip, tlock, nlock, WidMPI());
	  fflush(stdout);	  
#pragma omp atomic
	    ncps++;
#pragma omp master
	  {
	    if ((mbpt_reinit_ncps > 0 && ncps >= mbpt_reinit_ncps) ||
		(mbpt_reinit_mem > 0 && tmem >= mbpt_reinit_mem)) {
	      SetRadialCleanFlags();
	      ncps = 0;
	    }
	  }
	}
      }
      if (cpmeff) {
#pragma omp critical
	{
	  ttskip += tskip;
	  ttlock += tlock;
	  tnlock += nlock;
	  MPrintf(-1, "Time Skip/Lock: %12.5E %12.5E %12.5E %12.5E %ld\n", tskip, tlock, ttskip, ttlock, tnlock);
	  for (isym = 0; isym < MAX_SYMMETRIES; isym++) {
	    if (meff[isym] == NULL) continue;
	    q = meff[isym]->hsize;
	    for (k = 0; k < q; k++) {
	      if (meff[isym]->hab1[k] == NULL) continue;
	      for (i0 = 0; i0 < nhab1; i0++) {
		meff[isym]->hab1[k][i0] += imeff[isym]->hab1[k][i0];
		meff[isym]->hba1[k][i0] += imeff[isym]->hba1[k][i0];
	      }
	      for (i0 = 0; i0 < nhab; i0++) {
		meff[isym]->hab[k][i0] += imeff[isym]->hab[k][i0];
		meff[isym]->hba[k][i0] += imeff[isym]->hba[k][i0];
	      }
	      free(imeff[isym]->hab[k]);
	      free(imeff[isym]->hba[k]);
	      free(imeff[isym]->hab1[k]);
	      free(imeff[isym]->hba1[k]);	      
	    }
	    free(imeff[isym]->hab);
	    free(imeff[isym]->hba);
	    free(imeff[isym]->hab1);
	    free(imeff[isym]->hba1);
	    free(imeff[isym]);
	  }
	}	    
      } else {
	ttskip = tskip;
	ttlock = tlock;
	tnlock = nlock;
      }
    }

    MPrintf(-1, "MBPT Structure.\n");
    fflush(stdout);
    for (isym = 0; isym < MAX_SYMMETRIES; isym++) {
      if (meff[isym] == NULL) continue;
      if (meff[isym]->nbasis == 0) {
	k = 0;
	if (MyRankMPI() == 0) {
	  fwrite(&isym, sizeof(int), 1, f);
	  fwrite(&k, sizeof(int), 1, f);
	}
	continue;
      }
      heff = meff[isym]->heff;      
#if USE_MPI == 1
      if (NProcMPI() > 1) {
	for (i = 0; i < meff[isym]->hsize; i++) {
	  if (meff[isym]->hab1[i] != NULL) {
	    MPI_Allreduce(MPI_IN_PLACE, meff[isym]->hab[i], nhab, MPI_DOUBLE,
			  MPI_SUM, MPI_COMM_WORLD);
	    MPI_Allreduce(MPI_IN_PLACE, meff[isym]->hba[i], nhab, MPI_DOUBLE,
			  MPI_SUM, MPI_COMM_WORLD);
	    MPI_Allreduce(MPI_IN_PLACE, meff[isym]->hab1[i], nhab1, MPI_DOUBLE,
			  MPI_SUM, MPI_COMM_WORLD);
	    MPI_Allreduce(MPI_IN_PLACE, meff[isym]->hba1[i], nhab1, MPI_DOUBLE,
			  MPI_SUM, MPI_COMM_WORLD);
	  }
	}
	ttskip = TimeSkip();
	ttlock = TimeLock();
	MPI_Allreduce(MPI_IN_PLACE, &ttskip, 1, MPI_DOUBLE,
		      MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &ttlock, 1, MPI_DOUBLE,
		      MPI_SUM, MPI_COMM_WORLD);
      }
#endif
      h = GetHamilton(isym);
      AllocHamMem(h, meff[isym]->nbasis, meff[isym]->nbasis);
      h0 = meff[isym]->h0;
      h->pj = isym;
      h->n_basis = meff[isym]->nbasis;
      h->dim = h->n_basis;
      memcpy(h->basis, meff[isym]->basis, sizeof(int)*h->n_basis);
      h->hsize = h->dim*(h->dim+1)/2;
      k = ConstructHamilton(isym, nkg0, nkg, kg, 0, NULL, 1);
      h->heff = heff;
      DecodePJ(isym, &pp, &jj);
      if (MyRankMPI() != 0) continue;   
      fwrite(&isym, sizeof(int), 1, f);
      fwrite(&(h->dim), sizeof(int), 1, f);
      fflush(f);
      fflush(stdout);
      for (j = 0; j < h->dim; j++) {
	for (i = 0; i <= j; i++) {
	  k = j*(j+1)/2 + i;
	  if (meff[isym]->hab1[k] == NULL) {
	    a = h0[k];
	    m = j*h->dim + i;
	    heff[m] = a;
	    if (i < j) {
	      m = i*h->dim + j;
	      heff[m] = a;
	    }
	    b = 0.0;
	    c = 0.0;
	    q = -i-1;
	    k = -j-1;
	    fwrite(&q, sizeof(int), 1, f);
	    fwrite(&k, sizeof(int), 1, f);
	    fwrite(&a, sizeof(double), 1, f);
	    fwrite(&b, sizeof(double), 1, f);
	    fwrite(&c, sizeof(double), 1, f);
	    continue;
	  }  
	  hab1 = meff[isym]->hab1[k];
	  hba1 = meff[isym]->hba1[k];
	  hab = meff[isym]->hab[k];
	  hba = meff[isym]->hba[k];	
	  a = sqrt(jj+1.0);
	  for (i0 = 0; i0 < nhab1; i0++) {
	    hab1[i0] /= a;
	    hba1[i0] /= a;
	  }
	  for (i0 = 0; i0 < nhab; i0++) {
	    hab[i0] /= a;
	    hba[i0] /= a;
	  }
	  a = h0[k];
	  m = j*h->dim + i;
	  b = SumInterpH(n, ng, n2, ng2, hab, hab1, dw);
	  heff[m] = a+b;
	  if (i < j) {
	    m = i*h->dim + j;
	    c = SumInterpH(n, ng, n2, ng2, hba, hba1, dw);
	    heff[m] = a+c;
	  } else {
	    c = b;
	  }
	  fwrite(&i, sizeof(int), 1, f);
	  fwrite(&j, sizeof(int), 1, f);
	  fwrite(&a, sizeof(double), 1, f);
	  fwrite(&b, sizeof(double), 1, f);
	  fwrite(&c, sizeof(double), 1, f);
	  if (n3 != 2) {
	    fwrite(hab1, sizeof(double), nhab1, f);
	    if (i != j) {
	      fwrite(hba1, sizeof(double), nhab1, f);
	    }
	  }
	  if (n3 != 1) {
	    fwrite(hab, sizeof(double), nhab, f);
	    if (i != j) {
	      fwrite(hba, sizeof(double), nhab, f);
	    }
	  }
	}
      }
      fflush(f);
      if (DiagnolizeHamilton(h) < 0) {
	MPrintf(-1, "Diagnolizing Effective Hamiltonian Error\n");
	ierr = -1;
	goto ERROR;
      }
      AddToLevels(h, nkg0, kg);
      h->heff = NULL;
      tt1 = WallTime();
      dt = tt1-tt0;
      tt0 = tt1;
      MPrintf(-1, "Time = %12.5E, isym=%d\n", dt, isym);
      fflush(stdout);
    }
    if (MyRankMPI() == 0) {
      SortLevels(nlevels, -1, 0);
      SaveLevels(fn, nlevels, -1);
      fclose(f);
      
      tt1 = WallTime();
      dt = tt1 - tbg;
      tt0 = tt1;
      MPrintf(-1, "Total Time Structure= %12.5E %12.5E %12.5E %ld\n", dt, ttskip, ttlock, tnlock);
      fflush(stdout);
    }
  }
#if USE_MPI == 1
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  if (mbpt_tr.mktr > 0) {
    MPrintf(-1, "MBPT Transition.\n");
    fflush(stdout);
    ReinitRadial(2);
    if (MyRankMPI() == 0) {
      sprintf(tfn, "%s.tr", fn1);
      f = fopen(tfn, "w");
    }
    ncps = 0;
    ResetWidMPI();
#pragma omp parallel default(shared) private(n0,bra,ket,sbra,sket,bra1,ket1,bra2,ket2,sbra1,sket1,sbra2,sket2,cs,dt,dtt,k0,k1,c0,p0,c1,p1,m,bst0,kst0,m0,m1,ms0,ms1,q,q0,q1,k,mst,i0,i1,ct0,ct1,bst,kst,n1,bas0,bas1)
    {
      MBPT_TR *imtr;
      int cpmtr = 0;
#if CPMTR == 1
#if USE_MPI == 2
      if (mbpt_cpmtr && NProcMPI() > 1) {
	cpmtr = 1;
      }
#endif
#endif
      if (!cpmtr) imtr = mtr;
      else {
	InitTransitionMBPT(&imtr, n);
      }
      double ptt0, ptt1;
      cs = mbpt_cs;
      bas0 = mbpt_bas0;
      bas1 = mbpt_bas1;
      for (k0 = 0; k0 < nc; k0++) {
	c0 = cs[k0];
	for (k1 = 0; k1 < nc; k1++) {
	  c1 = cs[k1];
	  m = 0;      
	  bst0 = malloc(sizeof(int)*c0->n_csfs*c1->n_csfs);
	  kst0 = malloc(sizeof(int)*c0->n_csfs*c1->n_csfs);
	  for (m0 = 0; m0 < c0->n_csfs; m0++) {
	    ms0 = c0->symstate[m0];
	    UnpackSymState(ms0, &i0, &q0);	
	    DecodePJ(i0, &p0, &j0);
	    for (m1 = 0; m1 < c1->n_csfs; m1++) {
	      ms1 = c1->symstate[m1];
	      UnpackSymState(ms1, &i1, &q1);
	      DecodePJ(i1, &p1, &j1);
	      if (abs(j0 - j1) > 2*mbpt_tr.mktr) continue;
	      bst0[m] = m0;
	      kst0[m] = m1;
	      m++;
	    }
	  }
	  
	  if (m == 0) {
	    continue;
	  }
	  ptt0 = tt0;
	  /* mst pairs */
	  mst = m;
	  ct0 = c0;
	  ct1 = c1;
	  bst = bst0;
	  kst = kst0;
	  /* make sure ct0 and ct1 have the same set of shells */
	  n0 = PadStates(ct0, ct1, &bra, &ket, &sbra, &sket);
	  /* pointers 1 starts from 2nd virtual orb. */
	  /* pointers 2 starts from the real orb. */
	  bra1 = bra + 1;
	  ket1 = ket + 1;
	  bra2 = bra + 2;
	  ket2 = ket + 2;
	  sbra1 = sbra + 1;
	  sket1 = sket + 1;
	  sbra2 = sbra + 2;
	  sket2 = sket + 2;	
	  /* determine all real orbs */
	  int ks = 0;
	  int kd = 0;
	  for (k = 0; k < n0; k++) {	  
	    bas0[k] = OrbitalIndex(bra2[k].n, bra2[k].kappa, 0.0);
	    mbpt_bas0s[k] = 1000000;
	    mbpt_bas0d[k] = 1000000;
	    if(bra2[k].n <= mbpt_ne) {
	      int idx = IdxSD(bra2[k].n, bra2[k].kappa);
	      mbpt_bas0s[k] = mbpt_se[idx];
	      mbpt_bas0d[k] = mbpt_de[idx];
	      if (mbpt_bas0s[k] < 0) mbpt_bas0s[k] = 1000000;
	      if (mbpt_bas0d[k] < 0) mbpt_bas0d[k] = 1000000;
	    }
	  }
	  FreeIdxAry(&mbpt_ibas0, 2);
	  InitIdxAry(&mbpt_ibas0, n0, bas0);
	  /* determine all virtual orbs */
	  n1 = 0;
	  for (m = 0; m < nb; m++) {
	    k = IdxGet(&mbpt_ibas0, bas[m]);
	    if (k >= 0) continue;
	    bas1[n1] = bas[m];
	    n1++;
	  }
	  FreeIdxAry(&mbpt_ibas1, 2);
	  InitIdxAry(&mbpt_ibas1, n1, bas1);
	  mbptjp.nj = GetJpList(n1, bas1, mbptjp.jp);
	  /* 1-b 2-b term no virtual orb */
	  DeltaH12M0(imtr, n0, bra2, ket2, sbra2, sket2, mst, bst, kst,
		     ct0, ct1, &mbpt_ibas0, mbpt_bas0s, mbpt_bas0d,
		     &ing, nc, cs, 1);	    
	  /* 1-b 2-b term 1 virtual orb */
	  DeltaH12M1(imtr, n0+1, bra1, ket1, sbra1, sket1, mst, bst, kst,
		     ct0, ct1, &mbpt_ibas0, mbpt_bas0s, mbpt_bas0d,
		     &mbpt_ibas1, &ing, nc, cs, 1);
	  /* 1-b 1-b term no virtual orb */
	  DeltaH11M0(imtr, n0, bra2, ket2, sbra2, sket2, mst, bst, kst,
		     ct0, ct1, &mbpt_ibas0, mbpt_bas0s, &ing, nc, cs, 1);
	  /* 1-b 1-b term 1 virtual orb */
	  DeltaH11M1(imtr, n0+1, bra1, ket1, sbra1, sket1, mst, bst, kst, 
		     ct0, ct1, &mbpt_ibas0, mbpt_bas0s,
		     &mbpt_ibas1, &ing, nc, cs, 1);
	  free(bra);
	  free(ket);
	  free(sbra);
	  free(sket);
	  free(bst0);
	  free(kst0);
	  ptt1 = WallTime();
	  dt = ptt1-ptt0;
	  dtt = ptt1-tbg;
	  tt0 = ptt1;
	  double tmem = TotalSize();
	  MPrintf(0, "%3d %3d %3d %3d %3d %3d ... %12.5E %12.5E %12.5E\n", 
		  k0, k1, nc, mst, n0, n1, dt, dtt, tmem);
	  fflush(stdout);	  
#pragma omp atomic
	  ncps++;
#pragma omp master
	  {
	    if ((mbpt_reinit_ncps > 0 && ncps >= mbpt_reinit_ncps) ||
		(mbpt_reinit_mem > 0 && tmem >= mbpt_reinit_mem)) {
	      SetRadialCleanFlags();
	      ncps = 0;
	    }
	  }
	}
      }
      if (cpmtr) {
#pragma omp critical
	{
	  k = 2*mbpt_tr.mktr*MAX_SYMMETRIES;
	  for (j = 0; j < k; j++) {
	    if (mtr[j].nsym1 == 0) continue;
	    for (m = 0; m < mtr[j].nsym1; m++) {
	      mst = mtr[j].sym0->n_states*mtr[j].sym1[m]->n_states;
	      mst *= n*mbpt_tr.naw;
	      if (mst > 0) {
		for (q = 0; q < mst; q++) {
		  mtr[j].tma[m][q] += imtr[j].tma[m][q];
		  mtr[j].rma[m][q] += imtr[j].rma[m][q];
		}
	      }
	    }
	  }
	  FreeTransitionMBPT(imtr);
	}
      }
    }

    if (MyRankMPI() == 0) {
      fwrite(&mbpt_tr.mktr, sizeof(int), 1, f);
      fwrite(&mbpt_tr.naw, sizeof(int), 1, f);
      fwrite(&emin, sizeof(double), 1, f);
      fwrite(&emax, sizeof(double), 1, f);
    }
    k = 2*mbpt_tr.mktr*MAX_SYMMETRIES;
    for (j = 0; j < k; j++) {
      if (mtr[j].nsym1 == 0) continue;
      for (m = 0; m < mtr[j].nsym1; m++) {
	mst = mtr[j].sym0->n_states * mtr[j].sym1[m]->n_states;
	mst *= n * mbpt_tr.naw;	
	if (mst > 0) {
#if USE_MPI == 1
	  if (NProcMPI() > 1) {
	    MPI_Allreduce(MPI_IN_PLACE, mtr[j].tma[m], mst,
			  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	    MPI_Allreduce(MPI_IN_PLACE, mtr[j].rma[m], mst,
			  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	  }
#endif
	  if (MyRankMPI() == 0) {
	    fwrite(mtr[j].tma[m], sizeof(double), mst, f);
	    fwrite(mtr[j].rma[m], sizeof(double), mst, f);
	  }
	}
      }
    }		
    
    tt1 = WallTime();
    dt = tt1 - tbg;
    tt0 = tt1;
    MPrintf(-1, "Total Time Transition = %12.5E\n", dt);
    fflush(stdout);
    if (MyRankMPI() == 0) fclose(f);
  }
  
 ERROR:
  FreeIdxAry(&mbpt_ibas0, 2);
  FreeIdxAry(&mbpt_ibas1, 2);
  FreeIdxAry(&ing, 2);
  FreeIdxAry(&ing2, 2);
  ResetWidMPI();
#pragma omp parallel default(shared) private(cs)
  {
    cs = mbpt_cs;
    free(cs[nc]->shells);
    free(cs);
    free(mbpt_bas0);
    free(mbpt_bas0s);
    free(mbpt_bas0d);
    free(mbpt_bas1);
    free(mbptjp.jp);
    if (mbpt_nsplit) {
      FreeIdxAry(&mbptjp.ibs, 2);
    }
  }
  free(bas);
  free(dw);
  if (mbpt_nsplit) {
    for (i = 0; i < nb; i++) {
      free(mbpt_rij[i]);
    }
    free(mbpt_rij);
  }
  FreeEffMBPT(meff);
  FreeTransitionMBPT(mtr);
  return ierr;
}

int ReadMBPT(int nf, FILE *f[], MBPT_HAM *mbpt, int m) {
  int i, j, k, q, nb, nn2, n;
  double emin, emax;
  MBPT_TR *mtr;
    
  if (m == 0) {
    for (i = 0; i < nf; i++) {
      nb = fread(&(mbpt[i].n), sizeof(int), 1, f[i]);
      if (nb != 1) return -1;
      mbpt[i].ng = malloc(sizeof(int)*mbpt[i].n);
      nb = fread(mbpt[i].ng, sizeof(int), mbpt[i].n, f[i]);
      if (nb != mbpt[i].n) return -1;      
      nb = fread(&(mbpt[i].n2), sizeof(int), 1, f[i]);
      if (nb != 1) return -1;
      mbpt[i].ng2 = malloc(sizeof(int)*mbpt[i].n2);
      nb = fread(mbpt[i].ng2, sizeof(int), mbpt[i].n2, f[i]);
      if (nb != mbpt[i].n2) return -1;
      nb = fread(&(mbpt[i].n3), sizeof(int), 1, f[i]);
      if (nb != 1) return -1;
      mbpt[i].hab1 = malloc(sizeof(double)*mbpt[i].n*2);
      mbpt[i].hba1 = malloc(sizeof(double)*mbpt[i].n*2);
      mbpt[i].hab = malloc(sizeof(double)*mbpt[i].n*mbpt[i].n2*2);
      mbpt[i].hba = malloc(sizeof(double)*mbpt[i].n*mbpt[i].n2*2);
    }
    return 0;
  }

  if (m == 1) {
    for (i = 0; i < nf; i++) {
      nb = fread(&(mbpt[i].isym), sizeof(int), 1, f[i]);
      if (nb != 1) return -1;
      nb = fread(&(mbpt[i].dim), sizeof(int), 1, f[i]);
      if (nb != 1) return -1;
    }
    return 0;
  }
  
  if (m == 2) {
    for (i = 0; i < nf; i++) {
      if (mbpt[i].dim == 0) continue;
      nb = fread(&(mbpt[i].ibra), sizeof(int), 1, f[i]);
      if (nb != 1) return -1;
      nb = fread(&(mbpt[i].iket), sizeof(int), 1, f[i]);
      if (nb != 1) return -1;
      nb = fread(&(mbpt[i].a), sizeof(double), 1, f[i]);
      if (nb != 1) return -1;
      nb = fread(&(mbpt[i].b), sizeof(double), 1, f[i]);
      if (nb != 1) return -1;
      nb = fread(&(mbpt[i].c), sizeof(double), 1, f[i]);
      if (nb != 1) return -1;
      if (mbpt[i].ibra < 0 || mbpt[i].iket < 0) continue;
      if (mbpt[i].n3 != 2) {
	nb = fread(mbpt[i].hab1, sizeof(double), mbpt[i].n*2, f[i]);
	if (nb != mbpt[i].n*2) return -1;
	if (mbpt[i].ibra != mbpt[i].iket) {
	  nb = fread(mbpt[i].hba1, sizeof(double), mbpt[i].n*2, f[i]);
	  if (nb != mbpt[i].n*2) return -1;
	} else {
	  memcpy(mbpt[i].hba1, mbpt[i].hab1, sizeof(double)*mbpt[i].n*2);
	}
      } 
      if (mbpt[i].n3 != 1) {
	nn2 = 2*mbpt[i].n*mbpt[i].n2;
	nb = fread(mbpt[i].hab, sizeof(double), nn2, f[i]);
	if (nb != nn2) return -1;
	if (mbpt[i].ibra != mbpt[i].iket) {
	  nb = fread(mbpt[i].hba, sizeof(double), nn2, f[i]);
	  if (nb != nn2) return -1;
	} else {
	  memcpy(mbpt[i].hba, mbpt[i].hab, sizeof(double)*nn2);
	}
      }
    }
    return 0;
  }
  
  if (m == 3) {
    for (i = 0; i < nf; i++) {
      nb = fread(&mbpt_tr.mktr, sizeof(int), 1, f[i]);
      if (nb != 1) return -1;
      nb = fread(&mbpt_tr.naw, sizeof(int), 1, f[i]);
      if (nb != 1) return -1;
      nb = fread(&emin, sizeof(double), 1, f[i]);
      if (nb != 1) return -1;
      nb = fread(&emax, sizeof(double), 1, f[i]);
      if (nb != 1) return -1;
      InitTransitionMBPT(&(mbpt[i].mtr), mbpt[i].n);
      k = 2*mbpt_tr.mktr*MAX_SYMMETRIES;
      mtr = mbpt[i].mtr;
      for (j = 0; j < k; j++) {
	if (mtr[j].nsym1 == 0) continue;
	for (q = 0; q < mtr[j].nsym1; q++) {
	  n = mtr[j].sym0->n_states * mtr[j].sym1[q]->n_states;
	  n *= mbpt[i].n * mbpt_tr.naw;
	  nb = fread(mtr[j].tma[q], sizeof(double), n, f[i]);
	  if (nb != n) return -1;
	  nb = fread(mtr[j].rma[q], sizeof(double), n, f[i]);
	  if (nb != n) return -1;
	}
      }
    }    
    SetAWGridMBPT(emin, emax);
    return 0;
  }

  return 0;
}

void CombineMBPT(int nf, MBPT_HAM *mbpt, 
		 double *hab1, double *hba1,
		 double **hab, double **hba, 
		 double *nab1, double *nba1,
		 double **nab, double **nba,
		 int n0, int *ng0, int n, int *ng, int *n2, int **ng2) {
  int m, i, j, k, q, r, nn2;

  for (m = 0; m < nf; m++) {
    if (mbpt[m].ibra < 0 || mbpt[m].iket < 0) continue;
    if (mbpt[m].n3 != 2) {
      for (i = 0; i < mbpt[m].n; i++) {
	k = IBisect(mbpt[m].ng[i], n0, ng0);
	if (k < 0) continue;
	hab1[k] = mbpt[m].hab1[i];
	hba1[k] = mbpt[m].hba1[i];
	nab1[k] = mbpt[m].hab1[i+mbpt[m].n];
	nba1[k] = mbpt[m].hba1[i+mbpt[m].n];
      }
    }
    if (mbpt[m].n3 != 1) {
      nn2 = mbpt[m].n * mbpt[m].n2;
      for (i = 0; i < mbpt[m].n; i++) {	
	k = IBisect(mbpt[m].ng[i], n, ng);
	if (k < 0) continue;
	r = i*mbpt[m].n2;
	for (j = 0; j < mbpt[m].n2; j++) {
	  q = IBisect(mbpt[m].ng2[j], n2[k], ng2[k]);
	  if (q < 0) continue;
	  hab[k][q] = mbpt[m].hab[r+j];
	  hba[k][q] = mbpt[m].hba[r+j];
	  nab[k][q] = mbpt[m].hab[r+j + nn2];
	  nba[k][q] = mbpt[m].hba[r+j + nn2];
	}
      }
    }
  } 
}

void CombineTransitionMBPT(int nf, MBPT_HAM *mbpt, MBPT_TR *mtr, int n, int *ng) {
  int m, k, j, i, t, q, r, p, is0, is1;

  for (m = 0; m < nf; m++) {
    k = 2*mbpt_tr.mktr*MAX_SYMMETRIES;
    for (j = 0; j < k; j++) {
      if (mtr[j].nsym1 == 0) continue;
      for (q = 0; q < mtr[j].nsym1; q++) {
	for (is0 = 0; is0 < mtr[j].sym0->n_states; is0++) {
	  for (is1 = 0; is1 < mtr[j].sym1[q]->n_states; is1++) {
	    for (p = 0; p < mbpt_tr.naw; p++) {
	      t = (is0*mtr[j].sym1[q]->n_states+is1)*mbpt_tr.naw + p;
	      for (i = 0; i < mbpt[m].n; i++) {
		r = IBisect(mbpt[m].ng[i], n, ng);
		if (r < 0) continue;
		mtr[j].tma[q][t*n+r] = mbpt[m].mtr[j].tma[q][t*mbpt[m].n+i];
		mtr[j].rma[q][t*n+r] = mbpt[m].mtr[j].rma[q][t*mbpt[m].n+i];
	      }
	    }
	  }
	}
      }
    }
  }
}

void AdjustAngularZ(MBPT_TR *mtr) {
  SHAMILTON *h;
  ANGZ_DATUM *ad;
  ANGULAR_ZMIX **a;
  int nh, i0, i1, m0, m1, k0, k1;
  int i, k, iz, n, ns, p, q, p0, p1, j0, j1;

  SetMaxKMBPT(mbpt_tr.mktr);
  h = GetSHamilton(&nh);
  for (i0 = 0; i0 < nh; i0++) {
    DecodePJ(h[i0].pj, &p0, &j0);
    for (i1 = i0; i1 < nh; i1++) {
      DecodePJ(h[i1].pj, &p1, &j1);
      ns = AngularZMixStates(&ad, i0, i1);
      a = (ANGULAR_ZMIX **) ad->angz;
      ad->mk = malloc(sizeof(double *)*ns*2);
      n = mbpt_tr.mktr * mbpt_tr.naw;
      for (i = 0; i < 2*ns; i++) {
	ad->mk[i] = malloc(sizeof(double)*n);
	for (k = 0; k < n; k++) {
	  ad->mk[i][k] = 0.0;
	}
      }
      for (m0 = 0; m0 < h[i0].nbasis; m0++) {
	for (m1 = 0; m1 < h[i1].nbasis; m1++) {
	  iz = m0*h[i1].nbasis + m1;
	  for (k = 0; k < mbpt_tr.mktr; k++) {
	    n = mbpt_tr.naw * k;
	    k0 = h[i0].pj*2*mbpt_tr.mktr + 2*k;
	    if (IsEven(p0+p1+k)) k0++;
	    k1 = IBisect(h[i1].pj, mtr[k0].nsym1, mtr[k0].isym1);
	    if (k1 >= 0) {
	      p = m0*mtr[k0].sym1[k1]->n_states + m1;
	      p *= mbpt_tr.naw;
	      for (q = 0; q < mbpt_tr.naw; q++) {
		ad->mk[iz][n+q] += mtr[k0].tma[k1][p+q];
		if (i0 != i1) {
		  ad->mk[iz+ns][n+q] += mtr[k0].rma[k1][p+q];
		}
	      }
	    }
	    k0 = h[i1].pj*2*mbpt_tr.mktr + 2*k;
	    if (IsEven(p0+p1+k)) k0++;
	    k1 = IBisect(h[i0].pj, mtr[k0].nsym1, mtr[k0].isym1);
	    if (k1 >= 0) {
	      p = m1*mtr[k0].sym1[k1]->n_states + m0;
	      p *= mbpt_tr.naw;
	      for (q = 0; q < mbpt_tr.naw; q++) {
		ad->mk[iz][n+q] += mtr[k0].rma[k1][p+q];
		if (i0 != i1) {
		  ad->mk[iz+ns][n+q] += mtr[k0].tma[k1][p+q];
		}
	      }
	    }
	  }
	}
      }
    }
  }
}

void SaveTransitionMBPT(MBPT_TR *mtr) {
  char *fn;
  TFILE *f;
  LEVEL *lev1, *lev2;
  SYMMETRY *sym;
  STATE *st;
  TR_RECORD r;
  TR_HEADER tr_hdr;
  F_HEADER fhdr;
  double *awgrid, *rg, a, x, e, s, s0;
  int n, i, j, k, m, t, q, m1, m2, p;
  int i0, i1, p1, p2, j1, j2;

  if (MyRankMPI() != 0) return;
  fn = mbpt_tr.tfn;
  if (fn == NULL || mbpt_tr.nlow <= 0 || mbpt_tr.nup <= 0) return;  
  awgrid = mbpt_tr.awgrid;
  fhdr.type = DB_TR;
  strcpy(fhdr.symbol, GetAtomicSymbol());
  fhdr.atom = GetAtomicNumber();
  tr_hdr.nele = GetNumElectrons(0);
  n = GetNumLevels();
  f = OpenFile(fn, &fhdr);
  for (t = 1; t <= mbpt_tr.mktr; t++) {
    for (q = -1; q <= 1; q += 2) {      
      m = t*q;
      printf("Multipole %2d\n", m);
      tr_hdr.multipole = m;
      tr_hdr.gauge = GetTransitionGauge();
      tr_hdr.mode = 0;
      InitFile(f, &fhdr, &tr_hdr);
      ResetWidMPI();
#pragma omp parallel default(shared) private(i, j, k, sym, st, lev1, lev2, e, p1, j1, p2, j2, s0, i0, i1, m1, m2, a, p, rg, x, s, r)
      {
      for (j = 1; j < n; j++) {
	lev2 = GetLevel(j);
	k = lev2->pb;
	sym = GetSymmetry(lev2->pj);
	st = (STATE *) ArrayGet(&(sym->states), k);
	k = InGroups(st->kgroup, mbpt_tr.nup, mbpt_tr.up);
	if (k == 0) continue;
	DecodePJ(lev2->pj, &p2, &j2);
	for (i = 0; i < j; i++) {
	  lev1 = GetLevel(i);
	  k = lev1->pb;
	  sym = GetSymmetry(lev1->pj);
	  st = (STATE *) ArrayGet(&(sym->states), k);
	  k = InGroups(st->kgroup, mbpt_tr.nlow, mbpt_tr.low);
	  if (k == 0) continue;
	  int skip = SkipMPI();
	  if (skip) continue;
	  DecodePJ(lev1->pj, &p1, &j1);
	  e = 0.0;
	  k = TRMultipole(&s0, &e, m, i, j);
	  e *= FINE_STRUCTURE_CONST;
	  if (k != 0) continue;
	  s = 0;
	  i0 = lev1->pj*2*mbpt_tr.mktr + 2*(t-1);
	  if (q > 0) i0++;
	  i1 = IBisect(lev2->pj, mtr[i0].nsym1, mtr[i0].isym1);
	  if (i1 >= 0) {
	    for (m1 = 0; m1 < lev1->n_basis; m1++) {
	      for (m2 = 0; m2 < lev2->n_basis; m2++) {
		a = lev1->mixing[m1] * lev2->mixing[m2];
		p = lev1->basis[m1] * mtr[i0].sym1[i1]->n_states + lev2->basis[m2];
		p *= mbpt_tr.naw;	      
		rg = &(mtr[i0].tma[i1][p]);
		x = InterpolateMultipole(e, mbpt_tr.naw, awgrid, rg);
		if (tr_hdr.gauge == G_COULOMB && m < 0) {
		  x /= e;
		}
		s += a*x;
	      }
	    }
	  }
	  i0 = lev2->pj*2*mbpt_tr.mktr + 2*(t-1);
	  if (q > 0) i0++;
	  i1 = IBisect(lev1->pj, mtr[i0].nsym1, mtr[i0].isym1);
	  if (i1 >= 0) {
	    for (m1 = 0; m1 < lev1->n_basis; m1++) {
	      for (m2 = 0; m2 < lev2->n_basis; m2++) {
		a = lev1->mixing[m1] * lev2->mixing[m2];
		p = lev2->basis[m2] * mtr[i0].sym1[i1]->n_states + lev1->basis[m1];
		p *= mbpt_tr.naw;
		rg = &(mtr[i0].rma[i1][p]);
		x = InterpolateMultipole(e, mbpt_tr.naw, awgrid, rg);
		if (tr_hdr.gauge == G_COULOMB && m < 0) {
		  x /= e;
		}
		s += a*x;
	      }
	    }
	  }
	  r.lower = i;
	  r.upper = j;
	  r.strength = s0+s;
	  WriteTRRecord(f, &r, NULL);
	}
      }
      }
      DeinitFile(f, &fhdr);
    }
  }
  CloseFile(f, &fhdr);
  if (mbpt_tr.nlow > 0) {
    free(mbpt_tr.low);
  }
  if (mbpt_tr.nup > 0) {
    free(mbpt_tr.up);
  }
  mbpt_tr.nlow = 0;
  mbpt_tr.nup = 0;
  mbpt_tr.low = NULL;
  mbpt_tr.up = NULL;
}

int StructureReadMBPT(char *fn, char *fn2, int nf, char *fn1[],
		      int nkg, int *kg, int nkg0) {
  int ierr, m, i, j, k, k0, k1, nlevels, n2m;
  int isym, pp, jj, r, q, n, *ng, n0, *ng0, *n2, **ng2;
  double *dw, *z1, *z2, *x, *y, *z, *t, *heff, **hab, **hba;
  double **nab, **nba, *nab1, *nba1, *neff;
  double a, b, na, nb;
  HAMILTON *h;
  SYMMETRY *sym;
  MBPT_HAM *mbpt;
  MBPT_TR *mtr;
  char tfn1[1024];
  FILE **f1, *f2;

  if (MyRankMPI() != 0) return 0;
  ierr = 0;
  if (nkg0 <= 0 || nkg0 > nkg) nkg0 = nkg;
  mbpt = malloc(sizeof(MBPT_HAM)*nf);
  f1 = malloc(sizeof(FILE *)*nf);
  for (m = 0; m < nf; m++) {
    f1[m] = fopen(fn1[m], "r");
    if (f1[m] == NULL) {
      printf("cannot open file %s\n", fn1[m]);
      return -1;
    }
  }
  f2 = fopen(fn2, "w");
  if (f2 == NULL) {
    printf("cannot open file %s\n", fn2);
    return -1;
  }
  
  ierr = ReadMBPT(nf, f1, mbpt, 0);
  if (ierr < 0) {
    printf("cannot read MBPT Hamilton\n");
    return -1;
  }

  n = 0;
  n0 = 0;
  for (m = 0; m < nf; m++) {
    if (mbpt[m].n3 != 2) n0 += mbpt[m].n;
    if (mbpt[m].n3 != 1) n += mbpt[m].n;
  }    
  if (n0 > 0) ng0 = malloc(sizeof(int)*n0);
  if (n > 0) ng = malloc(sizeof(int)*n);
  i = 0;
  j = 0;
  for (m = 0; m < nf; m++) {
    if (mbpt[m].n3 != 1) {
      for (k = 0; k < mbpt[m].n; k++) {      
	ng[i++] = mbpt[m].ng[k];
      }
    }
    if (mbpt[m].n3 != 2) {
      for (k = 0; k < mbpt[m].n; k++) {
	ng0[j++] = mbpt[m].ng[k];
      }
    }
  }
  if (i > 0) {
    n = SortUnique(i, ng);
  }
  if (j > 0) {
    n0 = SortUnique(j, ng0);  
  }
  n2m = 0;
  if (n > 0) {
    n2 = malloc(sizeof(int)*n);
    ng2 = malloc(sizeof(int *)*n);    
    for (i = 0; i < n; i++) {
      n2[i] = 0;
      for (m = 0; m < nf; m++) {
	if (mbpt[m].n3 != 1) {
	  n2[i] += mbpt[m].n2;
	}
      }
      if (n2[i] > 0) {
	ng2[i] = malloc(sizeof(int)*n2[i]);
	j = 0;
	for (m = 0; m < nf; m++) {
	  if (mbpt[m].n3 != 1) {
	    for (k = 0; k < mbpt[m].n2; k++) {
	      ng2[i][j++] = mbpt[m].ng2[k];
	    }
	  }
	}
	if (j > 0) {
	  n2[i] = SortUnique(j, ng2[i]);
	}
      }
      if (n2[i] > n2m) n2m = n2[i];
    }
    hab = malloc(sizeof(double *)*n);
    hba = malloc(sizeof(double *)*n);
    nab = malloc(sizeof(double *)*n);
    nba = malloc(sizeof(double *)*n);
    for (i = 0; i < n; i++) {
      if (n2[i] > 0) {
	hab[i] = malloc(sizeof(double)*n2[i]);
	hba[i] = malloc(sizeof(double)*n2[i]);
	nab[i] = malloc(sizeof(double)*n2[i]);
	nba[i] = malloc(sizeof(double)*n2[i]);
      }
    }
  }
  n2m += n0;
  dw = malloc(sizeof(double)*n2m*7);
  z1 = dw;
  z2 = z1 + n2m;
  x = z2 + n2m;
  t = x + n2m;
  y = t + n2m;
  nab1 = y + n2m;
  nba1 = nab1 + n2m;

  printf("MBPT Structure.\n");
  fflush(stdout);
  nlevels = GetNumLevels();
  for (isym = 0; isym < MAX_SYMMETRIES; isym++) {
    k0 = ConstructHamilton(isym, nkg0, nkg, kg, 0, NULL, 101);
    if (k0 == -1) continue;
    h = GetHamilton(isym);
    sym = GetSymmetry(isym);
    DecodePJ(isym, &pp, &jj);
    ierr = ReadMBPT(nf, f1, mbpt, 1);
    if (ierr < 0) {
      printf("Error reading MBPT 1: %d\n", isym);
      goto ERROR;
    }
    k = 1;
    for (m = 0; m < nf; m++) {
      if (mbpt[m].dim == 0) {
	k = 0;
	break;
      }
    }
    if (k0 == 0) printf("sym: %3d %d %2d %d %3d\n", isym, pp, jj, k, h->dim);
    if (k == 0 || k0 < 0) {
      for (j = 0; j < h->dim; j++) {
	for (i = 0; i <= j; i++) {	  
	  ierr = ReadMBPT(nf, f1, mbpt, 2);
	}
      }
      continue;
    }
    heff = malloc(sizeof(double)*h->dim*h->dim);
    neff = malloc(sizeof(double)*h->dim*h->dim);
    h->heff = heff;
    for (j = 0; j < h->dim; j++) {
      for (i = 0; i <= j; i++) {
	k = j*(j+1)/2 + i;
	k0 = j*h->dim + i;
	k1 = i*h->dim + j;
	ierr = ReadMBPT(nf, f1, mbpt, 2);
	if (ierr < 0) {
	  printf("Error reading MBPT 2: %d %d %d\n", isym, i, j);
	  goto ERROR;	
	}
	heff[k0] = mbpt[0].a;
	heff[k1] = mbpt[0].a;
	neff[k0] = 0.0;
	neff[k1] = 0.0;
	if (mbpt[0].ibra < 0 || mbpt[0].iket < 0) goto OUT;
	CombineMBPT(nf, mbpt, z1, z2, hab, hba, nab1, nba1, nab, nba, 
		    n0, ng0, n, ng, n2, ng2);
	for (q = 0; q < n0; q++) {
	  fprintf(f2, "  %3d %3d %3d %3d %3d %12.5E %12.5E %12.5E %12.5E\n", 
		  isym, i, j, 0, ng0[q], z1[q], z2[q], nab1[q], nba1[q]);
	}
	for (r = 0; r < n; r++) {
	  for (q = 0; q < n2[r]; q++) {
	    fprintf(f2, "  %3d %3d %3d %3d %3d %12.5E %12.5E %12.5E %12.5E\n", 
		    isym, i, j, ng[r], ng[r]+ng2[r][q], 
		    hab[r][q], hba[r][q], nab[r][q], nba[r][q]);
	  }
	}
	fflush(f2);
	if (n0 > 0) {
	  for (q = 0; q < n0; q++) {
	    x[q] = ng0[q];
	  }
	  a = SumInterp1D(n0, z1, x, t, y);
	  heff[k0] += a;
	  na = SumInterp1D(n0, nab1, x, t, y);
	  neff[k0] += na;
	  b = a;
	  nb = na;
	  if (i < j) {
	    b = SumInterp1D(n0, z2, x, t, y);
	    heff[k1] += b;
	    nb = SumInterp1D(n0, nba1, x, t, y);
	    neff[k1] += nb;
	  }
	  fprintf(f2, "# %3d %3d %3d %3d %3d %12.5E %12.5E %15.8E %12.5E %12.5E\n",
		  isym, pp, jj, i, j, a, b, mbpt[0].a, na, nb);
	}
	if (n > 0) {
	  for (r = 0; r < n; r++) {
	    if (n2[r] > 0) {
	      for (q = 0; q < n2[r]; q++) {
		x[q] = ng[r] + ng2[r][q];
	      }
	      z1[r] = SumInterp1D(n2[r], hab[r], x, t, y);
	      if (i < j) z2[r] = SumInterp1D(n2[r], hba[r], x, t, y);
	      else z2[r] = z1[r];		
	      nab1[r] = SumInterp1D(n2[r], nab[r], x, t, y);
	      if (i < j) nba1[r] = SumInterp1D(n2[r], nba[r], x, t, y);
	    }
	  }
	  for (q = 0; q < n; q++) {
	    x[q] = ng[q];
	    fprintf(f2, "  %3d %3d %3d %3d %3d %12.5E %12.5E %12.5E %12.5E\n",
		    isym, i, j, -1, ng[q], z1[q], z2[q], nab1[q], nba1[q]);
	  }
	  fflush(f2);
	  a = SumInterp1D(n, z1, x, t, y);
	  heff[k0] += a;
	  na = SumInterp1D(n, nab1, x, t, y);
	  neff[k0] += na;
	  b = a;
	  nb = na;
	  if (i < j) {
	    b = SumInterp1D(n, z2, x, t, y);
	    heff[k1] += b;
	    nb = SumInterp1D(n, nba1, x, t, y);
	    neff[k1] += nb;
	  }
	  fprintf(f2, "# %3d %3d %3d %3d %3d %12.5E %12.5E %15.8E %12.5E %12.5E\n",
		  isym, pp, jj, i, j, a, b, mbpt[0].a, na, nb);
	}
      OUT:
	fprintf(f2, "# %3d %3d %3d %3d %3d %12.5E %12.5E %15.8E %12.5E %12.5E\n",
		isym, pp, jj, i, j,
		heff[k0]-mbpt[0].a, heff[k1]-mbpt[0].a, mbpt[0].a, 
		neff[k0], neff[k1]);
	fflush(f2);
      }
    }
    if (DiagnolizeHamilton(h) < 0) {
      printf("Diagnolizing Hamiltonian Error\n");
      ierr = -1;
      goto ERROR;
    }
    
    /* correct the normalization */
    y = h->mixing + h->dim;
    for (i = 0; i < h->dim; i++) {
      a = 0.0;
      for (j = 0; j < h->dim; j++) {
	for (k = 0; k < h->dim; k++) {
	  k0 = j*h->dim + k;
	  a += neff[k0] * y[j] * y[k];
	}
      }      
      b = sqrt(1.0 + a);
      fprintf(f2, "#NORM %3d %5d %12.5E %12.5E\n", isym, i, a, b);
      fflush(f2);
      for (j = 0; j < h->dim; j++) {
	y[j] /= b;
      }
      y += h->n_basis;
    }
    
    AddToLevels(h, nkg0, kg);
    free(heff);
    free(neff);
    h->heff = NULL;
  }

  SortLevels(nlevels, -1, 0);
  SaveLevels(fn, nlevels, -1);
  fclose(f2);
  for (m = 0; m < nf; m++) {
    fclose(f1[m]);    
  }

  if (mbpt_tr.mktr == 0) goto ERROR;
  for (m = 0; m < nf; m++) {
    sprintf(tfn1, "%s.tr", fn1[m]);
    f1[m] = fopen(tfn1, "r");
    if (f1[m] == NULL) {
      printf("no transition correction file %s\n", tfn1);
      goto ERROR;
    }
  }  
  ierr = ReadMBPT(nf, f1, mbpt, 3);
  if (ierr < 0) goto ERROR;  
  printf("MBPT Transition.\n");
  fflush(stdout);
  InitTransitionMBPT(&mtr, n0);
  CombineTransitionMBPT(nf, mbpt, mtr, n0, ng0);
  for (q = 0; q < n0; q++) {
    x[q] = ng0[q];
  }
  k = 2*mbpt_tr.mktr*MAX_SYMMETRIES;
  sprintf(tfn1, "%s.tr", fn2);
  f2 = fopen(tfn1, "w");
  for (j = 0; j < k; j++) {
    for (q = 0; q < mtr[j].nsym1; q++) {
      k0 = mtr[j].sym0->n_states * mtr[j].sym1[q]->n_states * mbpt_tr.naw;
      r = 0;
      for (m = 0; m < k0; m++) {
	for (k1 = 0; k1 < n0; k1++) {
	  z1[k1] = mtr[j].tma[q][m*n0+k1];
	  z2[k1] = mtr[j].rma[q][m*n0+k1];
	  fprintf(f2, "%5d %3d %5d %3d %5d %5d %5d %12.5E %12.5E\n", 
		  j, j/(2*mbpt_tr.mktr), q, mtr[j].isym1[q], m, k1, 
		  ng0[k1], z1[k1], z2[k1]);
	}
	mtr[j].tma[q][r] = SumInterp1D(n0, z1, x, t, y);
	mtr[j].rma[q][r] = SumInterp1D(n0, z2, x, t, y);
	fprintf(f2, "# %5d %3d %5d %3d %5d %12.5E %12.5E\n", 
		j, j/(2*mbpt_tr.mktr), q, mtr[j].isym1[q], m, 
		mtr[j].tma[q][r], mtr[j].rma[q][r]);
	r++;
      }
    }
  }
  fclose(f2);

  SetTransitionMode(0);
  SaveTransitionMBPT(mtr);
  AdjustAngularZ(mtr);

  for (m = 0; m < nf; m++) {
    FreeTransitionMBPT(mbpt[m].mtr);
  }
  FreeTransitionMBPT(mtr);

  for (m = 0; m < nf; m++) {
    fclose(f1[m]);
  }
  
 ERROR:
  for (m = 0; m < nf; m++) {
    free(mbpt[m].ng);
    free(mbpt[m].ng2);
    free(mbpt[m].hab1);
    free(mbpt[m].hba1);
    free(mbpt[m].hab);
    free(mbpt[m].hba);
  }
  free(f1);
  free(mbpt);
  free(dw);
  if (n > 0) {
    free(ng);
    for (i = 0; i < n; i++) {
      if (n2[i] > 0) {
	free(ng2[i]);
	free(hab[i]);
	free(hba[i]);
	free(nab[i]);
	free(nba[i]);
      }
    }
    free(ng2);
    free(hab);
    free(hba);
    free(nab);
    free(nba);
    free(n2);
  }
  if (n0 > 0) {
    free(ng0);
  }

  return ierr;
}
