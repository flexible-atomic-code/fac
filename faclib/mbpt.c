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

static double mbpt_freetr = 0;
static int mbpt_omp = 0;
static int mbpt_omp0 = 0;
static int mbpt_extra = 0;
static int mbpt_rand = 0;
static int mbpt_msort = 0;
static double mbpt_asort = 10.0;
static double mbpt_warn = 0.05;
static double mbpt_mwarn = 0.5;
static double mbpt_ignore = 50.0;
static double mbpt_warntr = 1.0;
static double mbpt_ignoretr = 10.0;
static int mbpt_ignorep = 0;
static long mbpt_ignoren = 0;
static double mbpt_angzc = 0.75;
static double mbpt_angzm = 0;
static int mbpt_adjaz = 0;
static int mbpt_savesum = 0;
static int mbpt_maxn = 0;
static int mbpt_maxm = 0;
static int mbpt_minn = 0;
static int mbpt_mini = 0;
static int mbpt_minn2 = 0;
static int mbpt_mini2 = 0;
static int mbpt_reinit_ncps = 0;
static double mbpt_reinit_mem = 0;
static int mbpt_nlev = 0;
static int *mbpt_ilev = NULL;
static double mbpt_mcut = 1e-4;
static double mbpt_mcut2 = 1e-1;
static double mbpt_mcut3 = 1e-1;
static double mbpt_mcut4 = 1.0;
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
static char mbpt_ccn[256] = "";
static ARRAY *mbpt_cca = NULL;
static double mbpt_wmix = 0.0025;
static double mbpt_nwmix = 0.0;
static struct {
  int nj;
  int *jp;
  IDXARY ibs;
  int *sa, *ma;
  double *sd1a, *sd2a, *c12a;
} mbptjp = {0, NULL};

static int mbpt_nmk0, mbpt_nmk1;
static ARRAY *mbpt_csary;

#pragma omp threadprivate(mbpt_cs, mbpt_cfg, mbpt_bas0, mbpt_bas0s, mbpt_bas0d, mbpt_bas1, mbpt_ibas0, mbpt_ibas1, mbptjp, mbpt_cca)
  
static TR_OPT mbpt_tr;

void PrintMBPTOptions(void) {
  printf("omp=%d\n", mbpt_omp0);
  printf("extra=%d\n", mbpt_extra);
  printf("rand=%d\n", mbpt_rand);
  printf("msort=%d\n", mbpt_msort);
  printf("asort=%g\n", mbpt_asort);
  printf("warn=%g/%g/%g\n", mbpt_warn, mbpt_mwarn, mbpt_warntr);
  printf("ignore=%g/%g\n", mbpt_ignore, mbpt_ignoretr);
  printf("ignorep=%d\n", mbpt_ignorep);
  printf("angzc=%g\n", mbpt_angzc);
  printf("angzm=%g\n", mbpt_angzm);
  printf("adjaz=%d\n", mbpt_adjaz);
  printf("savesum=%d\n", mbpt_savesum);
  printf("maxn=%d\n", mbpt_maxn);
  printf("maxm=%d\n", mbpt_maxm);
  printf("minn=%d\n", mbpt_minn);
  printf("mini=%d\n", mbpt_mini);
  printf("minn2=%d\n", mbpt_minn2);
  printf("mini2=%d\n", mbpt_mini2);
  printf("reinit_ncps=%d\n", mbpt_reinit_ncps);
  printf("reinit_meme=%g\n", mbpt_reinit_mem);
  printf("nlev=%d\n", mbpt_nlev);
  printf("mcut=%g %g %g %g\n", mbpt_mcut, mbpt_mcut2, mbpt_mcut3, mbpt_mcut4);
  printf("n3=%d\n", mbpt_n3);
  printf("nsplit=%d\n", mbpt_nsplit);
  printf("freetr=%g\n", mbpt_freetr);
  printf("ccn=%s\n", mbpt_ccn);
  printf("wmix=%g\n", mbpt_wmix);
  printf("nwmix=%g\n", mbpt_nwmix);
}

int SkipMPIM(int s) {
  if (s == mbpt_omp) {
    return SkipMPI();
  } else {
    return 0;
  }
}
    
void InitMBPT(void) {
  mbpt_tr.nktr = 0;
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

void TransitionMBPT(int nk, int *ks, int n) {
  if (mbpt_tr.nktr > 0) free(mbpt_tr.ks);
  if (nk <= 0) {
    mbpt_tr.nktr = 0;
    mbpt_tr.ks = NULL;
  }
  if (n < 0) {
    if (nk > 0) n = 3;
    else n = 0;
  }
  if (n > 0 && n != mbpt_tr.naw) {
    if (mbpt_tr.naw > 0) free(mbpt_tr.awgrid);
    mbpt_tr.awgrid = malloc(sizeof(double)*n);
  }
  if (nk > 0 && ks == NULL) {
    mbpt_tr.nktr = nk;
    mbpt_tr.ks = malloc(sizeof(int)*mbpt_tr.nktr);
    int i, j = 0;
    for (i = 1; i <= nk; i++) {
      mbpt_tr.ks[j++] = i;
    }
  } else {
    mbpt_tr.nktr = nk;
    mbpt_tr.ks = malloc(sizeof(int)*mbpt_tr.nktr);
    memcpy(mbpt_tr.ks, ks, sizeof(int)*mbpt_tr.nktr);
  }
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

void SetOptionMBPT(char *s, char *sp, int ip, double dp) {
  if (0 == strcmp(s, "mbpt:ignorep")) {
    mbpt_ignorep = ip;
    return;
  }
  if (0 == strcmp(s, "mbpt:warntr")) {
    mbpt_warntr = dp;
    return;
  }
  if (0 == strcmp(s, "mbpt:ignoretr")) {
    mbpt_ignoretr = dp;
    return;
  }
  if (0 == strcmp(s, "mbpt:warn")) {
    mbpt_warn = dp;
    return;
  }
  if (0 == strcmp(s, "mbpt:mwarn")) {
    mbpt_mwarn = dp;
    return;
  }
  if (0 == strcmp(s, "mbpt:ccn")) {
    strncpy(mbpt_ccn, sp, 256);
    return;
  }
  if (0 == strcmp(s, "mbpt:ignore")) {
    mbpt_ignore = dp;
    return;
  }
  if (0 == strcmp(s, "mbpt:angzc")) {
    mbpt_angzc = dp;
    return;
  }  
  if (0 == strcmp(s, "mbpt:angzm")) {
    mbpt_angzm = dp;
    return;
  }  
  if (0 == strcmp(s, "mbpt:adjaz")) {
    mbpt_adjaz = ip;
    return;
  }  
  if (0 == strcmp(s, "mbpt:rand")) {
    mbpt_rand = ip;
    return;
  }  
  if (0 == strcmp(s, "mbpt:msort")) {
    mbpt_msort = ip;
    return;
  }  
  if (0 == strcmp(s, "mbpt:asort")) {
    mbpt_asort = dp;
    return;
  }  
  if (0 == strcmp(s, "mbpt:n3")) {
    mbpt_n3 = ip;
    return;
  }  
  if (0 == strcmp(s, "mbpt:mcut2")) {
    mbpt_mcut2 = dp;
    return;
  }  
  if (0 == strcmp(s, "mbpt:mcut3")) {
    mbpt_mcut3 = dp;
    return;
  }  
  if (0 == strcmp(s, "mbpt:mcut4")) {
    mbpt_mcut4 = dp;
    return;
  }  
  if (0 == strcmp(s, "mbpt:wmix")) {
    mbpt_wmix = dp;
    return;
  }  
  if (0 == strcmp(s, "mbpt:nwmix")) {
    mbpt_nwmix = dp;
    return;
  }  
  if (0 == strcmp(s, "mbpt:omp")) {
    mbpt_omp0 = ip;
    return;
  }  
  if (0 == strcmp(s, "mbpt:extra")) {
    mbpt_extra = ip;
    return;
  }  
  if (0 == strcmp(s, "mbpt:savesum")) {
    mbpt_savesum = ip;
    return;
  }  
  if (0 == strcmp(s, "mbpt:minn")) {
    mbpt_minn = ip;
    return;
  }  
  if (0 == strcmp(s, "mbpt:reinit_ncps")) {
    mbpt_reinit_ncps = ip;
    return;
  }  
  if (0 == strcmp(s, "mbpt:reinit_mem")) {
    mbpt_reinit_mem = 1e9*ip;
    return;
  }
  if (0 == strcmp(s, "mbpt:freetr")) {
    mbpt_freetr = 1e9*dp;
  }
  return;
}

double AngZCutMBPT(void) {
  return mbpt_angzc;
}

void SetWarnMBPT(double f, double g) {
  mbpt_warn = f;
  mbpt_ignore = g;
}

void SetOptMBPT(int nr, int n3, double c, double d, double e, double f) {
  mbpt_rand = nr%10;
  mbpt_msort = nr/10;
  mbpt_n3 = n3;
  if (c >= 0) mbpt_mcut = c;
  else mbpt_mcut = 1e-3;
  if (d >= 0) mbpt_mcut2 = d;
  else mbpt_mcut2 = 1e-2;
  if (e >= 0) mbpt_mcut3 = e;
  else mbpt_mcut3 = 1e-1;
  if (f >= 0) mbpt_mcut4 = f;
  else mbpt_mcut4 = 1.0;
}

void SetExtraMBPT(int m) {
  int am = abs(m);
  mbpt_omp0 = am/1000000;
  am = am%1000000;
  mbpt_extra = am%10;
  mbpt_savesum = mbpt_extra/5;
  mbpt_extra = mbpt_extra%5;
  am /= 10;
  mbpt_minn = am%1000;
  am /= 1000;
  if (am == 0) {
    mbpt_reinit_ncps = 0;
    mbpt_reinit_mem = 0;
  } else if (m < 0) {
    mbpt_reinit_mem = 0;
    mbpt_reinit_ncps = am;
  } else {
    mbpt_reinit_ncps = 0;
    mbpt_reinit_mem = 1e9*am;
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
      mbpt_se[k] = nd;
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
    if (nk < 0) {      
      AllocHamMem(ha, -1, -1);
      AllocHamMem(ha, 0, 0);
      continue;
    }
    mb.nbasis = ha->dim;
    mb.isym = i;
    mb.basis = malloc(sizeof(int)*mb.nbasis);
    mb.ene = malloc(sizeof(double)*mb.nbasis);
    memcpy(mb.basis, ha->basis, sizeof(int)*mb.nbasis);
    memcpy(mb.ene, ha->hamilton, sizeof(double)*mb.nbasis);
    mb.bmax = mb.basis[mb.nbasis-1];
    ArrayAppend(&base, &mb, NULL);
    AllocHamMem(ha, -1, -1);
    AllocHamMem(ha, 0, 0);
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
	  ha = GetHamilton(i);
	  nk = ConstructHamiltonDiagonal(mbp->isym, 1, &kgp, 0);
	  if (nk < 0) {
	    AllocHamMem(ha, -1, -1);
	    AllocHamMem(ha, 0, 0);
	    continue;
	  }
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
	  AllocHamMem(ha, -1, -1);
	  AllocHamMem(ha, 0, 0);
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
  ArrayFreeLock(&ccfg, NULL);
  for (i = 0; i < base.dim; i++) {
    mbp = ArrayGet(&base, i);
    free(mbp->basis);
    free(mbp->ene);
  }
  ArrayFreeLock(&base, NULL);

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
	if (mbpt_maxm > 0 && mbpt_maxn > 0 && ng[i] > mbpt_maxn) continue;
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
      orb = GetOrbitalSolvedNoLock(i);
      QED1E(i, i);
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
    int i0, i1;
    if (c->shells[0].nq > 1) {
      if (c->shells[0].n > cs[nc]->nnrs) return -1;
      i0 = ShellToInt(c->shells[0].n, c->shells[0].kappa);
      i1 = i0;
    } else if (c->n_shells > 1) {
      if (c->shells[1].n > cs[nc]->nnrs) return -1;
      i0 = ShellToInt(c->shells[0].n, c->shells[0].kappa);
      i1 = ShellToInt(c->shells[1].n, c->shells[1].kappa);      
    } else {
      i0 = ShellToInt(c->shells[0].n, c->shells[0].kappa);
      i1 = mbpt_nmk1-1;
    }
    m = i0*mbpt_nmk1 + i1;
    ARRAY *csa = &mbpt_csary[m];
    int *ic;
    for (i = 0; i < csa->dim; i++) {
      ic = ArrayGet(csa, i);
      if (k != cs[*ic]->n_shells) continue;
      m = memcmp(c->shells, cs[*ic]->shells, sizeof(SHELL)*k);
      if (m == 0) return i;
    }
  }
  return -1;
}

double SumInterp1D0(int n, double *z, double *x, double *t, double *y, int m) {
  int i, k0, k;
  double a, b, r;
  r = 0.0;
  for (i = 0; i < n; i++) {
    r += z[i];
  }
  if (isnan(r) || isinf(r)) {
    MPrintf(-1, "invalid value SumInterp1D a: %g\n", r);
    //Abort(1);
  }
  if (1+r == 1) return 0.0;
  if (n == 1) return r;
  for (k0 = 1; k0 < n; k0++) {
    if (x[k0]-x[k0-1] > 1) break;
  }
  if (k0 == n) return r;
  spline_work(x, z, n, 1e31, 1e31, y, t);
  for (i = k0; i < n; i++) {
    for (a = x[i-1]+1.0; a < x[i]; a += 1.0) {
      k = splint(x, z, y, n, a, &b);
      if (k < 0) {
	MPrintf(-1, "SumInterp1D error: %d %d\n", k, m);
	Abort(1);
      }
      r += b;
    }    
  }
  return r;
}

double SumInterp1D(int n0, double *z, double *x, double *t, double *y, int m) {
  int i, k0, k1, nk, j, nw, n;
  double r, a, b, c, d, e, f, g, h;
  double p1, p2, p3, q1, q2, q3;
  double *dw0, *dw1;

  for (n = n0-1; n >= 0; n--) {
    if (z[n]) {
      n++;
      break;
    }
  }
  if (n == 0) return 0.0;
  r = 0.0;
  for (i = 0; i < n; i++) {
    r += z[i];
  }
  if (1+r == 1) return 0.0;
  if (n == 1) {
    if (isnan(r) || isinf(r)) {
      MPrintf(-1, "invalid value SumInterp1D a: %d %d %d %g\n", m, n, n0, r);
      for (i = 0; i < n0; i++) {
	MPrintf(-1, "%d %g %g\n", i, x[i], z[i]);
      }
      fflush(stdout);
      //Abort(1);
    }
    return r;
  }
  for (k0 = 1; k0 < n; k0++) {
    if (x[k0]-x[k0-1] > 1) break;
  }
  if (k0 == n) {
    k1 = n-3;
    if (k1 < 0) k1 = 0;
    for (i = k1; i < n; i++) {
      t[i] = log(x[i]);
      y[i] = log(fabs(z[i]));
      if (1+z[i] == 1) break;
    }
    if (i < n) nk = 1;
    else nk = n-k1;
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
  if (k1 == n) k1 = n-1;
  if (k1 < k0) k1 = k0;  
  for (i = k0; i < n; i++) {
    t[i] = log(x[i]);
  }
  dw0 = WorkSpace();
  nw = 2 + (int)(x[n-1]-x[0]);
  dw1 = dw0 + nw;
  if (k1 > k0) {
    nk = k1 - k0 + 1;
    j = 0;
    for (i = k0; i < k1; i++) {    
      for (a = x[i]+1; a < x[i+1]; a += 1.0) {
	dw0[j++] = log(a);
      }
    }
    if (j > 0) {
      UVIP3P(3, nk, t+k0, z+k0, j, dw0, dw1);
      for (i = 0; i < j; i++) {
	r += dw1[i];
      }
    }
  }
  for (i = k1; i < n; i++) {
    y[i] = log(fabs(z[i]));
  }
  nk = n - k1;
  if (nk > 1) {
    j = 0;
    for (i = k1+1; i < n; i++) {
      for (a = x[i-1]+1; a < x[i]; a += 1.0) {
	dw0[j++] = log(a);
      }
    }
    if (j > 0) {
      UVIP3P(3, nk, t+k1, y+k1, j, dw0, dw1);
      for (i = 0; i < j; i++) {
	b = exp(dw1[i]);
	if (z[n-1] < 0) b = -b;
	r += b;
      }
    }
  }
 END:
  h = 0.0;
  a = 0.0;
  b = 0.0;
  c = 0.0;
  if (mbpt_extra > 0) {
    if (nk > 2 && mbpt_extra <= 2) {
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
      if (c <= 0 && 2.0*log(x[i])*c+b < -1.1) {
	g = 1.0;
	d = x[i]+1.0;
	int iter = 0;
	while (g > EPS3) {
	  e = log(d);
	  f = exp(a + b*e + c*e*e);
	  if (z[i] < 0) f = -f;
	  h += f;
	  g = f*d/(-(2.0*e*c+b)-1.0);
	  g = fabs(g/r);
	  d += 1.0;
	  iter++;
	  if (iter >= 10000) {
	    MPrintf(-1, "maxiter reached in mbpt extrapolate: %d %d %d %d %g %g %g %g %g\n",
		    iter, m, i, nk, g, f, c, b, x[i]);
	    if (iter >= 10010) {
	      g = 0.0;
	      f = 0.0;
	      for (j = 0; j < n; j++) {
		MPrintf(-1, "%3d %12.5E %12.5E\n", j, x[j], z[j]);
	      }
	      MPrintf(-1, "# %12.5E %12.5E %12.5E %12.5E %12.5E\n", a, b, c, h, r);
	      break;
	    }
	  }
	}	
	g *= fabs(r);
	if (f < 0) g = -g;
	h += g-f;	
      }
    }
    if (h == 0 && nk > 2) {
      k0 = n-2;
      k1 = n-1;
      a = y[k1] - y[k0];
      b = t[k1] - t[k0];
      a = -a/b;      
      if (a > 1.1) {
	b = z[k1]*pow(x[k1]/(1.0+x[k1]), a);
	a -= 1.0;
	h = b*(1.0+x[k1])/a;
      }
      b = 0.0;
      c = 0.0;
    }
    if (mbpt_extra == 2 || mbpt_extra == 4) {
      for (i = 0; i < n; i++) {
	MPrintf(-1, "%3d %12.5E %12.5E\n", i, x[i], z[i]);
      }
      MPrintf(-1, "#a %4d %12.5E %12.5E %12.5E %12.5E %12.5E\n",
	      m, a, b, c, h, r);
    } else if (fabs(h) > 5*fabs(r)) {
      for (i = 0; i < n; i++) {
	MPrintf(-1, "%3d %12.5E %12.5E\n", i, x[i], z[i]);
      }
      MPrintf(-1, "#b %4d %12.5E %12.5E %12.5E %12.5E %12.5E\n",
	      m, a, b, c, h, r);
    }
    r += h;
  }
  if (isnan(r) || isinf(r)) {
    MPrintf(-1, "invalid value SumInterp1D b: %d %d %d %d %g %g\n",
	    m, n, n0, nk, r, h);
    for (i = 0; i < n0; i++) {
      MPrintf(-1, "%d %g %g\n", i, x[i], z[i]);
    }
    fflush(stdout);
    //Abort(1);
  }
  return r;
}

double SumInterpH(int n, int *ng, int n2, int *ng2, 
		  double *h, double *h1, double *w, int md) {
  double r, *p, *x, *y, *z;
  int m, i0, i1;

  i0 = n+n2;
  x = w + i0;
  y = x + i0;
  z = y + i0;
  if (n == 1) {
    r = h1[0];
  } else {
    for (i0 = 0; i0 < n; i0++) {
      x[i0] = ng[i0];
    }
    r = SumInterp1D(n, h1, x, w, y, md);
  }
  if (h) {
    p = h;
    for (i0 = 0; i0 < n; i0++) {
      for (i1 = 0; i1 < n2; i1++) {
	x[i1] = ng2[i1] + ng[i0];
      }
      z[i0] = SumInterp1D(n2, p, x, w, y, md+ng[i0]*10);
      p += n2;
    }  
    for (i0 = 0; i0 < n; i0++) {
      x[i0] = ng[i0];
    }
    r += SumInterp1D(n, z, x, w, y, -md);
  }

  return r;
}

void UnpackSymStateTR(MBPT_TR *mtr, int ms, int ip, int *s, int *m, int *is) {
  UnpackSymState(ms, s, m);
  *is = 2*mbpt_tr.nktr*(*s) + 2*ip;
  if (mtr[*is].idb[*s]) {
    *m = IdxGet(mtr[*is].idb[*s], *m);
  } else {
    *m = -5;
  }
}

void UnpackSymStateMBPT(MBPT_EFF **meff, int ms, int *s, int *m) {
  UnpackSymState(ms, s, m);
  if (meff[*s] && meff[*s]->idb) {
    *m = IdxGet(meff[*s]->idb, *m);
  } else {
    *m = -5;
  }
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
#pragma omp atomic
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
#pragma omp atomic
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
	  //ms0 = c0->symstate[q0];
	  //UnpackSymState(ms0, &s0, &m0);
	  //ms1 = c1->symstate[q1];
	  //UnpackSymState(ms1, &s1, &m1);	  
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
  int nn1 = orb->n;
  int ka1 = orb->kappa;
  ORBITAL *o1 = orb;
  orb = GetOrbital(ks2[3]);
  d1 += orb->energy + orb->qed;
  int nn2 = orb->n;
  int ka2 = orb->kappa;
  ORBITAL *o2 = orb;
  orb = GetOrbital(ks2[0]);
  d1 -= orb->energy + orb->qed;
  int nn3 = orb->n;
  int ka3 = orb->kappa;
  ORBITAL *o3 = orb;
  orb = GetOrbital(ks2[1]);
  d1 -= orb->energy + orb->qed;
  int nn4 = orb->n;
  int ka4 = orb->kappa;
  ORBITAL *o4 = orb;
  orb = GetOrbital(ks1[0]);
  d2 = orb->energy + orb->qed;
  int nn5 = orb->n;
  int ka5 = orb->kappa;
  ORBITAL *o5 = orb;
  orb = GetOrbital(ks1[1]);
  d2 += orb->energy + orb->qed;
  int nn6 = orb->n;
  int ka6 = orb->kappa;
  ORBITAL *o6 = orb;
  orb = GetOrbital(ks1[2]);
  d2 -= orb->energy + orb->qed;
  int nn7 = orb->n;
  int ka7 = orb->kappa;
  ORBITAL *o7 = orb;
  orb = GetOrbital(ks1[3]);
  d2 -= orb->energy + orb->qed;
  int nn8 = orb->n;
  int ka8 = orb->kappa;
  ORBITAL *o8 = orb;
  /*
  d1 = GetOrbital(ks2[2])->energy + GetOrbital(ks2[3])->energy;
  d1 -= GetOrbital(ks2[0])->energy + GetOrbital(ks2[1])->energy;
  d2 = GetOrbital(ks1[0])->energy + GetOrbital(ks1[1])->energy;
  d2 -= GetOrbital(ks1[2])->energy + GetOrbital(ks1[3])->energy;
  */
  int *sp, *mp;
  double *sd1p, *sd2p, *c12p;
  int warned = 0, ignored = 0;
  sd1p = mbptjp.sd1a;
  sd2p = mbptjp.sd2a;
  c12p = mbptjp.c12a;
  sp = mbptjp.sa;
  mp = mbptjp.ma;
  for (k = 0; k < mst; k++) {
    q0 = bst[k];
    q1 = kst[k];
    ms0 = c0->symstate[q0];
    UnpackSymStateMBPT(meff, ms0, &s0, &m0);
    ms1 = c1->symstate[q1];
    UnpackSymStateMBPT(meff, ms1, &s1, &m1);
    sp[k] = -1;
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
    if (m1 > m0) {
      m = m1*(m1+1)/2 + m0;
    } else {
      m = m0*(m0+1)/2 + m1;
    }
    mp[k] = m;
    sp[k] = s0;
    double c12 = c/(d1*d2);
    sd1 = c/d1;
    sd2 = c/d2;
    sd1p[k] = sd1;
    sd2p[k] = sd2;
    c12p[k] = c12;
    double sd1w = sqrt(fabs(sd1/d1));
    double sd2w = sqrt(fabs(sd2/d2));
    double sd1s = sd1w;
    double sd2s = sd2w;
    if (isnan(c12)) {
      sd1s = 1E50;
      sd2s = 1E50;
    } else if (isinf(c12)) {
      sd1s = 2E50;
      sd2s = 2E50;
    }
    if (meff[s0]->imbpt[m] > 1 && !warned) {     
      double wth = (mbpt_nwmix-meff[s0]->imbpt[m]);
      wth /= mbpt_nwmix - 2.0;
      wth = mbpt_warn + mbpt_mwarn*wth;
      if (wth > 0 && (sd1w > wth || sd2w > wth)) {
	int na[4], ka[4];
	na[0] = nn1;
	na[1] = nn3;
	na[2] = nn2;
	na[3] = nn4;
	ka[0] = ka1;
	ka[1] = ka3;
	ka[2] = ka2;
	ka[3] = ka4;
	CONFIG *cr = ExciteNRConfig(c1, 4, na, ka);
	char scr[2048] = "NULL";
	if (mbpt_cca) {
	  cr->sweight = d1;
	  cr->energy = d2;
	  cr->delta = sd1w;
	  cr->sth = sd2w;
	  cr->cth = sd1s;
	  cr->mde = sd2s;
	  cr->icfg = 2000+meff[s0]->imbpt[m];
	  cr->n_electrons = c0->igroup;
	  cr->n_csfs = c0->icfg;
	  cr->nnrs = c1->igroup;
	  cr->igroup = c1->icfg;
	  ArrayAppend(mbpt_cca, &cr, NULL);
	  warned = 1;
	} else {
	  if (cr != NULL) {
	    ConstructNRConfigName(scr, 2048, cr);
	    FreeConfigData(cr);
	    free(cr);
	  }
	  printf("large h22term warn:  %s %d %d %d %d %d %d %d %g %g %g %g %g\n",
		 scr, c0->igroup, c0->icfg, c1->igroup, c1->icfg, s0, m0, m1,
		 sd1w, sd2w, d1, d2, wth);
	  warned = 2;
	}      
      }
    }
    if (mbpt_ignore > 0) {// && (c0 != c1 || m0 != m1)) {
      if (sd1s > mbpt_ignore || sd2s > mbpt_ignore) {
	if (warned != 2 && (mbpt_ignorep&1)) {
	  int na[4], ka[4];
	  na[0] = nn1;
	  na[1] = nn3;
	  na[2] = nn2;
	  na[3] = nn4;
	  ka[0] = ka1;
	  ka[1] = ka3;
	  ka[2] = ka2;
	  ka[3] = ka4;
	  CONFIG *cr = ExciteConfig(c1, 4, na, ka);
	  char scr[2048] = "NULL";
	  if (cr != NULL) {
	    ConstructConfigName(scr, 2048, cr);
	    FreeConfigData(cr);
	    free(cr);
	  }	  
	  printf("large h22term ignore: %s %d %d %d %d %d %d %d %g %g %g %g %g %d\n",
		 scr,
		 c0->igroup, c0->icfg, c1->igroup, c1->icfg, s0, m0, m1,
		 sd1s, sd2s, d1, d2, mbpt_ignore, meff[s0]->imbpt[m]);
	}
	sd1 = 0;
	sd2 = 0;
	c12 = 0;
	ignored = 1;
#pragma omp atomic
	mbpt_ignoren++;
	break;
      }
    }
  }
  if (!ignored) {
    for (k = 0; k < mst; k++) {
      if (sp[k] < 0) continue;
      s0 = sp[k];
      m = mp[k];    
      if (i0 >= 0) {
	h1 = meff[s0]->hab[m];
	h2 = meff[s0]->hba[m];
	ng = meff[s0]->n * meff[s0]->n2;
	i1 = i0;
      } else {
	h1 = meff[s0]->hab1[m];
	h2 = meff[s0]->hba1[m];
	ng = meff[s0]->n;
	i1 = -(i0+1);
      }
#pragma omp atomic
      h1[i1] += sd1p[k];
#pragma omp atomic
      h2[i1] += sd2p[k];
      i1g = i1+ng;
#pragma omp atomic
      h1[i1g] += c12p[k];
#pragma omp atomic
      h2[i1g] += c12p[k];
    }
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
  int nn1 = orb->n;
  int ka1 = orb->kappa;
  orb = GetOrbital(ks[3]);
  d1 += orb->energy + orb->qed;
  int nn2 = orb->n;
  int ka2 = orb->kappa;
  orb = GetOrbital(ks[0]);
  d1 -= orb->energy + orb->qed;
  int nn3 = orb->n;
  int ka3 = orb->kappa;
  orb = GetOrbital(ks[1]);
  d1 -= orb->energy + orb->qed;
  int nn4 = orb->n;
  int ka4 = orb->kappa;
  orb = GetOrbital(k0);
  d2 = orb->energy + orb->qed;
  int nn5 = orb->n;
  int ka5 = orb->kappa;
  orb = GetOrbital(k1);
  d2 -= orb->energy + orb->qed;
  int nn6 = orb->n;
  int ka6 = orb->kappa;
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
  int *sp, *mp;
  double *sd1p, *sd2p, *c12p;
  int warned = 0, ignored = 0;
  sd1p = mbptjp.sd1a;
  sd2p = mbptjp.sd2a;
  c12p = mbptjp.c12a;
  sp = mbptjp.sa;
  mp = mbptjp.ma;
  for (k = 0; k < mst; k++) {
    sp[k] = -1;
    q0 = bst[k];
    q1 = kst[k];
    ms0 = c0->symstate[q0];
    UnpackSymStateMBPT(meff, ms0, &s0, &m0);
    ms1 = c1->symstate[q1];
    UnpackSymStateMBPT(meff, ms1, &s1, &m1);
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
    if (m0 <= m1) {
      m = m1*(m1+1)/2 + m0;
      mp[k] = -m;
    } else {
      m = m0*(m0+1)/2 + m1;
      mp[k] = m;
    }
    sp[k] = s0;
    double c12 = c/(d1*d2);
    sd = c/d1;
    se = c/d2;
    sd1p[k] = sd;
    sd2p[k] = se;
    c12p[k] = c12;
    double sd1w = sqrt(fabs(sd/d1));
    double sd2w = sqrt(fabs(se/d2));
    double sd1s = sd1w;
    double sd2s = sd2w;//c0->cth>0?c0->cth/fabs(d2):0;
    if (isnan(c12)) {
      sd1s = 1E50;
      sd2s = 1E50;
    } else if (isinf(c12)) {
      sd1s = 2E50;
      sd2s = 2E50;
    }
    //if (c0->icfg >= 0 || c1->icfg >= 0) {
    if (meff[s0]->imbpt[m] > 1 && !warned) {
      double wth = (mbpt_nwmix-meff[s0]->imbpt[m]);
      wth /= mbpt_nwmix - 2.0;
      wth = mbpt_warn + mbpt_mwarn*wth;
      if (wth > 0 && (sd1w > wth || sd2w > wth)) {
	int na[2], ka[2];
	na[0] = nn5;
	na[1] = nn6;
	ka[0] = ka5;
	ka[1] = ka6;
	CONFIG *cr = ExciteNRConfig(c0, 2, na, ka);
	char scr[2048] = "NULL";
	if (mbpt_cca) {
	  cr->sweight = d1;
	  cr->energy = d2;
	  cr->delta = sd1w;
	  cr->sth = sd2w;
	  cr->cth = sd1s;
	  cr->mde = sd2s;
	  cr->icfg = 1000 + meff[s0]->imbpt[m];
	  cr->n_electrons = c0->igroup;
	  cr->n_csfs = c0->icfg;
	  cr->nnrs = c1->igroup;
	  cr->igroup = c1->icfg;
	  ArrayAppend(mbpt_cca, &cr, NULL);
	  warned = 1;
	} else {
	  if (cr != NULL) {
	    ConstructNRConfigName(scr, 2048, cr);
	    FreeConfigData(cr);
	    free(cr);
	  }	
	  printf("large h12term warn: %s %d %d %d %d %d %d %d %g %g %g %g %g\n",
		 scr,
		 c0->igroup, c0->icfg, c1->igroup, c1->icfg, s0, m0, m1,
		 sd1s, sd2s, d1, d2, wth);
	  warned = 2;
	}
      }
    }
    if (mbpt_ignore > 0) {// && (c0 != c1 || m0 != m1)) {
      if (sd1s > mbpt_ignore || sd2s > mbpt_ignore) {
	if (warned != 2 && (mbpt_ignorep&1)) {
	  int na[2], ka[2];
	  na[0] = nn5;
	  na[1] = nn6;
	  ka[0] = ka5;
	  ka[1] = ka6;
	  CONFIG *cr = ExciteConfig(c0, 2, na, ka);
	  char scr[2048] = "NULL";
	  if (cr != NULL) {
	    ConstructConfigName(scr, 2048, cr);
	    FreeConfigData(cr);
	    free(cr);
	  }	
	  printf("large h12term ignore: %s %d %d %d %d %d %d %d %g %g %g %g %g %d\n",
		 scr,
		 c0->igroup, c0->icfg, c1->igroup, c1->icfg, s0, m0, m1,
		 sd1s, sd2s, d1, d2, mbpt_ignore, meff[s0]->imbpt[m]);
	}
	sd = 0;
	se = 0;
	c12 = 0;
	ignored = 1;
#pragma omp atomic
	mbpt_ignoren++;
	break;
      }
    }
  }
  if (!ignored) {
    for (k = 0; k < mst; k++) {
      if (sp[k] < 0) continue;
      s0 = sp[k];
      if (mp[k] <= 0) {
	m = -mp[k];
	h1 = meff[s0]->hab1[m];
	h2 = meff[s0]->hba1[m];
      } else {
	m = mp[k];
	h1 = meff[s0]->hba1[m];
	h2 = meff[s0]->hab1[m];
      }
      ng = meff[s0]->n;
#pragma omp atomic
      h1[i0] += sd1p[k];
#pragma omp atomic
      h2[i0] += sd2p[k];
      //c /= d1*d2;
      i0g = i0 + ng;
#pragma omp atomic
      h1[i0g] += c12p[k];
#pragma omp atomic
      h2[i0g] += c12p[k];
    }
  }
}

int TriangleTR(int j0, int j1) {
  int i;
  for (i = 0; i < mbpt_tr.nktr; i++) {
    if (Triangle(j0, j1, 2*mbpt_tr.ks[i])) break;
  }
  if (i == mbpt_tr.nktr) return 0;
  return 1;
}

void TR12Term(MBPT_TR *mtr, CONFIG *c0, CONFIG *c1,
	      int ns, SHELL *bra, SHELL *ket, 
	      SHELL_STATE *sbra, SHELL_STATE *sket,
	      int mst, int *bst, int *kst,
	      INTERACT_SHELL *s, int ph, int *ks, int k0, int k1,
	      FORMULA *fm, double **a, int i0, int ng) {  
  int kk, kk2, kmin, kmax, gauge, n, p1, j0, j1;
  int q0, q1, k, m0, m1, s0, s1, m, ms0, ms1, is0, is1, p, ip, p0;  
  double c, *r1, sd, se, yk[MKK], d2, rt;
  ORBITAL *orb;
  MBPT_PMA *pma;

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
  if (!TriangleTR(s[0].j, s[1].j)) return;
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
    for (ip = 0; ip < mbpt_tr.nktr; ip++) {
      p = mbpt_tr.ks[ip];
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
	  q0 = k*mbpt_tr.nktr + ip;
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
  int nn1 = orb->n;
  int ka1 = orb->kappa;
  orb = GetOrbital(ks[3]);
  d2 += orb->energy + orb->qed;
  int nn2 = orb->n;
  int ka2 = orb->kappa;
  orb = GetOrbital(ks[0]);
  d2 -= orb->energy + orb->qed;
  int nn3 = orb->n;
  int ka3 = orb->kappa;
  orb = GetOrbital(ks[1]);
  d2 -= orb->energy + orb->qed;
  int nn4 = orb->n;
  int ka4 = orb->kappa;
  double sd1s = c1->cth/fabs(d2);
  int warned = 0;
  if (c0->icfg >= 0 && c1->icfg >= 0) {
    if (mbpt_warntr > 0 && sd1s > mbpt_warntr) {
	char s1[16], s2[16], s3[16], s4[16];
	ShellString(nn1, ka1, -1, s1);
	ShellString(nn2, ka2, -1, s2);
	ShellString(nn3, ka3, -1, s3);
	ShellString(nn4, ka4, -1, s4);
	printf("large t12term warn: %s %s %s %s %d %d %d %d %g %g %g\n",
	       s1, s2, s3, s4,
	       c0->igroup, c0->icfg, c1->igroup, c1->icfg,
	       sd1s, d2, mbpt_warntr);
	warned = 1;
    }
  }
  if (mbpt_ignoretr > 0 && sd1s > mbpt_ignoretr) {
    if (!warned && (mbpt_ignorep&2)) {
      //char sc[64];
      char s1[16], s2[16], s3[16], s4[16];
      ShellString(nn1, ka1, -1, s1);
      ShellString(nn2, ka2, -1, s2);
      ShellString(nn3, ka3, -1, s3);
      ShellString(nn4, ka4, -1, s4);
      printf("large t12term ignore: %s %s %s %s %d %d %d %d %g %g %g\n",
	     s1, s2, s3, s4,
	     c0->igroup, c0->icfg, c1->igroup, c1->icfg,
	     sd1s, d2, mbpt_ignoretr);
    }
#pragma omp atomic
    mbpt_ignoren++;
    return;
  }

  gauge = GetTransitionGauge();
  ng -= mbpt_mini;
  for (kk = kmin; kk <= kmax; kk += 2) {
    kk2 = kk/2;
    SlaterTotal(&sd, &se, NULL, ks, kk, 0);  
    yk[kk2] = sd + se;
  }
  for (ip = 0; ip < mbpt_tr.nktr; ip++) {
    p = mbpt_tr.ks[ip];
    if (!Triangle(s[0].j, s[1].j, 2*p)) continue;
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
      UnpackSymStateTR(mtr, ms0, ip, &s0, &m0, &is0);
      ms1 = c1->symstate[q1];
      UnpackSymStateTR(mtr, ms1, ip, &s1, &m1, &is1);
      //is0 = 2*mbpt_tr.nktr*s0 + 2*ip;
      if (IsOdd(p0)) is0++;
      //is1 = IBisect(s1, mtr[is0].nsym1, mtr[is0].isym1);
      is1 = IdxGet(mtr[is0].ids, s1);
      if (is1 < 0) continue;
      int m01 = m0*mtr[is0].ndim1[is1] + m1;
      pma = mtr[is0].pma[is1][m01];
      if (pma == NULL) continue;
      q0 = k*mbpt_tr.nktr + ip;
      c = 0.0;
      for (kk = kmin; kk <= kmax; kk += 2) {
	kk2 = kk/2;
	c += a[kk2][q0]*yk[kk2];
      }
      for (p1 = 0; p1 < n; p1++) {
	q1 = p1*ng + i0;
	rt = c*r1[p1]/d2;
#pragma omp atomic
	pma->tma[q1] += rt;
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
  for (ip = 0; ip < mbpt_tr.nktr; ip++) {
    p = mbpt_tr.ks[ip];
    if (!Triangle(s[0].j, s[1].j, 2*p)) continue;
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
      UnpackSymStateTR(mtr, ms0, ip, &s0, &m0, &is0);
      DecodePJ(s0, NULL, &j0);
      ms1 = c1->symstate[q1];
      UnpackSymStateTR(mtr, ms1, ip, &s1, &m1, &is1);
      DecodePJ(s1, NULL, &j1);
      //is0 = 2*mbpt_tr.nktr*s0 + 2*ip;
      if (IsOdd(p0)) is0++;
      //is1 = IBisect(s1, mtr[is0].nsym1, mtr[is0].isym1);
      is1 = IdxGet(mtr[is0].ids, s1);
      if (is1 < 0) continue;
      pma = mtr[is0].pma[is1][m0*mtr[is0].ndim1[is1] + m1];
      if (pma == NULL) continue;
      q0 = k*mbpt_tr.nktr + ip;
      c = 0.0;
      for (kk = kmin; kk <= kmax; kk += 2) {
	kk2 = kk/2;
	c += a[kk2][q0]*yk[kk2];
      }
      if (IsOdd(abs(j0-j1+s[2].j-s[3].j+s[4].j-s[5].j+s[0].j-s[1].j)/2)) {
	c = -c;
      }
      for (p1 = 0; p1 < n; p1++) {
	q1 = p1*ng + i0;
	rt = c*r1[p1]/d2;
#pragma omp atomic
	pma->rma[q1] += rt;
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
  ORBITAL *orb0, *orb1, *orb2, *orb3;

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

  orb3 = GetOrbital(k3);
  d1 = orb3->energy + orb3->qed;
  orb2 = GetOrbital(k2);
  d1 -= orb2->energy + orb2->qed;
  orb0 = GetOrbital(k0);
  d2 = orb0->energy + orb0->qed;
  orb1 = GetOrbital(k1);
  d2 -= orb1->energy + orb1->qed;
  int *sp, *mp;
  double *sd1p, *sd2p, *c12p;
  int warned = 0, ignored = 0;
  sd1p = mbptjp.sd1a;
  sd2p = mbptjp.sd2a;
  c12p = mbptjp.c12a;
  sp = mbptjp.sa;
  mp = mbptjp.ma;
  for (k = 0; k < mst; k++) {
    q0 = bst[k];
    q1 = kst[k];
    sp[k] = -1;
    ms0 = c0->symstate[q0];
    UnpackSymStateMBPT(meff, ms0, &s0, &m0);
    ms1 = c1->symstate[q1];
    UnpackSymStateMBPT(meff, ms1, &s1, &m1);
    if (a[k] == 0) continue;
    if (m1 > m0) {
      m = m1*(m1+1)/2 + m0;
    } else {
      m = m0*(m0+1)/2 + m1;
    }
    sp[k] = s0;
    mp[k] = m;
    y = r1*r2*a[k];
    double y12 = y/(d1*d2);
    cd1 = y/d1;
    cd2 = y/d2;
    sd1p[k] = cd1;
    sd2p[k] = cd2;
    c12p[k] = y12;
    double sd1w = sqrt(fabs(cd1/d1));
    double sd2w = sqrt(fabs(cd2/d2));
    double sd1s = sd1w;//c1->cth>0?c1->cth/fabs(d1):0;
    double sd2s = sd2w;//c0->cth>0?c0->cth/fabs(d2):0;
    if (isnan(y12)) {
      sd1s = 1E50;
      sd2s = 1E50;
    } else if (isinf(y12)) {
      sd1s = 2E50;
      sd2s = 2E50;
    }
    if (meff[s0]->imbpt[m] > 1 && !warned) {
      double wth = (mbpt_nwmix-meff[s0]->imbpt[m]);
      wth /= mbpt_nwmix - 2.0;
      wth = mbpt_warn + mbpt_mwarn*wth;
      if (wth > 0 && (sd1w > wth || sd2w > wth)) {
	int na[2], ka[2];
	na[0] = orb0->n;
	na[1] = orb1->n;
	ka[0] = orb0->kappa;
	ka[1] = orb1->kappa;
	CONFIG *cr = ExciteNRConfig(c0, 2, na, ka);
	char scr[2048] = "NULL";
	if (mbpt_cca) {
	  cr->sweight = d1;
	  cr->energy = d2;
	  cr->delta = sd1w;
	  cr->sth = sd2w;
	  cr->cth = sd1s;
	  cr->mde = sd2s;
	  cr->icfg = meff[s0]->imbpt[m];
	  cr->n_electrons = c0->igroup;
	  cr->n_csfs = c0->icfg;
	  cr->nnrs = c1->igroup;
	  cr->igroup = c1->icfg;
	  ArrayAppend(mbpt_cca, &cr, NULL);
	  warned = 1;
	} else {
	  if (cr != NULL) {
	    ConstructNRConfigName(scr, 2048, cr);
	    FreeConfigData(cr);
	    free(cr);
	  }	
	  printf("large h11term warn: %s %d %d %d %d %d %d %d %g %g %g %g %g\n",
		 scr,
		 c0->igroup, c0->icfg, c1->igroup, c1->icfg, s0, m0, m1,
		 sd1s, sd2s, d1, d2, wth);
	  warned = 2;
	}
      }
    }
    if (mbpt_ignore > 0) {// && (c0 != c1 || m0 != m1)) {
      if (sd1s > mbpt_ignore || sd2s > mbpt_ignore) {
	if (warned != 2 && (mbpt_ignorep&1)) {
	  int na[2], ka[2];
	  na[0] = orb0->n;
	  na[1] = orb1->n;
	  ka[0] = orb0->kappa;
	  ka[1] = orb1->kappa;
	  CONFIG *cr = ExciteConfig(c0, 2, na, ka);
	  char scr[2048] = "NULL";
	  if (cr != NULL) {
	    ConstructConfigName(scr, 2048, cr);
	    FreeConfigData(cr);
	    free(cr);
	  }	
	  printf("large h11term ignore: %s %d %d %d %d %d %d %d %g %g %g %g %g %d\n",
		 scr,
		 c0->igroup, c0->icfg, c1->igroup, c1->icfg, s0, m0, m1,
		 sd1s, sd2s, d1, d2, mbpt_ignore, meff[s0]->imbpt[m]);
	}
	cd1 = 0;
	cd2 = 0;
	y12 = 0;
	ignored = 1;
#pragma omp atomic
	mbpt_ignoren++;
	break;
      }
    }
  }
  if (!ignored) {
    for (k = 0; k < mst; k++) {
      if (sp[k] < 0) continue;
      s0 = sp[k];
      m = mp[k];	
      h1 = meff[s0]->hab1[m];
      h2 = meff[s0]->hba1[m];
      ng = meff[s0]->n;
#pragma omp atomic
      h1[i0] += sd1p[k];
#pragma omp atomic
      h2[i0] += sd2p[k];
      i0g = i0 + ng;
#pragma omp atomic
      h1[i0g] += c12p[k];
#pragma omp atomic
      h2[i0g] += c12p[k];
    }
  }
}

void TR11Term(MBPT_TR *mtr, CONFIG *c0, CONFIG *c1, 
	      int ns, SHELL *bra, SHELL *ket,
	      SHELL_STATE *sbra, SHELL_STATE *sket, 
	      int mst, int *bst, int *kst,
	      INTERACT_SHELL *s, int ph, int k0, int k1, int k2, int k3,
	      FORMULA *fm, double *a, int i0, int ng) {
  double d2, *r1, r2, c, d;
  int q0, q1, k, m0, m1, m, ms0, ms1, s0, s1, p, ip, p0, n, is0, is1;
  int gauge, p1, j0, j1;
  ORBITAL *orb3, *orb2;
  MBPT_PMA *pma;
  
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
  if (!TriangleTR(s[0].j, s[1].j)) return;
  if (s[2].j != s[3].j) return;
  if (md <= 1) {
    FixJsZ(s, fm);
    for (ip = 0; ip < mbpt_tr.nktr; ip++) {
      p = mbpt_tr.ks[ip];
      fm->js[5] = 2*p;
      fm->js[6] = 0;
      fm->js[7] = 2*p;
      for (k = 0; k < mst; k++) {
	q0 = bst[k];
	q1 = kst[k];
	FixTotalJ(ns, sbra, bra, c0, q0);
	FixTotalJ(ns, sket, ket, c1, q1);	
	EvaluateTensor(ns, sbra, sket, s, 1, fm);
	q0 = k*mbpt_tr.nktr + ip;
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

  orb3 = GetOrbital(k3);
  d2 = orb3->energy + orb3->qed;
  orb2 = GetOrbital(k2);
  d2 -= orb2->energy + orb2->qed;
  double sd1w = sqrt(fabs(r2/d2));
  double sd1s = c1->cth>0?c1->cth/fabs(d2):0;
  int warned = 0;
  if (c0->icfg >= 0 && c1->icfg >= 0) {
    if (mbpt_warntr > 0 && sd1w > mbpt_warntr) {
      int na[2], ka[2];
      na[0] = orb2->n;
      na[1] = orb3->n;
      ka[0] = orb2->kappa;
      ka[1] = orb3->kappa;
      CONFIG *cr = ExciteNRConfig(c1, 2, na, ka);
      char scr[2048] = "NULL";
      if (cr != NULL) {
	ConstructNRConfigName(scr, 2048, cr);
	FreeConfigData(cr);
	free(cr);
      }
      printf("large t11term warn: %s %d %d %d %d %g %g %g\n",
	     scr, c0->igroup, c0->icfg, c1->igroup, c1->icfg,
	     sd1s, d2, mbpt_warntr);
      warned = 1;
    }
  }
  if (mbpt_ignoretr > 0 && sd1s > mbpt_ignoretr) {
    if (!warned && (mbpt_ignorep&2)) {
      int na[2], ka[2];
      na[0] = orb2->n;
      na[1] = orb3->n;
      ka[0] = orb2->kappa;
      ka[1] = orb3->kappa;
      CONFIG *cr = ExciteNRConfig(c1, 2, na, ka);
      char scr[2048] = "NULL";
      if (cr != NULL) {
	ConstructNRConfigName(scr, 2048, cr);
	FreeConfigData(cr);
	free(cr);
      }
      printf("large t11term ignore: %s %d %d %d %d %g %g %g\n",
	     scr, c0->igroup, c0->icfg, c1->igroup, c1->icfg,
		 sd1s, d2, mbpt_ignoretr);
    }
#pragma omp atomic
    mbpt_ignoren++;
    return;
  }
    
  /*
  d2 = GetOrbital(k3)->energy - GetOrbital(k2)->energy;
  */
  ng -= mbpt_mini;
  gauge = GetTransitionGauge();
  for (ip = 0; ip < mbpt_tr.nktr; ip++) {
    p = mbpt_tr.ks[ip];
    if (!Triangle(s[0].j, s[1].j, 2*p)) continue;
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
      UnpackSymStateTR(mtr, ms0, ip, &s0, &m0, &is0);
      ms1 = c1->symstate[q1];
      UnpackSymStateTR(mtr, ms1, ip, &s1, &m1, &is1);
      if (IsOdd(p0)) is0++;
      //is1 = IBisect(s1, mtr[is0].nsym1, mtr[is0].isym1);
      is1 = IdxGet(mtr[is0].ids, s1);
      if (is1 < 0) continue;
      int m01 = m0*mtr[is0].ndim1[is1] + m1;
      pma = mtr[is0].pma[is1][m01];
      if (pma == NULL) continue;
      q0 = k*mbpt_tr.nktr + ip;
      for (p1 = 0; p1 < n; p1++) {
	q1 = p1*ng + i0;
	d = a[q0]*r1[p1]*r2/d2;
#pragma omp atomic
	pma->tma[q1] += d;
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
      UnpackSymStateTR(mtr, ms0, ip, &s0, &m0, &is0);
      DecodePJ(s0, NULL, &j0);
      ms1 = c1->symstate[q1];
      UnpackSymStateTR(mtr, ms1, ip, &s1, &m1, &is1);
      DecodePJ(s1, NULL, &j1);
      if (IsOdd(p0)) is0++;
      //is1 = IBisect(s1, mtr[is0].nsym1, mtr[is0].isym1);
      is1 = IdxGet(mtr[is0].ids, s1);
      if (is1 < 0) continue;
      pma = mtr[is0].pma[is1][m0*mtr[is0].ndim1[is1] + m1];
      if (pma == NULL) continue;
      q0 = k*mbpt_tr.nktr + ip;
      if (IsOdd((abs(j0-j1+s[0].j-s[1].j)/2))) c = -r2;
      else c = r2;
      for (p1 = 0; p1 < n; p1++) {
	q1 = p1*ng + i0;	
	d = a[q0]*r1[p1]*c/d2;
#pragma omp atomic
	pma->rma[q1] += d;
      }
    }
  }
}

int IdxGetN1(IDXARY *ing, int m) {
  int i = IdxGet(ing, m);
  if (i < 0) return i;
  if (mbpt_mini <= 0) return i;
  i -= mbpt_mini;
  if (i < 0) i = 0;
  return i;
}

int IdxGetN2(IDXARY *ing, int m) {
  int i = IdxGet(ing, m);
  if (i < 0) return i;
  if (mbpt_mini2 <= 0) return i;
  i -= mbpt_mini2;
  if (i < 0) return -999;
  return i;
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
  i1 = IdxGetN1(ing, m1);
  if (i1 < 0) return;
  if (mbpt_maxn <= 0) {
    i2 = IdxGetN2(ing2, m2-m1);
    if (i2 == -999) i = -(i1+1);
    else {
      if (i2 < 0) return;
      i = i1*(ing2->n-mbpt_mini2) + i2;
    }
  } else {
    i = -(i1+1);
  }
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
	      if (SkipMPIM(1)) continue;
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
	      if (SkipMPIM(1)) continue;
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
		if (SkipMPIM(1)) continue;
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
		      i1 = IdxGetN1(ing, m1);
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
  i1 = IdxGetN1(ing, i1);
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
		  //int kket = k;
		  if (k >= 0) continue;
		  op[0] = ik;
		  op[1] = im;
		  om[0] = ia;
		  om[1] = ib;
		  k = CheckConfig(ns, bra, 2, op, 2, om, nc, cs);
		  //k = CheckConfig(ns, bra, 2, op, 2, om, 0, NULL);
		  if (k >= 0) {
		    continue;
		  }
		  //int kbra = k;
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
		  if (SkipMPIM(1)) continue;
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
      a[i] = malloc(sizeof(double)*mst*mbpt_tr.nktr);
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
	    if (SkipMPIM(1)) continue;
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
		    i1 = IdxGetN1(ing, m1);
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
  i1 = IdxGetN1(ing, i1);
  if (i1 < 0) return;
  k = MKK;
  for (i = 0; i < k; i++) {
    if (mode == 0) {
      a[i] = malloc(sizeof(double)*mst);
    } else {
      a[i] = malloc(sizeof(double)*mst*mbpt_tr.nktr);
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
	      if (SkipMPIM(1)) continue;
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
	      k = CheckConfig(ns, bra, 1, op, 1, om, nc, cs);
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
    a = malloc(sizeof(double)*mst*mbpt_tr.nktr);
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
	if (SkipMPIM(1)) continue;
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
	  i1 = IdxGetN1(ing, m1);
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
  i1 = IdxGetN1(ing, i1);
  if (i1 < 0) return;
  if (mode == 0) {
    a = malloc(sizeof(double)*mst);
  } else {
    a = malloc(sizeof(double)*mst*mbpt_tr.nktr);
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
	  if (SkipMPIM(1)) continue;
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
	  k = CheckConfig(ns, bra, 1, op, 1, om, nc, cs);
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
	  if (meff[i]->n2 > 0) {
	    free(meff[i]->hab[m]);
	    free(meff[i]->hba[m]);
	  }
	}
      }
      free(meff[i]->hab1);
      free(meff[i]->hba1);
      if (meff[i]->n2 > 0) {
	free(meff[i]->hab);
	free(meff[i]->hba);
      }
      free(meff[i]->h0);
      free(meff[i]->e0);
      free(meff[i]->heff);
      //free(meff[i]->basis);
      FreeIdxAry(meff[i]->idb, 0);
    }
    free(meff[i]);
  }
}

void InitTransitionMBPT(MBPT_TR **mtr0, IDXARY **idb) {
  MBPT_TR *mtr;
  HAMILTON *h, *h1;
  int i, ik, j, k, m, isym, k0, k1, m0, m1, ms0, mst;

  mtr = malloc(sizeof(MBPT_TR)*(2*mbpt_tr.nktr)*MAX_SYMMETRIES);
  *mtr0 = mtr;
  j = 0;
  for (isym = 0; isym < MAX_SYMMETRIES; isym++) {
    h = GetHamilton(isym);
    DecodePJ(isym, &k1, &k0);
    for (ik = 0; ik < mbpt_tr.nktr; ik++) {
      i = mbpt_tr.ks[ik];
      for (k = -1; k <= 1; k += 2) {
	mtr[j].isym0 = isym;
	mtr[j].m = k*i;
	mtr[j].nsym1 = 0;
	mtr[j].ndim0 = h->dim;
	mtr[j].idb = idb;
	if (h->dim > 0) {
	  mtr[j].nsym1 = 2*i + 1;
	  mtr[j].isym1 = malloc(sizeof(int)*mtr[j].nsym1);
	  mtr[j].ndim1 = malloc(sizeof(int)*mtr[j].nsym1);
	  mtr[j].pma = malloc(sizeof(MBPT_PMA **)*mtr[j].nsym1);
	  for (m = 0; m < mtr[j].nsym1; m++) {
	    mtr[j].ndim1[m] = 0;
	    mtr[j].pma[m] = NULL;
	  }
	  m = 0;
	  for (m0 = k0 - 2*i; m0 <= k0 + 2*i; m0 += 2) {
	    if (m0 < 0) {	      
	      continue;
	    }
	    m1 = k1;
	    if (k < 0 && IsOdd(i)) m1 = 1-m1;
	    if (k > 0 && IsEven(i)) m1 = 1-m1;
	    mtr[j].isym1[m] = IsEven(m1)? 2*m0 : (2*m0+1);
	    if (mtr[j].isym1[m] >= MAX_SYMMETRIES) continue;
	    h1 = GetHamilton(mtr[j].isym1[m]);
	    if (h1->dim > 0) {
	      mtr[j].ndim1[m] = h1->dim;
	      mst = h->dim * h1->dim;
	      mtr[j].pma[m] = malloc(sizeof(MBPT_PMA *)*mst);
	      for (ms0 = 0; ms0 < mst; ms0++) {
		mtr[j].pma[m][ms0] = NULL;
	      }
	      m++;
	    }
	  }
	  mtr[j].nsym1 = m;
	  if (mtr[j].nsym1 == 0) {
	    free(mtr[j].isym1);
	    free(mtr[j].ndim1);
	    free(mtr[j].pma);
	    mtr[j].ids = NULL;
	  } else {
	    mtr[j].ids = malloc(sizeof(IDXARY));
	    InitIdxAry(mtr[j].ids, mtr[j].nsym1, mtr[j].isym1);
	  }
	}
	j++;
      }
    }
  }
}

void FreeTransitionDataMBPT(MBPT_TR *mtr, int j) {
  int m, i, k;
  
  if (mtr[j].nsym1 == 0) return;
  for (m = 0; m < mtr[j].nsym1; m++) {
    k = mtr[j].ndim0 * mtr[j].ndim1[m];
    for (i = 0; i < k; i++) {
      if (mtr[j].pma[m][i] != NULL) {
	free(mtr[j].pma[m][i]->tma);
	free(mtr[j].pma[m][i]->rma);
      }
    }
    free(mtr[j].pma[m]);
  }
  FreeIdxAry(mtr[j].ids, 0);
  free(mtr[j].ndim1);
  free(mtr[j].pma);
  mtr[j].nsym1 = 0;
}

void FreeTransitionMBPT(MBPT_TR *mtr) {
  int j, n;

  if (mbpt_tr.nktr <= 0) return;
  n = 2*mbpt_tr.nktr*MAX_SYMMETRIES;
  for (j = 0; j < n; j++) {
    if (mtr[j].nsym1 == 0) continue;
    FreeTransitionDataMBPT(mtr, j);
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

int GetICP(int nt, CONFIG_PAIR *cp, int ncp, int icp, int *i0, int *i1) {
  int i, m, mb, mb1, mr, k;
  m = 0;
  for (i = 0; i < nt; i++) {
    m += cp[i].m;
  }
  if (ncp <= 0) {
    *i0 = 0;
    *i1 = nt;
    return m;
  }
  if (ncp >= nt) {
    *i0 = icp;
    *i1 = icp+1;
    return cp[icp].m;
  }
  mb = m/ncp;
  mr = m%ncp;
  printf("geticp: %d %d %d %d %d %d\n", icp, ncp, nt, m, mb, mr);
  k = 0;
  m = 0;
  *i0 = 0;
  mb1 = mb;
  if (k < mr) mb1++;
  for (i = 0; i < nt; i++) {
    m += cp[i].m;
    if (m >= mb1) {
      if (k == icp) {
	*i1 = i+1;
       return m;
      }
      k++;
      if (k >= mr) mb1 = mb;
      m = 0;
      *i0 = i+1;
    }      
  }
  *i1 = nt;
  return m;
}

int CompareMBPTCC(const void *p1, const void *p2) {
  CONFIG *c1, *c2;
  c1 = *((CONFIG **) p1);
  c2 = *((CONFIG **) p2);
  if (c1->n_shells != c2->n_shells) return c1->n_shells-c2->n_shells;
  int r = memcmp(c1->shells, c2->shells, sizeof(SHELL)*c1->n_shells);
  if (r) return r;
  double d2 = Max(c2->delta,c2->cth);
  double d1 = Max(c1->delta,c1->cth);
  if (d2 < d1) return -1;
  else if (d2 > d1) return 1;
  return 0;
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
int StructureMBPT1(char *fn, char *fn0, char *fn1,
		   int nkg, int *kg, int nk, int *nkm, 
		   int n, int *ng, int n2, int *ng2, int nkg00, int ncp,
		   int icpi, int icpf) {
  int *bas, *bas0, *bas1, nlevels, nb0, nb, n3, q, nhab, nhab1;
  int i, j, k, i0, i1, n0, n1, nr, nr2, *ngr, *ngr2;
  int ic, icp0, icp1, ncpt, icp, nkg0;
  int isym, ierr, nc, m, mks, *ks;
  int pp, jj, nmax, na, *ga, k0, k1, m0, m1, nmax1, mst, ncps;
  int p0, p1, j0, j1, j2, q0, q1, ms0, ms1, *bst, *kst, *bst0, *kst0;
  int pp0, pp1, jj0, jj1;
  char tfn[1024];
  SYMMETRY *sym;
  STATE *st;
  HAMILTON *h;
  SHELL *bra, *ket, *bra1, *ket1, *bra2, *ket2;
  SHELL_STATE *sbra, *sket, *sbra1, *sket1, *sbra2, *sket2;
  CONFIG **cs, *c0, *c1, *ct0, *ct1;
  CONFIG_GROUP *g;
  CONFIG_PAIR *cfgpair;
  ORBITAL *orb0, *orb1;
  MBPT_EFF *meff[MAX_SYMMETRIES];
  MBPT_TR *mtr;
  IDXARY *idb[MAX_SYMMETRIES];
  double a, b, c, *mix, *hab, *hba, emin, emax;
  double *h0, *heff, *neff, *hab1, *hba1, *dw, tt0, tt1, tbg, dt, dtt;
  FILE *f;
  IDXARY ing, ing2;
  ORBITAL *orb;
  int nkgp;
  int *kgp = NULL;
  int ip = 0, ncca = 0;
  
  if (mbpt_nwmix <= 0) mbpt_nwmix = 1.0/mbpt_wmix;
  if (nkg00 == 0) {
    printf("no valid configs for mbpt: %d %d\n", nkg, nkg00);
    return 0;
  }
  ing.n = ing.m = 0;
  ing2.n = ing2.m = 0;
  ierr = 0;
  n3 = mbpt_n3;
  int rh = 0;
  if (fn0 != NULL && nkg == 0) {
    ReadHamilton(fn0, &nkg0, &nkg, &kg, &nkgp, &kgp, -1);
    if (nkgp > 0) {
      i = nkg + nkgp;
      kg = realloc(kg, sizeof(int)*i);
      memcpy(kg+nkg, kgp, sizeof(int)*nkgp);
      nkg = i;
      free(kgp);
      kgp = NULL;
    }
    rh = 1;
  } else {
    nkg0 = abs(nkg00);
    if (nkg0 > nkg) nkg0 = nkg;
    nkgp = nkg-nkg0;    
  }
  printf("mbpt configs: %d %d %d\n", nkg, nkg0, nkgp);
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
#pragma omp master
	{
	  if (i >= nkg0) {
	    cs[k]->icfg = -(1+cs[k]->icfg);	  
	  }
	}
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
  }
  nmax = nmaxm;
  nmax1 = nmaxm1;
  cs = csm;
  mbpt_nmk0 = nmax*nmax;
  mbpt_nmk1 = nmax1*nmax1+1;
  int nmk2 = mbpt_nmk0*mbpt_nmk1;
  mbpt_csary = malloc(sizeof(ARRAY)*nmk2);
  for (i = 0; i < nmk2; i++) {
    ArrayInit(&mbpt_csary[i], sizeof(int), CONFIGS_BLOCK);
  }
  for (k = 0; k < nc; k++) {
    i = ShellToInt(cs[k]->shells[0].n, cs[k]->shells[0].kappa);
    if (cs[k]->shells[0].nq > 1) {
      j = i;
    } else if (cs[k]->n_shells > 1) {
      j = ShellToInt(cs[k]->shells[1].n, cs[k]->shells[1].kappa);
    } else {
      j = mbpt_nmk1-1;
    }
    m = i*mbpt_nmk1 + j;
    ArrayAppend(&mbpt_csary[m], &k, InitIntData);
  }      
    
  tt0 = WallTime();
  tbg = tt0;

  MPrintf(-1, "Construct Radial Basis, %d/%d...%11.4E %11.4E\n",
	  MyRankMPI(), NProcMPI(), TotalSize(), TotalArraySize());
  fflush(stdout);
  if (ng == NULL) {
    int n1min = n/1000;
    int n1max = n%1000;
    if (n1min < 1) n1min = 1;
    n = n1max-n1min+1;
    ng = malloc(sizeof(int)*n);
    for (i = 0; i < n; i++) {
      ng[i] = i+n1min;
    }
  } else {
    n = ConstructNGrid(n, &ng);
  }
  if (ng2 == NULL) {
    mbpt_maxn = n2%1000;
    mbpt_maxm = n2/1000;
    if (mbpt_maxm > 0) {
      n2 = mbpt_maxn - ng[0] + 1;
    } else {
      n2 = mbpt_maxn;
    }
    ng2 = malloc(sizeof(int)*n2);
    for (i = 0; i < n2; i++) {
      ng2[i] = i;
    }
  } else {
    n2 = ConstructNGrid(n2, &ng2);
    mbpt_maxn = 0;
    mbpt_maxm = 0;
    mbpt_minn2 = 0;
    mbpt_mini2 = 0;
    for (i = 1; i < n2; i++) {
      if (ng2[i]-ng2[i-1] > 1) {
	mbpt_minn2 = ng2[i-1];
	break;
      }
    }
    if (i == n2) {
      mbpt_maxn = ng2[i-1];
      mbpt_minn2 = 0;
    }
  }
  InitIdxAry(&ing, n, ng);
  InitIdxAry(&ing2, n2, ng2);
  if (mbpt_minn > 0) {
    mbpt_mini = IdxGet(&ing, mbpt_minn);
    if (mbpt_mini == -1) {
      mbpt_minn = 0;
      mbpt_mini = 0;
    } else if (mbpt_mini == -2) {      
      mbpt_mini = n-1;
    }
    nr = n - mbpt_mini;
    ngr = malloc(sizeof(int)*nr);
    for (i = 0; i < nr; i++) {
      ngr[i] = ng[i+mbpt_mini];
      if (ngr[i] < mbpt_minn) ngr[i] = mbpt_minn;
    }
  } else {
    nr = n;
    ngr = ng;
  }
  if (mbpt_minn2 > 0) {
    mbpt_mini2 = IdxGet(&ing2, mbpt_minn2);
    if (mbpt_mini2 < 0) {
      mbpt_minn2 = 0;
      mbpt_mini2 = 0;
    }
    nr2 = n2 - mbpt_mini2;
    ngr2 = malloc(sizeof(int)*nr2);
    for (i = 0; i < nr2; i++) {
      ngr2[i] = ng2[i+mbpt_mini2];
    }
  } else {
    nr2 = n2;
    ngr2 = ng2;
  }
  PrintMBPTOptions();
  nhab1 = nr*2;
  if (mbpt_maxn > 0) {
    nhab = 0;
  } else {
    nhab = nr*nr2*2;
  }
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

  tt1 = WallTime();
  dt = tt1-tt0;
  tt0 = tt1;
  MPrintf(-1, "Time = %12.5E %d %d %d %d %d %d %d %d %d %d\n",
	  dt, n, nr, n2, nr2, mbpt_maxm, mbpt_maxn,
	  mbpt_minn, mbpt_mini, mbpt_minn2, mbpt_mini2);
  fflush(stdout);
  if (nb < 0) return -1;

  MPrintf(-1, "CI Structure...%10.4E %10.4E %10.4E\n",
	  WallTime()-tbg, TotalSize(), TotalArraySize());
  fflush(stdout);
  nlevels = GetNumLevels();
  emax = -1E31;
  emin = 1E31;
		 
  if (nkgp > 0) {
    kgp = kg + nkg0;
    if (nkg00 < 0) {
      ip = 0;
    } else {
      ip = 1;
    }
  }
  char *fne = NULL;
  if (n3 < 0 && mbpt_tr.nktr <= 0) {
    fne = fn;
  }
  if (rh) {
    SolveStructure(fne, fn0, 0, NULL, 0, NULL, 0);
  } else {
    int *kg0 = kg;
    if (fne) {
      if (nkgp > 0) {
	kgp = malloc(sizeof(int)*nkgp);
	memcpy(kgp, kg+nkg0, sizeof(int)*nkgp);
      }
      kg0 = malloc(sizeof(int)*nkg0);
      memcpy(kg0, kg, sizeof(int)*nkg0);
    }
    SolveStructure(fne, fn0, nkg0, kg0, nkgp, kgp, ip);
  }
  if (fne) {
    for (isym = 0; isym < MAX_SYMMETRIES; isym++) {
      meff[isym] = NULL;
    }
    goto ERROR;
  }

  int npr, mr;
  mr = MPIRank(&npr);
  double *em0 = WorkSpace();
  double *em1 = em0 + npr;
  for (i = 0; i < npr; i++) {
    em0[i] = 1e31;
    em1[i] = -1e31;
  }
  ResetWidMPI();
#pragma omp parallel default(shared) private(isym, h, sym, ks, mks, m, mix, k, q, st, i, j, a, b, c, i0, mr)
  {
  mr = MPIRank(NULL);
  for (isym = 0; isym < MAX_SYMMETRIES; isym++) {    
    int skip = SkipMPI();
    if (skip) continue;
    meff[isym] = NULL;
    idb[isym] = NULL;
    h = GetHamilton(isym);
    if (h->dim <= 0) {
      continue;
    }
    sym = GetSymmetry(isym);
    meff[isym] = (MBPT_EFF *) malloc(sizeof(MBPT_EFF));
    meff[isym]->h0 = malloc(sizeof(double)*h->dsize);
    meff[isym]->e0 = malloc(sizeof(double)*h->dim);
    memcpy(meff[isym]->h0, h->hamilton, sizeof(double)*h->dsize);
    meff[isym]->basis = malloc(sizeof(int)*h->dim);
    meff[isym]->nbasis = h->dim;
    int dim2 = h->dim*h->dim;
    meff[isym]->hsize0 = dim2;
    meff[isym]->hsize = h->dsize;
    memcpy(meff[isym]->basis, h->basis, sizeof(int)*h->dim);
    meff[isym]->idb = malloc(sizeof(IDXARY));
    InitIdxAry(meff[isym]->idb, meff[isym]->nbasis, meff[isym]->basis);
    idb[isym] = meff[isym]->idb;
    meff[isym]->heff = malloc(sizeof(double)*dim2);
    meff[isym]->neff = malloc(sizeof(double)*dim2);
    if (h->n_basis > h->dim && fn0 == NULL) {
      meff[isym]->heff0 = malloc(sizeof(double)*dim2);
      double *ap0 = h->work+h->lwork+3*h->dsize2;
      memcpy(meff[isym]->heff0, ap0, sizeof(double)*dim2);
    } else {
      meff[isym]->heff0 = NULL;
    }
    meff[isym]->imbpt = malloc(sizeof(int)*h->dsize);
    meff[isym]->hab1 = malloc(sizeof(double *)*h->dsize);
    meff[isym]->hba1 = malloc(sizeof(double *)*h->dsize);
    meff[isym]->n = nr;
    if (nhab > 0) {
      meff[isym]->hab = malloc(sizeof(double *)*h->dsize);
      meff[isym]->hba = malloc(sizeof(double *)*h->dsize);
      meff[isym]->n2 = nr2;
    } else {
      meff[isym]->hab = NULL;
      meff[isym]->hba = NULL;
      meff[isym]->n2 = 0;
    }    

    ks = h->iwork;
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
	  if (h->mixing[k] < em0[mr]) em0[mr] = h->mixing[k];
	  if (h->mixing[k] > em1[mr]) em1[mr] = h->mixing[k];
	}
      } else {
	if (mbpt_nlev <= 0 || IBisect(m, mbpt_nlev, mbpt_ilev) >= 0) {
	  ks[mks++] = k;
	}
	m++;
	if (h->mixing[k] < em0[mr]) em0[mr] = h->mixing[k];
	if (h->mixing[k] > em1[mr]) em1[mr] = h->mixing[k];
      }
      mix += h->n_basis;
    }
    int ki=0, ke=0;
    int jig, iig;
    double *mix0 = h->mixing+h->dim;
    for (i = 0; i < h->dim; i++) {
      h->work[i] = 0;
      for (m = 0; m < mks; m++) {
	k = ks[m];
	mix = mix0 + h->n_basis*k;
	double am = fabs(mix[i]);
	if (h->work[i] < am) h->work[i] = am;
      }
    }
    for (j = 0; j < h->dim; j++) {
      if (nkg0 != nkg) {
	st = (STATE *) ArrayGet(&(sym->states), h->basis[j]);
	jig = InGroups(st->kgroup, nkg0, kg);
      } else {
	jig = 1;
      }
      for (i = 0; i <= j; i++) {
	if (nkg0 != nkg) {
	  st = (STATE *) ArrayGet(&(sym->states), h->basis[i]);
	  iig = InGroups(st->kgroup, nkg0, kg);
	} else {
	  iig = 1;
	}
	c = h->work[i]*h->work[j];
	k = j*(j+1)/2 + i;
	a = mbpt_mcut;
	if (iig || jig) a *= mbpt_mcut2;
	if (iig && jig) a *= mbpt_mcut4;
	if (i == j) a *= mbpt_mcut3;
	if (c >= a) {
	  if (c >= mbpt_wmix) {
	    meff[isym]->imbpt[k] = (int)(1+c/mbpt_wmix);
	    if (meff[isym]->imbpt[k] > mbpt_nwmix) {
	      meff[isym]->imbpt[k] = (int)mbpt_nwmix;
	    }
	  } else {
	    meff[isym]->imbpt[k] = 1;
	  }
	  meff[isym]->hab1[k] = malloc(sizeof(double)*nhab1);
	  meff[isym]->hba1[k] = malloc(sizeof(double)*nhab1);
	  if (nhab > 0) {
	    meff[isym]->hab[k] = malloc(sizeof(double)*nhab);
	    meff[isym]->hba[k] = malloc(sizeof(double)*nhab);
	  }
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
	  meff[isym]->imbpt[k] = 0;
	  ke++;
	}
      }
    }
    MPrintf(-1, "Eff Ham: %3d %3d %3d %3d %3d %10.3E %10.3E %12.5E %12.5E\n",
	    isym, h->dim, h->n_basis, ki, ke,
	    WallTime()-tbg, TotalSize(), em0[mr], em1[mr]);
  }
  }
  for (i = 0; i < npr; i++) {
    if (em0[i] < emin) emin = em0[i];
    if (em1[i] > emax) emax = em1[i];
  }
  if (mbpt_tr.nktr > 0) {
    double mem0 = TotalSize();
    emax = (emax-emin);
    emin = 0.1;
    emax *= FINE_STRUCTURE_CONST;
    emin *= FINE_STRUCTURE_CONST;
    SetAWGridMBPT(emin, emax);
    InitTransitionMBPT(&mtr, idb);
    ResetWidMPI();
#pragma omp parallel default(shared) private(isym, sym, h, i0, i1, mix, i, q, st, k, m)
    {
    for (isym = 0; isym < MAX_SYMMETRIES; isym++) {
      int skip = SkipMPI();
      if (skip) continue;
      h = GetHamilton(isym);
      if (h->dim <= 0) continue;
      sym = GetSymmetry(isym);
      int *ilev0 = h->iwork;
      int *ilev1 = h->iwork+h->dim;
      i0 = 0;
      i1 = 0;
      mix = h->mixing + h->dim;
      for (i = 0; i < h->dim; i++) {
	q = GetPrincipleBasis(mix, h->dim, NULL);
	st = (STATE *) ArrayGet(&(sym->states), h->basis[q]);
	if (mbpt_tr.nlow == 0) {
	  k = InGroups(st->kgroup, nkg0, kg);
	} else {
	  k = InGroups(st->kgroup, mbpt_tr.nlow, mbpt_tr.low);
	}
	if (k) {
	  ilev0[i0] = i;
	  i0++;
	}
	if (mbpt_tr.nup == 0) {
	  k = InGroups(st->kgroup, nkg0, kg);
	} else {
	  k = InGroups(st->kgroup, mbpt_tr.nup, mbpt_tr.up);
	}
	if (k) {
	  ilev1[i1] = i;
	  i1++;
	}
	mix += h->n_basis;
      }      
      double *mix0 = h->mixing + h->dim;
      double *wm0 = h->work;
      double *wm1 = h->work + h->dim;
      for (i = 0; i < h->dim; i++) {
	wm0[i] = 0;
	for (m = 0; m < i0; m++) {
	  k = ilev0[m];
	  mix = mix0 + h->n_basis*k;
	  double am = fabs(mix[i]);
	  if (wm0[i] < am) wm0[i] = am;
	}
	wm1[i] = 0;
	for (m = 0; m < i1; m++) {
	  k = ilev1[m];
	  mix = mix0 + h->n_basis*k;
	  double am = fabs(mix[i]);
	  if (wm1[i] < am) wm1[i] = am;
	}
      }
      MPrintf(-1, "TR Matrix: %3d %5d %5d %5d %10.3E %10.3E\n",
	      isym, h->dim, i0, i1, WallTime()-tbg, TotalSize());
    }
    }
    ResetWidMPI();
#pragma omp parallel default(shared) private(a, b, c, i0, i1, q0, q1, p0, p1, pp0, pp1, jj0, jj1, j1, m, i)
    {
    HAMILTON *h0, *h1;
    a = mbpt_mcut*mbpt_mcut2;
    for (i0 = 0; i0 < MAX_SYMMETRIES; i0++) {
      int skip = SkipMPI();
      if (skip) continue;
      int ki = 0, ke = 0;
      h0 = GetHamilton(i0);
      if (h0->dim <= 0) continue;
      q0 = i0*2*mbpt_tr.nktr;
      q1 = q0 + 2*mbpt_tr.nktr;
      DecodePJ(i0, &pp0, &jj0);
      for (i1 = 0; i1 < MAX_SYMMETRIES; i1++) {
	h1 = GetHamilton(i1);
	if (h1->dim <= 0) continue;
	DecodePJ(i1, &pp1, &jj1);
	for (p0 = 0; p0 < h0->dim; p0++) {
	  for (p1 = 0; p1 < h1->dim; p1++) {	    
	    c = h0->work[p0]*h1->work[p1+h1->dim];
	    b = h0->work[p0+h0->dim]*h1->work[p1];
	    if (b > c) c = b;
	    if (c >= a) {
	      for (i = q0; i < q1; i++) {
		j1 = IdxGet(mtr[i].ids, i1);
		if (j1 < 0) continue;
		if (!Triangle(jj0, jj1, 2*abs(mtr[i].m))) continue;
		MBPT_PMA *pma = mtr[i].pma[j1][p0*mtr[i].ndim1[j1]+p1];
		if (pma == NULL) {
		  pma = malloc(sizeof(MBPT_PMA));
		  pma->tma = malloc(sizeof(double)*nr*mbpt_tr.naw);
		  pma->rma = malloc(sizeof(double)*nr*mbpt_tr.naw);
		  for (m = 0; m < nr*mbpt_tr.naw; m++) {
		    pma->tma[m] = 0;
		    pma->rma[m] = 0;
		  }
		  mtr[i].pma[j1][p0*mtr[i].ndim1[j1]+p1] = pma;
		  ki++;
		  //printf("pma0: %d %d %d %d %d %d %d %d %d %d %g %g\n", i0, i1, h0[i0]->pj, h1[i1]->pj, ilev0[i0], ilev1[i1], p0, p1, h1[i1]->dim, mtr[i].ndim1[j1], c, a);
		}		
	      }
	    } else {
	      ke++;
	    }
	  }
	}
      }
      double mem1 = TotalSize();
      MPrintf(-1, "TR Mem: %3d %5d %10.3E %10.3E %10.3E %10.3E %10.3E %d %d\n",
	      i0, h0->dim, WallTime()-tbg, mem0, mem1, mem1-mem0, a, ki, ke);
    }
    }
  }
  
  if (fn0 == NULL) {
    for (isym = 0; isym < MAX_SYMMETRIES; isym++) {
      h = GetHamilton(isym);
      if (h->dim <= 0) continue;
      AllocHamMem(h, -1, -1);
      AllocHamMem(h, 0, 0);
    }
  }

  if (n3 >= 0) {
    MPrintf(-1, "Construct Effective Hamiltonian %d %12.5E %12.5E\n",
	    nc, WallTime()-tbg, TotalSize());
    fflush(stdout);
    for (k0 = 0; k0 < nc; k0++) {
      c0 = cs[k0];
      a = ZerothEnergyConfig(c0);
      b = ZerothResidualConfig(c0);
      for (m0 = 0; m0 < c0->n_csfs; m0++) {
	ms0 = c0->symstate[m0];
	UnpackSymStateMBPT(meff, ms0, &i0, &q0);
	if (q0 < 0) continue;
	if (meff[i0] && meff[i0]->nbasis > 0) {
	  meff[i0]->e0[q0] = a + 0.5*b;
	}
      }
    }
    ncpt = 0;
    for (k0 = 0; k0 < nc; k0++) {
      c0 = cs[k0];      
      p0 = ConfigParity(c0);
      for (k1 = k0; k1 < nc; k1++) {
	c1 = cs[k1];
	p1 = ConfigParity(c1);
	if (p0 != p1) continue;
	m = 0;      
	for (m0 = 0; m0 < c0->n_csfs; m0++) {
	  ms0 = c0->symstate[m0];
	  UnpackSymStateMBPT(meff, ms0, &i0, &q0);
	  if (q0 < 0) continue;
	  if (c0 == c1) q = m0;
	  else q = 0;
	  for (m1 = q; m1 < c1->n_csfs; m1++) {
	    ms1 = c1->symstate[m1];
	    UnpackSymStateMBPT(meff, ms1, &i1, &q1);
	    if (q1 < 0) continue;
	    if (i0 != i1) continue;
	    if (q0 <= q1) {
	      k = q1*(q1+1)/2 + q0;
	    } else {
	      k = q0*(q0+1)/2 + q1;
	    }
	    if (meff[i0] && meff[i0]->nbasis > 0 && meff[i0]->imbpt[k]) {
	      m++;
	      break;
	    }
	  }
	}
	if (m == 0) continue;	
	ncpt++;
      }
    }    
    if (ncp > ncpt) {
      printf("structure ncp > ncpt: %d %d\n", ncp, ncpt);
      ncp = ncpt;
    }
    if (ncp > 1) {
      if (icpf == -1) icpf = icpi;
      else if (icpf == -2) icpf = ncp-1;
      else if (icpf < icpi) icpf = icpi;
      else if (icpf >= ncp) icpf = ncp-1;
      if (icpi < 0 || icpi >= ncp) {
	printf("invalid config splits for structure: %d %d %d\n",
	       icpi, icpf, ncp);
	Abort(1);
      }
    } else {
      icpi = 0;
      icpf = 0;
    }
    cfgpair = malloc(sizeof(CONFIG_PAIR)*ncpt);
    RANDIDX *rid = NULL;
    if (mbpt_rand) rid = RandList(ncpt);
    ic = 0;
    mst = 0;
    for (k0 = 0; k0 < nc; k0++) {
      c0 = cs[k0];      
      p0 = ConfigParity(c0);
      for (k1 = k0; k1 < nc; k1++) {
	c1 = cs[k1];
	p1 = ConfigParity(c1);
	if (p0 != p1) continue;
	m = 0;      
	for (m0 = 0; m0 < c0->n_csfs; m0++) {
	  ms0 = c0->symstate[m0];
	  UnpackSymStateMBPT(meff, ms0, &i0, &q0);
	  if (q0 < 0) continue;
	  if (c0 == c1) q = m0;
	  else q = 0;
	  for (m1 = q; m1 < c1->n_csfs; m1++) {
	    ms1 = c1->symstate[m1];
	    UnpackSymStateMBPT(meff, ms1, &i1, &q1);
	    if (q1 < 0) continue;
	    if (i0 != i1) continue;
	    if (q0 <= q1) {
	      k = q1*(q1+1)/2 + q0;
	    } else {
	      k = q0*(q0+1)/2 + q1;
	    }
	    if (meff[i0] && meff[i0]->nbasis > 0 && meff[i0]->imbpt[k]) {
	      m++;
	    }
	  }
	}
	if (m == 0) continue;
	mst += m;
	int icr = ic;
	if (rid) icr = rid[ic].i;
	cfgpair[icr].k0 = k0;
	cfgpair[icr].k1 = k1;
	cfgpair[icr].m = m;
	ic++;
      }
    }
    if (rid) free(rid);
    for (icp = icpi; icp <= icpf; icp++) {
      if (MyRankMPI() == 0) {
	char *pc = fn1;
	while(*pc && *pc != '%') pc++;
	if (*pc) {
	  sprintf(tfn, fn1, icp);
	} else if (icpf > icpi) {
	  sprintf(tfn, "%s%02d", fn1, icp);
	} else {
	  sprintf(tfn, "%s", fn1);
	}
	f = fopen(tfn, "w");
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
	  if (ngr != ng) free(ngr);
	  if (ngr2 != ng2) free(ngr2);
	  free(ng);
	  free(ng2);
	  for (i = 0; i < nmk2; i++) {
	    ArrayFree(&mbpt_csary[i], NULL);
	  }
	  free(mbpt_csary);
	  return -1;
	}
	if (mbpt_savesum == 0) {
	  fwrite(&nr, sizeof(int), 1, f);
	  fwrite(ngr, sizeof(int), nr, f);
	  if (nhab > 0) {
	    fwrite(&nr2, sizeof(int), 1, f);
	    fwrite(ngr2, sizeof(int), nr2, f);
	  } else {
	    fwrite(&nhab, sizeof(int), 1, f);
	  }
	} else {
	  m = 1;
	  fwrite(&m, sizeof(int), 1, f);
	  fwrite(&m, sizeof(int), 1, f);
	  m = 0;
	  fwrite(&m, sizeof(int), 1, f);
	}
	fwrite(&n3, sizeof(int), 1, f);
      }
      for (isym = 0; isym < MAX_SYMMETRIES; isym++) {
	h = GetHamilton(isym);
	if (h->dim <= 0) continue;
	if (meff[isym] == NULL) continue;
	if (meff[isym]->nbasis <= 0) continue;
	for (j = 0; j < h->dim; j++) {
	  for (i = 0; i <= j; i++) {
	    k = j*(j+1)/2 + i;
	    if (meff[isym]->imbpt[k] < 0) {
	      meff[isym]->imbpt[k] = -meff[isym]->imbpt[k];
	    }
	    if (meff[isym]->imbpt[k]) {
	      for (i0 = 0; i0 < nhab1; i0++) {
		meff[isym]->hab1[k][i0] = 0.0;
		meff[isym]->hba1[k][i0] = 0.0;
	      }
	      for (i0 = 0; i0 < nhab; i0++) {
		meff[isym]->hab[k][i0] = 0.0;
		meff[isym]->hba[k][i0] = 0.0;
	      }
	    }
	  }
	}
      }
      int tmst = GetICP(ncpt, cfgpair, ncp, icp, &icp0, &icp1);
      MPrintf(-1, "MBPT structure cfgpair: %d %d %d %d %d %d %d %d %d %12.5E %12.5E\n",
	      icp, icpi, icpf, ncp, icp0, icp1, ncpt, tmst, mst,
	      WallTime()-tbg, TotalSize());
      for (ic = 0; ic < ncpt; ic++) {
	if (ic >= icp0 && ic < icp1) continue;
	k0 = cfgpair[ic].k0;
	k1 = cfgpair[ic].k1;
	c0 = cs[k0];
	c1 = cs[k1];
	m = 0;
	for (m0 = 0; m0 < c0->n_csfs; m0++) {
	  ms0 = c0->symstate[m0];
	  UnpackSymStateMBPT(meff, ms0, &i0, &q0);
	  if (q0 < 0) continue;
	  if (c0 == c1) q = m0;
	  else q = 0;
	  for (m1 = q; m1 < c1->n_csfs; m1++) {
	    ms1 = c1->symstate[m1];
	    UnpackSymStateMBPT(meff, ms1, &i1, &q1);
	    if (q1 < 0) continue;
	    if (i0 != i1) continue;
	    if (q0 <= q1) {
	      k = q1*(q1+1)/2 + q0;
	    } else {
	      k = q0*(q0+1)/2 + q1;
	    }
	    if (meff[i0] && meff[i0]->nbasis > 0 && meff[i0]->imbpt[k] > 0) {
	      meff[i0]->imbpt[k] = -meff[i0]->imbpt[k];
	    }
	  }
	}
      }    
      double ttskip = 0, ttlock=0;
      long long tnlock = 0;
      mbpt_ignoren = 0;
      if (mbpt_omp0 == 1 || mbpt_omp0 == 3) {
	mbpt_omp = 1;
      } else {
	mbpt_omp = 0;
      }
      double *dm = NULL;
      int *im = NULL;      
      if (mbpt_msort) {
	dm = malloc(sizeof(double)*(icp1-icp0));
	im = malloc(sizeof(int)*(icp1-icp0));
	double amst = tmst/npr;
	for (ic = icp0; ic < icp1; ic++) {
	  dm[ic-icp0] = amst/cfgpair[ic].m;
	  if (cfgpair[ic].k0 != cfgpair[ic].k1) {
	    dm[ic-icp0] *= mbpt_asort;
	  }
	}
	ArgSort(icp1-icp0, dm, im);
      }
      MPrintf(-1, "MBPT structure beg: %12.5E %12.5E\n",
	      WallTime()-tbg, TotalSize());
      ResetWidMPI();
#pragma omp parallel default(shared) private(isym,n0,bra,ket,sbra,sket,bra1,ket1,bra2,ket2,sbra1,sket1,sbra2,sket2,cs,dt,dtt,k0,k1,c0,p0,c1,p1,m,bst0,kst0,m0,m1,ms0,ms1,q,q0,q1,k,mst,i0,i1,ct0,ct1,bst,kst,n1,bas0,bas1,ic,ncps)
      {
      int cmst = 0;
      ncps = 0;
      double ptt0=0, ptt1=0, tskip=0, tlock=0;
      long long nlock=0;      
      cs = mbpt_cs;
      bas0 = mbpt_bas0;
      bas1 = mbpt_bas1;
      mbpt_ibas0.n = mbpt_ibas0.m = 0;
      mbpt_ibas1.n = mbpt_ibas1.m = 0;
      if (mbpt_ccn[0]) {
	mbpt_cca = malloc(sizeof(ARRAY));
	ArrayInit(mbpt_cca, sizeof(CONFIG *), 4096);
      }
      int ic0;
      for (ic0 = icp0; ic0 < icp1; ic0++) {
	if (SkipMPIM(0)) continue;
	if (im) {
	  ic = icp0 + im[ic0-icp0];
	} else {
	  ic = ic0;
	}
	k0 = cfgpair[ic].k0;
	k1 = cfgpair[ic].k1;
	c0 = cs[k0];
	c1 = cs[k1];
	bst0 = malloc(sizeof(int)*cfgpair[ic].m);
	kst0 = malloc(sizeof(int)*cfgpair[ic].m);
	m = 0;
	for (m0 = 0; m0 < c0->n_csfs; m0++) {
	  ms0 = c0->symstate[m0];
	  UnpackSymStateMBPT(meff, ms0, &i0, &q0);
	  if (q0 < 0) continue;
	  if (c0 == c1) q = m0;
	  else q = 0;
	  for (m1 = q; m1 < c1->n_csfs; m1++) {
	    ms1 = c1->symstate[m1];
	    UnpackSymStateMBPT(meff, ms1, &i1, &q1);
	    if (q1 < 0) continue;
	    if (i0 != i1) continue;
	    if (q0 <= q1) {
	      k = q1*(q1+1)/2 + q0;
	    } else {
	      k = q0*(q0+1)/2 + q1;
	    }
	    if (meff[i0] && meff[i0]->nbasis > 0 && meff[i0]->imbpt[k] > 0) {
	      bst0[m] = m0;
	      kst0[m] = m1;
	      m++;
	    }
	  }
	}
	ptt0 = tt0;
	/* mst pairs */
	mst = m;
	cmst += m;
	/* if q0 <= q1 for the 1st pair, so are for the rest pairs */
	ms0 = c0->symstate[bst0[0]];
	ms1 = c1->symstate[kst0[0]];
	UnpackSymStateMBPT(meff, ms0, &i0, &q0);
	UnpackSymStateMBPT(meff, ms1, &i1, &q1);
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
	mbptjp.sa = malloc(sizeof(int)*mst);
	mbptjp.ma = malloc(sizeof(int)*mst);
	mbptjp.sd1a = malloc(sizeof(double)*mst);
	mbptjp.sd2a = malloc(sizeof(double)*mst);
	mbptjp.c12a = malloc(sizeof(double)*mst);
	if (n3 != 2) {
	  /* 1-b 2-b term no virtual orb */
	  DeltaH12M0(meff, n0, bra2, ket2, sbra2, sket2, mst, bst, kst,
		     ct0, ct1, &mbpt_ibas0, mbpt_bas0s, mbpt_bas0d,
		     &ing, nc, cs, 0);	      
	  /* 1-b 2-b term 1 virtual orb */
	  DeltaH12M1(meff, n0+1, bra1, ket1, sbra1, sket1, mst, bst, kst,
		     ct0, ct1, &mbpt_ibas0, mbpt_bas0s, mbpt_bas0d,
		     &mbpt_ibas1, &ing, nc, cs, 0);	  
	  /* 2-b 1-b term no virtual orb */
	  DeltaH12M0(meff, n0, ket2, bra2, sket2, sbra2, mst, kst, bst,
		     ct1, ct0, &mbpt_ibas0, mbpt_bas0s, mbpt_bas0d,
		     &ing, nc, cs, 0);
	  /* 2-b 1-b term 1 virtual orb */
	  DeltaH12M1(meff, n0+1, ket1, bra1, sket1, sbra1, mst, kst, bst,
		     ct1, ct0, &mbpt_ibas0, mbpt_bas0s, mbpt_bas0d,
		     &mbpt_ibas1, &ing, nc, cs, 0);	  
	  /* 1-b 1-b term no virtual orb */
	  DeltaH11M0(meff, n0, bra2, ket2, sbra2, sket2, mst, bst, kst,
		     ct0, ct1, &mbpt_ibas0, mbpt_bas0s, &ing, nc, cs, 0);  
	  /* 1-b 1-b term 1 virtual orb */	    
	  DeltaH11M1(meff, n0+1, bra1, ket1, sbra1, sket1, mst, bst, kst, 
		     ct0, ct1, &mbpt_ibas0, mbpt_bas0s,
		     &mbpt_ibas1, &ing, nc, cs, 0);
	  /* 2-b 2-b term no virtual */
	  DeltaH22M0(meff, n0, bra2, ket2, sbra2, sket2, mst, bst, kst,
		     ct0, ct1, &mbpt_ibas0,
		     mbpt_bas0s, mbpt_bas0d, &ing, nc, cs);
	  /* 2-b 2-b term 1 virtual */
	  DeltaH22M1(meff, n0+1, bra1, ket1, sbra1, sket1, mst, bst, kst,
		     ct0, ct1, &mbpt_ibas0, mbpt_bas0s, mbpt_bas0d,
		     &mbpt_ibas1, &ing, nc, cs);
	}
      	
	if (n3 != 1) {
	  /* 2-b 2-b term 2 virtual */
	  DeltaH22M2(meff, n0+2, bra, ket, sbra, sket, mst, bst, kst,
		     ct0, ct1, &mbpt_ibas0, mbpt_bas0s, mbpt_bas0d,
		     &mbpt_ibas1, &ing, &ing2, nc, cs);
	}
	double tmemf = TotalSize();
	double amemf = TotalArraySize();
	free(bra);
	free(ket);
	free(sbra);
	free(sket);
	free(bst0);
	free(kst0);
	free(mbptjp.sa);
	free(mbptjp.ma);
	free(mbptjp.sd1a);
	free(mbptjp.sd2a);
	free(mbptjp.c12a);
	ptt1 = WallTime();
	dt = ptt1-ptt0;
	dtt = ptt1-tbg;
	tt0 = ptt1;
	tskip = TimeSkip();
	tlock = TimeLock();
	nlock = NumLock();
	double tmem = TotalSize();
	double amem = TotalArraySize();
	MPrintf(mbpt_omp?0:-1, "%4d %4d %4d %4d %3d %3d %3d %3d %3d %3d ... %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %ld %ld %d\n", 
		ic0, ic, icp0, icp1, k0, k1, nc, mst, cmst, tmst,
		dt, dtt, tmem, amem, tmemf-tmem, amemf-amem,
		tskip, tlock, nlock, mbpt_ignoren, mbpt_cca?mbpt_cca->dim:0);
	fflush(stdout);	  
#pragma omp master
	{
	  ncps++;
	  if ((mbpt_reinit_ncps > 0 && ncps >= mbpt_reinit_ncps) ||
	      (mbpt_reinit_mem > 0 && tmem >= mbpt_reinit_mem)) {
	    SetRadialCleanFlags();
	    ncps = 0;
	  }
	}
      }
      ttskip = tskip;
      ttlock = tlock;
      tnlock = nlock;
      if (mbpt_ccn[0]) {
#pragma omp atomic
	ncca += mbpt_cca->dim;
      }
      }
      if (im) {
	free(dm);
	free(im);
      }
      MPrintf(-1, "MBPT Structure ... %12.5E %12.5E %12.5E %ld\n",
	      WallTime()-tbg, TotalSize(), TotalArraySize(), mbpt_ignoren);
      fflush(stdout);
      for (isym = 0; isym < MAX_SYMMETRIES; isym++) {      
	if (meff[isym] == NULL || meff[isym]->nbasis == 0) {
	  k = 0;
	  if (MyRankMPI() == 0) {
	    fwrite(&isym, sizeof(int), 1, f);
	    fwrite(&k, sizeof(int), 1, f);
	  }
	  continue;
	}
	double wt0 = WallTime();
	heff = meff[isym]->heff;
	neff = meff[isym]->neff;
#if USE_MPI == 1
	if (NProcMPI() > 1) {
	  for (i = 0; i < meff[isym]->hsize; i++) {
	    if (meff[isym]->hab1[i] != NULL) {
	      if (nhab > 0) {
		MPI_Allreduce(MPI_IN_PLACE, meff[isym]->hab[i],
			      nhab, MPI_DOUBLE,
			      MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, meff[isym]->hba[i],
			      nhab, MPI_DOUBLE,
			      MPI_SUM, MPI_COMM_WORLD);
	      }
	      MPI_Allreduce(MPI_IN_PLACE, meff[isym]->hab1[i],
			    nhab1, MPI_DOUBLE,
			    MPI_SUM, MPI_COMM_WORLD);
	      MPI_Allreduce(MPI_IN_PLACE, meff[isym]->hba1[i],
			    nhab1, MPI_DOUBLE,
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
	if (fn0 == NULL) {
	  AllocHamMem(h, meff[isym]->nbasis, meff[isym]->nbasis);
	  h->pj = isym;
	  h->n_basis = meff[isym]->nbasis;
	  h->dim = h->n_basis;
	  memcpy(h->basis, meff[isym]->basis, sizeof(int)*h->n_basis);
	  h->hsize = h->dim*(h->dim+1)/2;
	}
	k = ConstructHamilton(isym, nkg0, nkg0, kg, 0, NULL, 1);
	h0 = meff[isym]->h0;
	h->heff = heff;
	DecodePJ(isym, &pp, &jj);
	if (MyRankMPI() != 0) continue;   
	fwrite(&isym, sizeof(int), 1, f);
	if (meff[isym]->heff0) {
	  k = -meff[isym]->hsize0;
	  fwrite(&k, sizeof(int), 1, f);
	  fwrite(meff[isym]->heff0, sizeof(double), meff[isym]->hsize0, f);
	}
	fwrite(&(h->dim), sizeof(int), 1, f);
	fwrite(meff[isym]->basis, sizeof(int), h->dim, f);
	fflush(f);
	fflush(stdout);
	double *wb, *wc, *wbn, *wcn;
	int iw = 0, nw = 0;
	for (j = 0; j < h->dim; j++) {
	  for (i = 0; i <= j; i++) {
	    k = j*(j+1)/2 + i;
	    if (meff[isym]->imbpt[k] > 0) nw++;
	  }
	}
	if (nw > 0) {
	  wb = malloc(sizeof(double)*nw);
	  wc = malloc(sizeof(double)*nw);
	  wbn = malloc(sizeof(double)*nw);
	  wcn = malloc(sizeof(double)*nw);
	} else {
	  wb = NULL;
	  wc = NULL;
	  wbn = NULL;
	  wcn = NULL;
	}
	ResetWidMPI();      
#pragma omp parallel default(shared) private(dw, i, j, k, a, b, c, m, hab1, hba1, hab, hba, i0, iw)
	{
	  dw = malloc(sizeof(double)*(nr+nr2)*4);
	  //MPrintf(-1, "suminterph0: %d %d %x\n", isym, nw, dw);
	  iw = 0;
	  for (j = 0; j < h->dim; j++) {
	    for (i = 0; i <= j; i++) {
	      k = j*(j+1)/2 + i;
	      if (meff[isym]->imbpt[k] <= 0) {
		a = h0[k];
		m = j*h->dim + i;
		heff[m] = a;
		neff[m] = 0.0;
		if (meff[isym]->heff0) {
		  heff[m] += meff[isym]->heff0[m];
		}
		if (i < j) {
		  m = i*h->dim + j;
		  heff[m] = a;
		  neff[m] = 0.0;
		  if (meff[isym]->heff0) {
		    heff[m] += meff[isym]->heff0[m];
		  }
		}
		continue;
	      }
	      int skip = SkipMPI();
	      if (skip) {
		iw++;
		continue;
	      }
	      hab1 = meff[isym]->hab1[k];
	      hba1 = meff[isym]->hba1[k];
	      double *nab, *nba, *nab1, *nba1;
	      nab1 = hab1 + nhab1/2;
	      nba1 = hba1 + nhab1/2;
	      if (nhab > 0) {
		hab = meff[isym]->hab[k];
		hba = meff[isym]->hba[k];
		nab = hab + nhab/2;
		nba = hba + nhab/2;
	      } else {
		hab = NULL;
		hba = NULL;
		nab = NULL;
		nba = NULL;
	      }
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
	      b = SumInterpH(nr, ngr, nr2, ngr2, hab, hab1, dw, 1);
	      wbn[iw] = SumInterpH(nr, ngr, nr2, ngr2, nab, nab1, dw, 2);
	      neff[m] = wbn[iw];
	      heff[m] = a+b;
	      if (meff[isym]->heff0) {
		heff[m] += meff[isym]->heff0[m];
	      }
	      if (i < j) {
		m = i*h->dim + j;
		c = SumInterpH(nr, ngr, nr2, ngr2, hba, hba1, dw, 3);
		wcn[iw] = SumInterpH(nr, ngr, nr2, ngr2, nba, nba1, dw, 4);
		neff[m] = wcn[iw];
		heff[m] = a+c;
		if (meff[isym]->heff0) {
		  heff[m] += meff[isym]->heff0[m];
		}
	      } else {
		c = b;
		wcn[iw] = wbn[iw];
	      }
	      wb[iw] = b;
	      wc[iw] = c;
	      iw++;
	    }
	  }
	  //MPrintf(-1, "suminterph1: %d %d %x\n", isym, nw, dw);
	  free(dw);
	}
	
	iw = 0;
	for (j = 0; j < h->dim; j++) {
	  for (i = 0; i <= j; i++) {
	    k = j*(j+1)/2 + i;
	    if (meff[isym]->imbpt[k] <= 0) {
	      q = -i-1;
	      fwrite(&q, sizeof(int), 1, f);
	      q = -j-1;
	      fwrite(&q, sizeof(int), 1, f);
	      a = h0[k];
	      fwrite(&a, sizeof(double), 1, f);
	      a = 0.0;
	      fwrite(&a, sizeof(double), 1, f);
	      fwrite(&a, sizeof(double), 1, f);	    
	    } else {
	      fwrite(&i, sizeof(int), 1, f);
	      fwrite(&j, sizeof(int), 1, f);
	      a = h0[k];
	      fwrite(&a, sizeof(double), 1, f);
	      fwrite(&wb[iw], sizeof(double), 1, f);
	      fwrite(&wc[iw], sizeof(double), 1, f);
	      if (mbpt_savesum == 0) {
		hab1 = meff[isym]->hab1[k];
		hba1 = meff[isym]->hba1[k];
		if (nhab > 0) {
		  hab = meff[isym]->hab[k];
		  hba = meff[isym]->hba[k];
		} else {
		  hab = NULL;
		  hba = NULL;
		}
		fwrite(hab1, sizeof(double), nhab1, f);
		if (i != j) {
		  fwrite(hba1, sizeof(double), nhab1, f);
		}
		if (n3 != 1 && nhab > 0) {
		  fwrite(hab, sizeof(double), nhab, f);
		  if (i != j) {
		    fwrite(hba, sizeof(double), nhab, f);
		  }
		}
	      } else {
		fwrite(&wb[iw], sizeof(double), 1, f);
		fwrite(&wbn[iw], sizeof(double), 1, f);
		if (i != j) {
		  fwrite(&wc[iw], sizeof(double), 1, f);
		  fwrite(&wcn[iw], sizeof(double), 1,f );
		}	      
	      }
	      iw++;
	    }
	  }
	}
	if (nw > 0) {
	  free(wb);
	  free(wc);
	  free(wbn);
	  free(wcn);
	}
	fflush(f);
	double wt1 = WallTime();
	MPrintf(-1, "SumInterpH: %d %d %d %10.3E\n", isym, h->dim, nw, wt1-wt0);
      }
      if (icpf == icpi && fn != NULL && strlen(fn) > 0) {
	ResetWidMPI();   
#pragma omp parallel default(shared) private(isym, h, i, j, k, k0, a, c)
	{
	  for (isym = 0; isym < MAX_SYMMETRIES; isym++) { 
	    if (meff[isym] == NULL) continue;
	    int skip = SkipMPI();
	    if (skip) continue;
	    h = GetHamilton(isym);
	    if (h->dim <= 0) continue;
	    double wt1 = WallTime();
	    if (DiagnolizeHamilton(h) < 0) {
	      MPrintf(-1, "Diagnolizing Effective Hamiltonian Error\n");
	      ierr = -1;
	      continue;
	    }
	    // correct the normalization
	    double *ym = h->mixing + h->dim;
	    for (i = 0; i < h->dim; i++) {
	      a = 0.0;
	      for (j = 0; j < h->dim; j++) {
		for (k = 0; k < h->dim; k++) {
		  k0 = j*h->dim + k;
		  a += meff[isym]->neff[k0] * ym[j] * ym[k];
		}
	      }	      
	      if (a > 0 && a < 1) {
		c = sqrt(1.0 + a);
		for (j = 0; j < h->n_basis; j++) {
		  ym[j] /= c;
		}
	      }
	      ym += h->n_basis;
	    }
	    
	    AddToLevels(h, nkg0, kg);
	    double wt2 = WallTime();
	    MPrintf(-1, "Diag Ham: %3d %5d %5d %10.4E %10.4E %10.4E %10.4E\n",
		    isym, h->dim, h->n_basis, wt2-wt1, wt2-tbg,
		    TotalSize(), TotalArraySize());
	    fflush(stdout);
	  }
	}
	if (MyRankMPI() == 0) {
	  SortLevels(nlevels, -1, 0);
	  SaveLevels(fn, nlevels, -1);
	}
	for (isym = 0; isym < MAX_SYMMETRIES; isym++) {
	  h = GetHamilton(isym);
	  if (h->dim <= 0) continue;
	  h->heff = NULL;
	  AllocHamMem(h, -1, -1);
	  AllocHamMem(h, 0, 0);
	}
	if (MyRankMPI() == 0) {
	  tt1 = WallTime();
	  dt = tt1 - tbg;
	  tt0 = tt1;
	  MPrintf(-1, "End MBPT Structure: %12.5E %12.5E %12.5E %12.5E %12.5E %ld\n", dt, TotalSize(), TotalArraySize(), ttskip, ttlock, tnlock);
	  fflush(stdout);
	}
      }
      fclose(f);
    }
    free(cfgpair);
  }
#if USE_MPI == 1
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  if (mbpt_tr.nktr > 0) {
    double ts = TotalSize();
    MPrintf(-1, "MBPT Transition: %12.5E %12.5E\n", WallTime()-tbg, ts);
    fflush(stdout);
    if (mbpt_freetr > 0 && ts > mbpt_freetr) {
      ReinitRadial(2);
      MPrintf(-1, "MBPT Transition ReinitRadial: %12.5E %12.5E\n", WallTime()-tbg, TotalSize());
    }
    if (MyRankMPI() == 0) {
      sprintf(tfn, "%s.tr", fn1);
      f = fopen(tfn, "w");
    }
    CONFIG_PAIR *cfgpair0 = malloc(sizeof(CONFIG_PAIR)*nc*nc);
    int *ncpti = malloc(sizeof(int)*npr);
    for (i = 0; i < npr; i++) ncpti[i] = 0;
    ResetWidMPI();
#pragma omp parallel default(shared) private(k0, c0, k1, c1, m, m0, ms0, i0, q0, p0, j0, pp0, pp1, m1, ms1, i1, q1, p1, j1, i, j, ic, mr)
    {
    mr = MPIRank(NULL);
    int w = 0;
    for (k0 = 0; k0 < nc; k0++) {
      c0 = cs[k0];
      for (k1 = 0; k1 < nc; k1++) {
	c1 = cs[k1];
	if (SkipWMPI(w++)) continue;
	ic = k0 + k1*nc;
	cfgpair0[ic].m = 0;
	m = 0;      
	for (m0 = 0; m0 < c0->n_csfs; m0++) {
	  ms0 = c0->symstate[m0];
	  UnpackSymStateMBPT(meff, ms0, &i0, &q0);
	  if (q0 < 0) continue;
	  DecodePJ(i0, &p0, &j0);
	  pp0 = i0*2*mbpt_tr.nktr;
	  pp1 = pp0 + 2*mbpt_tr.nktr;
	  for (m1 = 0; m1 < c1->n_csfs; m1++) {
	    ms1 = c1->symstate[m1];
	    UnpackSymStateMBPT(meff, ms1, &i1, &q1);
	    if (q1 < 0) continue;
	    DecodePJ(i1, &p1, &j1);
	    for (i = pp0; i < pp1; i++) {
	      j = IdxGet(mtr[i].ids, i1);
	      if (j < 0) continue;
	      if (mtr[i].pma[j][q0*mtr[i].ndim1[j]+q1] == NULL) continue;
	      if (!Triangle(j0, j1, 2*abs(mtr[i].m))) continue;
	      m++;
	      break;
	    }
	  }
	}
	if (m == 0) continue;
	cfgpair0[ic].k0 = k0;
	cfgpair0[ic].k1 = k1;
	cfgpair0[ic].m = m;
	ncpti[mr]++;
      }
    }
    }
    ncpt = 0;
    for (i = 0; i < npr; i++) {
      ncpt += ncpti[i];
    }
    MPrintf(-1, "MBPT ncpt: %d %12.5E %12.5E\n",
	    ncpt, WallTime()-tbg, TotalSize());
    if (ncp > ncpt) {
      printf("transition ncp > ncpt: %d %d\n", ncp, ncpt);
      ncp = ncpt;
    }
    if (ncp > 1 && (icp < 0 || icp >= ncp)) {
      printf("invalid config splits for transition: %d %d\n", icp, ncp);
      Abort(1);
    }    
    cfgpair = malloc(sizeof(CONFIG_PAIR)*ncpt);
    RANDIDX *rid = NULL;
    if (mbpt_rand) rid = RandList(ncpt);
    ic = 0;
    mst = 0;
    for (k0 = 0; k0 < nc; k0++) {
      for (k1 = 0; k1 < nc; k1++) {
	i = k0 + k1*nc;
	m = cfgpair0[i].m;
	if (m <= 0) continue;
	mst += m;
	int icr = ic;
	if (rid) icr = rid[ic].i;
	cfgpair[icr].k0 = cfgpair0[i].k0;
	cfgpair[icr].k1 = cfgpair0[i].k1;
	cfgpair[icr].m = m;
	ic++;
      }
    }
    free(rid);
    free(cfgpair0);
    free(ncpti);
    icp = icpi;
    int tmst = GetICP(ncpt, cfgpair, ncp, icp, &icp0, &icp1);
    MPrintf(-1, "MBPT transition cfgpair: %d %d %d %d %d %d %d %12.5E %12.5E\n",
	    icp, ncp, icp0, icp1, ncpt, tmst, mst, WallTime()-tbg, TotalSize());
    ncps = 0;
    mbpt_ignoren = 0;
    if (mbpt_omp0 == 2 || mbpt_omp0 == 3) {
      mbpt_omp = 1;
    } else {
      mbpt_omp = 0;
    }
    double *dm = NULL;
    int *im = NULL;      
    if (mbpt_msort) {
      dm = malloc(sizeof(double)*(icp1-icp0));
      im = malloc(sizeof(int)*(icp1-icp0));
      double amst = tmst/npr;
      for (ic = icp0; ic < icp1; ic++) {
	dm[ic-icp0] = amst/cfgpair[ic].m;
	if (cfgpair[ic].k0 != cfgpair[ic].k1) {
	  dm[ic-icp0] *= mbpt_asort;
	}
      }
      ArgSort(icp1-icp0, dm, im);
    }
    MPrintf(-1, "MBPT transition beg: %12.5E %12.5E\n", WallTime()-tbg, TotalSize());
    ResetWidMPI();
#pragma omp parallel default(shared) private(n0,bra,ket,sbra,sket,bra1,ket1,bra2,ket2,sbra1,sket1,sbra2,sket2,cs,dt,dtt,k0,k1,c0,p0,c1,p1,m,bst0,kst0,m0,m1,ms0,ms1,q,q0,q1,k,mst,i0,i1,j0,j1,ct0,ct1,bst,kst,n1,bas0,bas1,ic,pp0,pp1,i,j)
    {
      int cmst = 0;
      double ptt0, ptt1;
      cs = mbpt_cs;
      bas0 = mbpt_bas0;
      bas1 = mbpt_bas1;
      int ic0;
      for (ic0 = icp0; ic0 < icp1; ic0++) {
	int skip = SkipMPIM(0);
	if (skip) continue;
	if (im) {
	  ic = icp0 + im[ic0-icp0];
	} else {
	  ic = ic0;
	}
	k0 = cfgpair[ic].k0;
	k1 = cfgpair[ic].k1;
	c0 = cs[k0];
	c1 = cs[k1];
	bst0 = malloc(sizeof(int)*cfgpair[ic].m);
	kst0 = malloc(sizeof(int)*cfgpair[ic].m);
	m = 0;
	for (m0 = 0; m0 < c0->n_csfs; m0++) {
	  ms0 = c0->symstate[m0];
	  UnpackSymStateMBPT(meff, ms0, &i0, &q0);
	  if (q0 < 0) continue;
	  DecodePJ(i0, &p0, &j0);
	  pp0 = i0*2*mbpt_tr.nktr;
	  pp1 = pp0 + 2*mbpt_tr.nktr;
	  for (m1 = 0; m1 < c1->n_csfs; m1++) {
	    ms1 = c1->symstate[m1];
	    UnpackSymStateMBPT(meff, ms1, &i1, &q1);
	    if (q1 < 0) continue;
	    DecodePJ(i1, &p1, &j1);
	    for (i = pp0; i < pp1; i++) {
	      j = IdxGet(mtr[i].ids, i1);
	      if (j < 0) continue;
	      if (mtr[i].pma[j][q0*mtr[i].ndim1[j]+q1] == NULL) continue;
	      if (!Triangle(j0, j1, 2*abs(mtr[i].m))) continue;
	      bst0[m] = m0;
	      kst0[m] = m1;
	      m++;
	      break;
	    }
	  }
	}
	ptt0 = tt0;
	/* mst pairs */
	mst = m;
	cmst += m;
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
	DeltaH12M0(mtr, n0, bra2, ket2, sbra2, sket2, mst, bst, kst,
		   ct0, ct1, &mbpt_ibas0, mbpt_bas0s, mbpt_bas0d,
		   &ing, nc, cs, 1);	    
	/* 1-b 2-b term 1 virtual orb */
	DeltaH12M1(mtr, n0+1, bra1, ket1, sbra1, sket1, mst, bst, kst,
		   ct0, ct1, &mbpt_ibas0, mbpt_bas0s, mbpt_bas0d,
		   &mbpt_ibas1, &ing, nc, cs, 1);
	/* 1-b 1-b term no virtual orb */
	DeltaH11M0(mtr, n0, bra2, ket2, sbra2, sket2, mst, bst, kst,
		   ct0, ct1, &mbpt_ibas0, mbpt_bas0s, &ing, nc, cs, 1);
	/* 1-b 1-b term 1 virtual orb */
	DeltaH11M1(mtr, n0+1, bra1, ket1, sbra1, sket1, mst, bst, kst, 
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
	MPrintf(mbpt_omp?0:-1, "%4d %4d %4d %4d %3d %3d %3d %3d %3d %3d ... %12.5E %12.5E %12.5E %ld\n", 
		ic0, ic, icp0, icp1, k0, k1, nc, mst, cmst, tmst,
		dt, dtt, tmem, mbpt_ignoren);
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
    if (im) {
      free(dm);
      free(im);
    }
    free(cfgpair);
    k = 2*mbpt_tr.nktr*MAX_SYMMETRIES;
#if USE_MPI == 1
    if (NProcMPI() > 1) {
      for (j = 0; j < k; j++) {
	for (i = 0; i < mtr[j].nsym1; i++) {
	  for (m0 = 0; m0 < mtr[j].ndim0; m0++) {
	    for (m1 = 0; m1 < mtr[j].ndim1[i]; m1++) {
	      p0 = m0*mtr[j].ndim1[i] + m1;
	      if (mtr[j].pma[i][p0] == NULL) continue;	      
	      m = meff[mtr[j].isym0]->n*mbpt_tr.naw;	      
	      MPI_Allreduce(MPI_IN_PLACE, mtr[j].pma[i][p0]->tma, m,
			    MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	      MPI_Allreduce(MPI_IN_PLACE, mtr[j].pma[i][p0]->tma, m,
			    MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	    }
	  }
	}
      }
    }
#endif
    if (MyRankMPI() == 0) {
      fwrite(&mbpt_tr.nktr, sizeof(int), 1, f);
      fwrite(&mbpt_tr.naw, sizeof(int), 1, f);
      fwrite(&emin, sizeof(double), 1, f);
      fwrite(&emax, sizeof(double), 1, f);
      for (j = 0; j < k; j++) {
	fwrite(&mtr[j].m, sizeof(int), 1, f);
	fwrite(&mtr[j].nsym1, sizeof(int), 1, f);
	fwrite(&mtr[j].isym0, sizeof(int), 1, f);
	fwrite(&mtr[j].ndim0, sizeof(int), 1, f);
	if (mtr[j].nsym1 > 0) {
	  fwrite(mtr[j].isym1, sizeof(int), mtr[j].nsym1, f);
	  fwrite(mtr[j].ndim1, sizeof(int), mtr[j].nsym1, f);
	}
      }
      if (mbpt_savesum == 0) {
	for (j = 0; j < k; j++) {
	  if (mtr[j].nsym1 > 0) {
	    for (i = 0; i < mtr[j].nsym1; i++) {
	      for (m0 = 0; m0 < mtr[j].ndim0; m0++) {
		for (m1 = 0; m1 < mtr[j].ndim1[i]; m1++) {
		  p0 = m0*mtr[j].ndim1[i] + m1;
		  if (mtr[j].pma[i][p0] == NULL) {
		    m = 0;
		  } else {
		    m = meff[mtr[j].isym0]->n*mbpt_tr.naw;
		  }
		  fwrite(&m, sizeof(int), 1, f);
		  if (m > 0) {
		    fwrite(mtr[j].pma[i][p0]->tma, sizeof(double), m, f);
		    fwrite(mtr[j].pma[i][p0]->rma, sizeof(double), m, f);
		  }
		}
	      }
	    }
	  }	    
	}
      }
    }
    ResetWidMPI();
#pragma omp parallel default(shared) private(q, i, j, m0, m1, p0, q0, q1)
    {
      double *x = malloc(sizeof(double)*nr*5);
      double *t = x + nr;
      double *y = t + nr;
      double *zz1 = y + nr;
      double *zz2 = zz1 + nr;
      for (q = 0; q < nr; q++) {
	x[q] = ngr[q];
      }
      for (j = 0; j < k; j++) {
	double iwt0 = WallTime();
	for (i = 0; i < mtr[j].nsym1; i++) {
	  for (m0 = 0; m0 < mtr[j].ndim0; m0++) {
	    for (m1 = 0; m1 < mtr[j].ndim1[i]; m1++) {
	      int skip = SkipMPI();
	      if (skip) continue;
	      p0 = m0*mtr[j].ndim1[i] + m1;
	      if (mtr[j].pma[i][p0] == NULL) continue;
	      for (q0 = 0; q0 < mbpt_tr.naw; q0++) {
		for (q = 0; q < nr; q++) {
		  q1 = q0*nr + q;
		  zz1[q] = mtr[j].pma[i][p0]->tma[q1];
		  zz2[q] = mtr[j].pma[i][p0]->rma[q1];
		}
		mtr[j].pma[i][p0]->tma[q0] = SumInterp1D(nr, zz1, x, t, y, 5);
		mtr[j].pma[i][p0]->rma[q0] = SumInterp1D(nr, zz2, x, t, y, 6);
	      }
	    }
	  }
	}
	double iwt1 = WallTime();
	if (mtr[j].nsym1 > 0) {
	  MPrintf(0, "SumInterpTR: %4d %4d %4d %4d %12.5E\n",
		  j, mtr[j].isym0, mtr[j].m, mtr[j].nsym1, iwt1-iwt0);
	}
      }
      free(x);
    }
    tt1 = WallTime();
    dt = tt1-tbg;
    MPrintf(-1, "write MBPT transition: %12.5E %12.5E %ld\n",
	    dt, TotalSize(), mbpt_ignoren);
    if (mbpt_savesum && MyRankMPI() == 0) {
      for (j = 0; j < k; j++) {
	for (i = 0; i < mtr[j].nsym1; i++) {
	  for (m0 = 0; m0 < mtr[j].ndim0; m0++) {
	    for (m1 = 0; m1 < mtr[j].ndim1[i]; m1++) {
	      p0 = m0*mtr[j].ndim1[i] + m1;
	      if (mtr[j].pma[i][p0] == NULL) {
		m = 0;
	      } else {
		m = mbpt_tr.naw;
	      }
	      fwrite(&m, sizeof(int), 1, f);
	      if (m > 0) {
		fwrite(mtr[j].pma[i][p0]->tma, sizeof(double), m, f);
		fwrite(mtr[j].pma[i][p0]->rma, sizeof(double), m, f);
	      }
	    }
	  }
	}
      }
    }
    if (mbpt_tr.nlow > 0 && mbpt_tr.nup > 0) {
      tt1 = WallTime();
      dt = tt1-tbg;
      MPrintf(-1, "SaveTransitionMBPT: %12.5E %12.5E\n", dt, TotalSize());
      SetTransitionMode(0);
      SaveTransitionMBPT(mtr);
      printf("AdjustAngZ ... %12.5E %12.5E\n", WallTime()-tbg, TotalSize());
      AdjustAngularZ(mtr);
      printf("End TransitionMBPT ... %12.5E %12.5E\n",
	     WallTime()-tbg, TotalSize());
    }
    
    tt1 = WallTime();
    dt = tt1 - tbg;
    tt0 = tt1;
    MPrintf(-1, "Total Time Transition = %12.5E %12.5E %12.5E\n",
	    dt, TotalSize(), TotalArraySize());
    fflush(stdout);
    if (MyRankMPI() == 0) fclose(f);
  }
  
 ERROR:
  FreeIdxAry(&mbpt_ibas0, 2);
  FreeIdxAry(&mbpt_ibas1, 2);
  FreeIdxAry(&ing, 2);
  FreeIdxAry(&ing2, 2);
  CONFIG **cca = NULL;
  int icca = 0;
  if (ncca > 0) {
    cca = malloc(sizeof(CONFIG *)*ncca);
  }
  ResetWidMPI();
#pragma omp parallel default(shared) private(cs, i)
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
    if (mbpt_ccn[0]) {
#pragma omp critical
      {
	for (i = 0; i < mbpt_cca->dim; i++) {
	  CONFIG **c = ArrayGet(mbpt_cca, i);
	  cca[icca] = *c;
	  icca++;
	}
	ArrayFree(mbpt_cca, NULL);
	free(mbpt_cca);
	mbpt_cca = NULL;
      }
    }
  }
  if (ncca > 0) {
    FILE *fc = fopen(mbpt_ccn, "w");
    if (fc == NULL) {
      printf("cannot open mbpt_ccn: %s\n", mbpt_ccn);
    } else {
      qsort(cca, ncca, sizeof(CONFIG *), CompareMBPTCC);
      char scr[2048];
      icca = 0;
      for (i = 0; i < ncca; i++) {
	if (i > 0 && 0==CompareCfgPointer(&cca[i], &cca[i-1])) continue;
	ConstructNRConfigName(scr, 2048, cca[i]);
	fprintf(fc, "%s # %d %d %d %d %d %d %d %g %g %g %g %g %g\n",
		scr, icca, i, cca[i]->icfg,
		cca[i]->n_electrons, cca[i]->n_csfs,
		cca[i]->nnrs, cca[i]->igroup,
		cca[i]->sweight, cca[i]->energy,
		cca[i]->delta, cca[i]->sth, cca[i]->cth, cca[i]->mde);
	icca++;
      }
      for (i = 0; i < ncca; i++) {
	cca[i]->n_csfs = 0;
	cca[i]->nnrs = 0;
	FreeConfigData(cca[i]);
	free(cca[i]);
      }
      fclose(fc);
      free(cca);
    }
  }    
  free(bas);
  if (mbpt_nsplit) {
    for (i = 0; i < nb; i++) {
      free(mbpt_rij[i]);
    }
    free(mbpt_rij);
  }
  FreeEffMBPT(meff);
  FreeTransitionMBPT(mtr);
  if (ngr != ng) free(ngr);
  if (ngr2 != ng2) free(ngr2);
  free(ng);
  free(ng2);
  if (rh && nkg > 0) {
    free(kg);
  }
  for (i = 0; i < nmk2; i++) {
    ArrayFree(&mbpt_csary[i], NULL);
  }
  free(mbpt_csary);
  return ierr;
}

int ReadMBPT(int nf, FILE *f[], MBPT_HAM *mbpt, int m) {
  int i, j, k, q, nb, nn2, n, mst, m0, m1, p0;
  double emin, emax;
  MBPT_TR *mtr;
    
  if (m == 0) {
    for (i = 0; i < nf; i++) {
      nb = fread(&(mbpt[i].n), sizeof(int), 1, f[i]);
      mbpt[i].ng = malloc(sizeof(int)*mbpt[i].n);
      nb = fread(mbpt[i].ng, sizeof(int), mbpt[i].n, f[i]);
      nb = fread(&(mbpt[i].n2), sizeof(int), 1, f[i]);
      if (mbpt[i].n2 > 0) {
	mbpt[i].ng2 = malloc(sizeof(int)*mbpt[i].n2);
	nb = fread(mbpt[i].ng2, sizeof(int), mbpt[i].n2, f[i]);
      }
      nb = fread(&(mbpt[i].n3), sizeof(int), 1, f[i]);
      mbpt[i].hab1 = malloc(sizeof(double)*mbpt[i].n*2);
      mbpt[i].hba1 = malloc(sizeof(double)*mbpt[i].n*2);
      if (mbpt[i].n2 > 0) {
	mbpt[i].hab = malloc(sizeof(double)*mbpt[i].n*mbpt[i].n2*2);
	mbpt[i].hba = malloc(sizeof(double)*mbpt[i].n*mbpt[i].n2*2);
      }
      //printf("r0: %d %d %d %d %d\n", i, mbpt[i].n, mbpt[i].n2, mbpt[i].n3, mbpt[i].ng[0]);
    }
    return 0;
  }

  if (m == 1) {
    for (i = 0; i < nf; i++) {
      nb = fread(&(mbpt[i].isym), sizeof(int), 1, f[i]);
      //if (nb != 1) return -1;
      //printf("isym: %d %d\n", nb, mbpt[i].isym);
      nb = fread(&(mbpt[i].dim), sizeof(int), 1, f[i]);
      //printf("dim: %d %d\n", nb, mbpt[i].dim);
      if (mbpt[i].dim < 0) {
	mbpt[i].hsize0 = -mbpt[i].dim;
	mbpt[i].heff0 = malloc(sizeof(double)*mbpt[i].hsize0);
	nb = fread(mbpt[i].heff0, sizeof(double), mbpt[i].hsize0, f[i]);
	nb = fread(&(mbpt[i].dim), sizeof(int), 1, f[i]);
      } else {
	mbpt[i].hsize0 = 0;
	mbpt[i].heff0 = NULL;
      }
      if (mbpt[i].dim > 0) {
	mbpt[i].basis = malloc(sizeof(int)*mbpt[i].dim);
	nb = fread(mbpt[i].basis, sizeof(int), mbpt[i].dim, f[i]);
      }      
    }
    return 0;
  }
  
  if (m == 2) {
    for (i = 0; i < nf; i++) {
      //printf("r2: %d\n", mbpt[i].dim);
      if (mbpt[i].dim == 0) continue;
      nb = fread(&(mbpt[i].ibra), sizeof(int), 1, f[i]);
      //if (nb != 1) return -1;
      nb = fread(&(mbpt[i].iket), sizeof(int), 1, f[i]);
      //if (nb != 1) return -1;
      nb = fread(&(mbpt[i].a), sizeof(double), 1, f[i]);
      //if (nb != 1) return -1;
      nb = fread(&(mbpt[i].b), sizeof(double), 1, f[i]);
      //if (nb != 1) return -1;
      nb = fread(&(mbpt[i].c), sizeof(double), 1, f[i]);
      //if (nb != 1) return -1;
      //printf("r1: %d %d %d %g %g %g\n", i, mbpt[i].ibra, mbpt[i].iket, mbpt[i].a, mbpt[i].b, mbpt[i].c);
      if (mbpt[i].ibra < 0 || mbpt[i].iket < 0) {
	continue;
      }
      nb = fread(mbpt[i].hab1, sizeof(double), mbpt[i].n*2, f[i]);
      //if (nb != mbpt[i].n*2) return -1;
      if (mbpt[i].ibra != mbpt[i].iket) {
	nb = fread(mbpt[i].hba1, sizeof(double), mbpt[i].n*2, f[i]);
	//if (nb != mbpt[i].n*2) return -1;
      } else {
	memcpy(mbpt[i].hba1, mbpt[i].hab1, sizeof(double)*mbpt[i].n*2);
      }
      if (mbpt[i].n3 != 1 && mbpt[i].n2 > 0) {
	nn2 = 2*mbpt[i].n*mbpt[i].n2;
	nb = fread(mbpt[i].hab, sizeof(double), nn2, f[i]);
	//if (nb != nn2) return -1;
	if (mbpt[i].ibra != mbpt[i].iket) {
	  nb = fread(mbpt[i].hba, sizeof(double), nn2, f[i]);
	  //if (nb != nn2) return -1;
	} else {
	  memcpy(mbpt[i].hba, mbpt[i].hab, sizeof(double)*nn2);
	}
      }
    }
    return 0;
  }
  
  if (m == 3) {
    for (i = 0; i < nf; i++) {
      nb = fread(&mbpt_tr.nktr, sizeof(int), 1, f[i]);
      nb = fread(&mbpt_tr.naw, sizeof(int), 1, f[i]);
      nb = fread(&emin, sizeof(double), 1, f[i]);
      nb = fread(&emax, sizeof(double), 1, f[i]);
      k = 2*mbpt_tr.nktr*MAX_SYMMETRIES;
      mtr = malloc(sizeof(MBPT_TR)*k);
      for (j = 0; j < k; j++) {
	nb = fread(&mtr[j].m, sizeof(int), 1, f[i]);
	nb = fread(&mtr[j].nsym1, sizeof(int), 1, f[i]);
	nb = fread(&mtr[j].isym0, sizeof(int), 1, f[i]);
	nb = fread(&mtr[j].ndim0, sizeof(int), 1, f[i]);
	if (mtr[j].nsym1 > 0) {
	  mtr[j].isym1 = malloc(sizeof(int)*mtr[j].nsym1);
	  mtr[j].ndim1 = malloc(sizeof(int)*mtr[j].nsym1);
	  mtr[j].pma = malloc(sizeof(MBPT_PMA **)*mtr[j].nsym1);
	  for (q = 0; q < mtr[j].nsym1; q++) {
	    mtr[j].pma[q] = NULL;
	  }
	  nb = fread(mtr[j].isym1, sizeof(int), mtr[j].nsym1, f[i]);
	  nb = fread(mtr[j].ndim1, sizeof(int), mtr[j].nsym1, f[i]);	  
	  mtr[j].ids = malloc(sizeof(IDXARY));
	  InitIdxAry(mtr[j].ids, mtr[j].nsym1, mtr[j].isym1);
	} else {
	  mtr[j].ids = NULL;
	  mtr[j].isym1 = NULL;
	  mtr[j].ndim1 = NULL;
	  mtr[j].pma = NULL;
	}
      }
      mbpt[i].mtr = mtr;
    }
    SetAWGridMBPT(emin, emax);
    return 0;
  }

  if (m < 0) {
    j = -m-1;
    for (i = 0; i < nf; i++) {
      mtr = mbpt[i].mtr;
      for (q = 0; q < mtr[j].nsym1; q++) {
	mst = mtr[j].ndim0 * mtr[j].ndim1[q];
	mtr[j].pma[q] = malloc(sizeof(MBPT_PMA *)*mst);
	for (m0 = 0; m0 < mtr[j].ndim0; m0++) {
	  for (m1 = 0; m1 < mtr[j].ndim1[q]; m1++) {
	    p0 = m0*mtr[j].ndim1[q] + m1;
	    fread(&n, sizeof(int), 1, f[i]);
	    if (n == 0) {
	      mtr[j].pma[q][p0] = NULL;
	    } else {
	      mtr[j].pma[q][p0] = malloc(sizeof(MBPT_PMA));
	      mtr[j].pma[q][p0]->tma = malloc(sizeof(double)*n);
	      mtr[j].pma[q][p0]->rma = malloc(sizeof(double)*n);
	      fread(mtr[j].pma[q][p0]->tma, sizeof(double), n, f[i]);
	      fread(mtr[j].pma[q][p0]->rma, sizeof(double), n, f[i]);
	    }
	    //printf("rmtr: %d %d %d %d %d %d %d %x\n", j, i, q, m0, m1, p0, n, mtr[j].pma[q][p0]);
	  }
	}
      }
    }
    return 0;
  }
  return 0;
}

void CombineMBPT0(int nf, MBPT_HAM *mbpt, 
		  double *hab1, double *hba1,
		  double **hab, double **hba, 
		  double *nab1, double *nba1,
		  double **nab, double **nba,
		  int n0, int *ng0,
		  int n, int *ng,
		  int *n2, int **ng2) {
  ResetWidMPI();   
#pragma omp parallel default(shared)
  {
  int m, i, j, k, q, r, nn2;
#pragma omp for schedule(static)
  for (m = 0; m < nf; m++) {
    if (mbpt[m].ibra < 0 || mbpt[m].iket < 0) continue;
    for (i = 0; i < mbpt[m].n; i++) {
      //k = IdxGet(ing0, mbpt[m].ng[i]);
      k = IBisect(mbpt[m].ng[i], n0, ng0);
      if (k < 0) continue;
      hab1[k] = mbpt[m].hab1[i];
      hba1[k] = mbpt[m].hba1[i];
      nab1[k] = mbpt[m].hab1[i+mbpt[m].n];
	nba1[k] = mbpt[m].hba1[i+mbpt[m].n];
    }
    if (mbpt[m].n3 != 1) {
      nn2 = mbpt[m].n * mbpt[m].n2;
      for (i = 0; i < mbpt[m].n; i++) {
	//k = IdxGet(ing, mbpt[m].ng[i]);
	k = IBisect(mbpt[m].ng[i], n, ng);
	if (k < 0) continue;
	r = i*mbpt[m].n2;
	for (j = 0; j < mbpt[m].n2; j++) {
	  //q = IdxGet(&ing2[k], mbpt[m].ng2[j]);
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
}

void CombineMBPT(int nf, MBPT_HAM *mbpt, 
		 double *hab1, double *hba1,
		 double **hab, double **hba, 
		 double *nab1, double *nba1,
		 double **nab, double **nba,
		 IDXARY *ing0, IDXARY *ing, IDXARY *ing2) {
  int m, i, j, k, q, r, nn2;
  m = 0;
  for (i = 0; i < ing0->n; i++) {
    hab1[i] = 0;
    hba1[i] = 0;
    nab1[i] = 0;
    nba1[i] = 0;
  }
  if (mbpt[m].n3 != 1 && mbpt[m].n2 > 0) {
    for (i = 0; i < ing->n; i++) {
      for (j = 0; j < ing2[i].n; j++) {
	hab[i][j] = 0;
	hba[i][j] = 0;
	nab[i][j] = 0;
	nba[i][j] = 0;
      }
    }
  }
  ResetWidMPI();   
#pragma omp parallel default(shared) private(m,i,j,k,q,r,nn2)
  {
  int w = 0;
  for (m = 0; m < nf; m++) {
    if (SkipWMPI(w++)) continue;
    if (mbpt[m].ibra < 0 || mbpt[m].iket < 0) continue;
    for (i = 0; i < mbpt[m].n; i++) {
      k = IdxGet(ing0, mbpt[m].ng[i]);
      //k = IBisect(mbpt[m].ng[i], n0, ng0);
      if (k < 0) continue;
#pragma omp atomic
      hab1[k] += mbpt[m].hab1[i];
#pragma omp atomic
      hba1[k] += mbpt[m].hba1[i];
#pragma omp atomic
      nab1[k] += mbpt[m].hab1[i+mbpt[m].n];
#pragma omp atomic
      nba1[k] += mbpt[m].hba1[i+mbpt[m].n];
    }
    if (mbpt[m].n3 != 1 && mbpt[m].n2 > 0) {
      nn2 = mbpt[m].n * mbpt[m].n2;
      for (i = 0; i < mbpt[m].n; i++) {
	k = IdxGet(ing, mbpt[m].ng[i]);
	//k = IBisect(mbpt[m].ng[i], n, ng);
	if (k < 0) continue;
	r = i*mbpt[m].n2;
	for (j = 0; j < mbpt[m].n2; j++) {
	  q = IdxGet(&ing2[k], mbpt[m].ng2[j]);
	  //q = IBisect(mbpt[m].ng2[j], n2[k], ng2[k]);
	  if (q < 0) continue;
#pragma omp atomic
	  hab[k][q] += mbpt[m].hab[r+j];
#pragma omp atomic
	  hba[k][q] += mbpt[m].hba[r+j];
#pragma omp atomic
	  nab[k][q] += mbpt[m].hab[r+j + nn2];
#pragma omp atomic
	  nba[k][q] += mbpt[m].hba[r+j + nn2];
	}
      }
    }
  }
  }
}

void CombineTransitionMBPT(int nf, MBPT_HAM *mbpt, MBPT_TR *mtr,
			   int j, IDXARY *ing) {
  int m, i, t, q, r, p, is0, is1, n;

  m = 0;
  n = ing->n;
  for (q = 0; q < mtr[j].nsym1; q++) {
    for (is0 = 0; is0 < mtr[j].ndim0; is0++) {
      for (is1 = 0; is1 < mtr[j].ndim1[q]; is1++) {
	t = is0*mtr[j].ndim1[q]+is1;
	if (mbpt[0].mtr[j].pma[q][t]) {
	  mtr[j].pma[q][t] = malloc(sizeof(MBPT_PMA));
	  mtr[j].pma[q][t]->tma = malloc(sizeof(double)*n*mbpt_tr.naw);
	  mtr[j].pma[q][t]->rma = malloc(sizeof(double)*n*mbpt_tr.naw);
	  for (p = 0; p < mbpt_tr.naw; p++) {
	    for (i = 0; i < mbpt[m].n; i++) {
	      r = IdxGet(ing, mbpt[m].ng[i]);
	      if (r < 0) continue;
	      mtr[j].pma[q][t]->tma[p*n+r] = 0;
	      mtr[j].pma[q][t]->rma[p*n+r] = 0;
	    }
	  }
	}
      }
    }
  }
  ResetWidMPI();   
#pragma omp parallel default(shared) private(m,i,t,q,r,p,is0,is1,n)
  {
    n = ing->n;
    int w = 0;
    for (m = 0; m < nf; m++) {
      for (q = 0; q < mtr[j].nsym1; q++) {
	for (is0 = 0; is0 < mtr[j].ndim0; is0++) {
	  for (is1 = 0; is1 < mtr[j].ndim1[q]; is1++) {
	    if (SkipWMPI(w++)) continue;
	    t = is0*mtr[j].ndim1[q]+is1;
	    if (mtr[j].pma[q][t] == NULL) continue;
	    for (p = 0; p < mbpt_tr.naw; p++) {
	      for (i = 0; i < mbpt[m].n; i++) {
		r = IdxGet(ing, mbpt[m].ng[i]);
		if (r < 0) continue;
#pragma omp atomic
		mtr[j].pma[q][t]->tma[p*n+r] += mbpt[m].mtr[j].pma[q][t]->tma[p*mbpt[m].n+i];
#pragma omp atomic
		mtr[j].pma[q][t]->rma[p*n+r] += mbpt[m].mtr[j].pma[q][t]->rma[p*mbpt[m].n+i];
	      }
	    }
	  }
	}
      }
    }
  }
  for (m = 0; m < nf; m++) {
    FreeTransitionDataMBPT(mbpt[m].mtr, j);
  }
}

void AdjustAngularZ(MBPT_TR *mtr) {
  SHAMILTON *h;
  ANGZ_DATUM *ad;
  //ANGZ_ARY **a;
  int nh, i0, i1, m0, m1, k0, k1, mktr;
  int i, k, iz, n, ns, p, q, p0, p1, j0, j1, ik;

  if (mbpt_adjaz <= 0) return;
  mktr = 0;
  for (i = 0; i < mbpt_tr.nktr; i++) {
    if (mbpt_tr.ks[i] > mktr) mktr = mbpt_tr.ks[i];
  }
  mktr = Min(mktr, mbpt_adjaz);
  SetMaxKMBPT(mktr);
  h = GetSHamilton(&nh);
  ResetWidMPI();
#pragma omp parallel default(shared) private(i0, i1, p0, j0, p1, j1, ns, n, i, k, m0, m1, iz, ik, k0, k1, p, q, ad)
  {
  for (i0 = 0; i0 < nh; i0++) {
    DecodePJ(h[i0].pj, &p0, &j0);
    for (i1 = i0; i1 < nh; i1++) {
      DecodePJ(h[i1].pj, &p1, &j1);
      int skip = SkipMPI();
      if (skip) continue;
      ns = AngularZMixStates(&ad, i0, i1);
      if (ns <= 0) continue;
      ns = ad->nd;
      //a = (ANGZ_ARY **) ad->angz;
      ad->mk = malloc(sizeof(double *)*ns*2);
      n = mktr * mbpt_tr.naw;
      for (i = 0; i < 2*ns; i++) {
	ad->mk[i] = malloc(sizeof(double)*n);
	for (k = 0; k < n; k++) {
	  ad->mk[i][k] = 0.0;
	}
      }
      for (m0 = 0; m0 < h[i0].nlevs; m0++) {
	for (m1 = 0; m1 < h[i1].nlevs; m1++) {
	  iz = m0*h[i1].nlevs + m1;
	  for (ik = 0; ik < mbpt_tr.nktr; ik++) {
	    k = mbpt_tr.ks[ik]-1;
	    if (k >= mktr) continue;
	    n = mbpt_tr.naw * k;
	    k0 = h[i0].pj*2*mbpt_tr.nktr + 2*ik;
	    if (IsEven(p0+p1+k)) k0++;
	    //k1 = IBisect(h[i1].pj, mtr[k0].nsym1, mtr[k0].isym1);
	    k1 = IdxGet(mtr[k0].ids, h[i1].pj);
	    if (k1 >= 0) {
	      p = m0*mtr[k0].ndim1[k1] + m1;
	      if (mtr[k0].pma[k1][p]) {
		for (q = 0; q < mbpt_tr.naw; q++) {
		  ad->mk[iz][n+q] += mtr[k0].pma[k1][p]->tma[q];
		  if (i0 != i1) {
		    ad->mk[iz+ns][n+q] += mtr[k0].pma[k1][p]->rma[q];
		  }
		}
	      }
	    }
	    k0 = h[i1].pj*2*mbpt_tr.nktr + 2*ik;
	    if (IsEven(p0+p1+k)) k0++;
	    //k1 = IBisect(h[i0].pj, mtr[k0].nsym1, mtr[k0].isym1);
	    k1 = IdxGet(mtr[k0].ids, h[i0].pj);
	    if (k1 >= 0) {
	      p = m1*mtr[k0].ndim1[k1] + m0;
	      if (mtr[k0].pma[k1][p]) {
		for (q = 0; q < mbpt_tr.naw; q++) {
		  ad->mk[iz][n+q] += mtr[k0].pma[k1][p]->rma[q];
		  if (i0 != i1) {
		    ad->mk[iz+ns][n+q] += mtr[k0].pma[k1][p]->tma[q];
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
  int n, i, j, k, m, t, ti, q, m1, m2, p;
  int i0, i1, p1, p2, j1, j2, q1, q2;

  if (MyRankMPI() != 0) return;
  double wt0 = WallTime();
  fn = mbpt_tr.tfn;
  if (fn == NULL || mbpt_tr.nlow <= 0 || mbpt_tr.nup <= 0) return;  
  awgrid = mbpt_tr.awgrid;
  fhdr.type = DB_TR;
  strcpy(fhdr.symbol, GetAtomicSymbol());
  fhdr.atom = GetAtomicNumber();
  tr_hdr.nele = GetNumElectrons(0);
  n = GetNumLevels();
  f = OpenFile(fn, &fhdr);
  for (ti = 0; ti < mbpt_tr.nktr; ti++) {
    t = mbpt_tr.ks[ti];
    for (q = -1; q <= 1; q += 2) {      
      m = t*q;
      printf("multipole: %2d %10.3E %10.3E\n", m, WallTime()-wt0, TotalSize());
      tr_hdr.multipole = m;
      tr_hdr.gauge = GetTransitionGauge();
      tr_hdr.mode = 0;
      InitFile(f, &fhdr, &tr_hdr);
      ResetWidMPI();
#pragma omp parallel default(shared) private(i, j, k, sym, st, lev1, lev2, e, p1, j1, p2, j2, s0, i0, i1, m1, m2, a, p, rg, x, s, r, q1, q2)
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
	  i0 = lev1->pj*2*mbpt_tr.nktr + 2*ti;
	  if (q > 0) i0++;
	  i1 = IdxGet(mtr[i0].ids, lev2->pj);
	  //i1 = IBisect(lev2->pj, mtr[i0].nsym1, mtr[i0].isym1);
	  if (i1 >= 0) {
	    for (m1 = 0; m1 < lev1->n_basis; m1++) {
	      for (m2 = 0; m2 < lev2->n_basis; m2++) {
		a = lev1->mixing[m1] * lev2->mixing[m2];
		q1 = IdxGet(mtr[i0].idb[lev1->pj], lev1->basis[m1]);
		if (q1 < 0) continue;
		q2 = IdxGet(mtr[i0].idb[lev2->pj], lev2->basis[m2]);
		if (q2 < 0) continue;
		p = q1 * mtr[i0].ndim1[i1] + q2;
		if (mtr[i0].pma[i1][p] == NULL) continue;
		rg = mtr[i0].pma[i1][p]->tma;
		x = InterpolateMultipole(e, mbpt_tr.naw, awgrid, rg);
		if (tr_hdr.gauge == G_COULOMB && m < 0) {
		  x /= e;
		}
		s += a*x;
	      }
	    }
	  }
	  double sa = s;
	  i0 = lev2->pj*2*mbpt_tr.nktr + 2*ti;
	  if (q > 0) i0++;
	  i1 = IdxGet(mtr[i0].ids, lev1->pj);
	  //i1 = IBisect(lev1->pj, mtr[i0].nsym1, mtr[i0].isym1);
	  if (i1 >= 0) {
	    for (m1 = 0; m1 < lev1->n_basis; m1++) {
	      for (m2 = 0; m2 < lev2->n_basis; m2++) {
		a = lev1->mixing[m1] * lev2->mixing[m2];
		q1 = IdxGet(mtr[i0].idb[lev1->pj], lev1->basis[m1]);
		if (q1 < 0) continue;
		q2 = IdxGet(mtr[i0].idb[lev2->pj], lev2->basis[m2]);
		if (q2 < 0) continue;
		p = q2 * mtr[i0].ndim1[i1] + q1;
		if (mtr[i0].pma[i1][p] == NULL) continue;
		rg = mtr[i0].pma[i1][p]->rma;
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
	  a = s/s0;
	  a *= a;
	  if (a < mbpt_angzc) {
	    r.strength = s0+s;
	  } else {
	    r.strength = s0;	    
	  }
	  WriteTRRecord(f, &r, NULL);
	}
      }
      }
      if (mbpt_angzm > 0 && TotalSize() >= mbpt_angzm*1e9) {
	MPrintf(0, "CleanAngZArray: %2d %12.5E %12.5E\n",
		m, WallTime()-wt0, TotalSize());
	CleanAngZArray();
      }
      DeinitFile(f, &fhdr);
      fflush(stdout);
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
  double *x, *y, *z, *t, *zz1, *zz2, *ym;
  double **z1[MAX_SYMMETRIES], **z2[MAX_SYMMETRIES];
  double ***hab[MAX_SYMMETRIES], ***hba[MAX_SYMMETRIES];
  double ***nab[MAX_SYMMETRIES], ***nba[MAX_SYMMETRIES];
  double **nab1[MAX_SYMMETRIES], **nba1[MAX_SYMMETRIES];
  double *heff[MAX_SYMMETRIES], *neff[MAX_SYMMETRIES];
  double a, b, na, nb;
  HAMILTON *h;
  SYMMETRY *sym;
  MBPT_HAM *mbpt;
  MBPT_TR *mtr;
  char tfn1[1024];
  FILE **f1, *f2;
  IDXARY ing, ing0, *ing2;
  int mr, nr;
  int nkgp = 0;
  int *kgp = NULL;
  int rh = 0;
  IDXARY *idb[MAX_SYMMETRIES];

  mr = MPIRank(&nr);
  if (mr != 0) return 0;
  PrintMBPTOptions();
  double wt0 = WallTime();
  if (fn2 && strlen(fn2) > 0) {
    ReadHamilton(fn2, &nkg0, &nkg, &kg, &nkgp, &kgp, 1);
    if (nkgp > 0) {
      i = nkg + nkgp;
      kg = realloc(kg, sizeof(int)*i);
      memcpy(kg+nkg, kgp, sizeof(int)*nkgp);
      nkg = i;
      free(kgp);
      kgp = NULL;
    }
    rh = 1;
  }
  ierr = 0;
  if (nkg0 <= 0 || nkg0 > nkg) nkg0 = nkg;
  printf("mbpt configs: %d %d %d\n", nkg, nkg0, nkgp);
  mbpt = malloc(sizeof(MBPT_HAM)*nf);
  f1 = malloc(sizeof(FILE *)*nf);
  for (m = 0; m < nf; m++) {
    f1[m] = fopen(fn1[m], "r");
    if (f1[m] == NULL) {
      printf("cannot open file %s\n", fn1[m]);
      return -1;
    }
  }
  ierr = ReadMBPT(nf, f1, mbpt, 0);
  if (ierr < 0) {
    printf("cannot read MBPT Hamilton\n");
    return -1;
  }
  n = 0;
  n0 = 0;
  for (m = 0; m < nf; m++) {
    n0 += mbpt[m].n;
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
    for (k = 0; k < mbpt[m].n; k++) {
      ng0[j++] = mbpt[m].ng[k];
    }
  }
  if (i > 0) {
    n = SortUnique(i, ng);
  }
  if (j > 0) {
    n0 = SortUnique(j, ng0);  
  }
  n2m = 0;
  n2 = NULL;
  ng2 = NULL;
  ing2 = NULL;
  if (n > 0 && mbpt[0].n2 > 0) {
    n2 = malloc(sizeof(int)*n);
    ng2 = malloc(sizeof(int *)*n);
    ing2 = malloc(sizeof(IDXARY)*n);
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
      InitIdxAry(&ing2[i], n2[i], ng2[i]);	
    }
  }
  InitIdxAry(&ing, n, ng);
  InitIdxAry(&ing0, n0, ng0);
  n2m += n0;
  double wt1 = WallTime();
  printf("MBPT Structure ... %g %d %d %d\n", wt1-wt0, n0, n, n2m);
  fflush(stdout);
  nlevels = GetNumLevels();
  for (isym = 0; isym < MAX_SYMMETRIES; isym++) {
    wt0 = WallTime();
    //printf("rmbpt: %d\n", isym);
    ierr = ReadMBPT(nf, f1, mbpt, 1);
    if (mbpt[0].dim < 0) exit(1);
    idb[isym] = NULL;
    hab[isym] = NULL;
    hba[isym] = NULL;
    nab[isym] = NULL;
    nba[isym] = NULL;
    nab1[isym] = NULL;
    nba1[isym] = NULL;
    z1[isym] = NULL;
    z2[isym] = NULL;
    heff[isym] = NULL;
    neff[isym] = NULL;
    h = GetHamilton(isym);
    sym = GetSymmetry(isym);
    DecodePJ(isym, &pp, &jj);      
    k = 1;
    for (m = 0; m < nf; m++) {
      if (mbpt[m].dim == 0) {
	k = 0;
	break;
      }
    }
    if (k == 0) {
      for (j = 0; j < h->dim; j++) {
	for (i = 0; i <= j; i++) {	  
	  ierr = ReadMBPT(nf, f1, mbpt, 2);
	}
      }
      continue;
    }
    if (rh == 0) {
      AllocHamMem(h, mbpt[0].dim, mbpt[0].dim);
      h->pj = isym;
      h->n_basis = mbpt[0].dim;
      h->dim = mbpt[0].dim;
      h->hsize = h->dim*(h->dim+1)/2;
      memcpy(h->basis, mbpt[0].basis, sizeof(int)*h->n_basis);      
      k = ConstructHamilton(isym, 0, 0, NULL, 0, NULL, 1);
    }
    int *bs = malloc(sizeof(int)*h->dim);
    memcpy(bs, h->basis, sizeof(int)*h->dim);
    idb[isym] = malloc(sizeof(IDXARY));
    InitIdxAry(idb[isym], h->dim, bs);
    int dim2 = h->dim*h->dim;
    heff[isym] = malloc(sizeof(double)*dim2);
    neff[isym] = malloc(sizeof(double)*dim2);
    dim2 = h->dim*(h->dim+1)/2;
    z1[isym] = malloc(sizeof(double *)*dim2);
    z2[isym] = malloc(sizeof(double *)*dim2);
    nab1[isym] = malloc(sizeof(double *)*dim2);
    nba1[isym] = malloc(sizeof(double *)*dim2);
    for (i = 0; i < dim2; i++) {
      z1[isym][i] = NULL;
      z2[isym][i] = NULL;
      nab1[isym][i] = NULL;
      nba1[isym][i] = NULL;
    }
    if (n > 0 && mbpt[0].n2 > 0) {
      hab[isym] = malloc(sizeof(double **)*dim2);
      hba[isym] = malloc(sizeof(double **)*dim2);
      nab[isym] = malloc(sizeof(double **)*dim2);
      nba[isym] = malloc(sizeof(double **)*dim2);
      for (i = 0; i < dim2; i++) {
	hab[isym][i] = NULL;
	hba[isym][i] = NULL;
	nab[isym][i] = NULL;
	nba[isym][i] = NULL;
      }
    }
    for (j = 0; j < h->dim; j++) {
      for (i = 0; i <= j; i++) {
	k = j*(j+1)/2 + i;
	k0 = j*h->dim + i;
	k1 = i*h->dim + j;
	ierr = ReadMBPT(nf, f1, mbpt, 2);
	if (ierr < 0) {
	  MPrintf(-1, "Error reading MBPT 2: %d %d %d\n", isym, i, j);
	  Abort(1);
	}
	heff[isym][k0] = mbpt[0].a;
	heff[isym][k1] = mbpt[0].a;
	if (mbpt[0].heff0) {
	  heff[isym][k0] += mbpt[0].heff0[k0];
	  if (i != j) {
	    heff[isym][k1] += mbpt[0].heff0[k1];
	  }
	}
	neff[isym][k0] = 0.0;
	neff[isym][k1] = 0.0;
	for (m = 0; m < nf; m++) {
	  if (mbpt[m].ibra >= 0 && mbpt[m].iket >= 0) break;
	}	    
	if (m == nf) {
	  continue;
	}
	z1[isym][k] = malloc(sizeof(double)*n0);
	z2[isym][k] = malloc(sizeof(double)*n0);
	nab1[isym][k] = malloc(sizeof(double)*n0);
	nba1[isym][k] = malloc(sizeof(double)*n0);
	if (n > 0 && mbpt[0].n2 > 0) {	  
	  hab[isym][k] = malloc(sizeof(double *)*n);
	  hba[isym][k] = malloc(sizeof(double *)*n);
	  nab[isym][k] = malloc(sizeof(double *)*n);
	  nba[isym][k] = malloc(sizeof(double *)*n);
	  for (m = 0; m < n; m++) {
	    if (n2[m] > 0) {
	      hab[isym][k][m] = malloc(sizeof(double)*n2[m]);
	      hba[isym][k][m] = malloc(sizeof(double)*n2[m]);
	      nab[isym][k][m] = malloc(sizeof(double)*n2[m]);
	      nba[isym][k][m] = malloc(sizeof(double)*n2[m]);
	    }
	  }	
	}
	double **habi, **hbai, **nabi, **nbai;
	if (hab[isym] == NULL) {
	  habi = NULL;
	  hbai = NULL;
	  nabi = NULL;
	  nbai = NULL;
	} else {
	  habi = hab[isym][k];
	  hbai = hba[isym][k];
	  nabi = nab[isym][k];
	  nbai = nba[isym][k];
	}
	CombineMBPT(nf, mbpt, z1[isym][k], z2[isym][k],
		    habi, hbai,
		    nab1[isym][k], nba1[isym][k],
		    nabi, nbai, 
		    &ing0, &ing, ing2);
      }
    }
    wt1 = WallTime();
    MPrintf(0, "ReadMBPT: %d %d %d %d %d %g\n",
	    isym, h->dim, mbpt[0].dim, n, mbpt[0].n2, wt1-wt0);
    wt0 = wt1;
    ResetWidMPI();
#pragma omp parallel default(shared) private(i, j, k, k0, k1, x, t, y, a, b, na, nb, r, q)
    {
    x = malloc(sizeof(double)*n2m*3);  
    t = x + n2m;
    y = t + n2m;
    for (j = 0; j < h->dim; j++) {
      for (i = 0; i <= j; i++) {
	int skip = SkipMPI();
	if (skip) continue;
	k = j*(j+1)/2 + i;
	k0 = j*h->dim + i;
	k1 = i*h->dim + j;
	if (z1[isym][k] == NULL) continue;
	if (n0 > 0) {
	  for (q = 0; q < n0; q++) {
	    x[q] = ng0[q];
	  }
	  a = SumInterp1D(n0, z1[isym][k], x, t, y, 1);
	  heff[isym][k0] += a;
	  na = SumInterp1D(n0, nab1[isym][k], x, t, y, 2);
	  neff[isym][k0] += na;
	  b = a;
	  nb = na;
	  if (i < j) {
	    b = SumInterp1D(n0, z2[isym][k], x, t, y, 3);
	    heff[isym][k1] += b;
	    nb = SumInterp1D(n0, nba1[isym][k], x, t, y, 4);
	    neff[isym][k1] += nb;
	  }
	}
	if (n > 0 && hab[isym]) {
	  for (r = 0; r < n; r++) {
	    if (n2[r] > 0) {
	      for (q = 0; q < n2[r]; q++) {
		x[q] = ng[r] + ng2[r][q];
	      }
	      z1[isym][k][r] = SumInterp1D(n2[r], hab[isym][k][r],
					   x, t, y, ng[r]*10+1);
	      if (i < j) z2[isym][k][r] = SumInterp1D(n2[r], hba[isym][k][r],
						      x, t, y, ng[r]*10+2);
	      else z2[isym][k][r] = z1[isym][k][r];		
	      nab1[isym][k][r] = SumInterp1D(n2[r], nab[isym][k][r],
					     x, t, y, ng[r]*10+3);
	      if (i < j) nba1[isym][k][r] = SumInterp1D(n2[r], nba[isym][k][r],
							x, t, y, ng[r]*10+4);
	    }
	  }	  
	  for (q = 0; q < n; q++) {
	    x[q] = ng[q];
	  }
	  a = SumInterp1D(n, z1[isym][k], x, t, y, -1);
	  heff[isym][k0] += a;
	  na = SumInterp1D(n, nab1[isym][k], x, t, y, -2);
	  neff[isym][k0] += na;
	  b = a;
	  nb = na;
	  if (i < j) {
	    b = SumInterp1D(n, z2[isym][k], x, t, y, -3);
	    heff[isym][k1] += b;
	    nb = SumInterp1D(n, nba1[isym][k], x, t, y, -4);
	    neff[isym][k1] += nb;
	  }
	}
      OUT:
	continue;
      }
    }
    free(x);
    }
    for (k = 0; k < dim2; k++) {
      if (z1[isym][k]) {
	free(z1[isym][k]);
	free(z2[isym][k]);
	free(nab1[isym][k]);
	free(nba1[isym][k]);
	if (hab[isym]) {
	  for (m = 0; m < n; m++) {
	    if (n2[m] > 0) {
	      free(hab[isym][k][m]);
	      free(hba[isym][k][m]);
	      free(nab[isym][k][m]);
	      free(nba[isym][k][m]);
	    }
	  }
	  free(hab[isym][k]);
	  free(hba[isym][k]);
	  free(nab[isym][k]);
	  free(nba[isym][k]);
	}
      }
    }
    free(z1[isym]);
    free(z2[isym]);
    free(nab1[isym]);
    free(nba1[isym]);
    if (hab[isym]) {
      free(hab[isym]);
      free(hba[isym]);
      free(nab[isym]);
      free(nba[isym]);
    }    
    for (m = 0; m < nf; m++) {      
      free(mbpt[m].basis);
      mbpt[m].basis = NULL;
      if (mbpt[m].heff0) {
	free(mbpt[m].heff0);
	mbpt[m].heff0 = NULL;
	mbpt[m].hsize0 = 0;
      }
    }
    wt1 = WallTime();
    printf("SumInterp: %d %d %d %g\n", isym, h->dim, mbpt[0].dim, wt1-wt0);
  }
  
  ResetWidMPI();
#pragma omp parallel default(shared) private(isym, h, pp, jj, wt0, wt1, i, j, k, k0, a, b, ym, m)
  {
  for (isym = 0; isym < MAX_SYMMETRIES; isym++) {
    if (heff[isym] == NULL) continue;
    wt0 = WallTime();
    int skip = SkipMPI();
    if (skip) continue;
    h = GetHamilton(isym);
    DecodePJ(isym, &pp, &jj);
    h->heff = heff[isym];
    if (DiagnolizeHamilton(h) < 0) {
      printf("Diagnolizing Hamiltonian Error\n");
      ierr = -1;
      break;
    }
    
    /* correct the normalization */
    ym = h->mixing + h->dim;
    for (i = 0; i < h->dim; i++) {
      a = 0.0;
      for (j = 0; j < h->dim; j++) {
	for (k = 0; k < h->dim; k++) {
	  k0 = j*h->dim + k;
	  a += neff[isym][k0] * ym[j] * ym[k];
	}
      }
      if (a > 0 && a < 1) {
	b = sqrt(1.0 + a);
	for (j = 0; j < h->n_basis; j++) {
	  ym[j] /= b;
	}
      }
      ym += h->n_basis;
    }
    
    AddToLevels(h, nkg0, kg);
    wt1 = WallTime();
    double ts = TotalSize();
    MPrintf(-1, "Diag Ham: %3d %5d %5d %10.3E %10.3E\n",
	    isym, h->dim, h->n_basis, wt1-wt0, ts);
    fflush(stdout);    
    h->heff = NULL;
    free(heff[isym]);
    free(neff[isym]);
  }
  }
  SortLevels(nlevels, -1, 0);
  SaveLevels(fn, nlevels, -1);
  for (m = 0; m < nf; m++) {
    fclose(f1[m]);    
  }
  if (mbpt_tr.nktr > 0) {
    InitTransitionMBPT(&mtr, idb);
  }
  for (isym = 0; isym < MAX_SYMMETRIES; isym++) {
    AllocHamMem(h, -1, -1);
    AllocHamMem(h, 0, 0);
  }
  printf("End MBPT Structure: %g %g\n", WallTime()-wt0, TotalSize());
  if (mbpt_tr.nktr == 0) goto ERROR;
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
  printf("MBPT Transition ... %g\n", WallTime()-wt0);
  k = 2*mbpt_tr.nktr*MAX_SYMMETRIES;
  for (j = 0; j < k; j++) {
    ierr = ReadMBPT(nf, f1, mbpt, -(j+1));
    if (mbpt[0].mtr[j].nsym1 == 0) {
      continue;
    }
    CombineTransitionMBPT(nf, mbpt, mtr, j, &ing0);
  }
  printf("InterpTransitionMBPT ... %g\n", WallTime()-wt0);  
  ResetWidMPI();
  int m0, m1, p0, q0, q1;
#pragma omp parallel default(shared) private(q, i, j, m0, m1, p0, q0, q1)
  {
    double *x = malloc(sizeof(double)*n2m*5);
    double *t = x + n2m;
    double *y = t + n2m;
    double *zz1 = y + n2m;
    double *zz2 = zz1 + n2m;
    for (q = 0; q < n0; q++) {
      x[q] = ng0[q];
    }
    for (j = 0; j < k; j++) {
      double iwt0 = WallTime();
      for (i = 0; i < mtr[j].nsym1; i++) {
	for (m0 = 0; m0 < mtr[j].ndim0; m0++) {
	  for (m1 = 0; m1 < mtr[j].ndim1[i]; m1++) {
	    int skip = SkipMPI();
	    if (skip) continue;
	    p0 = m0*mtr[j].ndim1[i] + m1;
	    if (mtr[j].pma[i][p0] == NULL) continue;
	    for (q0 = 0; q0 < mbpt_tr.naw; q0++) {
	      for (q = 0; q < n0; q++) {
		q1 = q0*n0 + q;
		zz1[q] = mtr[j].pma[i][p0]->tma[q1];
		zz2[q] = mtr[j].pma[i][p0]->rma[q1];
	      }
	      mtr[j].pma[i][p0]->tma[q0] = SumInterp1D(n0, zz1, x, t, y, 5);
	      mtr[j].pma[i][p0]->rma[q0] = SumInterp1D(n0, zz2, x, t, y, 6);
	    }
	  }
	}
      }
      double iwt1 = WallTime();
      if (mtr[j].nsym1 > 0) {
	MPrintf(0, "SumInterpTR: %4d %4d %4d %4d %12.5E\n",
		j, mtr[j].isym0, mtr[j].m, mtr[j].nsym1, iwt1-iwt0);
      }
    }
    free(x);
  }
  SetTransitionMode(0);
  printf("SaveTransitionMBPT ... %12.5E %12.5E\n", WallTime()-wt0, TotalSize());
  fflush(stdout);
  SaveTransitionMBPT(mtr);
  printf("AdjustAngZ ... %12.5E %12.5E\n", WallTime()-wt0, TotalSize());
  fflush(stdout);
  AdjustAngularZ(mtr);
  printf("End TransitionMBPT ... %12.5E %12.5E\n", WallTime()-wt0, TotalSize());
  fflush(stdout);

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
    free(mbpt[m].hab1);
    free(mbpt[m].hba1);
    if (mbpt[m].n2 > 0) {
      free(mbpt[m].ng2);
      free(mbpt[m].hab);
      free(mbpt[m].hba);
    }
  }
  free(f1);
  free(mbpt);
  if (n > 0) {
    free(ng);
    FreeIdxAry(&ing, 2);
    if (n2) {
      for (i = 0; i < n; i++) {
	if (n2[i] > 0) {
	  free(ng2[i]);
	  FreeIdxAry(&ing2[i], 2);
	}
      }
      free(ng2);      
      free(n2);
      free(ing2);
    }
  }
  if (n0 > 0) {
    free(ng0);  
    FreeIdxAry(&ing0, 2);
  }
  if (rh && nkg > 0) {
    free(kg);
  }
  for (isym = 0; isym < MAX_SYMMETRIES; isym++) {
    if (idb[isym]) {
      FreeIdxAry(idb[isym], 0);
    }
  }
  return ierr;
}
