
#include "mbpt.h"
#include "cf77.h"

static char *rcsid="$Id: mbpt.c,v 1.1 2005/07/18 15:39:43 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

static int mbpt_extra = 0;
static int mbpt_isym = -1;
static int mbpt_ilev = -1;
static double mbpt_mcut = EPS3;
static int mbpt_n3 = 0;

void SetOptMBPT(int n3, double c) {
  mbpt_n3 = 0;
  mbpt_mcut = c;
}

void SetExtraMBPT(int m) {
  mbpt_extra = m;
}

void SetSymMBPT(int p, int j, int i) {
  if (p < 0 || j < 0) mbpt_isym = -1;
  else mbpt_isym = IsEven(p)?(2*j):(2*j+1);
  mbpt_ilev = i;
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

void BaseConfig(int n, int *kg, int *nbc, CONFIG **bc, 
		int *nbc1, CONFIG **bc1, int *nb, int **bk, int sr) {
  int nmax, nsm, ncs, i, j, km, k, t, m, p, jp, nq;
  int mcs, m1e, m2e, mcs1, *ik, em[2], b1e, b2e;
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

  mcs = ncs*nsm*(1+nsm);
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
	em[1] = ShellToInt(c0[i]->shells[jp].n, c0[i]->shells[jp].kappa);
	ConfigChangeNE(&m, c, k, ik, -2, em);
      }
    }
  }
  m2e = m;
  
  mcs1 = m1e*nsm + m2e*nsm*nsm;
  *bc1 = malloc(sizeof(CONFIG)*mcs1);
  c1 = *bc1;
  m = 0;
  c = *bc;
  for (i = 0; i < m1e; i++) {
    for (j = 0; j < k; j++) {
      if (ik[j] >= 0) ik[j] = 0;
    }
    for (j = 0; j < c[i].n_shells; j++) {
      t = ShellToInt(c[i].shells[j].n, c[i].shells[j].kappa);
      ik[t] = c[i].shells[j].nq;
    }
    for (j = 0; j < c[i].n_shells; j++) {
      em[0] = ShellToInt(c[i].shells[j].n, c[i].shells[j].kappa);
      ConfigChangeNE(&m, c1, k, ik, 1, em);
    }
  }
  c = c + m1e;
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
	em[1] = ShellToInt(c[i].shells[jp].n, c[i].shells[jp].kappa);
	ConfigChangeNE(&m, c1, k, ik, 2, em);
      }
    }
  }
  *nbc = m1e + m2e;  
  c = *bc;
  for (i = 0; i < *nbc; i++) {
    RemoveEmpty(c+i);
    if (i < m1e) c[i].n_electrons = nq-1;
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

  if (sr == 0) {
    for (i = 0; i < ncs; i++) {
      ConstructConfigName(sc, 2000, c0[i]);
      printf("0E: %4d   %s\n", i, sc);
    }
    c = *bc;
    for (i = 0; i < m1e; i++) {
      ConstructConfigName(sc, 2000, c+i);
      printf("1E: %4d   %s\n", i, sc);
    }
    c += m1e;
    for (i = 0; i < m2e; i++) {
      ConstructConfigName(sc, 2000, c+i);
      printf("2E: %4d   %s\n", i, sc);
    }
    c = *bc1;
    for (i = 0; i < *nbc1; i++) {
      ConstructConfigName(sc, 2000, c+i);
      printf("EX: %4d   %s\n", i, sc);
    }
    fflush(stdout);
  }
  free(c0);
}

int ConstructNGrid(int n, int **g) {
  int i, *p, *q, n0, n1, nt, m, di, dm;

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
  
int StructureMBPT0(char *fn, char *fn1, int n, int *s0, int kmax,
		   int n1, int *nm, int n2, int *nmp, char *gn0) {
  CONFIG_GROUP *g1, *g0;
  CONFIG *c1, *bc, *bc1;
  SYMMETRY *sym;
  STATE *st;
  LEVEL *lev;
  HAMILTON *ha;
  char gn[GROUP_NAME_LEN] = "_@nb@_";
  int nbc, nbc1, k;
  int nele0, nele1, nk;
  int i, p, np, nq, inp, inq, k0, k1, q;
  int jp, jq, kap, kaq, kp, kq, kp2, kq2;
  int ic, kgp, nb, *bk;
  int ncc, ncc0, ncc1, ncc2;
  int *bs1, nbs1, *bs0, nbs0, *s;
  double a1, a2, d, *e1;
  CORR_CONFIG cc, *ccp, *ccp1;
  ARRAY ccfg;
  int *icg, ncg, icg0, icg1, rg, np0, ig0;
  double a, b, xnq, ynq, *tq, tnq;
  double **de, **deq;
  double dnq[15];
  double *ham, **dh1, **dh2, *mix;
  int m, nlev, nlev1, ilev, t, kb, r, dim;  
  typedef struct _MBPT_BASE_ {
    int isym;
    int nbasis, nb2, bmax;
    int *basis;
    double *ene;
    double **dh1, **dh2;
    double **dy1, **dy2;
  } MBPT_BASE;
  MBPT_BASE mb, *mbp;
  ARRAY base;
  char fn2[256];
  FILE *f, *fp;
  double t0, t1, t2;
  int sr, nr;
#ifdef USE_MPI
  int *ics0, *ics1;
  double *de0, *de1;
  MPI_Comm_rank(MPI_COMM_WORLD, &sr);
  MPI_Comm_size(MPI_COMM_WORLD, &nr);
  printf("RANK: %2d, TOTAL: %2d\n", sr, nr);
#else
  sr = 0;
  nr = 1;
#endif
  
  t0 = clock();
  t0 /= CLOCKS_PER_SEC;

  if (gn0 == NULL) gn0 = gn;
  n1 = ConstructNGrid(n1, &nm);
  n2 = ConstructNGrid(n2, &nmp);

  ha = GetHamilton();
  BaseConfig(n, s0, &nbc, &bc, &nbc1, &bc1, &nb, &bk, sr);

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
	  nq = np + nmp[inq];
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
  ArrayInit(&base, sizeof(MBPT_BASE), 256);
  k = nbc;
  icg = malloc(sizeof(int)*(n1*n2*k+2));

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
    ccp = ArrayAppend(&ccfg, &cc, NULL);
  } 
  icg[t++] = ccfg.dim;      
  for (p = 0; p < nbc; p++) {
    c1 = bc + p;
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
	    ccp = ArrayAppend(&ccfg, &cc, NULL);
	  }
	}
	icg[t++] = ccfg.dim;
      }
    } else if (nele1 == nele0-2) {
      for (inp = 0; inp < n1; inp++) {
	np = nm[inp];
	for (inq = 0; inq < n2; inq++) {
	  nq = np + nmp[inq];
	  for (kp = 0; kp <= kmax; kp++) {
	    if (kp >= np) break;
	    kp2 = 2*kp;
	    for (jp = kp2-1; jp <= kp2+1; jp += 2) {
	      if (jp < 0) continue;
	      kap = GetKappaFromJL(jp, kp2);
	      k0 = OrbitalIndex(np, kap, 0);
	      k1 = IBisect(k0, nb, bk);
	      if (k1 >= 0) continue;
	      for (kq = 0; kq <= kmax; kq++) {
		if (kq >= nq) break;
		kq2 = 2*kq;
		for (jq = kq2-1; jq <= kq2+1; jq += 2) {
		  if (jq < 0) continue;
		  kaq = GetKappaFromJL(jq, kq2);
		  if (np == nq && kaq > kap) continue;
		  k0 = OrbitalIndex(nq, kaq, 0);
		  k1 = IBisect(k0, nb, bk);
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
		  ccp = ArrayAppend(&ccfg, &cc, NULL);
		}
	      }
	    }
	  }
	  icg[t++] = ccfg.dim;
	}
      }
    }
  }

  ncg = t-1;
  r = ncg/nr;
  icg0 = sr*r;
  if (sr < nr-1) {
    icg1 = (sr+1)*r;
  } else {
    icg1 = ncg;
  }
  printf("RANK: %2d, %d %d %d %d\n", sr, ncg, icg0, icg1, n1*n2*k+1);

  for (i = 0; i < MAX_SYMMETRIES; i++) {
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

  nlev = GetNumLevels();
  s = s0;
  for (i = 0; i < base.dim; i++) {    
    mbp = ArrayGet(&base, i);
    nk = ConstructHamilton(mbp->isym, n, n, s, 0, NULL);
    printf("RANK: %2d, dim: %3d %3d %6d\n", sr, i, mbp->isym, ha->dim);
    fflush(stdout);
    DiagnolizeHamilton();
    mix = ha->mixing + ha->dim;
    q = ha->n_basis;
    nbs0 = q*q;	
    mbp->dh1 = malloc(sizeof(double *)*nbs0);
    mbp->dh2 = malloc(sizeof(double *)*nbs0);
    for (t = 0; t < nbs0; t++) {
      mbp->dh1[t] = malloc(sizeof(double)*n1);
      mbp->dh2[t] = malloc(sizeof(double)*n1*n2);
      for (r = 0; r < n1; r++) {
	mbp->dh1[t][r] = 0.0;
      } 
      for (r = 0; r < n1*n2; r++) {
	mbp->dh2[t][r] = 0.0;
      }
    }
    AddToLevels(n, s);
    mbp->nb2 = nbs0;
  }

  if (sr == 0) {
    SaveLevels(fn, nlev, -1);
    sprintf(fn2, "%s.basis", fn1);
    GetBasisTable(fn2);
  }

  t1 = clock();
  t1 /= CLOCKS_PER_SEC;
  ncc2 = 0;
  for (rg = icg0; rg < icg1; rg++) {
    printf("RANK: %2d, %d %d\n", sr, rg, ncg);
    fflush(stdout);
    ncc0 = icg[rg];
    ncc = 0;
    for (p = icg[rg]; p < icg[rg+1]; p++) {
      ccp = ArrayGet(&ccfg, p);
      kgp = GroupIndex(gn0);
      AddMBPTConfig(kgp, ccp);
      ncc++;
      if (ncc == ccfg.block || p == icg[rg+1]-1) {
	t2 = clock();
	t2 /= CLOCKS_PER_SEC;
	printf("RANK: %2d, %5d %7d %7d %7d %10.3E %10.3E\n", sr, ncc, p+1, 
	       icg[rg+1], ccfg.dim, t2-t1, t2-t0);
	fflush(stdout);
	t1 = t2;
	for (i = 0; i < base.dim; i++) {
	  mbp = ArrayGet(&base, i);
	  nk = ConstructHamiltonDiagonal(mbp->isym, 1, &kgp, 0);
	  if (nk < 0) continue;
	  nbs1 = ha->dim;
	  bs1 = ha->basis;
	  e1 = ha->hamilton;	
	  printf("RANK: %2d, sym: %3d %3d %3d %6d\n", 
		 sr, i, mbp->isym, mbp->nbasis, nbs1);
	  fflush(stdout);
	  sym = GetSymmetry(mbp->isym);
	  nbs0 = mbp->nbasis;
	  bs0 = mbp->basis;
	  ham = malloc(sizeof(double)*nbs0*nbs1);
	  m = 0;
	  for (q = 0; q < nbs1; q++) {
	    st = ArrayGet(&(sym->states), bs1[q]);
	    r = st->kcfg + ncc0;
	    ccp1 = ArrayGet(&ccfg, r);
	    k1 = bs1[q];
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
	    dh1 = mbp->dh1;
	    dh2 = mbp->dh2;
	    for (r = 0; r < nbs0; r++) {
	      m = q*nbs0 + r;
	      a1 = ham[m];
	      if (1.0+a1 == 1.0) continue;
	      for (t = 0; t <= r; t++) {
		m = q*nbs0 + t;
		a2 = ham[m];
		if (1.0+a2 == 1.0) continue;
		a = a1*a2;
		d = a/(mbp->ene[t] - e1[q]);
		ic = r*nbs0 + t;
		if (ccp1->kq == 0) {
		  dh1[ic][ccp1->inp] += d;
		} else {
		  dh2[ic][(ccp1->inq)*n1 + ccp1->inp] += d;
		}
		if (r != t) {
		  d = a/(mbp->ene[r] - e1[q]);
		  ic = t*nbs0 + r;
		  if (ccp1->kq == 0) {
		    dh1[ic][ccp1->inp] += d;
		  } else {
		    dh2[ic][(ccp1->inq)*n1 + ccp1->inp] += d;
		  }
		}
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
    ncc2 += icg[rg+1] - icg[rg];
    if (ncc2 > 10*ccfg.dim) {
      ReinitRadial(1);
      ncc2 = 0;
    }
  }

#ifdef USE_MPI
  for (i = 0; i < base.dim; i++) {
    mbp = ArrayGet(&base, i);
    m = mbp->nb2*(n1 + n1*n2);
    de1 = malloc(sizeof(double)*m);
    de0 = malloc(sizeof(double)*m);
    p = 0;
    for (t = 0; t < mbp.nb2; t++) {
      for (inp = 0; inp < n1; inp++) {
	de1[p++] = mbp->dh1[t][inp];
	for (inq = 0; inq < n2; inq++) {
	  de1[p++] = mbp->dh2[t][inq*n1 + inp];
	}
      }
    }
    MPI_Reduce(de1, de0, m, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (sr == 0) {
      p = 0;
      for (t = 0; t < mbp.nb2; t++) {
	for (inp = 0; inp < n1; inp++) {
	  mbp->dh1[t][inp] = de0[p++];
	  for (inq = 0; inq < n2; inq++) {
	    mbp->dh2[t][inq*n1 + inp] = de0[p++];
	  }
	}
      }
    }
    free(de1);
    free(de0);
  }
#endif
     
  if (sr == 0) {
    nlev1 = GetNumLevels() - nlev;
    de = malloc(sizeof(double *)*nlev1);
    deq = malloc(sizeof(double *)*nlev1);
    for (i = 0; i < nlev1; i++) {
      de[i] = malloc(sizeof(double)*n1);
      deq[i] = malloc(sizeof(double)*n1*n2);
      for (t = 0; t < n1; t++) {
	de[i][t] = 0;
      }
      for (t = 0; t < n1*n2; t++) {
	deq[i][t] = 0;
      }
    }
    for (i = 0; i < nlev1; i++) {    
      ilev = i + nlev;
      lev = GetLevel(ilev);
      for (p = 0; p < base.dim; p++) {
	mbp = ArrayGet(&base, p);
	if (mbp->isym == lev->pj) break;
      }
      for (p = 0; p < lev->n_basis; p++) {
	r = IBisect(lev->basis[p], mbp->nbasis, mbp->basis);
	if (r < 0) continue;
	a1 = lev->mixing[p];
	for (q = 0; q < lev->n_basis; q++) {
	  t = IBisect(lev->basis[q], mbp->nbasis, mbp->basis);
	  if (t < 0) continue;
	  a2 = lev->mixing[q];
	  a = a1*a2;
	  m = r*mbp->nbasis + t;
	  for (inp = 0; inp < n1; inp++) {
	    de[ilev][inp] += a*(mbp->dh1[m][inp]);
	  }
	  for (inp = 0; inp < n1*n2; inp++) {
	    deq[ilev][inp] += a*(mbp->dh2[m][inp]);
	  }
	}
      }
    }

    f = fopen(fn1, "w");
    sprintf(fn2, "%s.detail", fn1);
    fp = fopen(fn2, "w");

    for (i = 0; i < nlev1; i++) {
      ilev = nlev + i;
      for (t = 0; t < n1; t++) {
	de[i][t] *= HARTREE_EV;
      }
      for (t = 0; t < n1*n2; t++) {
	deq[i][t] *= HARTREE_EV;
      }

      b = 0.0;
      a1 = 0.0;
      for (inp = 0; inp < n1; inp++) {	
	np = nm[inp];
	tq = deq[i] + inp*n2;
	for (inq = 0; inq < n2; inq++) {
	  nq = np + nmp[inq];
	  dnq[inq] = nq;
	  fprintf(fp, "%3d %3d %3d %12.5E %12.5E\n", 
		  ilev, np, nq, tq[inq], de[i][inp]);
	}
	for (inq = 0; inq < n2; inq++) {
	  if (tq[inq] < 0) break;
	}
	t = inq;
	for (; inq < n2; inq++) {
	  if (tq[inq] >= 0) break;
	}
	tnq = 0.0;
	if (inq == n2 && n2-t > 2) {
	  for (inq = 0; inq < t; inq++) {
	    tnq += tq[inq];
	  }
	  for (inq = t; inq < n2; inq++) {
	    tq[inq] = log(-tq[inq]);
	    dnq[inq] = log(dnq[inq]);
	  }
	  for (nq = np+nmp[t]; nq <= np + nmp[n2-1]; nq++) {
	    xnq = log(nq);
	    UVIP3P(3, n2-t, dnq+t, tq+t, 1, &xnq, &ynq);
	    tnq += -exp(ynq);
	  }
	  a = -(tq[n2-1] - tq[n2-2])/(dnq[n2-1]-dnq[n2-2]);
	  if (a > 1.0) {
	    tnq += -exp(tq[n2-1])*(np+nmp[n2-1])/(a-1.0);
	  }
	} else {
	  tnq = tq[0];
	  for (nq = np+1; nq <= np + nmp[n2-1]; nq++) {
	    xnq = nq;
	    UVIP3P(3, n2-1, dnq+1, tq+1, 1, &xnq, &ynq);
	    tnq += ynq;
	  }
	}
	de[i][inp] += tnq;
	a = de[i][inp];
	if (a) {
	  b += a;
	  fprintf(f, "%4d %5d %12.5E\n", np, ilev, a);
	}
      }
      fprintf(f, "\n#SUM %5d %12.5E\n\n", ilev, b);
      fprintf(fp, "\n");
    }
    fclose(f);
    fclose(fp);

    for (i = 0; i < nlev1; i++) {    
      free(de[i]);
      free(deq[i]);
    } 
    free(de);
    free(deq);
  }

 DONE:
  free(bc);
  free(bc1);
  free(bk);
  free(icg);
  ArrayFree(&ccfg, NULL);
  for (i = 0; i < base.dim; i++) {
    mbp = ArrayGet(&base, i);
    free(mbp->basis);
    free(mbp->ene);
    for (t = 0; t < mbp->nb2; t++) {
      free(mbp->dh1[t]);
      free(mbp->dh2[t]);
    }
    free(mbp->dh1);
    free(mbp->dh2);
  }
  ArrayFree(&base, NULL);

  t2 = clock();
  t2 /= CLOCKS_PER_SEC;
  printf("RANK: %2d, Total Time: %10.3E\n", sr, t2-t0);
  return 0;
}

int RadialBasisMBPT(int kmax, int n, int *ng, int **bas) {
  int nb, k, j, k2, i, m, ka;
  
  nb = n*(kmax+1)*2;
  (*bas) = malloc(sizeof(int)*nb);
  m = 0;
  for (k = 0; k <= kmax; k++) {
    k2 = 2*k;
    for (j = k2-1; j <= k2+1; j += 2) {
      if (j < 0) continue;
      ka = GetKappaFromJL(j, k2);
      for (i = 0; i < n; i++) {
	if (ng[i] <= k) continue;
	(*bas)[m] = OrbitalIndex(ng[i], ka, 0);
	/*
	printf("%2d %2d %2d %2d\n", (*bas)[m], ng[i], k, j);
	*/
	m++;
      }
    }
  }
  nb = m;
  (*bas) = ReallocNew(*bas, nb*sizeof(int));

  return nb;
}

int PadStates(SYMMETRY *sym, int si, int sj, SHELL **bra, SHELL **ket,
	      SHELL_STATE **sbra, SHELL_STATE **sket) {
  int i, j, m, k, ns;
  STATE *st1, *st2;
  SHELL *bra1, *ket1;
  SHELL_STATE *sbra1, *sket1, *s1, *s2;
  CONFIG *c1, *c2;

  st1 = ArrayGet(&(sym->states), si);
  st2 = ArrayGet(&(sym->states), sj);
  c1 = GetConfig(st1);
  c2 = GetConfig(st2);
  s1 = c1->csfs + st1->kstate;
  s2 = c2->csfs + st2->kstate;

  ns = c1->n_shells + c2->n_shells + 2;
  *bra = malloc(sizeof(SHELL)*ns);
  *ket = malloc(sizeof(SHELL)*ns);
  *sbra = malloc(sizeof(SHELL_STATE)*ns);
  *sket = malloc(sizeof(SHELL_STATE)*ns);
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
  if (!r) return 0.0;

  for (k0 = 1; k0 < n; k0++) {
    if (x[k0]-x[k0-1] > 1) break;
  }
  if (k0 == n) goto END;
  k0--;
  
  for (k1 = n-1; k1 >= k0; k1--) {
    if (z[n-1] > 0) {
      if (z[k1] <= 0) {
	k1++;
	break;
      }
    } else {
      if (z[k1] >= 0) {
	k1++;
	break;
      }
    }
  }
  if (k1 < 0) k1 = 0;

  for (i = Min(k0, k1); i < n; i++) {
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
  if (mbpt_extra) {
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
      if (c <= 0 && 2.0*log(x[i])*c+b < -1.1) {
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

double H22Term(int ns, SHELL_STATE *sbra, SHELL_STATE *sket,
	       INTERACT_SHELL *s, int *ks1, int *ks2,
	       FORMULA *fm, double *a, int md) {
  int m, kk1, kk2, kmin1, kmin2, kmax1, kmax2;
  int mkk1, mkk2, mkk;
  double c, y, sd1, sd2, se1, se2;
  double a1[MKK], a2[MKK];

  if (md <= 0) {
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
  if (kmax1 < kmin1) return 0.0;
  if (kmax2 < kmin2) return 0.0;
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
	EvaluateTensor(ns, sbra, sket, s, 1, fm);
	a[mkk+mkk2] = fm->coeff;
      }
    }
  }
  /* 
  ** these conditions must be placed after the above block,
  ** so that the recouping coeff. are calculated the first
  ** time even these parity conditions are not satisfied.
  */
  if (IsOdd((s[0].kl+s[1].kl+s[2].kl+s[3].kl)/2)) return 0.0;
  if (IsOdd((s[4].kl+s[5].kl+s[6].kl+s[7].kl)/2)) return 0.0;
  for (kk1 = kmin1; kk1 <= kmax1; kk1 += 2) {
    m = 0;
    mkk1 = kk1/2;
    mkk = mkk1*MKK;
    for (kk2 = kmin2; kk2 <= kmax2; kk2 += 2) {
      mkk2 = kk2/2;
      if (a[mkk+mkk2]) {
	m = 1;
	break;
      }
    }
    m = 1;
    if (m) {
      SlaterTotal(&sd1, &se1, NULL, ks1, kk1, 0);
      a1[mkk1] = sd1 + se1;
    }
  }
  for (kk2 = kmin2; kk2 <= kmax2; kk2 += 2) {
    m = 0;
    mkk2 = kk2/2;
    for (kk1 = kmin1; kk1 <= kmin1; kk1 += 2) {
      mkk1 = kk1/2;
      mkk = mkk1*MKK+mkk2;
      if (a[mkk]) {
	m = 1;
	break;
      }
    }
    m = 1;
    if (m) {
      SlaterTotal(&sd2, &se2, NULL, ks2, kk2, 0);
      a2[mkk2] = sd2 + se2;
    }
  }
  c = 0.0;
  for (kk1 = kmin1; kk1 <= kmax1; kk1 += 2) {
    mkk1 = kk1/2;
    mkk = mkk1*MKK;
    for (kk2 = kmin2; kk2 <= kmax2; kk2 += 2) {
      mkk2 = kk2/2;
      y = a[mkk+mkk2];
      if (fabs(y) < EPS10) continue;
      y /= sqrt((kk1+1.0)*(kk2+1.0));
      if (IsOdd((kk1+kk2)/2)) y = -y;      
      c += y*a1[mkk1]*a2[mkk2];
    }
  }
  return c;
}

double H12Term(int ns, SHELL_STATE *sbra, SHELL_STATE *sket,
	       INTERACT_SHELL *s, int *ks, int k0, int k1,
	       FORMULA *fm, double *a, int md) {  
  int kk, kmin, kmax;
  double c, c1, y, sd, se;

  if (md <= 0) {
    TriadsZ(1, 2, fm);
    RecoupleTensor(6, s, fm);
  }
  if (s[0].j != s[1].j) return 0.0;
  kmin = abs(s[2].j-s[3].j);
  kk = abs(s[4].j-s[5].j);
  kmin = Max(kmin, kk);
  kmax = s[2].j + s[3].j;
  kk = s[4].j + s[5].j;
  kmax = Min(kmax, kk);
  kk = 2*(MKK-1);
  kmax = Min(kmax, kk);
  if (kmax < kmin) return 0.0;
  if (md <= 1) {
    FixJsZ(s, fm);
    for (kk = kmin; kk <= kmax; kk += 2) {
      fm->js[7] = 0;
      fm->js[8] = kk;
      fm->js[9] = kk;
      fm->js[10] = 0;
      fm->js[11] = 0;
      EvaluateTensor(ns, sbra, sket, s, 1, fm);
      a[kk/2] = fm->coeff;
    }
  }

  /* 
  ** these conditions must be placed after the above block,
  ** so that the recouping coeff. are calculated the first
  ** time even these parity conditions are not satisfied.
  */
  if (IsOdd((s[0].kl+s[1].kl)/2)) return 0.0;
  if (IsOdd((s[2].kl+s[3].kl+s[4].kl+s[5].kl)/2)) return 0.0;
  c = 0.0;
  for (kk = kmin; kk <= kmax; kk += 2) {
    y = a[kk/2];
    if (fabs(y) < EPS10) continue;
    y /= sqrt(kk+1.0);
    if (IsOdd(kk/2)) y = -y;
    SlaterTotal(&sd, &se, NULL, ks, kk, 0);
    c += y*(sd + se);
  }
  if (fabs(c) < EPS10) return 0.0;
  ResidualPotential(&c1, k0, k1);
  c1 += QED1E(k0, k1);
  c1 *= sqrt(s[0].j+1.0);
  c *= c1;

  /* minus sign is from the definition of Z^k */
  return -c;
}

double H11Term(int ns, SHELL_STATE *sbra, SHELL_STATE *sket, 
	       INTERACT_SHELL *s, int k0, int k1, int k2, int k3,
	       FORMULA *fm, int md) {
  double y, c, c1, c2;

  if (md <= 0) {
    TriadsZ(1, 1, fm);
    RecoupleTensor(4, s, fm);
  }
  if (s[0].j != s[1].j) return 0.0;
  if (s[2].j != s[3].j) return 0.0;
  if (md <= 1) {
    FixJsZ(s, fm);
    fm->js[5] = 0;
    fm->js[6] = 0;
    fm->js[7] = 0;
    EvaluateTensor(ns, sbra, sket, s, 1, fm);
  }
  /* 
  ** these conditions must be placed after the above block,
  ** so that the recouping coeff. are calculated the first
  ** time even these parity conditions are not satisfied.
  */
  if (IsOdd((s[0].kl+s[1].kl)/2)) return 0.0;
  if (IsOdd((s[2].kl+s[3].kl)/2)) return 0.0;

  y = fm->coeff;
  if (fabs(y) < EPS10) return 0.0;
  ResidualPotential(&c1, k0, k1);
  c1 += QED1E(k0, k1);
  c1 *= sqrt(s[0].j+1.0);
  if (k3 != k0 || k2 != k1) {
    ResidualPotential(&c2, k2, k3);
    c2 += QED1E(k2, k3);
    c2 *= sqrt(s[3].j+1.0);
  } else {
    c2 = c1;
  }	  
  c = y*c1*c2; 
  return c;
}

void DeltaH22M2Loop(double *h1, double *h2, int ns, 
		    SHELL *bra, SHELL *ket, SHELL_STATE *sbra, SHELL_STATE *sket,
		    int n0, int *b0, int n1, int *b1, int n, int *ng, 
		    int n2, int *ng2, int ph,
		    int ia, int ib, int ic, int id, int ik, int im,
		    FORMULA *fm, double *a, int *j1, int *j2) {
  double c, d1, d2;
  int ip, iq;
  int ks1[4], ks2[4];
  int i, i1, i2, md, m1, m2;
  INTERACT_SHELL s[8];
  ORBITAL *o[8];
  
  ip = im;
  iq = ik;
  ks1[0] = b0[ia];
  ks1[1] = b0[ib];
  ks1[2] = b1[ik];
  ks1[3] = b1[im];
  ks2[0] = b1[iq];
  ks2[1] = b1[ip];
  ks2[2] = b0[ic];
  ks2[3] = b0[id];
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
      o[i] = GetOrbital(b1[ik]);
    } else if (i == 3 || i == 6) {
      o[i] = GetOrbital(b1[im]);
    } else {
      o[i] = GetOrbital(b0[s[i].index-2]);
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
  if (m1 > m2) {
    i = m1;
    m1 = m2;
    m2 = i;
  }  
  i1 = IBisect(m1, n, ng);
  if (i1 < 0) return;
  i2 = IBisect(m2-m1, n2, ng2);
  if (i2 < 0) return;
  if (s[1].j != *j1 || s[3].j != *j2) {
    md = (*j1 < 0 || *j2 < 0)?0:1;
    *j1 = s[1].j;
    *j2 = s[3].j;
  } else {
    md = 2;
  }
  c = H22Term(ns, sbra, sket, s, ks1, ks2, fm, a, md);
  if (fabs(c) < EPS10) return;
  if (IsOdd(ph)) c = -c;
  d1 = o[5]->energy+o[7]->energy-o[4]->energy-o[6]->energy;
  d2 = o[0]->energy+o[2]->energy-o[1]->energy-o[3]->energy;
  i = i1*n2 + i2;
  h1[i] += c/d1;
  h2[i] += c/d2;
  /*
  printf("22:2 %d %d %d %d %d %d %d %d %10.3E %10.3E %12.5E\n",
	 ks1[0],ks1[1],ks1[2],ks1[3],
	 ks2[0],ks2[1],ks2[2],ks2[3],d1,c,c/d1);
  */
}

void DeltaH22M2(double *h1, double *h2, int ns,
		SHELL *bra, SHELL *ket, SHELL_STATE *sbra, SHELL_STATE *sket,
		int n0, int *b0, int n1, int *b1, int n, int *ng, 
		int n2, int *ng2, int nc, CONFIG **cs) {
  int ia, ib, ic, id, ik, im, j1, j2;
  int m1, m2, nj, *jp, i, k;
  int op[4], om[4], ph;
  double a[MKK*MKK];
  ORBITAL *o;
  FORMULA fm;

  if (n1 <= 0) return;
  jp = malloc(sizeof(int)*(n1+1));
  jp[0] = 0;
  nj = 1;
  o = GetOrbital(b1[0]);
  j1 = GetJFromKappa(o->kappa);
  for (i = 1; i < n1; i++) {
    o = GetOrbital(b1[i]);
    j2 = GetJFromKappa(o->kappa);
    if (j2 != j1) {
      jp[nj] = i;
      nj++;
      j1 = j2;
    }
  }
  jp[nj] = n1;
  fm.js[0] = 0;
  for (ia = 0; ia < n0; ia++) {
    for (ib = 0; ib <= ia; ib++) {
      for (ic = 0; ic < n0; ic++) {
	for (id = 0; id <= ic; id++) {
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
	  if (ng2[0] == 0) {
	    j1 = -1;
	    j2 = -1;
	    for (im = 0; im < n1; im++) {
	      ik = im;
	      o = GetOrbital(b1[im]);
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
	      DeltaH22M2Loop(h1, h2, ns, bra, ket, sbra, sket, 
			     n0, b0, n1, b1, n, ng, n2, ng2, ph,
			     ia, ib, ic, id, ik, im, &fm, a, &j1, &j2);
	    }
	  }
	  j1 = -1;
	  j2 = -1;
	  for (m1 = 0; m1 < nj; m1++) {
	    for (m2 = 0; m2 <= m1; m2++) {
	      for (im = jp[m1]; im < jp[m1+1]; im++) {
		o = GetOrbital(b1[im]);
		ket[0].n = o->n;
		ket[0].kappa = o->kappa;
		ket[0].nq = 0;
		for (ik = jp[m2]; ik < jp[m2+1]; ik++) {
		  if (im <= ik) continue;
		  o = GetOrbital(b1[ik]);
		  ket[1].n = o->n;
		  ket[1].kappa = o->kappa;
		  ket[1].nq = 0;
		  if (ket[0].n <= cs[nc]->n_csfs && 
		      ket[1].n <= cs[nc]->n_csfs) {
		    om[0] = id+2;
		    om[1] = ic+2;
		    op[0] = 0;
		    op[1] = 1;
		    k = CheckConfig(ns, ket, 2, op, 2, om, nc, cs);
		    if (k >= 0) continue;
		  }
		  DeltaH22M2Loop(h1, h2, ns, bra, ket, sbra, sket, 
				 n0, b0, n1, b1, n, ng, n2, ng2, ph,
				 ia, ib, ic, id, ik, im, &fm, a, &j1, &j2);
		}
	      }
	    }
	  }
	}
      }
    }
  }
  free(jp);
}

void DeltaH22M1(double *h1, double *h2, int ns,
		SHELL *bra, SHELL *ket, SHELL_STATE *sbra, SHELL_STATE *sket,
		int n0, int *b0, int n1, int *b1, int n, int *ng,
		int nc, CONFIG **cs) {
  int ia, ib, ic, id, ik, im, ip, iq;
  int op[4], om[4], ph, k, ks1[4], ks2[4];
  int i, m1, i1, j1, md;
  double c, d1, d2, a[MKK*MKK];
  INTERACT_SHELL s[8];
  ORBITAL *o[8];
  FORMULA fm;

  fm.js[0] = 0;
  if (n1 <= 0) return;
  for (ia = 0; ia < n0; ia++) {
    for (ib = 0; ib <= ia; ib++) {
      for (ic = 0; ic < n0; ic++) {
	for (id = 0; id <= ic; id++) {
	  for (ik = 0; ik < n0; ik++) {
	    for (iq = 0; iq < n0; iq++) {
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
	      j1 = -1;
	      for (ip = 0; ip < n1; ip++) {
		o[1] = GetOrbital(b1[ip]);
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
		ks1[0] = b0[ia];
		ks1[1] = b0[ib];
		ks1[2] = b1[im];
		ks1[3] = b0[ik];
		ks2[0] = b0[iq];
		ks2[1] = b1[ip];
		ks2[2] = b0[ic];
		ks2[3] = b0[id];
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
		    o[i] = GetOrbital(b1[ip]);
		    if (i == 1) {
		      m1 = o[1]->n;
		      i1 = IBisect(m1, n, ng);	    
		      if (i1 < 0) break;
		    }
		  } else {
		    o[i] = GetOrbital(b0[s[i].index-1]);
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
		if (s[1].j != j1) {
		  md = j1<0?0:1;
		  j1 = s[1].j;
		} else {
		  md = 2;
		}
		c = H22Term(ns, sbra, sket, s, ks1, ks2, &fm, a, md);
		if (fabs(c) < EPS10) continue;
		if (IsOdd(ph)) c = -c;
		
		d1 = o[5]->energy+o[7]->energy-o[4]->energy-o[6]->energy;
		d2 = o[0]->energy+o[2]->energy-o[1]->energy-o[3]->energy;
		h1[i1] += c/d1;
		h2[i1] += c/d2;
		/*
		printf("22:1 %d %d %d %d %d %d %d %d %10.3E %10.3E %12.5E\n", 
		       ks1[0],ks1[1],ks1[2],ks1[3],
		       ks2[0],ks2[1],ks2[2],ks2[3],d1,c,c/d1);
		*/
	      }
	    }
	  }
	}
      }
    }
  }
}
		  
void DeltaH22M0(double *h1, double *h2, int ns,
		SHELL *bra, SHELL *ket, SHELL_STATE *sbra, SHELL_STATE *sket,
		int n0, int *b0, int n, int *ng, int nc, CONFIG **cs) {
  int ia, ib, ic, id, ik, im, ip, iq, i, i1;
  int op[4], om[4], ph, k, ks1[4], ks2[4];
  double c, d1, d2, a[MKK*MKK];
  INTERACT_SHELL s[8];
  ORBITAL *o[8];  
  FORMULA fm;

  i1 = bra[0].n;
  i1 = IBisect(i1, n, ng);
  if (i1 < 0) return;
  fm.js[0] = 0;
  for (ia = 0; ia < n0; ia++) {
    for (ib = 0; ib <= ia; ib++) {
      for (ic = 0; ic < n0; ic++) {
	for (id = 0; id <= ic; id++) {
	  for (ik = 0; ik < n0; ik++) {
	    for (im = 0; im <= ik; im++) {
	      for (ip = 0; ip < n0; ip++) {
		for (iq = 0; iq <= ip; iq++) {
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
		  ks1[0] = b0[ia];
		  ks1[1] = b0[ib];
		  ks1[2] = b0[ik];
		  ks1[3] = b0[im];
		  ks2[0] = b0[ip];
		  ks2[1] = b0[iq];
		  ks2[2] = b0[ic];
		  ks2[3] = b0[id];
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
		    o[i] = GetOrbital(b0[s[i].index]);
		    if (IsEven(i)) s[i].n = 1;
		    else s[i].n = -1;
		    s[i].kappa = o[i]->kappa;
		    GetJLFromKappa(s[i].kappa, &(s[i].j), &(s[i].kl));
		    s[i].nq_bra = bra[s[i].index].nq;
		    s[i].nq_ket = ket[s[i].index].nq;
		    s[i].index = ns-s[i].index-1;
		  }
		  c = H22Term(ns, sbra, sket, s, ks1, ks2, &fm, a, 0);
		  if (fabs(c) < EPS10) continue;
		  if (IsOdd(ph)) c = -c;
		  d1 = o[5]->energy+o[7]->energy-o[4]->energy-o[6]->energy;
		  d2 = o[0]->energy+o[2]->energy-o[1]->energy-o[3]->energy;
		  h1[i1] += c/d1;
		  h2[i1] += c/d2;
		  /*
		  printf("22:0 %d %d %d %d %d %d %d %d %10.3E %10.3E %12.5E\n", 
			 ks1[0],ks1[1],ks1[2],ks1[3],
			 ks2[0],ks2[1],ks2[2],ks2[3],d1,c,c/d1);
		  */
		}
	      }
	    }
	  }
	}
      }
    }
  }
}

void DeltaH12M1(double *h1, double *h2, int ns,
		SHELL *bra, SHELL *ket, SHELL_STATE *sbra, SHELL_STATE *sket,
		int n0, int *b0, int n1, int *b1, int n, int *ng,
		int nc, CONFIG **cs) {
  int ia, ib, ic, id, ik, im, j1, md;
  int m1, i1, k, i, op[3], om[3], ks[4], ph;
  double c, d1, d2, a[MKK];
  INTERACT_SHELL s[6];
  ORBITAL *o[6];
  FORMULA fm;

  fm.js[0] = 0;
  if (n1 <= 0) return;
  for (ia = 0; ia < n0; ia++) {
    if (bra[ia+1].nq == 0) continue;
    for (ib = 0; ib < n0; ib++) {
      for (ic = 0; ic <= ib; ic++) {
	for (id = 0; id < n0; id++) {	  
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
	  j1 = -1;
	  for (ik = 0; ik < n1; ik++) {
	    o[1] = GetOrbital(b1[ik]);
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
	    ks[0] = b0[id];
	    ks[1] = b1[im];
	    ks[2] = b0[ib];
	    ks[3] = b0[ic];
	    s[0].index = ia+1;
	    s[1].index = 0;
	    s[2].index = id+1;
	    s[3].index = ib+1;
	    s[4].index = 0;
	    s[5].index = ic+1;
	    i1 = 0;
	    for (i = 0; i < 6; i++) {
	      if (i == 1 || i == 4) {		
		o[i] = GetOrbital(b1[ik]);
		if (i == 1) {
		  m1 = o[1]->n;
		  i1 = IBisect(m1, n, ng);	    
		  if (i1 < 0) break;
		}
	      } else {
		o[i] = GetOrbital(b0[s[i].index-1]);
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
	    if (s[1].j != j1) {
	      md = j1<0?0:1;
	      j1 = s[1].j;
	    } else {
	      md = 2;
	    }
	    c = H12Term(ns, sbra, sket, s, ks, b0[ia], b1[ik], &fm, a, md);
	    if (fabs(c) < EPS10) continue;
	    if (IsOdd(ph)) c = -c;
	    d1 = o[3]->energy + o[5]->energy - o[2]->energy - o[4]->energy;
	    d2 = o[0]->energy - o[1]->energy;
	    h1[i1] += c/d1;
	    h2[i1] += c/d2;
	    /*
	    printf("12:1 %d %d %d %d %d %d %10.3E %10.3E %12.5E %10.3E %12.5E\n", 
		   b0[ia],b1[ik],ks[0],ks[1],ks[2],ks[3],d1,c,c/d1,d2,c/d2);
	    */
	  }
	}
      }
    }
  }
}

void DeltaH12M0(double *h1, double *h2, int ns,
		SHELL *bra, SHELL *ket, SHELL_STATE *sbra, SHELL_STATE *sket,
		int n0, int *b0, int n, int *ng, int nc, CONFIG **cs) {
  int ia, ib, ic, id, ik, im;
  int i1, k, i, op[3], om[3], ks[4], ph;
  double c, d1, d2, a[MKK];
  INTERACT_SHELL s[6];
  ORBITAL *o[6];
  FORMULA fm;
  
  i1 = bra[0].n;
  i1 = IBisect(i1, n, ng);
  if (i1 < 0) return;
  fm.js[0] = 0;
  for (ia = 0; ia < n0; ia++) {
    if (bra[ia].nq == 0) continue;
    for (ib = 0; ib < n0; ib++) {
      if (ket[ib].nq == 0) continue;
      for (ic = 0; ic <= ib; ic++) {
	if (ket[ic].nq == 0) continue;
	if (ic == ib && ket[ic].nq == 1) continue;
	for (id = 0; id < n0; id++) {
	  for (ik = 0; ik < n0; ik++) {
	    for (im = 0; im <= id; im++) {
	      if (ik == ia) continue;	      
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
	      ks[0] = b0[id];
	      ks[1] = b0[im];
	      ks[2] = b0[ib];
	      ks[3] = b0[ic];
	      s[0].index = ia;
	      s[1].index = ik;
	      s[2].index = id;
	      s[3].index = ib;
	      s[4].index = im;
	      s[5].index = ic;
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
		o[i] = GetOrbital(b0[s[i].index]);
		if (IsEven(i)) s[i].n = 1;
		else s[i].n = -1;
		s[i].kappa = o[i]->kappa;
		GetJLFromKappa(s[i].kappa, &(s[i].j), &(s[i].kl));
		s[i].nq_bra = bra[s[i].index].nq;
		s[i].nq_ket = ket[s[i].index].nq;
		s[i].index = ns-s[i].index-1;
	      }
	      c = H12Term(ns, sbra, sket, s, ks, b0[ia], b0[ik], &fm, a, 0);
	      if (fabs(c) < EPS10) continue;
	      if (IsOdd(ph)) c = -c;
	      d1 = o[3]->energy + o[5]->energy - o[2]->energy - o[4]->energy;
	      d2 = o[0]->energy - o[1]->energy;
	      h1[i1] += c/d1;
	      h2[i1] += c/d2;
	      /*
	      printf("12:0 %d %d %d %d %d %d %10.3E %10.3E %12.5E %10.3E %12.5E\n", 
		     b0[ia],b0[ik],ks[0],ks[1],ks[2],ks[3],d1,c,c/d1,d2,c/d2);
	      */
	    }
	  }
	}
      }
    }
  }
}

void DeltaH11M0(double *h1, double *h2, int ns, 
		SHELL *bra, SHELL *ket, SHELL_STATE *sbra, SHELL_STATE *sket,
		int n0, int *b0, int n, int *ng,
		int nc, CONFIG **cs) {
  int ia, ib, ik, im, k, k0, k1, k2, k3, i, i1;
  int op[2], om[2], ph;
  double y, d1, d2;
  INTERACT_SHELL s[4];
  ORBITAL *o0, *o1, *o2, *o3;
  FORMULA fm;
   
  i1 = bra[0].n;
  i1 = IBisect(i1, n, ng);
  if (i1 < 0) return;
  fm.js[0] = 0;
  for (ia = 0; ia < n0; ia++) {
    if (bra[ia].nq == 0) continue;
    for (ib = 0; ib < n0; ib++) { 
      if (ket[ib].nq == 0) continue;
      for (ik = 0; ik < n0; ik++) {
	for (im = 0; im < n0; im++) {
	  if (ia == ik) {
	    continue;
	  }
	  if (ib == im) {
	    continue;
	  }
	  op[0] = ia;
	  op[1] = im;
	  om[0] = ib;
	  om[1] = ik;
	  ph = CheckInteraction(ns, bra, ket, 2, op, 2, om);
	  if (ph < 0) continue;
	  op[0] = im;
	  om[0] = ib;
	  k = CheckConfig(ns, ket, 1, op, 1, om, nc, cs);	  
	  if (k >= 0) continue;
	  op[0] = ik;
	  om[0] = ia;
	  k = CheckConfig(ns, bra, 1, op, 1, om, 0, NULL);
	  if (k >= 0) continue;
	  k0 = b0[ia];
	  k1 = b0[ik];
	  k2 = b0[im];
	  k3 = b0[ib];
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
	  y = H11Term(ns, sbra, sket, s, k0, k1, k2, k3, &fm, 0);
	  if (fabs(y) < EPS10) continue;
	  if (IsOdd(ph)) y = -y;
	  d1 = o3->energy - o2->energy;
	  d2 = o0->energy - o1->energy;
	  h1[i1] += y/d1;
	  h2[i1] += y/d2;
	  /*
	  printf("11:0 %d %d %d %d %10.3E %10.3E %12.5E\n", k0,k1,k2,k3,d1,y,y/d1);
	  */
	}
      }
    }
  }
}
  

void DeltaH11M1(double *h1, double *h2, int ns, 
		SHELL *bra, SHELL *ket, SHELL_STATE *sbra, SHELL_STATE *sket,
		int n0, int *b0, int n1, int *b1, int n, int *ng,
		int nc, CONFIG **cs) {
  int ia, ib, ik, im, k, k0, k1, k2, k3, i, i1, m1;
  int op[2], om[2], ph, j1, md;
  double y, d1, d2;
  INTERACT_SHELL s[4];
  ORBITAL *o0, *o1, *o2, *o3;
  FORMULA fm;

  fm.js[0] = 0;
  if (n1 <= 0) return;
  for (ia = 0; ia < n0; ia++) {
    if (bra[ia+1].nq == 0) continue;
    for (ib = 0; ib < n0; ib++) { 
      if (ket[ib+1].nq == 0) continue;
      op[0] = ia;
      om[0] = ib;
      ph = CheckInteraction(ns-1, bra+1, ket+1, 1, op, 1, om);
      if (ph < 0) continue;
      j1 = -1;
      for (ik = 0; ik < n1; ik++) {
	o1 = GetOrbital(b1[ik]);
	ket[0].n = o1->n;
	ket[0].kappa = o1->kappa;
	ket[0].nq = 0;
	if (ket[0].n <= cs[nc]->n_csfs) {
	  om[0] = ib+1;
	  op[0] = 0;
	  k = CheckConfig(ns, ket, 1, op, 1, om, nc, cs);
	  if (k >= 0) continue;
	}
	im = ik;
	k0 = b0[ia];
	k1 = b1[ik];
	k2 = b1[im];
	k3 = b0[ib];
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
	i1 = IBisect(m1, n, ng);
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
	if (s[1].j != j1) {
	  md = j1<0?0:1;
	  j1 = s[1].j;
	} else {
	  md = 2;
	}
	y = H11Term(ns, sbra, sket, s, k0, k1, k2, k3, &fm, md);
	if (fabs(y) < EPS10) continue;
	if (IsOdd(ph)) y = -y;
	d1 = o3->energy - o2->energy;
	d2 = o0->energy - o1->energy;
	h1[i1] += y/d1;
	h2[i1] += y/d2;
	/*
	printf("11:1 %d %d %d %d %10.3E %10.3E %12.5E\n", k0,k1,k2,k3,d1,y,y/d1);
	*/
      }
    }
  }
}
  
int StructureMBPT1(char *fn, char *fn1, int nkg, int *kg, 
		   int kmax, int n, int *ng, int n2, int *ng2, int nkg0) {
  int *bas, *bas0, *bas1, *bas2, nlevels, nb, n3, q;
  int i, j, k, i0, i1, n0, n1, isym, ierr, nc, m;
  int pp, jj, nmax, na, *ga;
  SYMMETRY *sym;
  STATE *st;
  HAMILTON *h;
  SHELL *bra, *ket, *bra1, *ket1, *bra2, *ket2;
  SHELL_STATE *sbra, *sket, *sbra1, *sket1, *sbra2, *sket2;
  CONFIG **cs, cfg;
  CONFIG_GROUP *g;
  double a, b, c, *mix, *hab, *hba;
  double *h0, *heff, *hab1, *hba1, *dw;
  FILE *f;
					  
  ierr = 0;
  n3 = mbpt_n3;
  if (nkg0 <= 0 || nkg0 > nkg) nkg0 = nkg;

  n = ConstructNGrid(n, &ng);
  n2 = ConstructNGrid(n2, &ng2);
  na = n*n2+n;
  ga = malloc(sizeof(int)*na);
  k = 0;
  for (i = 0; i < n; i++) {
    ga[k++] = ng[i];
    for (j = 0; j < n2; j++) {
      if (ng2[j] == 0) continue;
      ga[k++] = ng[i]+ng2[j];
    }
  }
  na = SortUnique(k, ga);
  nb = RadialBasisMBPT(kmax, na, ga, &bas);
  free(ga);
  if (nb < 0) return -1;

  f = fopen(fn1, "w");
  if (f == NULL) {
    printf("cannot open file %s\n", fn1);
    free(bas);
    return -1;
  }
  fwrite(&n, sizeof(int), 1, f);
  fwrite(ng, sizeof(int), n, f);
  fwrite(&n2, sizeof(int), 1, f);
  fwrite(ng2, sizeof(int), n2, f);
  fwrite(&n3, sizeof(int), 1, f);

  dw = malloc(sizeof(double)*(n+n2)*4);
  hab = malloc(sizeof(double)*n*n2);
  hba = malloc(sizeof(double)*n*n2);
  hab1 = malloc(sizeof(double)*n);
  hba1 = malloc(sizeof(double)*n);

  bas0 = malloc(sizeof(int)*nb);
  bas1 = malloc(sizeof(int)*nb);
  bas2 = malloc(sizeof(int)*nb);

  nc = 0;
  for (i = 0; i < nkg; i++) {
    g = GetGroup(kg[i]);
    nc += g->n_cfgs;
  }
  cs = malloc(sizeof(CONFIG *)*(nc+2));
  k = 0;
  i1 = 0;
  nmax = 0;
  for (i = 0; i < nkg; i++) {
    g = GetGroup(kg[i]);    
    for (j = 0; j < g->n_cfgs; j++) {
      cs[k] = GetConfigFromGroup(i, j);
      if (cs[k]->n_shells > i1) i1 = cs[k]->n_shells;
      if (cs[k]->shells[0].n > nmax) nmax = cs[k]->shells[0].n;
      k++;
    }
  }
  /* use cfg.n_csfs to store the maximum n-value of the configuration in kg */
  cfg.n_csfs = nmax;
  cs[k] = &cfg;
  cs[k]->n_shells = i1;
  cs[k]->shells = malloc(sizeof(SHELL)*(i1+2));
  nlevels = GetNumLevels();
  for (isym = 0; isym < MAX_SYMMETRIES; isym++) {
    /*
    if (isym == 5) {
      STATE *st;
      CONFIG *pcf;
      double y1[10],y2[10],y3[10], d1, d2;
      k = ConstructHamiltonDiagonal(isym, nkg, kg, 0);
      h = GetHamilton();
      sym = GetSymmetry(isym);
      hab = malloc(sizeof(double)*h->dim);
      hba = malloc(sizeof(double)*h->dim);
      for (i = 2; i < sym->n_states; i++) {
	HamiltonElement1E2E(isym, 0, i, &a, &b);
	hab[i] = a+b;
	HamiltonElement1E2E(isym, 1, i, &a, &b);
	hba[i] = a+b;
      }
      for (i = 0; i < 10; i++) {
	y1[i] = 0.0;
	y2[i] = 0.0;
	y3[i] = 0.0;
      }
      for (i = 2; i < sym->n_states; i++) {
	st = ArrayGet(&(sym->states), i);
	pcf = GetConfig(st);
	k = pcf->shells[0].n;
	d1 = hab[i]*hab[i];
	d2 = hba[i]*hba[i];
	a = hab[i]*hba[i];
	b = h->hamilton[1] - h->hamilton[i];
	c = h->hamilton[0] - h->hamilton[i];
	printf("%3d %2d %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %d %d %d\n", 
	       i, k, d1, d2, a, c, b, d1/c, d2/b, a/b,
	       st->kgroup,st->kcfg,st->kstate);
	y1[k] += d1/c;
	y2[k] += d2/b;
	y3[k] += a/b;
      }
      for (i = 1; i < 10; i++) {
	printf("n: %2d %15.8E %15.8E %15.8E\n", i, y1[i], y2[i], y3[i]);
      }
      goto ERROR;
    } else {
      continue;
    }    
    */
    k = ConstructHamilton(isym, nkg0, nkg, kg, 0, NULL);
    if (k < 0) continue;
    if (mbpt_isym >= 0 && mbpt_isym != isym) {
      k = 0;
      fwrite(&isym, sizeof(int), 1, f);
      fwrite(&k, sizeof(int), 1, f);
      continue;
    }
    h = GetHamilton();
    h0 = malloc(sizeof(double)*h->hsize);
    memcpy(h0, h->hamilton, sizeof(double)*h->hsize);
    if (DiagnolizeHamilton() < 0) {
      printf("Diagnolizing Hamiltonian Error\n");
      ierr = -1;
      goto ERROR;
    }    
    h->heff = malloc(sizeof(double)*h->dim*h->dim);    
    heff = h->heff;
    sym = GetSymmetry(isym);
    DecodePJ(isym, &pp, &jj);
    fwrite(&isym, sizeof(int), 1, f);
    fwrite(&(h->dim), sizeof(int), 1, f);
    fflush(f);
    for (j = 0; j < h->dim; j++) {
      for (i = 0; i <= j; i++) {
	/*
	if (isym != 2 || i != 1 || j != 6) continue;
	printf("isym: %d %d %d %d %d\n", isym, i, j, h->basis[i], h->basis[j]);
	*/
	c = 0;
	m = -1;
	mix = h->mixing + h->dim;
	for (k = 0; k < h->dim; k++) {
	  if (nkg0 != nkg) {
	    q = GetPrincipleBasis(mix, h->dim, NULL);
	    st = (STATE *) ArrayGet(&(sym->states), h->basis[q]);
	    if (!InGroups(st->kgroup, nkg0, kg)) {
	      mix += h->n_basis;
	      continue;
	    }
	  }
	  m++;
	  if (mbpt_ilev >= 0 && m != mbpt_ilev) {
	    mix += h->n_basis;
	    continue;
	  }
	  a = fabs(mix[i]*mix[j]);
	  if (a > c) c = a;
	  mix += h->n_basis;
	}
	if (i == j) a = mbpt_mcut;
	else a = 5.0*mbpt_mcut;
	if (a > 0.1) a = 0.1;
	printf("%3d %d %2d %3d %3d %3d %10.3E %10.3E\n", 
	       isym, pp, jj, h->dim, i, j, c, a);
	fflush(stdout);
	if (c < a) {
	  k = j*(j+1)/2 + i;
	  a = h0[k];
	  k = j*h->dim + i;
	  heff[k] = a;
	  if (i < j) {
	    k = i*h->dim + j;
	    heff[k] = a;
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
	n0 = PadStates(sym, h->basis[i], h->basis[j], 
		       &bra, &ket, &sbra, &sket);	
	bra1 = bra + 1;
	ket1 = ket + 1;
	bra2 = bra + 2;
	ket2 = ket + 2;
	sbra1 = sbra + 1;
	sket1 = sket + 1;
	sbra2 = sbra + 2;
	sket2 = sket + 2;	
	for (k = 0; k < n0; k++) {	  
	  bas0[k] = OrbitalIndex(bra2[k].n, bra2[k].kappa, 0.0);
	  bas2[k] = bas0[k];
	}
	qsort(bas2, n0, sizeof(int), CompareInt);
	n1 = 0;
	for (m = 0; m < nb; m++) {
	  k = IBisect(bas[m], n0, bas2);
	  if (k >= 0) continue;
	  bas1[n1] = bas[m];
	  n1++;
	}
	for (i0 = 0; i0 < n; i0++) {
	  hab1[i0] = 0.0;
	  hba1[i0] = 0.0;
	}
	i1 = n*n2;
	for (i0 = 0; i0 < i1; i0++) {
	  hab[i0] = 0.0;
	  hba[i0] = 0.0;
	}
	if (n3 != 2) {
	  /* make sure that DeltaH12## routines must be called first
	     because of the factor of 2 multiplication */	
	  DeltaH12M0(hab1, hba1, n0, bra2, ket2, sbra2, sket2,
		     n0, bas0, n, ng, nc, cs);	
	  DeltaH12M1(hab1, hba1, n0+1, bra1, ket1, sbra1, sket1,
		     n0, bas0, n1, bas1, n, ng, nc, cs);	
	  if (i < j) {
	    DeltaH12M0(hba1, hab1, n0, ket2, bra2, sket2, sbra2,
		       n0, bas0, n, ng, nc, cs);	  
	    DeltaH12M1(hba1, hab1, n0+1, ket1, bra1, sket1, sbra1,
		       n0, bas0, n1, bas1, n, ng, nc, cs);
	  } else {
	    for (i0 = 0; i0 < n; i0++) {
	      hab1[i0] *= 2;
	      hba1[i0] *= 2;
	    }
	  }	
	  
	  DeltaH11M0(hab1, hba1, n0, bra2, ket2, sbra2, sket2, 
		     n0, bas0, n, ng, nc, cs);	
	  DeltaH11M1(hab1, hba1, n0+1, bra1, ket1, sbra1, sket1, 
		     n0, bas0, n1, bas1, n, ng, nc, cs);
	  
	  DeltaH22M0(hab1, hba1, n0, bra2, ket2, sbra2, sket2,
		     n0, bas0, n, ng, nc, cs);	
	  DeltaH22M1(hab1, hba1, n0+1, bra1, ket1, sbra1, sket1,
		     n0, bas0, n1, bas1, n, ng, nc, cs);
	}      
	if (n3 != 1) {
	  DeltaH22M2(hab, hba, n0+2, bra, ket, sbra, sket,
		     n0, bas0, n1, bas1, n, ng, n2, ng2, nc, cs);	
	}
	a = sqrt(jj+1.0);
	for (i0 = 0; i0 < n; i0++) {
	  hab1[i0] /= a;
	  hba1[i0] /= a;
	}
	for (i0 = 0; i0 < i1; i0++) {
	  hab[i0] /= a;
	  hba[i0] /= a;
	}
	k = j*(j+1)/2 + i;
	a = h0[k];
	k = j*h->dim + i;
	b = SumInterpH(n, ng, n2, ng2, hab, hab1, dw);
	heff[k] = a + b;
	if (i < j) {
	  k = i*h->dim + j;
	  c = SumInterpH(n, ng, n2, ng2, hba, hba1, dw);
	  heff[k] = a + c;
	} else {
	  c = b;
	}
	fwrite(&i, sizeof(int), 1, f);
	fwrite(&j, sizeof(int), 1, f);
	fwrite(&a, sizeof(double), 1, f);
	fwrite(&b, sizeof(double), 1, f);
	fwrite(&c, sizeof(double), 1, f);
	if (n3 != 2) {
	  fwrite(hab1, sizeof(double), n, f);
	  if (i != j) {
	    fwrite(hba1, sizeof(double), n, f);
	  }
	}
	if (n3 != 1) {
	  fwrite(hab, sizeof(double), i1, f);
	  if (i != j) {
	    fwrite(hba, sizeof(double), i1, f);
	  }
	}
	fflush(f);
	free(bra);
	free(ket);
	free(sbra);
	free(sket);
      }
    }
    if (DiagnolizeHamilton() < 0) {
      printf("Diagnolizing Hamiltonian Error\n");
      ierr = -1;
      goto ERROR;
    }
    AddToLevels(nkg0, kg);
    free(h->heff);
    h->heff = NULL;
    free(h0);
  }

  SortLevels(nlevels, -1);
  SaveLevels(fn, nlevels, -1);

 ERROR:
  fclose(f);
  free(bas);
  free(bas0);
  free(bas1);
  free(bas2);
  free(cs[nc]->shells);
  free(cs);
  free(dw);
  free(hab);
  free(hba);
  free(hab1);
  free(hba1);

  return ierr;
}

int ReadMBPT(int nf, FILE *f[], MBPT_HAM *mbpt, int m) {
  int i, nb, nn2;

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
      mbpt[i].hab1 = malloc(sizeof(double)*mbpt[i].n);
      mbpt[i].hba1 = malloc(sizeof(double)*mbpt[i].n);
      mbpt[i].hab = malloc(sizeof(double)*mbpt[i].n*mbpt[i].n2);
      mbpt[i].hba = malloc(sizeof(double)*mbpt[i].n*mbpt[i].n2);
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
	nb = fread(mbpt[i].hab1, sizeof(double), mbpt[i].n, f[i]);
	if (nb != mbpt[i].n) return -1;
	if (mbpt[i].ibra != mbpt[i].iket) {
	  nb = fread(mbpt[i].hba1, sizeof(double), mbpt[i].n, f[i]);
	  if (nb != mbpt[i].n) return -1;
	} else {
	  memcpy(mbpt[i].hba1, mbpt[i].hab1, sizeof(double)*mbpt[i].n);
	}
      } 
      if (mbpt[i].n3 != 1) {
	nn2 = mbpt[i].n*mbpt[i].n2;
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
  
  return 0;
}

void CombineMBPT(int nf, MBPT_HAM *mbpt, double *hab1, double *hba1,
		 double **hab, double **hba, int n0, int *ng0,
		 int n, int *ng, int *n2, int **ng2) {
  int m, i, j, k, q, r;

  for (m = 0; m < nf; m++) {
    if (mbpt[m].ibra < 0 || mbpt[m].iket < 0) continue;
    if (mbpt[m].n3 != 2) {
      for (i = 0; i < mbpt[m].n; i++) {
	k = IBisect(mbpt[m].ng[i], n0, ng0);
	if (k < 0) continue;
	hab1[k] = mbpt[m].hab1[i];
	hba1[k] = mbpt[m].hba1[i];
      }
    }
    if (mbpt[m].n3 != 1) {
      for (i = 0; i < mbpt[m].n; i++) {
	k = IBisect(mbpt[m].ng[i], n, ng);
	if (k < 0) continue;
	r = i*mbpt[m].n2;
	for (j = 0; j < mbpt[m].n2; j++) {
	  q = IBisect(mbpt[m].ng2[j], n2[k], ng2[k]);
	  if (q < 0) continue;
	  hab[k][q] = mbpt[m].hab[r+j];
	  hba[k][q] = mbpt[m].hba[r+j];
	}
      }
    }
  }
}

int StructureReadMBPT(char *fn, char *fn2, int nf, char *fn1[],
		      int nkg, int *kg, int nkg0) {
  int ierr, m, i, j, k, k0, k1, nlevels;
  int isym, pp, jj, r, q, n, *ng, n0, *ng0, *n2, **ng2;
  double *dw, *z1, *z2, *x, *y, *z, *t, *heff, **hab, **hba;
  HAMILTON *h;
  SYMMETRY *sym;
  MBPT_HAM *mbpt;
  FILE **f1, *f2;
  
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
  if (ierr < 0) return -1;

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
    }
    hab = malloc(sizeof(double *)*n);
    hba = malloc(sizeof(double *)*n);
    for (i = 0; i < n; i++) {
      if (n2[i] > 0) {
	hab[i] = malloc(sizeof(double)*n2[i]);
	hba[i] = malloc(sizeof(double)*n2[i]);
      }
    }
  }
  dw = malloc(sizeof(double)*(n0+n)*5);
  z1 = dw;
  z2 = z1 + n+n0;
  x = z2 + n+n0;
  t = x + n+n0;
  y = t + n+n0;

  nlevels = GetNumLevels();
  for (isym = 0; isym < MAX_SYMMETRIES; isym++) {
    k = ConstructHamiltonDiagonal(isym, nkg0, kg, -1);
    if (k < 0) continue;
    k = ConstructHamiltonDiagonal(isym, nkg, kg, -1);
    h = GetHamilton();
    h->heff = malloc(sizeof(double)*h->dim*h->dim);
    heff = h->heff;
    sym = GetSymmetry(isym);
    DecodePJ(isym, &pp, &jj);
    ierr = ReadMBPT(nf, f1, mbpt, 1);
    if (ierr < 0) {
      printf("Error reading MBPT %d\n", isym);
      goto ERROR;
    }
    if (mbpt[0].dim == 0) continue;
    if (mbpt_isym >= 0 && mbpt_isym != isym) {
      for (j = 0; j < h->dim; j++) {
	for (i = 0; i <= j; i++) {	  
	  ierr = ReadMBPT(nf, f1, mbpt, 2);
	}
      }
      continue;
    }
    for (j = 0; j < h->dim; j++) {
      for (i = 0; i <= j; i++) {
	k = j*(j+1)/2 + i;
	k0 = j*h->dim + i;
	k1 = i*h->dim + j;
	ierr = ReadMBPT(nf, f1, mbpt, 2);
	if (ierr < 0) {
	  printf("Error reading MBPT %d %d %d\n", isym, i, j);
	  goto ERROR;	
	}
	heff[k0] = mbpt[0].a;
	heff[k1] = mbpt[0].a;
	if (mbpt[0].ibra < 0 || mbpt[0].iket < 0) goto OUT;
	CombineMBPT(nf, mbpt, z1, z2, hab, hba, n0, ng0, n, ng, n2, ng2);
	for (q = 0; q < n0; q++) {
	  fprintf(f2, "  %3d %3d %3d %3d %3d %12.5E %12.5E\n", 
		  isym, i, j, 0, ng0[q], z1[q], z2[q]);
	}
	for (r = 0; r < n; r++) {
	  for (q = 0; q < n2[r]; q++) {
	    fprintf(f2, "  %3d %3d %3d %3d %3d %12.5E %12.5E\n", 
		    isym, i, j, ng[r], ng[r]+ng2[r][q], 
		    hab[r][q], hba[r][q]);
	  }
	}
	if (n0 > 0) {
	  for (q = 0; q < n0; q++) {
	    x[q] = ng0[q];
	  }
	  heff[k0] += SumInterp1D(n0, z1, x, t, y);
	  if (i < j) heff[k1] += SumInterp1D(n0, z2, x, t, y);
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
	    }
	  }
	  for (q = 0; q < n; q++) {
	    x[q] = ng[q];
	    fprintf(f2, "  %3d %3d %3d %3d %3d %12.5E %12.5E\n",
		    isym, i, j, -1, ng[q], z1[q], z2[q]);
	  }
	  heff[k0] += SumInterp1D(n, z1, x, t, y);
	  if (i < j) heff[k1] += SumInterp1D(n, z2, x, t, y);
	}
      OUT:
	fprintf(f2, "# %3d %3d %3d %3d %3d %12.5E %12.5E %15.8E\n",
		isym, pp, jj, i, j,
		heff[k0]-mbpt[0].a, heff[k1]-mbpt[0].a, mbpt[0].a);
      }
    }
    if (DiagnolizeHamilton() < 0) {
      printf("Diagnolizing Hamiltonian Error\n");
      ierr = -1;
      goto ERROR;
    }
    AddToLevels(nkg0, kg);
    free(h->heff);
    h->heff = NULL;
  }
  SortLevels(nlevels, -1);
  SaveLevels(fn, nlevels, -1);
  
 ERROR:
  fclose(f2);
  for (m = 0; m < nf; m++) {
    fclose(f1[m]);    
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
      }
    }
    free(ng2);
    free(hab);
    free(hba);
    free(n2);
  }
  if (n0 > 0) {
    free(ng0);
  }

  return ierr;
}
  
int StructureMBPT(char *fn, char *fn1, int n, int *s0, int kmax, 
		  int n1, int *nm, int n2, int *nmp, int n0, char *gn0) {
  if (gn0 == NULL) {
    StructureMBPT1(fn, fn1, n, s0, kmax, n1, nm, n2, nmp, n0);
  } else {
    StructureMBPT0(fn, fn1, n, s0, kmax, n1, nm, n2, nmp, gn0);
  }
  
  return 0;
}
