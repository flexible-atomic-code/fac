#include "recombination.h"
#include "time.h"

static char *rcsid="$Id: recombination.c,v 1.25 2001/10/12 03:15:01 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

#define NPARAMS 5

static int egrid_type = -1;
static int usr_egrid_type = -1;
static int output_format = 0;

static int interpolate_egrid;
static int n_egrid = 0;
static double egrid[MAXNE];
static double log_egrid[MAXNE];
static double sigma[MAXNE];
static double xegrid[MAXNTE][MAXNE];
static double log_xegrid[MAXNTE][MAXNE];

static int n_usr = 0;
static double usr_egrid[MAXNUSR];
static double log_usr[MAXNUSR];
static double xusr[MAXNUSR];
static double log_xusr[MAXNUSR];

static int n_tegrid = 0;
static double tegrid[MAXNTE];
static double log_te[MAXNTE];

static MULTI *pk_array;
static MULTI *qk_array;
static n_qk;

static struct {
  int n_spec;
  int n_frozen;
  int n_max;
  int max_kl;
  int kl_interp;
  int nkl0;
  int nkl;
  int pw_limits[2];
  int kl[MAXNKL+1];
  int kappa0[(MAXNKL+1)*2];
} pw_scratch = {6, 6, 100, MAXNKL, 10, 0, 0, {0, MAXKL}};

double ai_cut = EPS8;

static REC_COMPLEX rec_complex[MAX_COMPLEX];
int n_complex = 0;

void uvip3p_(int *np, int *ndp, double *x, double *y, 
	     int *n, double *xi, double *yi);


int SetAICut(double c) {
  ai_cut = c;
}

int SetRRTEGrid(int n, double emin, double emax) {
  n_tegrid = SetTEGrid(tegrid, log_te, n, emin, emax);
  return n_tegrid;
}

int SetRRTEGridDetail(int n, double *x) {
  n_tegrid = SetTEGridDetail(tegrid, log_te, n, x);
  return n_tegrid;
}

int SetUsrPEGridType(int type) {
  if (type >= 0) usr_egrid_type = type;
  return 0;
}

int SetPEGridDetail(int n, double *xg) {
  n_egrid = SetEGridDetail(egrid, log_egrid, n, xg);
  return n_egrid;
}

int SetPEGrid(int n, double emin, double emax, double eth) {
  n_egrid = SetEGrid(egrid, log_egrid, n, emin, emax, eth);
  return n_egrid;
}

int SetUsrPEGridDetail(int n, double *xg) {
  if (n > MAXNUSR) {
    printf("Max # of grid points reached \n");
    return -1;
  }
  n_usr = SetEGridDetail(usr_egrid, log_usr, n, xg);
  return n_usr;
}
  						      
int SetUsrPEGrid(int n, double emin, double emax, double eth) {
  if (n > MAXNUSR) {
    printf("Max # of grid points reached \n");
    return -1;
  }
  n_usr = SetEGrid(usr_egrid, log_usr, n, emin, emax, eth);
  return n_usr;
}

int AddRecPW(int n, int step) {
  int i, i2, kl, kl2;
  for (i = pw_scratch.nkl0; i < n+pw_scratch.nkl0; i++) {
    i2 = i*2;
    kl = pw_scratch.kl[i-1] + step;
    kl2 = kl*2;
    pw_scratch.kl[i] = kl;
    pw_scratch.kappa0[i2] = GetKappaFromJL(kl2-1, kl2);
    pw_scratch.kappa0[i2+1] = GetKappaFromJL(kl2+1, kl2);
    if (kl > pw_scratch.max_kl) break;
  }
  pw_scratch.nkl0 = i;
}

int SetRecSpectator(int n_max, int n_frozen, int n_spec) {
  pw_scratch.n_max = n_max;
  pw_scratch.n_frozen = n_frozen;
  pw_scratch.n_spec = n_spec;
}

int SetRecPWOptions(int kl_interp, int max_kl) {
  int k, j, m;

  pw_scratch.kl_interp = kl_interp;
  pw_scratch.max_kl = max_kl;
  pw_scratch.kl[0] = 0;
  pw_scratch.nkl0 = 1;
  pw_scratch.kappa0[0] = 0;
  pw_scratch.kappa0[1] = -1;

  AddRecPW(pw_scratch.kl_interp, 1);
  k = 2;
  j = 2;
  m = pw_scratch.kl[pw_scratch.nkl0-1];
  while (m+k <= max_kl) {
    AddRecPW(j, k);
    m = pw_scratch.kl[pw_scratch.nkl0-1];
    k *= 2;
  }
  pw_scratch.nkl = pw_scratch.nkl0;
}

int SetRecPWLimits(int min, int max) {
  pw_scratch.pw_limits[0] = min;
  pw_scratch.pw_limits[1] = max;
  return 0;
}

int ConstructRecGroupName(char *rgn, char *gn, int n) {
  sprintf(rgn, "+%d__%s", n, gn);
  return 0;
}

int IsRecombinedGroup(int i) {
  char *s;
  int n;
  s = GetGroup(i)->name;
  if (s[0] != '+') return 0;
  n = strtol(s, NULL, 10);
  return n;
}

int RecStates(int n, int k, int *kg) {
  int i, j, m, nsym, nlevels, ncfgs, kg0, t;
  ARRAY *clist;
  CONFIG *rcfg, *c;
  SHELL ns;
  CONFIG_GROUP *g;
  char *gn, rgn[GROUP_NAME_LEN];
  int nm;

  nm = 0;
  for (i = 0; i < k; i++) {
    g = GetGroup(kg[i]);
    clist = &(g->cfg_list);
    for (t = 0; t < g->n_cfgs; t++) {
      c = (CONFIG *) ArrayGet(clist, t);
      if (c->shells[0].n > nm) nm = c->shells[0].n;
    }
  }
  if (n < nm) return 0;

  nm++;

  if (pw_scratch.n_spec > nm) nm = pw_scratch.n_spec;

  if (n >= nm) {
    i = RecStatesFrozen(n, k, kg);
    return 0;
  }

  ns.n = n;
  ns.nq = 1;
  ncfgs = 0;
  for (i = 0; i < k; i++) {
    kg0 = kg[i];
    g = GetGroup(kg0);
    gn = g->name;
    ConstructRecGroupName(rgn, gn, n);
    if ((kg[i] = GroupExists(rgn)) >= 0) continue;
    kg[i] = AddGroup(rgn);
    if (kg[i] < 0) {
      printf("Can not add more Groups\n");
      abort();
    }
    ArrayAppend(rec_complex[n_complex].rg, kg+i);
    clist = &(g->cfg_list);
    for (t = 0; t < g->n_cfgs; t++) {
      c = (CONFIG *) ArrayGet(clist, t);
      for (j = 1; j < 2*pw_scratch.nkl0; j++) {
	if (pw_scratch.kl[j/2] >= n) break;
	ns.kappa = pw_scratch.kappa0[j];
	m = CompareShell(c->shells, &ns);
	if (m > 0) continue;
	if (m == 0) {
	  if (ShellClosed(c->shells)) continue;
	  else {
	    rcfg = malloc(sizeof(CONFIG));
	    rcfg->n_shells = c->n_shells;
	    rcfg->shells = malloc(sizeof(SHELL)*rcfg->n_shells);
	    memcpy(rcfg->shells, c->shells, sizeof(SHELL)*rcfg->n_shells);
	    rcfg->shells[0].nq += 1;
	  }
	} else {
	  rcfg = malloc(sizeof(CONFIG));
	  rcfg->n_shells = c->n_shells + 1;
	  rcfg->shells = malloc(sizeof(SHELL)*rcfg->n_shells);
	  memcpy(rcfg->shells+1, c->shells, sizeof(SHELL)*c->n_shells);
	  memcpy(rcfg->shells, &ns, sizeof(SHELL));
	}
	
	if (Couple(rcfg) < 0) return -3;
	if (AddConfigToList(kg[i], rcfg) < 0) return -4;
	ncfgs++;
      }
    }
  }
  if (ncfgs == 0) return -1;
  nlevels = GetNumLevels();
  nsym = MAX_SYMMETRIES;
  rec_complex[n_complex].n = n;
  rec_complex[n_complex].s0 = nlevels;
  for (i = 0; i < nsym; i++) {
    m = ConstructHamilton(i, k, kg, 0, NULL);
    if (m < 0) continue;
    j = DiagnolizeHamilton();
    if (j < 0) return -1;
    AddToLevels();
  }
  rec_complex[n_complex].s1 = GetNumLevels()-1;
  n_complex++;
  SortLevels(nlevels, -1);
  
  return 0;
}

int RecStatesFrozen(int n, int k, int *kg) {
  int i, j, m, nlevels, nsym, nstates, ko;
  int kl2, j1, j2, p1, p, jmin, jmax, tj;
  HAMILTON *h;
  LEVEL *lev;
  STATE *s;
  SYMMETRY *sym;
  int i0, i1, t, nt;

  nstates = 0;

  i0 = 0;
  if (n_complex == 0) {
    nt = 1;
    i1 = GetNumLevels();
    rec_complex[0].s0 = i1;
  } else {
    nt = n_complex;
  }
  for (t = 0; t < nt; t++) {
    i1 = rec_complex[t].s0;
    for (i = i0; i < i1; i++) {
      lev = GetLevel(i);
      m = lev->major_component;
      sym = GetSymmetry(lev->pj);
      s = (STATE *) ArrayGet(&(sym->states), m);
      if (!InGroups(s->kgroup, k, kg)) {
	continue;
      }
      j1 = lev->pj;
      DecodePJ(j1, &p1, &j1);
      
      for (j = 1; j < 2*pw_scratch.nkl0; j++) {
	kl2 = pw_scratch.kl[j/2];
	p = p1 + kl2;
	if (kl2 >= n) break;
	ko = OrbitalIndex(n, pw_scratch.kappa0[j], 0.0);
	j2 = GetJFromKappa(pw_scratch.kappa0[j]);
	jmin = abs(j2 - j1);
	jmax = j2 + j1;
	for (tj = jmin; tj <= jmax; tj += 2) {
	  AddStateToSymmetry(-(i+1), ko, tj, p, tj);
	  nstates++;
	}
      }
    }  
    i0 = rec_complex[t].s1+1;
  }

  if (nstates == 0) return -1;
  nsym = MAX_SYMMETRIES;
  nlevels = GetNumLevels();
  rec_complex[n_complex].n = n;
  rec_complex[n_complex].s0 = nlevels;
  for (i = 0; i < nsym; i++) {
    if (n >= pw_scratch.n_frozen) {
      i0 = 0;
      for (t = 0; t < nt; t++) {
	i1 = rec_complex[t].s0;
	for (j = i0; j < i1; j++) {
	  lev = GetLevel(j);
	  m = lev->major_component;
	  sym = GetSymmetry(lev->pj);
	  s = (STATE *) ArrayGet(&(sym->states), m);
	  if (!InGroups(s->kgroup, k, kg)) continue;
	  m = ConstructHamiltonFrozen(i, j, NULL, n);
	  if (m < 0) continue;
	  if (DiagnolizeHamilton() < 0) return -2;
	  AddToLevels();
	}
	i0 = rec_complex[t].s1+1;
      }
    } else {
      m = ConstructHamiltonFrozen(i, k, kg, n);
      if (m < 0) continue;
      if (DiagnolizeHamilton() < 0) return -2;
      AddToLevels();
    }
  }
  rec_complex[n_complex].s1 = GetNumLevels()-1;
  n_complex++;
  SortLevels(nlevels, -1);
  return 0;
}

double *RRRadialQkTable(int k0, int k1, int m) {
  int index[3], k;
  double **p, *qk;
  double rint[MAXNE];
  double r0, r1;
  ORBITAL *orb;
  int kappa0, jb0, klb0;
  int kappaf, jf, klf, kf;
  int jfmin, jfmax;
  int ite, ie;
  double eb, aw, e, pref;
  int mode, gauge, c;

  k = 2*abs(m);
  index[0] = k0;
  index[1] = k1;
  if (m >= 0) {
    index[2] = k;
  } else {
    index[2] = k-1;
  }

  p = (double **) MultiSet(qk_array, index, NULL);
  if (*p) {
    return *p;
  }

  orb = GetOrbital(k0);
  kappa0 = orb->kappa;
  GetJLFromKappa(kappa0, &jb0, &klb0);
  if (k1 != k0) {
    orb = GetOrbital(k1);
    if (GetJFromKappa(orb->kappa) != jb0) return NULL;
  }

  gauge = GetTransitionGauge();
  mode = GetTransitionMode();

  *p = (double *) malloc(sizeof(double)*n_tegrid*n_qk);
  qk = *p;

  /* the factor 2 comes from the conitinuum norm */
  pref = 2.0/((k+1.0)*(jb0+1.0));

  for (ite = 0; ite < n_tegrid; ite++) {
    if (mode == M_FR || m == 1 || ite == 0) {
      for (ie = 0; ie < n_egrid; ie++) {
	rint[ie] = 0.0;
      }
      eb = tegrid[ite];
      aw = FINE_STRUCTURE_CONST * eb;
      jfmin = jb0 - k;
      jfmax = jb0 + k;
      for (jf = jfmin; jf <= jfmax; jf += 2) {
	for (c = -1; c <= 1; c += 2) {
	  klf = jf + c;
	  if (jf <= 0 ||
	      klf < 0 ||
	      (m < 0 && IsOdd((klb0+klf+k)/2)) ||
	      (m > 0 && IsEven((klb0+klf+k)/2))) {
	    continue;
	  }
	  kappaf = GetKappaFromJL(jf, klf);
	  for (ie = 0; ie < n_egrid; ie++) {
	    e = egrid[ie];
	    kf = OrbitalIndex(0, kappaf, e);
	    if (mode == M_NR && m != 1) {
	      r0 = MultipoleRadialNR(m, kf, k0, gauge);
	      if (k1 == k0) {
		r1 = r0;
	      } else {
		r1 = MultipoleRadialNR(m, kf, k1, gauge);
	      }
	    } else {
	      r0 = MultipoleRadialFR(aw, m, kf, k0, gauge);
	      if (k1 == k0) {
		r1 = r0;
	      } else {
		r1 = MultipoleRadialFR(aw, m, kf, k1, gauge);
	      }
	    }	  
	    rint[ie] += r0*r1;
	  }  
	}
      }
      c = k-2;
      if (gauge == G_COULOMB && mode == M_FR && m < 0) {
	c -= 2;
      }
      for (ie = 0; ie < n_egrid; ie++) {
	e = eb + egrid[ie];
	rint[ie] *= e*pref;
	if (c) {
	  aw = FINE_STRUCTURE_CONST*e;
	  rint[ie] *= pow(aw, c);
	}
      }
    }
    if (interpolate_egrid) {
      SVDFit(NPARAMS, qk, NULL, 1E-3, n_egrid, xegrid[ite], log_xegrid[ite],
	     rint, sigma, RRRadialQkBasis);
    } else {
      for (ie = 0; ie < n_qk; ie++) {
	qk[ie] = rint[ie];
      }
    }
    qk += n_qk;
  }

  return *p;
}
  
int RRRadialQk(double *rqc, double te, int k0, int k1, int m) {
  int i, np, nd;
  int j, k;
  double *rqe, rq[MAXNTE];
 
  rqe = RRRadialQkTable(k0, k1, m);
  if (rqe == NULL) return -1;

  if (n_tegrid == 1) {
    for (i = 0; i < n_qk; i++) {
      rqc[i] = rqe[i];
    }
  } else {
    np = 3;
    nd = 1;
    for (i = 0; i < n_qk; i++) {
      j = i;
      for (k = 0; k < n_tegrid; k++) {
	rq[k] = rqe[j];
	j += n_qk;
      }
      uvip3p_(&np, &n_tegrid, tegrid, rq, &nd, &te, &rqc[i]);
    }
  }
  
  return 0;
}

void RRRadialQkBasis(int npar, double *yb, double x, double logx) {
  int i;
  double t;

  t = 1.0/x;
  yb[0] = t*t;
  t = sqrt(t);
  for (i = 1; i < npar; i++) {
    yb[i] = yb[i-1]*t;
  }
}

int BoundFreeOS(double *rqc, double *eb, int rec, int f, int m) {
  LEVEL *lev1, *lev2;
  ANGULAR_ZFB *ang;
  ORBITAL *orb1;
  int nz, ie, k;
  double rq[NPARAMS], a;
  int j1, j2;
  int i, j;

  lev1 = GetLevel(rec);
  lev2 = GetLevel(f);
  j1 = lev1->pj;
  j2 = lev2->pj;
  DecodePJ(j1, NULL, &j1);
  DecodePJ(j2, NULL, &j2);

  *eb = (lev2->energy - lev1->energy);

  nz = AngularZFreeBound(&ang, f, rec);
  if (nz <= 0) return -1;

  for (ie = 0; ie < n_qk; ie++) rqc[ie] = 0.0;
  for (i = 0; i < nz; i++) {
    for (j = 0; j <= i; j++) {
      k = RRRadialQk(rq, *eb, ang[i].kb, ang[j].kb, m);
      if (k < 0) continue;
      a = ang[i].coeff*ang[j].coeff;    
      if (j != i) {
	a *= 2;
      } 
      for (ie = 0; ie < n_qk; ie++) {	
	rqc[ie] += a * rq[ie];
      }
    }
  }

  free(ang);

  return 0;
}

int AIRate(double *rate, double *e, int rec, int f) {  
  LEVEL *lev1, *lev2;
  ANGULAR_ZxZMIX *ang;
  ANGULAR_ZFB *zfb;
  STATE *st;
  int k, nz, nzfb, ik, i, j1, j2, ij, kappaf, ip;
  int jf, k0, k1, kb, njf, nkappaf, klf, jmin, jmax;
  double *p, y2[MAXNE], r, s, log_e;
  double *ai_pk, ai_pk0[MAXNE];
  int np, nt;

  *rate = 0.0;
  lev1 = GetLevel(rec);
  lev2 = GetLevel(f);
  
  *e = lev1->energy - lev2->energy;
  if (*e <= 0.0) return -1;
  log_e = log(*e);

  i = lev1->major_component;

  j1 = lev1->pj;
  j2 = lev2->pj;

  st = (STATE *) ArrayGet(&(GetSymmetry(j1)->states), i);
  if (st->kgroup < 0) {
    k = GetOrbital(st->kcfg)->kappa;
  } else {
    k = (GetConfig(st)->shells[0]).kappa;
  }
  k = GetLFromKappa(k);
  k = k/2;
  if (k < pw_scratch.pw_limits[0] || k > pw_scratch.pw_limits[1]) return -1;

  DecodePJ(j1, NULL, &j1);
  DecodePJ(j2, NULL, &j2);

  jmin = abs(j1-j2);
  jmax = j1+j2;
  njf = (jmax - jmin)/2 + 1;
  nkappaf = njf*2;
  p = malloc(sizeof(double)*nkappaf);
  for (ip = 0; ip < nkappaf; ip++) p[ip] = 0.0;

  nz = AngularZxZFreeBound(&ang, f, rec);
  np = 3;
  nt = 1;
  if (nz > 0) {
    for (i = 0; i < nz; i++) {
      jf = ang[i].k0;
      kb = ang[i].k1;
      k0 = ang[i].k2;
      k1 = ang[i].k3;
      ij = (jf - jmin);
      for (ik = -1; ik <= 1; ik += 2) {
	klf = jf + ik;  
	kappaf = GetKappaFromJL(jf, klf);
	AIRadialPk(&ai_pk, k0, k1, kb, kappaf, ang[i].k);
	if (n_egrid > 1) {
	  uvip3p_(&np, &n_egrid, log_egrid, ai_pk, &nt, &log_e, &s);
	} else {
	  s = ai_pk[0];
	}
	ip = (ik == -1)? ij:(ij+1);
	p[ip] += s*ang[i].coeff;
      }
    }    
    free(ang);
  }
  
  nzfb = AngularZFreeBound(&zfb, f, rec);
  if (nzfb > 0) {
    for (i = 0; i < nzfb; i++) {
      kb = zfb[i].kb;
      jf = GetOrbital(kb)->kappa;
      jf = GetJFromKappa(jf);
      ij = jf - jmin;
      for (ik = -1; ik <= 1; ik += 2) {
	klf = jf + ik;
	kappaf = GetKappaFromJL(jf, klf);
	AIRadial1E(ai_pk0, kb, kappaf);
	if (n_egrid > 1) {
	  uvip3p_(&np, &n_egrid, log_egrid, ai_pk0, &nt, &log_e, &s);
	} else {
	  s = ai_pk0[0];
	}
	ip = (ik == -1)?ij:(ij+1);
	if (j2 > j1) {
	  if (IsEven(ij/2)) s = -s;
	} else {
	  if (IsOdd(ij/2)) s = -s;
	}
	p[ip] += s*zfb[i].coeff;
      }
    }
    free(zfb);
  }

  if (nz <= 0 && nzfb <= 0) return -1;

  r = 0.0;
  for (i = 0; i < nkappaf; i++) {
    r += p[i]*p[i];
  }
  /* the prefactor 4.0 includes the factor 2 from the continuum norm,
     otherwize, it should have been 2.0 */
  *rate = 4.0*r/(j1+1.0);

  free(p);

  return 0;
}

int AIRadial1E(double *ai_pk, int kb, int kappaf) {
  int kf;
  int i;

  for (i = 0; i < n_egrid; i++) {
    kf = OrbitalIndex(0, kappaf, egrid[i]);
    ResidualPotential(ai_pk+i, kf, kb);
  }
  return 0;
}  

int AIRadialPk(double **ai_pk, int k0, int k1, int kb, int kappaf, int k) {
  int i, j, kf;
  int ks[4];
  double e, sd, se;
  double **p;
  int index[5];

  if (kappaf > 0) {
    index[0] = 2*kappaf-1;
  } else {
    index[0] = -2*kappaf-2;
  }
  index[1] = kb;
  index[2] = k0;
  index[3] = k1;
  index[4] = k/2;

  p = (double **) MultiSet(pk_array, index, NULL);
  if (*p) {
    *ai_pk = *p;
    return 0;
  } 
 
  (*p) = (double *) malloc(sizeof(double)*n_egrid);
  *ai_pk = *p;
  for (i = 0; i < n_egrid; i++) {
    e = egrid[i];
    kf = OrbitalIndex(0, kappaf, e);
    ks[0] = k0;
    ks[1] = kf;
    ks[2] = k1;
    ks[3] = kb;
    SlaterTotal(&sd, &se, NULL, ks, k, 0);
    (*ai_pk)[i] = sd+se;
  }

  return 0;
}

int SaveRecRR(int nlow, int *low, int nup, int *up, 
	      char *fn, int m) {
  int i, j, k, n, ie, jp;
  char t;
  FILE *f;
  double rqc[NPARAMS], rqb[NPARAMS];
  double phi, rr, ee, eph, eb;
  double s, trr[MAXNUSR], tpi[MAXNUSR];
  int j1, j2, nlow0;
  LEVEL *lev1, *lev2;
  double e, emin, emax;
  double awmin, awmax;

  f = fopen(fn, "w");

  if (!f) return -1;

  emin = 1E10;
  emax = 1E-10;
  k = 0;
  for (i = 0; i < nlow; i++) {
    lev1 = GetLevel(low[i]);
    for (j = 0; j < nup; j++) {
      lev2 = GetLevel(up[j]);
      e = lev2->energy - lev1->energy;
      if (e > 0) k++;
      if (e < emin && e > 0) emin = e;
      if (e > emax) emax = e;
    }
  }
  if (k == 0) {
    printf("No recombination can occur\n");
    return 0;
  }
  
  e = 2.0*(emax - emin)/(emin+emax);
  if (tegrid[0] < 0.0) {
    if (e < 0.1) {
      SetRRTEGrid(1, 0.5*(emin+emax), emax);
    } else if (e < 0.5) {
      SetRRTEGrid(2, emin, emax);
    } else {
      if (k == 2) n_tegrid = 2;
      SetRRTEGrid(n_tegrid, emin, emax);
    }
  }
  
  if (m == 1 || GetTransitionMode() == M_FR) {
    FreeMultipoleArray();
    awmin = emin * FINE_STRUCTURE_CONST;
    awmax = emax * FINE_STRUCTURE_CONST;
    if (e < 0.3) {
      SetAWGrid(1, 0.5*(awmin+awmax), awmax);
    } else if (e < 1.0) {
      SetAWGrid(2, awmin, awmax);
    } else {
      SetAWGrid(3, awmin, awmax);
    }
  }  

  e = 0.5*(emin + emax);
  emin = 0.02*e;
  emax = 8.0*e;
  egrid_type = 1;
  if (usr_egrid_type < 0) usr_egrid_type = 1;
  interpolate_egrid = 1;

  if (n_usr > 0 && usr_egrid[0] < 0.0) {
    SetUsrPEGrid(n_usr, emin, emax, e);
  }  
  
  if (n_egrid == 0) {
    n_egrid = 6;
  }
  if (egrid[0] < 0.0) {
    if (n_usr > 0 && n_usr <= NPARAMS && egrid_type == usr_egrid_type) {
      SetPEGridDetail(n_usr, usr_egrid);
      interpolate_egrid = 0;
    } else {
      SetPEGrid(n_egrid, emin, emax, e);
    }
  }

  if (interpolate_egrid) {
    for (j = 0; j < n_egrid; j++) {
      for (i = 0; i < n_tegrid; i++) {
	xegrid[i][j] = egrid[j]/tegrid[i];
	if (egrid_type == 1) xegrid[i][j] += 1.0;
	log_xegrid[i][j] = log(xegrid[i][j]);
      }
      sigma[j] = 1.0/sqrt(xegrid[0][j]);
    }
    n_qk = NPARAMS;
  } else {
    n_qk = n_egrid;
  }

  if (m < 0) t = 'E';
  else t = 'M';
  fprintf(f, " Multipole %c%d\n\n", t, abs(m));

  fprintf(f,  " IEGRID:\n ");
  for (i = 0; i < n_tegrid; i++) {
    fprintf(f, "%10.4E ", tegrid[i]*HARTREE_EV);
  }
  fprintf(f, "\n");

  fprintf(f, " PEGRID:\n ");
  for (i = 0; i < n_egrid; i++) {
    fprintf(f, "%10.4E ", egrid[i]*HARTREE_EV);
  }
  fprintf(f, "\n\n");
  
  if (interpolate_egrid) {
    fprintf(f, " Strength = A/u^2 + B/u^2.5 + C/u^3 + D/u^3.5 + E/u^4\n");
  }
  if (output_format >= 0 && n_usr > 0) {
    if (usr_egrid_type == 0) fprintf(f, " Photon UsrEGrid\n");
    else fprintf(f, " Electron UsrEGrid\n");
  }
  fprintf(f, "\n");

  fprintf(f, "Free  2J   Bound 2J   E_Binding\n\n");
  for (i = 0; i < nup; i++) {
    j1 = LevelTotalJ(up[i]);    
    for (ie = 0; ie < n_usr; ie++) {
      trr[ie] = 0.0;
      tpi[ie] = 0.0;
    }
    nlow0 = 0;
    for (j = 0; j < nlow; j++) {
      j2 = LevelTotalJ(low[j]);
      k = BoundFreeOS(rqc, &eb, low[j], up[i], m);
      if (k < 0) continue;
      nlow0++;
      fprintf(f, "%-5d %-2d   %-5d %-2d   %10.4E\n",
	      up[i], j1, low[j], j2, eb*HARTREE_EV);
      if (interpolate_egrid) {
	if (output_format >= 0) {
	  for (ie = 0; ie < n_usr; ie++) {
	    xusr[ie] = usr_egrid[ie]/eb;
	    if (usr_egrid_type == 1) xusr[ie] += 1.0;
	    log_xusr[ie] = log(xusr[ie]);
	  }
	}
	for (ie = 0; ie < NPARAMS; ie++) {
	  fprintf(f, "%11.4E  ", rqc[ie]);
	}
	fprintf(f, "\n\n");
      }
      if (output_format >= 0 && n_usr > 0) {
	for (ie = 0; ie < n_usr; ie++) {
	  if (usr_egrid_type == 0) {
	    eph = usr_egrid[ie];
	    ee = eph - eb;
	  } else {
	    ee = usr_egrid[ie];
	    eph = eb + ee;
	  }
	  if (interpolate_egrid) {
	    RRRadialQkBasis(NPARAMS, rqb, xusr[ie], log_xusr[ie]);
	    s = 0.0;
	    for (jp = 0; jp < NPARAMS; jp++) s += rqc[jp]*rqb[jp];
	  } else {
	    s = rqc[ie];
	  }
	  phi = 2.0*PI*FINE_STRUCTURE_CONST*s*AREA_AU20;
	  rr = phi * pow(FINE_STRUCTURE_CONST*eph, 2) / (2.0*ee);
	  trr[ie] += rr;
	  tpi[ie] += phi;
	  fprintf(f, "%-10.3E %-10.3E %-10.3E %-10.3E\n", 
		usr_egrid[ie]*HARTREE_EV, rr, phi, s);
	}
	fprintf(f, "\n");
      }
    }
    if (output_format >= 0 && n_usr > 0 && nlow0 > 1) {
      fprintf(f, "==================================\n");
      fprintf(f, "%-5d   \tTotal\n", up[i]);
      for (ie = 0; ie < n_usr; ie++) {
	fprintf(f, "%-10.3E %-10.3E %-10.3E\n", 
		usr_egrid[ie]*HARTREE_EV, trr[ie], tpi[ie]);
      }
    }
    fprintf(f, "\n");
  }
  
  fclose(f);
  
  return 0;
}

int SaveDR(int nf, int *f, int na, int *a, int nb, int *b, int ng, int *g, 
	   char *fna, char *fnt, int channel) {
  int i, j, k, n, kb, j1, j2, m, do_transition;
  LEVEL *lev1, *lev2;
#ifdef PERFORM_STATISTICS
  clock_t start, stop;
  clock_t t_ai, t_rd, tt_ai, tt_rd;
  STRUCT_TIMING st_start, st_stop;
  ANGULAR_TIMING angt;
  RECOUPLE_TIMING recouplet;
  RAD_TIMING radt;
#endif
  double emin, emax;
  double *ea, *sa, tai, sd, e0; 
  double *et, *sr, *rd, elow, trd, trd1, tr_cut;
  char t;
  FILE *fa, *ft;

  tr_cut = GetTransitionCut();

  emin = 1E10;
  emax = 1E-16;
  k = 0;
  for (i = 0; i < na; i++) {
    lev1 = GetLevel(a[i]);
    for (j = 0; j < nf; j++) {
      lev2 = GetLevel(f[j]);
      e0 = lev1->energy - lev2->energy;
      if (e0 > 0) k++;
      if (e0 < emin && e0 > 0) emin = e0;
      if (e0 > emax) emax = e0;
    }
  }
  if (k == 0) {
    printf("No recombination channels\n");
    return 0;
  }
  
  if (n_egrid == 0) {
    n_egrid = 3;
  }
  if (egrid[0] < 0.0) {
    e0 = 2.0*(emax-emin)/(emax+emin);
    if (e0 < 0.1) {
      SetPEGrid(1, 0.5*(emin+emax), emax, 0.0);
    } else if (e0 < 0.5) {
      SetPEGrid(2, emin, emax, 0.0);
    } else {
      if (k == 2) n_egrid = 2;
      SetPEGrid(n_egrid, emin, emax, 0.0);
    }
  }

  if (GetTransitionMode() == M_FR) {
    emin = 1E10;
    emax = 1E-16;
    k = 0;
    for (i = 0; i < na; i++) {
      lev1 = GetLevel(a[i]);
      for (j = 0; j < nb; j++) {
	lev2 = GetLevel(b[j]);
	e0 = lev1->energy - lev2->energy;
	if (e0 > 0) k++;
	if (e0 < emin && e0 > 0) emin = e0;
	if (e0 > emax) emax = e0;
      }
    }
    if (k == 0) {
      printf("No decay routes\n");
      return 0;
    }    

    emin *= FINE_STRUCTURE_CONST;
    emax *= FINE_STRUCTURE_CONST;
    e0 = 2.0*(emax-emin)/(emin+emax);
    
    FreeMultipoleArray();
    if (e0 < 0.1) {
      SetAWGrid(1, 0.5*(emin+emax), emax);
    } else if (e0 < 1.0) {
      SetAWGrid(2, emin, emax);
    } else {
      SetAWGrid(3, emin, emax);
    }
  }
    
  if (nf <= 0 || na <= 0 || nb <= 0 || ng <= 0) return -1;

  fa = fopen(fna, "w");
  if (!fa) return -1;
  ft = fopen(fnt, "w");
  if (!ft) return -1;

  fprintf(fa, "DR Channel: %-4d  EGRID: ", channel);
  for (i = 0; i < n_egrid; i++) {
    fprintf(fa, "%10.4E ", egrid[i]*HARTREE_EV);
  }
  fprintf(fa, "\n\n");

  fprintf(fa, 
	  "Bound  2J\tFree   2J\tEe         AI(AU)    CS(10^-20)\n");

  
  fprintf(ft, "up     2J\tlow     2J\tDelta_E    M  gf        A(AU)\n");

  sa = malloc(sizeof(double)*nf);
  ea = malloc(sizeof(double)*nf);
  e0 = GetLevel(0)->energy;
  rd = malloc(sizeof(double)*nb);
  sr = malloc(sizeof(double)*nb);
  et = malloc(sizeof(double)*nb);

#ifdef PERFORM_STATISTICS
  start = clock();
  tt_ai = 0;
  tt_rd = 0;
#endif

  for (i = 0; i < na; i++) {
#ifdef PERFOR_STATISTICS
    GetStructTiming(&st_start);
#endif
    j1 = LevelTotalJ(a[i]);
    tai = 0.0;
    for (j = 0; j < nf; j++) {
      k = AIRate(sa+j, ea+j, a[i], f[j]);
      if (k < 0) continue;
      tai += sa[j];
    }
    if (tai < 1E-30) continue;
    do_transition = 0;
    for (k = 0; k < ng; k++) {
      for (j = 0; j < nf; j++) {
	if (g[k] == f[j]) {
	  if (sa[j] > ai_cut*tai) do_transition = 1;
	}
      }
    }
    for (j = 0; j < nf; j++) {
      if (sa[j] < ai_cut*tai) continue;
      j2  = LevelTotalJ(f[j]);
      sd = 0.5*(j1+1.0)*PI*PI*sa[j]/(ea[j]*(j2+1.0));
      sd *= AREA_AU20*HARTREE_EV;    
      ea[j] *= HARTREE_EV;
      fprintf(fa, "%-6d %-2d\t%-6d %-2d\t%10.4E %9.3E %9.3E\n",
	      a[i], j1, f[j], j2, ea[j], sa[j], sd);
    }

#ifdef PERFORM_STATISTICS
    stop = clock();
    t_ai = stop-start;
    start = stop;
    tt_ai += t_ai;
#endif

    fprintf(fa, "%-6d   \tTotal   \t           %9.3E\n\n", a[i], tai);
    
    if (!do_transition) continue;

    trd = 0.0;
    trd1 = 0.0;
    for (j = 0; j < nb; j++) {
      rd[j] = 0.0;
      elow = GetLevel(b[j])->energy;
      m = -1;
      et[j] = 0.0;
      k = OscillatorStrength(sr+j, et+j, m, b[j], a[i]);
      if (k != 0) continue;
      if (sr[j] < 1E-30) continue;
      rd[j] = 2*pow((FINE_STRUCTURE_CONST*et[j]),2)*FINE_STRUCTURE_CONST;
      rd[j] *= sr[j]/(j1+1.0);
      trd += rd[j];
      if (elow < e0) trd1 += rd[j];
    }
    if (trd < 1E-30) continue;
    for (j = 0; j < nb; j++) {
      if (rd[j] < (tr_cut*trd)) continue;
      et[j] *= HARTREE_EV;      
      if (m < 0) t = 'E';
      else t = 'M';
      elow = GetLevel(b[j])->energy;
      if (elow < e0) k = -b[j];
      else k = b[j];
      j2 = LevelTotalJ(b[j]);
      fprintf(ft, "%-6d %-2d\t%-7d %-2d\t%10.4E %c%d %9.3E %9.3E\n", 
	      a[i], j1, k, j2, et[j], t, abs(m), sr[j], rd[j]);
    }	
   
    fprintf(ft, 
	    "%-6d   \tTotal     \t              %9.3E %9.3E\n\n", 
	    a[i], trd, trd1);

#ifdef PERFORM_STATISTICS
    stop = clock();
    t_rd = stop-start;
    start = stop;
    tt_rd += t_rd;
#endif

#ifdef PERFORM_STATISTICS
  GetStructTiming(&st_stop);
  
  fprintf(ft, "Time in AI: %6.1E, RD: %6.1E\n", 
	  ((double) tt_ai)/CLOCKS_PER_SEC, ((double) tt_rd)/CLOCKS_PER_SEC);
  
  fprintf(ft, "AngZMix: %6.1E, AngZFB: %6.1E, AngZxZFB: %6.1E, SetH: %6.1E DiagH: %6.1E\n",
	  ((double) (st_stop.angz_mix))/CLOCKS_PER_SEC,
	  ((double) (st_stop.angz_fb))/CLOCKS_PER_SEC,
	  ((double) (st_stop.angzxz_fb))/CLOCKS_PER_SEC,
	  ((double) (st_stop.set_ham))/CLOCKS_PER_SEC,
	  ((double) (st_stop.diag_ham))/CLOCKS_PER_SEC);
  fprintf(ft, "AngZS: %6.1E, AngZFBS: %6.1E, AngZxZFBS: %6.1E, AddZ: %6.1E, AddZxZ: %6.1E\n",
	  ((double) (st_stop.angz_states))/CLOCKS_PER_SEC,
	  ((double) (st_stop.angzfb_states))/CLOCKS_PER_SEC,
	  ((double) (st_stop.angzxzfb_states))/CLOCKS_PER_SEC,
	  ((double) (st_stop.add_angz))/CLOCKS_PER_SEC,
	  ((double) (st_stop.add_angzxz))/CLOCKS_PER_SEC);

  GetAngularTiming(&angt);
  fprintf(ft, "W3J: %6.1E, W6J: %6.1E, W9J: %6.1E\n", 
	  ((double)angt.w3j)/CLOCKS_PER_SEC, 
	  ((double)angt.w6j)/CLOCKS_PER_SEC, 
	  ((double)angt.w9j)/CLOCKS_PER_SEC);
  GetRecoupleTiming(&recouplet);
  fprintf(ft, "AngZ: %6.1E, AngZxZ: %6.1E, Interact: %6.1E\n",
	  ((double)recouplet.angz)/CLOCKS_PER_SEC,
	  ((double)recouplet.angzxz)/CLOCKS_PER_SEC,
	  ((double)recouplet.interact)/CLOCKS_PER_SEC);
  GetRadTiming(&radt);
  fprintf(ft, "Dirac: %6.1E, 1E: %6.1E, 2E: %6.1E\n", 
	  ((double)radt.dirac)/CLOCKS_PER_SEC, 
	  ((double)radt.radial_1e)/CLOCKS_PER_SEC,
	  ((double)radt.radial_2e)/CLOCKS_PER_SEC);
  fprintf(ft, "\n");
#endif /* PERFORM_STATISTICS */

  }
  
  FreeRecPk();
  FreeSlaterArray();
  FreeMultipoleArray();

  free(sa);
  free(ea);
  free(sr);
  free(rd);
  free(et);
  fclose(fa);
  fclose(ft);
  return 0;
}

      
int SaveAI(int nlow, int *low, int nup, int *up, char *fn, int channel) {
  int i, j, k, n, kb, j1, j2;
  LEVEL *lev1, *lev2;
  double emin, emax;
  double *e, *s, tai, a, sdr;
  FILE *f;


  emin = 1E10;
  emax = 1E-10;
  k = 0;
  for (i = 0; i < nlow; i++) {
    lev1 = GetLevel(low[i]);
    for (j = 0; j < nup; j++) {
      lev2 = GetLevel(up[j]);
      a = lev1->energy - lev2->energy;
      if (a > 0) k++;
      if (a < emin && a > 0) emin = a;
      if (a > emax) emax = a;
    }
  }
  if (k == 0) {
    printf("No Recombination channels\n");
    return -1;
  }
  if (n_egrid == 0) {
    n_egrid = 3;
  }

  if (egrid[0] < 0.0) {
    a = 2.0*(emax-emin)/(emax+emin);
    if (a < 0.1) {
      SetPEGrid(1, 0.5*(emin+emax), emax, 0.0);
    } else if (a < 0.5) {
      SetPEGrid(2, emin, emax, 0.0);
    } else {
      if (k == 2) n_egrid = 2;
      SetPEGrid(n_egrid, emin, emax, 0.0);
    }
  }
 
  f = fopen(fn, "w");
  if (!f) return -1;
  fprintf(f, "DR Channel: %-4d  EGRID: ", channel);
  for (i = 0; i < n_egrid; i++) {
    fprintf(f, "%10.4E ", egrid[i]*HARTREE_EV);
  }
  fprintf(f, "\n\n");

  fprintf(f, "Bound  2J\tFree   2J\tEe         AI(AU)    CS(10^-20)\n");

  if (nup <= 0 || nlow <= 0) return -1;
  s = malloc(sizeof(double)*nup);
  e = malloc(sizeof(double)*nup);
  for (i = 0; i < nlow; i++) {
    j1 = LevelTotalJ(low[i]);
    tai = 0.0;
    for (j = 0; j < nup; j++) {
      k = AIRate(s+j, e+j, low[i], up[j]);
      if (k < 0) continue;
      tai += s[j];
    }
    if (tai < 1E-30) continue;
    for (j = 0; j < nup; j++) {
      if (s[j] < ai_cut*tai) continue;
      j2  = LevelTotalJ(up[j]);
      sdr = 0.5*(j1+1.0)*PI*PI*s[j]/(e[j]*(j2+1.0)) * AREA_AU20*HARTREE_EV;
      e[j] *= HARTREE_EV;
      fprintf(f, "%-6d %-2d\t%-6d %-2d\t%10.4E %9.3E %9.3E\n",
	      low[i], j1, up[j], j2, e[j], s[j], sdr);
    }
    fprintf(f, "%-6d   \tTotal    \t           %9.3E\n\n", low[i], tai);
  }

  FreeRecPk();
  FreeSlaterArray();
  FreeMultipoleArray();

  free(s);
  free(e);
  fclose(f);
  return 0;
}

int DROpen(int n, int *nlev, int **ops) {
  int i, j, n0, old_n;
  LEVEL *lev;
  double e0, e, z;

  e0 = GetLevel(0)->energy;
  z = GetResidualZ();
  z = z*z/2.0;

  (*ops) = malloc(sizeof(int)*n);
  j = 0;
  old_n = -1;
  for (i = 0; i < n; i++) {
    lev = GetLevel(nlev[i]);
    e = lev->energy - e0;
    if (e < EPS16) continue;
    n0 = ceil(sqrt(z/e));
    if (n0 != old_n) {
      (*ops)[j] = n0;
      j++;
      old_n = n0;
    }
  }

  return j;
}

static void _FreeRecPk(void *p) {
  double *dp;
  dp = *((double **) p);
  free(dp);
}

int FreeRecPk() {
  if (pk_array->array == NULL) return 0;
  MultiFreeData(pk_array->array, pk_array->ndim, _FreeRecPk);
  return 0;
}

int FreeRecQk() {
  if (qk_array->array == NULL) return 0;
  MultiFreeData(qk_array->array, qk_array->ndim, _FreeRecPk);
  return 0;
}

int FreeRecAngZ() {
  int i, j, *k;
  for (i = 0; i < n_complex; i++) {
    for (j = 0; j < (rec_complex[i].rg)->dim; j++) {
      k = (int *)ArrayGet(rec_complex[i].rg, j);  
      if (k) FreeAngZ(*k, -1);
    }
  }
}

int InitRecombination() {
  int blocks[5] = {5, 10, 5, 10, 5};
  int ndim;
  int i;
  
  ndim = 5;
  pk_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(pk_array, sizeof(double *), ndim, blocks);
  
  ndim = 3;
  qk_array = (MULTI *) malloc(sizeof(MULTI));
  blocks[0] = 10;
  blocks[1] = 10;
  blocks[2] = 4;
  MultiInit(qk_array, sizeof(double *), ndim, blocks);  

  n_complex = 0;
  for (i = 0; i < MAX_COMPLEX; i++) {
    rec_complex[i].rg = (ARRAY *) malloc(sizeof(ARRAY));
    ArrayInit(rec_complex[i].rg, sizeof(int), 64);
  }
  n_egrid = 0;
  egrid[0] = -1.0;
  n_usr = 0;
  usr_egrid[0] = -1.0;
  
  n_tegrid = 0;
  tegrid[0] = -1.0;

  output_format = 0;

  SetRecPWOptions(12, 12);

}







