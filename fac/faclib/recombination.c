#include "recombination.h"
#include "time.h"
#include "cf77.h"

static char *rcsid="$Id: recombination.c,v 1.73 2003/12/05 06:24:51 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

static int qk_mode;
static double qk_fit_tolerance;

static int egrid_type = -1;
static int usr_egrid_type = -1;

static int n_egrid = 0;
static double egrid[MAXNE];
static double log_egrid[MAXNE];
static double xegrid[MAXNE];
static double log_xegrid[MAXNE];
static double egrid_min;
static double egrid_max;
static int egrid_limits_type = 0;

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

#define MAXAIM 1024
#define NPARAMS 3
static ARRAY *hyd_qk_array;

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
} pw_scratch = {RECNSPEC, RECNFROZEN, 
		RECNMAX, RECLMAX, RECLMAX,
		0, 0, {0, RECLMAX}};

double ai_cut = AICUT;

static REC_COMPLEX rec_complex[MAX_COMPLEX];
int n_complex = 0;

int SetAICut(double c) {
  ai_cut = c;
  return 0;
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

int SetPEGridLimits(double min, double max, int type) {
  if (min <= 0) egrid_min = 0.05;
  else egrid_min = min;
  if (max <= 0) egrid_max = 8.0;
  else egrid_max = max;
  egrid_limits_type = type;
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
  return 0;
}

int SetRecSpectator(int n_max, int n_frozen, int n_spec) {
  if (n_max > 0) pw_scratch.n_max = n_max;
  if (n_frozen > 0) pw_scratch.n_frozen = n_frozen;
  if (n_spec > 0) pw_scratch.n_spec = n_spec;
  return 0;
}

int SetRecQkMode(int m, double tol) {
  if (m == QK_DEFAULT) qk_mode = QK_FIT;
  else qk_mode = m;
  if (tol > 0.0) qk_fit_tolerance = tol;
  return 0;
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
  return 0;
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

int RecStates(int n, int k, int *kg, char *fn) {
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
    i = RecStatesFrozen(n, k, kg, fn);
    return i;
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
      exit(1);
    }
    ArrayAppend(rec_complex[n_complex].rg, kg+i, NULL);
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
    AddToLevels(0, NULL);
  }
  rec_complex[n_complex].s1 = GetNumLevels()-1;
  n_complex++;
  SortLevels(nlevels, -1);
  SaveLevels(fn, nlevels, -1);

  return 0;
}

int RecStatesFrozen(int n, int k, int *kg, char *fn) {
  int i, j, m, nlevels, nsym, nstates, ko;
  int kl2, j1, j2, p1, p, jmin, jmax, tj;
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
      m = lev->pb;
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
	  m = lev->pb;
	  sym = GetSymmetry(lev->pj);
	  s = (STATE *) ArrayGet(&(sym->states), m);
	  if (!InGroups(s->kgroup, k, kg)) continue;
	  m = ConstructHamiltonFrozen(i, j, NULL, n);
	  if (m < 0) continue;
	  if (DiagnolizeHamilton() < 0) return -2;
	  AddToLevels(0, NULL);
	}
	i0 = rec_complex[t].s1+1;
      }
    } else {
      m = ConstructHamiltonFrozen(i, k, kg, n);
      if (m < 0) continue;
      if (DiagnolizeHamilton() < 0) return -2;
      AddToLevels(0, NULL);
    }
  }
  rec_complex[n_complex].s1 = GetNumLevels()-1;
  n_complex++;
  SortLevels(nlevels, -1);
  SaveLevels(fn, nlevels, -1);

  return 0;
}

void RRRadialQkHydrogenicParams(int np, double *p, 
				double z0, int n, int kl) {
#define NNE 12
  double **qk;
  double *t;
  double z, am, d, eth;
  double pc[200*NNE];
  double pnc[NNE], e[NNE];
  static double x[NNE], logx[NNE];
  int ipcp;
  static int iopt = 2;
  int ipvt[NPARAMS];
  int ierr;
  int lwa=5*NPARAMS+NNE;
  double wa[5*NPARAMS+NNE];
  double fvec[NNE], fjac[NNE*NPARAMS];
  double tol;
  int i, j, k, ne;

  qk = (double **) ArraySet(hyd_qk_array, n, NULL, InitPointerData);
  if (*qk == NULL) {
    tol = 1E-4;
    *qk = (double *) malloc(sizeof(double)*n*np);
    ipcp = 0;
    z = 1.0;
    am = 100.0;
    ne = NNE;
    if (iopt == 2) {
      x[0] = 1.01;
      logx[0] = log(x[0]);
      logx[ne-1] = log(80.0);
      d = (logx[ne-1] - logx[0])/(ne-1.0);
      for (i = 1; i < ne; i++) {
	logx[i] = logx[i-1] + d;
	x[i] = exp(logx[i]);
      }
      i = 200;
      PIXZ1(z, am, i, n, e, pc, NULL, pnc, 
	    ne, n, iopt, ipcp);
      iopt = 0;
    }
    
    eth = 0.5/(n*n);
    for (i = 0; i < ne; i++) {
      e[i] = (x[i]-1.0)*eth;
    }
    PIXZ1(z, am, ne, n, e, pc, NULL, pnc, 
	  ne, n, iopt, ipcp);
    t = *qk;    
    for (i = 0; i < n; i++) {
      k = i;
      for (j = 0; j < ne; j++) {
	pnc[j] = pc[k]/x[j];
	k += n;
      }
      t[0] = pnc[0];
      t[1] = 3.0*(i+1);
      t[2] = 1.0;
      ierr = NLSQFit(np, t, tol, ipvt, fvec, fjac, ne, wa, lwa, 
		     ne, x, logx, pnc, pnc, RRRadialQkFromFit, &i);
      t += np;
    }
  }

  t = (*qk) + kl*np;
  p[0] = t[0]/(z0*z0);
  for (i = 1; i < np; i++) {
    p[i] = t[i];
  }
#undef NNE
}  

void RRRadialQkFromFit(int np, double *p, int n, double *x, double *logx, 
		       double *y, double *dy, int ndy, void *extra) {
  int kl, i, k;
  double a, b, c, d, e, f, g, h;

  kl = *((int *) extra);
  if (ndy <= 0) {
    for (i = 0; i < n; i++) {
      a = 4.5 + kl;
      b = 0.5*p[1];
      c = sqrt(x[i]);
      d = p[2] + c;
      e = p[2] + 1;
      f = logx[i]*(b-a);
      h = log(e/d);
      g = h*p[1];
      y[i] = p[0] * exp(f+g);
    }
  } else {
    for (i = 0; i < n; i++) {
      k = i;
      a = 4.5 + kl;
      b = 0.5*p[1];
      c = sqrt(x[i]);
      d = p[2] + c;
      e = p[2] + 1;
      f = logx[i]*(b-a);
      h = log(e/d);
      g = h*p[1];
      a = exp(f+g);
      dy[k] = a;
      k += ndy;
      a = p[0]*a;
      dy[k] = 0.5*a*logx[i] + a*h;
      if (np == 3) {
	k += ndy;
	dy[k] = a*p[1]*(1.0/e - 1.0/d);
      }
    }
  }
}

double PICrossH(double z, int n0, int kl0, double e, int os) {
  double hp[NPARAMS], eth;
  double x, logx, r;

  eth = 0.5*z*z/(n0*n0);
  if (e < eth) return 0.0;
  RRRadialQkHydrogenicParams(NPARAMS, hp, z, n0, kl0);
  x = e/eth;
  logx = log(x);
  RRRadialQkFromFit(NPARAMS, hp, 1, &x, &logx, &r, NULL, 0, &kl0);
  r *= x;
  if (os) return r;
  else {
    r *= 2.0*PI*FINE_STRUCTURE_CONST*AREA_AU20;
    return r;
  }
}

double RRCrossH(double z, int n0, int kl0, double e) {
  double hp[NPARAMS], eth;
  double x, logx, r;
  
  eth = 0.5*z*z/(n0*n0);
  RRRadialQkHydrogenicParams(NPARAMS, hp, z, n0, kl0);
  x = 1.0 + e/eth;
  logx = log(x);
  RRRadialQkFromFit(NPARAMS, hp, 1, &x, &logx, &r, NULL, 0, &kl0);
  r *= x;
  r *= 2.0*PI*FINE_STRUCTURE_CONST*AREA_AU20;
  r *= 2.0*(2.0*kl0 + 1.0);
  x = FINE_STRUCTURE_CONST*(e+eth);
  x = x*x;
  r *= x/(2.0*e);
  
  return r;
}
  
int RRRadialQkTable(double *qr, int k0, int k1, int m) {
  int index[3], k, nqk;
  double **p, *qk, tq[MAXNE];
  double r0, r1, tq0[MAXNE];
  ORBITAL *orb;
  int kappa0, jb0, klb02, klb0;
  int kappaf, jf, klf, kf;
  int jfmin, jfmax;
  int ite, ie, i;
  double eb, aw, e, pref;
  int mode, gauge;
  int nh, klh;
  double hparams[NPARAMS];

  orb = GetOrbital(k0);
  kappa0 = orb->kappa;
  GetJLFromKappa(kappa0, &jb0, &klb02);
  klb0 = klb02/2;
  eb = -orb->energy;

  for (ie = 0; ie < n_egrid; ie++) {
    xegrid[ie] = 1.0 + egrid[ie]/eb;
    log_xegrid[ie] = log(xegrid[ie]);
  }

  GetHydrogenicNL(&nh, &klh, NULL, NULL);
  if (m == -1) {
    r0 = GetResidualZ();
    RRRadialQkHydrogenicParams(NPARAMS, hparams, r0, orb->n, klb0);
    RRRadialQkFromFit(NPARAMS, hparams, n_egrid, 
		      xegrid, log_xegrid, tq0, NULL, 0, &klb0);
    if (klb0 > klh || orb->n > nh) {
      if (k0 != k1) return -1;
      for (ie = 0; ie < n_egrid; ie++) {
	qr[ie] = tq0[ie]/eb;
      }
      k = n_egrid;
      for (ite = 1; ite < n_tegrid; ite++) {
	for (ie = 0; ie < n_egrid; ie++) {
	  qr[k] = qr[ie];
	  k++;
	}
      }
      return 0;
    }
  }
  
  k = 2*abs(m);
  index[0] = k0;
  index[1] = k1;
  if (m >= 0) {
    index[2] = k;
  } else {
    index[2] = k-1;
  }

  nqk = n_tegrid*n_egrid;
  p = (double **) MultiSet(qk_array, index, NULL, InitPointerData);
  if (*p) {
    for (i = 0; i < nqk; i++) {
      qr[i] = (*p)[i];
    }
    return 0;
  }

  gauge = GetTransitionGauge();
  mode = GetTransitionMode();

  *p = (double *) malloc(sizeof(double)*nqk);
  
  qk = *p;
  /* the factor 2 comes from the conitinuum norm */
  pref = 2.0/((k+1.0)*(jb0+1.0));
  
  for (ite = 0; ite < n_tegrid; ite++) {
    for (ie = 0; ie < n_egrid; ie++) {
      tq[ie] = 0.0;
    }
    aw = FINE_STRUCTURE_CONST * tegrid[ite];
    jfmin = jb0 - k;
    jfmax = jb0 + k;
    for (jf = jfmin; jf <= jfmax; jf += 2) {
      for (klf = jf-1; klf <= jf+1; klf += 2) {
	if (jf <= 0 ||
	    klf < 0 ||
	    (m < 0 && IsOdd((klb02+klf+k)/2)) ||
	    (m > 0 && IsEven((klb02+klf+k)/2))) {
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
	  tq[ie] += r0*r1;
	}  
      }
    }

    for (ie = 0; ie < n_egrid; ie++) {
      qk[ie] = tq[ie]*pref;
    }
    qk += n_egrid;
  }
  
  for (i = 0; i < nqk; i++) {
    qr[i] = (*p)[i];
  }
  return 0;
}
  
int RRRadialQk(double *rqc, double te, int k0, int k1, int m) {
  int i, j, np, nd, k;
  double rq[MAXNTE];
  double x0, rqe[MAXNTE*MAXNE];

  i = RRRadialQkTable(rqe, k0, k1, m);
  if (i < 0) return -1;

  if (n_tegrid == 1) {
    for (i = 0; i < n_egrid; i++) {
      rqc[i] = rqe[i];
    }
  } else {
    nd = 1;
    np = 3;
    x0 = te;
    for (i = 0; i < n_egrid; i++) {
      j = i;
      for (k = 0; k < n_tegrid; k++) {
	rq[k] = rqe[j];
	j += n_egrid;
      }
      UVIP3P(np, n_tegrid, tegrid, rq, nd, &x0, &rqc[i]);
    }
  }
  return 0;
}

int BoundFreeOS(double *rqu, double *rqc, double *eb, 
		int rec, int f, int m) {
  LEVEL *lev1, *lev2;
  ANGULAR_ZFB *ang;
  ORBITAL *orb;
  int nz, ie, k;
  double a, b, d, amax, eb0, z;
  double rq[MAXNE], tq[MAXNE];
  int j1, j2;
  int i, j, c;
  int gauge, mode;
  int nkl, nq;
  int kb, kbp, jb, klb, jbp;

  lev1 = GetLevel(rec);
  lev2 = GetLevel(f);
  j1 = lev1->pj;
  j2 = lev2->pj;
  DecodePJ(j1, NULL, &j1);
  DecodePJ(j2, NULL, &j2);

  *eb = (lev2->energy - lev1->energy);
  if (*eb <= 0.0) return -1;

  nz = AngularZFreeBound(&ang, f, rec);
  if (nz <= 0) return -1;

  gauge = GetTransitionGauge();
  mode = GetTransitionMode();
  c = 2*abs(m) - 2;
  if (gauge == G_COULOMB && mode == M_FR && m < 0) {
    c -= 2;
  } 
  for (ie = 0; ie < n_egrid; ie++) {
    tq[ie] = 0.0;
  }
  amax = 0.0;
  for (i = 0; i < nz; i++) {
    kb = ang[i].kb;
    orb = GetOrbital(kb);
    jbp = orb->kappa;
    GetJLFromKappa(jbp, &jb, &klb);
    klb /= 2;
    for (j = 0; j <= i; j++) {
      kbp = ang[j].kb;
      jbp = GetOrbital(kbp)->kappa;
      jbp = GetJFromKappa(jbp);
      if (jbp != jb) continue;
      k = RRRadialQk(rq, *eb, kb, kbp, m);
      if (k < 0) continue;
      a = ang[i].coeff*ang[j].coeff;
      if (j != i) {
	a *= 2;
      } else {
	if (a > amax) {
	  nkl = klb;
	  nq = orb->n; 
	  amax = a;
	  eb0 = -(orb->energy);
	}
      }
      for (ie = 0; ie < n_egrid; ie++) {
	tq[ie] += a*rq[ie];
      }
    }
  }
  if (qk_mode == QK_FIT) {
    z = GetResidualZ();
    RRRadialQkHydrogenicParams(NPARAMS, rqc, z, nq, nkl);
    for (ie = 0; ie < n_egrid; ie++) {
      xegrid[ie] = 1.0 + egrid[ie]/eb0;
      log_xegrid[ie] = log(xegrid[ie]);
    }

    for (ie = n_egrid-2; ie > 2; ie--) {
      a = log(tq[ie+1]/tq[ie]);
      b = xegrid[ie+1]/xegrid[ie];
      d = (sqrt(xegrid[ie]) + rqc[2])/(sqrt(xegrid[ie+1]) + rqc[2]);
      b = log(b);
      d = log(d);
      z = (a + (4.5+nkl)*b)/(0.5*b+d);
      if (a < 0 && z > 0) {
	rqc[1] = z;
	break;
      }
    }
    RRRadialQkFromFit(NPARAMS, rqc, n_egrid, xegrid, log_xegrid, 
		      rq, NULL, 0, &nkl);
    ie++;
    a = eb0*tq[ie]/rq[ie];
    rqc[0] *= a;
    rqc[3] = eb0;
    for (ie++; ie < n_egrid; ie++) {
      tq[ie] = a*(rq[ie]/eb0);
    }
    for (ie = 0; ie < n_egrid; ie++) {
      a = (*eb) + egrid[ie];
      rqu[ie] = tq[ie]*a;
    }
  } else {
    for (ie = 0; ie < n_egrid; ie++) { 
      a = *eb + egrid[ie];
      tq[ie] *= a;
      if (c) {
	a *= FINE_STRUCTURE_CONST;
	tq[ie] *= pow(a, c);
      }
    }
    if (qk_mode == QK_INTERPOLATE) {
      for (ie = 0; ie < n_egrid; ie++) {
	tq[ie] = log(tq[ie]);
      }
      k = 3;
      UVIP3P(k, n_egrid, log_egrid, tq, n_usr, log_usr, rqu);
      for (ie = 0; ie < n_usr; ie++) {
	rqu[ie] = exp(rqu[ie]);
      }
    } else {
      for (ie = 0; ie < n_usr; ie++) {
	rqu[ie] = tq[ie];
      }
    }
  }      
 
  free(ang);

  return nkl;
}

int AutoionizeRate(double *rate, double *e, int rec, int f, int msub) {  
  LEVEL *lev1, *lev2;
  ANGULAR_ZxZMIX *ang;
  ANGULAR_ZFB *zfb;
  STATE *st;
  int k, nz, nzfb, ik, i, j1, j2, ij, kappaf, ip;
  int jf, k0, k1, kb, njf, nkappaf, klf, jmin, jmax;
  double *p, r, s, log_e, a;
  double *ai_pk, ai_pk0[MAXNE];
  int np, nt, m1, m2, m;
  int kappafp, jfp, klfp, dkl;

  *rate = 0.0;
  lev1 = GetLevel(rec);
  lev2 = GetLevel(f);
  
  *e = lev1->energy - lev2->energy;
  if (*e <= 0.0) return -1;
  log_e = log(*e);

  i = lev1->pb;

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
	  UVIP3P(np, n_egrid, log_egrid, ai_pk, nt, &log_e, &s);
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
	  UVIP3P(np, n_egrid, log_egrid, ai_pk0, nt, &log_e, &s);
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

  if (!msub) {
    r = 0.0;
    for (i = 0; i < nkappaf; i++) {
      r += p[i]*p[i];
    }
    /* the prefactor 4.0 includes the factor 2 from the continuum norm,
       otherwize, it should have been 2.0 */
    *rate = 4.0*r/(j1+1.0);
    free(p);    
    return 0;
  } else {
    k = 0;
    for (m1 = -j1; m1 <= 0; m1 += 2) {
      for (m2 = -j2; m2 <= j2; m2 += 2) {
	m = m1-m2;
	rate[k] = 0;
	rate[k+1] = 0;
	for (i = 0; i < nkappaf; i++) {
	  if (IsOdd(i)) {
	    ij = i-1;
	    klf = 1;
	  } else {
	    ij = i;
	    klf = -1;
	  }
	  jf = ij+jmin;
	  klf += jf;
	  kappaf = GetKappaFromJL(jf, klf);
	  s = W3j(j2, jf, j1, m2, m, -m1);
	  rate[k] += s*s*p[i]*p[i];
	  if (m != 1 && m != -1) continue;
	  for (ip = 0; ip < nkappaf; ip++) {
	    if (IsOdd(ip)) {
	      ij = ip-1;
	      klfp = 1;
	    } else {
	      ij = ip;
	      klfp = -1;
	    }
	    jfp = ij + jmin;
	    klfp += jfp;
	    if (ip == i) {
	      r = W3j(klf, 1, jf, 0, m, -m);
	      r = r*r*s*s*p[i]*p[i];
	      r *= (klf+1.0)*(jf+1.0);
	      rate[k+1] += r;
	    } else {
	      kappafp = GetKappaFromJL(jfp, klfp);
	      for (ik = 0; ik < n_egrid; ik++) {
		k0 = OrbitalIndex(0, kappaf, egrid[ik]);
		k1 = OrbitalIndex(0, kappafp, egrid[ik]);
		ai_pk0[ik] = GetPhaseShift(k0);
		ai_pk0[ik] -= GetPhaseShift(k1);
	      }	      
	      if (n_egrid > 1) {
		UVIP3P(np, n_egrid, log_egrid, ai_pk0, nt, &log_e, &a);
	      } else {
		a = ai_pk0[0];
	      }
	      r = W3j(klf, 1, jf, 0, m, -m);
	      r *= W3j(klfp, 1, jfp, 0, m, -m);
	      r *= W3j(j2, jfp, j1, m2, m, -m1);
	      r = r*s*p[i]*p[ip];
	      r *= sqrt((klf+1.0)*(klfp+1.0)*(jf+1.0)*(jfp+1.0));
	      r *= cos(a);
	      dkl = (jf + jfp)/2 + 1;
	      if (IsOdd(dkl)) r = -r;
	      rate[k+1] += r;
	    }
	  }
	}
	rate[k] *= 4.0;
	rate[k+1] *= 2.0*PI*PI/(*e);
	k += 2;
      }
    }
    free(p);    
    return k;
  }
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
  int i, kf;
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

  p = (double **) MultiSet(pk_array, index, NULL, InitPointerData);
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

int PrepRREGrids(double e, double emax0) { 
  double rmin, rmax;
  double emin, emax;
  int j;

  if (egrid_limits_type == 0) {
    rmin = egrid_min;
    rmax = egrid_max;
  } else {
    rmin = egrid_min/e;
    rmax = egrid_max/e;
  }
  emin = rmin*e;
  emax = rmax*e;
  if (emax < emax0) {
    emax = 50.0*e;
    if (emax > emax0) emax = emax0;
  }
  egrid_type = 1;
  if (usr_egrid_type < 0) usr_egrid_type = 1;

  if (qk_mode == QK_EXACT) {
    if (n_egrid <= 0) {
      if (n_usr <= 0) n_usr = 6;
      if (usr_egrid[0] < 0.0) {
	SetUsrPEGrid(n_usr, emin, emax, e);
	usr_egrid_type = 1;
      }  
      SetPEGridDetail(n_usr, usr_egrid);
    } else {
      if (egrid[0] < 0.0) {
	SetPEGrid(n_egrid, emin, emax, e);
	usr_egrid_type = 1;
      }
      SetUsrPEGridDetail(n_egrid, egrid);
    }
  } else {
    if (n_egrid == 0) {
      n_egrid = 6;
    }
    if (egrid[0] < 0.0) {
      SetPEGrid(n_egrid, emin, emax, e);
    }
    if (qk_mode == QK_INTERPOLATE) {
      if (n_usr <= 0) {
	SetUsrPEGridDetail(n_egrid, egrid);
	usr_egrid_type = 1;
      } else if (usr_egrid[0] < 0) {
	SetUsrPEGrid(n_usr, emin, emax, e);
	usr_egrid_type = 1;
      }
    } else if (qk_mode == QK_FIT) {
      SetUsrPEGridDetail(n_egrid, egrid);
      usr_egrid_type = 1;
    }
  }

  if (qk_mode == QK_INTERPOLATE) {
    for (j = 0; j < n_usr; j++) {
      log_usr[j] = usr_egrid[j];
      if (usr_egrid_type == 1) log_usr[j] += e;
      log_usr[j] = log(log_usr[j]);
    }
  }
  return 0;
}
  
int SaveRecRR(int nlow, int *low, int nup, int *up, 
	      char *fn, int m) {
  int i, j, k, ie, ip;
  FILE *f;
  double rqu[MAXNUSR], qc[NPARAMS+1];
  double eb;
  LEVEL *lev1, *lev2;
  RR_RECORD r;
  RR_HEADER rr_hdr;
  F_HEADER fhdr;
  double e, emin, emax, emax0;
  double awmin, awmax;
  int nq, nqk;
  ARRAY subte;
  int isub, n_tegrid0, n_egrid0, n_usr0;
  int te_set, e_set, usr_set;
  double c, e0, e1;

  if (m != -1 && GetTransitionGauge() != G_BABUSHKIN && qk_mode == QK_FIT) {
    printf("QK_FIT mode is only available to LENGTH form of E1 transitions\n");
    printf("Changing QK_FIT to QK_INTERPOLATE.\n");
    SetRecQkMode(QK_INTERPOLATE, -1.0);
  }

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
    return 0;
  }
  
  if (tegrid[0] < 0) {
    te_set = 0;
  } else {
    te_set = 1;
  }
  if (egrid[0] < 0) {
    e_set = 0;
  } else {
    e_set = 1;
  }
  if (usr_egrid[0] < 0) {
    usr_set = 0;
  } else {
    usr_set = 1;
  }
  n_tegrid0 = n_tegrid;
  n_egrid0 = n_egrid;
  n_usr0 = n_usr;

  if (egrid_limits_type == 0) {
    emax0 = 0.5*(emin + emax)*egrid_max;
  } else {
    emax0 = egrid_max;
  }
  ArrayInit(&subte, sizeof(double), 128);
  ArrayAppend(&subte, &emin, NULL);
  c = 1.0/TE_MIN_MAX;
  if (!e_set || !te_set) {
    e = c*emin;
    while (e < emax) {
      ArrayAppend(&subte, &e, NULL);
      e *= c;
    }
  }
  ArrayAppend(&subte, &emax, NULL);

  if (qk_mode == QK_FIT) {
    nqk = NPARAMS+1;
    r.params = (float *) malloc(sizeof(float)*nqk);
  } else {
    nqk = 0;
  }

  fhdr.type = DB_RR;
  strcpy(fhdr.symbol, GetAtomicSymbol());
  fhdr.atom = GetAtomicNumber();
  rr_hdr.nele = GetNumElectrons(low[0]);
  rr_hdr.qk_mode = qk_mode;
  rr_hdr.nparams = nqk;
  rr_hdr.multipole = m;
  f = OpenFile(fn, &fhdr);
  
  e0 = emin;
  for (isub = 1; isub < subte.dim; isub++) {
    e1 = *((double *) ArrayGet(&subte, isub));
    if (isub == subte.dim-1) e1 = e1*1.001;
    emin = e1;
    emax = e0;
    k = 0;
    for (i = 0; i < nlow; i++) {
      lev1 = GetLevel(low[i]);
      for (j = 0; j < nup; j++) {
	lev2 = GetLevel(up[j]);
	e = lev2->energy - lev1->energy;
	if (e < e0 || e >= e1) continue;
	if (e < emin) emin = e;
	if (e > emax) emax = e;
	k++;
      }
    }
    if (k == 0) {
      e0 = e1;
      continue;
    }
    emin = e0;
    emax = e1;  
    if (m == 1 || GetTransitionMode() == M_FR) {
      e = (emax - emin)/(0.5*(emin+emax));
      if (!te_set) {
	if (e < 0.1) {
	  SetRRTEGrid(1, 0.5*(emin+emax), emax);
	} else if (e < 0.5) {
	  SetRRTEGrid(2, emin, emax);
	} else {
	  if (k == 2) n_tegrid = 2;
	  else if (n_tegrid0 == 0) n_tegrid = 3;
	  SetRRTEGrid(n_tegrid, emin, emax);
	}
      }
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
    } else {
      SetRRTEGrid(1, 0.5*(emin+emax), emax);
    }
    
    n_egrid = n_egrid0;
    n_usr = n_usr0;
    if (!usr_set) usr_egrid[0] = -1.0;
    if (!e_set) egrid[0] = -1.0;
    e = 0.5*(emin + emax);
    PrepRREGrids(e, emax0);
    
    if (qk_mode == QK_FIT && n_egrid <= NPARAMS) {
      printf("n_egrid must > %d to use QK_FIT mode\n", NPARAMS);
      return -1;
    }
    rr_hdr.n_tegrid = n_tegrid;
    rr_hdr.tegrid = tegrid;
    rr_hdr.n_egrid = n_egrid;
    rr_hdr.egrid = egrid;
    rr_hdr.n_usr = n_usr;
    rr_hdr.usr_egrid = usr_egrid;
    rr_hdr.egrid_type = egrid_type;
    rr_hdr.usr_egrid_type = usr_egrid_type;
    
    r.strength = (float *) malloc(sizeof(float)*n_usr);
    
    InitFile(f, &fhdr, &rr_hdr);
    
    for (i = 0; i < nup; i++) {
      lev1 = GetLevel(up[i]);
      for (j = 0; j < nlow; j++) {
	lev2 = GetLevel(low[j]);
	e = lev1->energy - lev2->energy;
	if (e < e0 || e >= e1) continue;
	nq = BoundFreeOS(rqu, qc, &eb, low[j], up[i], m);
	if (nq < 0) continue;
	r.b = low[j];
	r.f = up[i];
	r.kl = nq;
	
	if (qk_mode == QK_FIT) {
	  for (ip = 0; ip < nqk; ip++) {
	    r.params[ip] = (float) qc[ip];
	  }
	}
	
	for (ie = 0; ie < n_usr; ie++) {
	  r.strength[ie] = (float) rqu[ie];
	}
	WriteRRRecord(f, &r);
      }
    }      

    DeinitFile(f, &fhdr);
    
    free(r.strength);
    ReinitRadial(1);
    FreeRecQk();
    FreeRecPk();
    
    e0 = e1;
  }

  if (qk_mode == QK_FIT) {
    free(r.params);
  }
      
  ReinitRecombination(1);

  ArrayFree(&subte, NULL);
  CloseFile(f, &fhdr);

  return 0;
}
      
int SaveAI(int nlow, int *low, int nup, int *up, char *fn, 
	   int channel, int msub) {
  int i, j, k, t;
  LEVEL *lev1, *lev2;
  AI_RECORD r;
  AIM_RECORD r1;
  AI_HEADER ai_hdr;
  AIM_HEADER ai_hdr1;
  F_HEADER fhdr;
  double emin, emax;
  double e, s, tai, a, s1[MAXAIM];
  float rt[MAXAIM];
  FILE *f;
  ARRAY subte;
  double c, e0, e1, b;
  int isub, n_egrid0;
  int e_set;

  if (nup <= 0 || nlow <= 0) return -1;

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
    return 0;
  }

  if (egrid[0] < 0) {
    e_set = 0;
  } else {
    e_set = 1;
  }

  n_egrid0 = n_egrid;

  ArrayInit(&subte, sizeof(double), 128);
  ArrayAppend(&subte, &emin, NULL);
  c = 1.0/TE_MIN_MAX;
  if (!e_set) {
    b = c*emin;
    while (b < emax) {
      ArrayAppend(&subte, &b, NULL);
      b *= c;
    }
  }
  ArrayAppend(&subte, &emax, NULL);
  
  if (!msub) {
    fhdr.type = DB_AI;
  } else {
    fhdr.type = DB_AIM;
  }
  strcpy(fhdr.symbol, GetAtomicSymbol());
  fhdr.atom = GetAtomicNumber();
  if (!msub) {
    ai_hdr.nele = GetNumElectrons(low[0]);
    ai_hdr.channel = channel;
  } else {
    ai_hdr1.nele = GetNumElectrons(low[0]);
    ai_hdr1.channel = channel;
  }
  f = OpenFile(fn, &fhdr);

  e0 = emin;
  for (isub = 1; isub < subte.dim; isub++) {
    e1 = *((double *) ArrayGet(&subte, isub));
    if (isub == subte.dim-1) e1 = e1*1.001;
    emin = e1;
    emax = e0;
    k = 0;
    for (i = 0; i < nlow; i++) {
      lev1 = GetLevel(low[i]);
      for (j = 0; j < nup; j++) {
	lev2 = GetLevel(up[j]);
	a = lev1->energy - lev2->energy;
	if (a < e0 || a >= e1) continue;
	if (a < emin) emin = a;
	if (a > emax) emax = a;
	k++;
      }
    }
    if (k == 0) {
      e0 = e1;
      continue;
    }
    if (!e_set) {
      a = (emax-emin)/(0.5*(emax+emin));
      if (a < 0.1) {
	a = 0.5*(emin+emax);
	SetPEGrid(1, a, a, 0.0);
      } else if (a < 0.4) {
	SetPEGrid(2, emin, emax, 0.0);
      } else if (a < 1.0) {
	if (k == 2) n_egrid = 2;
	else if (n_egrid0 == 0)	n_egrid = 3;
	SetPEGrid(n_egrid, emin, emax, 0.0);
      } else {
	if (k == 2) n_egrid = 2;
	else if (n_egrid0 == 0) n_egrid = 4;
	SetPEGrid(n_egrid, emin, emax, 0.0);
      }
    }      
    
    if (!msub) {
      ai_hdr.n_egrid = n_egrid;
      ai_hdr.egrid = egrid;
      InitFile(f, &fhdr, &ai_hdr);
    } else {
      ai_hdr1.n_egrid = n_egrid;
      ai_hdr1.egrid = egrid;
      InitFile(f, &fhdr, &ai_hdr1);
    }
    for (i = 0; i < nlow; i++) {
      lev1 = GetLevel(low[i]);
      for (j = 0; j < nup; j++) {
	lev2 = GetLevel(up[j]);
	e = lev1->energy - lev2->energy;
	if (e < e0 || e >= e1) continue;
	if (!msub) {
	  k = AutoionizeRate(&s, &e, low[i], up[j], msub);
	  if (k < 0) continue;
	  if (s < ai_cut) continue;
	  r.b = low[i];
	  r.f = up[j];
	  r.rate = s;
	  WriteAIRecord(f, &r);
	} else {
	  k = AutoionizeRate(s1, &e, low[i], up[j], msub);
	  if (k < 0) continue;
	  r1.rate = rt;
	  s = 0;
	  for (t = 0; t < k; t++) {
	    r1.rate[t] = s1[t];
	    s += s1[t];
	  }
	  if (s < ai_cut) continue;
	  r1.b = low[i];
	  r1.f = up[j];
	  r1.nsub = k;
	  WriteAIMRecord(f, &r1);
	}
      }
    }

    DeinitFile(f, &fhdr);

    ReinitRadial(1);
    FreeRecQk();
    FreeRecPk();

    e0 = e1;
  }

  ReinitRecombination(1);
  
  ArrayFree(&subte, NULL);
  CloseFile(f, &fhdr);

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

static void FreeRecPkData(void *p) {
  double *dp;
  dp = *((double **) p);
  free(dp);
  *((double **) p) = NULL;
}

int FreeRecPk(void) {
  if (pk_array->array == NULL) return 0;
  MultiFreeData(pk_array->array, pk_array->ndim, FreeRecPkData);
  return 0;
}

int FreeRecQk(void) {
  if (qk_array->array == NULL) return 0;
  MultiFreeData(qk_array->array, qk_array->ndim, FreeRecPkData);
  return 0;
}

int FreeRecAngZ(void) {
  int i, j, *k;
  for (i = 0; i < n_complex; i++) {
    for (j = 0; j < (rec_complex[i].rg)->dim; j++) {
      k = (int *)ArrayGet(rec_complex[i].rg, j);  
      if (k) FreeAngZ(*k, -1);
    }
  }
  return 0;
}

int InitRecombination(void) {
  int blocks[5];
  int ndim;
  int i;
  
  ndim = 5;
  for (i = 0; i < ndim; i++) blocks[i] = MULTI_BLOCK5;
  pk_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(pk_array, sizeof(double *), ndim, blocks);
  
  ndim = 3;
  for (i = 0; i < ndim; i++) blocks[i] = MULTI_BLOCK3;
  qk_array = (MULTI *) malloc(sizeof(MULTI));
  blocks[0] = 10;
  blocks[1] = 10;
  blocks[2] = 4;
  MultiInit(qk_array, sizeof(double *), ndim, blocks);  
  
  hyd_qk_array = (ARRAY *) malloc(sizeof(ARRAY));
  ArrayInit(hyd_qk_array, sizeof(double *), 32);
  
  n_complex = 0;
  for (i = 0; i < MAX_COMPLEX; i++) {
    rec_complex[i].rg = (ARRAY *) malloc(sizeof(ARRAY));
    ArrayInit(rec_complex[i].rg, sizeof(int), 64);
  }
  n_egrid = 0;
  egrid[0] = -1.0;
  SetPEGridLimits(-1, -1, 0);
  n_usr = 0;
  usr_egrid[0] = -1.0;
  
  n_tegrid = 0;
  tegrid[0] = -1.0;

  SetRecQkMode(QK_DEFAULT, 0.1);
  SetRecPWOptions(RECLMAX, RECLMAX);
  return 0;
}

int ReinitRecombination(int m) {
  int i;

  if (m < 0) return 0;

  FreeRecQk();
  FreeRecPk();

  n_egrid = 0;
  egrid[0] = -1.0;
  SetPEGridLimits(-1, -1, 0);
  n_usr = 0;
  usr_egrid[0] = -1.0;
  
  n_tegrid = 0;
  tegrid[0] = -1.0;

  SetRecQkMode(QK_DEFAULT, 0.1);
  SetRecPWOptions(RECLMAX, RECLMAX);

  if (m > 0) return 0;

  for (i = 0; i < n_complex; i++) {
    ArrayFree(rec_complex[i].rg, NULL);
  }
  n_complex = 0;
  
  return 0;
}
