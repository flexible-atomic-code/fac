#include "recombination.h"
#include "time.h"

static char *rcsid="$Id: recombination.c,v 1.38 2002/01/17 19:52:23 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

#define NPARAMS 3

static int qk_mode;
static double qk_fit_tolerance;

static int egrid_type = -1;
static int usr_egrid_type = -1;

static int n_egrid = 0;
static double egrid[MAXNE];
static double log_egrid[MAXNE];
static double xegrid[MAXNTE][MAXNE];
static double log_xegrid[MAXNTE][MAXNE];
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
} pw_scratch = {8, 8, 500, MAXNKL, 10, 0, 0, {0, MAXKL}};

double ai_cut = EPS8;

static REC_COMPLEX rec_complex[MAX_COMPLEX];
int n_complex = 0;

void uvip3p_(int *np, int *ndp, double *x, double *y, 
	     int *n, double *xi, double *yi);
void pixz1_(double *z, double *am, int *ne, int *nl, double *phe,
	    double *pc, double *pcp, double *pnc, int *nde, int *ndl,
	    int *iopt, int *ipcp);

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
  if (max <= 0) egrid_max = 10.0;
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
      exit(1);
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
      m = lev->basis[0];
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
	  m = lev->basis[0];
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
  SaveLevels(fn, nlevels, -1);

  return 0;
}

void RRRadialQkHydrogenicParams(int np, double *p, int n, int kl) {
#define NNE 10
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

  qk = (double **) ArraySet(hyd_qk_array, n, NULL);
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
      logx[ne-1] = log(50.0);
      d = (logx[ne-1] - logx[0])/(ne-1.0);
      for (i = 1; i < ne; i++) {
	logx[i] = logx[i-1] + d;
	x[i] = exp(logx[i]);
      }
      i = 200;
      pixz1_(&z, &am, &i, &n, e, pc, NULL, pnc, 
	     &ne, &n, &iopt, &ipcp);
      iopt = 0;
    }
    
    eth = 0.5/(n*n);
    for (i = 0; i < ne; i++) {
      e[i] = (x[i]-1.0)*eth;
    }
    pixz1_(&z, &am, &ne, &n, e, pc, NULL, pnc, 
	   &ne, &n, &iopt, &ipcp);
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
  z = GetResidualZ();
  p[0] = t[0]/(z*z);
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
      k += ndy;
      dy[k] = a*p[1]*(1.0/e - 1.0/d);
    }
  }
}

int RRRadialQkTable(double *qr, int k0, int k1, int m) {
  int index[3], k, nqk;
  double **p, *qk, tq[MAXNE], sig[MAXNE];
  double r0, r1;
  ORBITAL *orb;
  int kappa0, jb0, klb02, klb0;
  int kappaf, jf, klf, kf;
  int jfmin, jfmax;
  int ite, ie, i;
  double eb, aw, e, pref, chisq, tol;
  int mode, gauge;
  int ipvt[NPARAMS];
  int ierr;
  int lwa=5*NPARAMS+MAXNE;
  double wa[5*NPARAMS+MAXNE];
  double fvec[MAXNE], fjac[MAXNE*NPARAMS];
  int nh, klh;
  double hparams[NPARAMS];

  orb = GetOrbital(k0);
  kappa0 = orb->kappa;
  GetJLFromKappa(kappa0, &jb0, &klb02);
  klb0 = klb02/2;
  GetHydrogenicNL(&nh, &klh);
  if (m == -1) {
    if (qk_mode == QK_FIT || orb->n >= nh || klb0 >= klh) {
      RRRadialQkHydrogenicParams(NPARAMS, hparams, orb->n, klb0);
    }
    if (orb->n >= nh || klb0 >= klh) {
      if (k0 != k1) return -1;
      if (qk_mode != QK_FIT) {
	for (ite = 0; ite < n_tegrid; ite++) {
	  RRRadialQkFromFit(NPARAMS, hparams, n_egrid, 
			    xegrid[ite], log_xegrid[ite], qr, NULL, 0, &klb0);
	  for (ie = 0; ie < n_egrid; ie++) {
	    qr[ie] /= tegrid[ite];
	  }
	  qr += n_tegrid;
	}
      } else {
	k = 0;
	for (ite = 0; ite < n_tegrid; ite++) {
	  for (i = 0; i < NPARAMS; i++) {
	    qr[k] = hparams[i];
	    k++;
	  }
	}
	qr[k] = klb0;
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

  if (qk_mode == QK_FIT) {
    nqk = n_tegrid*NPARAMS+1;
  } else {
    nqk = n_tegrid*n_egrid;
  }
  p = (double **) MultiSet(qk_array, index, NULL);
  if (*p) {
    for (i = 0; i < nqk; i++) {
      qr[i] = (*p)[i];
    }
    return 0;
  }

  gauge = GetTransitionGauge();
  mode = GetTransitionMode();

  if (qk_mode == QK_FIT) {
    *p = (double *) malloc(sizeof(double)*nqk);
    (*p)[n_tegrid*NPARAMS] = klb0;
    tol = qk_fit_tolerance * 1E-3;
    tol = Max(tol, 1E-4);
  } else {
    *p = (double *) malloc(sizeof(double)*nqk);
  }
  qk = *p;
  /* the factor 2 comes from the conitinuum norm */
  pref = 2.0/((k+1.0)*(jb0+1.0));
  
  for (ite = 0; ite < n_tegrid; ite++) {
    eb = tegrid[ite];
    for (ie = 0; ie < n_egrid; ie++) {
      tq[ie] = 0.0;
    }
    aw = FINE_STRUCTURE_CONST * eb;
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

    if (qk_mode == QK_FIT) {
      for (i = 0; i < NPARAMS; i++) {
	qk[i] = hparams[i];
      }
      for (ie = 0; ie < n_egrid; ie++) {
	tq[ie] *= eb*pref;
	sig[ie] = tq[ie];
      }
      ierr = NLSQFit(NPARAMS, qk, tol, ipvt, fvec, fjac, MAXNE, wa, lwa, 
		     n_egrid, xegrid[ite], log_xegrid[ite], 
		     tq, sig, RRRadialQkFromFit, &klb0);
      chisq = 0.0;
      for (ie = 0; ie < n_egrid; ie++) {
	chisq += fvec[ie]*fvec[ie];
      }
      chisq = sqrt(chisq/n_egrid);
      if (ierr > 3 || chisq > qk_fit_tolerance) {
	for (i = 0; i < n_egrid; i++) {
	  if (xegrid[ite][i] > 4.0) {
	    sig[i] = tq[i]*1E3;
	  } else {
	    sig[i] = (xegrid[ite][i]/xegrid[ite][0])*tq[i];
	  }
	}
	for (i = 0; i < NPARAMS; i++) {
	  qk[i] = hparams[i];
	}
	ierr = NLSQFit(2, qk, tol, ipvt, fvec, fjac, MAXNE, 
		       wa, lwa, n_egrid, xegrid[ite], log_xegrid[ite], 
		       tq, sig, RRRadialQkFromFit, &klb0);
      }      
      qk += NPARAMS;
    } else {
      for (ie = 0; ie < n_egrid; ie++) {
	qk[ie] = tq[ie]*pref;
      }
      qk += n_egrid;
    }
  }
  
  for (i = 0; i < nqk; i++) {
    qr[i] = (*p)[i];
  }
  return 0;
}
  
int RRRadialQk(double *rqc, double te, int k0, int k1, int m) {
  int i, j, np, nd, k;
  double rq[MAXNTE];
  double x0, rqe[MAXNTE*MAXNE+1];

  i = RRRadialQkTable(rqe, k0, k1, m);
  if (i < 0) return -1;

  if (qk_mode == QK_FIT) {
    if (n_tegrid == 1) {
      k = rqe[NPARAMS];
      for (i = 0; i < NPARAMS; i++) {
	rqc[i] = rqe[i];
      }
      rqc[0] *= (te/tegrid[0]);
    } else {
      nd = 1;
      np = 3;
      x0 = te;
      for (i = 0; i < NPARAMS; i++) {
	j = i;
	for (k = 0; k < n_tegrid; k++) {
	  rq[k] = rqe[j];
	  j += NPARAMS;
	}
	uvip3p_(&np, &n_tegrid, tegrid, rq, &nd, &x0, &rqc[i]);
      }
      k = rqe[NPARAMS*n_tegrid];
    }
    rqc[NPARAMS] = k;
    return k;
  } else {
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
	uvip3p_(&np, &n_tegrid, tegrid, rq, &nd, &x0, &rqc[i]);
      }
    }
    return 0;
  }
}

int BoundFreeOS(double *rqu, int *nqkc, double **rqc, double *eb, 
		int rec, int f, int m) {
  LEVEL *lev1, *lev2;
  ANGULAR_ZFB *ang;
  int nz, ie, k;
  double a;
  double rq[MAXNE], tq[MAXNUSR];
  int j1, j2;
  int i, j, c;
  int gauge, mode;
  int nq, nqk;
  double *p;
  int kb, kbp, jb, jbp;

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
  if (qk_mode == QK_FIT) {
    for (i = 0; i < n_usr; i++) {
      xusr[i] = usr_egrid[i]/(*eb);
      if (usr_egrid_type == 1) xusr[i] += 1.0;
      log_xusr[i] = log(xusr[i]);
    }
    nqk = NPARAMS+1;
    for (i = 0; i < n_usr; i++) {
      rqu[i] = 0.0;
    }
    nq = 0;
    p = *rqc;
    for (i = 0; i < nz; i++) {
      kb = ang[i].kb;
      jb = GetOrbital(kb)->kappa;
      jb = GetJFromKappa(jb);
      for (j = 0; j <= i; j++) {
        kbp = ang[j].kb;
        jbp = GetOrbital(kbp)->kappa;
        jbp = GetJFromKappa(jbp);
	if (jbp != jb) continue;
	if (nq == *nqkc) {
	  *nqkc *= 2;
	  *rqc = (double *) realloc(*rqc, sizeof(double)*(*nqkc));
	  p = *rqc + nq*nqk;
	}
	k = RRRadialQk(p, *eb, kb, kbp, m);
	if (k < 0) continue;
	a = ang[i].coeff*ang[j].coeff;    
	if (j != i) {
	  a *= 2;
	} 
	p[0] *= a;
	RRRadialQkFromFit(NPARAMS, p, n_usr, xusr, log_xusr, 
			  tq, NULL, 0, (void *) &k);
	for (ie = 0; ie < n_usr; ie++) {
	  rqu[ie] += tq[ie];
	}
	nq++;
	p += nqk;
      }
    }
    for (ie = 0; ie < n_usr; ie++) {
      rqu[ie] = xusr[ie]*rqu[ie];
    }
  } else {
    gauge = GetTransitionGauge();
    mode = GetTransitionMode();
    c = 2*abs(m) - 2;
    if (gauge == G_COULOMB && mode == M_FR && m < 0) {
      c -= 2;
    } 
    for (ie = 0; ie < n_egrid; ie++) {
      tq[ie] = 0.0;
    }
    for (i = 0; i < nz; i++) {
      kb = ang[i].kb;
      jb = GetOrbital(kb)->kappa;
      jb = GetJFromKappa(jb);
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
	} 
	for (ie = 0; ie < n_egrid; ie++) {
	  tq[ie] += a*rq[ie];
	}
      }
    }
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
      uvip3p_(&k, &n_egrid, log_egrid, tq, &n_usr, log_usr, rqu);
      for (ie = 0; ie < n_usr; ie++) {
	rqu[ie] = exp(rqu[ie]);
      }
    } else if (qk_mode == QK_EXACT) {
      for (ie = 0; ie < n_usr; ie++) {
	rqu[ie] = tq[ie];
      }
    }      
  }
 
  free(ang);

  return nq;
}

int AutoionizeRate(double *rate, double *e, int rec, int f) {  
  LEVEL *lev1, *lev2;
  ANGULAR_ZxZMIX *ang;
  ANGULAR_ZFB *zfb;
  STATE *st;
  int k, nz, nzfb, ik, i, j1, j2, ij, kappaf, ip;
  int jf, k0, k1, kb, njf, nkappaf, klf, jmin, jmax;
  double *p, r, s, log_e;
  double *ai_pk, ai_pk0[MAXNE];
  int np, nt;

  *rate = 0.0;
  lev1 = GetLevel(rec);
  lev2 = GetLevel(f);
  
  *e = lev1->energy - lev2->energy;
  if (*e <= 0.0) return -1;
  log_e = log(*e);

  i = lev1->basis[0];

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

int PrepRREGrids(double e) { 
  double rmin, rmax;
  double emin, emax;
  int i, j;

  if (egrid_limits_type == 0) {
    rmin = egrid_min;
    rmax = egrid_max;
  } else {
    rmin = egrid_min/e;
    rmax = egrid_max/e;
  }
  emin = rmin*e;
  emax = rmax*e;
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

  for (j = 0; j < n_egrid; j++) {
    log_egrid[j] = log(egrid[j] + e);
    for (i = 0; i < n_tegrid; i++) {
      xegrid[i][j] = 1.0 + egrid[j]/tegrid[i];
      log_xegrid[i][j] = log(xegrid[i][j]);
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
  double rqu[MAXNUSR], *qc;
  double eb;
  LEVEL *lev1, *lev2;
  RR_RECORD r;
  RR_HEADER rr_hdr;
  F_HEADER fhdr;
  double e, emin, emax;
  double awmin, awmax;
  int nqkc, nq, nqk, nshells;

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
  if (emin < TE_MIN_MAX*emax) {
    emin = TE_MIN_MAX*emax;
  }

  if (k == 0) {
    return 0;
  }
  
  if (m == 1 || GetTransitionMode() == M_FR) {
    e = (emax - emin)/(0.5*(emin+emax));
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
  
  if (m != -1 && GetTransitionGauge() != G_BABUSHKIN && qk_mode == QK_FIT) {
    printf("QK_FIT mode is only available to LENGTH form of E1 transitions\n");
    printf("Changing QK_FIT to QK_INTERPOLATE.\n");
    SetRecQkMode(QK_INTERPOLATE, -1.0);
  }

  e = 0.5*(emin + emax);
  PrepRREGrids(e);

  if (qk_mode == QK_FIT && n_egrid <= NPARAMS) {
    printf("n_egrid must > %d to use QK_FIT mode\n", NPARAMS);
    return -1;
  }

  if (qk_mode == QK_FIT) {
    nqkc = 10;
    nqk = NPARAMS+1;
    qc = (double *) malloc(sizeof(double)*nqkc*nqk);
  } else {
    nqk = 0;
  }

  fhdr.type = DB_RR;
  strcpy(fhdr.symbol, GetAtomicSymbol());
  fhdr.atom = GetAtomicNumber();
  rr_hdr.nele = GetNumElectrons(low[0]);
  rr_hdr.qk_mode = qk_mode;
  rr_hdr.nparams = nqk;
  rr_hdr.n_tegrid = n_tegrid;
  rr_hdr.n_egrid = n_egrid;
  rr_hdr.egrid_type = egrid_type;
  rr_hdr.n_usr = n_usr;
  rr_hdr.usr_egrid_type = usr_egrid_type;
  rr_hdr.multipole = m;
  rr_hdr.tegrid = tegrid;
  rr_hdr.egrid = egrid;
  rr_hdr.usr_egrid = usr_egrid;

  f = InitFile(fn, &fhdr, &rr_hdr);

  nshells = 1;
  if (qk_mode == QK_FIT) {
    r.params = (float *) malloc(sizeof(float)*nqk);
  }
  r.strength = (float *) malloc(sizeof(float)*n_usr);
  
  for (i = 0; i < nup; i++) {
    for (j = 0; j < nlow; j++) {
      nq = BoundFreeOS(rqu, &nqkc, &qc, &eb, low[j], up[i], m);
      if (nq < 0) continue;
      r.b = low[j];
      r.f = up[i];
      r.nshells = nq;
      if (r.nshells > nshells && qk_mode == QK_FIT) {
	r.params = (float *) realloc(r.params, sizeof(float)*nqk*nq);
	nshells = r.nshells;
      }
      
      if (qk_mode == QK_FIT) {
	ip = 0;
	for (k = 0; k < nq; k++) {
	  for (ie = 0; ie < nqk; ie++) {
	    r.params[ip] = (float) qc[ip];
	    ip++;
	  }
	}
      }

      for (ie = 0; ie < n_usr; ie++) {
	r.strength[ie] = (float) rqu[ie];
      }
      WriteRRRecord(f, &r);
    }
  }
  
  if (qk_mode == QK_FIT) free(r.params);
  free(r.strength);
  CloseFile(f, &fhdr);
  free(qc);

  return 0;
}

int SaveDR(int nf, int *f, int na, int *a, int nb, int *b, int ng, int *g, 
	   char *fna, char *fnt, int channel) {
  int i, j, k, j1, m, do_transition;
  LEVEL *lev1, *lev2;
  TR_RECORD rt;
  AI_RECORD ra;
  TR_HEADER tr_hdr;
  AI_HEADER ai_hdr;
  F_HEADER fhdra, fhdrt;
#ifdef PERFORM_STATISTICS
  clock_t start, stop;
  clock_t t_ai, t_rd, tt_ai, tt_rd;
  STRUCT_TIMING st_start, st_stop;
  ANGULAR_TIMING angt;
  RECOUPLE_TIMING recouplet;
  RAD_TIMING radt;
#endif
  double emin, emax;
  double *ea, *sa, tai, e0; 
  double *et, *sr, *rd, elow, trd, trd1, tr_cut;
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
  if (emin < TE_MIN_MAX*emax) {
    emin = TE_MIN_MAX*emax;
  }

  if (k == 0) {
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

  
  fhdra.type = DB_AI;
  strcpy(fhdra.symbol, GetAtomicSymbol());
  fhdra.atom = GetAtomicNumber();  
  ai_hdr.nele = GetNumElectrons(a[0]);
  ai_hdr.channel = channel;
  ai_hdr.n_egrid = n_egrid;
  ai_hdr.egrid = egrid;  
  fa = InitFile(fna, &fhdra, &ai_hdr);

  fhdrt.type = DB_TR;
  strcpy(fhdrt.symbol, GetAtomicSymbol());
  fhdrt.atom = GetAtomicNumber();
  tr_hdr.nele = ai_hdr.nele;
  tr_hdr.multipole = -1;
  tr_hdr.gauge = GetTransitionGauge();
  tr_hdr.mode = GetTransitionMode();
  ft = InitFile(fnt, &fhdrt, &tr_hdr);

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
      k = AutoionizeRate(sa+j, ea+j, a[i], f[j]);
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
      ra.b = a[i];
      ra.f = f[j];
      ra.rate = sa[j];
      WriteAIRecord(fa, &ra);
    }

#ifdef PERFORM_STATISTICS
    stop = clock();
    t_ai = stop-start;
    start = stop;
    tt_ai += t_ai;
#endif
    
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
      rt.upper = a[i];
      rt.lower = b[j];
      rt.strength = sr[j];
      WriteTRRecord(ft, &rt);
    }
#ifdef PERFORM_STATISTICS
    stop = clock();
    t_rd = stop-start;
    start = stop;
    tt_rd += t_rd;
#endif
  }

  free(sa);
  free(ea);
  free(sr);
  free(rd);
  free(et);
  CloseFile(fa, &fhdra);
  CloseFile(ft, &fhdrt);

#ifdef PERFORM_STATISTICS
  GetStructTiming(&st_stop);
  
  fprintf(perform_log, "Time in AI: %6.1E, RD: %6.1E\n", 
	  ((double) tt_ai)/CLOCKS_PER_SEC, ((double) tt_rd)/CLOCKS_PER_SEC);
  
  fprintf(perform_log, "AngZMix: %6.1E, AngZFB: %6.1E, AngZxZFB: %6.1E, SetH: %6.1E DiagH: %6.1E\n",
	  ((double) (st_stop.angz_mix))/CLOCKS_PER_SEC,
	  ((double) (st_stop.angz_fb))/CLOCKS_PER_SEC,
	  ((double) (st_stop.angzxz_fb))/CLOCKS_PER_SEC,
	  ((double) (st_stop.set_ham))/CLOCKS_PER_SEC,
	  ((double) (st_stop.diag_ham))/CLOCKS_PER_SEC);
  fprintf(perform_log, "AngZS: %6.1E, AngZFBS: %6.1E, AngZxZFBS: %6.1E, AddZ: %6.1E, AddZxZ: %6.1E\n",
	  ((double) (st_stop.angz_states))/CLOCKS_PER_SEC,
	  ((double) (st_stop.angzfb_states))/CLOCKS_PER_SEC,
	  ((double) (st_stop.angzxzfb_states))/CLOCKS_PER_SEC,
	  ((double) (st_stop.add_angz))/CLOCKS_PER_SEC,
	  ((double) (st_stop.add_angzxz))/CLOCKS_PER_SEC);

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
  fprintf(perform_log, "Dirac: %6.1E, 1E: %6.1E, 2E: %6.1E\n", 
	  ((double)radt.dirac)/CLOCKS_PER_SEC, 
	  ((double)radt.radial_1e)/CLOCKS_PER_SEC,
	  ((double)radt.radial_2e)/CLOCKS_PER_SEC);
  fprintf(perform_log, "\n");
#endif /* PERFORM_STATISTICS */
  
  return 0;
}

      
int SaveAI(int nlow, int *low, int nup, int *up, char *fn, int channel) {
  int i, j, k;
  LEVEL *lev1, *lev2;
  AI_RECORD r;
  AI_HEADER ai_hdr;
  F_HEADER fhdr;
  double emin, emax;
  double *e, *s, tai, a;
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
  if (emin < TE_MIN_MAX*emax) {
    emin = TE_MIN_MAX*emax;
  }

  if (k == 0) {
    return 0;
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
  if (nup <= 0 || nlow <= 0) return -1;
 
  fhdr.type = DB_AI;
  strcpy(fhdr.symbol, GetAtomicSymbol());
  fhdr.atom = GetAtomicNumber();
  ai_hdr.nele = GetNumElectrons(low[0]);
  ai_hdr.channel = channel;
  ai_hdr.n_egrid = n_egrid;
  ai_hdr.egrid = egrid;

  f = InitFile(fn, &fhdr, &ai_hdr);

  s = malloc(sizeof(double)*nup);
  e = malloc(sizeof(double)*nup);
  for (i = 0; i < nlow; i++) {
    tai = 0.0;
    for (j = 0; j < nup; j++) {
      k = AutoionizeRate(s+j, e+j, low[i], up[j]);
      if (k < 0) continue;
      tai += s[j];
    }
    if (tai < 1E-30) continue;
    for (j = 0; j < nup; j++) {
      if (s[j] < ai_cut*tai) continue;
      r.b = low[i];
      r.f = up[j];
      r.rate = s[j];
      WriteAIRecord(f, &r);
    }
  }

  free(s);
  free(e);
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

static void _FreeRecPk(void *p) {
  double *dp;
  dp = *((double **) p);
  free(dp);
}

int FreeRecPk(void) {
  if (pk_array->array == NULL) return 0;
  MultiFreeData(pk_array->array, pk_array->ndim, _FreeRecPk);
  return 0;
}

int FreeRecQk(void) {
  if (qk_array->array == NULL) return 0;
  MultiFreeData(qk_array->array, qk_array->ndim, _FreeRecPk);
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
  
  hyd_qk_array = (ARRAY *) malloc(sizeof(ARRAY));
  ArrayInit(hyd_qk_array, sizeof(double *), 10);
  
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
  SetRecPWOptions(12, 12);
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
  SetRecPWOptions(12, 12);

  if (m > 0) return 0;

  for (i = 0; i < n_complex; i++) {
    ArrayFree(rec_complex[i].rg, NULL);
  }
  n_complex = 0;
  
  return 0;
}






