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

#include "recombination.h"
#include "time.h"
#include "cf77.h"

static char *rcsid="$Id$";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

static int qk_mode;
static double qk_fit_tolerance;

static int egrid_type = -1;
static int usr_egrid_type = -1;

static int n_egrid = 0;
static int uta_egrid = 0;
static double egrid[MAXRRNE];
static double log_egrid[MAXRRNE];
static double egrid_min;
static double egrid_max;
static int egrid_limits_type = 0;

static int n_usr = 0;
static double usr_egrid[MAXRRNUSR];
static double log_usr[MAXRRNUSR];
static double xusr[MAXRRNUSR];
static double log_xusr[MAXRRNUSR];

static int n_tegrid = 0;
static int uta_tegrid = 0;
static double tegrid[MAXNTE];
static double log_te[MAXNTE];

static int n_cxegrid = 0;
static int t_cxegrid = 0;
static double *cxegrid0 = NULL;
static double *cxegrid = NULL;
static int cxldist = 5;
static double cxldistmj = 1.0;
static int cxrotcouple = 1.0;

static MULTI *pk_array;
static MULTI *qk_array;

#define MAXAIM 4096
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

static double _ai_cut = AICUT;
static double _rr_cut = RRCUT;
static int _asym_ci = 1;
static int _progress_report = 0;

static REC_COMPLEX rec_complex[MAX_COMPLEX];
int n_complex = 0;

static int maxaicache = MAXAICACHE;
static AICACHE aicache = {0, 0};

void SetMaxAICache(int n) {
  if (n > 0) {
    maxaicache = n;
  } else {
    maxaicache = MAXAICACHE;
  }
  FreeAICache(0);
}

void AllocAICache(void) {
  aicache.low = malloc(sizeof(int)*maxaicache);
  aicache.up = malloc(sizeof(int)*maxaicache);
  aicache.e = malloc(sizeof(double)*maxaicache);
  aicache.nz = malloc(sizeof(int)*maxaicache);
  aicache.nzf = malloc(sizeof(int)*maxaicache);
  aicache.az = malloc(sizeof(ANGULAR_ZxZMIX *)*maxaicache);
  aicache.azf = malloc(sizeof(ANGULAR_ZFB *)*maxaicache);
  int i;
  for (i = 0; i < maxaicache; i++) {
    aicache.low[i] = -1;
    aicache.up[i] = -1;
    aicache.e[i] = 0;
    aicache.nz[i] = 0;
    aicache.nzf[i] = 0;
    aicache.az[i] = NULL;
    aicache.azf[i] = NULL;
  }
}

void FreeAICache(int m) {
  int i;
  for (i = 0; i < aicache.nc; i++) {
    if (aicache.nz[i] > 0) {
      free(aicache.az[i]);
      aicache.az[i] = NULL;
      aicache.nz[i] = 0;
    }
    if (aicache.nzf[i] > 0) {
      free(aicache.azf[i]);
      aicache.azf[i] = NULL;
      aicache.nzf[i] = 0;
    }
  }
  if (m == 0) {
    free(aicache.low);
    free(aicache.up);
    free(aicache.e);
    free(aicache.nz);
    free(aicache.nzf);
    free(aicache.az);
    free(aicache.azf);
  }
}

static void FreeRecPkData(void *p) {
  double *dp;
  dp = *((double **) p);
  free(dp);
  *((double **) p) = NULL;
}

int SetAICut(double c) {
  _ai_cut = c;
  return 0;
}

int SetRRTEGrid(int n, double emin, double emax) {
  if (n == -1) {
    n_tegrid = 1;
    uta_tegrid = 1;
    tegrid[0] = -1;
    return 1;
  }
  uta_tegrid = 0;
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
  int i;
  
  if (n == -1) {
    n_egrid = 1;
    uta_egrid = 1;
    egrid[0] = 0;
    return 1;
  }
  uta_egrid = 0;
  n_egrid = SetEGrid(egrid, log_egrid, n, emin, emax, eth);
  if (uta_tegrid) {
    for (i = 0; i < n_egrid; i++) {
      egrid[i] /= eth;
      log_egrid[i] = log(egrid[i]);
    }
  }
  return n_egrid;
}

int SetUsrPEGridDetail(int n, double *xg) {
  if (n > MAXRRNUSR) {
    printf("Max # of grid points reached \n");
    return -1;
  }
  n_usr = SetEGridDetail(usr_egrid, log_usr, n, xg);
  return n_usr;
}
  						      
int SetUsrPEGrid(int n, double emin, double emax, double eth) {
  int i;
  
  if (n > MAXRRNUSR) {
    printf("Max # of grid points reached \n");
    return -1;
  }
  n_usr = SetEGrid(usr_egrid, log_usr, n, emin, emax, eth);
  if (uta_tegrid) {
    for (i = 0; i < n_usr; i++) {
      usr_egrid[i] /= eth;
      log_usr[i] = log(usr_egrid[i]);
    }
  }
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
	int nmax = GetOrbNMax(ns.kappa, 0);
	if (nmax > 0 && n > nmax) continue;
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
  if (ncfgs == 0) return 0;
  nlevels = GetNumLevels();
  nsym = MAX_SYMMETRIES;
  rec_complex[n_complex].n = n;
  rec_complex[n_complex].s0 = nlevels;
  for (i = 0; i < nsym; i++) {
    HAMILTON *h = GetHamilton(i);
    m = ConstructHamilton(i, k, k, kg, 0, NULL, 111);
    if (m < 0) continue;
    j = DiagnolizeHamilton(h);
    if (j < 0) return -1;
    AddToLevels(h, 0, NULL);
  }
  rec_complex[n_complex].s1 = GetNumLevels()-1;
  n_complex++;
  SortLevels(nlevels, -1, 0);
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
	int nmax = GetOrbNMax(pw_scratch.kappa0[j], 0);
	if (nmax > 0 && n > nmax) continue;
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

  if (nstates == 0) return 0;
  nsym = MAX_SYMMETRIES;
  nlevels = GetNumLevels();
  rec_complex[n_complex].n = n;
  rec_complex[n_complex].s0 = nlevels;
  for (i = 0; i < nsym; i++) {
    HAMILTON *h = GetHamilton(i);
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
	  m = ConstructHamiltonFrozen(i, j, NULL, n, 0, NULL, 0, 0, -1, -1);
	  if (m < 0) continue;
	  if (DiagnolizeHamilton(h) < 0) return -2;
	  AddToLevels(h, 0, NULL);
	}
	i0 = rec_complex[t].s1+1;
      }
    } else {
      m = ConstructHamiltonFrozen(i, k, kg, n, 0, NULL, 0, 0, -1, -1);
      if (m < 0) continue;
      if (DiagnolizeHamilton(h) < 0) return -2;
      AddToLevels(h, 0, NULL);
    }
  }
  rec_complex[n_complex].s1 = GetNumLevels()-1;
  n_complex++;
  SortLevels(nlevels, -1, 0);
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
#pragma omp critical
  {
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

double GauntBF(int n, int k, double z, double te, double e) {
  double ryd = HARTREE_EV/2;
  double z2 = z*z;
  int n2 = n*n;
  double en = z2*ryd/n2;
  double ef = exp(en/te);
  double r = PICrossH(z, e, n, k, 0);
  r *= 1.26e-3*pow(e/ryd,3)*ef*ryd/te;
  return r;
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

int RRRadialMultipoleTable(double *qr, int k0, int k1, int m) {
  int index[3], k, nqk;
  double **p, *qk;
  int kappaf, jf, klf, kf;
  int ite, ie, i;
  double aw, e, pref, te;
  int gauge;

  klf = k1;
  if (IsOdd(klf)) {
    klf++;
    jf = klf-1;
  } else {
    jf = klf+1;
  }
  kappaf = GetKappaFromJL(jf, klf);

  k = 2*abs(m);
  index[0] = k0;
  index[1] = k1;
  if (m >= 0) {
    index[2] = k;
  } else {
    index[2] = k-1;
  }

  nqk = n_tegrid*n_egrid;
  LOCK *lock = NULL;
  int locked = 0;
  p = (double **) MultiSet(qk_array, index, NULL, &lock,
			   InitPointerData, FreeRecPkData);
  if (lock && !(*p)) {
    SetLock(lock);
    locked = 1;
  }
  if (*p) {
    for (i = 0; i < nqk; i++) {
      qr[i] = (*p)[i];
    }
    if (locked) ReleaseLock(lock);
    return 0;
  }

  gauge = GetTransitionGauge();

  double *pd = (double *) malloc(sizeof(double)*nqk);
  
  qk = pd;
  /* the factor 2 comes from the conitinuum norm */
  pref = sqrt(2.0);
  te = 0.0;
  if (uta_tegrid) {
    te = fabs(GetOrbital(k0)->energy);
    aw = FINE_STRUCTURE_CONST*te;
  }
  for (ite = 0; ite < n_tegrid; ite++) {
    if (!uta_tegrid) aw = FINE_STRUCTURE_CONST * tegrid[ite];        
    for (ie = 0; ie < n_egrid; ie++) {
      if (uta_tegrid) {
	e = egrid[ie]*te;
      } else {
	e = egrid[ie];
      }
      kf = OrbitalIndex(0, kappaf, e);
      qk[ie] = MultipoleRadialFR(aw, m, k0, kf, gauge);
      qk[ie] *= pref;
    }
    qk += n_egrid;
  }
  
  for (i = 0; i < nqk; i++) {
    qr[i] = pd[i];
  }
  *p = pd;
  if (locked) ReleaseLock(lock);
#pragma omp flush
  return 0;
}

void TopUpRRRadialQk(double *qr, int k0, int kp, int mx) {
  int ite, ie, k, mx1;
  double r, rp, b, te, e;

  if (mx <= 1) return;
  if (mx > 100) mx = 100;

  if (uta_tegrid) te = fabs(GetOrbital(k0)->energy);
  k = 0;
  for (ite = 0; ite < n_tegrid; ite++) {
    if (!uta_tegrid) te = tegrid[ite];
    for (ie = 0; ie < n_egrid; ie++) {
      mx1 = mx;
      /*
      r = AsymmetryPI(k0, k0, egrid[ie], tegrid[ite], 0, 0, &mx1,
		      0, &b, NULL, NULL);
      if (r > 0) {
	if (kp != k0) {
	  mx1 = mx;
	  rp = AsymmetryPI(kp, kp, egrid[ie], tegrid[ite], 0, 0, &mx1,
			   0, &b, NULL, NULL);
	  if (rp > 0) {
	    r = sqrt(r*rp);
	  }
	}
	qr[k] /= r;
      }
      */
      if (uta_tegrid) {
	e = egrid[ie]*te;
      } else {
	e = egrid[ie];
      }
      r = AsymmetryPI(k0, kp, e, te, 0, 0, &mx1,
		      0, &b, NULL, NULL);
      if (r > 0) qr[k] /= r;
      k++;
    }
  }
}

int RRRadialQkTable(double *qr, int k0, int k1, int m0) {
  int index[3], k, nqk;
  double **p, *qk, tq[MAXRRNE];
  double r0, r1, tq0[MAXRRNE];
  ORBITAL *orb;
  int kappa0, jb0, klb02, klb0;
  int kappaf, jf, klf, kf;
  int jfmin, jfmax;
  int ite, ie, i;
  double eb, aw, e, pref;
  int gauge;
  int nh, klh;
  double hparams[NPARAMS];
  double xegrid[MAXRRNE], log_xegrid[MAXRRNE];
 
  orb = GetOrbital(k0);
  kappa0 = orb->kappa;
  GetJLFromKappa(kappa0, &jb0, &klb02);
  klb0 = klb02/2;
  eb = -orb->energy;

  for (ie = 0; ie < n_egrid; ie++) {
    if (uta_tegrid) {
      xegrid[ie] = 1 + egrid[ie];
    } else {
      xegrid[ie] = 1.0 + egrid[ie]/eb;
    }
    log_xegrid[ie] = log(xegrid[ie]);
  }

  GetHydrogenicNL(&nh, &klh, NULL, NULL);
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
    TopUpRRRadialQk(qr, k0, k0, m0);
    return 0;
  }

  k = 2;
  index[0] = k0;
  index[1] = k1;
  
  LOCK *lock = NULL;
  int locked = 0;
  nqk = n_tegrid*n_egrid;
  p = (double **) MultiSet(qk_array, index, NULL, &lock,
			   InitPointerData, FreeRecPkData);
  if (lock && !(*p)) {
    SetLock(lock);
    locked = 1;
  }
  if (*p) {
    for (i = 0; i < nqk; i++) {
      qr[i] = (*p)[i];
    }
    if (locked) {      
      ReleaseLock(lock);
    }
    return 0;
  }

  gauge = GetTransitionGauge();

  double *pd = (double *) malloc(sizeof(double)*nqk);
  
  qk = pd;
  /* the factor 2 comes from the conitinuum norm */
  pref = 2.0/((k+1.0)*(jb0+1.0));
  if (uta_tegrid) {
    aw = FINE_STRUCTURE_CONST * eb;
  }
  for (ite = 0; ite < n_tegrid; ite++) {
    for (ie = 0; ie < n_egrid; ie++) {
      tq[ie] = 0.0;
    }
    if (!uta_tegrid) aw = FINE_STRUCTURE_CONST * tegrid[ite];
    jfmin = jb0 - k;
    jfmax = jb0 + k;
    for (jf = jfmin; jf <= jfmax; jf += 2) {
      for (klf = jf-1; klf <= jf+1; klf += 2) {
	if (jf <= 0 ||
	    klf < 0 ||
	    IsOdd((klb02+klf+k)/2)) {
	  continue;
	}
	kappaf = GetKappaFromJL(jf, klf);
	for (ie = 0; ie < n_egrid; ie++) {
	  if (uta_tegrid) {
	    e = egrid[ie]*eb;
	  } else {
	    e = egrid[ie];
	  }
	  kf = OrbitalIndex(0, kappaf, e);
	  r0 = MultipoleRadialFR(aw, -1, k0, kf, gauge);
	  if (k1 == k0) {
	    r1 = r0;
	  } else {
	    r1 = MultipoleRadialFR(aw, -1, k1, kf, gauge);
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
  TopUpRRRadialQk(pd, k0, k1, m0);
  for (i = 0; i < nqk; i++) {
    qr[i] = pd[i];
  }
  *p = pd;
  if (locked) ReleaseLock(lock);
#pragma omp flush
  return 0;
}
  
int RRRadialMultipole(double *rqc, double te, int k0, int k1, int m) {
  int i, j, np, nd, k;
  double rq[MAXNTE];
  double x0, rqe[MAXNTE*MAXRRNE];

  i = RRRadialMultipoleTable(rqe, k0, k1, m);
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
  
int RRRadialQk(double *rqc, double te, int k0, int k1, int m) {
  int i, j, np, nd, k;
  double rq[MAXNTE];
  double x0, rqe[MAXNTE*MAXRRNE];

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

int BoundFreeMultipole(FILE *fp, int rec, int f, int m) {
  LEVEL *lev1, *lev2;
  ANGULAR_ZFB *ang;
  ORBITAL *orb;
  int i, nz, ie, k, kb, j1, j2, n, p1, p2;
  int jmin, jmax, jt, jfmin, jfmax, jf, klf, kf, jb, klb;
  double rq[MAXRRNE], rqu[MAXRRNE], eb, a;

  lev1 = GetLevel(rec);
  lev2 = GetLevel(f);
  j1 = lev1->pj;
  j2 = lev2->pj;
  DecodePJ(j1, &p1, &j1);
  DecodePJ(j2, &p2, &j2);

  eb = (lev2->energy - lev1->energy);
  if (eb <= 0.0) return -1;
  
  eb = (lev2->energy - lev1->energy);
  if (eb <= 0.0) return -1;

  nz = AngularZFreeBound(&ang, f, rec);
  if (nz <= 0) return -1;

  k = abs(m)*2;
  jmax = j1 + k;
  jmin = abs(j1 - k);
  for (jt = jmin; jt <= jmax; jt += 2) {
    jfmin = abs(jt - j2);
    jfmax = jt + j2;
    for (jf = jfmin; jf <= jfmax; jf += 2) {
      for (klf = jf-1; klf <= jf+1; klf += 2) {    
	kf = klf;
	if (jf < klf) kf--;	
	for (ie = 0; ie < n_egrid; ie++) {
	  rqu[ie] = 0.0;
	}
	if (jf <= 0 ||
	    klf < 0 ||
	    (m < 0 && IsOdd(p1+p2+klf/2+k/2)) ||
	    (m > 0 && IsEven(p1+p2+klf/2+k/2))) {
	  continue;
	}
	for (i = 0; i < nz; i++) {
	  kb = ang[i].kb;
	  orb = GetOrbital(kb);
	  GetJLFromKappa(orb->kappa, &jb, &klb);
	  if ((m < 0 && IsOdd((klb+klf+k)/2)) ||
	      (m > 0 && IsEven((klb+klf+k)/2))) {
	    continue;
	  }
	  n = RRRadialMultipole(rq, eb, kb, kf, m);
	  if (n < 0) continue;
	  a = ang[i].coeff*sqrt(jt+1.0)*W6j(j2, jf, jt, k, j1, jb);
	  if (IsOdd((jt+j1-k)/2)) a = -a;
	  if (fabs(a) > EPS16) {
	    for (ie = 0; ie < n_egrid; ie++) {
	      rqu[ie] += a * rq[ie];
	    }
	  }
	}
	for (ie = 0; ie < n_egrid; ie++) {
	  fprintf(fp, "%5d %2d %5d %2d  %2d  %3d %3d %3d %12.5E %12.5E %12.5E\n",
		  rec, j1, f, j2, m, jt, klf, jf, eb, egrid[ie], rqu[ie]);
	}
      }
    }
  }

  fprintf(fp, "\n");
  free(ang);

  return 0;
}

int BoundFreeOSUTA1(double *rqu, double *rqc, double *eb, 
		    int rec, int f, int m0) {
  INTERACT_DATUM *idatum;
  LEVEL *lev1, *lev2;
  int j1, ns, q1, ie;
  ORBITAL *orb;
  double a, b, d, eb0, z;
  double rq[MAXRRNE], tq[MAXRRNE];
  double xegrid[MAXRRNE], log_xegrid[MAXRRNE];
  int nkl, nq, k;
  int klb, jb, kb;
  
  lev1 = GetLevel(rec);
  lev2 = GetLevel(f);
  
  *eb = (lev2->energy - lev1->energy);
  if (*eb <= 0.0) return -1;

  idatum = NULL;
  ns = GetInteract(&idatum, NULL, NULL, lev2->iham, lev1->iham,
		   lev2->pb, lev1->pb, 0, 0, 1);  
  if (ns <= 0) return -1;
  if (idatum->s[1].index < 0 || idatum->s[3].index >= 0) {
    free(idatum->bra);
    free(idatum);
    return -1;
  }

  j1 = idatum->s[1].j;
  q1 = idatum->s[1].nq_ket;
  kb = OrbitalIndex(idatum->s[1].n, idatum->s[1].kappa, 0.0);
  orb = GetOrbital(kb);
  eb0 = -(orb->energy);
  GetJLFromKappa(orb->kappa, &jb, &klb);
  klb /= 2;
  
  for (ie = 0; ie < n_egrid; ie++) {
    tq[ie] = 0.0;
  }
  
  k = RRRadialQk(rq, *eb, kb, kb, m0);
  nq = orb->n;
  nkl = klb;
  for (ie = 0; ie < n_egrid; ie++) {
    tq[ie] += (jb+1.0)*rq[ie];
  }

  if (qk_mode == QK_FIT) {
    z = GetResidualZ();
    RRRadialQkHydrogenicParams(NPARAMS, rqc, z, nq, nkl);
    for (ie = 0; ie < n_egrid; ie++) {
      if (uta_tegrid) {
	xegrid[ie] = 1 + egrid[ie];
      } else {
	xegrid[ie] = 1.0 + egrid[ie]/eb0;
      }
      log_xegrid[ie] = log(xegrid[ie]);
    }

    for (ie = n_egrid-2; ie > 2; ie--) {
      a = log(tq[ie+1]/tq[ie]);
      b = xegrid[ie+1]/xegrid[ie];
      d = (sqrt(xegrid[ie]) + rqc[2])/(sqrt(xegrid[ie+1]) + rqc[2]);
      b = log(b);
      d = log(d) + 0.5*b;
      if (d > 0.05) {
	z = (a + (4.5+nkl)*b)/d;
	if (a < 0 && z > 0) {
	  rqc[1] = z;
	  break;
	}
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
      if (uta_tegrid) {
	a = (*eb)*(1+egrid[ie]);
      } else {
	a = (*eb) + egrid[ie];
      }
      rqu[ie] = tq[ie]*a;
    }
  } else {
    for (ie = 0; ie < n_egrid; ie++) {
      if (uta_tegrid) {
	a = (*eb)*(1+egrid[ie]);
      } else {
	a = (*eb) + egrid[ie];
      }
      tq[ie] *= a;
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

  d = (lev1->ilev+1.0)*(q1/(j1+1.0));
  for (ie = 0; ie < n_usr; ie++) {
    rqu[ie] *= d;
  }
  rqc[0] *= d;

  free(idatum->bra);
  free(idatum);
  return nkl;
}

void BoundFreeOSFit(double *rqu, double *rqc,
		    double *tq, int nq, int nkl,
		    double eb, double eb0) {
  int ie;
  double a, b, d, z;
  double rq0[MAXRRNE];
  double xegrid[MAXRRNE], log_xegrid[MAXRRNE];
  
  z = GetResidualZ();
  RRRadialQkHydrogenicParams(NPARAMS, rqc, z, nq, nkl);
  for (ie = 0; ie < n_egrid; ie++) {
    if (uta_tegrid) {
      xegrid[ie] = 1 + egrid[ie];
    } else {
      xegrid[ie] = 1.0 + egrid[ie]/eb0;
    }
    log_xegrid[ie] = log(xegrid[ie]);
  }

  for (ie = n_egrid-2; ie > 2; ie--) {
    a = log(tq[ie+1]/tq[ie]);
    b = xegrid[ie+1]/xegrid[ie];
    d = (sqrt(xegrid[ie]) + rqc[2])/(sqrt(xegrid[ie+1]) + rqc[2]);
    b = log(b);
    d = log(d) + 0.5*b;
    if (d > 0.05) {
      z = (a + (4.5+nkl)*b)/d;
      if (a < 0 && z > 0) {
	rqc[1] = z;
	break;
      }
    }
  }
  RRRadialQkFromFit(NPARAMS, rqc, n_egrid, xegrid, log_xegrid, 
		    rq0, NULL, 0, &nkl);
  ie++;
  a = eb0*tq[ie]/rq0[ie];
  rqc[0] *= a;
  rqc[3] = eb0;
  for (ie++; ie < n_egrid; ie++) {
    tq[ie] = a*(rq0[ie]/eb0);
  }
  for (ie = 0; ie < n_egrid; ie++) {
    if (uta_tegrid) {
      a = eb*(1+egrid[ie]);
    } else {
      a = eb + egrid[ie];
    }
    rqu[ie] = tq[ie]*a;
  }
}

int BoundFreeOSUTAShell(int m0, int kgb, int kcb,
			int kb, double *rq, int *qb) {
  int k, i;
  CONFIG *c;
  ORBITAL *orb;
  double eb;

  orb = GetOrbital(kb);  
  c = GetConfigFromGroup(kgb, kcb);
  for (i = 0; i < c->n_shells; i++) {
    if (c->shells[i].n == orb->n && c->shells[i].kappa == orb->kappa) {
      break;
    }
  }
  *qb = 0.0;
  if (i == c->n_shells) return 0;
  *qb = c->shells[i].nq;
  k = FactorNR(c, orb->n, orb->kappa);
  if (k > 0) {
    *qb = k;
  }
  eb = -orb->energy;
  k = RRRadialQk(rq, eb, kb, kb, m0);

  return 0;
}

int BoundFreeOSUTA0(int m0, int kgf, int kgb, int kcf, int kcb,
		    double eb, double *rq, int *kb, int *qb) {
  INTERACT_DATUM *idatum;
  int ns, k, ie;
  CONFIG *c;
  
  idatum = NULL;  
  ns = GetInteract(&idatum, NULL, NULL, kgf, kgb, kcf, kcb, 0, 0, 1);  
  if (ns <= 0) return -1;
  if (idatum->s[1].index < 0 || idatum->s[3].index >= 0) {
    free(idatum->bra);
    free(idatum);
    return -1;
  }
  
  *qb = idatum->s[1].nq_ket;
  c = GetConfigFromGroup(kgb, kcb);
  k = FactorNR(c, idatum->s[1].n, idatum->s[1].kappa);
  if (k > 0) {
    *qb = k;
  }
  *kb = OrbitalIndex(idatum->s[1].n, idatum->s[1].kappa, 0.0);      
  k = RRRadialQk(rq, eb, *kb, *kb, m0);
  free(idatum->bra);
  free(idatum);
  return 0;
}

int BoundFreeOSUTA(double *rqu, double *rqc, double *eb, 
		   int rec, int f, int m0) {
  SYMMETRY *sym;
  STATE *s;
  CONFIG *cfg;
  LEVEL *lev1, *lev2;
  int t, j, p;
  int k1, q1, ns, ie;
  ORBITAL *orb;
  double a, eb0, wb, wm;
  double rq[MAXRRNE], tq[MAXRRNE];
  double xegrid[MAXRRNE], log_xegrid[MAXRRNE];
  int nkl, nq, k, r;
  int klb, jb, kb, qb;
  
  lev1 = GetLevel(rec);
  lev2 = GetLevel(f);
  if (lev1->n_basis > 0 && lev2->n_basis > 0) {
    return BoundFreeOS(rqu, rqc, eb, rec, f, m0);
  }
  *eb = (lev2->energy - lev1->energy);
  if (*eb <= 0.0) return -1;
    
  if (lev1->n_basis == 0 && lev2->n_basis == 0) {
    r = BoundFreeOSUTA0(m0, lev2->iham, lev1->iham, lev2->pb, lev1->pb,
			*eb, rq, &kb, &qb);    
    if (r < 0) return -1;
    wb = qb*(lev1->ilev+1.0);
    for (ie = 0; ie < n_egrid; ie++) {
      tq[ie] = rq[ie]*wb;
    }
  } else if (lev1->n_basis > 0) {
    sym = GetSymmetry(lev1->pj);
    DecodePJ(lev1->pj, &p, &j);
    wm = 0;
    for (ie = 0; ie < n_egrid; ie++) tq[ie] = 0.0;
    for (t = 0; t < lev1->n_basis; t++) {
      s = (STATE *) ArrayGet(&(sym->states), lev1->basis[t]);
      r = BoundFreeOSUTA0(m0, lev2->iham, s->kgroup, lev2->pb, s->kcfg,
			  *eb, rq, &k1, &q1);
      if (r < 0) continue;
      wb = (j+1.0)*q1*lev1->mixing[t]*lev1->mixing[t];
      if (wb > wm) {
	wm = wb;
	kb = k1;
	qb = q1;
      }
      for (ie = 0; ie < n_egrid; ie++) {
	tq[ie] += wb*rq[ie];
      }
    }
  } else if (lev2->n_basis > 0) {
    sym = GetSymmetry(lev2->pj);
    DecodePJ(lev2->pj, &p, &j);
    wm = 0;
    for (ie = 0; ie < n_egrid; ie++) tq[ie] = 0.0;
    for (t = 0; t < lev2->n_basis; t++) {
      s = (STATE *) ArrayGet(&(sym->states), lev2->basis[t]);
      cfg = GetConfig(s);
      r = BoundFreeOSUTA0(m0, s->kgroup, lev1->iham, s->kcfg, lev1->pb,
			  *eb, rq, &k1, &q1);
      if (r < 0) continue;
      orb = GetOrbital(k1);
      wb = (lev1->ilev+1.0)*q1*lev2->mixing[t]*lev2->mixing[t];
      wb *= (j+1.0)/fabs(cfg->sweight);
      if (wb > wm) {
	wm = wb;
	kb = k1;
	qb = q1;
      }
      for (ie = 0; ie < n_egrid; ie++) {
	tq[ie] += wb*rq[ie];
      }
    }
  }
  a = 0.0;
  for (ie = 0; ie < n_egrid; ie++) a += tq[ie];
  if (a < 1e-31) return -1;
  
  orb = GetOrbital(kb);
  eb0 = -orb->energy;
  GetJLFromKappa(orb->kappa, &jb, &klb);
  klb /= 2;
  nq = orb->n;
  nkl = klb;
  if (qk_mode == QK_FIT) {
    BoundFreeOSFit(rqu, rqc, tq, nq, nkl, *eb, eb0);
  } else {
    for (ie = 0; ie < n_egrid; ie++) {
      if (uta_tegrid) {
	a = (*eb)*(1+egrid[ie]);
      } else {
	a = (*eb) + egrid[ie];
      }
      tq[ie] *= a;
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

  return nq*1000+nkl;
}
    
int RecOccupation(int **nk, double **nq, double **dn,
		  double *eb, int rec, int f) {
  LEVEL *lev1, *lev2;
  ANGULAR_ZFB *ang;
  ORBITAL *orb;
  int nz, nmax, kmax, k, i, j, kb, kbp;
  double a, *qk;
  
  lev1 = GetLevel(rec);
  lev2 = GetLevel(f);
  *eb = (lev2->energy - lev1->energy);
  nz = AngularZFreeBound(&ang, f, rec);
  if (nz <= 0) return 0;
  nmax = 0;
  for (i = 0; i < nz; i++) {
    kb = ang[i].kb;
    orb = GetOrbital(kb);
    if (orb->n > nmax) nmax = orb->n;
  }
  kmax = nmax*nmax;
  qk = malloc(sizeof(double)*kmax);
  for (k = 0; k < kmax; k++) qk[k] = 0.0;
  for (i = 0; i < nz; i++) {
    kb = ang[i].kb;
    orb = GetOrbital(kb);
    k = ShellToInt(orb->n, orb->kappa);
    a = ang[i].coeff;
    qk[k] += a;
  }
  nmax = 0;
  for (k = 0; k < kmax; k++) {
    if (qk[k]) nmax++;
  }
  *nk = malloc(sizeof(int)*nmax);
  *nq = malloc(sizeof(double)*nmax);
  *dn = malloc(sizeof(double)*nmax);
  i = 0;
  for (k = 0; k < kmax; k++) {
    if (qk[k]) {
      IntToShell(k, &kb, &kbp);
      j = OrbitalIndex(kb, kbp, 0);
      orb = GetOrbital(j);
      (*dn)[i] = orb->dn;
      (*nk)[i] = kb;
      GetJLFromKappa(kbp, &j, &kb);
      (*nk)[i] = 100*(*nk)[i] + kb/2;
      if (j < kb) {
	(*nk)[i] = -(*nk)[i];
      }
      (*nq)[i] = qk[k];
      i++;
    }
  }
  free(qk);
  free(ang);
  return nmax;
}
    
int CXCross(CXTGT *cxt, double *cx, double *eb, int rec, int f,
	    int nmc, double **mc) {
  LEVEL *lev1, *lev2;
  ANGULAR_ZFB *ang;
  ORBITAL *orb;
  int nz, k, i, j, m, kb, kbp, j1, j2, mn, mk, ie;
  double a, c, z, am, e;
  double orx, ov12, ode, ovdr, olam, obeta;
  
  lev1 = GetLevel(rec);
  lev2 = GetLevel(f);
  j1 = lev1->pj;
  j2 = lev2->pj;
  DecodePJ(j1, NULL, &j1);
  DecodePJ(j2, NULL, &j2);
  z = GetResidualZ();
  if (lev1->ibase < 0 || lev1->ibase == f) m = 0;
  else m = 1;
  *eb = (lev2->energy - lev1->energy);
  if (*eb <= 0) return 0;
  nz = AngularZFreeBound(&ang, f, rec);
  if (nz <= 0) return 0;
  am = 0;
  mn = 0;
  mk = 0;
  for (ie = 0; ie < n_cxegrid; ie++) {
    cx[ie] = 0.0;
  }

  for (i = 0; i < nz; i++) {
    kb = ang[i].kb;
    orb = GetOrbital(kb);
    k = GetLFromKappa(orb->kappa)/2;
    for (j = 0; j <= i; j++) {
      kbp = ang[j].kb;
      if (kbp != kb) continue;
      a = ang[i].coeff*ang[j].coeff;
      if (a > am) {
	am = a;
	mn = orb->n;
	mk = k;
      }
      if (j != i) {
	a *= 2;
      }
      orx = -1;
      ov12 = -1;
      ovdr = -1;
      olam = -1;
      ode = -1;
      obeta = -1;
      a *= 1.0/(2*(2*k+1.0)*(j2+1.0));
      for (ie = 0; ie < n_cxegrid; ie++) {
	e = cxegrid[ie];
	c = LandauZenerCX(cxt, orb->n, k, m, z, e, *eb,
			  &orx, &ov12, &ovdr, &olam, &ode, &obeta);
	if (orb->n <= nmc && c > 0) {
	  c = log(c) + mc[orb->n-1][ie];
	  if (c <= -225) c = 0.0;
	  else c = exp(c);
	}
	cx[ie] += a*c;
      }
    }
  }
  free(ang);
  for (ie = 0; ie < n_cxegrid; ie++) {
    if (cx[ie] > 0) break;
  }
  if (ie >= n_cxegrid) return 0;
  return mn*100 + mk;
}

int BoundFreeOS(double *rqu, double *rqc, double *eb, 
		int rec, int f, int m0) {
  LEVEL *lev1, *lev2;
  ANGULAR_ZFB *ang;
  ORBITAL *orb;
  int nz, ie, k;
  double a, amax, eb0;
  double rq0[MAXRRNE], rq[MAXRRNE], tq[MAXRRNE];
  double xegrid[MAXRRNE], log_xegrid[MAXRRNE];
  int i, j;
  int nkl, nq;
  int kb, kbp, jb, klb, jbp;

  lev1 = GetLevel(rec);
  lev2 = GetLevel(f);

  *eb = (lev2->energy - lev1->energy);
  if (*eb <= 0.0) return -1;
  nz = AngularZFreeBound(&ang, f, rec);
  if (nz <= 0) return -1;
  
  for (ie = 0; ie < n_egrid; ie++) {
    tq[ie] = 0.0;
  }
  amax = 0.0;
  int kp0;
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
      k = RRRadialQk(rq, *eb, kb, kbp, m0);
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
	  kp0 = orb->kappa;
	  for (ie = 0; ie < n_egrid; ie++) rq0[ie] = rq[ie];
	}
      }
      for (ie = 0; ie < n_egrid; ie++) {
	tq[ie] += a*rq[ie];
      }
    }
  }
  if (qk_mode == QK_FIT) {
    BoundFreeOSFit(rqu, rqc, tq, nq, nkl, *eb, eb0);
  } else {
    for (ie = 0; ie < n_egrid; ie++) {
      if (uta_tegrid) {
	a = (*eb)*(1+egrid[ie]);
      } else {
	a = (*eb) + egrid[ie];
      }
      tq[ie] *= a;
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

  return nq*1000+nkl;
}

int AutoionizeRateUTA0(double *rate, double *e,
		       int kgf, int kgb, int kcf, int kcb) {
  INTERACT_DATUM *idatum;
  int j0, j1, jb, ns, q0, q1, qb;
  int k0, k1, kb, kmin, kmax, jmin, jmax;
  int jf, ik, klf, kappaf, k, np, nt, j, jm;
  double a, b, r, s, log_e, *ai_pk;
  int ni, tp, nleft, id, nd=2500, done[2500];
  CONFIG *cfg;
  
  *rate = 0.0;
  log_e = log(*e);
  
  idatum = NULL;
  ns = GetInteract(&idatum, NULL, NULL, kgf, kgb, kcf, kcb, 0, 0, 1);
  if (ns <= 0) return -1;
  if (idatum->s[3].index < 0) {
    free(idatum->bra);
    free(idatum);
    return -1;
  }

  kb = OrbitalIndex(idatum->s[1].n, idatum->s[1].kappa, 0.0);
  k0 = OrbitalIndex(idatum->s[2].n, idatum->s[2].kappa, 0.0);
  k1 = OrbitalIndex(idatum->s[3].n, idatum->s[3].kappa, 0.0);
  j0 = idatum->s[2].j;
  j1 = idatum->s[3].j;
  jb = idatum->s[1].j;
  q0 = idatum->s[2].nq_ket;
  q1 = idatum->s[3].nq_ket;
  qb = idatum->s[1].nq_ket;
  cfg = GetConfigFromGroup(kgb, kcb);
  tp = FactorNR(cfg, idatum->s[2].n, idatum->s[2].kappa);
  if (tp > 0) {
    q0 = tp;
  }
  tp = FactorNR(cfg, idatum->s[1].n, idatum->s[1].kappa);
  if (tp > 0) {
    qb = tp;
  }
  for (id = 0; id < nd; id++) done[id] = 0;
  if (idatum->s[1].index != idatum->s[3].index) {
    kmin = abs(j0-j1);
    kmax = j0 + j1;
    jmin = 1;
    jmax = j1+j0+jb;
    np = 3;
    nt = 1;
    r = 0.0;
    ni = 0;
    while (1) {
      nleft = 0;
      id = 0;
      for (jf = jmin; jf <= jmax; jf += 2) {
	for (ik = -1; ik <= 1; ik += 2) {
	  klf = jf + ik;
	  kappaf = GetKappaFromJL(jf, klf);
	  for (k = kmin; k <= kmax; k += 2) {
	    if (!Triangle(j0, j1, k) || !Triangle(jb, jf, k)) continue;
	    if (id >= nd) {
	      if (ni == 0) {
		tp = AIRadialPk(&ai_pk, *e, k0, k1, kb, kappaf, k, 0);
	      } else {
		continue;
	      }
	    } else {
	      if (done[id]) {
		id++;
		continue;
	      }
	      tp = AIRadialPk(&ai_pk, *e, k0, k1, kb, kappaf, k, 1);
	      if (tp == -9999) {
		id++;
		nleft++;
		continue;
	      }
	      id++;
	      done[id] = 1;
	    }
	    if (n_egrid > 1) {
	      UVIP3P(np, n_egrid, log_egrid, ai_pk, nt, &log_e, &s);
	    } else {
	      s = ai_pk[0];
	    }
	    s = s*s/(k + 1.0);
	    r += s;
	  }
	}
      }
      if (nleft == 0) break;
      ni++;
    }
    r *= 4.0*(q1/(j1+1.0))*(qb/(jb+1.0))*((j0+1.0-q0)/(j0+1.0));
  } else {
    jm = 2*j1;
    r = 0.0;
    np = 3;
    nt = 1;
    ni = 0;
    while (1) {
      nleft = 0;
      id = 0;      
      for (j = 0; j <= jm; j += 4) {
	jmin = abs(j-j0);
	jmax = j+j0;
	for (jf = jmin; jf <= jmax; jf += 2) {
	  for (ik = -1; ik <= 1; ik += 2) {
	    klf = jf + ik;
	    kappaf = GetKappaFromJL(jf, klf);
	    kmin = abs(j0-j1);
	    kmax = j0 + j1;
	    a = 0.0;
	    for (k = kmin; k <= kmax; k += 2) {
	      if (!Triangle(jb, jf, k)) continue;
	      b = W6j(j, jf, j0, k, j1, j1);
	      if (fabs(b) < EPS30) continue;
	      if (id >= nd) {
		if (ni == 0) {
		  tp = AIRadialPk(&ai_pk, *e, k0, k1, kb, kappaf, k, 0);
		} else {
		  continue;
		}
	      } else {
		if (done[id]) {
		  id++;
		  continue;
		}
		tp = AIRadialPk(&ai_pk, *e, k0, k1, kb, kappaf, k, 1);
		if (tp == -9999) {
		  id++;
		  nleft++;
		  continue;
		}
		id++;
		done[id] = 1;
	      }
	      if (n_egrid > 1) {
		UVIP3P(np, n_egrid, log_egrid, ai_pk, nt, &log_e, &s);
	      } else {
		s = ai_pk[0];
	      } 
	      a += b*s;
	    }
	    r += a*a*2.0*(j+1.0);
	  }
	}
      }
      if (nleft == 0) break;
      ni++;
    }
    r *= 4.0*(q1/(j1+1.0))*((qb-1.0)/jb)*((j0+1.0-q0)/(j0+1.0));
  }

  *rate = r;
  
  free(idatum->bra);
  free(idatum);
  return 0;
}

int AutoionizeRateUTA(double *rate, double *e, int rec, int f) {
  LEVEL *lev1, *lev2;
  SYMMETRY *sym;
  STATE *s;
  CONFIG *c;
  int t, j, p, r;
  double wb, rq;

  *rate = 0.0;
  lev1 = GetLevel(rec);
  lev2 = GetLevel(f);

  if (lev1->n_basis > 0 && lev2->n_basis > 0) {
    return AutoionizeRate(rate, e, rec, f, 0);
  }

  if (lev1->n_basis == 0 && lev2->n_basis == 0) {
    r = AutoionizeRateUTA0(rate, e, lev2->iham, lev1->iham,
			   lev2->pb, lev1->pb);
    return r;
  } else if (lev1->n_basis > 0) {
    sym = GetSymmetry(lev1->pj);
    DecodePJ(lev1->pj, &p, &j);
    for (t = 0; t < lev1->n_basis; t++) {
      s = (STATE *) ArrayGet(&(sym->states), lev1->basis[t]);
      rq = 0.0;
      r = AutoionizeRateUTA0(&rq, e, lev2->iham, s->kgroup,
			     lev2->pb, s->kcfg);
      if (r < 0) continue;
      c = GetConfig(s);
      wb = lev1->mixing[t]*lev1->mixing[t];
      wb *= (j+1.0)/fabs(c->sweight);
      *rate += wb*rq;
    }
  } else if (lev2->n_basis > 0) {
    sym = GetSymmetry(lev2->pj);
    DecodePJ(lev2->pj, &p, &j);
    for (t = 0; t < lev2->n_basis; t++) {
      s = (STATE *) ArrayGet(&(sym->states), lev2->basis[t]);
      rq = 0.0;
      r = AutoionizeRateUTA0(&rq, e, s->kgroup, lev1->iham,
			     s->kcfg, lev1->pb);
      if (r < 0) continue;
      c = GetConfig(s);
      wb = lev2->mixing[t]*lev2->mixing[t];
      wb *= (j+1.0)/fabs(c->sweight);
      *rate += wb*rq;
    }
  }
  
  return 0;
}

int AutoionizeRateUTA1(double *rate, double *e, int rec, int f) {
  INTERACT_DATUM *idatum;
  LEVEL *lev1, *lev2;
  int j0, j1, jb, ns, q0, q1, qb;
  int k0, k1, kb, kmin, kmax, jmin, jmax;
  int jf, ik, klf, kappaf, k, np, nt, j, jm;
  double a, b, r, s, log_e, *ai_pk;
  
  *rate = 0.0;
  lev1 = GetLevel(rec);
  lev2 = GetLevel(f);
  
  log_e = log(*e);
  
  idatum = NULL;
  ns = GetInteract(&idatum, NULL, NULL, lev2->iham, lev1->iham,
		   lev2->pb, lev1->pb, 0, 0, 1);  
  if (ns <= 0) return -1;
  if (idatum->s[3].index < 0) {
    free(idatum->bra);
    free(idatum);
    return -1;
  }

  kb = OrbitalIndex(idatum->s[1].n, idatum->s[1].kappa, 0.0);
  k0 = OrbitalIndex(idatum->s[2].n, idatum->s[2].kappa, 0.0);
  k1 = OrbitalIndex(idatum->s[3].n, idatum->s[3].kappa, 0.0);
  j0 = idatum->s[2].j;
  j1 = idatum->s[3].j;
  jb = idatum->s[1].j;
  q0 = idatum->s[2].nq_ket;
  q1 = idatum->s[3].nq_ket;
  qb = idatum->s[1].nq_ket;

  if (idatum->s[1].index != idatum->s[3].index) {
    kmin = abs(j0-j1);
    kmax = j0 + j1;
    jmin = 1;
    jmax = j1+j0+jb;
    np = 3;
    nt = 1;
    r = 0.0;
    for (jf = jmin; jf <= jmax; jf += 2) {
      for (ik = -1; ik <= 1; ik += 2) {
	klf = jf + ik;
	kappaf = GetKappaFromJL(jf, klf);
	for (k = kmin; k <= kmax; k += 2) {
	  if (!Triangle(j0, j1, k) || !Triangle(jb, jf, k)) continue;
	  AIRadialPk(&ai_pk, *e, k0, k1, kb, kappaf, k, 0);
	  if (n_egrid > 1) {
	    UVIP3P(np, n_egrid, log_egrid, ai_pk, nt, &log_e, &s);
	  } else {
	    s = ai_pk[0];
	  }
	  s = s*s/(k + 1.0);
	  r += s;
	}
      }
    }  
    r *= 4.0*(q1/(j1+1.0))*(qb/(jb+1.0))*((j0+1.0-q0)/(j0+1.0));
  } else {
    jm = 2*j1;
    r = 0.0;
    np = 3;
    nt = 1;
    for (j = 0; j <= jm; j += 4) {
      jmin = abs(j-j0);
      jmax = j+j0;
      for (jf = jmin; jf <= jmax; jf += 2) {
	for (ik = -1; ik <= 1; ik += 2) {
	  klf = jf + ik;
	  kappaf = GetKappaFromJL(jf, klf);
	  kmin = abs(j0-j1);
	  kmax = j0 + j1;
	  a = 0.0;
	  for (k = kmin; k <= kmax; k += 2) {
	    if (!Triangle(jb, jf, k)) continue;
	    b = W6j(j, jf, j0, k, j1, j1);
	    if (fabs(b) < EPS30) continue;
	    AIRadialPk(&ai_pk, *e, k0, k1, kb, kappaf, k, 0);
	    if (n_egrid > 1) {
	      UVIP3P(np, n_egrid, log_egrid, ai_pk, nt, &log_e, &s);
	    } else {
	      s = ai_pk[0];
	    } 
	    a += b*s;
	  }
	  r += a*a*2.0*(j+1.0);
	}
      }
    }
    r *= 4.0*(q1/(j1+1.0))*((qb-1.0)/jb)*((j0+1.0-q0)/(j0+1.0));
  }

  *rate = r;
  
  free(idatum->bra);
  free(idatum);
  return 0;
}

int AutoionizeRate(double *rate, double *e, int rec, int f, int msub) {  
  LEVEL *lev1, *lev2;
  ANGULAR_ZxZMIX *ang;
  ANGULAR_ZFB *zfb;
  STATE *st;
  int k, nz, nzfb, ik, i, j1, j2, ij, kappaf, ip;
  int jf, k0, k1, kb, njf, nkappaf, klf, jmin, jmax;
  double *p, r, s, log_e, a, pe;
  double *ai_pk, ai_pk0[MAXRRNE];
  int np, nt, m1, m2, m;
  int kappafp, jfp, klfp, dkl;

  /*
  nz = aicache.nz[ic];
  nzfb = aicache.nzf[ic];
  ang = aicache.az[ic];
  zfb = aicache.azf[ic];
  if (nz <= 0 && nzfb <= 0) return -1;
  */
  *rate = 0.0;
  lev1 = GetLevel(rec);
  lev2 = GetLevel(f);
  if (uta_egrid) pe = *e;
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
    //int iter = 0;
    short *ia;
    ia = malloc(sizeof(short)*nz);
    for (i = 0; i < nz; i++) {
      ia[i] = 0;
    }
    while (1) {
      //iter++;
      int nleft = 0;      
      for (i = 0; i < nz; i++) {
	if (ia[i] == 3) continue;
	jf = ang[i].k0;
	kb = ang[i].k1;
	k0 = ang[i].k2;
	k1 = ang[i].k3;
	kappafp = GetOrbital(kb)->kappa;
	klfp = GetLFromKappa(kappafp);
	kappafp = GetOrbital(k0)->kappa;
	klfp += GetLFromKappa(kappafp);
	kappafp = GetOrbital(k1)->kappa;
	klfp += GetLFromKappa(kappafp);      
	ij = (jf - jmin);
	for (ik = -1; ik <= 1; ik += 2) {
	  if (ik == -1) {
	    if (ia[i] == 1) continue;
	  } else {
	    if (ia[i] == 2) continue;
	  }
	  klf = jf + ik;  	
	  if (IsOdd((klfp+klf)/2)) {
	    if (ik == -1) ia[i] |= 1;
	    else ia[i] |= 2;
	    continue;
	  }
	  kappaf = GetKappaFromJL(jf, klf);
	  int type = AIRadialPk(&ai_pk, pe, k0, k1, kb, kappaf, ang[i].k, 1);
	  if (type == -9999) {
	    nleft++;
	    continue;
	  }
	  if (ik == -1) {
	    ia[i] |= 1;
	  } else {
	    ia[i] |= 2;
	  }
	  if (n_egrid > 1) {
	    UVIP3P(np, n_egrid, log_egrid, ai_pk, nt, &log_e, &s);
	  } else {
	    s = ai_pk[0];
	  }
	  ip = (ik == -1)? ij:(ij+1);
	  p[ip] += s*ang[i].coeff;
	}
      }
      if (nleft == 0) break;
    }
    free(ang);
    //aicache.az[ic] = NULL;
    //aicache.nz[ic] = 0;
    free(ia);
  }
  
  nzfb = AngularZFreeBound(&zfb, f, rec);
  if (nzfb > 0) {
    for (i = 0; i < nzfb; i++) {
      kb = zfb[i].kb;
      kappaf = GetOrbital(kb)->kappa;
      GetJLFromKappa(kappaf, &jf, &klf);
      ij = jf - jmin;
      ik = klf - jf;
      AIRadial1E(ai_pk0, kb, kappaf, *e);
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
    free(zfb);
    //aicache.azf[ic] = NULL;
    //aicache.nzf[ic] = 0;
  }

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
		if (!uta_egrid) pe = egrid[ik];
		k0 = OrbitalIndex(0, kappaf, pe);
		k1 = OrbitalIndex(0, kappafp, pe);
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

int AIRadial1E(double *ai_pk, int kb, int kappaf, double pe) {
  int kf;
  int i;

  for (i = 0; i < n_egrid; i++) {
    if (!uta_egrid) pe = egrid[i];
    kf = OrbitalIndex(0, kappaf, pe);
    ResidualPotential(ai_pk+i, kf, kb);
  }
  return 0;
}  

int AIRadialPk(double **ai_pk, double pe, int k0, int k1, int kb, int kappaf,
	       int k, int trylock) {
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
  LOCK *lock = NULL;
  int locked = 0;
  p = (double **) MultiSet(pk_array, index, NULL, &lock,
			   InitPointerData, FreeRecPkData);
  if (lock && !(*p)) {
    if (trylock) {
      if (0 == TryLock(lock)) {
	locked = 1;
      } else {
	return -9999;
      }
    } else {
      SetLock(lock);
      locked = 1;
    }
  }
  if (*p) {
    *ai_pk = *p;
    if (locked) ReleaseLock(lock);
    return 0;
  } 
  double *pd = (double *) malloc(sizeof(double)*n_egrid);
  *ai_pk = pd;
  if (uta_egrid) e = pe;
  for (i = 0; i < n_egrid; i++) {
    if (!uta_egrid) e = egrid[i];
    kf = OrbitalIndex(0, kappaf, e);
    ks[0] = k0;
    ks[1] = kf;
    ks[2] = k1;
    ks[3] = kb;
    SlaterTotal(&sd, &se, NULL, ks, k, 0);
    (*ai_pk)[i] = sd+se;
  }
  *p = pd;
  if (locked) ReleaseLock(lock);
#pragma omp flush
  return 0;
}

void ProcessAICache(int msub, int iuta, TFILE *f) {
  int ic, iz, ilow, iup, skip;

  if (!iuta) {
    ResetWidMPI();
#pragma omp parallel default(shared) private(ic, ilow, iup, skip)
    {
      for (ic = 0; ic < aicache.nc; ic++) {
	ilow = aicache.low[ic];
	iup = aicache.up[ic];
	skip = SkipMPI();
	if (skip) continue;
	aicache.nz[ic] = AngularZxZFreeBound(&aicache.az[ic], iup, ilow);
	aicache.nzf[ic] = AngularZFreeBound(&aicache.azf[ic], iup, ilow);
      }
    }
    /*
    ResetWidMPI();
#pragma omp parallel default(shared) private(ic, iz, skip)
    {
      int jf, kb, k0, k1, kappafp, klfp, ij, ik, type, kappaf, klf;
      double *ai_pk;
      for (ic = 0; ic < aicache.nc; ic++) {
	skip = SkipMPI();
	if (skip) continue;
	LEVEL *lev1 = GetLevel(aicache.low[ic]);
	LEVEL *lev2 = GetLevel(aicache.up[ic]);
	int j1 = lev1->pj;
	int j2 = lev2->pj;
	DecodePJ(j1, NULL, &j1);
	DecodePJ(j2, NULL, &j2);
	int jmin = abs(j1-j2);
	for (iz = 0; iz < aicache.nz[ic]; iz++) {
	  jf = aicache.az[ic][iz].k0;
	  kb = aicache.az[ic][iz].k1;
	  k0 = aicache.az[ic][iz].k2;
	  k1 = aicache.az[ic][iz].k3;
	  kappafp = GetOrbital(kb)->kappa;
	  klfp = GetLFromKappa(kappafp);
	  kappafp = GetOrbital(k0)->kappa;
	  klfp += GetLFromKappa(kappafp);
	  kappafp = GetOrbital(k1)->kappa;
	  klfp += GetLFromKappa(kappafp);
	  ij = jf-jmin;
	  for (ik = -1; ik <= 1; ik += 2) {
	    klf = jf + ik;
	    if (IsOdd((klfp+klf)/2)) continue;
	    kappaf = GetKappaFromJL(jf, klf);
	    type = AIRadialPk(&ai_pk, 0.0, k0, k1, kb, kappaf,
			      aicache.az[ic][iz].k, 1);
	  }
	}
      }
      MPrintf(-1, "tpk0: %g %g\n", tpk0, tpk1);
    }
    wt1=WallTime();
    printf("aip: %g %g %g\n", wt1-wt0, tpk0, tpk1);
    wt0=wt1;
    ResetWidMPI();
#pragma omp parallel default(shared) private(ic, iz, skip)
    {
      double ai_pk0[MAXNE];
      int kb, kappaf;
      for (ic = 0; ic < aicache.nc; ic++) {
	skip = SkipMPI();
	if (skip) continue;
	for (iz = 0; iz < aicache.nzf[ic]; iz++) {
	  kb = aicache.azf[ic][iz].kb;
	  kappaf = GetOrbital(kb)->kappa;
	  AIRadial1E(ai_pk0, kb, kappaf, 0.0);
	}
      }
    }
    wt1 = WallTime();
    printf("aie: %g\n", wt1-wt0);
    */
  }    
  ResetWidMPI();
#pragma omp parallel default(shared) private(ic, ilow, iup, skip)
  {
    float rt[MAXAIM];
    double e, s[MAXAIM], s1;
    int k, t;
    AI_RECORD r;
    AIM_RECORD r1;
    for (ic = 0; ic < aicache.nc; ic++) {
      skip = SkipMPI();
      if (skip) continue;
      ilow = aicache.low[ic];
      iup = aicache.up[ic];
      e = aicache.e[ic];
      if (iuta) {
	k = AutoionizeRateUTA(s, &e, ilow, iup);
      } else {
	k = AutoionizeRate(s, &e, ilow, iup, msub);
      }
      if (k < 0) continue;
      if (!msub) {
	if (s[0] <= _ai_cut) continue;
	r.b = ilow;
	r.f = iup;
	r.rate = s[0];
	WriteAIRecord(f, &r);
      } else {
	r1.rate = rt;
	s1 = 0;
	for (t = 0; t < k; t++) {
	  r1.rate[t] = s[t];
	  s1 += s[t];
	}
	if (s1 <= _ai_cut) continue;
	r1.b = ilow;
	r1.f = iup;
	r1.nsub = k;
	WriteAIMRecord(f, &r1);
      }
    }
  }
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

  if (n_egrid == 0) {
    n_egrid = 6;
  }
  if (egrid[0] < 0.0) {
    SetPEGrid(n_egrid, emin, emax, e);
  }
  if (n_usr <= 0) {
    SetUsrPEGridDetail(n_egrid, egrid);
    usr_egrid_type = 1;
  } else if (usr_egrid[0] < 0) {
    SetUsrPEGrid(n_usr, emin, emax, e);
    usr_egrid_type = 1;
  }

  if (qk_mode == QK_INTERPOLATE) {
    for (j = 0; j < n_egrid; j++) {
      log_egrid[j] = egrid[j];      
      if (egrid_type == 1) {
	if (uta_tegrid) {
	  log_egrid[j] += 1.0;
	} else {
	  log_egrid[j] += e;
	}
      }
      log_egrid[j] = log(log_egrid[j]);
    }
    for (j = 0; j < n_usr; j++) {
      log_usr[j] = usr_egrid[j];
      if (usr_egrid_type == 1) {
	if (uta_tegrid) {
	  log_usr[j] += 1.0;
	} else {
	  log_usr[j] += e;
	}
      }
      log_usr[j] = log(log_usr[j]);
    }
  }
  return 0;
}

int SaveRRMultipole(int nlow, int *low, int nup, int *up, char *fn, int m) {
  int i, j, k;
  FILE *f;
  LEVEL *lev1, *lev2;
  double e, emin, emax, emax0;
  double awmin, awmax;
  ARRAY subte;
  int isub, n_tegrid0, n_egrid0, n_usr0;
  int te_set, e_set, usr_set;
  double c, e0, e1;

  f = fopen(fn, "w");
  if (f == NULL) return -1;
  fprintf(f, "# This file contains 11 columns\n");
  fprintf(f, "#  1, Bound state index\n");
  fprintf(f, "#  2, 2*J_Bound\n");
  fprintf(f, "#  3, Ionized state index\n");
  fprintf(f, "#  4, 2*J_Ionized\n");
  fprintf(f, "#  5, Multipole type, -1=E1, 1=M1, -2=E2, 2=M2, ...\n");
  fprintf(f, "#  6, 2*J_total, coupled angular momentum of J_Ionized and J\n");
  fprintf(f, "#  7, 2*L, L is the orbital angular momentum of the photo-electron\n");
  fprintf(f, "#  8, 2*J, J is the total angular momentum of the photo-electron\n");
  fprintf(f, "#  9, E_th, Ionization threshold in Hartree\n");
  fprintf(f, "# 10, E_e, photo-electron energy in Hartree\n");
  fprintf(f, "# 11, Multipole matrix element in atomic unit\n");
  fprintf(f, "\n\n");
  
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
    if (e < 0.1) {
      SetAWGrid(1, 0.5*(awmin+awmax), awmax);
    } else if (e < 0.5) {
      SetAWGrid(2, awmin, awmax);
    } else {
      SetAWGrid(3, awmin, awmax);
    }
    
    n_egrid = n_egrid0;
    n_usr = n_usr0;
    if (!usr_set) usr_egrid[0] = -1.0;
    if (!e_set) egrid[0] = -1.0;
    e = 0.5*(emin + emax);
    PrepRREGrids(e, emax0);
    
    for (i = 0; i < nup; i++) {
      lev1 = GetLevel(up[i]);
      for (j = 0; j < nlow; j++) {
	lev2 = GetLevel(low[j]);
	e = lev1->energy - lev2->energy;
	if (e < e0 || e >= e1) continue;
	BoundFreeMultipole(f, low[j], up[i], m);
      }
    }
    
    ReinitRadial(1);
    FreeRecQk();
    FreeRecPk();
    
    e0 = e1;
  }
  
  fclose(f);
  ArrayFreeLock(&subte, NULL);
  ReinitRecombination(1);

  return 0;
}

int SaveRecOccupation(int nlow, int *low, int nup, int *up, char *fn) {
  int i, j, n;
  double eb;
  RO_RECORD r;
  RO_HEADER roc_hdr;
  F_HEADER fhdr;
  TFILE *f;

  fhdr.type = DB_RO;
  strcpy(fhdr.symbol, GetAtomicSymbol());
  fhdr.atom = GetAtomicNumber();
  roc_hdr.nele = GetNumElectrons(low[0]);
  f = OpenFile(fn, &fhdr);
  InitFile(f, &fhdr, &roc_hdr);
  ResetWidMPI();
#pragma omp parallel default(shared) private(r, i, j, eb, n)
  {
    for (j = 0; j < nup; j++) {
      for (i = 0; i < nlow; i++) {
	int skip = SkipMPI();
	if (skip) continue;
	n = RecOccupation(&r.nk, &r.nq, &r.dn, &eb, low[i], up[j]);
	if (n <= 0) continue;
	r.b = low[i];
	r.f = up[j];
	r.n = n;
	WriteRORecord(f, &r);	
	free(r.nk);
	free(r.nq);
	free(r.dn);
      }
    }
  }
  DeinitFile(f, &fhdr);
  CloseFile(f, &fhdr);
  return 0;
}

int SetCXEGrid(int n, double e0, double e1, double *eg, int ilog, int t) {
  double de;
  int i;
  
  if (n_cxegrid > 0) {
    free(cxegrid0);
    free(cxegrid);
    n_cxegrid = 0;    
    cxegrid = NULL;
  }
  n_cxegrid = n;
  t_cxegrid = t;
  cxegrid = malloc(sizeof(double)*n);
  cxegrid0 = malloc(sizeof(double)*n);
  if (eg != NULL) {
    for (i = 0; i < n; i++) {
      cxegrid0[i] = eg[i]/HARTREE_EV;
      cxegrid[i] = cxegrid0[i];
    }
    return 0;
  }
  e0 /= HARTREE_EV;
  e1 /= HARTREE_EV;
  if (ilog) {
    e0 = log(e0);
    e1 = log(e1);
  }
  de = (e1-e0)/(n-1);
  cxegrid0[0] = e0;
  for (i = 1; i < n; i++) {
    cxegrid0[i] = cxegrid0[i-1] + de;
  }
  if (ilog) {
    for (i = 0; i < n; i++) {
      cxegrid0[i] = exp(cxegrid0[i]);
    }
  }
  for (i = 0; i < n; i++) {
    cxegrid[i] = cxegrid0[i];
  }
  return 0;
}

int SaveCX(int nlow, int *low, int nup, int *up, char *fn) {
  int i, j, n;
  double eb, pmass;
  CX_RECORD r;
  CX_HEADER cx_hdr;
  F_HEADER fhdr;
  TFILE *f;
  CXTGT *cxt = GetCXTarget();

  if (cxt == NULL || cxt->z <= 0) {
    printf("CX target not setup\n");
    return -1;
  }
  if (t_cxegrid < 0 || cxegrid == NULL) {
    printf("CX egrid not setup: %d\n", t_cxegrid);
    return -1;
  }
  pmass = GetAtomicMass();
  for (i = 0; i < n_cxegrid; i++) {
    cxegrid[i] = cxegrid0[i];
    if (t_cxegrid > 0) cxegrid[i] /= pmass;
  }
  double z = GetResidualZ();
  int nmc;
  double **mc = LandauZenerMC(cxt, z, n_cxegrid, cxegrid, &nmc);
  fhdr.type = DB_CX;
  strcpy(fhdr.symbol, GetAtomicSymbol());
  fhdr.atom = GetAtomicNumber();
  cx_hdr.nele = GetNumElectrons(low[0]);
  cx_hdr.te0 = t_cxegrid;
  cx_hdr.ne0 = n_cxegrid;
  cx_hdr.e0 = cxegrid0;  
  memcpy(cx_hdr.tgts, cxt->symbol, 128);
  cx_hdr.tgtz = cxt->z;
  cx_hdr.tgtm = cxt->m;
  cx_hdr.tgta = cxt->a;
  cx_hdr.tgtb = cxt->b;
  cx_hdr.tgte = cxt->e;
  cx_hdr.tgtx = cxt->x;
  cx_hdr.ldist = cxldist;
  f = OpenFile(fn, &fhdr);
  InitFile(f, &fhdr, &cx_hdr);
  ResetWidMPI();
#pragma omp parallel default(shared) private(r, i, j, eb, n)
  {
    r.cx = malloc(sizeof(double)*n_cxegrid);
    for (j = 0; j < nup; j++) {
      for (i = 0; i < nlow; i++) {
	int skip = SkipMPI();
	if (skip) continue;
	n = CXCross(cxt, r.cx, &eb, low[i], up[j], nmc, mc);
	if (n <= 0) continue;
	r.b = low[i];
	r.f = up[j];
	r.vnl = n;
	WriteCXRecord(f, &r);
      }
    }
    free(r.cx);
  }
  DeinitFile(f, &fhdr);
  CloseFile(f, &fhdr);
  for (i = 0; i < nmc; i++) {
    free(mc[i]);
  }
  free(mc);
  return 0;
}

int SaveRecRR(int nlow, int *low, int nup, int *up, 
	      char *fn, int m0) {
  int i, j, k, ie, ip;
  TFILE *f;
  double rqu[MAXRRNUSR], qc[NPARAMS+1];
  double eb;
  LEVEL *lev1, *lev2;
  RR_RECORD r;
  RR_HEADER rr_hdr;
  F_HEADER fhdr;
  double e, emin, emax, emax0;
  double awmin, awmax;
  int nq, nqk, qkmo;
  ARRAY subte;
  int isub, n_tegrid0, n_egrid0, n_usr0;
  int te_set, e_set, usr_set, iuta, auta;
  double c, e0, e1;
  RANDIDX *rid0, *rid1;

  double tstart = WallTime();

  iuta = IsUTA();
  auta = iuta && UTAGrid();
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
      if (lev1->n_basis > 0 || lev2->n_basis > 0) auta = 0;
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
  if (!auta) {
    c = 1.0/TE_MIN_MAX;
    if (!e_set || !te_set) {
      e = c*emin;
      while (e < emax) {
	ArrayAppend(&subte, &e, NULL);
	e *= c;
      }
    }
  }
  ArrayAppend(&subte, &emax, NULL);

  fhdr.type = DB_RR;
  strcpy(fhdr.symbol, GetAtomicSymbol());
  fhdr.atom = GetAtomicNumber();
  rr_hdr.nele = GetNumElectrons(low[0]);
  rr_hdr.multipole = m0;
  f = OpenFile(fn, &fhdr);

  if (NProcMPI() > 1) {
    rid0 = RandList(nlow);
    rid1 = RandList(nup);
  } else {
    rid0 = NULL;
    rid1 = NULL;
  }
  int nproc = 0, *ntrans = NULL;
  if (_progress_report != 0) {
    ntrans = InitTransReport(&nproc);
  }
  e0 = emin*0.999;
  for (isub = 1; isub < subte.dim; isub++) {
    e1 = *((double *) ArrayGet(&subte, isub));
    if (isub == subte.dim-1) e1 = e1*1.001;
    if (auta) {
      emin = e0;
      emax = e1;
    } else {
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
    }
    e = (emax - emin)/(0.5*(emin+emax));
    if (auta) {
      SetRRTEGrid(-1, emin, emax);
    } else {
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
    }
    FreeMultipoleArray();
    awmin = emin * FINE_STRUCTURE_CONST;
    awmax = emax * FINE_STRUCTURE_CONST;
    if (auta) {
      SetAWGrid(-1, awmin, awmax);
    } else {
      if (e < 0.1) {
	SetAWGrid(1, 0.5*(awmin+awmax), awmax);
      } else if (e < 0.5) {
	SetAWGrid(2, awmin, awmax);
      } else {
	SetAWGrid(3, awmin, awmax);
      }
    }
    n_egrid = n_egrid0;
    n_usr = n_usr0;
    if (!usr_set) usr_egrid[0] = -1.0;
    if (!e_set) egrid[0] = -1.0;
    e = 0.5*(emin + emax);
    PrepRREGrids(e, emax0);
    
    qkmo = qk_mode;
    if (qk_mode == QK_FIT && n_egrid <= NPARAMS) {
      qk_mode = QK_EXACT;
    }
    if (qk_mode == QK_FIT) {
      nqk = NPARAMS+1;
    } else {
      nqk = 0;
    }
    rr_hdr.qk_mode = qk_mode;
    rr_hdr.nparams = nqk;
    rr_hdr.n_tegrid = n_tegrid;
    rr_hdr.tegrid = tegrid;
    rr_hdr.n_egrid = n_egrid;
    rr_hdr.egrid = egrid;
    rr_hdr.n_usr = n_usr;
    rr_hdr.usr_egrid = usr_egrid;
    rr_hdr.egrid_type = egrid_type;
    rr_hdr.usr_egrid_type = usr_egrid_type;
        
    InitFile(f, &fhdr, &rr_hdr);
    ResetWidMPI();
#pragma omp parallel default(shared) private(r, i, j, lev1, lev2, e, nq, ip, ie, rqu, qc, eb)
    {
    if (qk_mode == QK_FIT) {
      r.params = (float *) malloc(sizeof(float)*nqk);
    }
    r.strength = (float *) malloc(sizeof(float)*n_usr);
    int myrank = MPIRank(NULL);
    int ilow, iup;
    int ipr = 0;
    for (i = 0; i < nup; i++) {
      if (rid1) iup = up[rid1[i].i];
      else iup = up[i];
      lev1 = GetLevel(iup);
      for (j = 0; j < nlow; j++) {
	if (rid0) ilow = low[rid0[j].i];
	else ilow = low[j];
	lev2 = GetLevel(ilow);
	e = lev1->energy - lev2->energy;
	if (e < e0 || e >= e1) continue;
	int skip = SkipMPI();
	if (skip) continue;
	if (iuta) {
	  nq = BoundFreeOSUTA(rqu, qc, &eb, ilow, iup, m0);
	} else {
	  nq = BoundFreeOS(rqu, qc, &eb, ilow, iup, m0);
	}
	if (nq < 0) continue;
	r.b = ilow;
	r.f = iup;
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
	if (ntrans) {
	  ntrans[myrank]++;
	  if (_progress_report > 0) {
	    if (myrank == 0 && ntrans[0]%_progress_report == 0) {
	      PrintTransReport(nproc, tstart, ntrans, "RR", ipr++);
	    }
	  } else if (_progress_report < 0) {
	    if (ntrans[myrank]%(-_progress_report) == 0) {
	      PrintTransReport(-myrank, tstart, ntrans, "RR", ipr++);
	    }
	  }
	}
      }
    }
    if (qk_mode == QK_FIT) {
      free(r.params);
    }
    free(r.strength);
    }
    
    DeinitFile(f, &fhdr);    
    ReinitRadial(1);
    FreeRecQk();
    FreeRecPk();
    
    e0 = e1;
    qk_mode = qkmo;
  }
      
  ReinitRecombination(1);
  if (rid0) free(rid0);
  if (rid1) free(rid1);
  ArrayFreeLock(&subte, NULL);
  CloseFile(f, &fhdr);
  if (_progress_report != 0) {
    PrintTransReport(nproc, tstart, ntrans, "RR", -1);
  }

  return 0;
}

int SaveAI(int nlow, int *low, int nup, int *up, char *fn, 
	   double eref, int msub) {
#ifdef PERFORM_STATISTICS
  STRUCT_TIMING structt;
  ANGULAR_TIMING angt;
  RECOUPLE_TIMING recouplet;
  RAD_TIMING radt;
#endif
  int i, j, k, t;
  LEVEL *lev1, *lev2;
  AI_RECORD r;
  AIM_RECORD r1;
  AI_HEADER ai_hdr;
  AIM_HEADER ai_hdr1;
  F_HEADER fhdr;
  double emin, emax, e, a, s, s1[MAXAIM];
  float rt[MAXAIM];
  TFILE *f;
  ARRAY subte;
  double c, e0, e1, b;
  int isub, n_egrid0;
  int e_set, iuta, auta;
  RANDIDX *rid0, *rid1;

  double tstart = WallTime();
  
  iuta = IsUTA();
  auta = iuta && UTAGrid();
  if (iuta && msub) {
    printf("cannot call AITableMSub with UTA mode\n");
    return -1;
  }

  if (nup <= 0 || nlow <= 0) return -1;

  eref /= HARTREE_EV;
  emin = 1E10;
  emax = 1E-10;
  k = 0;
  for (i = 0; i < nlow; i++) {
    lev1 = GetLevel(low[i]);
    for (j = 0; j < nup; j++) {
      lev2 = GetLevel(up[j]);
      a = lev1->energy - lev2->energy;
      if (a < 0) a -= eref;
      if (a > 0) k++;
      if (a < emin && a > 0) emin = a;
      if (a > emax) emax = a;
      if (lev1->n_basis > 0 || lev2->n_basis > 0) auta = 0;
    }
  }
  if (k == 0) {
    return 0;
  }
  /*
  if (!iuta) {
    AllocAICache();
  }
  */
  if (egrid[0] < 0) {
    e_set = 0;
  } else {
    e_set = 1;
  }

  n_egrid0 = n_egrid;

  ArrayInit(&subte, sizeof(double), 128);
  ArrayAppend(&subte, &emin, NULL);
  if (!auta) {
    c = 1.0/TE_MIN_MAX;
    if (!e_set) {
      b = c*emin;
      while (b < emax) {
	ArrayAppend(&subte, &b, NULL);
	b *= c;
      }
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
    ai_hdr.emin = eref;
  } else {
    ai_hdr1.nele = GetNumElectrons(low[0]);
    ai_hdr1.emin = eref;
  }
  f = OpenFile(fn, &fhdr);

  if (NProcMPI() > 1) {
    rid0 = RandList(nlow);
    rid1 = RandList(nup);
  } else {
    rid0 = NULL;
    rid1 = NULL;
  }
  int nproc = 0;
  int *ntrans = NULL;
  if (_progress_report != 0) {
    ntrans = InitTransReport(&nproc);
  }
  e0 = emin*0.999;
  for (isub = 1; isub < subte.dim; isub++) {
    e1 = *((double *) ArrayGet(&subte, isub));
    if (isub == subte.dim-1) e1 = e1*1.001;
    if (auta) {
      emin = e0;
      emax = e1;
    } else {
      emin = e1;
      emax = e0;
      k = 0;
      for (i = 0; i < nlow; i++) {
	lev1 = GetLevel(low[i]);
	for (j = 0; j < nup; j++) {
	  lev2 = GetLevel(up[j]);
	  a = lev1->energy - lev2->energy;
	  if (a < 0) a -= eref;
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
    }
    if (auta) {
      SetPEGrid(-1, a, a, 0.0);
    } else {
      if (!e_set) {
	a = (emax-emin)/(0.5*(emax+emin));
	if (a < EPS3) {
	  a = 0.5*(emin+emax);
	  SetPEGrid(1, a, a, 0.0);
	} else if (a < 0.4) {
	  SetPEGrid(2, emin, emax, 0.0);
	} else if (a < 1.0) {
	  if (k == 2) n_egrid = 2;
	  else if (n_egrid0 == 0) n_egrid = 3;
	  SetPEGrid(n_egrid, emin, emax, 0.0);
	} else {
	  if (k == 2) n_egrid = 2;
	  else if (n_egrid0 == 0) n_egrid = 4;
	  SetPEGrid(n_egrid, emin, emax, 0.0);
	}
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
    ResetWidMPI();
#pragma omp parallel default(shared) private(i, j, lev1, lev2, e, k, s, r, r1, s1, t, rt)
    {
    int myrank = MPIRank(NULL);
    int ilow, iup;
    int ipr = 0;
    for (i = 0; i < nlow; i++) {
      if (rid0) ilow = low[rid0[i].i];
      else ilow = low[i];
      lev1 = GetLevel(ilow);
      for (j = 0; j < nup; j++) {
	if (rid1) iup = up[rid1[j].i];
	else iup = up[j];
	lev2 = GetLevel(iup);	
	int skip = SkipMPI();
	if (skip) continue;
	e = lev1->energy - lev2->energy;
	if (e < 0 && lev1->ibase != iup) e -= eref;
	if (e < e0 || e >= e1) continue;
	
	if (!msub) {
	  if (iuta) {
	    k = AutoionizeRateUTA(&s, &e, ilow, iup);
	  } else {
	    k = AutoionizeRate(&s, &e, ilow, iup, msub);
	  }
	  if (k < 0) continue;
	  if (s < _ai_cut) continue;
	  r.b = ilow;
	  r.f = iup;
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
	  if (s < _ai_cut) continue;
	  r1.b = ilow;
	  r1.f = iup;
	  r1.nsub = k;
	  WriteAIMRecord(f, &r1);
	}
	if (ntrans) {
	  ntrans[myrank]++;
	  if (_progress_report > 0) {
	    if ( myrank == 0 && ntrans[0]%_progress_report == 0) {
	      PrintTransReport(nproc, tstart, ntrans, "AI", ipr++);
	    }
	  } else if (_progress_report < 0) {
	    if (ntrans[myrank]%(-_progress_report) == 0) {
	      PrintTransReport(-myrank, tstart, ntrans, "AI", ipr++);
	    }
	  }
	} 
      }
    }
    }
    /*
      aicache.low[ic] = ilow;
      aicache.up[ic] = iup;
      aicache.e[ic] = e;
      ic++;
      if (ic == maxaicache) {
      aicache.nc = ic;
      ic = 0;
      ProcessAICache(msub, iuta, f);
      }
      }
      }
      if (ic > 0) {
      aicache.nc = ic;
      ProcessAICache(msub, iuta, f);
      }
    */
    DeinitFile(f, &fhdr);
    ReinitRadial(1);
    FreeRecQk();
    FreeRecPk();

    e0 = e1;
  }

  if (rid0) free(rid0);
  if (rid1) free(rid1);
  if (_progress_report != 0) {
    PrintTransReport(nproc, tstart, ntrans, "AI", -1);
  }
  ReinitRecombination(1);
  //FreeAICache(0);
  ArrayFreeLock(&subte, NULL);
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
#endif /* PERFORM_STATISTICS */

  return 0;
}

double AsymmetryPI(int k0, int k0p, double e, double te,
		   int mp, int mx0, int *mx1,
		   int m, double *b, double *pqa, double *pqa2) {
  ORBITAL *orb0, *orb0p;
  double **ak, **akp;
  int *nak, **kak, gauge, ip0, ip1, im, mm;
  int L, L2, i, p, q, j0, j1, j2, kl0, kl1, kl2, mx;
  int jmin, jmax, se, kappa, k, Lp, Lp2, ip, pp, q2;
  double aw, aw0, ph1, ph2, c, d, b0, b1, w3j1, w3j2;
  double bm, bp, bmc, bpc, bmp, bmi, bpi, bmt, bpt;

  orb0 = GetOrbital(k0);
  GetJLFromKappa(orb0->kappa, &j0, &kl0);
  if (te <= 0) {
    te = -orb0->energy;
  }
  aw = FINE_STRUCTURE_CONST*te;
  aw0 = FINE_STRUCTURE_CONST*(e + te);
  gauge = GetTransitionGauge();
  
  mx = *mx1;
  nak = malloc(sizeof(int)*mx);
  ak = malloc(sizeof(double *)*mx);
  if (k0p != k0) {
    akp = malloc(sizeof(double *)*mx);
  } else {
    akp = ak;
  }
  kak = malloc(sizeof(int *)*mx);
  for (i = 0; i < m*2+1; i++) {
    b[i] = 0.0;
  }
  b1 = 0.0;
  bm = 0.0;
  bmp = 0.0;
  bmt = 0.0;
  bpt = 0.0;
  for (im = mx0; im < mx; im++) {
    L = im/2 + 1;
    L2 = 2*L;
    jmin = abs(j0-L2);
    jmax = j0 + L2;
    nak[im] = 1 + ((jmax - jmin)/2);
    ak[im] = malloc(sizeof(double)*nak[im]);
    if (k0p != k0) {
      akp[im] = malloc(sizeof(double)*nak[im]);
    }
    kak[im] = malloc(sizeof(int)*nak[im]);
    j1 = jmin;
    if (IsEven(im)) {
      if (im > mx0+2) {
	bmp = fabs(bm);
      }
      bm = 0.0;
      se = IsEven(L+kl0/2);
    } else {
      se = IsOdd(L+kl0/2);
    }
    c = sqrt((2.0/(L2+1.0))/aw0);
    if (L < MMFR) c *= pow(aw0, L);
    for (p = 0; p < nak[im]; p++) {
      kl1 = j1 - 1;
      if (se) {
	if (IsOdd(kl1/2)) kl1 += 2;
      } else {
	if (IsEven(kl1/2)) kl1 += 2;
      }
      kappa = GetKappaFromJL(j1, kl1);
      kak[im][p] = kappa;
      k = OrbitalIndex(0, kappa, e);
      if (IsEven(im)) {
	if (mp <= 0) {
	  ak[im][p] = MultipoleRadialFR(aw, -L, k0, k, gauge);
	  if (k0p != k0) {
	    akp[im][p] = MultipoleRadialFR(aw, -L, k0p, k, gauge);
	  }
	  if (IsOdd((j0+j1)/2+(L2+kl0-kl1)/4)) {
	    ak[im][p] = -ak[im][p];
	    if (k0p != k0) akp[im][p] = -akp[im][p];
	  }
	} else {
	  ak[im][p] = 0.0;
	  if (k0p != k0) akp[im][p] = 0.0;
	}
      } else {
	if (mp >= 0) {
	  ak[im][p] = MultipoleRadialFR(aw, L, k0, k, gauge);
	  if (k0p != k0) {
	    akp[im][p] = MultipoleRadialFR(aw, L, k0p, k, gauge);
	  }
	  if (IsOdd((j0+j1)/2+(L2+2+kl0-kl1)/4)) {
	    ak[im][p] = -ak[im][p];
	    if (k0p != k0) akp[im][p] = -akp[im][p];
	  }
	} else {
	  ak[im][p] = 0.0;
	  if (k0p != k0) akp[im][p] = 0.0;
	}
      }
      ak[im][p] *= c;      
      if (k0p != k0) {
	akp[im][p] *= c;
      }
      b0 = ak[im][p]*akp[im][p];
      b[0] += b0;
      bm += b0;
      if (im == mx0) {
	b1 += b0;
      }
      j1 += 2;
    }
    if (m == 0) {
      if (_rr_cut > 0 && IsOdd(im) && bmp > 0) {	
	c = fabs(bm)/bmp;
	if (c < 1) {
	  c = fabs(bm)/(1-c);
	  if (c < _rr_cut*fabs(b[0])) {
	    mx = im+1;
	    break;
	  }
	}
      }
    } else {
      for (i = mx0; i <= im; i++) {
	L = i/2 + 1;
	L2 = 2*L;
	if (i < im) {
	  ip0 = im;
	  ip1 = im;
	} else {
	  ip0 = mx0;
	  ip1 = im;
	}
	for (p = 0; p < nak[i]; p++) {
	  GetJLFromKappa(kak[i][p], &j1, &kl1);
	  k = OrbitalIndex(0, kak[i][p], e);
	  ph1 = GetPhaseShift(k)-(PI/4)*kl1;
	  for (ip = ip0; ip <= ip1; ip++) {
	    Lp = ip/2 + 1;
	    Lp2 = 2*Lp;
	    for (pp = 0; pp < nak[ip]; pp++) {
	      GetJLFromKappa(kak[ip][pp], &j2, &kl2);
	      k = OrbitalIndex(0, kak[ip][pp], e);
	      ph2 = GetPhaseShift(k)-(PI/4)*kl2;
	      c = (j1+1.0)*(j2+1.0)*(kl1+1.0)*(kl2+1.0)*(L2+1.0)*(Lp2+1.0);
	      c = sqrt(c)*ak[i][p]*akp[ip][pp];
	      if (ph1 != ph2) c *= cos(ph1-ph2);
	      if (IsOdd((L2-Lp2+j1-j2+j0+1)/2)) c = -c;
	      for (q = 0; q < m; q++) {
		q2 = 2*q;
		d = c*(q2 + 1.0);
		d *= W3j(kl1, kl2, q2, 0, 0, 0);
		d *= W6j(kl1, kl2, q2, j2, j1, 1);
		d *= W6j(j1, j2, q2, Lp2, L2, j0);
		w3j1 = W3j(q2, Lp2, L2, 0, 2, -2);
		w3j2 = 0.0;
		b[q+1] += d*w3j1;
		if (q >= 2) {
		  if (IsOdd((kl2-Lp2-kl0)/2)) d = -d;
		  d *= exp(0.5*(ln_factorial[q-2]-ln_factorial[q+2]));
		  d *= q*(q-1.0);
		  w3j2 = W3j(q2, Lp2, L2, 4, -2, -2);
		  b[q+1+m] += d*w3j2;
		}
	      }	  
	    }
	  }
	}
      }
      if (_rr_cut > 0 && IsOdd(im)) {
	mm = 2*(1+(im-1)/2)+1;
	bmc = 0.0;
	for (q = 0; q < mm; q += 2) {
	  bmc += b[q+1]*pqa[q];
	}
	bmi = bmc - bmt;
	bmt = bmc;
	bpc = 0.0;
	for (q = 2; q < mm; q += 2) {
	  bpc += pqa2[q-2]*b[q+m+1]/(q*(q-1.0));
	}
	bpi = bpc - bpt;
	bpt = bpc;	
	if (bmp > 0) {
	  c = fabs(bm)/bmp;
	  //printf("im: %d %g %g %g %g %g %g %g %g %g %g %g\n", im, c, bmp, bm, bm/(1-c), b[0], bmt, bmi, bmi/(1-c), bpt, bpi, bpi/(1-c));
	  if (c < 1) {
	    c = 1-c;
	    ph1 = fabs(bm)/c;
	    d = fabs(bpi)/c;
	    c = fabs(bmi)/c;
	    if (ph1 < _rr_cut*fabs(b[0]) &&
		c < _rr_cut*fabs(bmt) &&
		d <= _rr_cut*fabs(bpt)) {
	      mx = im+1;
	      break;
	    }
	  }
	}
      }
    }
  }
  
  for (i = 1; i < 2*m+1; i++) {
    b[i] /= b[0];
  }  
  
  double rb = b1/b[0];
  c = 2.0*PI/(j0+1.0);
  b[0] *= c;
  for (i = mx0; i < mx; i++) {
    free(ak[i]);
    free(kak[i]);
    if (k0p != k0) free(akp[i]);
  }
  free(nak);
  free(ak);
  if (k0p != k0) free(akp);
  free(kak);
  
  *mx1 = mx;
  return rb;
}

int SaveAsymmetry(char *fn, char *s, int mx, double te) {
  ORBITAL *orb0;
  CONFIG *cfg;
  char *p, sp[16], js;
  int k, ns, i, j, q, ncfg, m, mlam, mxi, mx0, mp;
  int kappa, n, jj, kl, k0, k1, mx1[MAXRRNUSR], mm, mt;
  double **b, **bi, e0, e, emin, emax, a, phi;
  double phi90, phi1, phi2, bphi, rp;
  double *pqa, *pqa2, nu1, theta;
  int *ipqa, ierr, nudiff, mu1;
  FILE *f;
  TFILE *fr;
  F_HEADER fh;
  RO_HEADER h;
  RO_RECORD r;
  int swp;

  fr = OpenFileRO(s, &fh, &swp);
  if (fr == NULL) {
    if (strlen(s)<2 || !isdigit(s[0])) return 0;
    ns = StrSplit(s, ' ');
    if (ns <= 0) return 0;
  }
  te /= HARTREE_EV;
  mxi = mx;
  mx0 = mx/1000;
  mp = mx0/1000;
  if (mp > 1) mp = -1;
  mx0 = mx0%1000;
  mx = mx%1000;
  mlam = 1 + (mx-1)/2;
  m = 2*mlam+1;
  pqa = malloc(sizeof(double)*m);
  pqa2 = malloc(sizeof(double)*m);
  ipqa = malloc(sizeof(int)*m);
  theta = acos(0.0);
  nu1 = 0;
  nudiff = mlam*2;
  mu1 = 0;
  DXLEGF(nu1, nudiff, mu1, mu1, theta, 3, pqa, ipqa, &ierr);
  nu1 = 2;
  nudiff = mlam*2-2;
  mu1 = 2;
  DXLEGF(nu1, nudiff, mu1, mu1, theta, 3, pqa2, ipqa, &ierr);
  int m2p = 2*m+1;
  p = s;
  if (fr == NULL) {
    f = fopen(fn, "a");
  } else {
    f = fopen(fn, "w");
  }
  if (n_usr <= 0) {
    if (n_egrid <= 0) {
      n_usr = 6;
    } else {
      n_usr = n_egrid;
    }
  }
  b = malloc(sizeof(double *)*n_usr);
  if (fr) {
    bi = malloc(sizeof(double *)*n_usr);
  }    
  for (i = 0; i < n_usr; i++) {
    b[i] = malloc(sizeof(double)*(m2p+1));
    if (fr) {
      bi[i] = malloc(sizeof(double)*(m2p+1));
    }
  }
  int eset = egrid[0] >= 0;
  int uset = usr_egrid[0] >= 0;
  if (fr == NULL) {
    for (k = 0; k < ns; k++) {
      while (*p == ' ') p++;
      ncfg = GetConfigFromString(&cfg, p);
      for (j = ncfg-1; j >= 0; j--) {
	if (cfg[j].n_shells != 1) continue;
	n = (cfg[j].shells)[0].n;
	kappa = (cfg[j].shells)[0].kappa;
	free(cfg[j].shells);
	GetJLFromKappa(kappa, &jj, &kl);
	if (jj > kl) js = '+';
	else js = '-';
	SpecSymbol(sp, kl/2);
	k0 = OrbitalIndex(n, kappa, 0);
	if (te <= 0) {
	  orb0 = GetOrbital(k0);
	  e0 = -(orb0->energy);
	} else {
	  e0 = te;
	}
	SetAWGrid(1, e0*FINE_STRUCTURE_CONST, e0);
	if (!uset) {
	  if (eset) {
	    SetUsrPEGridDetail(n_egrid, egrid);
	  } else {
	    if (egrid_limits_type == 0) {
	      emin = egrid_min*e0;
	      emax = egrid_max*e0;
	    } else {
	      emin = egrid_min;
	      emax = egrid_max;
	    }
	    SetUsrPEGrid(n_usr, emin, emax, e0);
	  }
	}
	if (k == 0) {
	  for (i = 0; i < n_usr; i++) {
	    xusr[i] = 2.0*usr_egrid[i];
	    xusr[i] *= (1.0+0.5*FINE_STRUCTURE_CONST2*usr_egrid[i]);
	    xusr[i] = sqrt(xusr[i])/FINE_STRUCTURE_CONST;
	    xusr[i] *= HARTREE_EV;
	  }
	}
	double et0 = e0;
	mt = 0;
	for (i = 0; i < n_usr; i++) {
	  e = usr_egrid[i];
	  mx1[i] = mx;
	  b[i][m2p] = AsymmetryPI(k0, k0, e, et0, mp, mx0, &mx1[i],
				  m, b[i], pqa, pqa2);
	  if (mx1[i] > mt) mt = mx1[i];
	}
	mm = 2*(1 + (mt-1)/2)+1;
	e0 *= HARTREE_EV;
	fprintf(f, "#  %2s  %2d %2d\n", 
		GetAtomicSymbol(), (int)GetAtomicNumber(), (int)GetResidualZ());
	fprintf(f, "#  %d%s%c %d %d %d %12.5E  %d %d %d\n",
		n, sp, js, n, kl, jj, e0, n_usr, mxi, mt);
	
	for (i = 0; i < n_usr; i++) {
	  e = usr_egrid[i]*HARTREE_EV;
	  a = (e+e0)/xusr[i];
	  a = a*a;
	  phi = b[i][0]*AREA_AU20;
	  phi90 = 0.0;
	  for (q = 0; q < mm; q += 2) {
	    phi90 += b[i][q+1]*pqa[q];
	  }
	  phi1 = phi90;
	  phi2 = phi90;
	  for (q = 2; q < mm; q += 2) {
	    bphi = pqa2[q-2]*b[i][q+m+1]/(q*(q-1.0));
	    phi1 -= bphi; /* parallel, Phi=0 */
	    phi2 += bphi; /* perpendicular Phi=90 */
	  }
	  phi90 *= phi/(4.0*PI);
	  if (phi1) rp = phi2/phi1;
	  else rp = phi2>0?1e30:-1e30;
	  fprintf(f, "%12.5E %12.5E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %d\n",
		  e, e+e0, phi, phi90, a*phi, a*phi90, rp, b[i][m2p], mx1[i]/2);
	}      
	for (q = 0; q < mm; q++) {
	  for (i = 0; i < n_usr; i++) {
	    e = usr_egrid[i]*HARTREE_EV;
	    fprintf(f, "%12.5E %15.8E %15.8E\n",
		    e, b[i][q+1], b[i][q+1+m]);
	  }
	}
      }
      if (ncfg > 0) free(cfg);
      while (*p) p++;
      p++;
    }
  } else {
    int nb, ir, nr, ip, iq, iq0, vn, vl, vnq, vlq, jq;
    LEVEL *lev0=NULL, *lev1=NULL;
    for (nb = 0; nb < fh.nblocks; nb++) {
      nr = ReadROHeader(fr, &h, swp);
      for (ir = 0; ir < h.ntransitions; ir++) {
	nr = ReadRORecord(fr, &r, swp);	       	    
	for (i = 0; i < n_usr; i++) {
	  for (q = 0; q <= m2p; q++) {
	    b[i][q] = 0.0;
	  }
	}
	lev0 = GetLevel(r.b);
	lev1 = GetLevel(r.f);
	te = lev1->energy - lev0->energy;
	e0 = te*HARTREE_EV;
	SetAWGrid(1, te*FINE_STRUCTURE_CONST, te);
	if (!uset) {
	  if (eset) {
	    SetUsrPEGridDetail(n_egrid, egrid);
	  } else {
	    if (egrid_limits_type == 0) {
	      emin = egrid_min*te;
	      emax = egrid_max*te;
	    } else {
	      emin = egrid_min;
	      emax = egrid_max;
	    }
	    SetUsrPEGrid(n_usr, emin, emax, te);
	  }
	}
	for (i = 0; i < n_usr; i++) {
	  xusr[i] = 2.0*usr_egrid[i];
	  xusr[i] *= (1.0+0.5*FINE_STRUCTURE_CONST2*usr_egrid[i]);
	  xusr[i] = sqrt(xusr[i])/FINE_STRUCTURE_CONST;
	  xusr[i] *= HARTREE_EV;
	}
	fprintf(f, "#  %2s  %2d %2d\n", 
		GetAtomicSymbol(), (int)GetAtomicNumber(), (int)GetResidualZ());
	for (ip = 0; ip < r.n; ip++) {
	  vn = abs(r.nk[ip]);
	  vl = 2*(vn%100);
	  vn = vn/100;
	  if (r.nk[ip] < 0) {
	    j = vl-1;
	  } else {
	    j = vl+1;
	  }
	  k0 = OrbitalIndex(vn, GetKappaFromJL(j, vl), 0);
	  if (_asym_ci > 0) iq0 = 0;
	  else iq0 = ip;
	  for (iq = iq0; iq <= ip; iq++) {
	    if (iq == ip) {
	      k1 = k0;
	    } else {
	      vnq = abs(r.nk[iq]);
	      vlq = 2*(vnq%100);
	      vnq = vnq/100;
	      if (r.nk[iq] < 0) {
		jq = vlq-1;
	      } else {
		jq = vlq+1;
	      }
	      if (jq != j) continue;
	      k1 = OrbitalIndex(vnq, GetKappaFromJL(jq, vlq), 0);
	    }	      
	    double nq = r.nq[ip]*r.nq[iq];
	    if (ip != iq) nq = 2*nq;
	    mt = 0;
	    for (i = 0; i < n_usr; i++) {
	      e = usr_egrid[i];
	      mx1[i] = mx;
	      bi[i][m2p] = AsymmetryPI(k0, k1, e, te, mp, mx0, &mx1[i],
				       m, bi[i], pqa, pqa2);
	      b[i][0] += bi[i][0]*nq;
	      if (mx1[i] > mt) mt = mx1[i];
	      for (q = 1; q <= m2p; q++) {
		b[i][q] += bi[i][q]*bi[i][0]*nq;
	      }
	    }
	    mm = 2*(1 + (mt-1)/2)+1;
	  }
	}
	fprintf(f, "#  %6d %6d %12.5E  %d %d %d\n",
		r.b, r.f, e0, n_usr, mxi, mt);
	int pb, pf, jb, jf;
	DecodePJ(lev0->pj, &pb, &jb);
	DecodePJ(lev1->pj, &pf, &jf);
	double wb = 1.0+jb;
	double wf = 1.0+jf;
	for (i = 0; i < n_usr; i++) {
	  for (q = 1; q <= m2p; q++) {
	    b[i][q] /= b[i][0];
	  }
	  e = usr_egrid[i]*HARTREE_EV;
	  a = (e+e0)/xusr[i];
	  a = a*a;
	  phi = b[i][0]*AREA_AU20;
	  phi90 = 0.0;
	  for (q = 0; q < mm; q += 2) {
	    phi90 += b[i][q+1]*pqa[q];
	  }
	  phi1 = phi90;
	  phi2 = phi90;
	  for (q = 2; q < mm; q += 2) {
	    bphi = pqa2[q-2]*b[i][q+m+1]/(q*(q-1.0));
	    phi1 -= bphi; /* parallel, Phi=0 */
	    phi2 += bphi; /* perpendicular Phi=90 */
	  }
	  phi90 *= phi/(4.0*PI);
	  if (phi1) rp = phi2/phi1;
	  else rp = phi2>0?1e30:-1e30;	
	  fprintf(f, "%12.5E %12.5E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %d\n",
		  e, e+e0, phi/wb, phi90/wb,
		  a*phi/wf, a*phi90/wf, rp, b[i][m2p], mx1[i]/2);
	}      
	for (q = 0; q < mm; q++) {
	  for (i = 0; i < n_usr; i++) {
	    e = usr_egrid[i]*HARTREE_EV;
	    fprintf(f, "%12.5E %15.8E %15.8E\n",
		    e, b[i][q+1], b[i][q+1+m]);
	  }
	}
      }      	    
    }
  }
  
  fclose(f);
  if (fr) {
    FCLOSE(fr);
  }
  for (i = 0; i < n_usr; i++) {
    free(b[i]);
    if (fr) {
      free(bi[i]);
    }
  }
  free(b);
  if (fr) {
    free(bi);
  }
  free(pqa);
  free(pqa2);
  free(ipqa);
  
  ReinitRadial(1);
  ReinitRecombination(1);

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

int FreeRecPk(void) {
  MultiFreeData(pk_array, FreeRecPkData);
  return 0;
}

int FreeRecQk(void) {
  if (qk_array->array == NULL) return 0;
  MultiFreeData(qk_array, FreeRecPkData);
  return 0;
}

int InitRecombination(void) {
  int blocks[5];
  int ndim;
  int i;
  
  ndim = 5;
  for (i = 0; i < ndim; i++) blocks[i] = MULTI_BLOCK5;
  pk_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(pk_array, sizeof(double *), ndim, blocks, "rpk_array");
  
  ndim = 2;
  for (i = 0; i < ndim; i++) blocks[i] = MULTI_BLOCK3;
  qk_array = (MULTI *) malloc(sizeof(MULTI));
  blocks[0] = 16;
  blocks[1] = 16;
  MultiInit(qk_array, sizeof(double *), ndim, blocks, "rqk_array");  
  
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

  SetMaxAICache(-1);
  SetCXEGrid(36, 1e-2, 1e5, NULL, 1, 0);
  
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

  if (m > 0) return 0;

  SetRecQkMode(QK_DEFAULT, 0.1);
  SetRecPWOptions(RECLMAX, RECLMAX);

  for (i = 0; i < n_complex; i++) {
    ArrayFree(rec_complex[i].rg, NULL);
  }
  n_complex = 0;
  
  return 0;
}

void SetOptionRecombination(char *s, char *sp, int ip, double dp) {
  if (0 == strcmp(s, "recombination:cxldist")) {
    cxldist = ip;
    return;
  }
  if (0 == strcmp(s, "recombination:cxldistmj")) {
    cxldistmj = dp;
    return;
  }
  if (0 == strcmp(s, "recombination:cxrotcouple")) {
    cxrotcouple = dp;
    return;
  }
  if (0 == strcmp(s, "recombination:progress_report")) {
    _progress_report = ip;
    return;
  }
  if (0 == strcmp(s, "recombination:ai_cut")) {
    _ai_cut = dp;
    return;
  }
  if (0 == strcmp(s, "recombination:rr_cut")) {
    _rr_cut = dp;
    return;
  }
  if (0 == strcmp(s, "recombination:maxkl")) {
    pw_scratch.max_kl = ip;
    return;
  }
  if (0 == strcmp(s, "recombination:interpkl")) {
    pw_scratch.kl_interp = ip;
    return;
  }
  if (0 == strcmp(s, "recombination:asym_ci")) {
    _asym_ci = ip;
    return;
  }
}

double LandauZenerDE(CXTGT *cx, double z, double de, double r) {
  double vp, vr, vc, r2, rp, r0;
  if (cx->x > 0) r0 = cx->x/cx->b;
  else {
    r0 = pow(0.5*cx->a*z/(-cx->x*25.0), 0.25);
  }
  rp = r + r0;
  r2 = rp*rp;
  vp = -cx->a*z*z*0.5/(r2*r2);
  vr = 25.0*z*exp(-cx->b*r);
  vc = (z-1.0)/r;
  return vp+vr-vc+de;
}

double LandauZenerDEDR(CXTGT *cx, double z, double r, double *lam) {
  double r2 = r*r;
  double r0;
  if (cx->x > 0) r0 = cx->x/cx->b;
  else {
    r0 = pow(0.5*cx->a*z/(-cx->x*25.0), 0.25);
  }
  double rp = r + r0;
  double rp2 = rp*rp;
  double vdr = (z-1.0)/r2;
  double v1 = cx->a*z*z/(rp2*rp2);
  double v2 = 25.0*z*exp(-cx->b*r);
  double mu = 2.0*v1/rp - cx->b*v2;
  *lam = 0.5*v1 - v2;
  if (mu < 0) mu = 0.0;  
  if (*lam < 0) *lam = 0.0;
  if (r0 > 0) {
    mu *= 1.0-exp(-r/r0);
  }
  return vdr+mu;
}

double LandauZenerRX(CXTGT *cx, double z, double de, double r) {
  double dv, r0, r1;
  int i;
  dv = LandauZenerDE(cx, z, de, r);
  r0 = r;
  r1 = r;
  if (dv > 0) {
    i = 0;
    r0 = r;
    while (dv > 0) {
      r1 = r0;
      i++;
      if (i > 1000) {
	printf("LZ RX fail0: %g %g %g\n", r0, de, dv);
	return 0.0;
      }
      r0 /= 1.05;
      dv = LandauZenerDE(cx, z, de, r0);
    }
  } else if (dv < 0) {
    i = 0;
    r1 = r;
    while (dv < 0) {
      r0 = r1;
      i++;
      if (i > 1000) {
	printf("LZ RX fail1: %g %g %g\n", r1, de, dv);
	return 0.0;
      }
      r1 *= 1.05;
      dv = LandauZenerDE(cx, z, de, r1);
    }
  } else {
    return r;
  }
  i = 0;
  while (fabs(r1-r0)/(r0+r1) > 1e-4) {
    i++;
    if (i > 1000) {
      printf("LZ RX fail2: %g %g\n", r, dv);
      return 0.0;
    }
    r = 0.5*(r0+r1);
    dv = LandauZenerDE(cx, z, de, r);
    if (dv < 0) {
      r0 = r;
    } else if (dv > 0) {
      r1 = r;
    } else {
      break;
    }
  }
  r = 0.5*(r0+r1);
  return r;
}

double LandauZenerUX(CXTGT *cx, double z, double rx, int m) {
  double ux;
  double sz = sqrt(z);
  ux = exp(-1.324*rx*sqrt(2.0*cx->e)/sz)*(9.13/sz);
  if (m > 0) {
    double u2 = 0.5*exp(-0.4*rx);
    int i;  
    for (i = 0; i < m; i++) {
      ux *= u2;
    }
  }
  return ux;
}

double RotAlpha(double nzv, double b, double rx) {
  double x = nzv/b;
  return acos(b/rx)*sqrt(1 + x*x);
}

double LandauZenerRC(double rx, double nzv, double beta, int n) {
  if (n == 1) return 0.0;  
  double b0, b1, b;
  int np = 2*(n-1);
#define nbs 25
  double bsq[nbs], asq[nbs], qn[nbs];  
  const int nr = 6;
  double p, bx, x, x2, sa, sb, qni;
  int i, ir;
  double asp, a, db, bsp, bsp1, sigma;

  b1 = rx;
  bsp1 = 1.0;
  sigma = 0.0;
  for (ir = 0; ir < nr; ir++) {    
    b0 = nzv*0.1;
    if (b0 > b1) b0 = b1;
    asp = PI*(ir+1);
    a = RotAlpha(nzv, b0, rx);
    int i0 = 0;
    while (a < asp) {
      i0++;
      b1 = b0;
      b0 *= 0.5;
      a = RotAlpha(nzv, b0, rx);
      if (i0 > 1000) {
	printf("LZRC fail0: %d %g %g\n", ir, a, asp);
      }
    }
    int i1 = 0;
    while (fabs(b1-b0)/fabs(b0+b1) > 1e-8) {
      i1++;
      if (i1 > 1000) {
	printf("LZRC fail1: %g %g\n", b0, b1);
      }
      b = 0.5*(b0+b1);
      a = RotAlpha(nzv, b, rx);
      if (a > asp) {
	b0 = b;
      } else if (a < asp) {
	b1  = b;
      } else {
	break;
      }
    }
    bsp = 0.5*(b0+b1);
    b1 = bsp;
    bsp = bsp/rx;
    bsp = bsp*bsp;
    //printf("#bsp: %d %d %d %d %g %g %g\n", ir, n, i0, i1, bsp, a, asp);
    db = (bsp1-bsp)/(nbs-1);
    bsq[0] = bsp;
    asq[0] = PI*(ir+1);
    qn[0] = 0.0;
    bsq[nbs-1] = bsp1;
    asq[nbs-1] = PI*ir;
    qn[nbs-1] = 0.0;
    for (i = 1; i < nbs-1; i++) {
      bsq[i] = bsq[i-1] + db;
      b = sqrt(bsq[i])*rx;
      asq[i] = RotAlpha(nzv, b, rx);
      sa = sin(asq[i]);
      sa *= sa;
      sb = sin(atan(b/nzv));
      sb *= sb;
      qni = 1-pow(1-sa*sb, np);
      x2 = 1.0/(1-bsq[i]);
      x = sqrt(x2);
      bx = beta*x;
      if (bx < 1e-3) {
	p = bx;
      } else {
	p = 1-exp(-bx);
      }
      qn[i] = qni*p*p;
      //printf("%d %d %g %g %g %g %g %g %g %g\n", ir, i, bsq[i], b, asq[i], qni, bx, p, qn[i], PI*rx*rx*sigma);
    }
    sigma += Simpson(qn, 0, nbs-1)*db;
    bsp1 = bsp;
  }
  bsq[nbs-1] = bsq[0];
  asq[nbs-1] = asq[0];
  qn[nbs-1] = qn[0];
  bsq[0] = 0;
  asq[0] = 0;
  qn[0] = 0.0;
  db = bsp1/(nbs-1);  
  for (i = 1; i < nbs-1; i++) {
    bsq[i] = bsq[i-1] + db;
    b = sqrt(bsq[i])*rx;
    sa = 0.5;
    sb = sin(atan(b/nzv));
    sb *= sb;
    qni = 1-pow(1-sa*sb, np);
    x2 = 1.0/(1-bsq[i]);
    x = sqrt(x2);
    bx = beta*x;
    if (bx < 1e-3) {
      p = bx;
    } else {
      p = 1-exp(-bx);
    }
    qn[i] = qni*p*p;
    //printf("%d %d %g %g %g %g %g %g %g %g\n", ir, i, bsq[i], b, asq[i], qni, bx, p, qn[i], PI*rx*rx*sigma);
  }
  sigma += Simpson(qn, 0, nbs-1)*db;
  //printf("rotcoup: %d %g %g %g %g\n", n, rx, beta, nzv, PI*rx*rx*sigma);
#undef nbc
  return PI*rx*rx*sigma;
}

double **LandauZenerMC(CXTGT *cx, double z, int ne0, double *e0, int *nmc) {
  int n, i, j, n1, n2;
  double lam[1000], rx[1000], de[1000], beta[1000];
  double v12[1000], vdr[1000];
  double **pc, r;
  for (n = 1; n < 1000; n++) {
    lam[n] = -1;
    rx[n] = -1;
    de[n] = -1;
    v12[n] = -1;
    vdr[n] = -1;
    beta[n] = 0;
  }
  n2 = (int)(1+pow(z, 0.768));
  n1 = 1000;
  for (i = 0; i < ne0; i++) {
    for (n = 1; n < n1; n++) {
      r = LandauZenerCX(cx, n, -1, -1, z, e0[i], -1,
			rx+n, v12+n, vdr+n, lam+n, de+n, beta+n);
      if (i == 0 && n > n2 && rx[n] <= 0) {
	n1 = n-1;
	pc = malloc(sizeof(double *)*n1);
	for (j = 0; j < n1; j++) {
	  pc[j] = malloc(sizeof(double)*ne0);
	}
	break;
      }
    }
    pc[n1-1][i] = 0.0;
    for (n = n1-1; n >= 1; n--) {
      pc[n-1][i] = pc[n][i] - beta[n+1]*2.0;
    }
  }
  *nmc = n1;
  return pc;
}
		      
double LandauZenerCX(CXTGT *cx, int n, int k, int m0,
		     double z, double e0, double ei,
		     double *orx, double *ov12, double *ovdr,
		     double *olam, double *ode, double *obeta) {
  double de, beta, sigma, f3, lam, elam, rx, v, v12, vdr, wnk;
  double *mass = GetAtomicMassTable();
  int m;

  if (m0 < 0) m = 0;
  else m = m0;
  if (orx == NULL || *orx < 0) {
    if (ei < 0) {
      ei = -EnergyH(z, n, -1);
    }
    de = ei - cx->e;
    if (de < 0.0) {
      if (orx) *orx = 0;
      return 0.0;
    }
    rx = LandauZenerRX(cx, z, de, 3*n*n/z);
    if (rx <= 0) {
      if (orx) *orx = 0;
      return 0.0;
    }
    v12 = LandauZenerUX(cx, z, rx, m);
    vdr = LandauZenerDEDR(cx, z, rx, &lam);
    if (orx) {
      *ode = de;
      *ov12 = v12;
      *ovdr = vdr;
      *orx = rx;
      *olam = lam;
    }
  } else {
    de = *ode;
    rx = *orx;
    v12 = *ov12;
    vdr = *ovdr;
    lam = *olam;
  }
  wnk = 1.0;
  v = sqrt(2*e0/AMU);
  int iz = (int)(z-1);
  double rmass = mass[iz]*cx->m/(mass[iz]+cx->m);
  elam = lam/(e0*rmass);
  double slam = sqrt(1+elam);
  beta = TWO_PI*v12*v12/(v*vdr*slam);
  if (beta <= 0) {
    beta = 0.0;
  }
  *obeta = beta;
  if (m0 < 0) return 0.0;
  if (beta <= 0) return 0.0;
  double eb = exp(-beta);
  if (beta < 1e-3) {
    f3 = beta;
  } else {
    f3 = eb*(EXPINT(beta, 3) - eb*EXPINT(2*beta, 3));  
    if (f3 < 0) f3 = 0.0;
  }
  sigma = FOUR_PI*rx*rx*(1+elam)*f3;
  if (cxrotcouple > 0) {
    double sigrc = LandauZenerRC(rx*slam, 1.5*n/(z*v), beta, n);
    sigma += cxrotcouple*sigrc;
  }
  if (k >= 0 && k < n) {
    wnk = LandauZenerLD(ln_factorial, n, k, z, v*rx, cxldist, cxldistmj);
    sigma *= wnk;
  }
  return sigma;
}

void LandauZenerBareCX(char *fn, double z, int n0, double *e0, int ldm) {
  CXTGT *cx = GetCXTarget();
  FILE *f;

  f = NULL;
  if (fn != NULL) {
    f = fopen(fn, "w");
    if (f == NULL) {
      printf("cannot open file: %s\n", fn);
    }
  }
  if (f == NULL) f = stdout;
  double lam[1000], rx[1000], de[1000], beta[1000];
  double v12[1000], vdr[1000], cs[1000];
  char ss[8], su[8], snk[16];
  double rt, r, w, v;
  int n1, n2, i, j, k;
  for (i = 0; i < 1000; i++) {
    lam[i] = -1;
    rx[i] = -1;
    de[i] = -1;
    v12[i] = -1;
    vdr[i] = -1;
    beta[i] = 0;
    cs[i] = 0;
  }
  n1 = 1;
  n2 = (int)(1+z/sqrt(2*cx->e));
  r = e0[0]/HARTREE_EV;
  for (j = n2; j >= n1; j--) {
    cs[j] = LandauZenerCX(cx, j, -1, 0, z, r, -1.0,
			  rx+j, v12+j, vdr+j, lam+j, de+j, beta+j);
    if (rx[j] > 0 && beta[j] > 0 && cs[j] > 1e-10) break;
  }
  n2 = j;
  fprintf(f, "# z=%g\n", z);
  fprintf(f, "# tgts=%s\n", cx->symbol);
  fprintf(f, "# tgtz=%g\n", cx->z);
  fprintf(f, "# tgtm=%g\n", cx->m);
  fprintf(f, "# tgta=%g\n", cx->a);
  fprintf(f, "# tgtb=%g\n", cx->b);
  fprintf(f, "# tgte=%g\n", cx->e);
  fprintf(f, "# tgtx=%g\n", cx->x);
  fprintf(f, "# ldist=%d\n", ldm);
  fprintf(f, "# ldistmj=%g\n", cxldistmj);
  fprintf(f, "%-10s", "# (eV/u)");
  for (j = n2; j >= n1; j--) {
    if (ldm < 0) {
      fprintf(f, " n=%-6d", j);
    } else {
      for (k = j-1; k >=0; k--) {
	SpecSymbol(ss, k);
	StrUpper(su, ss, 8);
	sprintf(snk, "%d%s(2%s)", j, ss, su);
	fprintf(f, " %-8s", snk);
      }
    }
  }
  fprintf(f, " Total\n");
  for (i = 0; i < n0; i++) {
    e0[i] /= HARTREE_EV;
    v = sqrt(2*e0[i]/AMU);
    rt = 0.0;
    for (j = n2; j >= n1; j--) {
      cs[j] = LandauZenerCX(cx, j, -1, 0, z, e0[i], -1.0,
			    rx+j, v12+j, vdr+j, lam+j, de+j, beta+j);
    }
    if (i == 0) {
      fprintf(f, "%-10s", "# RX");
      for (j = n2; j >= n1; j--) {
	if (ldm < 0) {
	  fprintf(f, " %8.2E", rx[j]);
	} else {
	  for (k = j-1; k >= 0; k--) {
	    fprintf(f, " %8.2E", rx[j]);
	  }
	}
      }
      fprintf(f, "\n");
      fprintf(f, "%-10s", "# V12");
      for (j = n2; j >= n1; j--) {
	if (ldm < 0) {
	  fprintf(f, " %8.2E", v12[j]);
	} else {
	  for (k = j-1; k >= 0; k--) {
	    fprintf(f, " %8.2E", v12[j]);
	  }
	}
      }
      fprintf(f, "\n");
      fprintf(f, "%-10s", "# VDR");
      for (j = n2; j >= n1; j--) {
	if (ldm < 0) {
	  fprintf(f, " %8.2E", vdr[j]);
	} else {
	  for (k = j-1; k >= 0; k--) {
	    fprintf(f, " %8.2E", vdr[j]);
	  }
	}
      }
      fprintf(f, "\n");
      fprintf(f, "%-10s", "# LAM");
      for (j = n2; j >= n1; j--) {
	if (ldm < 0) {
	  fprintf(f, " %8.2E", lam[j]);
	} else {
	  for (k = j-1; k >= 0; k--) {
	    fprintf(f, " %8.2E", lam[j]);
	  }
	}
      }
      fprintf(f, "\n");
    }
    fprintf(f, "%10.4E", e0[i]*HARTREE_EV);
    double pc = 0.0;
    for (j = n2; j >= n1; j--) {
      if (cs[j] > 0) {
	cs[j] = log(cs[j]) + pc;
	if (cs[j] <= -225) cs[j] = 0;
	else cs[j] = exp(cs[j]);
      }
      pc -= beta[j]*2.0;
    }
    for (j = n2; j >= n1; j--) {
      r = cs[j]*AREA_AU20*1e-4;
      if (ldm < 0) {
	fprintf(f, " %8.2E", r);
      } else {
	for (k = j-1; k >= 0; k--) {
	  w = LandauZenerLD(ln_factorial, j, k, z, v*rx[j], ldm, cxldistmj);
	  fprintf(f, " %8.2E", r*w);
	}
      }
      rt += r;
    }
    fprintf(f, " %8.2E\n", rt);
  }
  fflush(f);
  if (f != stdout) fclose(f);
}
