#include "ionization.h"

static char *rcsid="$Id: ionization.c,v 1.18 2001/10/12 03:15:00 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

#define NINT 15
#define MAXNQK (MAXNE*(MAXNE+1)/2)
#define NPARAMS 4

static int output_format = 0;
static int egrid_type = -1;
static int usr_egrid_type = -1;
static int pw_type = -1;
static int qk_mode = QKDETAIL;

static double yegrid0[NINT];
static double log_yegrid0[NINT];

static int n_usr = 0;
static double usr_egrid[MAXNUSR];
static double log_usr[MAXNUSR];
static double xusr[MAXNUSR];
static double log_xusr[MAXNUSR];
static double qk_usr[MAXNUSR];

static int n_egrid = 0;
static double egrid[MAXNE];
static double log_egrid[MAXNE];
static double sigma[MAXNE];
static double xegrid[MAXNTE][MAXNE];
static double log_xegrid[MAXNTE][MAXNE];

static int n_tegrid = 0;
static double tegrid[MAXNTE];
static double log_te[MAXNTE];

static double _dwork[17*MAXNQK];
static int _iwork[39*MAXNQK];

void uvip3p_(int *np, int *ndp, double *x, double *y, 
	     int *n, double *xi, double *yi);
void sdbi3p_(int *md, int *ndp, double *xd, double *yd, 
	     double *zd, int *nip, double *xi, double *yi,
	     double *zi, int *ierr, double *wk, int *iwk);

static struct {
  int max_k;
  int qr;
  int max_kl;
  int max_kl_eject;
  int kl_cb;
  double tolerence;
  int nkl0;
  int nkl;
  int ns;
  double kl[MAXNKL+1];
  double log_kl[MAXNKL];
} pw_scratch = {12, 0, MAXKL, 8, 0, 1E-2, 0, 0};

static MULTI *qk_array;

int SetIEGrid(int n, double emin, double emax) {
  n_tegrid = SetTEGrid(tegrid, log_te, n, emin, emax);
  return n_tegrid;
}

int SetIEGridDetail(int n, double *x) {
  n_tegrid = SetTEGridDetail(tegrid, log_te, n, x);
  return n_tegrid;
}

int SetCIPWOptions(int qr, int max, int max_eject, int kl_cb, double tol) {
  pw_scratch.max_k = GetMaxRank();
  pw_scratch.qr = qr;
  if (max > MAXKL) {
    printf("The maximum partial wave reached in Ionization: %d\n", MAXKL);
    abort();
  }
  pw_scratch.max_kl = max;
  pw_scratch.max_kl_eject = max_eject;
  pw_scratch.kl_cb = kl_cb;
  pw_scratch.tolerence = tol;
  pw_scratch.nkl0 = 1;
  pw_scratch.kl[0] = 0;
  pw_scratch.log_kl[0] = -100.0;
  pw_scratch.nkl = 0;
  return 0;
}  

int SetCIEGridDetail(int n, double *xg) {
  n_egrid = SetEGridDetail(egrid, log_egrid, n, xg);
  return n_egrid;
}

int SetCIFormat(int m) {
  output_format = m;
  return m;
}

int SetUsrCIEGridType(int type) {
  if (type >= 0) usr_egrid_type = type;
  return 0;
}

int SetCIEGrid(int n, double emin, double emax, double eth) {
  n_egrid = SetEGrid(egrid, log_egrid, n, emin, emax, eth);
  return n_egrid;
}

int SetUsrCIEGridDetail(int n, double *xg) {
  if (n > MAXNUSR) {
    printf("Max # of grid points reached \n");
    return -1;
  }
  n_usr = SetEGridDetail(usr_egrid, log_usr, n, xg);
  return n_usr;
}

int SetUsrCIEGrid(int n, double emin, double emax, double eth) {
  if (n > MAXNUSR) {
    printf("Max # of grid points reached \n");
    return -1;
  }
  n_usr = SetEGrid(usr_egrid, log_usr, n, emin, emax, eth);
  return n_usr;
}

int SetCIPWGrid(int ns, int *n, int *step) {
  if (pw_scratch.nkl0 <= 0) SetCIPWOptions(0, 500, 8, 50, 5E-2);
  pw_scratch.nkl = SetPWGrid(&(pw_scratch.nkl0),
			     pw_scratch.kl,
			     pw_scratch.log_kl,
			     pw_scratch.max_kl,
			     &ns, n, step);
  pw_scratch.ns = ns;
  return 0;
}

int CIRadialQk(double *qk, int ie1, int ie2, int kb, int kbp, int k) {
  double pk[MAXNTE][MAXNKL];
  double y2[MAXNKL];
  double e1, e2, e0, te;
  ORBITAL *orb;
  int kappab, jb, klb, i, j, t, kl, klp;
  int js[4], ks[4];
  int jmin, jmax, ko2;
  int kf, kappaf, kf0, kappa0, kf1, kappa1;
  int kl0, kl0p, kl1, kl1p;
  int j0, j1, j1min, j1max;
  double z, z2, r, rp, sd, se, s;
  double qkjt, qkj, qklt, qkl, qkl0;
  double eps, a, b, h;
  int kl_max0, kl_max1, kl_max2, max0;
  int type, last_kl0, second_last_kl0;

  ko2 = k/2;

  for (i = 0; i < n_tegrid; i++) {
    for (j = 0; j < pw_scratch.nkl0; j++) {
      pk[i][j] = 0.0;
    }
  }
  for (i = 0; i < n_tegrid; i++) qk[i] = 0.0;
  
  e1 = egrid[ie1];
  r = GetRMax();
  z = GetResidualZ();
  z2 = z*z;
  t = r*sqrt(e1+2.0*z/r);
  kl_max0 = pw_scratch.max_kl;
  kl_max0 = Min(kl_max0, t);
  max0 = pw_scratch.ns;

  e2 = egrid[ie2];
  t = r*sqrt(e2+2.0*z/r);
  kl_max2 = pw_scratch.max_kl_eject;
  kl_max2 = Min(kl_max2, t);

  orb = GetOrbital(kb);
  kappab = orb->kappa;
  GetJLFromKappa(kappab, &jb, &klb);
  klb /= 2;
  jmin = abs(k - jb);
  jmax = k + jb;

  qkjt = 0.0;
  for (j = jmin; j <= jmax; j += 2) {
    if (qkjt && fabs(qkj/qkjt) < eps) break;
    qkj = 0.0;
    for (klp = j - 1; klp <= j + 1; klp += 2) {
      kappaf = GetKappaFromJL(j, klp);
      kl = klp/2;
      if (kl > kl_max2) break;
      if (IsEven(kl + klb + ko2)) {
	type = ko2;
      } else {
        type = -1;
      }
      if (kl < pw_scratch.qr) {
	js[2] = 0;
      } else {
	js[2] = j;
	if (IsOdd(kl)) {
	  if (kappaf < 0) kappaf = -kappaf - 1;
	} else {
	  if (kappaf > 0) kappaf = -kappaf - 1;
	}
      }
      kf = OrbitalIndex(0, kappaf, e2);	
      ks[2] = kf;  

      if (type >= 0 && type < CBMULTIPOLES) {
	kl_max1 = pw_scratch.kl_cb;
      } else {
	kl_max1 = kl_max0;
      }
      eps = pw_scratch.tolerence;
      if (type >= CBMULTIPOLES) {
	z = GetCoulombBetheAsymptotic(tegrid[0]+e2, e1);
      }

      last_kl0 = 0;
      second_last_kl0 = 0; 
      qklt = 0.0;
      for (t = 0; ; t++) {
        if (second_last_kl0) {
          last_kl0 = 1;
        } else {
	  kl0 = pw_scratch.kl[t];
          if (kl0 > max0) {
	    if (type < 0) {
	      rp *= qkl;
	      rp = rp/qklt;
	      if (rp < eps) last_kl0 = 1;
	    } else if (type >= CBMULTIPOLES) {
	      h = z*qkl;
	      h = h/(h+qklt);
	      if (h < eps) last_kl0 = 1;
	      else {
		rp = fabs(1.0 - rp/z);
		rp *= h;
		if (rp < eps) last_kl0 = 1;
	      }
	    } else {
	      z = (GetCoulombBethe(ie2, 0, ie1, type, 1))[t-1];
	      h = z*qkl;
	      h = h/(h+qklt);
	      if (h < eps) last_kl0 = 1;
	      else {
		z = (GetCoulombBethe(ie2, 0, ie1, type, 0))[t-1];
		z = z/(1.0-z);
		rp = fabs(1.0 - rp/z);
		rp *= h;
		if (rp < eps) last_kl0 = 1;
	      }
	    }
	    if (pw_scratch.kl[t+1] > kl_max1) {
	      second_last_kl0 = 1;
	    }  
	  } 
	}	  
	kl0p = 2*kl0;
	qkl0 = qkl;
        qkl = 0.0;
	for (j0 = abs(kl0p - 1); j0 <= kl0p + 1; j0 += 2) {
	  kappa0 = GetKappaFromJL(j0, kl0p);
	  if (kl0 < pw_scratch.qr) {
	    js[1] = 0;
	  } else {
	    js[1] = j0;
	    if (IsOdd(kl0)) {
	      if (kappa0 < 0) kappa0 = -kappa0 - 1;
	    } else {
	      if (kappa0 > 0) kappa0 = -kappa0 - 1;
	    }
	  }
	  j1min = abs(j0 - k);
	  j1max = j0 + k;
	  for (j1 = j1min; j1 <= j1max; j1 += 2) {
	    for (kl1p = j1 - 1; kl1p <= j1 + 1; kl1p += 2) {
	      kappa1 = GetKappaFromJL(j1, kl1p);
	      kl1 = kl1p/2;
	      if (kl1 < pw_scratch.qr) {
		js[3] = 0;
	      } else {
		js[3] = j1;
		if (IsOdd(kl1)) {
		  if (kappa1 < 0) kappa1 = -kappa1 - 1;
		} else {
		  if (kappa1 > 0) kappa1 = -kappa1 - 1;
		}
	      }
	      kf1 = OrbitalIndex(0, kappa1, e1);
	      ks[3] = kf1;
	      for (i = 0; i < n_tegrid; i++) {
		te = tegrid[i];		
		e0 = e1 + e2 + te;
		kf0 = OrbitalIndex(0, kappa0, e0);
		ks[1] = kf0;

		js[0] = 0;
		ks[0] = kb;
		if (kl1 >= pw_scratch.qr &&
		    kl0 >= pw_scratch.qr) {
		  SlaterTotal(&sd, &se, js, ks, k, -1);
		} else {
		  SlaterTotal(&sd, &se, js, ks, k, 1);
		} 
		r = sd + se;
                if (!r) break;

		if (kbp != kb) {
		  js[0] = 0;
		  ks[0] = kbp;
		  if (kl1 >= pw_scratch.qr &&
		      kl0 >= pw_scratch.qr) {
		    SlaterTotal(&sd, &se, js, ks, k, -1);
		  } else {
		    SlaterTotal(&sd, &se, js, ks, k, 1);
		  } 
		  rp = sd + se;
                  if (!rp) break;
		} else {
		  rp = r;
		}
		r = r*rp;
		if (i == 0) qkl += r;
		pk[i][t] += r;
	      }
	    }
	  }
	}
        if (t == 0) {
	  qklt += qkl;
	} else if (!last_kl0 && !second_last_kl0) {
	  if (qkl + 1.0 == 1.0) rp = 0.0;
	  else rp = qkl / qkl0;
	  h = kl0 - pw_scratch.kl[t-1];
	  if (h > 1) {
	    a = (1.0 - rp) * qkl0;
	    rp = pow(rp, 1.0/h);
	    rp = rp/(1.0 - rp);
	    a *= rp;
	    qklt += a;
	  } else {
	    qklt += qkl;
	    rp = rp/(1.0 - rp);
	  }
	} else if (last_kl0) {
	  if (type >= CBMULTIPOLES) {
	    for (i = 0; i < n_tegrid; i++) {
	      b = GetCoulombBetheAsymptotic(tegrid[i]+e2, e1);
	      r = qkl * b;
	      qk[i] += r;  
	    }
	  } else if (type >= 0) {
	    for (i = 0; i < n_tegrid; i++) {
	      b = (GetCoulombBethe(ie2, i, ie1, type, 1))[t];
	      r = qkl*b;
	      qk[i] += r;
	    }    
          } 
	  break;
	}
      }
      qkj += qklt;
    }  
    qkjt += qkj;
  }
  
  t = pw_scratch.nkl;
  for (i = 0; i < n_tegrid; i++) {
    spline(pw_scratch.log_kl, pk[i], t, 1E30, 1E30, y2);
    r = pk[i][0];
    for (j = 1; j < t; j++) {
      r += pk[i][j];
      kl0 = pw_scratch.kl[j-1];
      kl1 = pw_scratch.kl[j];
      for (kl = kl0+1; kl < kl1; kl++) {
	splint(pw_scratch.log_kl, pk[i], y2, t, LnInteger(kl), &s);
	r += s;
      }
    } 
    qk[i] += r;
  }

  return 0;
}

int CIRadialQkIntegrated(double *qke, double te, int kb, int kbp) {
  int np, i, j, k, nd;
  double *qk, qkc[MAXNTE];

  qk = CIRadialQkIntegratedTable(kb, kbp);
  if (qk == NULL) {
    return -1;
  }

  if (n_tegrid > 1) {
    nd = 1;
    np = 3;
    for (i = 0; i < NPARAMS; i++) {
      k = i;
      for (j = 0; j < n_tegrid; j++) {
	qkc[j] = qk[k];
	k += NPARAMS;
      }
      uvip3p_(&np, &n_tegrid, tegrid, qkc, &nd, &te, &qke[i]);
    }
  } else {
    for (i = 0; i < NPARAMS; i++) {
      qke[i] = qk[i];
    }
  }

  return 0;
}

void CIRadialQkIntegratedBasis(int npar, double *yb, double x, double log_x) {
  double y1, y2;

  y1 = 1.0/x;
  y2 = 1.0 - y1;
  yb[0] = log_x;
  yb[1] = y2*y2;
  yb[2] = log_x*y1 ;
  yb[3] = y2;
}
 
double *CIRadialQkIntegratedTable(int kb, int kbp) {
  int index[3], ie1, ie2, ite, i, j;
  double **p, *qkc, *qk, qkt[MAXNE];
  double xr[MAXNQK], yr[MAXNQK], zr[MAXNE];
  double x, y, integrand[NINT];
  int np, nd, ierr, nqk;

  index[0] = 1;
  index[1] = kb;
  index[2] = kbp;
  
  p = (double **) MultiSet(qk_array, index, NULL);
  if (*p) {
    return (*p);
  } else if (qk_mode == QKFIT) {
    return NULL;
  }    

  qk = CIRadialQkTable(kb, kbp);
  if (qk == NULL) return NULL;

  *p = (double *) malloc(sizeof(double)*n_tegrid*NPARAMS);
  qkc = *p;

  nd = 1;
  for (ite = 0; ite < n_tegrid; ite++) {
    for (ie1 = 0; ie1 < n_egrid; ie1++) {
      zr[ie1] = log(1.0 + egrid[ie1]/tegrid[ite]);
    }
    j = 0;
    for (ie1 = 0; ie1 < n_egrid; ie1++) {
      for (ie2 = 0; ie2 <= ie1; ie2++) {
	xr[j] = zr[ie1];
	yr[j] = zr[ie2];
	j++;
      }
    }
    nqk = j;
    np = 0;
    for (ie1 = 0; ie1 < n_egrid; ie1++) {
      for (i = 0; i < NINT; i++) {	
	x = log(1.0 + egrid[ie1]*(1.0-yegrid0[i])/tegrid[ite]);
	y = log(1.0 + egrid[ie1]*yegrid0[i]/tegrid[ite]);
	if (np == 0) np = 1;
	else np = 3;
	sdbi3p_(&np, &nqk, xr, yr, qk, &nd, &x, &y, &integrand[i],
		&ierr, _dwork, _iwork);
      }   
      x = Simpson(integrand, 0, NINT-1);
      x *= yegrid0[1]*egrid[ie1];
      qkt[ie1] = x;
    }
    
    SVDFit(NPARAMS, qkc, NULL, 1E-3, n_egrid, xegrid[ite], log_xegrid[ite],
	   qkt, sigma, CIRadialQkIntegratedBasis);
    qk += nqk;
    qkc += NPARAMS;
  }

  /*
  qkc = *p;
  for (j = 0; j < NPARAMS; j++) {
    i = j;
    for (ite = 0; ite < n_tegrid; ite++) {
      printf("%d %d %10.3E %10.3E\n", 
	     kb, kbp, tegrid[ite], qkc[i]);
      i += NPARAMS;
    }
    printf("\n\n");
  }
  printf("===============\n");
  */

  return (*p);
}

double *CIRadialQkTable(int kb, int kbp) {
  int index[3], ie, i, j, k, t;
  int nqk, ie1, ie2, ite;
  double **p, *qk, qi[MAXNTE];

  index[0] = 0;
  index[1] = kb;
  index[2] = kbp;

  p = (double **) MultiSet(qk_array, index, NULL);
  if (*p) {
    return (*p);
  } else if (qk_mode != QKDETAIL) {
    return NULL;
  }

  nqk = n_egrid*(n_egrid+1)/2;
  *p = (double *) malloc(sizeof(double)*n_tegrid*nqk);
  qk = *p;

  for (j = 0; j < nqk*n_tegrid; j++) qk[j] = 0.0;
  j = 0;
  for (ie1 = 0; ie1 < n_egrid; ie1++) {
    for (ie2 = 0; ie2 <= ie1; ie2++) {
      for (k = 0; k <= pw_scratch.max_k; k += 2) {	
	CIRadialQk(qi, ie1, ie2, kb, kbp, k);
	for (ite = 0; ite < n_tegrid; ite++) {
	  t = nqk*ite + j;
	  qi[ite] /= k + 1.0;
	  qk[t] += qi[ite];
	}
	if (k > 2 && fabs(qi[0]/qk[j]) < pw_scratch.tolerence) break;
      }
      j++;
    }
  }

  return qk;
}
  
int IonizeStrength(double *qkc, double *te, int b, int f) {
  int i, ip, j, ierr;
  LEVEL *lev1, *lev2;
  ANGULAR_ZFB *ang;
  double c, r, qke[NPARAMS];
  int nz, j0, j0p, kb, kbp;

  lev1 = GetLevel(b);
  lev2 = GetLevel(f);
  *te = lev2->energy - lev1->energy;
  if (*te <= 0) return -1;
  
  nz = AngularZFreeBound(&ang, f, b);
  if (nz <= 0) return -1;

  for (j = 0; j < NPARAMS; j++) {
    qkc[j] = 0.0;
  }

  for (i = 0; i < nz; i++) {
    kb = ang[i].kb;
    j0 = GetOrbital(kb)->kappa;
    j0 = GetJFromKappa(j0);
    for (ip = 0; ip <= i; ip++) {
      kbp = ang[ip].kb;
      j0p = GetOrbital(kbp)->kappa;
      j0p = GetJFromKappa(j0p);
      if (j0p != j0) continue;
      c = ang[i].coeff*ang[ip].coeff/(j0+1.0);
      if (ip != i) {
	c *= 2.0;
      }
      ierr = CIRadialQkIntegrated(qke, (*te), kb, kbp);
      if (ierr < 0) continue;
      for (j = 0; j < NPARAMS; j++) {
	qkc[j] += c * qke[j];
      }
    }
  }

  for (j = 0; j < NPARAMS; j++) {
    qkc[j] *= 16.0;
  }
  
  free(ang);

  return 0;
}

int SaveIonization(int nb, int *b, int nf, int *f, char *fn) {
  int i, j, k;
  int j1, j2, ie;
  FILE *file;
  LEVEL *lev1, *lev2;
  double delta, emin, emax, e, e0;
  double qkc[NPARAMS], qkb[NPARAMS], s;
  
  file = fopen(fn, "w");
  if (!file) return -1;

  if (n_tegrid == 0) {
    n_tegrid = 3;
  }

  emin = 1E10;
  emax = 1E-10;
  k = 0;
  for (i = 0; i < nb; i++) {
    lev1 = GetLevel(b[i]);
    for (j = 0; j < nf; j++) {
      lev2 = GetLevel(f[j]);
      e = lev2->energy - lev1->energy;
      if (e > 0) k++;
      if (e < emin && e > 0) emin = e;
      if (e > emax) emax = e;
    }
  }
  if (k == 0) {
    printf("No ionizations can occur\n");
    return 0;
  }
  if (tegrid[0] < 0.0) {
    e = 2.0*(emax-emin)/(emax+emin);
    if (e < 0.2) {
      SetIEGrid(1, 0.5*(emin+emax), emax);
    } else if (e < 0.5) {
      SetIEGrid(2, emin, emax);
    } else {
      if (k == 2) n_tegrid = 2;
      SetIEGrid(n_tegrid, emin, emax);
    }
  }

  e = 0.5*(emin + emax);
  emin = 0.05*e;
  emax = 8.0*e;
  egrid_type = 1;
  pw_type = 0;
  if (usr_egrid_type < 0) usr_egrid_type = 1;

  if (n_usr > 0 && usr_egrid[0] < 0.0) {
    SetUsrCIEGrid(n_usr, emin, emax, e);
  }    

  if (n_egrid == 0) {    
    n_egrid = 6;
  }
  if (egrid[0] < 0.0) {
    SetCIEGrid(n_egrid, emin, emax, e);
  }

  for (ie = 0; ie < n_egrid; ie++) {
    for (i = 0; i < n_tegrid; i++) {
      xegrid[i][ie] = egrid[ie]/tegrid[i];
      if (egrid_type == 1) xegrid[i][ie] += 1.0;
      log_xegrid[i][ie] = log(xegrid[i][ie]);
    }
    sigma[ie] = 1.0/sqrt(xegrid[0][ie]);
  }
  yegrid0[0] = 0.0;
  delta = 0.5/(NINT-1.0);
  for (i = 1; i < NINT; i++) {
    yegrid0[i] = yegrid0[i-1] + delta;
  }

  if (pw_scratch.nkl == 0) {
    SetCIPWGrid(0, NULL, NULL);
  }
  
  e = GetResidualZ();
  PrepCoulombBethe(n_egrid, n_tegrid, n_egrid, e, egrid, tegrid, egrid,
		   pw_scratch.nkl, pw_scratch.kl, egrid_type, pw_type, 0);
  
  fprintf(file, " IEGRID:\n   ");
  for (i = 0; i < n_tegrid; i++) {
    fprintf(file, "%10.4E ", tegrid[i]*HARTREE_EV);
  }
  fprintf(file, "\n");

  fprintf(file, " EGRID:\n   ");
  for (i = 0; i < n_egrid; i++) {
    fprintf(file, "%10.4E ", egrid[i]*HARTREE_EV);
  }
  fprintf(file, "\n\n");
  fprintf(file, " Strength = A*ln(u) + B*(1-1/u)^2 + C*ln(u)/u + D*(1-1/u)\n");
  if (output_format >= 0 && n_usr >= 0) {	
    if (usr_egrid_type == 0) fprintf(file, " Incident Electron UsrEGrid ");
    else fprintf(file, " Scattered Electron UsrEGrid ");
    fprintf(file, "\n");
  }
  fprintf(file, "\n");

  fprintf(file, "Bound 2J   Free  2J   Delta_E\n");

  for (i = 0; i < nb; i++) {
    j1 = LevelTotalJ(b[i]);
    for (j = 0; j < nf; j++) {
      j2 = LevelTotalJ(f[j]);
      k = IonizeStrength(qkc, &e, b[i], f[j]);
      if (k < 0) continue;
      fprintf(file, "%-5d %-2d   %-5d %-2d   %10.4E \n",
	      b[i], j1, f[j], j2, e*HARTREE_EV);
      for (k = 0; k < NPARAMS; k++) {
	fprintf(file, "%11.4E ", qkc[k]);
      }
      fprintf(file, "\n\n");
      if (output_format >= 0 && n_usr > 0) {	
	for (ie = 0; ie < n_usr; ie++) {
	  xusr[ie] = usr_egrid[ie]/e;
	  if (usr_egrid_type == 1) xusr[ie] += 1.0;
	  log_xusr[ie] = log(xusr[ie]);
	  CIRadialQkIntegratedBasis(4, qkb, xusr[ie], log_xusr[ie]);
	  s = 0;
	  for (k = 0; k < NPARAMS; k++) s += qkb[k]*qkc[k];
	  fprintf(file, "%-10.3E ", 
		  usr_egrid[ie]*HARTREE_EV);
	  if (output_format != 2) {
	    fprintf(file, "%-10.3E ", s);
	  } 
	  if (output_format != 1) {
	    e0 = usr_egrid[ie];
	    if (usr_egrid_type == 1) e0 += e;
	    s *= AREA_AU20/(2*e0*(j1+1.0));
	    fprintf(file, "%-10.3E ", s);
	  }
	  fprintf(file, "\n");
	}
	fprintf(file, "\n");
      }      
    }
  }

  fclose(file);
  return 0;
}

void _FreeIonizationQk(void *p) {
  double *dp;
  dp = *((double **) p);
  free(dp);
}

int FreeIonizationQk() {
  ARRAY *b;
  
  b = qk_array->array;
  if (b == NULL) return 0;
  MultiFreeData(b, qk_array->ndim, _FreeIonizationQk);
  return 0;
}

int InitIonization() {
  int blocks[3] = {2, 10, 10};
  int ndim = 3;
  
  qk_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(qk_array, sizeof(double *), ndim, blocks);

  n_egrid = 0;
  n_tegrid = 0;
  n_usr = 0;
  egrid[0] = -1.0;
  tegrid[0] = -1.0;
  usr_egrid[0] = -1.0;
  output_format = 0;

  return 0;
}

