#include "ionization.h"

#define NINT 31

static int output_format = 0;
static int egrid_type = -1;
static int usr_egrid_type = -1;
static int pw_type = -1;

static int n_usr = 0;
static double usr_egrid[MAXNUSR];
static double log_usr[MAXNUSR];
static int n_egrid = 0;
static double egrid[MAXNE];
static double log_egrid[MAXNE];
static int n_tegrid = 0;
static double tegrid[MAXNTE];
static double log_te[MAXNTE];
static double _dwork[17*MAXNE*(MAXNE+1)/2];
static int _iwork[39*MAXNE*(MAXNE+1)/2];

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

int SetCIEGridType(int utype, int etype) {
  if (utype >= 0) usr_egrid_type = utype;
  if (etype >= 0) egrid_type = etype;
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

double *CIRadialQk(int ie1, int ie2, int kb, int kbp, int k) {
  int index[4];
  double pk[MAXNTE][MAXNKL];
  double y2[MAXNKL];
  double **p, *qk, e1, e2, e0, te;
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
  index[0] = ie1;
  index[1] = ie2;
  index[2] = kb;
  index[3] = ko2;

  p = (double **) MultiSet(qk_array, index, NULL);
  if (*p) {
    return (*p);
  }

  for (i = 0; i < n_tegrid; i++) {
    for (j = 0; j < pw_scratch.nkl0; j++) {
      pk[i][j] = 0.0;
    }
  }
  (*p) = (double *) malloc(sizeof(double)*n_tegrid);
  qk = *p;
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
  return qk;
}

int IonizeStrength(double *s, double *te, int b, int f) {
  int i, ip, j, k, t, nt, np, ndp, md, ierr;
  LEVEL *lev1, *lev2;
  double c, r, r0, scale;
  int nz, j0, j0p, kb, kbp;
  double rq[MAXNE*(MAXNE+1)/2];
  double xrq[MAXNE*(MAXNE+1)/2];
  double yrq[MAXNE*(MAXNE+1)/2];
  double *qk, *x, x0; 
  double e, delta, ratio, e1, e2, ehalf;
  double xr[NINT], yr[NINT], integrand[NINT];
  ANGULAR_ZFB *ang;
  int ie1, ie2;

  lev1 = GetLevel(b);
  lev2 = GetLevel(f);
  *te = lev2->energy - lev1->energy;
  if (*te <= 0) return -1;
  
  nz = AngularZFreeBound(&ang, f, b);
  if (nz <= 0) return -1;

  i = 0;
  for(ie1 = 0; ie1 < n_egrid; ie1++) {
    for (ie2 = 0; ie2 <= ie1; ie2++) {
      r = log(1.0 + egrid[ie1]/(egrid[ie2]+(*te)));
      xrq[i] = r;
      r = log(1.0 + egrid[ie2]/(*te));
      yrq[i] = r;
      rq[i++] = 0.0;      
    }
  }
  ndp = i;

  np = 3;
  nt = 1;
  x = log_te;
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
      t = 0;
      for(ie1 = 0; ie1 < n_egrid; ie1++) {
	for (ie2 = 0; ie2 <= ie1; ie2++) {
	  r0 = 0.0;
	  for (k = 0; k <= pw_scratch.max_k; k += 2) {
	    qk = CIRadialQk(ie1, ie2, kb, kbp, k);
	    if (n_tegrid == 1) {
	      r = qk[0];
	    } else {
	      x0 = log(*te);
	      uvip3p_(&np, &n_tegrid, x, qk, &nt, &x0, &r);
	    }
	    r = r/(k+1.0);
	    r0 += r;
	    if (k > 2 && fabs(r/r0) < pw_scratch.tolerence) break;
	  }
	  rq[t++] += c*r0;
	}
      }      
    }
  }

  for (j = 0; j < n_usr; j++) {
    e = usr_egrid[j];
    if (usr_egrid_type == 0) {
      e -= (*te);
    }
    ehalf = e*0.5;
    yr[0] = 0;
    nt = NINT-1;
    yr[nt] = log(1.0 + ehalf/(*te));
    delta = (yr[nt] - yr[0])/(nt);
    for (t = 1; t < nt; t++) {
      yr[t] = yr[t-1] + delta;
    }
    nt = NINT;    
    for (t = 0; t < nt; t++) {
      e2 = (*te)*(exp(yr[t])-1.0);
      e1 = e - e2;
      xr[t] = log(1.0 + e1/(e2 + (*te)));
    }
    if (j == 0) md = 1;
    else md = 3;
    sdbi3p_(&md, &ndp, xrq, yrq, rq, &nt, xr, yr, 
	    integrand, &ierr, _dwork, _iwork);
    for (t = 0; t < nt; t++) {
      e2 = (*te)*(exp(yr[t])-1.0);
      integrand[t] *= e2 + (*te);
    }
    r0 = Simpson(integrand, 0, nt-1);
    r0 *= delta;
    s[j] = 16.0*r0;
  }
  
  free(ang);

  return 0;
}

int SaveIonization(int nb, int *b, int nf, int *f, char *fn) {
  int i, j, k;
  int j1, j2, ie;
  FILE *file;
  LEVEL *lev1, *lev2;
  double emin, emax, e, s[MAXNUSR];
  
  if (n_usr == 0) {
    printf("No ionization energy specified \n");
    return -1;
  }

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
    if (e < 0.1) {
      SetIEGrid(1, 0.5*(emin+emax), emax);
    } else if (e < 0.3) {
      SetIEGrid(2, emin, emax);
    } else {
      if (k == 2) n_tegrid = 2;
      SetIEGrid(n_tegrid, emin, emax);
    }
  }

  e = 0.5*(emin + emax);
  emin = 0.1*e;
  emax = 8.0*e;
  egrid_type = 1;
  pw_type = 0;
  if (usr_egrid_type < 0) usr_egrid_type = 1;

  if (n_usr == 0) {
    n_usr = 6;
  }
  if (usr_egrid[0] < 0.0) {
    if (n_egrid > n_usr && egrid_type == usr_egrid_type) {
      SetUsrCIEGridDetail(n_egrid, egrid);
    } else {
      SetUsrCIEGrid(n_usr, emin, emax, e);
    }
  }
  
  if (n_egrid == 0) {    
    n_egrid = 6;
  }
  if (egrid[0] < 0.0) {
    if (usr_egrid_type == 0) 
      emax = Min(usr_egrid[n_usr-1]-e, emax);
    else 
      emax = Min(usr_egrid[n_usr-1], emax);
    SetCIEGrid(n_egrid, emin, emax, e);
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
  if (usr_egrid_type == 0) fprintf(file, " Incident Electron UsrEGrid, ");
  else fprintf(file, " Scattered Electron UsrEGrid, ");
  if (output_format != 2) fprintf(file, "Ionization Strength, ");
  if (output_format != 1) fprintf(file, "Cross Section");
  fprintf(file, "\n\n");

  fprintf(file, "Bound 2J\tFree  2J\tDelta_E\n");

  for (i = 0; i < nb; i++) {
    j1 = LevelTotalJ(b[i]);
    for (j = 0; j < nf; j++) {
      j2 = LevelTotalJ(f[j]);
      k = IonizeStrength(s, &e, b[i], f[j]);
      if (k < 0) continue;
      fprintf(file, "%-5d %-2d\t%-5d %-2d\t%10.4E\n",
	      b[i], j1, f[j], j2, e*HARTREE_EV);
      for (ie = 0; ie < n_usr; ie++) {
	fprintf(file, "%-10.3E ", 
		usr_egrid[ie]*HARTREE_EV);
	if (output_format != 2) {
	  fprintf(file, "%-10.3E ", s[ie]);
	} 
	if (output_format != 1) {
	  s[ie] = AREA_AU20*(s[ie])/(2*(e+usr_egrid[ie])*(j1+1.0));
	  fprintf(file, "%-10.3E ", s[ie]);
	}
	fprintf(file, "\n");
      }
      fprintf(file, "\n");
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
  int blocks[4] = {8, 8, 10, 6};
  int ndim = 4;
  
  qk_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(qk_array, sizeof(double *), ndim, blocks);
  n_egrid = 0;
  n_tegrid = 0;
  n_usr = 0;
  egrid[0] = -1.0;
  tegrid[0] = -1.0;
  usr_egrid[0] = -1.0;
  return 0;
}

