#include "transition.h"
#include <time.h>

static char *rcsid="$Id: transition.c,v 1.28 2005/01/10 22:06:49 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

/* the following options controll the computation methods of OS.
   gauge = 1 coulomb gauge (velocity form)
         = 2 babushkin gauge (length form)

   mode = 0 use relativistic expression for radial integrals.
        = 1 use non-relativistic approximation.
   
   max_e, the maximum rank of electric multipoles.
   max_m, the maximum rank of magnetic multipoles.
*/
static struct {
  int gauge;
  int mode;
  int max_e;
  int max_m;
  double eps;
} transition_option = {DGAUGE, DMODE, ERANK, MRANK, TRCUT};

typedef struct {
  TR_RECORD r;
  TR_EXTRA rx;
  int ks[2];
} TR_DATUM;

int SetTransitionCut(double c) {
  transition_option.eps = c;
  return 0;
}

double GetTransitionCut(void) {
  return transition_option.eps;
}

void SetTransitionMode(int m) {
  transition_option.mode = m;
}

void SetTransitionGauge(int m) {
  transition_option.gauge = m;
}

void SetTransitionMaxE(int m) {
  transition_option.max_e = m;
}

void SetTransitionMaxM(int m) {
  transition_option.max_m = m;
}

void SetTransitionOptions(int gauge, int mode, 
			  int max_e, int max_m) {
  transition_option.gauge = gauge;
  transition_option.mode = mode;
  transition_option.max_e = max_e;
  transition_option.max_m = max_m;
}

int GetTransitionGauge(void) {
  return transition_option.gauge;
}

int GetTransitionMode(void) {
  return transition_option.mode;
}

int TRMultipoleUTA(double *strength, TR_EXTRA *rx, 
		   int m, int lower, int upper, int *ks) {
  int m2, ns, k0, k1, q1, q2;
  int p1, p2, j1, j2, ia, ib;
  LEVEL *lev1, *lev2;
  double r, aw;
  INTERACT_DATUM *idatum;
  
  *strength = 0.0;
  lev1 = GetLevel(lower);
  if (lev1 == NULL) return -1;
  lev2 = GetLevel(upper);
  if (lev2 == NULL) return -1;

  p1 = lev1->pj;
  p2 = lev2->pj;
  if (m > 0 && IsEven(p1+p2+m)) return -1;
  if (m < 0 && IsOdd(p1+p2+m)) return -1;
  
  idatum = NULL;
  ns = GetInteract(&idatum, NULL, NULL, lev1->iham, lev2->iham,
		   lev1->pb, lev2->pb, 0, 0, 0);
  if (ns <= 0) return -1;
  if (idatum->s[0].index < 0 || idatum->s[3].index >= 0) return -1;

  if (idatum->s[0].nq_bra > idatum->s[0].nq_ket) {
    ia = ns-1-idatum->s[0].index;
    ib = ns-1-idatum->s[1].index;
    j1 = idatum->s[0].j;
    j2 = idatum->s[1].j;
    q1 = idatum->s[0].nq_bra;
    q2 = idatum->s[1].nq_bra;
    k0 = OrbitalIndex(idatum->s[0].n, idatum->s[0].kappa, 0.0);
    k1 = OrbitalIndex(idatum->s[1].n, idatum->s[1].kappa, 0.0);
    PackNRShell(ks, idatum->s[0].n, idatum->s[0].kl, q1);
    if (idatum->s[0].kappa > 0) ks[0] |= 0x01000000;
    PackNRShell(ks+1, idatum->s[1].n, idatum->s[1].kl, q2);
    if (idatum->s[1].kappa > 0) ks[1] |= 0x01000000;
  } else {
    ia = ns-1-idatum->s[1].index;
    ib = ns-1-idatum->s[0].index;    
    j1 = idatum->s[1].j;
    j2 = idatum->s[0].j;
    q1 = idatum->s[1].nq_bra;
    q2 = idatum->s[0].nq_bra;
    k1 = OrbitalIndex(idatum->s[0].n, idatum->s[0].kappa, 0.0);
    k0 = OrbitalIndex(idatum->s[1].n, idatum->s[1].kappa, 0.0);
    PackNRShell(ks, idatum->s[1].n, idatum->s[1].kl, q1);
    if (idatum->s[1].kappa > 0)  ks[0] |= 0x01000000;
    PackNRShell(ks+1, idatum->s[0].n, idatum->s[0].kl, q2);
    if (idatum->s[0].kappa > 0)  ks[1] |= 0x01000000;
  }

  m2 = 2*abs(m);
  if (!Triangle(j1, j2, m2)) return -1;

  rx->energy = (lev2->energy - lev1->energy);
  rx->energy += ConfigEnergyShift(ns, idatum->bra, ia, ib, m2);
  rx->sdev = sqrt(ConfigEnergyVariance(ns, idatum->bra, ia, ib, m2));
  aw = FINE_STRUCTURE_CONST * rx->energy;
  if (aw < 0.0) return -1;
  
  if (transition_option.mode == M_NR && m != 1) {
    r = MultipoleRadialNR(m, k0, k1, transition_option.gauge);
  } else {
    r = MultipoleRadialFR(aw, m, k0, k1, transition_option.gauge);
  }

  *strength = sqrt((lev1->ilev+1.0)*q1*(j2+1.0-q2)/((j1+1.0)*(j2+1.0)))*r;

  return 0;
}

int TRMultipole(double *strength, double *energy,
		int m, int lower, int upper) {
  int m2;
  int p1, p2, j1, j2;
  LEVEL *lev1, *lev2;
  double s, r, aw;
  int nz, i;
  ANGULAR_ZMIX *ang;

  *strength = 0.0;
  lev1 = GetLevel(lower);
  if (lev1 == NULL) return -1;
  lev2 = GetLevel(upper);
  if (lev2 == NULL) return -1;
  if (*energy <= 0.0) {
    *energy = lev2->energy - lev1->energy;
  }
  if (*energy <= 0.0) return -1;

  DecodePJ(lev1->pj, &p1, &j1);
  DecodePJ(lev2->pj, &p2, &j2);

  m2 = 2*abs(m);

  if (!Triangle(j1, j2, m2)) return -1;
  if (m > 0 && IsEven(p1+p2+m)) return -1;
  if (m < 0 && IsOdd(p1+p2-m)) return -1;

  aw = FINE_STRUCTURE_CONST * (*energy);
  if (aw < 0.0) return -1;

  s = 0.0;
  
  nz = AngularZMix(&ang, lower, upper, m2, m2);
    
  if (nz <= 0) return -1;

  for (i = 0; i < nz; i++) {
    if (ang[i].k != m2) continue;
    if (transition_option.mode == M_NR && m != 1) {
      r = MultipoleRadialNR(m, ang[i].k0, ang[i].k1, 
			    transition_option.gauge);
    } else {
      r = MultipoleRadialFR(aw, m, ang[i].k0, ang[i].k1,
			    transition_option.gauge);
    }
    s += r * ang[i].coeff;
  }
  free(ang);	  
  
  *strength = s;

  return 0;
}

int GetLowestMultipole(int p1, int j1, int p2, int j2) {
  int m;

  if (j1 == 0 && j2 == 0) return 0;

  m = abs(j1-j2);
  if (IsOdd(m)) return 0;
  m = m/2;

  if (m == 0) m = 1;
  
  if (IsEven(p1+m+p2)) m = -m;

  return m;
}

static int CompareNRConfig(const void *p1, const void *p2) {
  CONFIG *c1, *c2;

  c1 = (CONFIG *) p1;
  c2 = (CONFIG *) p2;
  if (c1->nnrs > c2->nnrs) return 1;
  else if (c1->nnrs < c2->nnrs) return -1;
  else {
    return memcmp(c1->nrs, c2->nrs, sizeof(int)*c1->nnrs);
  }
}

static int CompareNRLevel(const void *p1, const void *p2) {
  int *i1, *i2;
  LEVEL *lev1, *lev2;
  CONFIG *c1, *c2;
  
  i1 = p1;
  i2 = p2;
  lev1 = GetLevel(*i1);
  lev2 = GetLevel(*i2);
  c1 = GetConfigFromGroup(lev1->iham, lev1->pb);
  c2 = GetConfigFromGroup(lev2->iham, lev2->pb);
  return CompareNRConfig(c1, c2);
}

static int CompareTRDatum(const void *p1, const void *p2) {
  TR_DATUM *r1, *r2;

  r1 = p1;
  r2 = p2;
  if ((r1->r).upper < (r2->r).upper) return -1;
  else if ((r1->r).upper > (r2->r).upper) return 1;
  else {
    if ((r1->r).lower < (r2->r).lower) return -1;
    else if ((r1->r).lower > (r2->r).lower) return 1;
    else return 0;
  }
}

int SaveTransition(int nlow, int *low, int nup, int *up, 
		   char *fn, int m) {
  int i, j, k, n, jup;
  FILE *f;
  LEVEL *lev1, *lev2;
  TR_RECORD r;
  TR_EXTRA rx;
  TR_HEADER tr_hdr;
  F_HEADER fhdr;
  double *s, *et, *a, trd, gf;
  double e0, emin, emax;
  int *alev;
  int ic0, ic1, nic0, nic1, *nc0, *nc1, j0, j1, ntr;
  int imin, imax, jmin, jmax, nrs0, nrs1, ir, ir0;
  double ep, em, wp, wm, w0, de, cp, cm;
  CONFIG *c0, *c1;
  TR_DATUM *rd;
  int mj = 0xFF000000, mn = 0xFFFFFF;

#ifdef PERFORM_STATISTICS
  STRUCT_TIMING structt;
  ANGULAR_TIMING angt;
  RECOUPLE_TIMING recouplet;
  RAD_TIMING radt;
#endif
  
  n = 0;
  if (nlow == 0 || nup == 0) {
    n = GetNumLevels();
    if (n <= 0) return -1;
    alev = malloc(sizeof(int)*n);
    if (!alev) return -1;
    
    for (i = 0; i < n; i++) alev[i] = i;

    if (nlow == 0) {
      nlow = n; 
      low = alev;
    }
    if (nup == 0) {
      nup = n;
      up = alev;
    }
  }
  
  if (m == 1 || transition_option.mode == M_FR) {
    k = 0;
    emin = 1E10;
    emax = 1E-10;
    for (i = 0; i < nlow; i++) {
      lev1 = GetLevel(low[i]);
      for (j = 0; j < nup; j++) {
	lev2 = GetLevel(up[j]);
	e0 = lev2->energy - lev1->energy;
	if (e0 > 0) k++;
	if (e0 < emin && e0 > 0) emin = e0;
	if (e0 > emax) emax = e0;
      }
    }
      
    if (k == 0) {
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
    
  if (nlow <= 0 || nup <= 0) return -1;
  
  fhdr.type = DB_TR;
  strcpy(fhdr.symbol, GetAtomicSymbol());
  fhdr.atom = GetAtomicNumber();
  tr_hdr.nele = GetNumElectrons(low[0]);
  tr_hdr.multipole = m;
  tr_hdr.gauge = GetTransitionGauge();
  if (m == 1) { /* always FR for M1 transitions */
    tr_hdr.mode = M_FR;
  } else {
    tr_hdr.mode = GetTransitionMode();
  }
  f = OpenFile(fn, &fhdr);
  InitFile(f, &fhdr, &tr_hdr);
    
  if (IsUTA()) {
    qsort(low, nlow, sizeof(int), CompareNRLevel);
    nc0 = malloc(sizeof(int)*nlow);    
    ic0 = 0;
    for (i = 0; i < nlow; i++) {
      lev1 = GetLevel(low[i]);
      c1 = GetConfigFromGroup(lev1->iham, lev1->pb);
      if (i > 0 && CompareNRConfig(c1, c0)) {
	nc0[ic0++] = i;
      }
      c0 = c1;
    }
    nc0[ic0] = nlow;
    nic0 = ic0+1;
    if (up != low) {
      qsort(up, nup, sizeof(int), CompareNRLevel);
      nc1 = malloc(sizeof(int)*nup);    
      ic1 = 0;
      for (i = 0; i < nup; i++) {
	lev1 = GetLevel(up[i]);
	c1 = GetConfigFromGroup(lev1->iham, lev1->pb);
	if (i > 0 && CompareNRConfig(c1, c0) != 0) {
	  nc1[ic1++] = i;
	}
	c0 = c1;
      }
      nc1[ic1] = nup;
      nic1 = ic1+1;
    } else {
      nc1 = nc0;
      nic1 = nic0;
    }
    imin = 0;
    for (ic0 = 0; ic0 < nic0; ic0++) {
      imax = nc0[ic0];
      jmin = 0;
      lev1 = GetLevel(low[imin]);
      c0 = GetConfigFromGroup(lev1->iham, lev1->pb);
      for (ic1 = 0; ic1 < nic1; ic1++) {
	jmax = nc1[ic1];
	lev2 = GetLevel(up[jmin]);
	c1 = GetConfigFromGroup(lev2->iham, lev2->pb);
	ir = 0;
	ntr = (jmax-jmin)*(imax-imin);
	rd = malloc(sizeof(TR_DATUM)*ntr);
	ep = 0.0;
	em = 0.0;
	e0 = 0.0;
	wp = 0.0;
	wm = 0.0;
	w0 = 0.0;
	ir0 = -1;
	for (i = imin; i < imax; i++) {
	  for (j = jmin; j < jmax; j++) {
	    k = TRMultipoleUTA(&gf, &(rd[ir].rx), m, low[i], up[j], rd[ir].ks);
	    if (k != 0) {
	      rd[ir].r.lower = -1;
	      rd[ir].r.upper = -1;
	      ir++;
	      continue;
	    }
	    ir0 = ir;
	    rd[ir].r.lower = low[i];
	    rd[ir].r.upper = up[j];
	    rd[ir].r.strength = gf;
	    rd[ir].rx.sci = 1.0;
	    if (m == -1) {
	      gf = OscillatorStrength(m, rd[ir].rx.energy, rd[ir].r.strength, NULL);
	      j0 = rd[ir].ks[0]&mj;
	      j1 = rd[ir].ks[1]&mj;
	      if (j0==0 && j1==0) {
		wp += gf;
		ep += gf*rd[ir].rx.energy;
	      } else if (j0 && j1) {
		wm += gf;
		em += gf*rd[ir].rx.energy;
	      }
	      e0 += gf*rd[ir].rx.energy;
	      w0 += gf;
	    }
	    ir++;
	  }
	}	
	if (wm > 0.0 && ir0 >= 0) {
	  nrs0 = 0;
	  nrs1 = 0;
	  for (i = 0; i < c0->nnrs; i++) {
	    if (c0->nrs[i]>>8 == (rd[ir0].ks[0]&mn)>>8) nrs0 = c0->nrs[i];
	    else if (c0->nrs[i]>>8 == (rd[ir0].ks[1]&mn)>>8) nrs1 = c1->nrs[i];
	  }
	  if (nrs0 == 0) {
	    nrs0 = rd[ir0].ks[0] & mn;
	  }
	  if (nrs1 == 0) {
	    nrs1 = rd[ir0].ks[1] & mn;
	  }
	  ep /= wp;
	  em /= wm;
	  e0 /= w0;
	  de = ConfigEnergyShiftCI(nrs0, nrs1);
	  cm = 1.0 + de/(e0 - ep);
	  cp = 1.0 + de/(e0 - em);
	  if (cm < EPS3) {
	    cp = (wp+wm)/wp;
	    cm = 0.0;
	  } else if (cp < EPS3) {
	    cm = (wp+wm)/wm;
	    cp = 0.0;
	  }
	  ir = 0;
	  for (i = imin; i < imax; i++) {
	    for (j = jmin; j < jmax; j++) {
	      if (rd[ir].r.lower < 0) {
		ir++;
		continue;
	      }
	      j0 = rd[ir].ks[0]&mj;
	      j1 = rd[ir].ks[1]&mj;
	      if (j0==0 && j1==0) {
		rd[ir].rx.sci = cp;
	      } else if (j0 && j1) {
		rd[ir].rx.sci = cm;
	      }
	      ir++;
	    }
	  }
	}
	qsort(rd, ntr, sizeof(TR_DATUM), CompareTRDatum);
	ir = 0;
	for (ir = 0; ir < ntr; ir++) {
	  if (rd[ir].r.lower < 0) {
	    continue;
	  }
	  WriteTRRecord(f, &(rd[ir].r), &(rd[ir].rx));
	}
	free(rd);
	jmin = jmax;
      }
      imin = imax;
    }    
    free(nc0);
    if (up != low) free(nc1);
  } else {
    a = malloc(sizeof(double)*nlow);
    s = malloc(sizeof(double)*nlow);
    et = malloc(sizeof(double)*nlow);
    for (j = 0; j < nup; j++) {
      jup = LevelTotalJ(up[j]);
      trd = 0.0;
      for (i = 0; i < nlow; i++) {
	a[i] = 0.0;
	et[i] = 0.0;
	k = TRMultipole(s+i, et+i, m, low[i], up[j]);
	if (k != 0) continue;
	gf = OscillatorStrength(m, et[i], s[i], &(a[i]));
	a[i] /= jup+1.0;
	trd += a[i];
      } 
      if (trd < 1E-30) continue;
      r.upper = up[j];
      for (i = 0; i < nlow; i++) {
	if (a[i] < (transition_option.eps * trd)) continue;
	r.lower = low[i];
	r.strength = s[i];
	WriteTRRecord(f, &r, NULL);
      }
    }
    
    free(a);
    free(s);
    free(et);
  }

  DeinitFile(f, &fhdr);
  CloseFile(f, &fhdr);
  ReinitRadial(1);
  if (n > 0) {
    free(alev);
  }

#ifdef PERFORM_STATISTICS
  GetStructTiming(&structt);
  fprintf(perform_log, "AngZMix: %6.1E, AngZFB: %6.1E, AngZxZFB: %6.1E, SetH: %6.1E DiagH: %6.1E\n",
	  ((double) (structt.angz_mix))/CLOCKS_PER_SEC,
	  ((double) (structt.angz_fb))/CLOCKS_PER_SEC,
	  ((double) (structt.angzxz_fb))/CLOCKS_PER_SEC,
	  ((double) (structt.set_ham))/CLOCKS_PER_SEC,
	  ((double) (structt.diag_ham))/CLOCKS_PER_SEC);
  fprintf(perform_log, "AngZS: %7d %7d %6.1E, %6.1E, AngZFBS: %6.1E, AngZxZFBS: %6.1E, AddZ: %6.1E, AddZxZ: %6.1E\n",
	  structt.n_angz_states, structt.n_angz_states_load,
	  ((double) (structt.angz_states))/CLOCKS_PER_SEC, 
	  ((double) (structt.angz_states_load))/CLOCKS_PER_SEC,
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
