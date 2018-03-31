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

#include "transition.h"
#include <time.h>

static char *rcsid="$Id$";
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
  double eps0;
  double eps;
} transition_option = {DGAUGE, DMODE, ERANK, MRANK, TRCUT0, TRCUT};

typedef struct {
  TR_RECORD r;
  TR_EXTRA rx;
  int ks[2];
} TR_DATUM;

int SetTransitionCut(double c0, double c) {
  if (c0 >= 0) {
    transition_option.eps0 = c0;
  } else {
    transition_option.eps0 = TRCUT0;
  }
  if (c >= 0) {
    transition_option.eps = c;
  } else {
    transition_option.eps = TRCUT;
  }
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
  if (idatum->s[0].index < 0 || idatum->s[3].index >= 0) {
    free(idatum->bra);
    free(idatum);
    return -1;
  }

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
  if (!Triangle(j1, j2, m2)) {
    free(idatum->bra);    
    free(idatum);
    return -1;
  }

  rx->energy = (lev2->energy - lev1->energy);

  rx->energy += ConfigEnergyShift(ns, idatum->bra, ia, ib, m2);
  rx->sdev = sqrt(ConfigEnergyVariance(ns, idatum->bra, ia, ib, m2));

  aw = FINE_STRUCTURE_CONST * rx->energy;
  if (aw < 0.0) {
    free(idatum->bra);
    free(idatum);
    return -1;
  }
  
  if (transition_option.mode == M_NR && m != 1) {
    r = MultipoleRadialNR(m, k0, k1, transition_option.gauge);
  } else {
    r = MultipoleRadialFR(aw, m, k0, k1, transition_option.gauge);
  }

  *strength = sqrt((lev1->ilev+1.0)*q1*(j2+1.0-q2)/((j1+1.0)*(j2+1.0)))*r;
  
  free(idatum->bra);
  free(idatum);
  return 0;
}

int TRMultipole(double *strength, double *energy,
		int m, int lower, int upper) {
  int m0, m1, m2;
  int p1, p2, j1, j2;
  LEVEL *lev1, *lev2;
  double s, r, a, aw, *mbk, tr;
  int nz, i, nmk;
  ANGULAR_ZMIX *ang;

  lev1 = GetLevel(lower);
  if (lev1 == NULL) return -1;
  lev2 = GetLevel(upper);
  if (lev2 == NULL) return -1;
  if (*energy <= 0.0) {
    *energy = lev2->energy - lev1->energy;
  }
  if (*energy <= 0.0) return -1;
  aw = FINE_STRUCTURE_CONST * (*energy);
  if (aw < 0.0) return -1;

  DecodePJ(lev1->pj, &p1, &j1);
  DecodePJ(lev2->pj, &p2, &j2);
  if (j1 == 0 && j2 == 0) return -1;

  if (m != 0) {
    m2 = 2*abs(m);    
    if (!Triangle(j1, j2, m2)) return -1;
    if (m > 0 && IsEven(p1+p2+m)) return -1;
    if (m < 0 && IsOdd(p1+p2-m)) return -1;    
    
    s = 0.0;
    nz = AngularZMix(&ang, lower, upper, m2, m2, &nmk, &mbk);
    if (nz <= 0 && nmk < m2/2) {
      if (nmk > 0) free(mbk);
      return -1;
    }
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
    if (nmk >= m2/2) {
      r = mbk[m2/2-1];
      if (transition_option.gauge == G_COULOMB && m < 0) {
	r /= aw;
      }
      a = r/s;
      a *= a;
      if (a < 0.75) {
	s += r;
      }
    }
    if (nz > 0) {
      free(ang);	
    }  
    if (nmk > 0) {
      free(mbk);
    }

    *strength = s;
  } else {
    m0 = abs(j1-j2);
    if (m0 == 0) m0 += 2;
    m1 = (j1+j2);
    if (m0 > transition_option.max_m && m0 > transition_option.max_e) {
      return -1;
    }
    tr = 0.0;
    nz = AngularZMix(&ang, lower, upper, m0, m1, &nmk, &mbk);
    for (m2 = m0; m2 <= m1; m2 += 2) {      
      m = m2/2;
      if (IsEven(p1+p2+m)) m = -m;
      s = 0.0;
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
      if (nmk >= m2/2) {
	r = mbk[m2/2-1];
	if (transition_option.gauge == G_COULOMB && m < 0) {
	  r /= aw;
	}
	a = r/s;
	a *= a;
	if (a < 0.75) {
	  s += r;
	}
      }
      r = OscillatorStrength(m, *energy, s, &a);
      tr += r;
      if (tr > 0 && r/tr < transition_option.eps0) break;
    }
    if (nz > 0) {
      free(ang);
    }
    if (nmk > 0) {
      free(mbk);
    }
    *strength = tr;
  }
  return 0;
}

int TRMultipoleEB(double *strength, double *energy, int m, int lower, int upper) {
  LEVEL *lev1, *lev2;
  LEVEL *plev1, *plev2;
  int i1, i2, j1, j2, p1, p2, k, m2;
  int ilev1, ilev2, mlev1, mlev2, q;
  double r, a, c;

  lev1 = GetEBLevel(lower);
  if (lev1 == NULL) return -1;
  lev2 = GetEBLevel(upper);
  if (lev2 == NULL) return -1;
  
  *energy = lev2->energy - lev1->energy;
  if (*energy <= 0.0) return -1;
  
  m2 = 2*abs(m);

  for (q = 0; q <= m2; q++) strength[q] = 0.0;
  for (i1 = 0; i1 < lev1->n_basis; i1++) {
    if (lev1->mixing[i1] == 0) continue;
    DecodeBasisEB(lev1->basis[i1], &ilev1, &mlev1);
    plev1 = GetLevel(ilev1);
    DecodePJ(plev1->pj, &p1, &j1);
    for (i2 = 0; i2 < lev2->n_basis; i2++) {
      if (lev2->mixing[i2] == 0) continue;
      c = lev1->mixing[i1]*lev2->mixing[i2];
      DecodeBasisEB(lev2->basis[i2], &ilev2, &mlev2);
      plev2 = GetLevel(ilev2);
      DecodePJ(plev2->pj, &p2, &j2);
      k = TRMultipole(&r, energy, m, ilev1, ilev2);
      if (k != 0) continue;
      a = W3j(j1, m2, j2, -mlev1, mlev1-mlev2, mlev2);
      q = (mlev1-mlev2)/2+abs(m);
      if (a == 0) continue;
      if (IsOdd((j1-mlev1)/2)) a = -a;
      strength[q] += c*r*a;
      /*
      printf("%d %d %d %d %2d %2d %2d %10.3E %10.3E %10.3E %10.3E %10.3E\n",
	     lower, upper, ilev1, ilev2, mlev1, mlev2, q-1,c, a, r, c*a*r, 
	     strength[q]);
      */
    }
  }

  return 0;
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
  
  i1 = (int *) p1;
  i2 = (int *) p2;
  lev1 = GetLevel(*i1);
  lev2 = GetLevel(*i2);
  c1 = GetConfigFromGroup(lev1->iham, lev1->pb);
  c2 = GetConfigFromGroup(lev2->iham, lev2->pb);
  return CompareNRConfig(c1, c2);
}

static int CompareTRDatum(const void *p1, const void *p2) {
  TR_DATUM *r1, *r2;

  r1 = (TR_DATUM *) p1;
  r2 = (TR_DATUM *) p2;
  if ((r1->r).upper < (r2->r).upper) return -1;
  else if ((r1->r).upper > (r2->r).upper) return 1;
  else {
    if ((r1->r).lower < (r2->r).lower) return -1;
    else if ((r1->r).lower > (r2->r).lower) return 1;
    else return 0;
  }
}

int SaveTransitionEB0(int nlow, int *low, int nup, int *up, 
		      char *fn, int m) {
  int k, i, j, nq;
  double emin, emax, e0, s[101], et;
  F_HEADER fhdr;
  TRF_HEADER tr_hdr;
  TRF_RECORD r;
  LEVEL *lev1, *lev2;
  TFILE *f;
  
  if (nlow <= 0 || nup <= 0) return -1;  
  
  if (m == 1 || transition_option.mode == M_FR) {
    k = 0;
    emin = 1E10;
    emax = 1E-10;
    for (i = 0; i < nlow; i++) {
      lev1 = GetEBLevel(low[i]);
      for (j = 0; j < nup; j++) {
	lev2 = GetEBLevel(up[j]);
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
    if (e0 < EPS3) {
      SetAWGrid(1, 0.5*(emin+emax), emax);
    } else if (e0 < 1.0) {
      SetAWGrid(2, emin, emax);
    } else {
      SetAWGrid(3, emin, emax);
    }
  }
  fhdr.type = DB_TRF;
  strcpy(fhdr.symbol, GetAtomicSymbol());
  fhdr.atom = GetAtomicNumber();  
  lev1 = GetEBLevel(low[0]);
  DecodeBasisEB(lev1->pb, &i, &j);  
  tr_hdr.nele = GetNumElectrons(i);
  tr_hdr.multipole = m;
  tr_hdr.gauge = GetTransitionGauge();
  if (m == 1) { /* always FR for M1 transitions */
    tr_hdr.mode = M_FR;
  } else {
    tr_hdr.mode = GetTransitionMode();
  }
  nq = 2*abs(m) + 1;
  GetFields(&tr_hdr.bfield, &tr_hdr.efield, &tr_hdr.fangle);
    
  f = OpenFile(fn, &fhdr);
  InitFile(f, &fhdr, &tr_hdr);

  ResetWidMPI();
#pragma omp parallel default(shared) private(i, j, k, e0, r, s, et)
  {
  r.strength = (float *) malloc(sizeof(float)*nq);
  for (j = 0; j < nup; j++) {
    for (i = 0; i < nlow; i++) {
      int skip = SkipMPI();
      if (skip) continue;
      k = TRMultipoleEB(s, &et, m, low[i], up[j]);
      if (k != 0) continue;
      e0 = 0.0;
      for (k = 0; k < nq; k++) {
	r.strength[k] = s[k];
	if (s[k]) e0 = s[k];
      }
      if (e0 == 0.0) continue;
      r.lower = low[i];
      r.upper = up[j];
      WriteTRFRecord(f, &r);
    }
  }
  free(r.strength);
  }
  DeinitFile(f, &fhdr);
  CloseFile(f, &fhdr);

  return 0;
}
      
int SaveTransition0(int nlow, int *low, int nup, int *up, 
		    char *fn, int m) {
  int i, j, k, jup;
  TFILE *f;
  LEVEL *lev1, *lev2;
  TR_RECORD r;
  TR_EXTRA rx;
  TR_HEADER tr_hdr;
  F_HEADER fhdr;
  double *s, *et, *a, trd, gf;
  double e0, emin, emax;
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
  
  if (nlow <= 0 || nup <= 0) return -1;
  
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
  if (e0 < EPS3) {
    SetAWGrid(1, 0.5*(emin+emax), emax);
  } else if (e0 < 1.0) {
    SetAWGrid(2, emin, emax);
  } else {
    SetAWGrid(3, emin, emax);
  }
  
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
    ResetWidMPI();
#pragma omp parallel default(shared) private(imin, imax, jmin, jmax, lev1, lev2, c0, c1, ir, ntr, rd, ep, em, e0, wp, wm, w0, i, j, ic0, ic1, k, ir0, gf, j0, j1, nrs0, nrs1, de, cm, cp)
    {
    imin = 0;
    for (ic0 = 0; ic0 < nic0; ic0++) {
      imax = nc0[ic0];
      jmin = 0;
      lev1 = GetLevel(low[imin]);
      c0 = GetConfigFromGroup(lev1->iham, lev1->pb);
      for (ic1 = 0; ic1 < nic1; ic1++) {
	int skip = SkipMPI();
	if (skip) continue;
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
	      gf = OscillatorStrength(m, rd[ir].rx.energy, 
				      rd[ir].r.strength, NULL);
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
    }
  } else {
    //PrepAngZStates(nlow, low, nup, up);
    ResetWidMPI();
#pragma omp parallel default(shared) private(a, s, et, j, jup, trd, i, k, gf, r)
    {
      a = malloc(sizeof(double)*nlow);
      s = malloc(sizeof(double)*nlow);
      et = malloc(sizeof(double)*nlow);
      for (j = 0; j < nup; j++) {
	int skip = SkipMPI();
	if (skip) continue;
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
	  if (a[i] <= 0 || a[i] < (transition_option.eps * trd)) continue;
	  r.lower = low[i];
	  r.strength = s[i];
	  WriteTRRecord(f, &r, NULL);
	}
      }      
      free(a);
      free(s);
      free(et);
    }
  }
  DeinitFile(f, &fhdr);
  CloseFile(f, &fhdr);

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

int OverlapLowUp(int nlow, int *low, int nup, int *up) {
  int i, j, n;
  int *lowinup, *upinlow, *icom;

  lowinup = (int *) malloc(sizeof(int)*nlow);
  upinlow = (int *) malloc(sizeof(int)*nup);

  for (i = 0; i < nlow; i++) {
    lowinup[i] = -1;
  }
  for (i = 0; i < nup; i++) {
    upinlow[i] = -1;
  }
  qsort(low, nlow, sizeof(int), CompareInt);
  if (up != low) {
    qsort(up, nup, sizeof(int), CompareInt);
  }
  for (i = 0; i < nlow; i++) {
    lowinup[i] = IBisect(low[i], nup, up);    
    if (lowinup[i] >= 0) {
      upinlow[lowinup[i]] = i;
    }
  }
  icom = (int *) malloc(sizeof(int)*nlow);
  n = 0;
  for (i = 0; i < nlow; i++) {
    if (lowinup[i] >= 0) icom[n++] = low[i];
  }
  j = 0;
  for (i = 0; i < nlow; i++) {
    if (lowinup[i] < 0) {
      low[j++] = low[i];
    }
  }
  for (i = 0; i < n; i++) {
    low[j++] = icom[i];
  }
  j = 0;
  for (i = 0; i < nup; i++) {
    if (upinlow[i] < 0) {
      up[j++] = up[i];
    }
  }
  for (i = 0; i < n; i++) {
    up[j++] = icom[i];
  }
  
  free(lowinup);
  free(upinlow);
  free(icom);
  /*
#if USE_MPI == 2
  int mr, nr;
  mr = MPIRank(&nr);
  if (nr > 1) {
    RandIntList(nup-n, up);
    RandIntList(nlow-n, low);
    if (n > 1) {
      RandIntList(n, up+nup-n);
      RandIntList(n, low+nlow-n);
    }
  }
#endif
  */
  return n;
}
  
int SaveTransition(int nlow, int *low, int nup, int *up,
		   char *fn, int m) {
  int n, *alev, i, nc;
  
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
  if (nlow <= 0 || nup <= 0) return -1;

  nc = OverlapLowUp(nlow, low, nup, up);
  SaveTransition0(nc, low+nlow-nc, nc, up+nup-nc, fn, m);
  SaveTransition0(nc, low+nlow-nc, nup-nc, up, fn, m);
  SaveTransition0(nup-nc, up, nc, low+nlow-nc, fn, m);
  SaveTransition0(nlow-nc, low, nup, up, fn, m);
  SaveTransition0(nup, up, nlow-nc, low, fn, m);

  if (n > 0) free(alev);
  ReinitRadial(1);

  return 0;
}
  
int GetLowUpEB(int *nlow, int **low, int *nup, int **up, 
	       int nlow0, int *low0, int nup0, int *up0) {  
  int i, j, ilev, mlev, n;
  LEVEL *lev;
 
  n = GetNumEBLevels();
  if (n == 0) return -1;

  *low = malloc(sizeof(int)*n);
  *up = malloc(sizeof(int)*n);
  *nlow = 0;
  *nup = 0;
  for (i = 0; i < n; i++) {
    lev = GetEBLevel(i);
    DecodeBasisEB(lev->pb, &ilev, &mlev);
    for (j = 0; j < nlow0; j++) {
      if (low0[j] == ilev) {
	(*low)[(*nlow)++] = i;	
	break;
      }      
    }
    for (j = 0; j < nup0; j++) {
      if (up0[j] == ilev) {
	(*up)[(*nup)++] = i;	
	break;
      }
    }
  }

  return 0;
}

int SaveTransitionEB(int nlow0, int *low0, int nup0, int *up0,
		     char *fn, int m) {
  int n, nlow, *low, nup, *up, nc;

  n = GetLowUpEB(&nlow, &low, &nup, &up, nlow0, low0, nup0, up0);
  if (n == -1) return 0;

  nc = OverlapLowUp(nlow, low, nup, up);
  SaveTransitionEB0(nc, low+nlow-nc, nc, up+nup-nc, fn, m);
  SaveTransitionEB0(nc, low+nlow-nc, nup-nc, up, fn, m);
  SaveTransitionEB0(nup-nc, up, nc, low+nlow-nc, fn, m);
  SaveTransitionEB0(nlow-nc, low, nup, up, fn, m);
  SaveTransitionEB0(nup, up, nlow-nc, low, fn, m);

  free(low);
  free(up);
  ReinitRadial(1);

  return 0;
}

int PolarizeCoeff(char *ifn, char *ofn, int i0, int i1) {
  TFILE *f1;
  FILE *f2;
  int n, i, t, tp, k, q, s, sp, m, mp, m2;
  double a, c, e;
  F_HEADER fh;
  TRF_HEADER h;
  TRF_RECORD r;
  EN_SRECORD *mem_en_table;
  int swp, mem_en_table_size;
  
  mem_en_table = GetMemENFTable(&mem_en_table_size);
  if (mem_en_table == NULL) {
    printf("Energy table has not been built in memory.\n");
    return -1;
  }

  f1 = OpenFileRO(ifn, &fh, &swp);
  if (f1 == NULL) {
    printf("cannot open file %s\n", ifn);
    return -1;
  }  

  if (fh.type != DB_TRF || fh.nblocks == 0) {
    printf("File %s is not of DB_TRF type\n", ifn);
    goto DONE;
  }
  
  if (strcmp(ofn, "-") == 0) {
    f2 = stdout;
  } else {
    f2 = fopen(ofn, "w");
  }
  if (f2 == NULL) {
    printf("cannot open file %s\n", ofn);
    goto DONE;
  }
    
  while (1) {
    n = ReadTRFHeader(f1, &h, swp);
    if (n == 0) break;
    m2 = 2*abs(h.multipole);
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadTRFRecord(f1, &r, swp, &h);
      if ((r.lower == i0 || i0 < 0) && (r.upper == i1 || i1 < 0)) {
	e = mem_en_table[r.upper].energy - mem_en_table[r.lower].energy;
	e = FINE_STRUCTURE_CONST*e;
	e = e*e*e;
	e *= RATE_AU/(4.0*PI);
	for (t = -1; t <= 1; t += 2) {
	  for (tp = -1; tp <= 1; tp += 2) {
	    for (k = 0; k <= m2; k++) {
	      for (q = -k; q <= k; q++) {
		c = 0.0;
		for (s = 0; s <= m2; s++) {
		  m = s - m2/2;
		  mp = m + q;
		  sp = mp + m2/2;
		  if (sp < 0 || sp > m2) continue;
		  a = r.strength[s]*r.strength[sp];
		  a *= W3j(m2, m2, 2*k, 2*m, -2*mp, 2*q);
		  a *= W3j(m2, m2, 2*k, 2*t, -2*tp, 2*(tp-t));
		  if (IsOdd(abs(mp-tp))) a = -a;
		  c += a;
		}
		c *= 2.0*k + 1.0;
		c *= e;
		fprintf(f2, "%4d %4d %2d %2d %2d %2d %2d %15.8E\n",
			r.lower, r.upper, t, tp, k, q, tp-t, c);
	      }
	    }
	  }
	}
      }
      free(r.strength);
    }
  }

 DONE:
  FCLOSE(f1);
  
  if (f2) {
    if (f2 != stdout) {
      fclose(f2);
    } else {
      fflush(f2);
    }
  }

  return 0;
}
  
