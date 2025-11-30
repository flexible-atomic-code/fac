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

#include <time.h>
#include "transition.h"
#include "mbpt.h"
#include "global.h"

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
  int k0, k1;
} TR_DATUM;

typedef struct {
  float strength;
  float energy;
  float sdev;
  float sci;
  int k0;
  int k1;
} TR_UTA;

static int _ufmin = 0;
static int _ufmax = -1;
static int _nuf_tr = 0;
static int _uumin = 0;
static int _uumax = -1;
static int _nuu_tr = 0;
static TR_UTA *_ufuu_tr = NULL;

static int _sfu_nmax = 0;
static int _sfu_ke = 0;
static double *_efu = NULL;

static int _tr_all = 0;
static double _tr_aw = -1;
static double _lower_emin = 0.0;
static double _lower_emax = 0.0;
static double _upper_emin = 0.0;
static double _upper_emax = 0.0;
static double _tr_emin = 0.0;
static double _tr_emax = 0.0;
static int _uta_e1 = 0;
static int _progress_report = 0;

int OutOfERange(double e0, double e1, double de) {
  if (_lower_emin > 0 && e0 < _lower_emin) return 1;
  if (_lower_emax > 0 && e0 > _lower_emax) return 1;
  if (_upper_emin > 0 && e1 < _upper_emin) return 1;
  if (_upper_emax > 0 && e1 > _upper_emax) return 1;
  if (_tr_emin > 0 && de < _tr_emin) return 1;
  if (_tr_emax > 0 && de > _tr_emax) return 1;
  return 0;
}

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

int TRMultipoleUTA0I(double *strength, TR_EXTRA *rx, INTERACT_DATUM *idatum,
		     int m0, int k0, int k1, int p1, int p2, int j1, int j2,
		     int ns, int *ks, int ia, int ib,
		     double te, int q1, int q2, double w0, double w1) {
  int m, m2;
  double aw, qw1, qw2, r;

  *strength = 0.0;
  if (m0 == 0) {
    m = abs(j1-j2)/2;
    if (m == 0) m += 1;
    if (IsEven(p1+p2+m)) m = -m;
  } else {
    m = m0;
  }
  m2 = 2*abs(m);
  if ((_uta_e1 && m != -1) || !Triangle(j1, j2, m2)) {
    return 0;
  }
  rx->energy = te;
  aw = FINE_STRUCTURE_CONST * rx->energy;
  if (aw < 0.0) {
    return 0;
  }
  
  if (ks) {
    rx->energy += ConfigEnergyShift(ns, idatum->bra, ia, ib, m2);
    r = ConfigEnergyVariance(ns, idatum->bra, ia, ib, m2);
    rx->sdev = sqrt(r);
  } else {
    rx->sdev = 0.0;
  }
  rx->sci = 1.0;
  
  if (transition_option.mode == M_NR && m != 1) {
    r = MultipoleRadialNR(m, k0, k1, transition_option.gauge);
  } else {
    r = MultipoleRadialFR(aw, m, k0, k1, transition_option.gauge);
  }
  /*
  qw1 = q1 * (j1+1.0)/w0;
  qw2 = q2 * (j2+1.0)/w1;
  *strength = sqrt(qw1*(j2+1.0-qw2)/((j1+1.0)*(j2+1.0)))*r;
  */
  *strength = sqrt(q1*(w1-q2)/(w0*w1))*r;
  if (m0 == 0) {
    *strength = OscillatorStrength(m, rx->energy, *strength, &r);
  }
  return m;
}

int TRMultipoleUTA0(double *strength, TR_EXTRA *rx, 
		    int m, int kg0, int kg1, int kc0, int kc1,
		    double te, int p1, int p2, int *ks, int *pk0, int *pk1) {
  int m2, ns, k0, k1, q1, w2, q2, m0, m1, u0, u1;
  int j1, j2, ia, ib, i0, i1, t0, t1, i, t, u;
  double r, eg, w0, w1;
  CONFIG *c, *c1;
  INTERACT_DATUM *idatum;
  INTERACT_SHELL *s;

  m0 = m;
  *pk0 = -1;
  *pk1 = -1;
  *strength = 0.0;
  rx->energy = 0.0;
  rx->sdev = 0.0;
  rx->sci = 1.0;
  c = NULL;
  if (_uta_e1) {
    c = GetConfigFromGroup(kg0, kc0);
    if (c->parity >= 0 && c1->parity >= 0) {
      if (c->parity == c1->parity) return 0;
    }
  }
  
  idatum = NULL;
  ns = GetInteract(&idatum, NULL, NULL, kg0, kg1, kc0, kc1, 0, 0, 0);
  if (ns <= 0) return 0;  
  if (idatum->s[0].index < 0 || idatum->s[3].index >= 0) {
    free(idatum->bra);
    free(idatum);
    return 0;
  }

  if (c == NULL) c = GetConfigFromGroup(kg0, kc0);
  
  s = idatum->s;
  k0 = -1;
  k1 = -1;
  
  if (s[0].nq_bra > s[0].nq_ket) {
    ia = ns-1-s[0].index;
    ib = ns-1-s[1].index;
    j1 = s[0].j;
    j2 = s[1].j;
    q1 = s[0].nq_bra;
    q2 = s[1].nq_bra;
    w2 = FactorNR(c, s[0].n, s[0].kappa);
    if (w2 > 0) q1 = w2;
    if (s[0].nr < 2 && s[1].nr < 2) {
      k0 = OrbitalIndex(s[0].n, s[0].kappa, 0.0);
      k1 = OrbitalIndex(s[1].n, s[1].kappa, 0.0);
      if (ks) {
	PackNRShell(ks, s[0].n, s[0].kl, q1);
	if (s[0].kappa > 0) ks[0] |= 0x01000000;
	PackNRShell(ks+1, s[1].n, s[1].kl, q2);
	if (s[1].kappa > 0) ks[1] |= 0x01000000;
      }
    } else {
      ks = NULL;
    }
  } else {
    ia = ns-1-s[1].index;
    ib = ns-1-s[0].index;    
    j1 = s[1].j;
    j2 = s[0].j;
    q1 = s[1].nq_bra;
    q2 = s[0].nq_bra;
    w2 = FactorNR(c, s[1].n, s[1].kappa);
    if (w2 > 0) q1 = w2;
    if (s[0].nr < 2 && s[1].nr < 2) {
      k1 = OrbitalIndex(s[0].n, s[0].kappa, 0.0);
      k0 = OrbitalIndex(s[1].n, s[1].kappa, 0.0);
      if (ks) {
	PackNRShell(ks, s[1].n, s[1].kl, q1);
	if (s[1].kappa > 0)  ks[0] |= 0x01000000;
	PackNRShell(ks+1, s[0].n, s[0].kl, q2);
	if (s[0].kappa > 0)  ks[1] |= 0x01000000;
      }
    } else {
      ks = NULL;
    }
  }

  if (s[0].nr == 2) {
    if (s[0].kappa >= 0) {
      i1 = -(1+s[0].kappa/2);
      i0 = -(i1+1);
      u0 = i0;
      w0 = 2*(s[0].kappa+1.0);
    } else {
      i0 = 0;
      i1 = 2*s[0].n-2;
      u0 = -1;
      w0 = 2*s[0].n*s[0].n;
    }
  } else {
    i0 = s[0].kappa;
    i1 = s[0].kappa;
    u0 = i0;
    w0 = 2*abs(s[0].kappa);
  }
  if (s[1].nr == 2) {
    if (s[1].kappa >= 0) {
      t1 = -(1+s[1].kappa/2);
      t0 = -(t1+1);
      u1 = t0;
      w1 = 2*(s[1].kappa+1.0);
    } else {
      t0 = 0;
      t1 = 2*s[1].n-2;
      u1 = -1;
      w1 = 2*s[1].n*s[1].n;
    }
  } else {
    t0 = s[1].kappa;
    t1 = s[1].kappa;
    u1 = t0;
    w1 = 2*abs(s[1].kappa);
  }
  if (s[0].nr < 2 && s[1].nr < 2) {
    if (p1 < 0) {
      p1 = s[0].kl/2;
    }
    if (p2 < 0) {
      p2 = s[1].kl/2;
    }
    m = TRMultipoleUTA0I(strength, rx, idatum, m0, k0, k1, p1, p2, j1, j2,
			 ns, ks, ia, ib, te, q1, q2, w0, w1);
  } else {
    m = 0;
    u = u1;
    for (i = abs(i0); i <= abs(i1); i++) {
      GetJLFromKappa(u0, &j1, &p1);
      k0 = OrbitalIndex(s[0].n, u0, 0.0);
      p1 /= 2;
      u1 = u;
      for (t = abs(t0); t <= abs(t1); t++) {
	GetJLFromKappa(u1, &j2, &p2);
	k1 = OrbitalIndex(s[1].n, u1, 0.0);
	p2 /= 2;
	m1 = TRMultipoleUTA0I(&r, rx, idatum, m0, k0, k1, p1, p2, j1, j2,
			      ns, ks, ia, ib, te, q1, q2, w0, w1);
	if (m1) {
	  *strength += r;
	  if (m == 0) m = m1;
	  else if (m < 0) {
	    if (m1 < 0 && m1 > m) m = m1;
	  } else {
	    if (m1 < 0) m = m1;
	    else if (m1 < m) m = m1;
	  }
	}
	if (u1 < 0) u1 = -u1;
	else u1 = -(u1+1);
      }
      if (u0 < 0) u0 = -u0;
      else u0 = -(u0+1);
    }    
  }
  free(idatum->bra);
  free(idatum);
  if (ks && m == -1) {
    *pk0 = k0;
    *pk1 = k1;
  }
  return m;
}

int TRMultipoleUTA(double *strength, TR_EXTRA *rx, 
		   int m, int lower, int upper, int *ks, int *pk0, int *pk1) {  
  LEVEL *lev1, *lev2;
  SYMMETRY *sym;
  STATE *s, *v;
  CONFIG *cfg;
  int kg0, kc0, kg1, kc1, iuu, iuf, u;
  int r, t, j, p, k, mk, m0, p1, p2, k0, k1;
  double wb, si, te, wm, eg, w1, w2, xe, xe2, xs, xc;

  m0 = m;
  *strength = 0.0;
  *pk0 = -1;
  *pk1 = -1;
  lev1 = GetLevel(lower);
  if (lev1 == NULL) return 0;
  lev2 = GetLevel(upper);
  if (lev2 == NULL) return 0;
  if (lev1->nele != lev2->nele) return 0;
  eg = EGroundIon(lev1->nele);
  if (OutOfERange(lev1->energy-eg, lev2->energy-eg,
		  lev2->energy-lev1->energy)) return 0;
  p1 = lev1->pj;
  p2 = lev2->pj;
  if (p1 >= 0 && p2 >= 0) {
    if (m > 0 && IsEven(p1+p2+m)) return 0;
    if (m < 0 && IsOdd(p1+p2+m)) return 0;
  }
  
  te = lev2->energy - lev1->energy;
  if (te <= 0) return 0;

  kg0 = -1;
  kg1 = -1;
  kc0 = -1;
  kc1 = -1;
  mk = 0;
  if (lev1->n_basis > 0 && lev2->n_basis > 0) {
    if (TransUTA() == 0) {
      k = TRMultipole(strength, &te, m, lower, upper, pk0, pk1);
      mk = (k == 0);
    } else {
      sym = GetSymmetry(lev1->pj);
      DecodePJ(lev1->pj, &p, &j);
      *strength = 0.0;
      wm = 0.0;
      for (t = 0; t < lev1->n_basis; t++) {
	s = (STATE *) ArrayGet(&(sym->states), lev1->basis[t]);
	w1 = lev1->mixing[t]*lev1->mixing[t];
	for (r = 0; r < lev2->n_basis; r++) {
	  v = (STATE *) ArrayGet(&(sym->states), lev2->basis[r]);
	  w2 = lev2->mixing[r]*lev2->mixing[r];
	  if (s->kgroup != kg0 || v->kgroup != kg1 ||
	      s->kcfg != kc0 || v->kcfg != kc1) {
	    k = TRMultipoleUTA0(&si, rx, m, s->kgroup, v->kgroup,
				s->kcfg, v->kcfg, te, p1, p2, NULL, &k0, &k1);
	    kg0 = s->kgroup;
	    kg1 = v->kgroup;
	    kc0 = s->kcfg;
	    kc1 = v->kcfg;
	  }
	  if (k == 0) continue;
	  wb = (j+1.0)*w1*w2;
	  if (wb > wm) {
	    wm = wb;
	    mk = k;
	    *pk0 = k0;
	    *pk1 = k1;
	  }
	  if (m0 != 0) si *= si;
	  *strength += wb*si;
	}
      }
      if (m0 != 0) *strength = sqrt(*strength);
    }
    rx->energy = te;
    rx->sdev = 0.0;
    rx->sci = 1.0;
    return mk;
  }
  CONFIG *cc = GetConfigFromGroup(lev2->iham, lev2->pb);
  if (lev1->n_basis == 0 && lev2->n_basis == 0) {
    k = TRMultipoleUTA0(strength, rx, m, lev1->iham, lev2->iham,
			lev1->pb, lev2->pb, te, p1, p2, ks, pk0, pk1);
    wb = lev1->ilev + 1.0;
    if ((lower >= 0 && upper >= 0)||(lower <= 0 && upper <= 0)) {
      if (m0 != 0) wb = sqrt(wb);
      *strength *= wb;
    }
    rx->sci = 1.0;
    return k;
  } else if (lev1->n_basis > 0) {
    sym = GetSymmetry(lev1->pj);
    DecodePJ(lev1->pj, &p, &j);
    *strength = 0.0;
    wm = 0.0;
    xe = 0.0;
    xe2 = 0.0;
    xs = 0.0;
    xc = 0.0;
    for (t = 0; t < lev1->n_basis; t++) {
      s = (STATE *) ArrayGet(&(sym->states), lev1->basis[t]);
      /*
      if (s->kgroup != kg0 || lev2->iham != kg1 ||
	  s->kcfg != kc0 || lev2->pb != kc1) {
	k = TRMultipoleUTA0(&si, rx, m, s->kgroup, lev2->iham,
			    s->kcfg, lev2->pb, te, p1, p2, NULL);
	kg0 = s->kgroup;
	kg1 = lev2->iham;
	kc0 = s->kcfg;
	kc1 = lev2->pb;
      }
      if (k == 0) continue;
      */
      k = -1;
      cfg = GetConfig(s);
      iuf = cfg->iulev - _ufmin;      
      iuu = upper - _uumin;
      u = iuu*_nuf_tr + iuf;
      si = _ufuu_tr[u].strength;
      wb = (j+1.0)*lev1->mixing[t]*lev1->mixing[t];
      if (wb > wm) {
	wm = wb;
	mk = k;
	*pk0 = _ufuu_tr[u].k0;
	*pk1 = _ufuu_tr[u].k1;
      }
      if (m0 != 0) si *= si;
      wb *= si;
      *strength += wb;      
      xe += wb*_ufuu_tr[u].energy;
      xe2 += wb*_ufuu_tr[u].energy*_ufuu_tr[u].energy;
      xs += wb*_ufuu_tr[u].sdev*_ufuu_tr[u].sdev;
      xc += wb*_ufuu_tr[u].sci;
    }
    wb = *strength;
    if (wb > 0) {
      if (m0 != 0) *strength = sqrt(*strength);
      xe /= wb;
      xe2 /= wb;
      xs /= wb;
      xc /= wb;
      xe2 = xe2 - xe*xe;
      if (xe2 < 0) xe2 = 0.0;
    }      
    rx->energy = xe;
    rx->sdev = sqrt(xs + xe2);
    rx->sci = xc;
  } else if (lev2->n_basis > 0) {
    sym = GetSymmetry(lev2->pj);
    DecodePJ(lev2->pj, &p, &j);
    *strength = 0.0;
    wm = 0.0;
    xe = 0.0;
    xe2 = 0.0;
    xs = 0.0;
    xc = 0.0;
    for (t = 0; t < lev2->n_basis; t++) {
      s = (STATE *) ArrayGet(&(sym->states), lev2->basis[t]);
      cfg = GetConfig(s);
      /*
      if (lev1->iham != kg0 || s->kgroup != kg1 ||
	  lev1->pb != kc0 || s->kcfg != kc1) {
	k = TRMultipoleUTA0(&si, rx, m, lev1->iham, s->kgroup,
			    lev1->pb, s->kcfg, te, p1, p2, NULL);
	kg0 = lev1->iham;
	kg1 = s->kgroup;
	kc0 = lev1->pb;
	kc1 = s->kcfg;
      }
      if (k == 0) continue;
      */
      k = -1;
      iuf = cfg->iulev - _ufmin;      
      iuu = lower - _uumin;
      u = iuu*_nuf_tr + iuf;
      si = _ufuu_tr[u].strength;
      wb = lev2->mixing[t]*lev2->mixing[t]*(lev1->ilev+1.0);
      wb *= (j+1.0)/fabs(cfg->sweight);
      if (wb > wm) {
	wm = wb;
	mk = k;
	*pk0 = _ufuu_tr[u].k0;
	*pk1 = _ufuu_tr[u].k1;
      }
      if (m0 != 0) si *= si;
      wb *= si;
      *strength += wb;      
      xe += wb*_ufuu_tr[u].energy;
      xe2 += wb*_ufuu_tr[u].energy*_ufuu_tr[u].energy;
      xs += wb*_ufuu_tr[u].sdev*_ufuu_tr[u].sdev;
      xc += wb*_ufuu_tr[u].sci;
    }
    wb = *strength;
    if (wb > 0) {
      if (m0 != 0) *strength = sqrt(*strength);
      xe /= wb;
      xe2 /= wb;
      xs /= wb;
      xc /= wb;
      xe2 = xe2 - xe*xe;
      if (xe2 < 0) xe2 = 0.0;
    }
    rx->energy = xe;
    rx->sdev = sqrt(xs + xe2);
    rx->sci = xc;
  }
  return mk;
}

// return multipoletype if the rate is present, as m may be passed in as 0
int TRMultipoleUTA1(double *strength, TR_EXTRA *rx, 
		    int m, int lower, int upper, int *ks) {
  int m2, ns, k0, k1, q1, q2, m0;
  int p1, p2, j1, j2, ia, ib;
  LEVEL *lev1, *lev2;
  double r, aw, eg;
  INTERACT_DATUM *idatum;

  m0 = m;
  *strength = 0.0;
  lev1 = GetLevel(lower);
  if (lev1 == NULL) return 0;
  lev2 = GetLevel(upper);
  if (lev2 == NULL) return 0;

  if (lev1->nele != lev2->nele) return 0;
  eg = EGroundIon(lev1->nele);
  if (OutOfERange(lev1->energy-eg, lev2->energy-eg,
		  lev2->energy-lev1->energy)) return 0;
  p1 = lev1->pj;
  p2 = lev2->pj;
  if (m > 0 && IsEven(p1+p2+m)) return 0;
  if (m < 0 && IsOdd(p1+p2+m)) return 0;
  
  idatum = NULL;
  ns = GetInteract(&idatum, NULL, NULL, lev1->iham, lev2->iham,
		   lev1->pb, lev2->pb, 0, 0, 0);
  if (ns <= 0) return 0;
  if (idatum->s[0].index < 0 || idatum->s[3].index >= 0) {
    free(idatum->bra);
    free(idatum);
    return 0;
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

  if (m == 0) {
    m = abs(j1-j2)/2;
    if (m == 0) m += 1;
    if (IsEven(p1+p2+m)) m = -m;
  }
  m2 = 2*abs(m);
  if (!Triangle(j1, j2, m2)) {
    free(idatum->bra);    
    free(idatum);
    return 0;
  }

  rx->energy = (lev2->energy - lev1->energy);
  rx->energy += ConfigEnergyShift(ns, idatum->bra, ia, ib, m2);
  rx->sdev = sqrt(ConfigEnergyVariance(ns, idatum->bra, ia, ib, m2));

  aw = FINE_STRUCTURE_CONST * rx->energy;
  if (aw < 0.0) {
    free(idatum->bra);
    free(idatum);
    return 0;
  }
  
  if (transition_option.mode == M_NR && m != 1) {
    r = MultipoleRadialNR(m, k0, k1, transition_option.gauge);
  } else {
    r = MultipoleRadialFR(aw, m, k0, k1, transition_option.gauge);
  }
  *strength = sqrt((lev1->ilev+1.0)*q1*(j2+1.0-q2)/((j1+1.0)*(j2+1.0)))*r;
  if (m0 == 0) {
    *strength = OscillatorStrength(m, rx->energy, *strength, &r);
  }
  free(idatum->bra);
  free(idatum);
  return m;
}

void AddSFU(int md, int k0, int k1, int ke, double e, double r) {
  if (_sfu_nmax > 0 && k0 >= 0 && k1 >= 0) {
    int n0, n1, i;
    n0 = GetOrbital(k0)->n;
    n1 = GetOrbital(k1)->n;
    if (n0 > 0 && n1 > 0 && n0 <= _sfu_nmax && n1 <= _sfu_nmax) {
      int nm1 = _sfu_nmax*_sfu_nmax;
      int nm2 = nm1*2;
      int nm3 = nm1*3;
      int nm4 = nm1*4;
      int nm5 = nm1*5;
      i = (n1-1)*_sfu_nmax + n0-1;
      if (md == 0) {
#pragma omp atomic
	_efu[i] += r*e;
#pragma omp atomic
	_efu[i+nm2] += r;
#pragma omp atomic
	_efu[i+nm4] += 1;
      } else {
#pragma omp atomic
	_efu[i+nm1] += r*e;
#pragma omp atomic
	_efu[i+nm3] += r;
#pragma omp atomic
	_efu[i+nm5] += 1;
      }
    }
  }
}

int TRMultipole(double *strength, double *energy,
		int m, int lower, int upper, int *pk0, int *pk1) {
  int m0, m1, m2;
  int p1, p2, j1, j2, n0, n1;
  LEVEL *lev1, *lev2;
  double s, r, a, aw, *mbk, tr, eg, wa, wm;
  int nz, i, nmk;
  ANGULAR_ZMIX *ang;

  *pk0 = -1;
  *pk1 = -1;
  lev1 = GetLevel(lower);
  if (lev1 == NULL) return -1;
  lev2 = GetLevel(upper);
  if (lev2 == NULL) return -1;
  
  if (lev1->nele != lev2->nele) return -1;
  eg = EGroundIon(lev1->nele);
  if (OutOfERange(lev1->energy-eg, lev2->energy-eg,
		  lev2->energy-lev1->energy)) return -1;
  
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
    n0 = -1;
    n1 = -1;
    for (i = 0; i < nz; i++) {
      if (ang[i].k != m2) continue;
      if (transition_option.mode == M_NR && m != 1) {
	r = MultipoleRadialNR(m, ang[i].k0, ang[i].k1, 
			      transition_option.gauge);
      } else {
	r = MultipoleRadialFR(aw, m, ang[i].k0, ang[i].k1,
			      transition_option.gauge);
      }
      wa = r*ang[i].coeff;
      s += wa;
      wa = fabs(wa);
      if (m == -1 && wa > wm) {
	wm = wa;
	n0 = ang[i].k0;
	n1 = ang[i].k1;
      }
    }
    if (nmk >= m2/2) {
      r = mbk[m2/2-1];
      if (transition_option.gauge == G_COULOMB && m < 0) {
	r /= aw;
      }
      a = r/s;
      a *= a;
      if (a < AngZCutMBPT()) {
	s += r;
      }
    }
    if (nz > 0) {
      free(ang);	
    }  
    if (nmk > 0) {
      free(mbk);
    }
    *pk0 = n0;
    *pk1 = n1;
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
    wm = 0.0;
    n0 = -1;
    n1 = -1;
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
	wa = r*ang[i].coeff;
	s += wa;
	wa = fabs(wa);
	if (m == -1 && wa > wm) {
	  wm = wa;
	  n0 = ang[i].k0;
	  n1 = ang[i].k1;
	}
      }
      if (nmk >= m2/2) {
	r = mbk[m2/2-1];
	if (transition_option.gauge == G_COULOMB && m < 0) {
	  r /= aw;
	}
	a = r/s;
	a *= a;
	if (a < AngZCutMBPT()) {
	  s += r;
	}
      }
      r = OscillatorStrength(m, *energy, s, &a);
      tr += r;
      if (tr > 0 && r/tr < transition_option.eps0) break;
    }
    *pk0 = n0;
    *pk1 = n1;
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
  int i1, i2, j1, j2, p1, p2, k, m2, k0, k1;
  int ilev1, ilev2, mlev1, mlev2, q;
  double r, a, c, eg;

  lev1 = GetEBLevel(lower);
  if (lev1 == NULL) return -1;
  lev2 = GetEBLevel(upper);
  if (lev2 == NULL) return -1;
  
  if (lev1->nele != lev2->nele) return -1;
  
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
      k = TRMultipole(&r, energy, m, ilev1, ilev2, &k0, &k1);
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
  SYMMETRY *sym;
  STATE *s;
  
  i1 = (int *) p1;
  i2 = (int *) p2;
  lev1 = GetLevel(*i1);
  lev2 = GetLevel(*i2);
  if (lev1->n_basis > 0) {
    sym = GetSymmetry(lev1->pj);
    s = (STATE *) ArrayGet(&(sym->states), lev1->pb);
    c1 = GetConfig(s);
  } else {
    c1 = GetConfigFromGroup(lev1->iham, lev1->pb);
  }
  if (lev2->n_basis > 0) {
    sym = GetSymmetry(lev2->pj);
    s = (STATE *) ArrayGet(&(sym->states), lev2->pb);
    c2 = GetConfig(s);
  } else {
    c2 = GetConfigFromGroup(lev2->iham, lev2->pb);
  }
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
  tr_hdr.nele = lev1->nele;
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
  SYMMETRY *sym;
  STATE *st;
  double *s, *et, *a, trd, gf;
  double e0, emin, emax;
  int iuta, auta, ic0, ic1, nic0, nic1, *nc0, *nc1, j0, j1, ntr;
  int imin, imax, jmin, jmax, nrs0, nrs1, ir, ir0, nele, k0, k1, nut;
  double ep, em, wp, wm, w0, de, cp, cm, eg;
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

  iuta = IsUTA();
  auta = iuta && UTAGrid();
  nele = GetNumElectrons(low[0]);
  eg = EGroundIon(nele);
  nut = GetNumULevels()+1;
  for (i = 0; i < nlow; i++) {
    lev1 = GetLevel(low[i]);
    for (j = 0; j < nup; j++) {
      lev2 = GetLevel(up[j]);
      if (lev1->nele != lev2->nele) continue;
      if (fn == NULL) {
	if (low[i] >= 0 && up[j] >= 0) continue;
      }
      e0 = lev2->energy - lev1->energy;
      if (OutOfERange(lev1->energy-eg, lev2->energy-eg, e0)) continue;
      if (e0 > 0) {
	k++;
	if (e0 < emin && e0 > 0) emin = e0;
	if (e0 > emax) emax = e0;
	if (lev1->n_basis > 0 || lev2->n_basis > 0) auta = 0;
      }
    }
  }
    
  if (k == 0) {
    return 0;
  }
    
  emin *= FINE_STRUCTURE_CONST;
  emax *= FINE_STRUCTURE_CONST;
  e0 = 2.0*(emax-emin)/(emin+emax);
    
  FreeMultipoleArray();
  if (auta) {
    SetAWGrid(-1, emin, emax);
  } else {
    if (_tr_aw >= 0) {
      SetAWGrid(1, _tr_aw, _tr_aw);
    } else {
      if (e0 < EPS3) {
	SetAWGrid(1, 0.5*(emin+emax), emax);
      } else if (e0 < 1.0) {
	SetAWGrid(2, emin, emax);
      } else {
	SetAWGrid(3, emin, emax);
      }
    }
  }
  if (fn) {
    fhdr.type = DB_TR;
    strcpy(fhdr.symbol, GetAtomicSymbol());
    fhdr.atom = GetAtomicNumber();
    tr_hdr.nele = nele;
    tr_hdr.multipole = m;
    tr_hdr.gauge = GetTransitionGauge();
    if (m == 1) { /* always FR for M1 transitions */
      tr_hdr.mode = M_FR;
    } else {
      tr_hdr.mode = GetTransitionMode();
    }
    f = OpenFile(fn, &fhdr);
    InitFile(f, &fhdr, &tr_hdr);
  }
  double tstart = WallTime();
  int nproc=0, *ntrans = NULL;
  if (_progress_report != 0) {
    ntrans = InitTransReport(&nproc);
  }
  if (iuta) {
    qsort(low, nlow, sizeof(int), CompareNRLevel);
    nc0 = malloc(sizeof(int)*nlow);    
    ic0 = 0;
    for (i = 0; i < nlow; i++) {
      lev1 = GetLevel(low[i]);
      if (lev1->n_basis > 0) {
	sym = GetSymmetry(lev1->pj);
	st = (STATE *) ArrayGet(&(sym->states), lev1->pb);
	c1 = GetConfig(st);
      } else {
	c1 = GetConfigFromGroup(lev1->iham, lev1->pb);
      }
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
	if (lev1->n_basis > 0) {
	  sym = GetSymmetry(lev1->pj);
	  st = (STATE *) ArrayGet(&(sym->states), lev1->pb);
	  c1 = GetConfig(st);
	} else {
	  c1 = GetConfigFromGroup(lev1->iham, lev1->pb);
	}
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
    int myrank = MPIRank(NULL);
    int ipr = 0;
    for (ic0 = 0; ic0 < nic0; ic0++) {
      imax = nc0[ic0];
      jmin = 0;
      lev1 = GetLevel(low[imin]);
      if (lev1->n_basis > 0) {
	sym = GetSymmetry(lev1->pj);
	st = (STATE *) ArrayGet(&(sym->states), lev1->pb);
	c0 = GetConfig(st);
      } else {
	c0 = GetConfigFromGroup(lev1->iham, lev1->pb);
      }
      for (ic1 = 0; ic1 < nic1; ic1++) {
	int skip = SkipMPI();
	jmax = nc1[ic1];
	if (skip) {
	  jmin = jmax;
	  continue;
	}
	lev2 = GetLevel(up[jmin]);	
	if (lev2->n_basis > 0) {
	  sym = GetSymmetry(lev2->pj);
	  st = (STATE *) ArrayGet(&(sym->states), lev2->pb);
	  c1 = GetConfig(st);
	} else {
	  c1 = GetConfigFromGroup(lev2->iham, lev2->pb);
	}
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
	    if (IsPreloadedTR(low[i], up[j], m)) {
	      rd[ir].r.lower = -nut;
	      rd[ir].r.upper = -nut;
	      ir++;
	      continue;
	    }
	    if (fn == NULL) {
	      if (low[i] >= 0 && up[j] >= 0) {
		rd[ir].r.lower = -nut;
		rd[ir].r.upper = -nut;
		ir++;
		continue;
	      }
	    }
	    k = TRMultipoleUTA(&gf, &(rd[ir].rx), m, low[i], up[j],
			       rd[ir].ks, &rd[ir].k0, &rd[ir].k1);
	    if (k == 0) {
	      rd[ir].r.lower = -nut;
	      rd[ir].r.upper = -nut;
	      ir++;
	      continue;
	    }
	    ir0 = ir;
	    rd[ir].r.lower = low[i];
	    rd[ir].r.upper = up[j];
	    rd[ir].r.strength = gf;
	    if (k == -1 && lev1->n_basis == 0 && lev2->n_basis == 0) {
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
	if (wp > 0 && wm > 0.0 && ir0 >= 0 && c0->nr < 1 && c1->nr < 1) {
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
	  if (e0 != ep) {
	    cm = 1.0 + de/(e0 - ep);
	  } else {
	    cm = 1.0;
	  }
	  if (e0 != em) {
	    cp = 1.0 + de/(e0 - em);
	  } else {
	    cp = 1.0;
	  }
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
	      if (rd[ir].r.lower == -1 && rd[ir].r.upper == -1) {
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
	  if (rd[ir].r.lower == -nut && rd[ir].r.upper == -nut) {
	    continue;
	  }
	  if (fn) {
	    if (_sfu_nmax > 0 && rd[ir].r.lower >= 0 && rd[ir].r.upper >= 0) {
	      if (_sfu_ke < nele) _sfu_ke = nele;
	      lev1 = GetLevel(rd[ir].r.lower);
	      lev2 = GetLevel(rd[ir].r.upper);
	      if (lev1->n_basis > 0 && lev2->n_basis > 0 &&
		  rd[ir].k0 >= 0 && rd[ir].k1 >= 0) {
		sym = GetSymmetry(lev2->pj);
		st = (STATE *) ArrayGet(&(sym->states), lev2->pb);
		c1 = GetConfig(st);
		ORBITAL *orb = GetOrbital(rd[ir].k1);
		int kl0, kl1;
		kl1 = GetLFromKappa(orb->kappa);
		kl0 = GetLFromKappa(c1->shells[0].kappa);
		if (orb->n == c1->shells[0].n && kl0 == kl1) {
		  w0 = rd[ir].r.strength;
		  e0 = lev2->energy - lev1->energy;
		  AddSFU(0, rd[ir].k0, rd[ir].k1, nele, e0, w0);
		}
	      }
	    }
	    WriteTRRecord(f, &(rd[ir].r), &(rd[ir].rx));
	  } else {
	    int iuf, iuu;
	    if (rd[ir].r.lower < 0 && rd[ir].r.upper >= 0) {
	      iuf = -(rd[ir].r.lower+1);
	      iuu = rd[ir].r.upper;
	    } else if (rd[ir].r.lower >= 0 && rd[ir].r.upper < 0) {
	      iuf = -(rd[ir].r.upper+1);
	      iuu = rd[ir].r.lower;
	    } else {
	      iuu = -1;
	      iuf = -1;
	    }
	    if (iuu >= 0 && iuf >= 0) {
	      iuu -= _uumin;
	      iuf -= _ufmin;
	      int u = iuu*_nuf_tr + iuf;
	      _ufuu_tr[u].strength = rd[ir].r.strength;
	      _ufuu_tr[u].energy = rd[ir].rx.energy;
	      _ufuu_tr[u].sdev = rd[ir].rx.sdev;
	      _ufuu_tr[u].sci = rd[ir].rx.sci;
	      _ufuu_tr[u].k0 = rd[ir].k0;
	      _ufuu_tr[u].k1 = rd[ir].k1;
	    } else if (_sfu_nmax > 0 &&
		       rd[ir].r.lower < 0 &&
		       rd[ir].r.upper < 0 &&
		       rd[ir].k0 >= 0 &&
		       rd[ir].k1 >= 0) {
	      if (_sfu_ke < nele) _sfu_ke = nele;
	      lev2 = GetLevel(rd[ir].r.upper);
	      c1 = GetConfigFromGroup(lev2->iham, lev2->pb);	      
	      ORBITAL *orb = GetOrbital(rd[ir].k1);
	      int kl1, kl0;
	      kl1 = GetLFromKappa(orb->kappa);
	      kl0 = GetLFromKappa(c1->shells[0].kappa);
	      if (c1->shells[0].n == orb->n && kl0 == kl1) {
		w0 = rd[ir].r.strength*rd[ir].rx.sci;
		AddSFU(1, rd[ir].k0, rd[ir].k1, nele, rd[ir].rx.energy, w0);
	      }
	    }
	  }
	  if (ntrans) {
	    ntrans[myrank]++;
	    if (_progress_report > 0) {
	      if (myrank == 0 && ntrans[0]%_progress_report == 0) {
		PrintTransReport(nproc, tstart, ntrans, "TR", ipr++);
	      }
	    } else if (_progress_report < 0) {
	      if (ntrans[myrank]%(-_progress_report) == 0) {
		PrintTransReport(-myrank, tstart, ntrans, "TR", ipr++);
	      }
	    }
	  }
	}
	free(rd);
	jmin = jmax;
      }
      imin = imax;
    }    
    }
    free(nc0);
    if (up != low) free(nc1);
  } else {
    //PrepAngZStates(nlow, low, nup, up);
    ResetWidMPI();
#pragma omp parallel default(shared) private(a, s, et, j, jup, trd, i, k, gf, r, k0, k1)
    {
      int myrank = MPIRank(NULL);
      int ipr = 0;
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
	  s[i] = 0.0;
	  if (IsPreloadedTR(low[i], up[j], m)) continue;
	  k = TRMultipole(s+i, et+i, m, low[i], up[j], &k0, &k1);
	  if (k != 0) continue;
	  gf = OscillatorStrength(m, et[i], s[i], &(a[i]));
	  a[i] /= jup+1.0;
	  trd += a[i];
	  if (ntrans) {
	    ntrans[myrank]++;
	    if (_progress_report > 0) {
	      if (myrank == 0 && ntrans[0]%_progress_report == 0) {
		PrintTransReport(nproc, tstart, ntrans, "TR", ipr++);
	      }
	    } else if (_progress_report < 0) {
	      if (ntrans[myrank]%(-_progress_report) == 0) {
		PrintTransReport(-myrank, tstart, ntrans, "TR", ipr++);
	      }
	    }
	  }
	}
	if (trd < 1E-30) continue;
	r.upper = up[j];
	for (i = 0; i < nlow; i++) {
	  if (!s[i]) continue;
	  if (!_tr_all && a[i] < (transition_option.eps * trd)) continue;
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
  if (fn) {
    DeinitFile(f, &fhdr);
    CloseFile(f, &fhdr);
  }
  if (_progress_report != 0) {
    PrintTransReport(nproc, tstart, ntrans, "TR", -1);
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

int SaveTransition1(int nlow, int *low, int nup, int *up,
		    char *fn, int m) {
  int nc;

  nc = OverlapLowUp(nlow, low, nup, up);
  if (nc > 0) {
    SaveTransition0(nc, low+nlow-nc, nc, up+nup-nc, fn, m);
    SaveTransition0(nc, low+nlow-nc, nup-nc, up, fn, m);    
    SaveTransition0(nup-nc, up, nc, low+nlow-nc, fn, m);
  }
  SaveTransition0(nlow-nc, low, nup, up, fn, m);
  SaveTransition0(nup, up, nlow-nc, low, fn, m);
  
  return 0;
}
   
int SaveTransition(int nlow, int *low, int nup, int *up,
		   char *fn, int m) {
  int n, *alev, i, r;  
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
  
  if (!IsUTA()) {
    r = SaveTransition1(nlow, low, nup, up, fn, m);
  } else {
    int *ufu0, nuf0, nuu0, ufmin0, ufmax0, uumin0, uumax0;  
    int *ufu1, nuf1, nuu1, ufmin1, ufmax1, uumin1, uumax1;  
    int nufu0, nufu1;
    nufu0 = GetULevelIndices(nlow, low, &ufu0, &nuf0, &nuu0,
			     &ufmin0, &ufmax0, &uumin0, &uumax0);
    if (up == low) {
      ufu1 = ufu0;
      nuf1 = nuf0;
      nuu1 = nuu0;
      ufmin1 = ufmin0;
      ufmax1 = ufmax0;
      uumin1 = uumin0;
      uumax1 = uumax0;
    } else {
      nufu1 = GetULevelIndices(nup, up, &ufu1, &nuf1, &nuu1,
			       &ufmin1, &ufmax1, &uumin1, &uumax1);
    }
    if (nuf0 > 0 && nuf1 > 0) {
      r = SaveTransition1(nuf0, ufu0, nuf1, ufu1, NULL, m);
    }
    if ((nuf0 == 0 && nuf1 == 0) || (nuu0 == 0 && nuu1 == 0)) {
      if (nufu1 > 0 && ufu1 != ufu0) free(ufu1);
      if (nufu0 > 0) free(ufu0);
      r = SaveTransition1(nlow, low, nup, up, fn, m);
    } else {
      _ufmin = GetNumULevels();
      _ufmax = 0;
      _uumin = GetNumLevels();
      _uumax = 0;
      if (ufmax0 >= ufmin0) {
	_ufmin = ufmin0;
	_ufmax = ufmax0;
      }
      if (ufmax1 >= ufmin1) {
	_ufmin = Min(_ufmin, ufmin1);
	_ufmax = Max(_ufmax, ufmax1);
      }
      if (uumax0 >= uumin0) {
	_uumin = uumin0;
	_uumax = uumax0;
      }
      if (uumax1 >= uumin1) {
	_uumin = Min(_uumin, uumin1);
	_uumax = Max(_uumax, uumax1);
      }
      _nuf_tr = _ufmax-_ufmin+1;
      _nuu_tr = _uumax-_uumin+1;
      _ufuu_tr = malloc(sizeof(TR_UTA)*_nuf_tr*_nuu_tr);
      int i;
      for (i = 0; i < _nuf_tr*_nuu_tr; i++) {
	_ufuu_tr[i].strength = 0;
	_ufuu_tr[i].energy = 0;
	_ufuu_tr[i].sdev = 0;
	_ufuu_tr[i].sci = 0;
      }
      r = SaveTransition1(nuf0, ufu0, nuu1, ufu1+nuf1, NULL, m);
      r = SaveTransition1(nuu0, ufu0+nuf0, nuf1, ufu1, NULL, m);
      if (nufu1 > 0 && ufu1 != ufu0) free(ufu1);
      if (nufu0 > 0) free(ufu0);
      r = SaveTransition1(nlow, low, nup, up, fn, m);
      _ufmin = 0;
      _ufmax = -1;
      _nuf_tr = 0;
      _uumin = 0;
      _uumax = -1;
      _nuu_tr = 0;
      free(_ufuu_tr);
      _ufuu_tr = NULL;
    }
  }
  if (n > 0) free(alev);

  ReinitRadial(1);

  return r;
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

void SetNMaxSFU(int n) {
  if (_sfu_nmax > 0) free(_efu);
  _sfu_nmax = n;
  int m = n*n*6;
  _efu = malloc(sizeof(double)*m);
  int i;
  for (i = 0; i < m; i++) _efu[i] = 0.0;
}

void PrintSFU(char *fn) {
  FILE *f;

  f = fopen(fn, "w");
  if (f == NULL) {
    printf("cannot open file: %s\n", fn);
    return;
  }

  int i, j, k, nm1, nm2, nm3, nm4, nm5;
  double e0, e1;
  fprintf(f, "# %3d %3d\n", _sfu_ke, _sfu_nmax);
  nm1 = _sfu_nmax*_sfu_nmax;
  nm2 = nm1*2;
  nm3 = nm1*3;
  nm4 = nm1*4;
  nm5 = nm1*5;
  for (i = 0; i < _sfu_nmax; i++) {
    for (j = 0; j < _sfu_nmax; j++) {
      k = j*_sfu_nmax + i;
      e0 = 0.0;
      e1 = 0.0;
      if (_efu[k+nm2] > 0) {
	e0 = (_efu[k]/_efu[k+nm2])*HARTREE_EV;
      }
      if (_efu[k+nm3] > 0) {
	e1 = (_efu[k+nm1]/_efu[k+nm3])*HARTREE_EV;
      }
      fprintf(f, "%5d %3d %12.5E %12.5E %12.5E %12.5E %12.5E %8d %8d\n",
	      i+1, j+1, e0, e1,
	      _efu[k+nm2], _efu[k+nm3], e0-e1,
	      (int)(_efu[k+nm4]+0.1), (int)(_efu[k+nm5]+0.1));
    }
  }
  fclose(f);
}

void SetOptionTransition(char *s, char *sp, int ip, double dp) {
  if (0 == strcmp(s, "transition:aw")) {
    _tr_aw = dp;
    return;
  }
  if (0 == strcmp(s, "transition:ls_all")) {
    _tr_all = ip;
    return;
  }
  if (0 == strcmp(s, "transition:progress_report")) {
    _progress_report = ip;
    return;
  }
  if (0 == strcmp(s, "transition:lower_emin")) {
    _lower_emin = dp/HARTREE_EV;
    return;
  }
  if (0 == strcmp(s, "transition:lower_emax")) {
    _lower_emax = dp/HARTREE_EV;
    return;
  }
  if (0 == strcmp(s, "transition:upper_emin")) {
    _upper_emin = dp/HARTREE_EV;
    return;
  }
  if (0 == strcmp(s, "transition:upper_emax")) {
    _upper_emax = dp/HARTREE_EV;
    return;
  }
  if (0 == strcmp(s, "transition:tr_emin")) {
    _tr_emin = dp/HARTREE_EV;
    return;
  }
  if (0 == strcmp(s, "transition:tr_emax")) {
    _tr_emax = dp/HARTREE_EV;
    return;
  }
  if (0 == strcmp(s, "transition:uta_e1")) {
    _uta_e1 = ip;
    return;
  }
  if (0 == strcmp(s, "transition:sfu_nmax")) {
    SetNMaxSFU(ip);
    return;
  }
  if (0 == strcmp(s, "transition:sfu_ofn")) {
    PrintSFU(sp);
    return;
  }
}
