#include "transition.h"
#include <time.h>

static char *rcsid="$Id: transition.c,v 1.13 2002/02/28 16:55:05 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

/* the following options controll the computation methods of OS.
   gauge = 1 coulomb gauge (velocity form)
         = 2 babushkin gauge (length form)

   mode = 0 use relativistic expression for radial integrals.
        = 1 use non-resltivistic approximation.
   
   max_e, the maximum rank of electric multipoles.
   max_m, the maximum rank of magnetic multipoles.
*/
static struct {
  int gauge;
  int mode;
  int max_e;
  int max_m;
  double eps;
} transition_option = {2, 1, 4, 4, EPS4};

int SetTransitionCut(double c) {
  transition_option.eps = c;
  return 0;
}

double GetTransitionCut(void) {
  return transition_option.eps;
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

int OscillatorStrength(double *strength, double *energy,
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

  if (m > transition_option.max_m || m < -transition_option.max_e) return 1;
  if (!Triangle(j1, j2, m2)) return -1;
  if (m > 0 && IsEven(p1+p2+m)) return -1;
  if (m < 0 && IsOdd(p1+p2-m)) return -1;

  aw = FINE_STRUCTURE_CONST * (*energy);
  if (aw < 0.0) return -1;

  s = 0.0;
  
  nz = AngularZMix(&ang, lower, upper, m2, m2);
  if (nz <= 0) return -1;
  for (i = 0; i < nz; i++) {
    if (transition_option.mode == M_NR && 
	!(m == 1 && ang[i].k0 != ang[i].k1)) {
      r = MultipoleRadialNR(m, ang[i].k0, ang[i].k1, 
			    transition_option.gauge);
    } else {
      r = MultipoleRadialFR(aw, m, ang[i].k0, ang[i].k1,
			    transition_option.gauge);
    }
    s += r * ang[i].coeff;
  }
 
  free(ang);	  
  
  *strength = s*s/(m2+1.0);
  *strength *= (*energy);
  m2 = m2 - 2;
  if (transition_option.gauge == G_COULOMB && 
      transition_option.mode == M_FR && m < 0) 
    m2 = m2 - 2; 

  if (m2) {
    (*strength) *= pow(aw, m2);
  }

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


int SaveTransition(int nlow, int *low, int nup, int *up, 
		   char *fn, int m) {
  int i, j, k, n, jup;
  FILE *f;
  LEVEL *lev1, *lev2;
  TR_RECORD r;
  TR_HEADER tr_hdr;
  F_HEADER fhdr;
  double *s, *et, *a, trd;
  double e0, emin, emax;
  int *alev;
#ifdef PERFORM_STATISTICS
  ARRAY_TIMING arrayt;
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
  tr_hdr.mode = GetTransitionMode();
  f = OpenFile(fn, &fhdr);
  InitFile(f, &fhdr, &tr_hdr);

  a = malloc(sizeof(double)*nlow);
  s = malloc(sizeof(double)*nlow);
  et = malloc(sizeof(double)*nlow);
  for (j = 0; j < nup; j++) {
    jup = LevelTotalJ(up[j]);
    trd = 0.0;
    for (i = 0; i < nlow; i++) {
      a[i] = 0.0;
      et[i] = 0.0;
      k = OscillatorStrength(s+i, et+i, m, low[i], up[j]);
      if (k != 0) continue;
      if (s[i] < 1E-30) continue;
      a[i] = 2*pow((FINE_STRUCTURE_CONST*et[i]),2)*FINE_STRUCTURE_CONST;
      a[i] *= s[i]/(jup+1.0);
      trd += a[i];
    } 
    if (trd < 1E-30) continue;
    r.upper = up[j];
    for (i = 0; i < nlow; i++) {
      if (a[i] < (transition_option.eps * trd)) continue;
      r.lower = low[i];
      r.strength = s[i];
      WriteTRRecord(f, &r);
    }
  }

  DeinitFile(f, &fhdr);
  CloseFile(f, &fhdr);

  free(a);
  free(s);
  free(et);
  if (n > 0) {
    free(alev);
  }

#ifdef PERFORM_STATISTICS
  GetArrayTiming(&arrayt);
  fprintf(perform_log, "Multi: %6.1E, Array: %6.1E\n",
	  ((double) (arrayt.multi))/CLOCKS_PER_SEC,
	  ((double) (arrayt.array))/CLOCKS_PER_SEC);
  GetStructTiming(&structt);
  fprintf(perform_log, "AngZMix: %6.1E, AngZFB: %6.1E, AngZxZFB: %6.1E, SetH: %6.1E DiagH: %6.1E\n",
	  ((double) (structt.angz_mix))/CLOCKS_PER_SEC,
	  ((double) (structt.angz_fb))/CLOCKS_PER_SEC,
	  ((double) (structt.angzxz_fb))/CLOCKS_PER_SEC,
	  ((double) (structt.set_ham))/CLOCKS_PER_SEC,
	  ((double) (structt.diag_ham))/CLOCKS_PER_SEC);
  fprintf(perform_log, "AngZS: %6.1E, %ld %ld, AngZFBS: %6.1E, AngZxZFBS: %6.1E, AddZ: %6.1E, AddZxZ: %6.1E\n",
	  ((double) (structt.angz_states))/CLOCKS_PER_SEC, 
	  structt.angz_states_load, structt.angz_states_calc,
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
