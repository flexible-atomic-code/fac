#include "transition.h"
#include <time.h>

static char *rcsid="$Id: transition.c,v 1.5 2001/09/14 16:32:06 mfgu Exp $";
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
}

double GetTransitionCut() {
  return transition_option.eps;
}

void SetTransitionOptions(int gauge, int mode, 
			  int max_e, int max_m) {
  transition_option.gauge = gauge;
  transition_option.mode = mode;
  transition_option.max_e = max_e;
  transition_option.max_m = max_m;
}

int GetTransitionGauge() {
  return transition_option.gauge;
}

int OscillatorStrength(double *strength, double *energy, 
			  int *multipole, int lower, int upper) {
  int m, m2, n;
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

  m = GetLowestMultipole(p1, j1, p2, j2);
  if (m > transition_option.max_m || m < -transition_option.max_e) return 1;
  if (*multipole != 0) {
    if (m != *multipole) {
      if (m == 1 && *multipole == -2) m = -2;
      else return -1;
    }
  } else {
    *multipole = m;
  }

  aw = FINE_STRUCTURE_CONST * (*energy);
  if (aw < 0.0) return -1;

  m2 = 2*abs(m);
  s = 0.0;
  
  nz = AngularZMix(&ang, lower, upper, m2, m2);
  if (nz <= 0) return -1;
  for (i = 0; i < nz; i++) {
    if (fabs(ang[i].coeff) > EPS10) {
      if (transition_option.mode) {
	r = MultipoleRadialNR(m, ang[i].k0, ang[i].k1, 
			      transition_option.gauge);
      } else {
	r = MultipoleRadial(aw, m, ang[i].k0, ang[i].k1,
			    transition_option.gauge);
      }
      s += r * ang[i].coeff;
    }
  }
 
  free(ang);	  
  
  *strength = s*s/(m2+1.0);
  *strength *= (*energy) * pow(aw, m2-2);
#if FAC_DEBUG
  fprintf(debug_log, "**%12.8lE %12.8lE \n", s, *strength);
#endif
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
		   char *fn, int multipole) {
  int i, j, k, n, m, jup, jlow;
  FILE *f;
  double *s, *et, *a, trd, trd1, elow, e0;
  char t;
  int *alev;

  f = fopen(fn, "w");
  if (!f) return -1;

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

  fprintf(f, "up     2J\tlow     2J\tDelta_E    M  gf        A(AU)\n");
  if (nlow <= 0 || nup <= 0) return -1;
  e0 = GetLevel(0)->energy;
  a = malloc(sizeof(double)*nlow);
  s = malloc(sizeof(double)*nlow);
  et = malloc(sizeof(double)*nlow);
  for (j = 0; j < nup; j++) {
    jup = LevelTotalJ(up[j]);
    trd = 0.0;
    trd1 = 0.0;
    for (i = 0; i < nlow; i++) {
      a[i] = 0.0;
      elow = GetLevel(low[i])->energy;
      m = multipole;
      et[i] = 0.0;
      k = OscillatorStrength(s+i, et+i, &m, low[i], up[j]);
      if (k != 0) continue;
      if (m == 0) continue;
      if (s[i] < 1E-30) continue;
      a[i] = 2*pow((FINE_STRUCTURE_CONST*et[i]),2)*FINE_STRUCTURE_CONST;
      a[i] *= s[i]/(jup+1.0);
      trd += a[i];
      if (elow < e0) trd1 += a[i];
    } 
    if (trd < 1E-30) continue;
    for (i = 0; i < nlow; i++) {
      if (a[i] < (transition_option.eps * trd)) continue;
      jlow = LevelTotalJ(low[i]);
      elow = GetLevel(low[i])->energy;
      if (elow < e0) k = -low[i];
      else k = low[i];
      et[i] *= HARTREE_EV;
      if (m < 0) t = 'E';
      else t = 'M';
      fprintf(f, "%-6d %-2d\t%-7d %-2d\t%10.4E %c%d %9.3E %9.3E\n", 
	      up[j], jup, k, jlow, et[i], t, abs(m), s[i], a[i]);
    }	
    fprintf(f, "%-6d   \tTotal     \t              %9.3E %9.3E\n\n", 
	    up[j], trd, trd1); 
  }
  free(a);
  free(s);
  free(et);
  fclose(f);
  if (n > 0) {
    free(alev);
  }
  return 0;
}
