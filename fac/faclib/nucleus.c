#include "nucleus.h"

static char *rcsid="$Id: nucleus.c,v 1.13 2004/03/11 00:26:05 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

static NUCLEUS atom;
static char _ename[N_ELEMENTS][3] = 
{"H", "He", "Li", "Be", "B", "C", "N", "O", "F",
 "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", 
 "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", 
 "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
 "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", 
 "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
 "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", 
 "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf",
 "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl",
 "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
 "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es",
 "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt"};

static double _emass[N_ELEMENTS] = 
{1, 4, 7, 9, 11, 12, 14, 16, 19, 20, 23, 24, 27, 28, 31, 32, 35, 40, 39, 
 40, 45, 48, 51, 52, 55, 56, 58, 59, 64, 65, 70, 73, 75, 79, 80, 84, 85,
 88, 89, 91, 93, 96, 98, 101, 103, 106, 108, 112, 115, 119, 122, 128, 127,
 131, 133, 137, 139, 140, 141, 144, 147, 150, 152, 157, 159, 162, 165, 167,
 169, 173, 175, 178, 181, 184, 186, 190, 190, 195, 197, 200, 204, 207, 209, 
 210, 210, 222, 223, 226, 227, 232, 231, 238, 237, 242, 243, 247, 247, 249, 
 254, 253, 256, 254, 257, 257, 260, 263, 262, 265, 266};


char *GetAtomicSymbolTable(void) {
  return (char *) _ename;
}

double *GetAtomicMassTable(void) {
  return _emass;
}

int SetAtom(char *s, double z, double mass, double rn) {
  int i;

  if (s == NULL) return -1;
  strncpy(atom.symbol, s, 2); 
  if (z <= 0 || mass <= 0) {
    for (i = 0; i < N_ELEMENTS; i++) {
      if (strncasecmp(_ename[i], s, 2) == 0) {
	if (z <= 0) atom.atomic_number = i+1;
	if (mass <= 0) atom.mass = _emass[i];
	break;
      }
    }
    if (i == N_ELEMENTS) return -1;
  }

  if (z > 0) {
    atom.atomic_number = z;
  } 
  if (mass > 0.0) {
    atom.mass = mass;
  }
  if (rn < 0.0) {
    atom.rn = 2.2677E-5 * pow(atom.mass, 1.0/3);
  } else {
    atom.rn = rn;
  }

  return 0;
}

double GetAtomicMass(void) {
  return atom.mass;
}

double GetAtomicNumber(void) {
  return atom.atomic_number;
}

char *GetAtomicSymbol(void) {
  return atom.symbol;
}

double GetAtomicEffectiveZ(double r) {
  double x, y;
  if (r > atom.rn) {
    return (double) atom.atomic_number;
  } else {
    x = r/atom.rn;
    y = 3.0 - x*x;
    y = x*y*0.5*(atom.atomic_number);
    return y;
  }
}

double GetAtomicR(void) {
  return atom.rn;
}
