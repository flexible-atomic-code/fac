#include "nucleus.h"

static NUCLEUS atom;

int SetAtom(char *s, double z, double mass) {
  if (s == NULL) return -1;
  strncpy(atom.symbol, s, 2); 
  atom.atomic_number = z;
  if (mass <= 0.0) atom.mass = 2.0*z;
  atom.rn = 2.2677E-5 * pow(mass, 1.0/3);
}

double GetAtomicNumber() {
  return atom.atomic_number;
}

char *GetAtomSymbol() {
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
