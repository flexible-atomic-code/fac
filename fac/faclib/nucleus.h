#ifndef _NUCLEUS_H_
#define _NUCLEUS_H_

#include "global.h"

#define N_ELEMENTS 109

typedef struct _NUCLEUS_ {
  char symbol[5];
  double atomic_number;
  double mass;
  double rn;
} NUCLEUS;


int SetAtom(char *s, double z, double mass);
char *GetAtomicSymbolTable();
double *GetAtomicMassTable();
double GetAtomicNumber();
double GetAtomicMass();
char *GetAtomicSymbol();
double GetAtomicEffectiveZ(double r);

#endif

