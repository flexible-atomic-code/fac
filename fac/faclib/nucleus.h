#ifndef _NUCLEUS_H_
#define _NUCLEUS_H_

#include "global.h"

typedef struct _NUCLEUS_ {
  char symbol[5];
  double atomic_number;
  double mass;
  double rn;
} NUCLEUS;


int SetAtom(char *s, double z, double mass);
char *GetAtomicSymbolTable(void);
double *GetAtomicMassTable(void);
double GetAtomicNumber(void);
double GetAtomicMass(void);
double GetAtomicR(void);
char *GetAtomicSymbol(void);
double GetAtomicEffectiveZ(double r);

#endif

