#ifndef _POLARIZATION_H_
#define _POLARIZATION_H_

typedef struct _MLEVEL_ {
  short j;
  short p;
  double energy;
  double dtotal;
  int ic;
} MLEVEL;

typedef struct _MTR_ {
  int multipole;
  int lower;
  int upper;
  int n;
  double rtotal;
  double *rates;
} MTR;

typedef struct _MCE_ {
  int lower;
  int upper;
  int n;
  double *rates;
} MCE;

int InitPolarization(void);
int SetMLevels(char *fn, char *tfn);
int SetMCERates(char *fn, double energy);
int PopulationTable(char *fn, double eden);
int PolarizationTable(char *fn);

#endif
