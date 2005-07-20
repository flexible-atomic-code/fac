#ifndef _MBPT_H_
#define _MBPT_H_

#include "structure.h"


typedef struct _CORR_CONFIG_ {
  CONFIG *c;
  int inp, inq, np, nq;
  int ig, ncs, kp, kq;
} CORR_CONFIG;

typedef struct _MBPT_HAM_ {
  int n, n2, n3;
  int *ng, *ng2;
  int isym, dim;
  int ibra, iket;
  double a, b, c;
  double *hab1, *hba1;
  double *hab, *hba;
} MBPT_HAM;

typedef struct _MBPT_EFF_ {
  int nbasis, *basis;
  double *h0, *heff;
  double **hab1, **hba1;
  double **hab, **hba;
} MBPT_EFF;

int StructureMBPT(char *fn, char *fn1, int n, int *s0, int kmax, 
		  int n1, int *nm, int n2, int *nmp, int n0, char *gn0);
int StructureReadMBPT(char *fn, char *fn2, int nf, char *fn1[], 
		      int nkg, int *kg, int nkg0);
void SetExtraMBPT(int m);
void SetOptMBPT(int n3, double c);
void SetSymMBPT(int p, int j, int i);

#endif
