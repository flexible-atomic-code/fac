#ifndef _RMATRIX_H_ 
#define _RMATRIX_H_

#include "global.h"
#include "structure.h"
#include "radial.h"

#define NBTERMS 5
#define NBFIT 6
typedef struct _RBASIS_ {
  int kmax, nbk, nkappa, nbuttle;
  int ib0, ib1;
  double rb0, rb1, bqp;
  int **basis;
  double **ebuttle, **cbuttle[NBTERMS];
  double **ek, **w0, **w1;
} RBASIS;

typedef struct _RMATRIX_ {
  int nts, nkappa, nchan, nchan0;
  int *chans, *ts;
  int ndim;
  int isym, p, j;
  double et0, *et, *ek, **w0, **w1;
  double *rmatrix[3];
  double energy;
} RMATRIX;

int InitRMatrix(void);
void ClearRMatrixBasis(void);
void ReadRMtraixBasis(char *fn);
void WriteRMatrixBasis(char *fn);
void RMatrixBoundary(double r0, double r1, double b);
int RMatrixBasis(char *fn, int kmax, int nb);
int IndexFromKappa(int k);
int KappaFromIndex(int i);
void RMatrixTargets(int nt, int *kt, int nc, int *kc);
int ReadRMatrixSurface(FILE *f, RMATRIX *rmx, int m);
void WriteRMatrixSurface(FILE *f, double **wik0, double **wik1, int m);
int RMatrixSurface(char *fn);
int RMatrix(double e, RMATRIX *rmx, int m);
void TestRMatrix(double e, int m, char *fn1, char *fn2, char *fn3);

#endif
